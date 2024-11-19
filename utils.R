#!/usr/bin/env Rscript

# opens all metadata tables
vars <- c()

for (file in list.files("data/sra-metadata", full.names = TRUE)) {
  filename <- stringr::str_extract(file, "(?<=metadata\\/).+(?=.txt)")
  vars <- append(vars, filename)
  assign(x = filename, value = readr::read_csv(file))
}


accessions <- purrr::map(vars,
                         ~  {
                           data <- eval(parse(text = .x))

                           runs <- data |> dplyr::select(BioProject, Organism, Run, LibraryLayout)

                           if ("genotype" %in% colnames(data)) {
                             runs$genotype <- data$genotype
                           }
                           if ("strain" %in% colnames(data)) {
                             runs$genotype <- data$strain
                           }

                           return(runs)
                         }
) |> purrr::list_rbind()

kt2440 <- accessions |>
  dplyr::filter(dplyr::if_else(is.na(genotype),
                               stringr::str_ends(Organism, "KT2440"),
                               genotype %in% c("KT2440", "wild type", "KT2440 wildtype")))

suffixes <- c(".fastq","_1.fastq", "_2.fastq")

# Defines a function to trim reads using Trimmomatic

getTrimmomaticCall <- function(BioProject, Run, LibraryLayout, outdir) {

  message(glue::glue("Generating Trimmomatic call for {Run} in Bioproject {BioProject}\n The library layout is {LibraryLayout}"))


  if (LibraryLayout == "PAIRED") {

    expression <- glue::glue(
      'java -jar software/Trimmomatic-0.39/trimmomatic-0.39.jar PE \\
                  	-threads 8 \\
                  	-trimlog {outdir}/trimlog/{Run}_trimlog.txt \\
                  	-baseout {outdir}/{Run}_trim.fastq \\
                  	data/{BioProject}/rawdata/{Run}_1.fastq data/{BioProject}/rawdata/{Run}_2.fastq \\
                  	ILLUMINACLIP:software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \\
                  	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50')
  } else {

    if (LibraryLayout == "SINGLE") {

      expression <- glue::glue('
                    java -jar software/Trimmomatic-0.39/trimmomatic-0.39.jar SE \\
                    -threads 8 \\
                    -trimlog {outdir}/trimlog/{Run}_trimlog.txt \\
                    data/{BioProject}/rawdata/{Run}.fastq \\
                    {outdir}/{Run}_trim.fastq \\
                    ILLUMINACLIP:software/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:6 \\
                    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:40')

    } else {
      message("Library Layout not detected. Skipping...")
      break()
    }
  }

  return(expression)
}


# Defines a function to run RSEM
getRSEMCall <- function(BioProject, Run, LibraryLayout, outDir, rsemReference) {

  if (LibraryLayout == "SINGLE") {

    message("running rsem-calculate-expression in single end mode")

    call <- glue::glue("rsem-calculate-expression --bowtie2 --no-bam-output --num-threads 8 \\
    data/{BioProject}/trimmed/{Run}_trim.fastq \\
    {rsemReference} \\
    {outDir}/{Run}")

  } else {
    if (LibraryLayout == "PAIRED") {

      message("running rsem-calculate-expression in paired end mode")

      call <- glue::glue(
      "rsem-calculate-expression --bowtie2 --no-bam-output --num-threads 8 --paired-end \\
       data/{BioProject}/trimmed/{Run}_trim_1P.fastq data/{BioProject}/trimmed/{Run}_trim_2P.fastq \\
       {rsemReference} \\
       {outDir}/{Run}"
      )

    } else {
      message("Library layout not detected!")
      break()
      }
  }

  return(call)
}


# Defines helper functions for dealing with DESeq2


prepareMetadata <- function(metadataTable, makeConditionsFrom = NULL, na.value = "") {

  if(is.null(makeConditionsFrom)) {
    message("No columns to collapse for conditions were chosen!")
    message("available columns in metadata file are:")
    message(stringr::str_flatten_comma(colnames(metadataTable), last = " and "))

    return(invisible())
  }

  out <- metadataTable |>
    dplyr::select(Run, dplyr::all_of(makeConditionsFrom)) |>
    dplyr::mutate(dplyr::across(dplyr::all_of(makeConditionsFrom), ~ tidyr::replace_na(.x,replace = na.value))) |>
    tidyr::unite(col = Condition, dplyr::all_of(makeConditionsFrom),sep = "_") |>
    dplyr::mutate(Condition = stringr::str_replace_all(Condition, pattern = "\\+", replacement = ""),
                  Condition = stringr::str_squish(Condition),
                  Condition = stringr::str_replace_all(Condition, pattern = " ", replacement = "_"),
                  Condition = stringr::str_replace_all(Condition, pattern = ",", replacement = "."),
                  Condition = stringr::str_to_lower(Condition),
                  Condition = forcats::as_factor(Condition))


  message("the new Condition column values are:")
  message(paste0(unique(out$Condition), "\n"))
  message("Use these for specifying levels with conditionLevelOrder()")

  return(out)
}

`conditionOrder<-` <- function(preparedMetadata, value) {
  levels(preparedMetadata[["Condition"]]) <- value

  preparedMetadata
}

runTxi <- function(metadata, dataset, prefix = "rsem/genbank") {

  md <- metadata

  files <- c()
  dl_files <- list.files(glue::glue("{prefix}/{dataset}"), full.names = TRUE, pattern = glue::glue(".genes.results"))

  for (run in unique(md$Run)) {
    filename <- glue::glue("{prefix}/{dataset}/{run}.genes.results")
    if( filename %in% dl_files) {
      names(filename) <- run
      files <- append(files, filename)
    }
  }

  txi.rsem <- tximport::tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)


  # identifies if there are transcripts with length == 0
  zeroLengthTranscripts <- which(txi.rsem$length == 0, arr.ind = TRUE)

  # if so, filter them out as suggested by michael love (https://support.bioconductor.org/p/132686/)
  if (length(zeroLengthTranscripts) != 0) {

    message("zero-length transcripts were detected in the dataset. Filtering...")
    transcriptNames <- rownames(zeroLengthTranscripts)

    # updates the matrices in txis.rsem to remove the culprits
    txi.rsem$abundance <- txi.rsem$abundance[-which(rownames(txi.rsem$abundance) %in% transcriptNames), ]
    txi.rsem$counts <- txi.rsem$counts[-which(rownames(txi.rsem$counts) %in% transcriptNames), ]
    txi.rsem$length <- txi.rsem$length[-which(rownames(txi.rsem$length) %in% transcriptNames), ]
  }

  # identifies if there are samples with an all 0 count value
  zeroCountSamples <- tibble::as_tibble(txi.rsem$counts) |>
    dplyr::summarise(dplyr::across(where(is.numeric), sum)) |>
    tidyr::pivot_longer(cols = everything(),
                        names_to = "Run",
                        values_to = "Total") |>
    dplyr::filter(Total == 0) |>
    dplyr::pull(Run)

  if(length(zeroCountSamples) > 0) {
    message(glue::glue("Runs {stringr::str_flatten_comma(zeroCountSamples, last = ' and ')} are constant, zero-count runs. Removing for DESeq2 compatibility..."))

    colnames <- colnames(txi.rsem$counts)

    txi.rsem$abundance <- txi.rsem$abundance[,!colnames%in% zeroCountSamples]
    txi.rsem$length <- txi.rsem$length[,! colnames %in% zeroCountSamples]
    txi.rsem$counts <- txi.rsem$counts[,!colnames %in% zeroCountSamples]
  }

  return(txi.rsem)

}

getDESeqResults <- function(ddsTxi, dataset, comparisons) {

  out <- NULL

  for (i in 1:nrow(comparisons)) {
    comp <- dplyr::slice(comparisons, i)
    contrastVector <- c(comp$factorName, comp$numeratorLevel, comp$denominatorLevel)

    res <- DESeq2::results(ddsTxi, contrast = contrastVector) |>
      tibble::as_tibble(rownames = NA) |>
      tibble::rownames_to_column(var = "transcript") |>
      dplyr::mutate(dataset = dataset,
                    condition = contrastVector[2],
                    control = contrastVector[3]) |>
      dplyr::relocate(dataset, condition, control, .before = transcript)

    out <- dplyr::bind_rows(out, res)
  }

  return(out)
}

readGEOmetadata <- function(series_matrix_path) {

  output <- readr::read_lines(series_matrix_path) |>
    stringr::str_subset(pattern = "\\!Sample_geo_accession|\\!Sample_description") |>
    stringr::str_remove_all(pattern = '\\"') |>
    stringr::str_remove_all(pattern = "!") |>
    stringr::str_split(pattern = "\\t|' '") |>
    dplyr::bind_cols() |>
    janitor::row_to_names(1)

  return(output)

}

processRNAseq <- function(data, folder.prefix = "rsem/genbank") {
  tempdataname <- data$dataset
  tempmeta <- data$metadata
  order <- data$conditionorder

  conditionOrder(tempmeta) <- order


  # obtaining DESeq2 results...
  txi <- runTxi(tempmeta, dataset = tempdataname, prefix = folder.prefix)

  tempmeta <- tempmeta |> dplyr::filter(Run %in% colnames(txi$counts))

  # calling DESeq2 on the tximport object
  ddsTxi <- DESeq2::DESeqDataSetFromTximport(txi, colData = tempmeta, design = ~ Condition) |>
    DESeq2::DESeq()

  all_comparisons <- data.frame(factorName = "Condition",
                                numeratorLevel = levels(tempmeta$Condition)[-1],
                                denominatorLevel = levels(tempmeta$Condition)[1])
  res <- getDESeqResults(ddsTxi, tempdataname, comparisons = all_comparisons) |>
    dplyr::mutate(Source = folder.prefix)

  return(res)
}
