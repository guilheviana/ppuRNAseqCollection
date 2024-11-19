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
