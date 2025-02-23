source("utils.R")
source("prepare-loci-info.R")

data <- list()
datasets <- kt2440 |> dplyr::pull(BioProject) |> unique()

# PRJNA338603 ----
dataset <- datasets[1]
metadata <- prepareMetadata(get(dataset),makeConditionsFrom = c("stress","time_point"), na.value = "0 min")

order <- c("none_0_min", "nacl_7_min", "nacl_60_min", "imipenem_7_min", "imipenem_60_min", "h2o2_7_min", "h2o2_60_min")

data <- append(data, list(list(dataset = dataset, metadata = metadata, conditionorder = order)))

# PRJNA341451 ----
dataset <- datasets[2]
metadata <- prepareMetadata(get(dataset), makeConditionsFrom = "treatment")
order <- c("glucose", "toluene_vapor_phase", "toluene_shock")

data <- append(data, list(list(dataset = dataset, metadata = metadata, conditionorder = order)))


# PRJNA450701 ----
dataset <- datasets[3]
metadata <- prepareMetadata(get(dataset), makeConditionsFrom = "Library Name") |>
  dplyr::mutate(Condition = stringr::str_extract(Condition, "^.{2}"))
order <- c("ck", "ls", "is", "hs")

data <- append(data, list(list(dataset = dataset, metadata = metadata, conditionorder = order)))


# PRJNA480068 ----
dataset <- datasets[4]
metadata <- prepareMetadata(get(dataset), makeConditionsFrom = "carbon_source")
order <- c("glucose", "myristic_acid")

data <- append(data, list(list(dataset = dataset, metadata = metadata, conditionorder = order)))


# PRJNA486005 ----
dataset <- datasets[5]
metadata <- prepareMetadata(get(dataset), makeConditionsFrom = "carbon_source")|>
  dplyr::mutate(Condition = stringr::str_extract(Condition, "(?<=in_).+"))
order <- c("glucose", "citrate", "ferulic_acid", "serine")

data <- append(data, list(list(dataset = dataset, metadata = metadata, conditionorder = order)))

# PRJNA533248 ---- <- will come back to this one later
dataset <- datasets[6]

series_metadata <- c("metadata/GSE129947_series_matrix.txt.gz") |>
  purrr::map_df(readGEOmetadata) |>
  dplyr::mutate(Condition = stringr::str_remove_all(Sample_description, "_[A-Z]$"),
                Condition = stringr::str_squish(Condition),
                Condition = stringr::str_to_lower(Condition))

metadata <- get(dataset) |>
  dplyr::left_join(series_metadata, by = c( "Sample Name" = "Sample_geo_accession")) |>
  dplyr::select(Run, Condition) |>
  dplyr::mutate(Condition = forcats::fct(Condition))

order <- c("s_0min", "s_15min", "p1_15min", "p3_15min", "p5_15min", "s_1h", "s_3h", "s_9h", "s_25h", "p1_25h", "p3_25h", "p5_25h")

data <- append(data, list(list(dataset = dataset, metadata = metadata, conditionorder = order)))


# PRJNA630234 ----

dataset <- datasets[7]
metadata <- prepareMetadata(get(dataset), makeConditionsFrom = "growth")|>
  dplyr::mutate(
    Condition = stringr::str_remove_all(Condition,"(in_the_)|(supplemented_with_)|(_m9)|(\\d{1}.*g/l_)"),
    Condition = stringr::str_squish(Condition))
  order <- c("glucose", "glucose_triethanolamine_acetate", "glucose_triethylamine_hydrogen_sulfate", "glucose_acetate")
conditionOrder(metadata) <- order

data <- append(data, list(list(dataset = dataset, metadata = metadata, conditionorder = order)))


# PRJNA796354 ----
dataset <- datasets[8]

# since this dataset does not include the sample conditions in the metadata
# we will retrieve them from the series matrix manually downloaded from GEO and located in the
# metadata folder in this project

# it is divided into 3 separate files, which we will group now
series_metadata <- c("metadata/GSE193493-GPL28478_series_matrix.txt.gz",
                     "metadata/GSE193493-GPL28971_series_matrix.txt.gz",
                     "metadata/GSE198395_series_matrix.txt.gz") |>
                  purrr::map(readGEOmetadata) |>
                  purrr::list_rbind() |>
                  dplyr::mutate(Condition = stringr::str_remove_all(Sample_description, "^Pseudomonas.+(with|medium)"),
                                Condition = stringr::str_remove_all(Condition, "supplemented"),
                                Condition = stringr::str_remove_all(Condition, "as\\s.+$"),
                                Condition = stringr::str_remove_all(Condition, "[Rr]eplicate.+$"),
                                Condition = stringr::str_remove_all(Condition, "of a mixture of"),
                                Condition = stringr::str_remove_all(Condition, "2.5 g/L "),
                                Condition = stringr::str_squish(Condition),
                                Condition = stringr::str_to_lower(Condition))

# in this case, we will merge the metadata with the dataset info
metadata <- get(dataset) |>
            dplyr::left_join(series_metadata, by = c( "Sample Name" = "Sample_geo_accession")) |>
            dplyr::select(Run, Condition) |>
            dplyr::mutate(Condition = stringr::str_replace_all(Condition, ", and ", "_"),
                          Condition = stringr::str_replace_all(Condition, " and ", "_"),
                          Condition = stringr::str_replace_all(Condition, ", ", "_"),
                          Condition = forcats::fct(Condition))
order <- c("glucose", "gluconate", "fructose", "ferulate", "coumarate", "glucose_gluconate",  "fructose_glucose", "coumarate_ferulate", "fructose_glucose_gluconate")

data <- append(data, list(list(dataset = dataset, metadata = metadata, conditionorder = order)))


# PRJNA885112 ----
dataset <- datasets[9]
metadata <- prepareMetadata(get(dataset), makeConditionsFrom = "condition") |>
  dplyr::mutate(Condition = stringr::str_replace(Condition, " \\+ ", "_"))
order <- c("lb", "lb_selenite")

data <- append(data, list(list(dataset = dataset, metadata = metadata, conditionorder = order)))

# running deseq2

genbank_rnaseq <- purrr::map(data, ~ processRNAseq(.x, folder.prefix = "rsem/genbank"))

genbank_processed <- genbank_rnaseq |>
  dplyr::bind_rows() |>
  dplyr::left_join(genbank_lookup)

refseq_rnaseq <- purrr::map(data, ~ processRNAseq(.x, folder.prefix = "rsem/refseq"))

refseq_processed <- refseq_rnaseq |>
  dplyr::bind_rows() |>
  dplyr::left_join(refseq_lookup, relationship = "many-to-many")

# saving output files
unique_conditions <- dplyr::bind_rows(genbank_processed, refseq_processed) |>
  dplyr::select(condition) |>
  dplyr::distinct()

dir.create("./output/")
if (!file.exists("./output/unique_conditions.tsv")) readr::write_tsv(unique_conditions, file = "./output/unique_conditions.tsv")

if (file.exists("./output/unique_conditions_MANUALFORMATTED.csv")) {
  # get the formatted conditions table (manually curated)
  cond_format <- readr::read_tsv("output/unique_conditions_MANUALFORMATTED.csv")

  genbank_processed <- genbank_processed |>
    dplyr::left_join(cond_format) |>
    dplyr::relocate(dataset, condition_formatted, locus_tag, baseMean, log2FoldChange,
                    lfcSE, stat, pvalue, padj, condition, control, transcript, Source)

  readr::write_tsv(genbank_processed, "../output/tables/genbank_rnaseq.tsv")
  saveRDS(genbank_processed, "genbank_processed.RDS")

  refseq_processed <- refseq_processed |>
    dplyr::left_join(cond_format) |>
    dplyr::relocate(dataset, condition_formatted, new_locus_tag, old_locus_tag, baseMean, log2FoldChange,
                    lfcSE, stat, pvalue, padj, condition, control, transcript, Source)

  readr::write_tsv(refseq_processed, "../output/tables/refseq_rnaseq.tsv")
  saveRDS(refseq_processed, "refseq_processed.RDS")

} else {
  message("manually formated condition names not found. Exiting...")
  break()
}
