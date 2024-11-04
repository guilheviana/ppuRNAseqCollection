source("utils.R")

data <- list()
datasets <- kt2440 |> dplyr::pull(BioProject) |> unique()

# PRJNA338603 ----
dataset <- datasets[1]
metadata <- prepareMetadata(get(dataset),makeConditionsFrom = c("stress","time_point"), na.value = "0 min")

order <- c("none_0 min", "NaCl_7 min", "NaCl_60 min", "Imipenem_7 min", "Imipenem_60 min", "H2O2_7 min", "H2O2_60 min")

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
order <- c("CK", "LS", "IS", "HS")

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
  dplyr::mutate(Condition = stringr::str_remove_all(Condition, "(in the )|(supplemented with)|( M9)|(\\d{1}.+g/L)"),
                Condition = stringr::str_squish(Condition))
  order <- c("glucose", "glucose_triethanolamine_acetate", "glucose_triethylamine_hydrogen_sulfate", "glucose_acetate")
conditionOrder(metadata) <- order

data <- append(data, list(list(dataset = dataset, metadata = metadata, conditionorder = order)))


# PRJNA796354 ----  <- will come back to this one later
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
order <- c("LB", "LB_selenite")

data <- append(data, list(list(dataset = dataset, metadata = metadata, conditionorder = order)))


all_deseq2 <- NULL
#
# for (ds in unique(data$dataset)) {
#   tempdata <- dplyr::filter(data, dataset == ds)
#   tempmeta <- tempdata$metadata
#   order <- tempdata$order
#
#   conditionOrder(tempmeta) <- order
#
#
#   # obtaining DESeq2 results...
#   txi <- runTxi(tempmeta, dataset = ds)
#
#   tempmeta <- tempmeta |> dplyr::filter(Run %in% colnames(txi$counts))
#
#   # calling DESeq2 on the tximport object
#   ddsTxi <- DESeq2::DESeqDataSetFromTximport(txi, colData = tempmeta, design = ~ Condition) |>
#             DESeq2::DESeq()
#
#   all_comparisons <- data.frame(factorName = "Condition",
#                                 numeratorLevel = levels(tempmeta$Condition)[-1],
#                                 denominatorLevel = levels(tempmeta$Condition)[1])
#
#   res <- getDESeqResults(ddsTxi, ds, comparisons = all_comparisons)
#
#   dplyr::bind_rows(all_deseq2, res)
# }
#
#
