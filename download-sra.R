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

             runs <- data |> dplyr::select(BioProject, Organism, Run)

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

for (i in 1:nrow(kt2440)) {

  BioProject <- kt2440 |> dplyr::slice(i) |> dplyr::pull(BioProject)
  Run <- kt2440 |> dplyr::slice(i) |> dplyr::pull(Run)
  possibleOutFiles <- glue::glue("{Run}{suffixes}")

  outdir <- glue::glue("data/{BioProject}")

  if (!dir.exists(outdir)) {dir.create(path = outdir)}

  message <- glue::glue("{BioProject}: retrieving run {Run} [{i}/{nrow(kt2440)}]")
  message(message)

  if (!any(possibleOutFiles %in% list.files(glue::glue("{outdir}/rawdata")))) {
    fasterq_call <- glue::glue("fasterq-dump -p {Run} -e 8 -O {outdir}/rawdata")
    system(fasterq_call)
  } else {
    message(glue::glue("{Run} already present in output directory. Skipping..."))
  }
}
