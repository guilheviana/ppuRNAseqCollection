#!/usr/bin/env Rscript

source("import-metadata.R")

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
