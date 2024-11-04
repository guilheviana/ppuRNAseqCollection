#!/usr/bin/env Rscript

# please activate rsem virtual environment before running this script!

for (i in 1:nrow(kt2440)) {

  # gathering run information
  BioProject <- kt2440 |> dplyr::slice(i) |> dplyr::pull(BioProject)
  Run <- kt2440 |> dplyr::slice(i) |> dplyr::pull(Run)

  # creating the output directories
  outDir <- glue::glue("rsem/{BioProject}")
  if (!dir.exists(outdir)) {dir.create(path = outDir, recursive = TRUE)}

  # determines library layout of the run
  files <- list.files(glue::glue("data/{BioProject}/trimmed"), pattern = Run)
  matches <- stringr::str_detect(files, pattern = glue::glue("{Run}_trim.fastq"))
  LibraryLayout <- dplyr::if_else(all(matches), "SINGLE", "PAIRED")

  call <- getRSEMCall(BioProject = BioProject, Run = Run, LibraryLayout = LibraryLayout, outDir = outDir)
  message(call)
}

