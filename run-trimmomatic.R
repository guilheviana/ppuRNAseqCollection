#!/usr/bin/env Rscript

# Obs: This script assumes that Trimmomatic (version 0.39) is downloaded and
# located in a directory called "software/trimmomatic"
# Trimmomatic can be downloaded at http://www.usadellab.org/cms/?page=trimmomatic

# In a very rare case (PRJNA450701) there are paired-end and single-end reads in the same
# run, without this difference being specified in the bioproject metadata file
source("utils.R")

for (i in 1:nrow(kt2440)) {

  # gathering run information
  BioProject <- kt2440 |> dplyr::slice(i) |> dplyr::pull(BioProject)
  Run <- kt2440 |> dplyr::slice(i) |> dplyr::pull(Run)

  # creating the output directories
  outdir <- glue::glue("data/{BioProject}/trimmed")
  if (!dir.exists(outdir)) {dir.create(path = outdir)}
  if (!dir.exists(glue::glue("{outdir}/trimlog"))) {dir.create(path = glue::glue("{outdir}/trimlog"))}

  # determines library layout of the run
  files <- list.files(glue::glue("data/{BioProject}/rawdata/"), pattern = Run)
  matches <- stringr::str_detect(files, pattern = glue::glue("{Run}.fastq"))
  LibraryLayout <- dplyr::if_else(all(matches), "SINGLE", "PAIRED")

  # trimmomatic function call
  call <- getTrimmomaticCall(BioProject = BioProject, Run = Run, LibraryLayout = LibraryLayout, outdir = outdir)
  system(call)
}
