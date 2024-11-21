#!/usr/bin/env Rscript

source("utils.R")

# please activate rsem virtual environment before running this script!

# # checks if the reference rsem files exist, if not create them
# if (!file.exists("software/rsem/refseq.1.bt2")){
#   system("rsem-prepare-reference --bowtie2 software/GCF_000007565.2_ASM756v2_cds_from_genomic.fna software/rsem/refseq")
# }
# if (!file.exists("software/rsem/genbank.1.bt2")){
#   system("rsem-prepare-reference --bowtie2 software/GCA_000007565.2_ASM756v2_cds_from_genomic.fna software/rsem/genbank")
# }
# #

prefixes <- c("genbank", "refseq")

for (prefix in prefixes) {
  for (i in 1:nrow(kt2440)) {
    # gathering run information
    BioProject <- kt2440 |> dplyr::slice(i) |> dplyr::pull(BioProject)
    Run <- kt2440 |> dplyr::slice(i) |> dplyr::pull(Run)

    # creating the output directories
    outDir <- glue::glue("rsem/{prefix}/{BioProject}")
    if (!dir.exists(outDir)) {dir.create(path = outDir, recursive = TRUE)}

    # determines library layout of the run
    files <- list.files(glue::glue("data/{BioProject}/trimmed"), pattern = Run)
    matches <- stringr::str_detect(files, pattern = glue::glue("{Run}_trim.fastq"))
    LibraryLayout <- dplyr::if_else(all(matches), "SINGLE", "PAIRED")

    # preparing the rsem reference
    reference <- glue::glue("software/rsem/{prefix}")

    call <- getRSEMCall(BioProject = BioProject, Run = Run, LibraryLayout = LibraryLayout, outDir = outDir, rsemReference = reference)
    system(call)
  }
}

