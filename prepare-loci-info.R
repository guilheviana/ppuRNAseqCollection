# make the transcript names correspond to the loci names

#genbank
genbank_fa <- read.delim("GCA_000007565.2_ASM756v2_cds_from_genomic.fna.gz", FALSE) |>
  dplyr::filter(stringr::str_starts(V1, pattern = ">")) |> tibble::deframe()

genbank_lookup <- tibble::tibble(
transcript = stringr::str_extract(genbank_fa, "(?<=>).*?(?=\\s)"),
locus_tag = stringr::str_extract(genbank_fa, "(?<=locus_tag=).+?(?=])")
)

#refseq
#
# rsem substitutes the protein id for the locus tag
refseq_fa <- read.delim("GCF_000007565.2_ASM756v2_cds_from_genomic.fna.gz", FALSE) |>
  dplyr::filter(stringr::str_starts(V1, pattern = ">")) |> tibble::deframe()

refseq_lookup <- tibble::tibble(
  og_name = stringr::str_extract(refseq_fa, "(?<=>).*?(?=\\s)"),
  locus_tag = stringr::str_extract(refseq_fa, "(?<=locus_tag=).+?(?=])"),
  protein_id = stringr::str_extract(refseq_fa, "(?<=protein_id=).+?(?=])")
) |>
  dplyr::mutate(new_name = stringr::str_replace(og_name, protein_id, locus_tag))
