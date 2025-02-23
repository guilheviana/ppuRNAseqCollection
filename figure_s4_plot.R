source("./plot_utils.R")
library(ggplot2)

genbank <- readRDS("genbank_processed.RDS")
refseq <- readRDS("refseq_processed.RDS")
load("affected_oldformat") # loads the vector with the genome type annotation

# figure with the statistically significant genes

genbank_categorized <- genbank |>
  dplyr::mutate(shifted_gene = dplyr::if_else(locus_tag %in% affected_loci.old_codes, TRUE, FALSE),
                significant = dplyr::if_else(padj < 0.05 & abs(log2FoldChange) >= 2, TRUE, FALSE))


refseq_categorized <- refseq |>
  dplyr::mutate(shifted_gene = dplyr::if_else(old_locus_tag %in% affected_loci.old_codes, TRUE, FALSE),
                significant = dplyr::if_else(padj < 0.05 & abs(log2FoldChange) >= 2, TRUE, FALSE))


conditions <- c("coumarate", "glucose_triethylamine_hydrogen_sulfate", "s_1h", "myristic_acid")


for(cond in conditions) {

  a <- findDEGs(condition = "coumarate")
  vennDiag <- vennPlot(a)

  View(a)

  # ggsave(vennDiag,
  #        filename = glue::glue("Fig S4. - VennDiagram - {cond}.svg"),
  #        device = svg, path = "../output/panels/",
  #        width = 6,
  #        height = 6,
  #        units = "in",
  #        dpi = 300)
}


