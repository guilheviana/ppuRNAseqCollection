source("./plot_utils.R")
library(ggplot2)

genbank <- readRDS("genbank_processed.RDS")
refseq <- readRDS("refseq_processed.RDS")
load("affected_oldformat") # loads the vector with the genome type annotation

conditions <- c("coumarate", "glucose_triethylamine_hydrogen_sulfate", "s_1h", "myristic_acid")

  for (i in conditions) {
  plot <- scatterPlot(cond.choose = i)
  ggsave(glue::glue("figure 3 - {i}.svg"), plot = plot, device = svg, width = 6, height = 6, path = "../output/panels/", dpi = 300)
}
