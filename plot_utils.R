scatterPlot <- function(cond.choose) {

  gb_rnaseq <- genbank |> dplyr::filter(condition == cond.choose) |>
    dplyr::select(log2FoldChange, locus_tag, condition_formatted)

  rs_rnaseq <- refseq |> dplyr::filter(condition == cond.choose) |>
    dplyr::select(log2FoldChange, old_locus_tag, condition_formatted)

  joined_rnaseq <- dplyr::inner_join(gb_rnaseq, rs_rnaseq,
                                     by = c("locus_tag" = "old_locus_tag", "condition_formatted"),
                                     suffix = c(".genbank", ".refseq")) |>
    dplyr::mutate(
      affected = dplyr::if_else(locus_tag %in% affected_loci.old_codes, TRUE, FALSE),
      delta_lfc = abs(log2FoldChange.genbank - log2FoldChange.refseq))

  # retrieving the outliers that are unshifted
  unshifted_large_variations <- joined_rnaseq |>
    dplyr::filter(delta_lfc >= 2.5)

  caption <- unique(joined_rnaseq$condition_formatted)

  ggplot(joined_rnaseq, aes(x = log2FoldChange.genbank, y = log2FoldChange.refseq, color = affected)) +
    geom_vline(xintercept = 0, color = "grey70") +
    geom_hline(yintercept = 0, color = "grey70") +
    geom_point(size = 2, alpha = 1) +
   ggrepel::geom_text_repel(data = unshifted_large_variations, aes(label = locus_tag), color = "black", force_pull = 0.3, size = 4.5) +
    scale_color_manual(name = "",
                       labels = c("Same position in both annotations", "Shifted genes"),
                       values = c( "TRUE" = "steelblue", "FALSE" = "grey"))+
    guides(color = guide_legend(override.aes = list(size=5))) +
    labs(subtitle = caption, y = "RefSeq Log2 Fold Change", x = "GenBank Log2 Fold Change") +
    theme_light(base_size = 16) +
    theme(
      panel.border = element_rect(color = NA),
      panel.grid.minor = element_blank(),

      plot.subtitle = ggtext::element_markdown(hjust = 0.5),

      legend.position = "inside",
      legend.background = element_blank(),
      legend.position.inside = c(0.3, 0.9),
      legend.title = element_blank(),
      legend.text = element_text(color = "black", size = 14),

      axis.ticks = element_blank(),
      axis.text = element_text(color = "black", size = 16)
    )
}

vennPlot <- function(DEGs_table) {

  condition <- DEGs_table |> dplyr::pull(condition) |> unique()

  gb_exclusive <- DEGs_table |> dplyr::filter(parameter == "genbank_exc") |> dplyr::pull(loci)
  rs_exclusive <- DEGs_table |> dplyr::filter(parameter == "refseq_exc") |> dplyr::pull(loci)
  intersection <- DEGs_table |> dplyr::filter(parameter == "common") |> dplyr::pull(loci)

  if(any(is.na(gb_exclusive))) {gb_exclusive <- "No exclusive DEGs"}
  if(any(is.na(rs_exclusive))) {rs_exclusive <- "No exclusive DEGs"}

  if(length(intersection) <= 10){
    intersection_label <- stringr::str_flatten(intersection, collapse = "\n")
  } else {
    intersection_label <- glue::glue("{length(intersection)} shared DEGs")
  }

  vennInfo <- data.frame(x = c(0.65, 1.5, 2.35),
                         y = rep(0, 3),
                         label = c(stringr::str_flatten(gb_exclusive, collapse = "\n"),
                                   intersection_label,
                                   stringr::str_flatten(rs_exclusive, collapse = "\n")),
                         title = c("GenBank", NA, "RefSeq"),
                         x_title = c(0.4, NA, 2.5))

  circles <- data.frame(x0 = c(1,2), y0 = c(0,0), r = c(1,1), groups = c("1", "2"))

  custom_colors <- c("#8cba9d", "steelblue")

  vennDiag <- ggplot(vennInfo) +
    ggforce::geom_circle(data = circles,
                         aes(x0 = x0, y0 = y0, r = r, fill = groups),
                         radius = 3, color = NA, alpha = 0.7, show.legend = FALSE) +
    geom_text(data = vennInfo,
              aes(x = x, y = y, label = label), size = 4) +
    geom_text(data = vennInfo,
              aes(x = x_title, y = y + 1, label = title), size = 5) +
    scale_fill_manual(values = custom_colors) +
    labs(title = condition) +
    coord_fixed(clip = "off") +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = 20, hjust = 0.5, face = "bold", margin = margin(0, 0, 0.5, 0, "cm"))
    )

  return(vennDiag)
}

findDEGs <- function(condition,
                     dataset.refseq = refseq_categorized,
                     dataset.genbank = genbank_categorized) {

  filter <- condition

  filter.format <- dplyr::bind_rows(dataset.genbank, dataset.refseq) |>
                   dplyr::filter(condition == filter) |>
                   dplyr::pull(condition_formatted) |>
                   unique()

  g <- dataset.genbank |>
    dplyr::filter(significant, condition == filter) |>
    dplyr::pull(locus_tag) |> sort() |> unique()
  r <- dataset.refseq |>
    dplyr::filter(significant, condition == filter) |>
    dplyr::pull(old_locus_tag) |> sort() |> unique()

  intersection <- intersect(g, r)
  genbank_exc <- setdiff(g, r)
  refseq_exc <- setdiff(r,g)

  if(identical(intersection, character(0))) {
    intersection <- NA
  }
  if(identical(genbank_exc, character(0))) {
    genbank_exc <- NA
  }
  if(identical(refseq_exc, character(0))) {
    refseq_exc <- NA
  }


  output <- tibble::tibble(parameter = "common", loci = intersection, condition = filter.format) |>
    dplyr::bind_rows(
      tibble::tibble(parameter = "genbank_exc", loci = genbank_exc, condition = filter.format)
    ) |>
    dplyr::bind_rows(
      tibble::tibble(parameter = "refseq_exc", loci = refseq_exc, condition = filter.format)
    )
  return(output)
}
