#!/usr/bin/env Rscript
#
# ESR1 locus: chromatin accessibility across the six groups
#
# Run from the analysis directory:
#
#   Rscript atacSeq/custom/ESR1/plotESR1Locus.R
#
# Reads the 10mNorm bigWigs, the peak atlas counts and the edgeR results
# of atacSeq/diff. Writes to atacSeq/custom/ESR1:
#
#   Proj_18143__ESR1_Tracks.pdf  locus overview (group means, then every
#                                sample), then one zoom per ESR1 TSS
#                                cluster
#   Proj_18143__ESR1_Peaks.pdf   TMM log2 CPM by group of the ESR1
#                                promoter peaks and the differential peaks
#   Proj_18143__ESR1_Peaks.xlsx  every atlas peak in the region
#

cat("## START: plotESR1Locus.R\n")

suppress({
  library(GenomicRanges)
  library(rtracklayer)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(edgeR)
  library(patchwork)
  library(openxlsx)
  library(fs)
  library(tidyverse)
})

PROJ <- "Proj_18143"
GENE <- "ESR1"
FLANK <- 50e3      # bp shown on each side of the gene
TSS_GAP <- 5e3     # TSSs closer than this share one zoom window
ZOOM_PAD <- 3e3    # bp shown on each side of a TSS cluster
PROMOTER <- 1e3    # a peak within this of a TSS is a promoter peak
FDR_CUT <- 0.05
N_BINS <- 1500     # most coverage bins across a panel
MIN_BIN <- 10      # smallest bin in bp
PEAKS_PER_PAGE <- 6
PEAKS_PER_ROW <- 3
GROUP_ORDER <- c("EV", "ED1", "ED6", "BRAF", "BD1", "BD6")

OUT_DIR <- "atacSeq/custom/ESR1"
ATAC_SDIR <- "ATAC-seq"

source(path(ATAC_SDIR, "R", "palettes.R"))
GROUP_COLORS <- set_names(PAL_CAT[seq_along(GROUP_ORDER)], GROUP_ORDER)
FEATURE_COLOR <- "black"   # gene models and atlas peaks

################################################################################
# Samples
################################################################################

manifest <- read_csv("sampleManifest.csv", show_col_types = FALSE) |>
  mutate(Group = fct(Group, levels = GROUP_ORDER),
         BigWig = path("atacSeq/bigwig",
                       str_glue("{MapID}_postProcess.shifted.10mNorm.bw"))) |>
  arrange(Group, SampleID)

if (!all(file_exists(manifest$BigWig))) {
  rlang::abort(str_c("missing bigWig: ",
                     str_c(manifest$BigWig[!file_exists(manifest$BigWig)],
                           collapse = " ")))
}

################################################################################
# Gene models
################################################################################

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

gene_id <- AnnotationDbi::mapIds(org.Hs.eg.db, GENE, "ENTREZID", "SYMBOL")
gene_tx <- GenomicFeatures::transcriptsBy(txdb, "gene")[[gene_id]]
gene_range <- range(gene_tx)
# The peak plot keeps to GENE: its span plus the promoter window at the
# 5' end
gene_span <- resize(gene_range, width(gene_range) + PROMOTER, fix = "end")
strand(gene_range) <- "*"
region <- gene_range + FLANK
chrom <- as.character(seqnames(region))

tx <- GenomicFeatures::transcripts(txdb,
                                   columns = c("tx_id", "tx_name", "gene_id")) |>
  subsetByOverlaps(region)

tx_tbl <- tibble(tx_id = as.character(tx$tx_id),
                 tx_name = tx$tx_name,
                 gene_id = map_chr(as.list(tx$gene_id), \(x) x[1]),
                 Start = start(tx),
                 End = end(tx),
                 Strand = as.character(strand(tx))) |>
  filter(!is.na(gene_id))

symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, unique(tx_tbl$gene_id),
                                 "SYMBOL", "ENTREZID")
tx_tbl <- tx_tbl |>
  mutate(Symbol = coalesce(unname(symbols[gene_id]), gene_id),
         TSS = if_else(Strand == "+", Start, End))

exon_tbl <- GenomicFeatures::exonsBy(txdb, "tx")[tx_tbl$tx_id] |>
  as.data.frame() |>
  as_tibble() |>
  dplyr::select(tx_id = group_name, Start = start, End = end)

# Gene model track: GENE only, one row per distinct exon structure, named
# by its first transcript
gene_rows <- exon_tbl |>
  semi_join(tx_tbl |> filter(Symbol == GENE), by = "tx_id") |>
  arrange(tx_id, Start) |>
  summarize(Key = str_c(Start, "-", End, collapse = ","), .by = tx_id) |>
  left_join(tx_tbl, by = "tx_id") |>
  arrange(Start, End, tx_name) |>
  distinct(Key, .keep_all = TRUE) |>
  mutate(Row = str_c(Symbol, " ", tx_name))

gene_exons <- exon_tbl |>
  inner_join(gene_rows |> dplyr::select(tx_id, Row), by = "tx_id")

model_rows <- gene_rows |>
  dplyr::select(Row, Start, End, Strand)
model_exons <- gene_exons |>
  dplyr::select(Row, Start, End)

# Zoom windows: the GENE TSSs, clustered
gene_tss <- gene_tx |> resize(1, fix = "start")
strand(gene_tss) <- "*"
tss_clusters <- GenomicRanges::reduce(gene_tss, min.gapwidth = TSS_GAP)
tss_clusters$nTx <- countOverlaps(tss_clusters, gene_tss)
zooms <- tss_clusters + ZOOM_PAD

################################################################################
# Atlas peaks, normalized counts and edgeR results in the region
################################################################################

raw_counts <- read_tsv("atacSeq/atlas/peaks_raw_fcCounts.txt", comment = "#",
                       show_col_types = FALSE)

count_matrix <- raw_counts |>
  dplyr::select(PeakNo = Geneid, ends_with(".bam")) |>
  as.data.frame() |>
  column_to_rownames("PeakNo") |>
  as.matrix()
colnames(count_matrix) <- colnames(count_matrix) |>
  path_file() |>
  str_remove("_postProcess.bam$")
count_matrix <- count_matrix[, manifest$MapID]
colnames(count_matrix) <- manifest$SampleID

# TMM over all samples, as in diffAnalysisPairwise.R before its
# per-comparison step
dge <- DGEList(counts = count_matrix, group = manifest$Group)
dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
eff_lib_size <- dge$samples$lib.size * dge$samples$norm.factors

peaks <- raw_counts |>
  dplyr::select(PeakNo = Geneid, Chr, Start, End, Length) |>
  filter(Chr == chrom, End >= start(region), Start <= end(region)) |>
  arrange(Start)

log_cpm <- cpm(count_matrix[peaks$PeakNo, , drop = FALSE],
               lib.size = eff_lib_size, log = TRUE)

diff_res <- readRDS(path("atacSeq/diff", str_glue("{PROJ}__DiffPeaksEdgeR_V5.rds")))
comparisons <- diff_res$comparisons |>
  map_dfr(\(cmp) tibble(Tag = cmp$tag, Label = cmp$label))

diff_tbl <- diff_res$comparisons |>
  map_dfr(\(cmp) {
    cmp$table |>
      filter(PeakNo %in% peaks$PeakNo) |>
      mutate(Tag = cmp$tag, Label = cmp$label)
  }) |>
  mutate(Label = fct(Label, levels = comparisons$Label))

diff_sig <- diff_tbl |>
  filter(FDR < FDR_CUT) |>
  left_join(peaks, by = "PeakNo") |>
  mutate(Direction = fct(if_else(logFC > 0, "up", "down"),
                         levels = c("down", "up")))

################################################################################
# Peak locations
################################################################################

peaks_gr <- GRanges(peaks$Chr, IRanges::IRanges(peaks$Start, peaks$End))

#' Name the genes whose features overlap each peak
#'
#' @param features GRanges with a Symbol column
#' @param what Character feature name appended to each symbol
#' @return Tibble of peak index i and Label
overlap_labels <- function(features, what) {
  hits <- findOverlaps(peaks_gr, features, ignore.strand = TRUE)
  tibble(i = queryHits(hits), Symbol = features$Symbol[subjectHits(hits)]) |>
    distinct() |>
    summarize(Label = str_c(sort(Symbol), " ", what, collapse = ", "),
              .by = i)
}

tss_window_gr <- GRanges(chrom,
                         IRanges::IRanges(tx_tbl$TSS - PROMOTER,
                                          tx_tbl$TSS + PROMOTER),
                         Symbol = tx_tbl$Symbol)
exon_gr <- exon_tbl |>
  left_join(tx_tbl |> dplyr::select(tx_id, Symbol), by = "tx_id") |>
  with(GRanges(chrom, IRanges::IRanges(Start, End), Symbol = Symbol))
body_gr <- GRanges(chrom, IRanges::IRanges(tx_tbl$Start, tx_tbl$End),
                   Symbol = tx_tbl$Symbol)

gene_tss_pos <- start(gene_tss)
dist_col <- str_glue("DistTo{GENE}TSS")

peaks <- peaks |>
  mutate(i = row_number(),
         Center = round((Start + End) / 2)) |>
  left_join(overlap_labels(tss_window_gr, "promoter") |>
              dplyr::rename(Promoter = Label), by = "i") |>
  left_join(overlap_labels(exon_gr, "exon") |>
              dplyr::rename(Exon = Label), by = "i") |>
  left_join(overlap_labels(body_gr, "intron") |>
              dplyr::rename(Intron = Label), by = "i") |>
  mutate(Location = coalesce(Promoter, Exon, Intron, "intergenic"),
         "{dist_col}" := map_dbl(Center, \(x) {
           d <- x - gene_tss_pos
           d[which.min(abs(d))]
         })) |>
  dplyr::select(-i, -Promoter, -Exon, -Intron)

################################################################################
# Coverage
################################################################################

#' Mean 10mNorm coverage of every sample in bins across a window
#'
#' @param window GRanges of one range
#' @return Tibble of SampleID, Group, Start, End, Score
sample_coverage <- function(window) {
  n_bins <- min(N_BINS, floor(width(window) / MIN_BIN))
  manifest |>
    dplyr::select(SampleID, Group, BigWig) |>
    pmap_dfr(\(SampleID, Group, BigWig) {
      bins <- summary(BigWigFile(as.character(BigWig)), which = window, size = n_bins,
                      type = "mean", defaultValue = 0)[[1]]
      tibble(SampleID = SampleID, Group = Group,
             Start = start(bins), End = end(bins), Score = bins$score)
    })
}

#' Mean over the replicates of each group
#'
#' @param cov Tibble from sample_coverage()
#' @return Tibble of Group, Start, End, Score
group_coverage <- function(cov) {
  cov |>
    summarize(Score = mean(Score), .by = c(Group, Start, End))
}

################################################################################
# Plots
################################################################################

theme_track <- theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(1, 4, 1, 4))

#' x scale and limits shared by every track of a window
#'
#' @param window GRanges of one range
#' @return List of ggplot2 components
x_window <- function(window) {
  labels <- if (width(window) > 100e3) {
    scales::label_number(scale = 1e-6, accuracy = 0.01, suffix = " Mb")
  } else {
    scales::label_comma()
  }
  list(scale_x_continuous(labels = labels),
       coord_cartesian(xlim = c(start(window), end(window)), expand = FALSE))
}

#' Coverage track, one facet per row of Track
#'
#' @param cov Tibble of Track, Group, Start, End, Score
#' @param window GRanges of one range
#' @return ggplot
plot_coverage <- function(cov, window) {
  tss <- gene_tss_pos[gene_tss_pos >= start(window) &
                        gene_tss_pos <= end(window)]
  ggplot(cov) +
    geom_vline(xintercept = tss, colour = REF_LINE, linewidth = 0.2,
               linetype = "dotted") +
    geom_rect(aes(xmin = Start, xmax = End + 1, ymin = 0, ymax = Score,
                  fill = Group)) +
    facet_grid(Track ~ ., switch = "y") +
    scale_fill_manual(values = GROUP_COLORS, guide = "none") +
    scale_y_continuous(n.breaks = 3, expand = expansion(mult = c(0, 0.05))) +
    x_window(window) +
    theme_track +
    theme(strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9),
          panel.spacing.y = unit(1, "pt"),
          panel.border = element_rect(colour = TILE_BORDER, fill = NA))
}

#' Rectangles of features, widened to be visible at the window's scale
#'
#' @param feats Tibble with Start, End
#' @param window GRanges of one range
#' @return feats with XMin and XMax
widen <- function(feats, window) {
  min_width <- width(window) / 400
  feats |>
    mutate(Mid = (Start + End) / 2,
           Half = pmax(End - Start + 1, min_width) / 2,
           XMin = Mid - Half,
           XMax = Mid + Half)
}

#' Atlas peak track
#'
#' @param window GRanges of one range
#' @return ggplot
plot_peaks <- function(window) {
  ggplot(widen(peaks, window)) +
    geom_rect(aes(xmin = XMin, xmax = XMax, ymin = 0.2, ymax = 0.8),
              fill = FEATURE_COLOR) +
    scale_y_continuous(breaks = 0.5, labels = "Atlas peaks",
                       limits = c(0, 1)) +
    x_window(window) +
    theme_track
}

#' Differential peaks track: one row per comparison
#'
#' @param window GRanges of one range
#' @return ggplot
plot_diff <- function(window) {
  labels <- levels(diff_tbl$Label)
  sig <- widen(diff_sig, window) |>
    mutate(Y = length(labels) + 1 - as.integer(Label))
  ggplot(sig) +
    geom_hline(yintercept = seq_along(labels), colour = NA_GREY,
               linewidth = 0.2) +
    geom_rect(aes(xmin = XMin, xmax = XMax, ymin = Y - 0.35, ymax = Y + 0.35,
                  fill = Direction)) +
    scale_fill_manual(values = SIDE_COLORS, drop = FALSE,
                      name = str_glue("FDR < {FDR_CUT}"),
                      labels = c(down = "lower in 1st group",
                                 up = "higher in 1st group")) +
    scale_y_continuous(breaks = seq_along(labels), labels = rev(labels),
                       limits = c(0.5, length(labels) + 0.5)) +
    x_window(window) +
    theme_track +
    theme(legend.position = "right")
}

#' Gene model track, the last of a page, with the x axis
#'
#' @param window GRanges of one range
#' @return ggplot
plot_genes <- function(window) {
  rows <- model_rows |>
    filter(End >= start(window), Start <= end(window)) |>
    mutate(Y = n() + 1 - row_number())
  exons <- model_exons |>
    inner_join(rows |> dplyr::select(Row, Y), by = "Row")
  step <- width(window) / 40
  arrows <- rows |>
    reframe(X = seq(max(Start, start(window)), min(End, end(window)),
                    by = step),
            Dir = if_else(Strand == "+", 1, -1),
            .by = c(Row, Y))

  ggplot(rows) +
    geom_segment(aes(x = Start, xend = End, y = Y, yend = Y),
                 colour = FEATURE_COLOR, linewidth = 0.3) +
    geom_segment(data = arrows,
                 aes(x = X, xend = X + Dir * step * 0.2, y = Y, yend = Y),
                 colour = FEATURE_COLOR, linewidth = 0.3,
                 arrow = arrow(length = unit(1.2, "mm"))) +
    geom_rect(data = exons,
              aes(xmin = Start, xmax = End + 1, ymin = Y - 0.3, ymax = Y + 0.3),
              fill = FEATURE_COLOR) +
    scale_y_continuous(breaks = rows$Y, labels = rows$Row,
                       limits = c(0.5, nrow(rows) + 0.5)) +
    x_window(window) +
    labs(x = str_glue("{chrom} position (hg38)")) +
    theme_track +
    theme(axis.title.x = element_text(),
          axis.text.x = element_text(),
          axis.ticks.x = element_line())
}

#' One page: coverage, atlas peaks, differential peaks and gene models
#'
#' @param cov Tibble of Track, Group, Start, End, Score
#' @param window GRanges of one range
#' @param title Character page title
#' @param subtitle Character page subtitle
#' @return patchwork
plot_locus <- function(cov, window, title, subtitle) {
  n_tracks <- n_distinct(cov$Track)
  n_genes <- model_rows |>
    filter(End >= start(window), Start <= end(window)) |>
    nrow()
  wrap_plots(plot_coverage(cov, window), plot_peaks(window),
             plot_diff(window), plot_genes(window), ncol = 1) +
    plot_layout(heights = c(n_tracks, 0.4, 0.3 * nrow(comparisons),
                            0.3 * n_genes + 0.6),
                guides = "collect") +
    plot_annotation(title = title, subtitle = subtitle,
                    theme = theme(plot.title = element_text(size = 12),
                                  plot.subtitle = element_text(size = 9)))
}

#' Format a window as chr:start-end
#'
#' @param window GRanges of one range
#' @return Character
locus_str <- function(window) {
  str_glue("{seqnames(window)}:{scales::comma(start(window))}-",
           "{scales::comma(end(window))}")
}

cov_note <- str_c("Coverage: 10mNorm bigWig (reads per 10M), mean per bin. ",
                  "Dotted lines: ", GENE, " TSSs. Differential peaks: edgeR, ",
                  "FDR < ", FDR_CUT, "; 1st group is the first named.")

region_cov <- sample_coverage(region)

pages <- list(
  plot_locus(group_coverage(region_cov) |> mutate(Track = Group),
             region,
             str_glue("{GENE} +/- {FLANK / 1e3} kb, {locus_str(region)}: ",
                      "group means"),
             cov_note),
  plot_locus(region_cov |>
               mutate(Track = fct(SampleID, levels = manifest$SampleID)),
             region,
             str_glue("{GENE} +/- {FLANK / 1e3} kb, {locus_str(region)}: ",
                      "every sample"),
             cov_note)
)

for (k in seq_along(zooms)) {
  window <- zooms[k]
  pages <- c(pages, list(plot_locus(
    group_coverage(sample_coverage(window)) |> mutate(Track = Group),
    window,
    str_glue("{GENE} TSS cluster {k} of {length(zooms)} ",
             "({tss_clusters$nTx[k]} transcripts), {locus_str(window)}: ",
             "group means"),
    cov_note)))
}

dir_create(OUT_DIR)
tracks_pdf <- path(OUT_DIR, str_glue("{PROJ}__{GENE}_Tracks.pdf"))
pdf(tracks_pdf, width = 11, height = 8.5)
walk(pages, print)
dev.off()

################################################################################
# Peak accessibility by group
################################################################################

sample_log_cpm <- log_cpm |>
  as.data.frame() |>
  rownames_to_column("PeakNo") |>
  pivot_longer(-PeakNo, names_to = "SampleID", values_to = "log2CPM") |>
  left_join(manifest |> dplyr::select(SampleID, Group), by = "SampleID")

# Peaks in gene_span only; the flanks are in the tracks and the xlsx.
# GENE promoter peaks first, then the rest; each by its smallest FDR
# over the comparisons
min_fdr <- diff_tbl |>
  slice_min(FDR, n = 1, with_ties = FALSE, by = PeakNo) |>
  dplyr::select(PeakNo, MinFDR = FDR, MinFDRLabel = Label)

dot_peaks <- peaks |>
  mutate(IsPromoter = str_detect(Location, str_glue("\\b{GENE} promoter"))) |>
  filter(End >= start(gene_span), Start <= end(gene_span),
         IsPromoter | PeakNo %in% diff_sig$PeakNo) |>
  left_join(diff_sig |> count(PeakNo, name = "nSig"), by = "PeakNo") |>
  left_join(min_fdr, by = "PeakNo") |>
  arrange(desc(IsPromoter), MinFDR) |>
  mutate(nSig = replace_na(nSig, 0),
         FDRText = if_else(is.na(MinFDR), "not tested",
                           str_glue("min FDR {signif(MinFDR, 2)} ",
                                    "({MinFDRLabel})")),
         Facet = str_glue("{PeakNo}  {Chr}:{scales::comma(Start)}-",
                          "{scales::comma(End)}\n{Location}\n",
                          "{FDRText}\nFDR < {FDR_CUT} in {nSig} of ",
                          "{nrow(comparisons)}"),
         Facet = fct_inorder(Facet))

dot_data <- sample_log_cpm |>
  inner_join(dot_peaks |> dplyr::select(PeakNo, Facet), by = "PeakNo")

dot_pages <- dot_peaks |>
  mutate(Page = (row_number() - 1) %/% PEAKS_PER_PAGE) %>%
  split(.$Page)

peaks_pdf <- path(OUT_DIR, str_glue("{PROJ}__{GENE}_Peaks.pdf"))
pdf(peaks_pdf, width = 11, height = 8.5)
walk(dot_pages, \(page) {
  p <- dot_data |>
    filter(PeakNo %in% page$PeakNo) |>
    mutate(Facet = fct_drop(Facet)) |>
    ggplot(aes(Group, log2CPM, colour = Group)) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.6,
                 linewidth = 0.3, colour = REF_LINE) +
    geom_point(size = 2.5) +
    facet_wrap(~Facet, ncol = PEAKS_PER_ROW, scales = "free_y") +
    scale_colour_manual(values = GROUP_COLORS, guide = "none") +
    labs(x = NULL, y = "log2 CPM (TMM, all samples)",
         title = str_glue("{GENE} peaks ({locus_str(gene_span)}): promoter ",
                          "peaks and peaks with FDR < {FDR_CUT} in any ",
                          "comparison"),
         subtitle = str_glue("Points: samples. Bars: group mean. Order: ",
                             "{GENE} promoter peaks first, then by smallest ",
                             "FDR in any comparison.")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = REF_LINE, fill = NA,
                                      linewidth = 0.2),
          strip.text = element_text(size = 10, hjust = 0))
  print(p)
})
dev.off()

################################################################################
# Peak table
################################################################################

mean_cols <- sample_log_cpm |>
  summarize(log2CPM = mean(log2CPM), .by = c(PeakNo, Group)) |>
  arrange(Group) |>
  pivot_wider(names_from = Group, values_from = log2CPM,
              names_prefix = "MeanLog2CPM_")

stat_cols <- diff_tbl |>
  dplyr::select(PeakNo, Tag, logFC, FDR) |>
  mutate(Tag = fct(Tag, levels = comparisons$Tag)) |>
  arrange(Tag) |>
  pivot_wider(names_from = Tag, values_from = c(logFC, FDR),
              names_vary = "slowest")

sample_cols <- log_cpm |>
  as.data.frame() |>
  rename_with(\(x) str_c("Log2CPM_", x)) |>
  rownames_to_column("PeakNo")

peak_table <- peaks |>
  dplyr::select(PeakNo, Chr, Start, End, Length, Location,
                all_of(dist_col)) |>
  left_join(mean_cols, by = "PeakNo") |>
  left_join(stat_cols, by = "PeakNo") |>
  left_join(sample_cols, by = "PeakNo") |>
  mutate(across(starts_with(c("MeanLog2CPM_", "Log2CPM_", "logFC_")),
                \(x) round(x, 3)),
         across(starts_with("FDR_"), \(x) signif(x, 3)))

info <- tibble(
  Item = c("Region", "Gene", "Peaks", "Start, End", "Location",
           dist_col, "MeanLog2CPM_<group>", "logFC_<comparison>",
           "FDR_<comparison>", "Log2CPM_<sample>"),
  Description = c(
    str_glue("{locus_str(region)} (hg38): {GENE} +/- {FLANK / 1e3} kb"),
    str_glue("{GENE}, Entrez {gene_id}, TxDb.Hsapiens.UCSC.hg38.knownGene"),
    "Every peak of the merged atlas (atacSeq/atlas) in the region",
    "1-based, inclusive (SAF)",
    str_glue("Gene promoter (a TSS within {PROMOTER} bp), else exon, else ",
             "intron, else intergenic; all genes in the region"),
    str_glue("Peak center minus the nearest {GENE} TSS, bp"),
    "Mean over the group's samples of Log2CPM",
    str_c("edgeR logFC of the comparison (atacSeq/diff); positive is higher ",
          "in the first group named. Empty: the peak was not tested ",
          "(filterByExpr)"),
    "edgeR FDR of the comparison (BH)",
    str_c("edgeR log2 CPM (prior count 2) with TMM factors computed over ",
          "all samples; the comparisons use per-comparison TMM factors")
  )
)

xlsx_file <- path(OUT_DIR, str_glue("{PROJ}__{GENE}_Peaks.xlsx"))
write.xlsx(list(Peaks = peak_table, Info = info), xlsx_file,
           firstRow = TRUE, colWidths = "auto", overwrite = TRUE)

cat("\nRegion:", locus_str(region), "\n")
cat("Atlas peaks in region:", nrow(peaks), "\n")
cat(str_glue("Peaks in the peak plot ({locus_str(gene_span)}): ",
             "{nrow(dot_peaks)}"), "\n")
cat(str_glue("{GENE} TSS clusters (zoom pages): {length(zooms)}"), "\n")
cat(str_glue("Peaks with FDR < {FDR_CUT} per comparison:"), "\n")
diff_sig |>
  count(Label, Direction, .drop = FALSE) |>
  pivot_wider(names_from = Direction, values_from = n, values_fill = 0) |>
  print()
cat("\nWrote:", tracks_pdf, peaks_pdf, xlsx_file, sep = "\n  ")

cat("\n## END: plotESR1Locus.R\n")
