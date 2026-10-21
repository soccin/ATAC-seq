suppressPackageStartupMessages(library(optparse))

usage <- "
  %prog [options] [RUNTAG]

Pathway enrichment (preranked GSEA, fgsea) of the differential peaks
of R/diffAnalysisPairwise.R against MSigDB gene sets (msigdbr).

Run from the analysis directory, after diffAnalysisPairwise.R with the
same RUNTAG (default: none). Its input is that script's
<PROJ>_<RUNTAG>_DiffPeaksEdgeR_V5.rds, which holds every tested peak
of every comparison. The genome is GENOME= in 00.RUNSTATUS.txt; the
annotation, the MSigDB species and the default collections come from
lib/genomes/GENOMES.tsv.

For each comparison every tested peak gets the nearest gene
(ChIPseeker), peaks too far from that gene's TSS (--mode) are dropped,
and each gene gets the statistic sign(logFC) * -log10(PValue) of its
peak with the largest absolute value. The genes ranked by it are
tested with fgseaMultilevel against each collection; padj is
Benjamini-Hochberg within the collection. A positive NES means the set
is more accessible in Group2 of the comparison.

COLLECTIONS is a comma separated list of COLL or COLL:SUBCOLL
(msigdbr's gs_collection and gs_subcollection, e.g. H, C2:CP:REACTOME,
C5:GO:BP; mouse: MH, M2:CP:REACTOME, M5:GO:BP). A subcollection also
matches the ones below it: C2:CP is CP and every CP:*. An item
starting with - removes those sets from every other item, e.g.
-C2:CP:KEGG_LEGACY. Each item without - is one fgsea run.

OUTPUT FILES (in the current directory):
  <PROJ>_<RUNTAG>_PathwayGSEA.xlsx  Summary, one sheet per comparison,
                                    Parameters
  <PROJ>_<RUNTAG>_PathwayGSEA.pdf   NES of the top sets, per comparison
                                    and collection
  <PROJ>_<RUNTAG>_PathwayGSEA.rds   gene ranks and full results
  With --out-tag T the names end _T_PathwayGSEA.*"

option_list <- list(
  make_option(c("-c", "--collections"), type = "character", default = NULL,
              help = "MSigDB collections [default: MSIGDB_SETS of GENOMES.tsv]"),
  make_option(c("-m", "--mode"), type = "character", default = "promoter",
              help = "promoter or distance [default: %default]"),
  make_option(c("-w", "--tss-window"), type = "integer", default = 1000L,
              help = "promoter mode: keep peaks with |distanceToTSS| <= this [default: %default]"),
  make_option(c("-d", "--max-distance"), type = "integer", default = 50000L,
              help = "distance mode: keep peaks with |distanceToTSS| <= this [default: %default]"),
  make_option("--min-size", type = "integer", default = 15L,
              help = "smallest gene set tested, in ranked genes [default: %default]"),
  make_option("--max-size", type = "integer", default = 500L,
              help = "largest gene set tested, in ranked genes [default: %default]"),
  make_option(c("-o", "--out-tag"), type = "character", default = "",
              help = "tag added to the output names, e.g. for a second --mode [default: none]")
)
parsed <- parse_args(OptionParser(usage = usage, option_list = option_list),
                     positional_arguments = c(0, 1),
                     convert_hyphens_to_underscores = TRUE)
opts <- parsed$options
RUNTAG <- if (len(parsed$args) == 1) parsed$args[1] else ""

if (!opts$mode %in% c("promoter", "distance")) {
  stop("--mode must be promoter or distance, not ", opts$mode, call. = FALSE)
}

cat("## START: pathwayGSEA.R\n")

cat("Loading libraries ...")
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicRanges)
  library(fgsea)
  library(msigdbr)
  library(tidyverse)
  library(fs)
  library(openxlsx)
})
cat(" done\n")

# The ATAC-seq checkout is the parent of this script's directory
SDIR <- commandArgs(trailingOnly = FALSE) |>
  str_subset("^--file=") |>
  str_remove("^--file=") |>
  normalizePath() |>
  dirname() |>
  dirname()
source(file.path(SDIR, "R", "genomeInfo.R"))

GSEA_SEED <- 42
FDR_CUT <- 0.05
N_TOP <- 15

# Genome of the run, its annotation and its MSigDB, from GENOMES.tsv
GENOME <- runGenome()
cat("GENOME =", GENOME, "\n")
txdb_pkg <- genomeInfo(GENOME, "TXDB", SDIR)
txdb <- getExportedValue(txdb_pkg, txdb_pkg)
annoDb <- genomeInfo(GENOME, "ORGDB", SDIR)
add_chr <- genomeInfo(GENOME, "ADD_CHR", SDIR) == "yes"
mito_contig <- genomeInfo(GENOME, "MITO", SDIR)
species <- genomeInfo(GENOME, "SPECIES", SDIR)
msigdb_db <- genomeInfo(GENOME, "MSIGDB_DB", SDIR)
collections_spec <- opts$collections %||% genomeInfo(GENOME, "MSIGDB_SETS", SDIR)

# Project for output naming, as in diffAnalysisPairwise.R: the last
# directory of the working directory path named Proj_* or B-<digits>,
# else ATAC
project_id <- getwd() |>
  str_split_1("/") |>
  str_subset("^Proj_|^B-\\d+") |>
  tail(1)
if (len(project_id) == 0) {
  project_id <- "ATAC"
}

input_rds <- cc(project_id, RUNTAG, "DiffPeaksEdgeR_V5.rds")
if (!file_exists(input_rds)) {
  stop(str_glue(
    "no {input_rds} in {getwd()}; run R/diffAnalysisPairwise.R with ",
    "RUNTAG [{RUNTAG}] first"
  ), call. = FALSE)
}
diff_results <- readRDS(input_rds)
if (diff_results$genome != GENOME) {
  stop(str_glue("{input_rds} is for {diff_results$genome}, the run is {GENOME}"),
       call. = FALSE)
}

out_prefix <- c(RUNTAG, opts$out_tag) |>
  keep(nzchar) |>
  str_c(collapse = "_")
output_xlsx <- cc(project_id, out_prefix, "PathwayGSEA.xlsx")
output_pdf <- cc(project_id, out_prefix, "PathwayGSEA.pdf")
output_rds <- cc(project_id, out_prefix, "PathwayGSEA.rds")

#' Parse the --collections list
#'
#' @param spec Character, comma separated COLL[:SUBCOLL] items, an item
#'   starting with - removes sets
#' @return Tibble with item, collection, subcollection, exclude
parse_collections <- function(spec) {
  tibble(item = str_trim(str_split_1(spec, ","))) |>
    filter(item != "") |>
    mutate(exclude = str_starts(item, "-"),
           item = str_remove(item, "^-"),
           collection = str_extract(item, "^[^:]+"),
           subcollection = str_remove(item, "^[^:]+:?"))
}

#' Select the gene sets of one collection item
#'
#' An empty subcollection selects the whole collection; otherwise the
#' subcollection itself and every one below it (CP matches CP:REACTOME).
#'
#' @param gene_sets Tibble from msigdbr
#' @param subcollection Character subcollection, may be ""
#' @return Logical vector over the rows of gene_sets
in_subcollection <- function(gene_sets, subcollection) {
  subcollection == "" |
    gene_sets$gs_subcollection == subcollection |
    str_starts(gene_sets$gs_subcollection, str_glue("{subcollection}:"))
}

collections <- parse_collections(collections_spec)
if (!any(!collections$exclude)) {
  stop("--collections names no collection to test: ", collections_spec,
       call. = FALSE)
}

# Every item must name gene sets of this MSigDB; an unknown one stops
# the run with the list of valid names
known <- msigdbr_collections(db_species = msigdb_db)
msigdb_version <- unique(known$db_version)
for (i in seq_len(nrow(collections))) {
  item <- collections[i, ]
  hits <- known |>
    filter(gs_collection == item$collection,
           in_subcollection(known, item$subcollection))
  if (nrow(hits) == 0) {
    valid <- known |>
      mutate(name = if_else(gs_subcollection == "", gs_collection,
                            str_c(gs_collection, ":", gs_subcollection))) |>
      pull(name)
    stop(str_glue(
      "unknown MSigDB collection [{item$item}] for {msigdb_db} ",
      "({msigdb_version}); valid: {str_c(valid, collapse = ' ')}"
    ), call. = FALSE)
  }
}

cat("MSigDB =", msigdb_version, "\n")
cat("Collections =", collections_spec, "\n")

# Gene sets as Entrez IDs, the gene IDs of the UCSC knownGene TxDb
gene_sets_all <- collections$collection |>
  unique() |>
  map(\(coll) {
    msigdbr(db_species = msigdb_db, species = species, collection = coll)
  }) |>
  bind_rows() |>
  distinct(gs_collection, gs_subcollection, gs_name, ncbi_gene)

excluded <- collections |> filter(exclude)
for (i in seq_len(nrow(excluded))) {
  gene_sets_all <- gene_sets_all |>
    filter(!(gs_collection == excluded$collection[i] &
               in_subcollection(gene_sets_all, excluded$subcollection[i])))
}

tested_collections <- collections |> filter(!exclude)
gene_set_lists <- map(seq_len(nrow(tested_collections)), \(i) {
  sets <- gene_sets_all |>
    filter(gs_collection == tested_collections$collection[i],
           in_subcollection(gene_sets_all, tested_collections$subcollection[i]))
  list(
    pathways = split(sets$ncbi_gene, sets$gs_name),
    subcollection = sets |> distinct(gs_name, gs_subcollection)
  )
})
names(gene_set_lists) <- tested_collections$item

# Nearest gene of every tested peak, once for all comparisons. The SAF
# Start is 1-based, as GRanges. Mitochondrial peaks have no gene to
# rank
tested_peaks <- diff_results$comparisons |>
  map("table") |>
  map("PeakNo") |>
  unlist() |>
  unique()

peaks <- diff_results$peaks |>
  filter(PeakNo %in% tested_peaks, Chr != mito_contig)
if (add_chr) {
  peaks <- peaks |> mutate(Chr = str_c("chr", Chr))
}
peaks_gr <- peaks |>
  dplyr::select(PeakNo, Chr, Start, End) |>
  makeGRangesFromDataFrame(seqnames.field = "Chr", keep.extra.columns = TRUE)

cat("Annotating", format(len(peaks_gr), big.mark = ","), "peaks ...\n")
max_distance <- if (opts$mode == "promoter") opts$tss_window else opts$max_distance
peak_genes <- annotatePeak(peaks_gr, TxDb = txdb,
                           tssRegion = c(-opts$tss_window, opts$tss_window),
                           annoDb = annoDb, verbose = FALSE) |>
  as.data.frame() |>
  tibble() |>
  filter(!is.na(geneId), abs(distanceToTSS) <= max_distance) |>
  dplyr::select(PeakNo, Chr = seqnames, Start = start, End = end, geneId,
                SYMBOL, distanceToTSS)
cat("Peaks within", max_distance, "bp of the TSS of their nearest gene:",
    format(nrow(peak_genes), big.mark = ","), "\n")

#' Rank the genes of one comparison
#'
#' Each gene gets the statistic sign(logFC) * -log10(PValue) of its
#' peak with the largest absolute value (ties: the peak nearest the
#' TSS). A p-value of 0 (an underflow) is set to the smallest nonzero
#' one, as on the volcano plot.
#'
#' @param results Tibble of every tested peak (PeakNo, logFC, PValue, ...)
#' @return Tibble with one row per gene, highest stat first
rank_genes <- function(results) {
  p_floor <- min(results$PValue[results$PValue > 0])
  results |>
    inner_join(peak_genes, by = "PeakNo") |>
    mutate(stat = sign(logFC) * -log10(pmax(PValue, p_floor))) |>
    arrange(desc(abs(stat)), abs(distanceToTSS)) |>
    distinct(geneId, .keep_all = TRUE) |>
    arrange(desc(stat)) |>
    dplyr::select(geneId, SYMBOL, stat, PeakNo, Chr, Start, End,
                  distanceToTSS, logFC, logCPM, PValue, FDR)
}

#' Run fgsea for one comparison over every collection
#'
#' @param comp One comparison of the input RDS (tag, label, group1,
#'   group2, table)
#' @return List with the ranks and the results table
run_gsea <- function(comp) {
  cat("========================================================\n")
  cat(comp$label, "\n")

  ranks_tbl <- rank_genes(comp$table)
  ranks <- set_names(ranks_tbl$stat, ranks_tbl$geneId)
  symbols <- set_names(ranks_tbl$SYMBOL, ranks_tbl$geneId)
  cat("Genes ranked =", format(len(ranks), big.mark = ","), "\n")

  results <- imap(gene_set_lists, \(sets, coll) {
    set.seed(GSEA_SEED)
    res <- fgseaMultilevel(sets$pathways, ranks, minSize = opts$min_size,
                           maxSize = opts$max_size, eps = 0, nproc = 1)
    cat("  ", coll, ":", nrow(res), "sets,",
        sum(res$padj < FDR_CUT, na.rm = TRUE), "with padj <", FDR_CUT, "\n")
    res |>
      tibble() |>
      mutate(Collection = coll,
             leadingEdge = map_chr(leadingEdge, \(ids) {
               str_c(coalesce(unname(symbols[ids]), ids), collapse = ",")
             })) |>
      left_join(sets$subcollection, by = c("pathway" = "gs_name"))
  }) |>
    bind_rows() |>
    mutate(Collection = factor(Collection, levels = names(gene_set_lists))) |>
    arrange(Collection, pval) |>
    mutate(Collection = as.character(Collection)) |>
    dplyr::select(Collection, Subcollection = gs_subcollection, pathway, NES,
                  pval, padj, log2err, ES, size, leadingEdge)

  list(tag = comp$tag, label = comp$label, group1 = comp$group1,
       group2 = comp$group2, ranks = ranks_tbl, results = results)
}

gsea_list <- map(diff_results$comparisons, run_gsea)
cat("========================================================\n\n")

# Summary: per comparison and collection, sets tested and significant
summary_stats <- gsea_list |>
  map(\(g) {
    g$results |>
      group_by(Collection) |>
      summarize(SetsTested = n(),
                SigPositiveNES = sum(padj < FDR_CUT & NES > 0, na.rm = TRUE),
                SigNegativeNES = sum(padj < FDR_CUT & NES < 0, na.rm = TRUE),
                .groups = "drop") |>
      mutate(Comparison = g$tag, GenesRanked = nrow(g$ranks), .before = 1)
  }) |>
  bind_rows() |>
  mutate(Collection = factor(Collection, levels = names(gene_set_lists))) |>
  arrange(match(Comparison, map_chr(gsea_list, "tag")), Collection) |>
  mutate(Collection = as.character(Collection))

parameters <- tribble(
  ~Parameter, ~Value,
  "Genome", GENOME,
  "Input", input_rds,
  "Mode", opts$mode,
  "Max distance to TSS (bp)", as.character(max_distance),
  "Rank statistic", "sign(logFC) * -log10(PValue), the peak with the largest |stat| per gene",
  "Collections", collections_spec,
  "MSigDB", msigdb_version,
  "msigdbr", as.character(packageVersion("msigdbr")),
  "fgsea", as.character(packageVersion("fgsea")),
  "ChIPseeker", as.character(packageVersion("ChIPseeker")),
  "TxDb", txdb_pkg,
  "Min set size", as.character(opts$min_size),
  "Max set size", as.character(opts$max_size),
  "Seed", as.character(GSEA_SEED),
  "padj", "Benjamini-Hochberg within each collection and comparison",
  "NES sign", "positive: more accessible in Group2 (the second group of the comparison)"
)

# Excel: Summary (150% zoom, as in diffAnalysisPairwise.R), one sheet
# per comparison, Parameters
result_sheets <- map(gsea_list, "results")
names(result_sheets) <- map_chr(gsea_list, "tag") |> substr(1, 31)
xlsx_sheets <- c(list(Summary = summary_stats), result_sheets,
                 list(Parameters = parameters))

wb <- createWorkbook()
iwalk(xlsx_sheets, \(sheet_data, sheet_name) {
  addWorksheet(wb, sheet_name, zoom = if (sheet_name == "Summary") 150 else 100)
  writeData(wb, sheet_name, sheet_data)
})
setColWidths(wb, "Summary", cols = seq_len(ncol(summary_stats)),
             widths = "auto")
setColWidths(wb, "Parameters", cols = 1:2, widths = "auto")
saveWorkbook(wb, output_xlsx, overwrite = TRUE)

saveRDS(list(parameters = parameters, gsea = gsea_list), output_rds)

#' NES bar chart of the top sets of one comparison and collection
#'
#' The N_TOP sets with the smallest p-value on each side of 0. A bar is
#' colored by the group its set is more accessible in when padj <
#' FDR_CUT, grey otherwise, with the colors of the MA and volcano plots.
#'
#' @param res Tibble of fgsea results of one collection
#' @param g One element of gsea_list
#' @param coll Character collection item
#' @return ggplot, or NULL if no set has an NES
plot_top_sets <- function(res, g, coll) {
  if (all(is.na(res$NES))) {
    return(NULL)
  }
  # fgseaMultilevel gives NA (NES, pval, padj) for a set whose p-value
  # it cannot estimate; such sets stay in the tables but are not drawn
  top <- res |>
    filter(!is.na(NES)) |>
    group_by(Positive = NES > 0) |>
    slice_min(pval, n = N_TOP, with_ties = FALSE) |>
    ungroup() |>
    mutate(Status = case_when(padj >= FDR_CUT ~ "ns",
                              NES > 0 ~ "up",
                              TRUE ~ "down"),
           Status = factor(Status, levels = c("ns", "down", "up")),
           Set = fct_reorder(pathway, NES))

  status_colors <- c(ns = "#bdbdbd", down = "#1f78b4", up = "#e31a1c")
  status_labels <- c(ns = "Not significant",
                     down = str_glue("Higher in {g$group1}"),
                     up = str_glue("Higher in {g$group2}"))
  n_sig <- sum(res$padj < FDR_CUT, na.rm = TRUE)
  subtitle <- str_glue("{nrow(res)} sets tested, {n_sig} with padj < ",
                       "{FDR_CUT}; top {N_TOP} by p-value on each side")

  ggplot(top, aes(NES, Set, fill = Status)) +
    theme_light(base_size = 11) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = 0, color = "grey40") +
    scale_fill_manual(values = status_colors, labels = status_labels,
                      name = str_glue("padj < {FDR_CUT}"), drop = FALSE) +
    scale_y_discrete(labels = \(x) str_trunc(x, 70)) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank()) +
    labs(title = str_glue("{g$label}: {coll}"),
         subtitle = subtitle,
         x = str_glue("NES (positive: more accessible in {g$group2})"),
         y = NULL)
}

pdf(output_pdf, width = 11, height = 8.5)
for (g in gsea_list) {
  for (coll in names(gene_set_lists)) {
    p <- plot_top_sets(g$results |> filter(Collection == coll), g, coll)
    if (!is.null(p)) {
      print(p)
    }
  }
}
dev.off()

cat("Wrote", output_xlsx, output_pdf, output_rds, "\n")
cat("## END: pathwayGSEA.R\n")
