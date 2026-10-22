suppressPackageStartupMessages(library(optparse))

usage <- "
  %prog [options] [RUNTAG]

Pathway enrichment (preranked GSEA, fgsea) of the differential peaks
of R/diffAnalysisPairwise.R against MSigDB gene sets (msigdbr).

Run from the analysis directory, after diffAnalysisPairwise.R with the
same RUNTAG (default: none). Its input is that script's
atacSeq/diff/<PROJ>_<RUNTAG>_DiffPeaksEdgeR_V5.rds, which holds every tested peak
of every comparison. The genome is GENOME= in 00.RUNSTATUS.txt; the
annotation, the MSigDB species and the default collections come from
lib/genomes/GENOMES.tsv.

Each peak is assigned to genes by their TSSs (transcript starts of the
TxDb): in promoter mode to every gene with a TSS within --tss-window,
in distance mode to the gene with the nearest TSS (all of them on a
tie) within --max-distance; a peak can count for more than one gene.
For each comparison each gene gets the statistic
sign(logFC) * -log10(PValue) of its tested peak with the largest
absolute value. The genes ranked by it are
tested with fgseaMultilevel against each collection; padj is
Benjamini-Hochberg within the collection. A positive NES means the set
is more accessible in Group2 of the comparison.

COLLECTIONS is a comma separated list of COLL or COLL:SUBCOLL
(msigdbr's gs_collection and gs_subcollection, e.g. H, C2:CP:REACTOME,
C5:GO:BP; mouse: MH, M2:CP:REACTOME, M5:GO:BP). A subcollection also
matches the ones below it: C2:CP is CP and every CP:*. An item
starting with - removes those sets from every other item, e.g.
-C2:CP:KEGG_LEGACY. Each item without - is one fgsea run.

The significant sets (padj < 0.05) of each collection are collapsed
to independent ones (fgsea collapsePathways); the plots show only
independent sets, the xlsx every tested set.

OUTPUT FILES (in atacSeq/diff):
  <PROJ>_<RUNTAG>_PathwayGSEA.xlsx  Summary, one sheet per comparison
                                    (every tested set), Parameters
  <PROJ>_<RUNTAG>_PathwayGSEA.pdf   NES bar chart per comparison and
                                    collection; no page where no set
                                    is significant
  <PROJ>_<RUNTAG>_PathwayGSEA_Dot.pdf
                                    the same sets as a dot plot
  <PROJ>_<RUNTAG>_PathwayGSEA_Heatmap.pdf
                                    NES across comparisons, per
                                    collection (two or more comparisons)
  <PROJ>_<RUNTAG>_PathwayGSEA_Enrichment.pdf
                                    running enrichment score of the
                                    top sets, and a heatmap of each
                                    one's leading-edge genes
  <PROJ>_<RUNTAG>_PathwayGSEA.rds   gene ranks and full results
  With --out-tag T the names end _T_PathwayGSEA*"

option_list <- list(
  make_option(c("-c", "--collections"), type = "character", default = NULL,
              help = "MSigDB collections [default: MSIGDB_SETS of GENOMES.tsv]"),
  make_option(c("-m", "--mode"), type = "character", default = "promoter",
              help = "promoter or distance [default: %default]"),
  make_option(c("-w", "--tss-window"), type = "integer", default = 1000L,
              help = "promoter mode: a peak counts for every gene with a TSS within this many bp [default: %default]"),
  make_option(c("-d", "--max-distance"), type = "integer", default = 50000L,
              help = "distance mode: a peak counts for the gene(s) with the nearest TSS within this many bp [default: %default]"),
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
source(file.path(SDIR, "R", "palettes.R"))

GSEA_SEED <- 42
FDR_CUT <- 0.05
N_TOP <- 15
# A bar chart page has room for at least this many bars, so that one or
# two sets are not drawn as bars a page high
MIN_BAR_SLOTS <- 10
# Running-score plots and leading-edge heatmaps: sets per side of NES 0
# for each comparison and collection, and leading-edge genes drawn
N_ENRICH <- 3
N_LEADING_EDGE <- 50
# NES heatmap across comparisons: sets per collection
N_HEATMAP <- 50

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

out_dir <- "atacSeq/diff"
input_rds <- path(out_dir, cc(project_id, RUNTAG, "DiffPeaksEdgeR_V5.rds"))
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
output_xlsx <- path(out_dir, cc(project_id, out_prefix, "PathwayGSEA.xlsx"))
output_pdf <- path(out_dir, cc(project_id, out_prefix, "PathwayGSEA.pdf"))
output_dot_pdf <- path(out_dir, cc(project_id, out_prefix, "PathwayGSEA_Dot.pdf"))
output_heatmap_pdf <- path(out_dir, cc(project_id, out_prefix, "PathwayGSEA_Heatmap.pdf"))
output_enrich_pdf <- path(out_dir, cc(project_id, out_prefix, "PathwayGSEA_Enrichment.pdf"))
output_rds <- path(out_dir, cc(project_id, out_prefix, "PathwayGSEA.rds"))

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

# The genes of every peak of the atlas, once for all comparisons. A
# gene's TSSs are the starts of its transcripts in the TxDb (gene IDs
# are Entrez IDs, as in msigdbr). In promoter mode a peak belongs to
# every gene with a TSS within --tss-window of it, so a promoter shared
# by two genes counts for both; in distance mode to the gene with the
# nearest TSS, or to all of them on a tie, within --max-distance. No
# gene is chosen over another: ChIPseeker's annotatePeak, used before,
# gave a peak one gene and broke ties in a way that depended on the
# other peaks it was given, so a comparison's results changed with the
# other comparisons of the run, and 5-7% of the ranked genes (most often
# a protein-coding gene and a noncoding RNA at one promoter) depended on
# the tie. The SAF Start is 1-based, as GRanges. Mitochondrial peaks
# have no gene to rank
peaks <- diff_results$peaks |>
  filter(Chr != mito_contig)
if (add_chr) {
  peaks <- peaks |> mutate(Chr = str_c("chr", Chr))
}
peaks_gr <- peaks |>
  dplyr::select(PeakNo, Chr, Start, End) |>
  makeGRangesFromDataFrame(seqnames.field = "Chr", keep.extra.columns = TRUE)

tss_gr <- GenomicFeatures::transcriptsBy(txdb, by = "gene") |>
  unlist() |>
  resize(width = 1, fix = "start")
tss_gr$geneId <- names(tss_gr)
names(tss_gr) <- NULL

cat("Assigning", format(len(peaks_gr), big.mark = ","), "peaks to genes ...\n")
max_distance <- if (opts$mode == "promoter") opts$tss_window else opts$max_distance
# distance() is 0 for a TSS in or next to the peak, and <= d exactly for
# the hits of findOverlaps(maxgap = d)
hits <- if (opts$mode == "promoter") {
  findOverlaps(peaks_gr, tss_gr, maxgap = max_distance, ignore.strand = TRUE)
} else {
  nearest(peaks_gr, tss_gr, select = "all", ignore.strand = TRUE)
}
peak_genes <- tibble(
  PeakNo = peaks_gr$PeakNo[queryHits(hits)],
  geneId = tss_gr$geneId[subjectHits(hits)],
  distanceToTSS = distance(peaks_gr[queryHits(hits)],
                           tss_gr[subjectHits(hits)], ignore.strand = TRUE)
) |>
  filter(distanceToTSS <= max_distance) |>
  group_by(PeakNo, geneId) |>
  summarize(distanceToTSS = min(distanceToTSS), .groups = "drop") |>
  left_join(peaks |> dplyr::select(PeakNo, Chr, Start, End), by = "PeakNo")

org_db <- getExportedValue(annoDb, annoDb)
gene_symbols <- suppressMessages(AnnotationDbi::mapIds(
  org_db, keys = unique(peak_genes$geneId), column = "SYMBOL",
  keytype = "ENTREZID"))
peak_genes <- peak_genes |>
  mutate(SYMBOL = unname(gene_symbols[geneId]))
cat("Peaks within", max_distance, "bp of a TSS",
    if (opts$mode == "distance") "(nearest gene)", ":",
    format(n_distinct(peak_genes$PeakNo), big.mark = ","), "for",
    format(n_distinct(peak_genes$geneId), big.mark = ","), "genes\n")

#' Rank the genes of one comparison
#'
#' Each gene gets the statistic sign(logFC) * -log10(PValue) of its
#' peak with the largest absolute value (ties: the peak nearest the
#' TSS). A peak can be the peak of more than one gene. A p-value of 0
#' (an underflow) is set to the smallest nonzero one, as on the volcano
#' plot.
#'
#' @param results Tibble of every tested peak (PeakNo, logFC, PValue, ...)
#' @return Tibble with one row per gene, highest stat first
rank_genes <- function(results) {
  p_floor <- min(results$PValue[results$PValue > 0])
  results |>
    inner_join(peak_genes, by = "PeakNo", relationship = "one-to-many") |>
    mutate(stat = sign(logFC) * -log10(pmax(PValue, p_floor))) |>
    arrange(desc(abs(stat)), abs(distanceToTSS)) |>
    distinct(geneId, .keep_all = TRUE) |>
    arrange(desc(stat)) |>
    dplyr::select(geneId, SYMBOL, stat, PeakNo, Chr, Start, End,
                  distanceToTSS, logFC, logCPM, PValue, FDR)
}

#' Independent sets among the significant sets of one collection
#'
#' fgsea's collapsePathways on the sets with padj < FDR_CUT, smallest
#' p-value first: a set whose enrichment is no longer significant (p >
#' 0.05) given a more significant set is redundant and collapsed into
#' that set. It permutes, so it runs with the seed.
#'
#' @param res Tibble of fgsea results of one collection
#' @param pathways Named list of the collection's gene sets (Entrez IDs)
#' @param ranks Named vector of the rank statistics
#' @return Tibble with pathway and CollapsedInto (NA for an independent
#'   set), one row per significant set
collapse_sets <- function(res, pathways, ranks) {
  sig <- res |>
    filter(!is.na(padj), padj < FDR_CUT) |>
    arrange(pval)
  if (nrow(sig) < 2) {
    return(tibble(pathway = sig$pathway, CollapsedInto = NA_character_))
  }
  set.seed(GSEA_SEED)
  collapsed <- collapsePathways(data.table::as.data.table(sig), pathways,
                                ranks)
  tibble(pathway = names(collapsed$parentPathways),
         CollapsedInto = as.character(unname(collapsed$parentPathways)))
}

#' Run fgsea for one comparison over every collection
#'
#' @param comp One comparison of the input RDS (tag, label, group1,
#'   group2, table)
#' @return List with the ranks and the results table; leadingEdgeIds
#'   holds the leading edge as Entrez IDs, for the plots
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
                           maxSize = opts$max_size, eps = 0, nproc = 1) |>
      tibble()
    collapsed <- collapse_sets(res, sets$pathways, ranks)
    cat("  ", coll, ":", nrow(res), "sets,", nrow(collapsed), "with padj <",
        str_c(FDR_CUT, ","), sum(is.na(collapsed$CollapsedInto)),
        "independent\n")
    res |>
      left_join(collapsed, by = "pathway") |>
      mutate(Collection = coll,
             leadingEdgeIds = leadingEdge,
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
                  pval, padj, log2err, ES, size, CollapsedInto, leadingEdge,
                  leadingEdgeIds)

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
                SigIndependent = sum(padj < FDR_CUT & is.na(CollapsedInto),
                                     na.rm = TRUE),
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
  "Peak to gene", if (opts$mode == "promoter") "every gene with a TSS within the max distance" else "the gene(s) with the nearest TSS, within the max distance",
  "GenomicRanges", as.character(packageVersion("GenomicRanges")),
  "TxDb", txdb_pkg,
  "Min set size", as.character(opts$min_size),
  "Max set size", as.character(opts$max_size),
  "Seed", as.character(GSEA_SEED),
  "padj", "Benjamini-Hochberg within each collection and comparison",
  "Comparison sheets", as.character(str_glue("every tested set, NOT filtered for padj < {FDR_CUT}")),
  "Collapsing", as.character(str_glue("fgsea collapsePathways on the sets with padj < {FDR_CUT}; CollapsedInto names the set that explains a redundant one; the plots show only independent sets")),
  "NES sign", "positive: more accessible in Group2 (the second group of the comparison)"
)

# Excel: Summary (150% zoom, as in diffAnalysisPairwise.R), one sheet
# per comparison, Parameters
result_sheets <- map(gsea_list, \(g) g$results |> dplyr::select(-leadingEdgeIds))
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

# The plots show only significant, independent sets: padj < FDR_CUT and
# not collapsed into another set. fgseaMultilevel gives NA (NES, pval,
# padj) for a set whose p-value it cannot estimate; such sets stay in
# the tables but are never drawn.

#' The significant, independent sets of one collection to plot
#'
#' @param res Tibble of fgsea results of one collection
#' @param n Maximum sets per side of NES 0, smallest p-value first
#' @return Tibble of at most 2 n sets
top_sets <- function(res, n) {
  res |>
    filter(!is.na(padj), padj < FDR_CUT, is.na(CollapsedInto)) |>
    group_by(Positive = NES > 0) |>
    slice_min(pval, n = n, with_ties = FALSE) |>
    ungroup()
}

#' Subtitle with the counts of tested, significant and independent sets
#'
#' @param res Tibble of fgsea results of one collection
#' @param g One element of gsea_list
#' @return Character
counts_subtitle <- function(res, g) {
  sig <- res |> filter(!is.na(padj), padj < FDR_CUT)
  ind <- sig |> filter(is.na(CollapsedInto))
  str_glue("{nrow(res)} sets tested, {nrow(sig)} with padj < {FDR_CUT}, ",
           "{nrow(ind)} independent ({sum(ind$NES > 0)} higher in ",
           "{g$group2}, {sum(ind$NES < 0)} in {g$group1})")
}

#' Set axis of the bar and dot charts
#'
#' ggplot's default discrete expansion is 0.6 on each side; fewer than
#' MIN_BAR_SLOTS sets get the missing slots as space below them, so a
#' bar is never thicker than with MIN_BAR_SLOTS bars.
#'
#' @param n_sets Number of sets drawn
#' @return ggplot scale
set_axis <- function(n_sets) {
  n_pad <- max(0, MIN_BAR_SLOTS - n_sets)
  scale_y_discrete(labels = \(x) str_trunc(x, 70),
                   expand = expansion(add = c(0.6 + n_pad, 0.6)))
}

#' NES bar chart of the independent significant sets of one comparison
#' and collection
#'
#' At most N_TOP sets with the smallest p-value on each side of 0. A bar
#' is colored by the group its set is more accessible in
#' (SIDE_COLORS).
#'
#' @param res Tibble of fgsea results of one collection
#' @param g One element of gsea_list
#' @param coll Character collection item
#' @return ggplot, or NULL if no set has padj < FDR_CUT
plot_bar_sets <- function(res, g, coll) {
  top <- top_sets(res, N_TOP)
  if (nrow(top) == 0) {
    return(NULL)
  }
  top <- top |>
    mutate(Status = if_else(NES > 0, "up", "down"),
           Status = factor(Status, levels = c("down", "up")),
           Set = fct_reorder(pathway, NES))

  status_labels <- c(down = str_glue("Higher in {g$group1}"),
                     up = str_glue("Higher in {g$group2}"))

  ggplot(top, aes(NES, Set, fill = Status)) +
    theme_light(base_size = 11) +
    # show.legend: keys for both groups even when one has no bar
    geom_col(width = 0.7, show.legend = TRUE) +
    geom_vline(xintercept = 0, color = REF_LINE) +
    scale_fill_manual(values = SIDE_COLORS, labels = status_labels,
                      name = str_glue("padj < {FDR_CUT}"), drop = FALSE) +
    set_axis(nrow(top)) +
    theme(legend.position = "bottom", panel.grid.major.y = element_blank()) +
    labs(title = str_glue("{g$label}: {coll}"),
         subtitle = str_glue("{counts_subtitle(res, g)}\n",
                             "At most {N_TOP} per side drawn"),
         x = str_glue("NES (positive: more accessible in {g$group2})"),
         y = NULL)
}

#' NES dot plot of the same sets as plot_bar_sets
#'
#' The dot size is the number of ranked genes in the set, the fill
#' -log10(padj) on PAL_QC, clipped to its 5th and 95th percentiles. The
#' dots have an outline, so the palest still shows on white.
#'
#' @param res Tibble of fgsea results of one collection
#' @param g One element of gsea_list
#' @param coll Character collection item
#' @return ggplot, or NULL if no set has padj < FDR_CUT
plot_dot_sets <- function(res, g, coll) {
  top <- top_sets(res, N_TOP)
  if (nrow(top) == 0) {
    return(NULL)
  }
  top <- top |> mutate(Set = fct_reorder(pathway, NES))

  ggplot(top, aes(NES, Set)) +
    theme_light(base_size = 11) +
    geom_vline(xintercept = 0, color = REF_LINE) +
    geom_point(aes(size = size, fill = -log10(padj)), shape = 21,
               color = POINT_OUTLINE, stroke = 0.3) +
    scale_sequential(-log10(top$padj), PAL_QC, name = "-log10(padj)") +
    scale_size(name = "Genes in set", range = c(2, 8)) +
    set_axis(nrow(top)) +
    theme(legend.position = "right") +
    labs(title = str_glue("{g$label}: {coll}"),
         subtitle = str_glue("{counts_subtitle(res, g)}\n",
                             "At most {N_TOP} per side drawn; color ",
                             "clipped to the 5th-95th percentile"),
         x = str_glue("NES (positive: more accessible in {g$group2})"),
         y = NULL)
}

#' NES of one collection's sets across all comparisons
#'
#' Rows: the sets significant and independent in at least one
#' comparison, at most N_HEATMAP with the smallest padj, clustered on
#' their NES. A cell is colored by the NES where the set has padj <
#' FDR_CUT in that comparison (PAL_SIGNED, symmetric at +/- the 95th
#' percentile of |NES|), grey otherwise or where it was not tested.
#'
#' @param coll Character collection item
#' @return ggplot, or NULL if no set is significant in any comparison
plot_nes_heatmap <- function(coll) {
  comp_labels <- map_chr(gsea_list, "label")
  all_res <- gsea_list |>
    map(\(g) g$results |> filter(Collection == coll) |>
          mutate(Comparison = g$label)) |>
    bind_rows()

  candidates <- all_res |>
    filter(!is.na(padj), padj < FDR_CUT, is.na(CollapsedInto)) |>
    group_by(pathway) |>
    summarize(minPadj = min(padj), .groups = "drop")
  if (nrow(candidates) == 0) {
    return(NULL)
  }
  rows <- candidates |>
    slice_min(minPadj, n = N_HEATMAP, with_ties = FALSE)

  cells <- all_res |>
    filter(pathway %in% rows$pathway) |>
    mutate(NES_sig = if_else(!is.na(padj) & padj < FDR_CUT, NES, NA_real_),
           Comparison = factor(Comparison, levels = comp_labels)) |>
    dplyr::select(pathway, Comparison, NES_sig) |>
    complete(pathway, Comparison)

  # Cluster the rows on the NES, 0 where not significant
  nes_wide <- cells |>
    mutate(NES_sig = coalesce(NES_sig, 0)) |>
    pivot_wider(names_from = Comparison, values_from = NES_sig)
  nes_mat <- as.matrix(nes_wide[, -1])
  rownames(nes_mat) <- nes_wide$pathway
  row_order <- if (nrow(nes_mat) > 1) {
    rownames(nes_mat)[hclust(dist(nes_mat))$order]
  } else {
    rownames(nes_mat)
  }

  shown <- if (nrow(candidates) > N_HEATMAP) {
    str_glue(", the {N_HEATMAP} with the smallest padj shown")
  } else {
    ""
  }
  ggplot(cells, aes(Comparison, factor(pathway, levels = row_order),
                    fill = NES_sig)) +
    theme_minimal(base_size = 9) +
    geom_tile(color = TILE_BORDER) +
    scale_signed(cells$NES_sig, name = "NES") +
    scale_y_discrete(labels = \(x) str_trunc(x, 70)) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank()) +
    labs(title = str_glue("NES across comparisons: {coll}"),
         subtitle = str_glue(
           "{nrow(candidates)} sets with padj < {FDR_CUT} and independent ",
           "in at least one comparison{shown}.\nGrey: padj >= {FDR_CUT} ",
           "or not tested. Positive NES: more accessible in the first ",
           "group of the comparison. Color clipped at the 95th ",
           "percentile of |NES|."),
         x = NULL, y = NULL)
}

#' Running enrichment score of one set
#'
#' fgsea's plotEnrichment, redrawn from plotEnrichmentData with the
#' colors of R/palettes.R: the curve in the color of the set's side
#' (SIDE_COLORS), the maximum and minimum dashed.
#'
#' @param set Tibble row of fgsea results
#' @param stats Named vector of the rank statistics
#' @param coll Character collection item
#' @return ggplot
plot_running_score <- function(set, stats, coll) {
  pd <- plotEnrichmentData(gene_set_lists[[coll]]$pathways[[set$pathway]],
                           stats)
  side_color <- SIDE_COLORS[[if_else(set$NES > 0, "up", "down")]]
  ggplot(pd$curve, aes(rank, ES)) +
    theme_light(base_size = 9) +
    geom_hline(yintercept = c(pd$posES, pd$negES), color = REF_LINE,
               linetype = "dashed") +
    geom_hline(yintercept = 0, color = REF_LINE) +
    geom_segment(data = pd$ticks,
                 aes(x = rank, xend = rank, y = -pd$spreadES / 16,
                     yend = pd$spreadES / 16),
                 linewidth = 0.2) +
    geom_line(color = side_color) +
    labs(title = str_trunc(set$pathway, 50),
         subtitle = str_glue("NES {round(set$NES, 2)}, ",
                             "padj {signif(set$padj, 2)}"),
         x = "Gene rank", y = "Enrichment score") +
    theme(plot.title = element_text(size = 9),
          plot.subtitle = element_text(size = 8))
}

#' Running-score plots of the top sets of one comparison and collection
#'
#' At most N_ENRICH sets per side of NES 0, one row per side with sets:
#' first the sets higher in group2, then those higher in group1. A page
#' always has room for two rows, so the plots keep their size.
#'
#' @param top Tibble from top_sets
#' @param g One element of gsea_list
#' @param coll Character collection item
#' @return cowplot page
plot_running_page <- function(top, g, coll) {
  stats <- set_names(g$ranks$stat, g$ranks$geneId)
  sides <- tibble(Positive = c(TRUE, FALSE),
                  Group = c(g$group2, g$group1)) |>
    filter(Positive %in% top$Positive)
  row_plots <- map(sides$Positive, \(positive) {
    sets <- top |> filter(Positive == positive) |> arrange(pval)
    plots <- map(seq_len(nrow(sets)), \(i) {
      plot_running_score(sets[i, ], stats, coll)
    })
    c(plots, rep(list(NULL), N_ENRICH - len(plots)))
  })
  if (len(row_plots) == 1) {
    row_plots <- c(row_plots, list(rep(list(NULL), N_ENRICH)))
  }
  rows_text <- if (nrow(sides) == 2) {
    str_glue("top row higher in {g$group2}, bottom row higher in {g$group1}")
  } else {
    str_glue("higher in {sides$Group}")
  }
  title <- cowplot::ggdraw() +
    cowplot::draw_label(
      str_glue("{g$label}: {coll}. Running enrichment score; {rows_text}"),
      size = 12, hjust = 0, x = 0.01)
  grid_plots <- cowplot::plot_grid(plotlist = list_flatten(row_plots),
                                   nrow = 2, ncol = N_ENRICH)
  cowplot::plot_grid(title, grid_plots, ncol = 1, rel_heights = c(0.06, 1))
}

#' Heatmap of the leading-edge genes of one set
#'
#' Each row is a leading-edge gene, by the peak that gave it its rank,
#' at most N_LEADING_EDGE with the largest |rank statistic| (largest at
#' the top). Each peak's logCPM of the comparison's samples, the values
#' edgeR tested, z-scored across the samples.
#'
#' @param set Tibble row of fgsea results
#' @param g One element of gsea_list
#' @param comp The matching comparison of the input RDS (log_cpm,
#'   samples)
#' @param coll Character collection item
#' @return ggplot
plot_leading_edge <- function(set, g, comp, coll) {
  le_genes <- g$ranks |>
    filter(geneId %in% set$leadingEdgeIds[[1]]) |>
    mutate(Label = coalesce(SYMBOL, geneId))
  shown <- le_genes |>
    slice_max(abs(stat), n = N_LEADING_EDGE, with_ties = FALSE) |>
    arrange(abs(stat))

  samples <- comp$samples |>
    mutate(Group = factor(Group, levels = c(g$group1, g$group2))) |>
    arrange(Group)
  # Rows by gene: two genes at one promoter share their peak, and show
  # the same values
  log_cpm <- comp$log_cpm[shown$PeakNo, samples$SampleID, drop = FALSE]
  rownames(log_cpm) <- shown$geneId
  z <- t(scale(t(log_cpm)))
  z[is.nan(z)] <- 0

  cells <- as_tibble(z, rownames = "geneId") |>
    pivot_longer(-geneId, names_to = "SampleID", values_to = "z") |>
    left_join(samples, by = "SampleID") |>
    mutate(geneId = factor(geneId, levels = shown$geneId),
           SampleID = factor(SampleID, levels = samples$SampleID))

  ggplot(cells, aes(SampleID, geneId, fill = z)) +
    theme_minimal(base_size = 9) +
    geom_tile() +
    facet_grid(~Group, scales = "free_x", space = "free_x") +
    scale_signed(cells$z, name = "logCPM\nz-score") +
    scale_y_discrete(labels = set_names(shown$Label, shown$geneId)) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank()) +
    labs(title = str_glue("{g$label}: {str_trunc(set$pathway, 70)}"),
         subtitle = str_glue(
           "{coll}; NES {round(set$NES, 2)}, padj {signif(set$padj, 2)}. ",
           "Leading edge: {nrow(le_genes)} genes, ",
           "{nrow(shown)} with the largest |rank statistic| shown.\n",
           "Each row is the gene's ranked peak; its logCPM is z-scored ",
           "across the samples; color clipped at the 95th percentile ",
           "of |z|."),
         x = NULL, y = NULL)
}

pdf(output_pdf, width = 11, height = 8.5)
for (g in gsea_list) {
  for (coll in names(gene_set_lists)) {
    p <- plot_bar_sets(g$results |> filter(Collection == coll), g, coll)
    if (!is.null(p)) {
      print(p)
    }
  }
}
dev.off()

pdf(output_dot_pdf, width = 11, height = 8.5)
for (g in gsea_list) {
  for (coll in names(gene_set_lists)) {
    p <- plot_dot_sets(g$results |> filter(Collection == coll), g, coll)
    if (!is.null(p)) {
      print(p)
    }
  }
}
dev.off()

# The heatmap compares comparisons; with one there is nothing to compare.
# A heatmap PDF of an earlier run with the same names is removed, so it
# is not delivered with these results
if (len(gsea_list) > 1) {
  pdf(output_heatmap_pdf, width = 11, height = 8.5)
  for (coll in names(gene_set_lists)) {
    p <- plot_nes_heatmap(coll)
    if (!is.null(p)) {
      print(p)
    }
  }
  dev.off()
} else {
  cat("One comparison: no NES heatmap across comparisons\n")
  if (file_exists(output_heatmap_pdf)) {
    file_delete(output_heatmap_pdf)
    cat("Removed", output_heatmap_pdf, "of an earlier run\n")
  }
}

# The leading-edge heatmaps need each comparison's logCPM, which
# diffAnalysisPairwise.R saves in its RDS since the GSEA plots were
# added; with an older RDS they are left out
has_log_cpm <- diff_results$comparisons |>
  map_lgl(\(comp) !is.null(comp$log_cpm)) |>
  all()
if (!has_log_cpm) {
  cat("WARNING:", input_rds, "has no logCPM; no leading-edge heatmaps.",
      "Rerun R/diffAnalysisPairwise.R to get them\n")
}

pdf(output_enrich_pdf, width = 11, height = 8.5)
for (i in seq_along(gsea_list)) {
  g <- gsea_list[[i]]
  comp <- diff_results$comparisons[[i]]
  for (coll in names(gene_set_lists)) {
    top <- top_sets(g$results |> filter(Collection == coll), N_ENRICH)
    if (nrow(top) == 0) {
      next
    }
    print(plot_running_page(top, g, coll))
    if (has_log_cpm) {
      top <- top |> arrange(desc(Positive), pval)
      for (j in seq_len(nrow(top))) {
        print(plot_leading_edge(top[j, ], g, comp, coll))
      }
    }
  }
}
dev.off()

cat("Wrote", output_xlsx, output_pdf, output_dot_pdf,
    if (len(gsea_list) > 1) output_heatmap_pdf, output_enrich_pdf,
    output_rds, "\n")
cat("## END: pathwayGSEA.R\n")
