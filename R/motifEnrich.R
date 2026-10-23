suppressPackageStartupMessages(library(optparse))

usage <- "
  %prog [options] [RUNTAG]

Motif enrichment (monaLisa, JASPAR 2024) in the differential peaks of
R/diffAnalysisPairwise.R.

Run from the analysis directory, after diffAnalysisPairwise.R with the
same RUNTAG (default: none), as a job through motifEnrich.sh at the
top of this repo; motifEnrich.sh -h gives the bsub line and how to
choose --width, --min-peaks and the job size. Its input is
diffAnalysisPairwise.R's atacSeq/diff/<PROJ>_<RUNTAG>_DiffPeaksEdgeR_V5.rds, which holds every tested peak
of every comparison, and the per-sample MACS2 summits in atacSeq/macs.
The genome is GENOME= in 00.RUNSTATUS.txt; its FASTA comes from
lib/genomes/GENOMES.tsv. The motifs are the JASPAR 2024 CORE
vertebrate non-redundant profiles of lib/motifs.

Each peak of the atlas is represented by a fixed-width region,
--width bp centered on its strongest summit: of the MACS2 summits of
every sample inside the peak, the one with the largest -log10 p. For
each comparison two sets of peaks are tested: up (FDR < 0.05,
logFC > 0, more accessible in Group2) and down (FDR < 0.05,
logFC < 0, more accessible in Group1), the significance of
diffAnalysisPairwise.R. The background of both is the comparison's
other tested peaks (FDR >= 0.05). A set with fewer than --min-peaks
peaks is not tested.

monaLisa's calcBinnedMotifEnrR weights the background sequences to the
G+C content and k-mer composition of the set, then tests each motif
for more set sequences with a hit than expected from the weighted
background (Fisher's exact test, one-sided). padj is Benjamini-Hochberg
within the set (comparison and direction).

OUTPUT FILES (in atacSeq/diff):
  <PROJ>_<RUNTAG>_MotifEnrichment.xlsx  Summary, one sheet per comparison
                                        (every tested motif), Parameters
  <PROJ>_<RUNTAG>_MotifEnrichment.pdf   heatmap of the top motifs of each
                                        comparison, and of all
                                        comparisons together; no page
                                        where no motif is significant
  <PROJ>_<RUNTAG>_MotifEnrichment.rds   regions and full results"

option_list <- list(
  make_option(c("-w", "--width"), type = "integer", default = 501L,
              help = "width of the region around each peak's strongest summit, odd [default: %default]"),
  make_option(c("-m", "--min-peaks"), type = "integer", default = 50L,
              help = "smallest set of peaks tested [default: %default]"),
  make_option(c("-c", "--cores"), type = "integer",
              default = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")),
              help = "cores for the motif scan [default: SLURM_CPUS_PER_TASK, else 1]")
)
parsed <- parse_args(OptionParser(usage = usage, option_list = option_list),
                     positional_arguments = c(0, 1),
                     convert_hyphens_to_underscores = TRUE)
opts <- parsed$options
RUNTAG <- if (len(parsed$args) == 1) parsed$args[1] else ""

if (is.na(opts$width) || opts$width < 51 || opts$width %% 2 != 1) {
  stop("--width must be an odd number of at least 51, not ", opts$width,
       call. = FALSE)
}
if (is.na(opts$min_peaks) || opts$min_peaks < 1) {
  stop("--min-peaks must be at least 1, not ", opts$min_peaks, call. = FALSE)
}
if (is.na(opts$cores) || opts$cores < 1) {
  stop("--cores must be at least 1, not ", opts$cores, call. = FALSE)
}

cat("## START: motifEnrich.R\n")

cat("Loading libraries ...")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(Rsamtools)
  library(TFBSTools)
  library(monaLisa)
  library(BiocParallel)
  library(SummarizedExperiment)
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

FDR_CUT <- 0.05
# monaLisa's defaults: a motif hit scores at least MIN_SCORE (log2
# odds, bits); a sequence with more than MAX_FRAC_N N is dropped; the
# background is weighted for k-mers up to MAX_KMER
MIN_SCORE <- 10
MAX_FRAC_N <- 0.7
MAX_KMER <- 3L
# TFBSTools toPWM default: pseudocount added to each count of a PFM
PWM_PSEUDOCOUNT <- 0.8
# Per-comparison heatmap: motifs per set; across comparisons: motifs
N_TOP <- 20
N_HEATMAP <- 50

JASPAR_FILE <- path(SDIR, "lib", "motifs",
                    "JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt")
JASPAR_NAME <- "JASPAR 2024 CORE vertebrates non-redundant"

# Genome of the run and its FASTA, from GENOMES.tsv
GENOME <- runGenome()
cat("GENOME =", GENOME, "\n")
fasta_file <- genomeInfo(GENOME, "FASTA", SDIR)
genome_file <- path(SDIR, "lib", "genomes", str_glue("{GENOME}.genome"))

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

output_xlsx <- path(out_dir, cc(project_id, RUNTAG, "MotifEnrichment.xlsx"))
output_pdf <- path(out_dir, cc(project_id, RUNTAG, "MotifEnrichment.pdf"))
output_rds <- path(out_dir, cc(project_id, RUNTAG, "MotifEnrichment.rds"))

# The FASTA must hold every contig of the genome file with the same
# length, or the regions would be read from another assembly
contigs <- read_tsv(genome_file, col_names = c("Chr", "ContigLength"),
                    col_types = "ci")
fai_file <- str_c(fasta_file, ".fai")
if (!file_exists(fai_file)) {
  stop(str_glue("no FASTA index {fai_file}"), call. = FALSE)
}
fai <- read_tsv(fai_file, col_names = c("Chr", "FastaLength"),
                col_types = "ci---")
bad_contigs <- contigs |>
  left_join(fai, by = "Chr") |>
  filter(is.na(FastaLength) | FastaLength != ContigLength)
if (nrow(bad_contigs) > 0) {
  stop(str_glue(
    "{fasta_file} does not match {genome_file}: ",
    "{str_c(bad_contigs$Chr, collapse = ' ')} missing or of another length"
  ), call. = FALSE)
}
cat("FASTA =", fasta_file, "\n")

# The peaks tested in at least one comparison
tested_peaks <- diff_results$comparisons |>
  map(\(comp) comp$table$PeakNo) |>
  unlist() |>
  unique()
peaks <- diff_results$peaks |>
  filter(PeakNo %in% tested_peaks)
# The SAF Start is 1-based, as GRanges
peaks_gr <- peaks |>
  dplyr::select(PeakNo, Chr, Start, End) |>
  makeGRangesFromDataFrame(seqnames.field = "Chr", keep.extra.columns = TRUE)

# The summits of every sample's blacklist-filtered peaks, the peaks the
# atlas was merged from. Column 5 is MACS2's -log10 p at the summit;
# the BED start is 0-based
summit_files <- dir_ls("atacSeq/macs", recurse = TRUE,
                       glob = "*_summits.blacklistFiltered.bed")
if (len(summit_files) == 0) {
  stop("no *_summits.blacklistFiltered.bed in atacSeq/macs", call. = FALSE)
}
cat("Reading", len(summit_files), "summit files ...")
summits <- summit_files |>
  map(\(f) read_tsv(f, col_names = c("Chr", "Start0", "End", "Name", "Score"),
                    col_types = "cii-d")) |>
  bind_rows()
cat(" done,", format(nrow(summits), big.mark = ","), "summits\n")
summits_gr <- GRanges(summits$Chr, IRanges(summits$Start0 + 1, width = 1),
                      Score = summits$Score)

# Each peak's strongest summit; ties go to the leftmost. A peak with no
# summit (the merge removes the blacklisted part of a merged peak, which
# can leave a piece without one) is centered on its middle
hits <- findOverlaps(peaks_gr, summits_gr)
strongest <- tibble(PeakNo = peaks_gr$PeakNo[queryHits(hits)],
                    Summit = start(summits_gr)[subjectHits(hits)],
                    Score = summits_gr$Score[subjectHits(hits)]) |>
  arrange(PeakNo, desc(Score), Summit) |>
  distinct(PeakNo, .keep_all = TRUE)
regions <- peaks |>
  left_join(strongest, by = "PeakNo") |>
  mutate(HasSummit = !is.na(Summit),
         Summit = coalesce(Summit, as.integer((Start + End) %/% 2)))
n_no_summit <- sum(!regions$HasSummit)
cat("Peaks:", format(nrow(regions), big.mark = ","), "tested,",
    n_no_summit, "without a summit (centered)\n")

# The region: --width bp centered on the summit, moved inside the contig
# where it would run past an end, so every region has the same width
half_width <- (opts$width - 1) %/% 2
regions <- regions |>
  left_join(contigs, by = "Chr") |>
  mutate(RegionStart = pmax(1L, Summit - half_width),
         RegionStart = pmin(RegionStart, ContigLength - opts$width + 1L),
         RegionEnd = RegionStart + opts$width - 1L) |>
  dplyr::select(-ContigLength)

cat("Reading", format(nrow(regions), big.mark = ","), "sequences ...")
regions_gr <- GRanges(regions$Chr,
                      IRanges(regions$RegionStart, regions$RegionEnd))
seqs <- getSeq(FaFile(fasta_file), regions_gr)
names(seqs) <- regions$PeakNo
cat(" done\n")

# monaLisa drops a sequence with more than MAX_FRAC_N N; drop it here,
# so the set and background sizes reported are those tested
frac_n <- as.vector(letterFrequency(seqs, "N")) / opts$width
regions <- regions |> mutate(FracN = frac_n)
n_many_n <- sum(frac_n > MAX_FRAC_N)
seqs <- seqs[frac_n <= MAX_FRAC_N]
cat("Sequences with more than", MAX_FRAC_N, "N, dropped:", n_many_n, "\n")

# Motifs: the PFMs as log2 odds PWMs against a uniform background,
# named by their JASPAR ID (the names are not unique across versions)
pfms <- readJASPARMatrix(JASPAR_FILE, matrixClass = "PFM")
pwms <- toPWM(pfms, pseudocounts = PWM_PSEUDOCOUNT)
names(pwms) <- ID(pwms)
cat("Motifs:", len(pwms), "from", JASPAR_FILE, "\n")

bpparam <- if (opts$cores > 1) MulticoreParam(opts$cores) else SerialParam()
cat("Cores =", opts$cores, "\n")

#' Motif enrichment of the up and down sets of one comparison
#'
#' The sets with at least --min-peaks sequences and the background are
#' bins of one calcBinnedMotifEnrR call with the background as its zero
#' bin, so the sequences are scanned once and each set is tested
#' against the background alone. The call also tests the background
#' against itself; that column is dropped, and padj is computed here
#' within each set.
#'
#' @param comp One comparison of the input RDS (tag, label, group1,
#'   group2, table)
#' @return List with the set sizes and statuses, and the results
#'   table (one row per set and motif)
run_motifs <- function(comp) {
  cat("========================================================\n")
  cat(comp$label, "\n")
  t0 <- Sys.time()

  peak_sets <- comp$table |>
    filter(PeakNo %in% names(seqs)) |>
    mutate(Set = case_when(
      FDR < FDR_CUT & logFC > 0 ~ "up",
      FDR < FDR_CUT & logFC < 0 ~ "down",
      .default = "background"
    ))

  sets <- tibble(Direction = c("up", "down"),
                 HigherIn = c(comp$group2, comp$group1)) |>
    mutate(Peaks = map_int(Direction, \(d) sum(peak_sets$Set == d)),
           Background = sum(peak_sets$Set == "background"),
           Tested = Peaks >= opts$min_peaks,
           Status = if_else(Tested, "tested",
                            str_c("not tested: fewer than ", opts$min_peaks,
                                  " peaks")))
  walk(seq_len(nrow(sets)), \(i) {
    cat("  ", sets$Direction[i], "(higher in", str_c(sets$HigherIn[i], "):"),
        sets$Peaks[i], "peaks,", sets$Status[i], "\n")
  })
  cat("   background:", sets$Background[1], "peaks\n")

  tested <- sets |> filter(Tested)
  if (nrow(tested) == 0) {
    return(list(tag = comp$tag, label = comp$label, group1 = comp$group1,
                group2 = comp$group2, sets = sets, results = NULL))
  }

  peak_sets <- peak_sets |>
    filter(Set %in% c(tested$Direction, "background"))
  bins <- factor(peak_sets$Set, levels = c(tested$Direction, "background")) |>
    setZeroBin("background")
  se <- calcBinnedMotifEnrR(seqs[as.character(peak_sets$PeakNo)], bins, pwms,
                            background = "zeroBin", test = "fisher",
                            maxFracN = MAX_FRAC_N, maxKmerSize = MAX_KMER,
                            min.score = MIN_SCORE, BPPARAM = bpparam)

  motifs <- rowData(se) |>
    as_tibble() |>
    dplyr::select(MotifID = motif.id, Motif = motif.name,
                  MotifGC = motif.percentGC)
  results <- map(tested$Direction, \(d) {
    col_data <- colData(se)[d, ]
    motifs |>
      mutate(Direction = d,
             log2enr = assay(se, "log2enr")[, d],
             pval = 10^-assay(se, "negLog10P")[, d],
             padj = p.adjust(pval, method = "BH"),
             FgWithHit = 100 * assay(se, "sumForegroundWgtWithHits")[, d] /
               col_data$totalWgtForeground,
             BgWithHit = 100 * assay(se, "sumBackgroundWgtWithHits")[, d] /
               col_data$totalWgtBackground)
  }) |>
    bind_rows() |>
    left_join(sets |> dplyr::select(Direction, HigherIn), by = "Direction") |>
    mutate(Direction = factor(Direction, levels = c("up", "down"))) |>
    arrange(Direction, pval) |>
    mutate(Direction = as.character(Direction)) |>
    dplyr::select(Direction, HigherIn, MotifID, Motif, log2enr, pval, padj,
                  FgWithHit, BgWithHit, MotifGC)

  sets <- sets |>
    left_join(results |>
                group_by(Direction) |>
                summarize(MotifsSig = sum(padj < FDR_CUT, na.rm = TRUE)),
              by = "Direction")
  walk(seq_len(nrow(tested)), \(i) {
    d <- tested$Direction[i]
    cat("  ", d, ":", sets$MotifsSig[sets$Direction == d], "motifs with padj <",
        FDR_CUT, "\n")
  })
  cat("   time", format(round(Sys.time() - t0, 1)), "\n")

  list(tag = comp$tag, label = comp$label, group1 = comp$group1,
       group2 = comp$group2, sets = sets, results = results)
}

motif_list <- map(diff_results$comparisons, run_motifs)
cat("========================================================\n\n")

# Summary: per comparison and direction, the set, its background and
# the motifs with padj < FDR_CUT
summary_stats <- motif_list |>
  map(\(m) {
    m$sets |>
      mutate(Comparison = m$tag, .before = 1) |>
      mutate(MotifsTested = if_else(Tested, len(pwms), 0L),
             MotifsSig = if (has_name(m$sets, "MotifsSig")) coalesce(MotifsSig, 0L) else 0L) |>
      dplyr::select(Comparison, Direction, HigherIn, Peaks, Background,
                    Status, MotifsTested, MotifsSig)
  }) |>
  bind_rows() |>
  dplyr::rename(!!str_glue("MotifsPadj<{FDR_CUT}") := MotifsSig)

parameters <- tribble(
  ~Parameter, ~Value,
  "Genome", GENOME,
  "Input", as.character(input_rds),
  "FASTA", fasta_file,
  "Motifs", str_glue("{JASPAR_NAME}, {len(pwms)} profiles ({path_file(JASPAR_FILE)})"),
  "PWM", str_glue("TFBSTools toPWM, log2 odds against a uniform background, pseudocount {PWM_PSEUDOCOUNT}"),
  "Region", str_glue("{opts$width} bp centered on the peak's strongest summit (largest MACS2 -log10 p of any sample's summits in the peak)"),
  "Peaks without a summit", str_glue("{n_no_summit}, centered on the middle of the peak"),
  "Sequences dropped", str_glue("{n_many_n}, more than {MAX_FRAC_N} N"),
  "Up set", str_glue("FDR < {FDR_CUT} and logFC > 0: more accessible in Group2"),
  "Down set", str_glue("FDR < {FDR_CUT} and logFC < 0: more accessible in Group1"),
  "Background", str_glue("the other tested peaks of the comparison, FDR >= {FDR_CUT}"),
  "Min peaks", as.character(opts$min_peaks),
  "Method", "monaLisa calcBinnedMotifEnrR, background = zeroBin (the background bin)",
  "Composition correction", str_glue("background weighted to the set's G+C content, then to its k-mers up to {MAX_KMER}"),
  "Motif hit", str_glue("score >= {MIN_SCORE} (log2 odds), either strand (monaLisa findMotifHits, matchPWM); a sequence counts once per motif"),
  "Test", "Fisher's exact test, one-sided (more set sequences with a hit than expected)",
  "padj", "Benjamini-Hochberg within each set (comparison and direction)",
  "log2enr", "log2 enrichment of sequences with a hit, set over weighted background (monaLisa, pseudocount 8)",
  "FgWithHit, BgWithHit", "percent of set and of (weighted) background sequences with a hit",
  "MotifGC", "percent G+C of the motif",
  "Comparison sheets", str_glue("every tested motif, NOT filtered for padj < {FDR_CUT}"),
  "monaLisa", as.character(packageVersion("monaLisa")),
  "TFBSTools", as.character(packageVersion("TFBSTools")),
  "Rsamtools", as.character(packageVersion("Rsamtools"))
)

# Excel: Summary (150% zoom, as in diffAnalysisPairwise.R), one sheet
# per comparison with motifs tested, Parameters
with_results <- motif_list |> keep(\(m) !is.null(m$results))
result_sheets <- map(with_results, "results")
names(result_sheets) <- map_chr(with_results, "tag") |> substr(1, 31)
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

saveRDS(list(parameters = parameters, regions = regions, motifs = motif_list),
        output_rds)

#' Heatmap of log2 enrichment, cells with padj < FDR_CUT
#'
#' Rows are motifs, clustered on log2enr (0 where not significant);
#' columns are sets. A cell is colored by log2enr where the motif has
#' padj < FDR_CUT in that set (PAL_EXPR, clipped to the 5th and 95th
#' percentiles), grey otherwise.
#'
#' @param res Tibble of results with Column (factor), MotifID, Motif,
#'   log2enr and padj, for the motifs drawn
#' @param title,subtitle Character
#' @return ggplot
plot_motif_heatmap <- function(res, title, subtitle) {
  cells <- res |>
    mutate(Label = str_glue("{Motif} ({MotifID})"),
           log2enr_sig = if_else(!is.na(padj) & padj < FDR_CUT, log2enr,
                                 NA_real_)) |>
    dplyr::select(Label, Column, log2enr_sig) |>
    complete(Label, Column)

  wide <- cells |>
    mutate(log2enr_sig = coalesce(log2enr_sig, 0)) |>
    pivot_wider(names_from = Column, values_from = log2enr_sig)
  mat <- as.matrix(wide[, -1])
  rownames(mat) <- wide$Label
  row_order <- if (nrow(mat) > 1) {
    rownames(mat)[hclust(dist(mat))$order]
  } else {
    rownames(mat)
  }

  ggplot(cells, aes(Column, factor(Label, levels = row_order),
                    fill = log2enr_sig)) +
    theme_minimal(base_size = 9) +
    geom_tile(color = TILE_BORDER) +
    scale_sequential(cells$log2enr_sig, PAL_EXPR, name = "log2\nenrichment") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank()) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL)
}

#' Heatmap of the top motifs of one comparison
#'
#' The N_TOP motifs with the smallest p-value among those with padj <
#' FDR_CUT of each tested set; one column per set.
#'
#' @param m One element of motif_list
#' @return ggplot, or NULL if no motif has padj < FDR_CUT
plot_comparison <- function(m) {
  if (is.null(m$results)) {
    return(NULL)
  }
  top <- m$results |>
    filter(!is.na(padj), padj < FDR_CUT) |>
    group_by(Direction) |>
    slice_min(pval, n = N_TOP, with_ties = FALSE) |>
    ungroup()
  if (nrow(top) == 0) {
    return(NULL)
  }
  tested <- m$sets |> filter(Tested)
  columns <- str_glue("Higher in {tested$HigherIn}\n({tested$Direction}, ",
                      "{tested$Peaks} peaks)") |>
    set_names(tested$Direction)
  res <- m$results |>
    filter(MotifID %in% top$MotifID) |>
    mutate(Column = factor(columns[Direction], levels = columns))
  plot_motif_heatmap(
    res,
    title = str_glue("{m$label}: motif enrichment"),
    subtitle = str_glue(
      "{counts_subtitle(m)}; background: {tested$Background[1]} peaks ",
      "with FDR >= {FDR_CUT}\n",
      "At most {N_TOP} motifs per set with the smallest p-value drawn.\n",
      "Grey: padj >= {FDR_CUT}. Color clipped to the 5th-95th percentile.")
  )
}

#' Subtitle with the motifs significant in each set of one comparison
#'
#' @param m One element of motif_list
#' @return Character
counts_subtitle <- function(m) {
  m$sets |>
    mutate(Text = if_else(
      Tested,
      as.character(str_glue("{Direction}: {MotifsSig} of {len(pwms)} motifs ",
                            "with padj < {FDR_CUT}")),
      as.character(str_glue("{Direction}: {Peaks} peaks, not tested")))) |>
    pull(Text) |>
    str_c(collapse = "; ")
}

#' Heatmap of the top motifs across every tested set
#'
#' Rows: the motifs with padj < FDR_CUT in at least one set, at most
#' N_HEATMAP with the smallest padj.
#'
#' @return ggplot, or NULL if no motif is significant in any set
plot_all <- function() {
  all_res <- with_results |>
    map(\(m) m$results |> mutate(Column = str_glue("{m$label}: {Direction}"))) |>
    bind_rows()
  candidates <- all_res |>
    filter(!is.na(padj), padj < FDR_CUT) |>
    group_by(MotifID) |>
    summarize(minPadj = min(padj), .groups = "drop")
  if (nrow(candidates) == 0) {
    return(NULL)
  }
  rows <- candidates |>
    arrange(minPadj, MotifID) |>
    head(N_HEATMAP)
  columns <- unique(all_res$Column)
  res <- all_res |>
    filter(MotifID %in% rows$MotifID) |>
    mutate(Column = factor(Column, levels = columns))
  shown <- if (nrow(candidates) > N_HEATMAP) {
    str_glue(", the {N_HEATMAP} with the smallest padj shown")
  } else {
    ""
  }
  plot_motif_heatmap(
    res,
    title = "Motif enrichment across comparisons",
    subtitle = str_glue(
      "{nrow(candidates)} motifs with padj < {FDR_CUT} in at least one ",
      "set{shown}.\nup: more accessible in the first group of the ",
      "comparison, down: in the second.\nGrey: padj >= {FDR_CUT} or not ",
      "tested. Color clipped to the 5th-95th percentile.")
  )
}

pdf(output_pdf, width = 8.5, height = 11)
for (m in motif_list) {
  p <- plot_comparison(m)
  if (!is.null(p)) {
    print(p)
  }
}
n_sets <- sum(map_int(with_results, \(m) sum(m$sets$Tested)))
if (n_sets > 1) {
  p <- plot_all()
  if (!is.null(p)) {
    print(p)
  }
}
dev.off()

cat("Wrote", output_xlsx, output_pdf, output_rds, "\n")
cat("## END: motifEnrich.R\n")
