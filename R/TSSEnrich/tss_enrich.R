#!/usr/bin/env Rscript
# TSS Enrichment analysis - R port of ENCODE encode_task_tss_enrich.py
# Outputs: {prefix}.tss_enrich.qc, {prefix}.tss_enrich.png, {prefix}.large_tss_enrich.png
#
# IMPORTANT - BAM INPUT REQUIREMENT:
#   Duplicate reads MUST be physically REMOVED from the BAM before running
#   this script. A BAM where duplicates have only been MARKED (e.g., the
#   output of Picard MarkDuplicates with default settings, or samtools markdup)
#   but not removed will produce INFLATED, INCORRECT scores because all
#   duplicate reads will be counted as signal.
#
#   Duplicate reads are reads arising from PCR amplification of the same
#   original DNA fragment. They are not independent observations of chromatin
#   accessibility and must be excluded.
#
#   To produce a correct input BAM:
#     1. Run MarkDuplicates (Picard/GATK) or equivalent to FLAG duplicates
#        (SAM FLAG bit 0x400 = 1024).
#     2. Then filter them OUT:
#          samtools view -F 1024 -b marked.bam > nodup.bam
#          samtools index nodup.bam
#
#   The "--nodup-bam" argument name reflects this requirement: the BAM must
#   have duplicates removed (not-duplicated = nodup), not merely marked.

suppressPackageStartupMessages({
    library(optparse)
    library(GenomicRanges)
    library(Rsamtools)
    library(ggplot2)
    library(dplyr)
})

# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------
option_list <- list(
    make_option("--nodup-bam",    type="character", help="BAM file with duplicate reads PHYSICALLY REMOVED (not merely marked). Use samtools view -F 1024 after MarkDuplicates to produce this file."),
    make_option("--chrsz",        type="character", help="2-col chromosome sizes file"),
    make_option("--tss",          type="character", help="TSS BED file (BED6, 1-bp regions)"),
    make_option("--out-dir",      type="character", default=NULL, help="Output directory (default: out/tssEnrich/<sampleId> from BAM @RG SM tag)"),
    make_option("--read-len",     type="integer",   help="Read length (integer)"),
    make_option("--read-len-log", type="character", help="File containing read length")
)
opts <- parse_args(OptionParser(option_list=option_list))

bam_file     <- opts[["nodup-bam"]]
chrsz        <- opts[["chrsz"]]
tss_file     <- opts[["tss"]]
out_dir      <- opts[["out-dir"]]
read_len     <- opts[["read-len"]]
read_len_log <- opts[["read-len-log"]]

if (!is.null(read_len_log))
    read_len <- as.integer(trimws(readLines(read_len_log, n=1)))
if (is.null(read_len))
    stop("Must supply --read-len or --read-len-log")

if (is.null(out_dir)) {
    hdr     <- system2("samtools", c("view", "-H", bam_file), stdout=TRUE)
    rg_lines <- grep("^@RG", hdr, value=TRUE)
    sm_vals  <- regmatches(rg_lines, regexpr("(?<=SM:)[^\t]+", rg_lines, perl=TRUE))
    sm_vals  <- unique(sm_vals)
    if (length(sm_vals) == 0) stop("No SM tag found in BAM @RG header; use --out-dir")
    if (length(sm_vals) > 1) warning("Multiple SM values found; using first: ", sm_vals[1])
    sampleId <- sm_vals[1]
    out_dir  <- file.path("out", "tssEnrich", sampleId)
    message("--out-dir not specified; using: ", out_dir)
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

prefix <- file.path(out_dir, sub("\\.bam$", "", basename(bam_file)))
message("Output prefix: ", prefix)

# ---------------------------------------------------------------------------
# Parameters (matching Python defaults)
# ---------------------------------------------------------------------------
BP_EDGE       <- 2000L          # bp on each side of TSS
BINS          <- 400L           # number of output bins
WIN           <- 2L*BP_EDGE+1L  # window size = 4001 bp
SHIFT         <- -read_len/2    # = -50.5 for read_len=101
NUM_EDGE_BINS <- 10L            # int(100 / (2*bp_edge/bins))

# metaseq local_coverage padded window:
#   pad_start = win_start - shift_width  (= win_start + 50.5 for shift=-50.5)
#   pad_stop  = win_stop  + shift_width  (= win_stop  - 50.5)
# These are truncated (toward zero) in Python's C int assignment.
PAD_START_OFFSET <- trunc(-SHIFT)    # +50
PAD_STOP_OFFSET  <- trunc(SHIFT)     # -51 (trunc(-50.5) = -50... wait)
# trunc(-50.5) in R = -50 (toward zero)
# so pad_stop = win_stop + trunc(SHIFT) = win_stop + trunc(-50.5) = win_stop - 50
# Let me compute directly:
# Python: pad_stop = win_stop + shift_width = win_stop + (-50.5)
#         stored as C int: trunc(win_stop - 50.5) = win_stop - 51 for positive values
#         Actually: (long)(1000.0 - 50.5) = (long)(949.5) = 949  [truncation toward zero]
#         So: trunc(win_stop - 50.5) = win_stop - 51 when (win_stop - 50.5) is not integer
#         win_stop is always integer, so win_stop - 50.5 is always .5 fractional
#         trunc(win_stop - 50.5) = win_stop - 51  [because floor of x.5 = x for positive]
# Hmm: trunc(949.5) = 949 = 1000 - 51? No: 1000 - 50.5 = 949.5, trunc = 949 = 1000 - 51. Yes!
PAD_DELTA <- trunc(abs(SHIFT))  # = 50 (trunc of 50.5)
# pad_start = win_start + PAD_DELTA     (win_start - SHIFT = win_start + 50.5, trunc = win_start+50)
# pad_stop  = win_stop  - PAD_DELTA - 1 (win_stop + SHIFT = win_stop - 50.5, trunc = win_stop-51)

# Interpolation grid: np.linspace(0, WIN-1, BINS) equivalent
x_orig   <- seq_len(WIN) - 1L       # 0, 1, ..., 4000
x_interp <- seq(0, WIN-1, length.out=BINS)

# ---------------------------------------------------------------------------
# Load chromosome sizes and TSS BED
# ---------------------------------------------------------------------------
message("Loading TSS file...")
chrom_sizes <- read.table(chrsz, header=FALSE, sep="\t",
                          col.names=c("chrom","size"), stringsAsFactors=FALSE)

tss_df <- read.table(tss_file, header=FALSE, sep="\t", stringsAsFactors=FALSE,
                     col.names=c("chrom","start","end","name","score","strand")) %>%
    left_join(chrom_sizes, by="chrom") %>%
    filter(!is.na(size)) %>%
    # slop: expand by BP_EDGE on each side, clip to chrom bounds
    mutate(start = pmax(0L, start - BP_EDGE),
           end   = pmin(size, end + BP_EDGE)) %>%
    # keep only full-width windows
    filter((end - start) == WIN) %>%
    select(chrom, start, end, name, score, strand)

message(sprintf("Using %d TSS windows", nrow(tss_df)))

# ---------------------------------------------------------------------------
# Ensure BAM is indexed
# ---------------------------------------------------------------------------
if (!file.exists(paste0(bam_file, ".bai"))) {
    message("Indexing BAM...")
    indexBam(bam_file)
}

# ---------------------------------------------------------------------------
# GUARD: Abort if duplicate reads are present in the BAM.
#
# Duplicate reads (PCR copies of the same original DNA fragment) must be
# PHYSICALLY REMOVED before running this analysis. If they are only MARKED
# (SAM FLAG bit 0x400 = 1024 set but reads still present in the file), they
# will be counted as signal and the TSS enrichment score will be inflated
# and incorrect.
#
# This check reads the idxstats index to get total mapped reads, then uses
# flagstat to count how many carry the duplicate flag. If any are found the
# script aborts with instructions for how to fix the input.
# ---------------------------------------------------------------------------
message("Checking BAM for duplicate reads...")
flagstat_out <- system2("samtools", c("flagstat", bam_file),
                        stdout=TRUE, stderr=FALSE)
dup_line <- grep("^[0-9]+ \\+ [0-9]+ duplicate", flagstat_out, value=TRUE)
if (length(dup_line) > 0L) {
    n_dups <- as.integer(sub("^([0-9]+) .*", "\\1", dup_line))
    if (n_dups > 0L) {
        stop(sprintf(
            paste0(
                "\n\n",
                "ERROR: BAM file contains %d reads with the duplicate flag set (SAM FLAG 0x400).\n",
                "\n",
                "This script requires that duplicate reads be PHYSICALLY REMOVED from the BAM,\n",
                "not merely marked. Marked-but-present duplicates are counted as signal and\n",
                "will produce an inflated, incorrect TSS enrichment score.\n",
                "\n",
                "To fix, filter out duplicate-flagged reads before running this script:\n",
                "  samtools view -F 1024 -b %s > nodup.bam\n",
                "  samtools index nodup.bam\n",
                "Then rerun with --nodup-bam nodup.bam\n"
            ),
            n_dups, bam_file
        ))
    }
}
message("BAM duplicate check passed (0 duplicate-flagged reads present).")

# ---------------------------------------------------------------------------
# Build signal matrix: N_TSS x BINS
#
# Exact port of metaseq _local_coverage() for BAM:
#   1. Fetch reads from padded window [win_start + PAD_DELTA, win_stop - PAD_DELTA - 1]
#      (inset, not outset, because shift_width is negative)
#   2. For each CIGAR M/=/X segment [seg_start, seg_end):
#      - + strand: shift by SHIFT; - strand: shift by -SHIFT
#      - Clip to window and accumulate into bp-resolution profile (length WIN)
#   3. Downsample 4001-bp profile to BINS via linear interpolation (np.interp)
#   4. Reverse profile for minus-strand TSS windows
# ---------------------------------------------------------------------------
message("Building signal matrix...")

# Helper: parse all M/=/X segments from a vector of CIGAR strings.
# Returns a data.frame with columns: read_idx, seg_start (0-based), seg_len
# Only handles ops that consume reference: M, D, N, =, X (I/S/H/P do not).
parse_cigars <- function(pos0, cigars) {
    # pos0: 0-based leftmost positions (integer vector, length n_reads)
    # cigars: character vector of CIGAR strings
    # Returns list of (seg_start, seg_len, read_idx) for M/=/X ops only
    n <- length(cigars)
    all_starts <- integer(0)
    all_lens   <- integer(0)
    all_ridx   <- integer(0)

    for (j in seq_len(n)) {
        cig <- cigars[j]
        # Fast path: pure NM (e.g. "101M") — most common ATAC case
        m <- regmatches(cig, regexpr("^([0-9]+)M$", cig))
        if (length(m) == 1L && nchar(m) > 0L) {
            all_starts <- c(all_starts, pos0[j])
            all_lens   <- c(all_lens,   as.integer(sub("M","",m)))
            all_ridx   <- c(all_ridx,   j)
            next
        }
        # General CIGAR parsing
        lens_raw <- as.integer(regmatches(cig, gregexpr("[0-9]+", cig))[[1]])
        ops_raw  <- regmatches(cig, gregexpr("[A-Z]",  cig))[[1]]
        cur <- pos0[j]
        for (k in seq_along(ops_raw)) {
            op  <- ops_raw[k]
            len <- lens_raw[k]
            if (op == "M" || op == "=" || op == "X") {
                all_starts <- c(all_starts, cur)
                all_lens   <- c(all_lens,   len)
                all_ridx   <- c(all_ridx,   j)
            }
            if (op == "M" || op == "D" || op == "N" || op == "=" || op == "X")
                cur <- cur + len
        }
    }
    list(start=all_starts, len=all_lens, ridx=all_ridx)
}

# Vectorized profile accumulation for a set of segments.
# Uses the diff-of-cumsum (start/stop sentinel) trick — fully vectorized, no R loop.
accumulate_segments <- function(seg_starts, seg_lens, ws, WIN) {
    profile <- numeric(WIN)
    if (length(seg_starts) == 0L) return(profile)
    si <- pmax(0L,  as.integer(seg_starts)            - ws)
    ei <- pmin(WIN, as.integer(seg_starts + seg_lens) - ws)
    keep <- si < WIN & ei > 0L & si < ei
    si <- si[keep] + 1L   # convert to 1-based R index for starts
    ei <- ei[keep] + 1L   # one past end (also 1-based)
    if (length(si) == 0L) return(profile)
    delta <- numeric(WIN + 1L)
    # tabulate replaces the loop: adds 1 at each start, subtracts 1 at each end
    delta <- delta + tabulate(si,   nbins=WIN+1L)
    delta <- delta - tabulate(ei,   nbins=WIN+1L)
    cumsum(delta)[seq_len(WIN)]
}

CHUNK_SIZE <- 500L
n_tss      <- nrow(tss_df)
bam_mat    <- matrix(0.0, nrow=n_tss, ncol=BINS)

# Pre-compute padded fetch coordinates for all windows
pad_starts <- pmax(0L, tss_df$start + PAD_DELTA)
pad_stops  <- tss_df$end - PAD_DELTA - 1L
valid_win  <- pad_stops > pad_starts

for (ci in seq_len(ceiling(n_tss / CHUNK_SIZE))) {
    idx <- ((ci-1L)*CHUNK_SIZE + 1L):min(ci*CHUNK_SIZE, n_tss)

    # Batch fetch: one scanBam call for all windows in chunk
    chunk_valid <- idx[valid_win[idx]]
    if (length(chunk_valid) == 0L) next

    chunk_gr <- GRanges(
        seqnames = tss_df$chrom[chunk_valid],
        ranges   = IRanges(pad_starts[chunk_valid]+1L, pad_stops[chunk_valid])
    )
    param    <- ScanBamParam(which=chunk_gr, what=c("pos","cigar","flag"))
    bam_data <- scanBam(bam_file, param=param)

    for (i in seq_along(chunk_valid)) {
        ti   <- chunk_valid[i]
        ws   <- tss_df$start[ti]
        wstr <- tss_df$strand[ti]

        rd <- bam_data[[i]]
        if (length(rd$pos) == 0L) next

        pos0     <- rd$pos - 1L
        is_minus <- bitwAnd(rd$flag, 16L) != 0L

        segs_p <- parse_cigars(pos0[!is_minus], rd$cigar[!is_minus])
        segs_m <- parse_cigars(pos0[ is_minus], rd$cigar[ is_minus])

        prof_p <- accumulate_segments(trunc(segs_p$start + SHIFT), segs_p$len, ws, WIN)
        prof_m <- accumulate_segments(trunc(segs_m$start - SHIFT), segs_m$len, ws, WIN)

        profile <- prof_p + prof_m

        binned <- approx(x=x_orig, y=profile, xout=x_interp,
                         method="linear", rule=2)$y
        if (wstr == "-") binned <- rev(binned)
        bam_mat[ti, ] <- binned
    }

    if (ci %% 10L == 0L)
        message(sprintf("  chunk %d / %d", ci, ceiling(n_tss/CHUNK_SIZE)))
}

message("Matrix built.")

# ---------------------------------------------------------------------------
# Greenleaf normalization
# ---------------------------------------------------------------------------
bin_means <- colMeans(bam_mat)
avg_noise <- mean(c(bin_means[1:NUM_EDGE_BINS],
                    bin_means[(BINS-NUM_EDGE_BINS+1L):BINS]))
bam_mat   <- bam_mat / avg_noise
bin_means <- colMeans(bam_mat)
tss_score <- max(bin_means)

message(sprintf("TSS Enrichment Score: %.6f", tss_score))

# ---------------------------------------------------------------------------
# Write QC file
# ---------------------------------------------------------------------------
qc_file <- paste0(prefix, ".tss_enrich.qc")
writeLines(as.character(tss_score), qc_file)
message("Written: ", qc_file)

# ---------------------------------------------------------------------------
# Simple line plot  (matching Python's simple tss_enrich.png)
# ---------------------------------------------------------------------------
x_pos   <- seq(-BP_EDGE, BP_EDGE, length.out=BINS)
plot_df <- data.frame(pos=x_pos, signal=bin_means)

p <- ggplot(plot_df, aes(x=pos, y=signal)) +
    geom_line(color="red") +
    geom_vline(xintercept=0, linetype="dotted", color="black") +
    labs(x="Distance from TSS (bp)", y="TSS Enrichment") +
    theme_classic()

tss_plot_file <- paste0(prefix, ".tss_enrich.png")
ggsave(tss_plot_file, plot=p, width=6, height=4, dpi=150)
message("Written: ", tss_plot_file)

# ---------------------------------------------------------------------------
# Large heatmap plot  (matching Python's large_tss_enrich.png)
# Rows sorted ascending by row mean (weakest at top), aggregate line overlaid
# ---------------------------------------------------------------------------
message("Generating large TSS plot...")

row_means  <- rowMeans(bam_mat)
sort_order <- order(row_means)          # ascending
mat_sorted <- bam_mat[sort_order, ]

all_vals   <- as.vector(mat_sorted)
vmin_val   <- quantile(all_vals, 0.05)
vmax_pct   <- 0.99
if (quantile(all_vals, vmax_pct) == 0.0) vmax_pct <- 1.0
vmax_val   <- quantile(all_vals, vmax_pct)

mat_clipped <- pmin(pmax(mat_sorted, vmin_val), vmax_val)
n_rows      <- nrow(mat_clipped)

hm_df <- data.frame(
    bin = rep(x_pos, each=n_rows),
    row = rep(seq_len(n_rows), times=BINS),
    val = as.vector(mat_clipped)
)

line_df <- data.frame(pos=x_pos, signal=bin_means)
scale_y  <- n_rows / max(bin_means, 1)

p_large <- ggplot() +
    geom_raster(data=hm_df, aes(x=bin, y=row, fill=val)) +
    scale_fill_gradient(low="white", high="black", name="Signal") +
    geom_line(data=line_df, aes(x=pos, y=signal*scale_y),
              color="black", linewidth=0.5) +
    geom_vline(xintercept=0, linetype="dotted", color="grey50") +
    labs(x="Distance from TSS (bp)", y="TSS (sorted by strength)") +
    theme_classic() +
    theme(legend.position="right")

large_plot_file <- paste0(prefix, ".large_tss_enrich.png")
ggsave(large_plot_file, plot=p_large, width=5, height=10, dpi=150)
message("Written: ", large_plot_file)

message("All done.")
