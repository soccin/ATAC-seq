args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
    cat("
================================================================================
  TSS Enrichment Synthetic Data Test (Ground Truth Validation)
================================================================================

USAGE:
  Rscript syntheticTSS.R GENOME [OUT_PREFIX]

ARGUMENTS:
  GENOME      Genome: hg19, b38, or mouse
  OUT_PREFIX  Output file prefix (default: syntheticTSS)

OUTPUT:
  <OUT_PREFIX>_profile.pdf    Profile plot for synthetic data
  <OUT_PREFIX>_scores.txt     Score comparison table

TEST DESIGN:
  - Take 500 + strand TSSs from chr1
  - Place SIGNAL_READS reads per TSS at exactly the TSS position (+/-JITTER bp)
  - Place BACKGROUND_READS read per TSS uniformly across the +-HALFWIN bp window
  - Expected enrichment: SIGNAL_READS / BACKGROUND_READS * (window / signal_window)
  - Both methods should show a sharp peak at position 0
  - This is the only test with a known ground-truth expected score
================================================================================
")
    quit()
}

GENOME     <- args[1]
OUT_PREFIX <- if (length(args) >= 2) args[2] else "syntheticTSS"

# ---- Parameters ----------------------------------------------------------------

HALFWIN          <- 2000L
BINSIZE          <- 10L
BG_BP            <- 200L
SCORE_BP         <- 200L

N_TSS            <- 500L    # number of TSSs to use from chr1
SIGNAL_READS     <- 50L     # reads per TSS placed at the TSS
JITTER           <- 5L      # place signal reads uniformly within +-JITTER bp of TSS
BACKGROUND_READS <- 1L      # reads per TSS placed uniformly across the window
READ_LENGTH      <- 50L     # synthetic read length

# ---- Dependencies --------------------------------------------------------------

suppressPackageStartupMessages({
    library(rtracklayer)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(ggplot2)
    library(ATACseqQC)
})

# ---- TxDb selection ------------------------------------------------------------

getTxDb <- function(genome) {
    if (genome == "hg19") {
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    } else if (genome == "b38") {
        TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    } else if (genome == "mouse") {
        TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    } else {
        stop("Unknown genome: ", genome, ". Use hg19, b38, or mouse.")
    }
}

# ---- TSS windows (from computeTSSEnrichment.R) ---------------------------------

getTSSWindows <- function(txdb, halfwin=2000L, binsize=10L) {
    tss <- unique(resize(transcripts(txdb), width=1L, fix="start"))
    std_chrs <- paste0("chr", c(1:22, "X", "Y"))
    tss <- keepSeqlevels(tss, intersect(std_chrs, seqlevels(tss)), pruning.mode="coarse")
    winsize <- 2L * halfwin + binsize
    win <- trim(resize(tss, width=winsize, fix="center"))
    win[width(win) == winsize]
}

# ---- Coverage profile (from computeTSSEnrichment.R) ----------------------------

computeProfile <- function(reads_gr, windows, halfwin=2000L, binsize=10L) {
    cov <- coverage(reads_gr)

    winsize <- 2L * halfwin + binsize
    nbins   <- winsize %/% binsize

    bin_sum <- numeric(nbins)
    n_total <- 0L

    chroms <- intersect(names(cov), unique(as.character(seqnames(windows))))

    for (chr in chroms) {
        cv <- cov[[chr]]
        w  <- windows[seqnames(windows) == chr]
        chr_len <- length(cv)
        w <- w[end(w) <= chr_len]

        for (s in c("+", "-")) {
            ws <- w[strand(w) == s]
            if (length(ws) == 0L) next

            v       <- Views(cv, start=start(ws), end=end(ws))
            pos_sum <- as.numeric(Reduce("+", as(v, "IntegerList")))

            if (length(pos_sum) != winsize) next

            if (s == "-") pos_sum <- rev(pos_sum)

            bin_sum <- bin_sum + colSums(matrix(pos_sum, nrow=binsize))
            n_total <- n_total + length(ws)
        }
    }

    profile   <- bin_sum / n_total
    positions <- seq(-halfwin, halfwin, by=binsize)  # 401 bins, TSS at center (0)

    list(profile=profile, positions=positions, n=n_total)
}

# ---- Normalize and score -------------------------------------------------------

normProfile <- function(res, bg_bp=200L, score_bp=200L) {
    binsize <- res$positions[2] - res$positions[1]
    bg_n    <- round(bg_bp / binsize)
    nbins   <- length(res$profile)

    bg   <- mean(c(res$profile[seq_len(bg_n)], res$profile[(nbins - bg_n + 1L):nbins]))
    norm <- res$profile / bg

    score <- max(norm[abs(res$positions) <= score_bp])
    list(norm=norm, positions=res$positions, score=score, n=res$n)
}

# ---- Load TSSs -----------------------------------------------------------------

cat("Genome:", GENOME, "\n")
txdb    <- getTxDb(GENOME)
windows <- getTSSWindows(txdb, HALFWIN, BINSIZE)

# Use N_TSS TSSs from chr1 (both strands).
# Including both strands is required: if only + strand TSSs have signal reads,
# nearby - strand TSS windows see those reads reversed at negative positions,
# creating a spurious peak upstream of the TSS.
tss_all  <- unique(resize(transcripts(txdb), width=1L, fix="start"))
tss_chr1 <- tss_all[seqnames(tss_all) == "chr1"]
tss_chr1 <- sort(tss_chr1)

# start(tss) is the 5' end (TSS position) for both strands after resize(fix="start").
# Require enough room for a full window on both sides.
tss_valid <- tss_chr1[start(tss_chr1) > HALFWIN + READ_LENGTH]

if (length(tss_valid) < N_TSS) {
    cat(sprintf("WARNING: only %d valid chr1 TSSs available (requested %d). Using all.\n",
        length(tss_valid), N_TSS))
    N_TSS <- length(tss_valid)
}

tss_use <- tss_valid[seq_len(N_TSS)]
cat(sprintf("Using %d TSSs from chr1 (%d + strand, %d - strand)\n",
    N_TSS,
    sum(strand(tss_use) == "+"),
    sum(strand(tss_use) == "-")))

# ---- Compute expected score analytically ---------------------------------------
# Reads are centered at the TSS (start = TSS - READ_LENGTH/2 + offset).
# All SIGNAL_READS reads cover the TSS bin (10 bp centered at -5 or +5 bp).
# The bin with maximum coverage gets contributions from all SIGNAL_READS reads,
# each contributing READ_LENGTH / binsize positions, but the bin is only 10 bp wide.
# For reads centered at the TSS, peak bin coverage per TSS ~ SIGNAL_READS reads.
# Background: BACKGROUND_READS reads per TSS spread uniformly over 2*HALFWIN bp.
#   Expected coverage per 10 bp bin = BACKGROUND_READS * READ_LENGTH / (2*HALFWIN/binsize)
#                                   = BACKGROUND_READS * READ_LENGTH * BINSIZE / (2*HALFWIN)
# Expected enrichment = peak_bin_coverage / background_bin_coverage
bg_per_bin    <- BACKGROUND_READS * READ_LENGTH * BINSIZE / (2L * HALFWIN + BINSIZE)
peak_per_bin  <- SIGNAL_READS     # all reads cover the TSS bin
expected_score <- peak_per_bin / bg_per_bin
cat(sprintf("Expected enrichment: %.1fx\n", expected_score))
cat(sprintf("  (signal: %d reads/TSS at TSS, binsize: %d bp, background: %d read/TSS over +-%d bp)\n",
    SIGNAL_READS, BINSIZE, BACKGROUND_READS, HALFWIN))

# ---- Generate synthetic reads --------------------------------------------------

set.seed(42L)

tss_pos <- start(tss_use)   # 5' end (TSS) for both strands after resize(fix="start")
chr_str <- as.character(seqnames(tss_use))

# Signal reads: centered at each TSS within +-JITTER bp.
# Centering (not starting at TSS) is critical: reads that start at the TSS
# extend only downstream, creating a peak 25 bp downstream, not at the TSS.
# start = TSS - READ_LENGTH/2 + offset places the read midpoint at TSS + offset.
center_offset   <- READ_LENGTH %/% 2L
signal_offsets  <- sample(-JITTER:JITTER, N_TSS * SIGNAL_READS, replace=TRUE)
signal_starts   <- rep(tss_pos, each=SIGNAL_READS) - center_offset + signal_offsets
signal_gr <- GRanges(
    seqnames = rep(chr_str, each=SIGNAL_READS),
    ranges   = IRanges(start=pmax(1L, signal_starts), width=READ_LENGTH),
    strand   = "*"
)

# Background reads: uniformly across +-HALFWIN window
bg_offsets  <- sample(-HALFWIN:(HALFWIN - READ_LENGTH), N_TSS * BACKGROUND_READS, replace=TRUE)
bg_starts   <- rep(tss_pos, each=BACKGROUND_READS) + bg_offsets
bg_gr <- GRanges(
    seqnames = rep(chr_str, each=BACKGROUND_READS),
    ranges   = IRanges(start=pmax(1L, bg_starts), width=READ_LENGTH),
    strand   = "*"
)

synthetic_gr <- c(signal_gr, bg_gr)
synthetic_gr <- sort(synthetic_gr)

cat(sprintf("Synthetic reads: %d signal + %d background = %d total\n",
    length(signal_gr), length(bg_gr), length(synthetic_gr)))

# ---- Run BED method on synthetic reads -----------------------------------------

cat("\nRunning BED method on synthetic data...\n")
res_syn    <- computeProfile(synthetic_gr, windows, HALFWIN, BINSIZE)
normed_syn <- normProfile(res_syn, BG_BP, SCORE_BP)

cat(sprintf("BED method score:  %.4f  (expected: ~%.1f)\n", normed_syn$score, expected_score))
cat(sprintf("Peak position:     %.0f bp (expected: ~-5 to +5 bp)\n",
    normed_syn$positions[which.max(normed_syn$norm)]))

# ---- Run ATACseqQC on synthetic data -------------------------------------------
# ATACseqQC requires a BAM. Write synthetic reads to a temporary BED and use
# the BED method only (ATACseqQC BAM-based run omitted here -- include separately
# if a BAM writer is available).
cat("\nNote: ATACseqQC requires a BAM file; synthetic BED-only run uses BED method.\n")
cat("To run ATACseqQC on synthetic data, write the GRanges to a BAM with rtracklayer\n")
cat("or Rsamtools and pass it to testTSS.R.\n")

# ---- Write synthetic BED for inspection / external tools -----------------------

bed_file <- paste0(OUT_PREFIX, "_synthetic.bed.gz")
export(synthetic_gr, bed_file, format="BED")
cat("Synthetic reads written to:", bed_file, "\n")

# ---- Plot ----------------------------------------------------------------------

df_syn <- data.frame(pos=normed_syn$positions, enrichment=normed_syn$norm)

peak_pos <- normed_syn$positions[which.max(normed_syn$norm)]

p_syn <- ggplot(df_syn, aes(pos, enrichment)) +
    geom_line(color="steelblue", linewidth=0.7) +
    geom_vline(xintercept=0, linetype="dashed", color="grey50", linewidth=0.5) +
    geom_vline(xintercept=peak_pos, linetype="dotted", color="tomato", linewidth=0.8) +
    annotate("text", x=Inf, y=Inf,
        label=sprintf(
            "Observed score: %.2f\nExpected score: ~%.1f\nPeak at: %.0f bp\n(n = %d TSSs)",
            normed_syn$score, expected_score, peak_pos, normed_syn$n),
        hjust=1.1, vjust=1.5, size=3.5) +
    labs(
        title="Synthetic TSS Test (Ground Truth)",
        subtitle=sprintf(
            "%d signal reads/TSS at +-5 bp; %d background read/TSS over +-%d bp",
            SIGNAL_READS, BACKGROUND_READS, HALFWIN),
        x="Distance to TSS (bp)", y="Normalized enrichment") +
    theme_light(base_size=14)

pdf_file <- paste0(OUT_PREFIX, "_profile.pdf")
pdf(pdf_file, width=8, height=5)
print(p_syn)
invisible(dev.off())
cat("Plot written to:", pdf_file, "\n")

# ---- Score table ---------------------------------------------------------------

score_file <- paste0(OUT_PREFIX, "_scores.txt")
lines <- c(
    paste("Test", "Score", "Expected", "PeakPos_bp", sep="\t"),
    paste("BED_method",
          sprintf("%.4f", normed_syn$score),
          sprintf("~%.1f",  expected_score),
          sprintf("%.0f",   peak_pos),
          sep="\t")
)
writeLines(lines, score_file)
cat("Scores written to:", score_file, "\n")
