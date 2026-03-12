args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    cat("
================================================================================
  ATAC-seq TSS Enrichment Score and Profile
================================================================================

USAGE:
  Rscript computeTSSEnrichment.R GENOME BED_FILE [SAMPLE_NAME] [OUT_PREFIX]

ARGUMENTS:
  GENOME       Genome: hg19, b38, or mouse
  BED_FILE     Tn5-shifted reads BED file (can be gzipped)
  SAMPLE_NAME  Label for plots (default: BED_FILE basename)
  OUT_PREFIX   Output file prefix (default: BED_FILE without extension)

OUTPUT:
  <OUT_PREFIX>_tss_enrich.csv  CSV with header: SampleName, Score, N_TSS
  <OUT_PREFIX>_TSSprofile.pdf  Profile plot

DETAILS:
  Computes the ENCODE-style TSS enrichment score:
    - Aggregate coverage in a +-2000bp window around each TSS
    - Normalize by mean coverage in the outermost 200bp flanks
    - Score = max normalized signal within +-200bp of TSS
  Center shift applied per-read: each read shifted by width(read) %/% 2 so that
  its midpoint lands at the Tn5 cut site (handles variable read lengths).
================================================================================
")
    quit()
}

GENOME      <- args[1]
BED_FILE    <- args[2]
SAMPLE_NAME <- if (length(args) >= 3) args[3] else basename(BED_FILE)
OUT_PREFIX  <- if (length(args) >= 4) args[4] else sub("\\.bed(\\.gz)?$", "", BED_FILE)

HALFWIN  <- 2000L   # bp each side of TSS
BINSIZE  <- 10L     # bp per bin for profile
BG_BP    <- 200L    # outermost Nbp used for background normalization
SCORE_BP <- 200L    # score = max signal within +-SCORE_BP of TSS

# ---- Dependencies ------------------------------------------------------------

suppressPackageStartupMessages({
    library(rtracklayer)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(ggplot2)
})

# ---- TxDb selection ----------------------------------------------------------

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

# ---- TSS windows -------------------------------------------------------------
# resize(transcripts, width=1, fix="start") is strand-aware:
#   + strand: anchors at leftmost coordinate (5' end = TSS)
#   - strand: anchors at rightmost coordinate (5' end = TSS)
# Then expand to a 2*HALFWIN bp window centered on the TSS.

getTSSWindows <- function(txdb, halfwin=2000L, binsize=10L) {
    tss <- unique(resize(transcripts(txdb), width=1L, fix="start"))
    std_chrs <- paste0("chr", c(1:22, "X", "Y"))
    tss <- keepSeqlevels(tss, intersect(std_chrs, seqlevels(tss)), pruning.mode="coarse")
    # Use 2*halfwin + binsize (4010 bp) so there are an odd number of bins (401)
    # with one bin centered exactly at the TSS (position 0).
    winsize <- 2L * halfwin + binsize
    win <- trim(resize(tss, width=winsize, fix="center"))
    win[width(win) == winsize]   # drop windows clipped at chromosome ends
}

# ---- Coverage profile --------------------------------------------------------
# For each chromosome, sum coverage across all TSS windows using Views.
# Negative-strand windows are reversed so position 1 = most upstream.
# Accumulation via Reduce("+") stays at C level -- no large matrix in memory.

computeProfile <- function(bedFile, windows, halfwin=2000L, binsize=10L) {
    cat("Loading reads:", bedFile, "\n")
    reads <- import(bedFile, format="BED")

    # Center reads on cut sites (ENCODE shift_width = -read_len/2).
    # Tn5-shifted BED: + strand reads START at cut site, - strand reads END at cut site.
    # Shift each read by its own half-width so its midpoint lands at the cut site.
    # Per-read centering handles variable read lengths (trimmed reads differ in length).
    plus_idx  <- which(as.character(strand(reads)) == "+")
    minus_idx <- which(as.character(strand(reads)) == "-")
    other_idx <- which(as.character(strand(reads)) == "*")
    reads[plus_idx]  <- shift(reads[plus_idx],  -(width(reads[plus_idx])  %/% 2L))
    reads[minus_idx] <- shift(reads[minus_idx], +(width(reads[minus_idx]) %/% 2L))
    reads[other_idx] <- shift(reads[other_idx], -(width(reads[other_idx]) %/% 2L))

    cov <- coverage(reads)

    winsize <- 2L * halfwin + binsize  # 4010 bp -> 401 bins centered at -2000,-1990,...,0,...,2000
    nbins   <- winsize %/% binsize

    bin_sum <- numeric(nbins)
    n_total <- 0L

    chroms <- intersect(names(cov), unique(as.character(seqnames(windows))))

    cat("Aggregating coverage over", length(windows), "TSS windows...\n")

    for (chr in chroms) {
        cv <- cov[[chr]]
        w  <- windows[seqnames(windows) == chr]

        # Filter to windows fully within the coverage Rle.
        # Rle length = position of last read, which may be shorter than the
        # chromosome length. Out-of-bounds Views return integer(0), and
        # bin_sum + numeric(0) = numeric(0), silently wiping the accumulator.
        chr_len <- length(cv)
        w <- w[end(w) <= chr_len]

        for (s in c("+", "-")) {
            ws <- w[strand(w) == s]
            if (length(ws) == 0L) next

            v       <- Views(cv, start=start(ws), end=end(ws))
            # as(v, "IntegerList") correctly expands each RleView into a
            # full integer vector of length winsize. viewApply(v, as.integer)
            # does NOT work -- it applies as.integer() to the RleViews object
            # itself, returning a length-1 value per view.
            pos_sum <- as.numeric(Reduce("+", as(v, "IntegerList")))

            if (length(pos_sum) != winsize) next   # safety: skip malformed sums

            if (s == "-") pos_sum <- rev(pos_sum)

            # Bin: matrix is (binsize x nbins), colSums gives per-bin totals
            bin_sum <- bin_sum + colSums(matrix(pos_sum, nrow=binsize))
            n_total <- n_total + length(ws)
        }
    }

    profile   <- bin_sum / n_total
    # Bin centers: -2000, -1990, ..., 0, ..., 1990, 2000 (401 bins, odd, TSS at 0)
    positions <- seq(-halfwin, halfwin, by=binsize)

    list(profile=profile, positions=positions, n=n_total)
}

# ---- Normalize and score -----------------------------------------------------

normProfile <- function(res, bg_bp=200L, score_bp=200L) {
    binsize <- res$positions[2] - res$positions[1]
    bg_n    <- round(bg_bp / binsize)
    nbins   <- length(res$profile)

    bg   <- mean(c(res$profile[seq_len(bg_n)], res$profile[(nbins - bg_n + 1L):nbins]))
    norm <- res$profile / bg

    score <- max(norm[abs(res$positions) <= score_bp])
    list(norm=norm, positions=res$positions, score=score, n=res$n)
}

# ---- Plot --------------------------------------------------------------------

plotProfile <- function(normed, sample_name) {
    df <- data.frame(pos=normed$positions, enrichment=normed$norm)
    ggplot(df, aes(pos, enrichment)) +
        geom_line(color="steelblue", linewidth=0.7) +
        geom_vline(xintercept=0, linetype="dashed", color="grey50", linewidth=0.5) +
        annotate("text", x=Inf, y=Inf,
            label=sprintf("TSS enrichment: %.2f\nn = %s TSSs",
                normed$score, formatC(normed$n, format="d", big.mark=",")),
            hjust=1.1, vjust=1.5, size=3.5) +
        labs(title=sample_name, x="Distance to TSS (bp)", y="Normalized enrichment") +
        theme_light(base_size=14)
}

# ---- Main --------------------------------------------------------------------

cat("Genome:", GENOME, "\n")
cat("Sample:", SAMPLE_NAME, "\n")
cat("Center shift: per-read (width %/% 2), handles variable read lengths\n")

txdb    <- getTxDb(GENOME)
windows <- getTSSWindows(txdb, HALFWIN, BINSIZE)

res    <- computeProfile(BED_FILE, windows, HALFWIN, BINSIZE)
normed <- normProfile(res, BG_BP, SCORE_BP)

cat(sprintf("TSS enrichment score: %.4f  (n = %s TSSs)\n",
    normed$score, formatC(normed$n, format="d", big.mark=",")))

score_file <- paste0(OUT_PREFIX, "_tss_enrich.csv")
write.csv(
    data.frame(SampleName=SAMPLE_NAME, Score=round(normed$score, 4), N_TSS=normed$n),
    score_file, row.names=FALSE
)
cat("Score:", score_file, "\n")

pdf_file <- paste0(OUT_PREFIX, "_tss_enrich_profile.pdf")
pdf(pdf_file, width=8, height=5)
print(plotProfile(normed, SAMPLE_NAME))
invisible(dev.off())
cat("Plot:", pdf_file, "\n")
