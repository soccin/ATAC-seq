args <- commandArgs(trailingOnly=TRUE)

DEBUG <- "-d" %in% args
args  <- args[args != "-d"]

if (length(args) < 2) {
    cat("
================================================================================
  ENCODE-exact TSS Enrichment Score and Profile
================================================================================

USAGE:
  Rscript encodeTSS.R [-d] GENOME BED_FILE [SAMPLE_NAME] [OUT_PREFIX]

ARGUMENTS:
  -d           Debug: dump full normalized enrichment profile to CSV
  GENOME       Genome: hg19, b38, or mouse
  BED_FILE     Tn5-shifted reads BED file (can be gzipped)
  SAMPLE_NAME  Label for plots (default: BED_FILE basename)
  OUT_PREFIX   Output file prefix (default: BED_FILE without extension)

OUTPUT:
  <OUT_PREFIX>_tss_enrich.csv          CSV with header: SampleName, Score, N_TSS
  <OUT_PREFIX>_tss_enrich_profile.pdf  Profile plot
  <OUT_PREFIX>_tss_enrich_debug.csv    (with -d) Position and NormalizedEnrichment

DETAILS:
  Direct R translation of encode_task_tss_enrich.py make_tss_plot().
  Source: github.com/ENCODE-DCC/atac-seq-pipeline

  Parameters matching ENCODE defaults:
    bins=401, bp_edge=2000   -> 401 bins of 10 bp over +-2000 bp, bin centered at 0
    greenleaf_norm=True
    num_edge_bins = 10       -> background from outermost 100 bp each side
    score = max over all 401 bins (no +-200 bp restriction)

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

# ENCODE hardcoded defaults (encode_task_tss_enrich.py make_tss_plot())
HALFWIN       <- 2000L  # bp_edge
BINSIZE       <- 10L    # (2*bp_edge) / bins = 4000/400
NUM_EDGE_BINS <- 10L    # int(100 / (2*bp_edge/bins)) = int(100/10)

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
# Same as computeTSSEnrichment.R: strand-aware resize to get exact TSS position,
# then expand to 2*HALFWIN window. Drop windows clipped at chromosome ends.

getTSSWindows <- function(txdb, halfwin=2000L, binsize=10L) {
    tss <- unique(resize(transcripts(txdb), width=1L, fix="start"))
    std_chrs <- paste0("chr", c(1:22, "X", "Y"))
    tss <- keepSeqlevels(tss, intersect(std_chrs, seqlevels(tss)), pruning.mode="coarse")
    winsize <- 2L * halfwin + binsize
    win <- trim(resize(tss, width=winsize, fix="center"))
    win[width(win) == winsize]
}

# ---- Coverage profile --------------------------------------------------------
# Matches ENCODE stranded=True: negative-strand windows reversed so pos 1 = TSS.
# Same implementation as computeTSSEnrichment.R computeProfile().

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

    winsize <- 2L * halfwin + binsize  # 4010 bp -> 401 bins centered at -2000,...,0,...,2000
    nbins   <- winsize %/% binsize

    bin_sum <- numeric(nbins)
    n_total <- 0L

    chroms <- intersect(names(cov), unique(as.character(seqnames(windows))))

    cat("Aggregating coverage over", length(windows), "TSS windows...\n")

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

# ---- Normalize and score (ENCODE Greenleaf method) ---------------------------
# encode_task_tss_enrich.py:
#   num_edge_bins = int(100 / (2*bp_edge/bins))            # = 10
#   avg_noise = (sum(bin_means[:10]) + sum(bin_means[-10:])) / (2*10)
#   bam_array /= avg_noise
#   tss_point_val = max(bam_array.mean(axis=0))            # max over all bins

normProfile <- function(res, num_edge_bins=10L) {
    nbins <- length(res$profile)

    # Background: mean of outermost num_edge_bins on each side
    # Python: (sum(bin_means[:10]) + sum(bin_means[-10:])) / (2*10)
    edge_vals <- c(res$profile[seq_len(num_edge_bins)],
                   res$profile[(nbins - num_edge_bins + 1L):nbins])
    avg_noise <- mean(edge_vals)

    norm <- res$profile / avg_noise

    # Score: max over all 401 bins (no +-200 bp restriction)
    score <- max(norm)

    list(norm=norm, positions=res$positions, score=score, n=res$n)
}

# ---- Plot --------------------------------------------------------------------

plotProfile <- function(normed, sample_name) {
    df <- data.frame(pos=normed$positions, enrichment=normed$norm)
    ggplot(df, aes(pos, enrichment)) +
        geom_line(color="steelblue", linewidth=0.7) +
        geom_vline(xintercept=0, linetype="dashed", color="grey50", linewidth=0.5) +
        annotate("text", x=Inf, y=Inf,
            label=sprintf("ENCODE TSS enrichment: %.2f\nn = %s TSSs",
                normed$score, formatC(normed$n, format="d", big.mark=",")),
            hjust=1.1, vjust=1.5, size=3.5) +
        labs(title=sample_name,
             subtitle="ENCODE-exact: background=100bp edges, score=max over +-2000bp",
             x="Distance to TSS (bp)", y="Normalized enrichment") +
        theme_light(base_size=14)
}

# ---- Main --------------------------------------------------------------------

cat("Genome:", GENOME, "\n")
cat("Sample:", SAMPLE_NAME, "\n")
cat("Center shift: per-read (width %/% 2), handles variable read lengths\n")
cat("Method: ENCODE-exact (background=100bp, score=max all bins)\n")

txdb    <- getTxDb(GENOME)
windows <- getTSSWindows(txdb, HALFWIN, BINSIZE)

res    <- computeProfile(BED_FILE, windows, HALFWIN, BINSIZE)
normed <- normProfile(res, NUM_EDGE_BINS)

cat(sprintf("ENCODE TSS enrichment score: %.4f  (n = %s TSSs)\n",
    normed$score, formatC(normed$n, format="d", big.mark=",")))

score_file <- paste0(OUT_PREFIX, "_tss_enrich.csv")
write.csv(
    data.frame(SampleName=SAMPLE_NAME, Score=round(normed$score, 4), N_TSS=normed$n),
    score_file, row.names=FALSE
)
cat("Score:", score_file, "\n")

if (DEBUG) {
    debug_file <- paste0(OUT_PREFIX, "_tss_enrich_debug.csv")
    write.csv(
        data.frame(Position_bp=normed$positions, NormalizedEnrichment=normed$norm),
        debug_file, row.names=FALSE
    )
    cat("Debug profile:", debug_file, "\n")
}

pdf_file <- paste0(OUT_PREFIX, "_tss_enrich_profile.pdf")
pdf(pdf_file, width=8, height=5)
print(plotProfile(normed, SAMPLE_NAME))
invisible(dev.off())
cat("Plot:", pdf_file, "\n")
