args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 5) {
    cat("
================================================================================
  TSS Enrichment Validation: Profile Comparison and Parameter Matching
================================================================================

USAGE:
  Rscript validateTSS.R GENOME SHIFTED_BED UNSHIFTED_BED BAM_FILE SAMPLE_NAME [OUT_PREFIX]

ARGUMENTS:
  GENOME        Genome: hg19, b38, or mouse
  SHIFTED_BED   Tn5-shifted reads BED file (can be gzipped)
  UNSHIFTED_BED Unshifted reads BED file (can be gzipped)
                  Create from BAM: samtools view -b BAM | bedtools bamtobed -i - | gzip > unshifted.bed.gz
  BAM_FILE      Name-sorted BAM for ATACseqQC (samtools sort -n input.bam -o namesort.bam)
  SAMPLE_NAME   Label for plots
  OUT_PREFIX    Output file prefix (default: SAMPLE_NAME)

OUTPUT:
  <OUT_PREFIX>_validate.pdf   Multi-panel PDF with all tests

TESTS:
  Test 1: Profile shape -- both methods on same sample (visual agreement check)
  Test 2: Shift sensitivity -- shifted vs unshifted profiles (validates Tn5 shift)
  Test 4: Parameter-matched ATACseqQC vs our method (checks parameterization)
================================================================================
")
    quit()
}

GENOME        <- args[1]
SHIFTED_BED   <- args[2]
UNSHIFTED_BED <- args[3]
BAM_FILE      <- args[4]
SAMPLE_NAME   <- args[5]
OUT_PREFIX    <- if (length(args) >= 6) args[6] else SAMPLE_NAME

# ---- Parameters matching computeTSSEnrichment.R --------------------------------

HALFWIN  <- 2000L
BINSIZE  <- 10L
BG_BP    <- 200L
SCORE_BP <- 200L

# ---- Dependencies --------------------------------------------------------------

suppressPackageStartupMessages({
    library(rtracklayer)
    library(GenomicRanges)
    library(GenomicFeatures)
    library(ggplot2)
    library(ATACseqQC)
})

# ---- TxDb selection (from computeTSSEnrichment.R) ------------------------------

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

getTSSWindows <- function(txdb, halfwin=2000L) {
    tss <- unique(resize(transcripts(txdb), width=1L, fix="start"))
    std_chrs <- paste0("chr", c(1:22, "X", "Y"))
    tss <- keepSeqlevels(tss, intersect(std_chrs, seqlevels(tss)), pruning.mode="coarse")
    win <- trim(resize(tss, width=2L * halfwin, fix="center"))
    win[width(win) == 2L * halfwin]
}

# ---- Coverage profile (from computeTSSEnrichment.R) ----------------------------

computeProfile <- function(bedFile, windows, halfwin=2000L, binsize=10L) {
    cat("Loading reads:", bedFile, "\n")
    reads <- import(bedFile, format="BED")
    cov   <- coverage(reads)

    winsize <- 2L * halfwin
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
    positions <- seq(-halfwin + binsize / 2L, halfwin - binsize / 2L, by=binsize)

    list(profile=profile, positions=positions, n=n_total)
}

# ---- Normalize and score (from computeTSSEnrichment.R) -------------------------

normProfile <- function(res, bg_bp=200L, score_bp=200L) {
    binsize <- res$positions[2] - res$positions[1]
    bg_n    <- round(bg_bp / binsize)
    nbins   <- length(res$profile)

    bg   <- mean(c(res$profile[seq_len(bg_n)], res$profile[(nbins - bg_n + 1L):nbins]))
    norm <- res$profile / bg

    score <- max(norm[abs(res$positions) <= score_bp])
    list(norm=norm, positions=res$positions, score=score, n=res$n)
}

# ---- Setup ---------------------------------------------------------------------

cat("Genome:       ", GENOME, "\n")
cat("Sample:       ", SAMPLE_NAME, "\n")
cat("Shifted BED:  ", SHIFTED_BED, "\n")
cat("Unshifted BED:", UNSHIFTED_BED, "\n")
cat("BAM file:     ", BAM_FILE, "\n\n")

txdb    <- getTxDb(GENOME)
windows <- getTSSWindows(txdb, HALFWIN)
txs     <- transcripts(txdb)

# ---- Test 1: Profile shape comparison ------------------------------------------
# Both methods applied to the same sample (shifted BED vs ATACseqQC from BAM).
# Qualitative check: both profiles must peak at position 0.

cat("=== Test 1: Profile shape comparison ===\n")

cat("Running BED method (shifted)...\n")
res1    <- computeProfile(SHIFTED_BED, windows, HALFWIN, BINSIZE)
normed1 <- normProfile(res1, BG_BP, SCORE_BP)
cat(sprintf("  BED method score: %.4f  (n = %s TSSs)\n",
    normed1$score, formatC(normed1$n, format="d", big.mark=",")))

cat("Running ATACseqQC (default parameters)...\n")
bamReads <- readBamFile(BAM_FILE, asMates=TRUE)
shifted_bam <- shiftGAlignmentsList(bamReads)
tsse_default <- TSSEscore(shifted_bam, txs)
cat(sprintf("  ATACseqQC score (default): %.4f\n", tsse_default$TSSEscore))

# ATACseqQC tsse$values: 20 windows of 100 bp covering -1000 to +1000 bp
# Center positions: seq(-950, 950, by=100)
ataqc_pos <- seq(-950, 950, by=100)
ataqc_val <- tsse_default$values

df_bed <- data.frame(
    pos         = normed1$positions,
    enrichment  = normed1$norm,
    method      = "BED (10 bp bins, +-2000 bp)"
)
# Scale ATACseqQC values to share the same y-axis: normalize to its background
# ATACseqQC already returns normalized values, so use them directly
ataqc_bg <- mean(c(head(ataqc_val, 2), tail(ataqc_val, 2)))
ataqc_norm <- if (ataqc_bg > 0) ataqc_val / ataqc_bg else ataqc_val
df_ataqc <- data.frame(
    pos         = ataqc_pos,
    enrichment  = ataqc_norm,
    method      = "ATACseqQC (100 bp bins, +-1000 bp)"
)

p_test1 <- ggplot(df_bed, aes(pos, enrichment)) +
    geom_line(aes(color=method), linewidth=0.7) +
    geom_line(data=df_ataqc, aes(color=method), linewidth=1.2) +
    geom_vline(xintercept=0, linetype="dashed", color="grey50", linewidth=0.5) +
    annotate("text", x=Inf, y=Inf,
        label=sprintf("BED score: %.2f\nATACseqQC score: %.2f",
            normed1$score, tsse_default$TSSEscore),
        hjust=1.1, vjust=1.5, size=3.5) +
    scale_color_manual(values=c("BED (10 bp bins, +-2000 bp)"="steelblue",
                                "ATACseqQC (100 bp bins, +-1000 bp)"="tomato")) +
    labs(title=sprintf("%s -- Test 1: Profile shape comparison", SAMPLE_NAME),
        subtitle="Both profiles should peak at position 0",
        x="Distance to TSS (bp)", y="Normalized enrichment", color=NULL) +
    theme_light(base_size=12) +
    theme(legend.position="bottom")

# ---- Test 2: Shift sensitivity -------------------------------------------------
# BED method applied to shifted vs unshifted reads.
# Unshifted peak should be at ~-4 bp; shifted peak should be at 0.

cat("\n=== Test 2: Shift sensitivity ===\n")

same_file <- (SHIFTED_BED == UNSHIFTED_BED)
if (same_file) {
    cat("  NOTE: SHIFTED_BED == UNSHIFTED_BED -- Test 2 shows a single profile.\n")
    cat("  To enable Test 2, provide an unshifted BED:\n")
    cat("    samtools view -b BAM | bedtools bamtobed -i - | gzip > unshifted.bed.gz\n")
}

cat("Running BED method (unshifted)...\n")
res2    <- computeProfile(UNSHIFTED_BED, windows, HALFWIN, BINSIZE)
normed2 <- normProfile(res2, BG_BP, SCORE_BP)
cat(sprintf("  BED method score (unshifted): %.4f  (n = %s TSSs)\n",
    normed2$score, formatC(normed2$n, format="d", big.mark=",")))

df_shifted   <- data.frame(pos=normed1$positions, enrichment=normed1$norm,
                            reads="Shifted (+4/-5 bp)")
df_unshifted <- data.frame(pos=normed2$positions, enrichment=normed2$norm,
                            reads="Unshifted")

subtitle2 <- if (same_file) {
    "Provide separate unshifted BED to see the ~4 bp peak offset"
} else {
    "Shifted peak at 0; unshifted peak at ~-4 bp (+ strand bias)"
}

p_test2 <- ggplot(df_shifted, aes(pos, enrichment)) +
    geom_line(aes(color=reads), linewidth=0.7) +
    geom_line(data=df_unshifted, aes(color=reads), linewidth=0.7) +
    geom_vline(xintercept=0, linetype="dashed", color="grey50", linewidth=0.5) +
    geom_vline(xintercept=-4, linetype="dotted", color="grey70", linewidth=0.5) +
    annotate("text", x=Inf, y=Inf,
        label=sprintf("Shifted score: %.2f\nUnshifted score: %.2f",
            normed1$score, normed2$score),
        hjust=1.1, vjust=1.5, size=3.5) +
    scale_color_manual(values=c("Shifted (+4/-5 bp)"="steelblue",
                                "Unshifted"="darkorange")) +
    labs(title=sprintf("%s -- Test 2: Shift sensitivity", SAMPLE_NAME),
        subtitle=subtitle2,
        x="Distance to TSS (bp)", y="Normalized enrichment", color=NULL) +
    theme_light(base_size=12) +
    theme(legend.position="bottom")

# ---- Test 4: Parameter-matched ATACseqQC run -----------------------------------
# Run TSSEscore with parameters matching our BED method.
# If scores converge, the difference is purely parameterization.

cat("\n=== Test 4: Parameter-matched ATACseqQC ===\n")

tsse_matched <- tryCatch({
    cat("Running ATACseqQC with matched parameters (upstream=2000, width=10, step=10)...\n")
    res <- TSSEscore(shifted_bam, txs,
        upstream=2000L, downstream=2000L,
        endSize=200L,
        width=10L, step=10L)
    cat(sprintf("  ATACseqQC score (matched): %.4f\n", res$TSSEscore))
    res
}, error=function(e) {
    cat("  ATACseqQC with matched parameters failed:", conditionMessage(e), "\n")
    cat("  Check if TSSEscore supports these arguments in your installed version.\n")
    NULL
})

# Build Test 4 plot
if (!is.null(tsse_matched)) {
    # matched parameters: 400 bins of 10 bp across +-2000 bp
    n_matched <- length(tsse_matched$values)
    matched_pos <- seq(-2000 + 5, 2000 - 5, length.out=n_matched)
    matched_bg <- mean(c(head(tsse_matched$values, 20), tail(tsse_matched$values, 20)))
    matched_norm <- if (matched_bg > 0) tsse_matched$values / matched_bg else tsse_matched$values

    df_matched <- data.frame(pos=matched_pos, enrichment=matched_norm,
                              method="ATACseqQC (matched: 10 bp bins, +-2000 bp)")
    df_ours    <- data.frame(pos=normed1$positions, enrichment=normed1$norm,
                              method="BED method (10 bp bins, +-2000 bp)")

    p_test4 <- ggplot(df_ours, aes(pos, enrichment)) +
        geom_line(aes(color=method), linewidth=0.7) +
        geom_line(data=df_matched, aes(color=method), linewidth=0.7) +
        geom_vline(xintercept=0, linetype="dashed", color="grey50", linewidth=0.5) +
        annotate("text", x=Inf, y=Inf,
            label=sprintf("BED score: %.2f\nATACseqQC matched: %.2f",
                normed1$score, tsse_matched$TSSEscore),
            hjust=1.1, vjust=1.5, size=3.5) +
        scale_color_manual(values=c(
            "BED method (10 bp bins, +-2000 bp)"="steelblue",
            "ATACseqQC (matched: 10 bp bins, +-2000 bp)"="forestgreen")) +
        labs(title=sprintf("%s -- Test 4: Parameter-matched ATACseqQC", SAMPLE_NAME),
            subtitle="If scores converge, difference was parameterization only",
            x="Distance to TSS (bp)", y="Normalized enrichment", color=NULL) +
        theme_light(base_size=12) +
        theme(legend.position="bottom")
} else {
    p_test4 <- ggplot() +
        annotate("text", x=0.5, y=0.5,
            label="Test 4 failed: ATACseqQC does not support matched parameters\nin this version. See console output.",
            hjust=0.5, vjust=0.5, size=5) +
        theme_void() +
        labs(title=sprintf("%s -- Test 4: Parameter-matched ATACseqQC (FAILED)", SAMPLE_NAME))
}

# ---- Score summary table -------------------------------------------------------

cat("\n=== Score Summary ===\n")
cat(sprintf("  BED method (shifted, 10 bp bins, +-2000 bp):      %.4f\n", normed1$score))
cat(sprintf("  BED method (unshifted, 10 bp bins, +-2000 bp):    %.4f\n", normed2$score))
cat(sprintf("  ATACseqQC (default: 100 bp bins, +-1000 bp):      %.4f\n", tsse_default$TSSEscore))
if (!is.null(tsse_matched)) {
    cat(sprintf("  ATACseqQC (matched: 10 bp bins, +-2000 bp):       %.4f\n", tsse_matched$TSSEscore))
}

# ---- Write output --------------------------------------------------------------

pdf_file <- paste0(OUT_PREFIX, "_validate.pdf")
pdf(pdf_file, width=10, height=6)
print(p_test1)
print(p_test2)
print(p_test4)
invisible(dev.off())
cat("\nPlots written to:", pdf_file, "\n")

score_file <- paste0(OUT_PREFIX, "_validate_scores.txt")
lines <- c(
    paste("Sample", "Method", "Score", sep="\t"),
    paste(SAMPLE_NAME, "BED_shifted",      sprintf("%.4f", normed1$score),         sep="\t"),
    paste(SAMPLE_NAME, "BED_unshifted",    sprintf("%.4f", normed2$score),         sep="\t"),
    paste(SAMPLE_NAME, "ATACseqQC_default",sprintf("%.4f", tsse_default$TSSEscore),sep="\t")
)
if (!is.null(tsse_matched)) {
    lines <- c(lines,
        paste(SAMPLE_NAME, "ATACseqQC_matched", sprintf("%.4f", tsse_matched$TSSEscore), sep="\t"))
}
writeLines(lines, score_file)
cat("Scores written to:", score_file, "\n")
