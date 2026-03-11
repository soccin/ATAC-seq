args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) {
    cat("\n   Usage: Rscript testTSS.R NAME_SORTED_BAM\n\n")
    cat("   BAM must be sorted by name: samtools sort -n input.bam -o input_namesort.bam\n\n")
    quit()
}

BAM_FILE <- args[1]

suppressPackageStartupMessages({
    library(ATACseqQC)
    library(GenomicFeatures)
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
})

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txs  <- transcripts(txdb)

cat("BAM file:", BAM_FILE, "\n")
cat("Transcripts:", length(txs), "\n")

bamReads <- readBamFile(BAM_FILE, asMates=TRUE)
shifted  <- shiftGAlignmentsList(bamReads)
tsse     <- TSSEscore(shifted, txs,
    upstream=2000L, downstream=2000L,
    endSize=200L,
    width=10L, step=10L)

cat("ATACseqQC TSS enrichment score:", tsse$TSSEscore, "\n")

out      <- data.frame(Sample=sub("\\.bam$", "", basename(BAM_FILE)), TSSEscore=tsse$TSSEscore)
csv_file <- sub("\\.bam$", "_TSSEscore_matched.csv", BAM_FILE)
write.csv(out, csv_file, row.names=FALSE)
cat("Score written to:", csv_file, "\n")
