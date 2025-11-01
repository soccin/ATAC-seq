args <- commandArgs(trailing = TRUE)
if (len(args) < 3) {
  cat("
================================================================================
  Pairwise Differential ChIP-seq Peak Analysis using edgeR
================================================================================

USAGE:
  Rscript diffAnalysisPairwise.R GENOME MANIFEST COMPARISONS [RUNTAG]

REQUIRED ARGUMENTS:
  GENOME        Genome assembly: 'human' or 'mouse'
  MANIFEST      Sample manifest CSV file (with header row)
  COMPARISONS   Comparison pairs CSV file (NO header row)

OPTIONAL ARGUMENTS:
  RUNTAG        Tag to append to output filenames (default: none)

OUTPUT FILES:
  - DiffPeaksEdgeR_V3.xlsx : Excel file with differential peaks
  - DiffPeaks_V3.pdf       : PDF with PCA and MA/volcano plots

--------------------------------------------------------------------------------
MANIFEST FILE FORMAT (with header):
  Three columns: MapID, SampleID, Group

  Example:
    MapID,SampleID,Group
    761-wt,761-wt,wt
    810-wt,810-wt,wt
    876-ko,876-ko,ko
    978-ko,978-ko,ko

  See example: R/sampleManifest.csv

--------------------------------------------------------------------------------
COMPARISONS FILE FORMAT (NO header):
  Two columns: Group1, Group2
  Each row defines one pairwise comparison

  Example:
    wt,ko

  This compares Group2 vs Group1 (ko vs wt)
  Sign convention: Group2 - Group1
  Positive logFC means higher in Group2 (ko)
  Negative logFC means higher in Group1 (wt)

  See example: R/comparisions.csv

--------------------------------------------------------------------------------
EXAMPLE USAGE:
  Rscript diffAnalysisPairwise.R mouse R/sampleManifest.csv R/comparisions.csv
  Rscript diffAnalysisPairwise.R mouse manifest.csv comps.csv run1

================================================================================

")
  quit()
}

#' Fix sample names based on naming convention
#'
#' Detects whether samples use PEmap or BIC naming and applies the appropriate
#' transformation to match the sample manifest.
#'
#' @param sample_names Character vector of raw sample names
#' @return Character vector of cleaned sample names matching manifest
fix_sample_names <- function(sample_names) {
  if (grepl("___MD", sample_names[1])) {
    fix_sample_names_pemap(sample_names)
  } else {
    fix_sample_names_bic(sample_names)
  }
}

#' Fix sample names using BIC naming convention
#'
#' Removes processing suffixes and path prefixes, then maps to manifest IDs.
#'
#' @param sample_names Character vector of BIC-format sample names
#' @return Character vector of cleaned sample names
fix_sample_names_bic <- function(sample_names) {
  sample_names |>
    gsub("_postProcess.*", "", x = _) |>
    gsub(".*_s_", "s_", x = _) |>
    (\(x) sampRename[x])() |>
    unname()
}

#' Fix sample names using PEmap naming convention
#'
#' Extracts basename and removes PEmap suffixes, then maps to manifest IDs.
#'
#' @param sample_names Character vector of PEmap-format sample names
#' @return Character vector of cleaned sample names
fix_sample_names_pemap <- function(sample_names) {
  sample_names |>
    basename() |>
    gsub("___.*", "", x = _) |>
    (\(x) sampRename[x])() |>
    unname()
}

#' Create reverse log transformation for ggplot2 volcano plots
#'
#' Transforms p-values so smaller values appear higher on the y-axis.
#'
#' @param base Numeric base for logarithm (default: e)
#' @return A transformation object for ggplot2 scales
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(
    paste0("reverselog-", format(base)),
    trans,
    inv,
    log_breaks(base = base),
    domain = c(1e-100, Inf)
  )
}

#' Create PNG file with Cairo graphics device
#'
#' @param filename Character path to output file
#' @param width Numeric width in inches (default: 14)
#' @param height Numeric height in inches (default: 8.5)
#' @param pointsize Numeric point size for text (default: 12)
#' @param res Numeric resolution in DPI (default: 150)
png_cairo <- function(filename, width = 14, height = 8.5, pointsize = 12,
                      res = 150) {
  png(
    filename,
    type = "cairo",
    units = "in",
    width = width,
    height = height,
    pointsize = pointsize,
    res = res
  )
}

#' Merge PNG files into a single PDF using ImageMagick
#'
#' @param file_spec Character file specification with printf-style formatting
merge_pngs <- function(file_spec) {
  file_re <- gsub("_%\\d+d", ".*", file_spec)
  pdf_file <- gsub("_%\\d+d.*", ".pdf", file_spec)
  system2(
    "convert",
    c(sort(dir_ls(regex = file_re)), pdf_file),
    stderr = cc("stderr", "mergePNGs", "convert", DATE())
  )
}

cat("Loading libraries ...")
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(AnnotationDbi)
  library(patchwork)
  library(scales)
  library(edgeR)
  library(ggrepel)
  library(ggsci)
  library(tidyverse)
  library(fs)
  library(openxlsx)
  library(ggrastr)
})
cat(" done\n")

# Locate feature counts summary file
# Works with both ATACseq and ChIPSeq pipeline structures
peak_counts_file <- fs::dir_ls(regex = "peaks_raw_fcCounts.txt.summary")

if (len(peak_counts_file) == 0) {
  peak_counts_file <- fs::dir_ls("out/macs",
                                 regex = "peaks_raw_fcCounts.txt.summary")
}

feature_summary <- read_tsv(peak_counts_file)

# Parse command line arguments
GENOME <- args[1]
MANIFEST_FILE <- args[2]
COMPARISON_FILE <- args[3]
RUNTAG <- if (len(args) == 4) args[4] else ""

# Set genome-specific annotation databases
if (GENOME == "human") {
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  annoDb <- "org.Hs.eg.db"
} else if (GENOME == "mouse") {
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  annoDb <- "org.Mm.eg.db"
} else {
  cat("\n\tUnknown GENOME", GENOME, "\n")
  cat("\tValid genomes: human, mouse\n\n")
  quit()
}

# Load sample manifest and create renaming lookup
manifest <- read_csv(MANIFEST_FILE) |> arrange(SampleID)

sampRename <- manifest$SampleID
names(sampRename) <- manifest$MapID

# Process feature counts summary for QC
feature_summary <- feature_summary |>
  gather(Sample, Count, -Status) |>
  mutate(Sample = fix_sample_names(Sample)) |>
  mutate(Status = gsub("_.*", "", Status)) |>
  group_by(Sample, Status) |>
  summarize(Counts = sum(Count), .groups = "drop") |>
  mutate(Status = ifelse(Status == "Assigned", "InPeaks", "Outside"))

# Load raw feature counts matrix
raw_counts <- read_tsv(gsub(".txt.summary$", ".txt", peak_counts_file),
                       comment = "#")

# Extract peak annotation
peak.annote <- raw_counts |>
  select(PeakNo = Geneid, Chr, Start, End, Strand, Length)

# Remove excluded samples from manifest
manifest <- manifest |> filter(!grepl("^EXC", Group))

# Create count matrix for edgeR
count_matrix <- raw_counts |>
  select(PeakNo = Geneid, matches(".bam$")) |>
  data.frame(check.names = FALSE) |>
  column_to_rownames("PeakNo")

colnames(count_matrix) <- fix_sample_names(colnames(count_matrix))

# Remove columns not in manifest (NA values from failed sample name mapping)
count_matrix <- count_matrix[, !is.na(colnames(count_matrix))]

# Validate all count matrix columns are in manifest
if (!all(colnames(count_matrix) %in% manifest$SampleID)) {
  cat("\nERROR: Count matrix columns don't match manifest SampleIDs\n")
  cat("Check sample naming conventions and manifest file\n\n")
  rlang::abort("Count matrix validation failed")
}

# Prepare edgeR DGEList object
count_matrix <- count_matrix[, manifest$SampleID]
group <- factor(manifest$Group)
dge <- DGEList(counts = count_matrix, group = group)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# Perform PCA for quality control
pca_result <- prcomp(cpm(dge, log = TRUE), scale = FALSE)
pca_data <- pca_result$rotation |>
  data.frame() |>
  rownames_to_column("SampleID") |>
  left_join(manifest, by = "SampleID")

# Define color palette for plots
plot_colors <- c(pal_uchicago("default")(9), pals::cols25())

# PCA plot without labels
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = Group, label = SampleID)) +
  theme_light(base_size = 16) +
  geom_point(size = 4, alpha = 0.6) +
  scale_color_manual(values = plot_colors)

# PCA plot with sample labels
pca_plot_labeled <- ggplot(pca_data,
                           aes(PC1, PC2, color = Group, label = SampleID)) +
  theme_light(base_size = 16) +
  geom_point(size = 2) +
  scale_color_manual(values = plot_colors) +
  geom_label_repel(
    color = "black",
    max.overlaps = Inf,
    min.segment.length = 0,
    size = 3,
    force = 16,
    max.time = 10,
    max.iter = 100000
  )

# Set up design matrix and estimate dispersion
design <- model.matrix(~0 + group)
dge <- estimateDisp(dge, design)

#' Perform differential peak analysis using edgeR LRT
#'
#' Subsets DGE object to comparison groups, re-filters peaks, performs
#' likelihood ratio test, and annotates significant peaks.
#'
#' @param dge_object DGEList object containing normalized counts
#' @param design Design matrix from full dataset
#' @param contrast Contrast vector defining comparison
#' @param fdr_cut FDR cutoff for significance (default: 0.05)
#' @return List with results table, comparison name, MA plot, and volcano plot
do_differential_analysis <- function(dge_object, design, contrast,
                                     fdr_cut = 0.05) {

  cat("Num peaks =", nrow(dge_object), "\n")
  cat("\nRefiltering peaks for this comparison\n\n")

  # Extract samples involved in this comparison
  comp_sample_idx <- which(rowSums(design[, contrast != 0]) > 0) |> unname()

  # Subset and re-filter data for comparison-specific analysis
  comp_group <- factor(group[comp_sample_idx])
  comp_counts <- getCounts(dge_object)[, comp_sample_idx]
  comp_design <- model.matrix(~0 + comp_group)

  dge_comp <- DGEList(counts = comp_counts, group = comp_group)
  keep <- filterByExpr(dge_comp, design = comp_design)
  dge_comp <- dge_comp[keep, , keep.lib.sizes = FALSE]
  dge_comp <- calcNormFactors(dge_comp)
  dge_comp <- estimateDisp(dge_comp, comp_design)

  # Perform likelihood ratio test
  comp_contrast <- contrast[contrast != 0]
  fit <- glmFit(dge_comp, comp_design)
  lrt <- glmLRT(fit, contrast = comp_contrast)
  model <- lrt

  # Create comparison tag for labeling
  comparison_string <- model$comparison
  comp_tag <- paste0(sort(strsplit(comparison_string, " ")[[1]],
                          decreasing = TRUE), collapse = "") |>
    gsub("-1\\*", "-", x = _) |>
    gsub("^1\\*", "", x = _) |>
    gsub("gComp", "", x = _)

  # Extract and process results table
  results_table <- model$table |>
    data.frame() |>
    rownames_to_column("PeakNo") |>
    tibble() |>
    mutate(FDR = p.adjust(PValue)) |>
    arrange(FDR, PValue) |>
    mutate(PValue.mod = ifelse(PValue < .Machine$double.eps^2,
                                .Machine$double.eps^2, PValue))

  # Apply IHW (Independent Hypothesis Weighting) for multiple testing
  library(IHW)
  ihw_result <- ihw(PValue ~ Length,
                    dat = results_table |> left_join(peak.annote, by = "PeakNo"),
                    alpha = 0.25)
  cat("IHW rejections =", rejections(ihw_result), "\n")

  # Create MA plot
  max_logfc <- max(abs(results_table$logFC))

  ma_plot <- ggplot(arrange(results_table, desc(PValue)),
                    aes(logCPM, logFC, color = FDR < fdr_cut)) +
    theme_light(base_size = 16) +
    rasterize(geom_point()) +
    scale_color_manual(values = c("#7f7f7f33", "#e31a1c")) +
    ggtitle(comp_tag) +
    scale_y_continuous(limits = c(-1, 1) * max_logfc)

  # Create volcano plot (PValue.mod clipped for readability)
  volcano_plot <- ggplot(arrange(results_table, desc(PValue)),
                         aes(logFC, PValue.mod, color = FDR < fdr_cut)) +
    theme_light(base_size = 16) +
    rasterize(geom_point()) +
    scale_y_continuous(trans = reverselog_trans(10)) +
    scale_x_continuous(limits = c(-1, 1) * max_logfc) +
    scale_color_manual(values = c("#7f7f7f33", "#e31a1c")) +
    ggtitle(comp_tag)

  # Filter for significant peaks (exclude mitochondrial)
  sig_results <- results_table |>
    filter(FDR < fdr_cut) |>
    left_join(peak.annote, by = "PeakNo") |>
    select(-matches("^F$|^LR")) |>
    filter(Chr != "MT")

  # Prepare significant peaks for ChIPseeker annotation
  sig_peaks_bed <- sig_results |>
    select(V1 = Chr, V2 = Start, V3 = End, V4 = PeakNo, V5 = PValue) |>
    mutate(V5 = -10 * log10(V5))

  # Annotate significant peaks with nearest genes
  if (nrow(sig_results) > 0) {
    # Add "chr" prefix if missing (genome-dependent)
    if (substr(sig_results$Chr[1], 1, 3) != "chr") {
      sig_peaks_bed <- sig_peaks_bed |>
        mutate(V1 = paste0("chr", V1)) |>
        data.frame()
    } else {
      sig_peaks_bed <- sig_peaks_bed |> data.frame()
    }

    # Convert to GRanges and annotate with TSS regions
    sig_peaks_granges <- ChIPseeker:::peakDF2GRanges(sig_peaks_bed)
    peak_annotations <- annotatePeak(sig_peaks_granges,
                                     TxDb = txdb,
                                     tssRegion = c(-5000, 5000),
                                     annoDb = annoDb)

    # Extract gene information
    gene_info <- as.data.frame(peak_annotations) |>
      tibble() |>
      dplyr::select(PeakNo = V4, SYMBOL, GENENAME, annotation, distanceToTSS)

    final_table <- sig_results |> left_join(gene_info, by = "PeakNo")
  } else {
    cat("\n   No significant peaks for", comp_tag, "at FDR", fdr_cut, "\n\n")
    final_table <- sig_results
  }

  list(
    tbl = final_table,
    comparison = comp_tag,
    p.ma = ma_plot,
    p.vc = volcano_plot,
    model = model
  )
}

# Extract project number from working directory for output naming
working_dir_parts <- strsplit(getwd(), "/")[[1]]
project_id <- grep("^Proj_|^B-\\d+", working_dir_parts, value = TRUE)

# Load comparisons file and extract group names from design matrix
group_names <- gsub("group", "", colnames(design))
comparisons <- read_csv(COMPARISON_FILE, col_names = FALSE)

# Run differential analysis for each comparison
results_list <- list()

for (comparison_pair in transpose(comparisons)) {
  # Create contrast vector: Group2 - Group1
  contrast <- as.numeric(group_names == comparison_pair$X2) -
              as.numeric(group_names == comparison_pair$X1)

  cat("========================================================\n")
  cat(str(comparison_pair))

  # Validate contrast is valid (must have exactly 2 groups)
  if (!(sum(contrast != 0) > 0 & sum(contrast) == 0)) {
    cat("\nERROR: Invalid contrast specification\n")
    cat("Check comparison file format\n\n")
    stop("Invalid contrast")
  }

  results_list[[len(results_list) + 1]] <-
    do_differential_analysis(dge, design, contrast)
}
cat("========================================================\n")
cat("========================================================\n\n")

# Compile results tables with comparison names
result_tables <- map(results_list, "tbl")
names(result_tables) <- map(results_list, "comparison") |>
  unlist() |>
  substr(1, 31)

# Create summary statistics
summary_stats <- map(result_tables, nrow) |>
  bind_rows() |>
  gather(Comparison, NumSig)

# Write results to Excel file
output_xlsx <- cc(project_id, RUNTAG, "DiffPeaksEdgeR_V3.xlsx")
write.xlsx(
  c(list(Summary = summary_stats), result_tables),
  output_xlsx
)

# Generate PDF with PCA and differential analysis plots
output_pdf <- cc(project_id, RUNTAG, "DiffPeaks_V3.pdf")
pdf(output_pdf, width = 11, height = 8.5)

print(pca_plot)
print(pca_plot_labeled)

for (i in seq(len(results_list))) {
  print(results_list[[i]]$p.ma)
  print(results_list[[i]]$p.vc)
}

dev.off()
