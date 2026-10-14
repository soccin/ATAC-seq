# BIC ATAC-seq Pipeline (v2.0.0)

## Introduction

The results are delivered in one folder, `atacSeq`. `macs` holds the
peaks called for each sample. `atlas` holds the merged peak atlas of
the project and the raw count matrix over it, which is the input to
any normalization or differential analysis. `bigwig` holds a
normalized signal track per sample for a genome browser. `metrics`
holds the QC metrics and plots. Projects with replicates have a
`diff` folder with the differential accessibility analysis. `docs`
holds this document, `output.pdf`, and the methods document,
`methods.pdf`, which describes how each file was made.

In the names below, `<sample>` is the sample ID of the BAM and `<no>`
is the project number (for example `18143`). Coordinates are 0-based
with an exclusive end (BED convention) unless a file says otherwise.
Peaks, the atlas and the signal tracks cover the autosomes and the X
and Y chromosomes only; mitochondrial, unplaced and other contigs are
left out.

## Directory Structure

```
atacSeq
├── atlas
│   ├── macsPeaksMerged.saf
│   ├── peaks_raw_fcCounts.txt
│   └── peaks_raw_fcCounts.txt.summary
├── bigwig
│   └── <sample>_postProcess.shifted.10mNorm.bw
├── macs
│   └── <sample>_postProcess.shifted
│       ├── <sample>_postProcess.shifted_peaks.narrowPeak
│       ├── <sample>_postProcess.shifted_summits.bed
│       ├── <sample>_postProcess.shifted_peaks.xls
│       ├── <sample>_postProcess.shifted_peaks.blacklistFiltered.narrowPeak
│       └── <sample>_postProcess.shifted_summits.blacklistFiltered.bed
├── metrics
│   ├── Proj_<no>__postInsDistribution.pdf
│   ├── Proj_<no>__ATACSeqQC.pdf
│   ├── Proj_<no>__ATACSeqDepth.csv
│   ├── <sample>_postProcess___INS.txt
│   ├── <sample>_postProcess___INS.pdf
│   ├── <sample>_postProcess___MITO.txt
│   ├── <sample>_postProcess.tss_enrich.csv
│   ├── <sample>_postProcess.tss_enrich.pdf
│   └── <sample>_postProcess.large_tss_enrich.png
├── docs
│   ├── output.pdf
│   └── methods.pdf
└── diff                      (projects with replicates only)
    ├── Proj_<no>__DiffPeaksEdgeR_V4.xlsx
    └── Proj_<no>__DiffPeaks_V4.pdf
```

The filtered BAM files (`<sample>_postProcess.bam`) are not delivered.
If you need them, please contact us.

## Peaks (`macs`)

One folder per sample with the output of MACS2 and two files derived
from it. The folder and file names carry the suffix
`_postProcess.shifted` because the peaks were called on the filtered,
Tn5-shifted reads.

**MACS2 output.** These three files are MACS2's own, delivered
unchanged:

- `*_peaks.narrowPeak`: every peak MACS2 called at p < 0.01, in
  ENCODE narrowPeak format (BED6+4). The columns are chromosome,
  start, end, peak name, score, strand (always `.`), fold enrichment
  at the summit, -log10 p-value, -log10 q-value, and the position of
  the summit relative to the peak start. A peak with several summits
  appears once per summit, with a letter suffix on its name
  (`_peak_12a`, `_peak_12b`).
- `*_summits.bed`: one line per summit, the single base of highest
  signal in each peak, with the peak name and its -log10 p-value.
- `*_peaks.xls`: the same peaks as a tab-delimited text table (not an
  Excel file) with 1-based coordinates, peak length, summit position,
  pileup, -log10 p-value, fold enrichment and -log10 q-value. The
  header lines record the MACS2 command and settings.

**Filtered peaks.** The same peaks after the ENCODE post-processing
steps (see `methods.pdf`): at most 300,000 peaks per sample, the most
significant by p-value, clipped to the chromosome ends, with every
peak that overlaps the ENCODE blacklist removed and the score capped
at 1000. Peaks keep their MACS2 names, so a peak can be looked up in
the three files above.

- `*_peaks.blacklistFiltered.narrowPeak`: the filtered peaks, same
  columns as `*_peaks.narrowPeak`, sorted by p-value.
- `*_summits.blacklistFiltered.bed`: the summits of the filtered
  peaks.

Use the `blacklistFiltered` files for any per-sample peak analysis.
They are the peaks the atlas was built from. The per-sample peak
counts in these files are the ones to read against ENCODE's data
standards.

## Peak Atlas and Counts (`atlas`)

### Atlas (`macsPeaksMerged.saf`)

The peak atlas: the filtered peaks of all samples pooled and merged,
with peaks within 500 bp of each other joined into one region. The
atlas is one set of regions shared by every sample, so that the
samples can be compared region by region. It is in SAF format, as read
by featureCounts, with no header line and five tab-separated columns:

| Column | Contents |
| --- | --- |
| `GeneID` | region name, `Peak_1`, `Peak_2`, ... in genome order |
| `Chr` | chromosome |
| `Start` | start, 1-based |
| `End` | end, inclusive |
| `Strand` | always `+` |

Note that SAF coordinates are 1-based and inclusive, unlike the BED
and narrowPeak files. To convert a region to BED, subtract 1 from
`Start`.

Because the atlas pools and merges the peaks of all samples, its
number of regions grows with the number of samples and is not a
per-sample peak count.

### Raw counts (`peaks_raw_fcCounts.txt`)

The number of fragments of each sample in each atlas region, from
featureCounts. The first line is a comment with the featureCounts
command; the second is the header. The columns are the five atlas
columns (`Geneid`, `Chr`, `Start`, `End`, `Strand`), the region
`Length`, and then one column per sample, named by the path of its
filtered BAM (`out/<sample>/<sample>_postProcess.bam`).

A fragment is a read pair; it is counted once in every region it
overlaps by at least one base, so a fragment that spans two adjacent
regions is counted in both. The counts are raw, not normalized. They
are the input to the DESeq2 scale factors, the QC plots and the
differential analysis, and are what to use for any analysis that
needs counts.

### Count summary (`peaks_raw_fcCounts.txt.summary`)

For each sample, the number of fragments assigned to the atlas
(`Assigned`) and the number not assigned, by reason (the
`Unassigned_*` rows; `Unassigned_NoFeatures` is fragments that
overlap no atlas region). The sum of the column is the number of
fragments in the filtered BAM.

## Signal Tracks (`bigwig`)

`<sample>_postProcess.shifted.10mNorm.bw` is the coverage of the
filtered, Tn5-shifted reads of the sample, scaled to 10 million reads
(`10mNorm`), in bigWig format. The tracks can be loaded directly into
IGV or the UCSC Genome Browser. Because every sample is scaled to the
same number of reads, the tracks can be compared by eye across
samples. The value at a base is the number of shifted reads covering
it, multiplied by 10 million and divided by the sample's total number
of shifted reads.

## Metrics (`metrics`)

### QC report (`Proj_<no>__ATACSeqQC.pdf`)

Six pages built from the count matrix, with samples colored by
group:

1. The percentage of each sample's fragments that fall in the atlas
   (`InPeaks`) against those outside it. This is a signal-to-noise
   measure; a sample with a much lower percentage than the others may
   have failed. It is computed against the pooled atlas, and its
   denominator includes mitochondrial fragments, so it is not ENCODE's
   per-sample FRiP and should not be read against ENCODE's FRiP
   standard.
2. The same as counts, in millions of fragments.
3. Sequencing depth in peaks: each sample's fragments in the atlas,
   its ratio to the median sample, and the share of the peaks used
   for the PCA in which it has fewer than 5 fragments. A sample below
   half the median is marked `LOW DEPTH`. Such a sample has noisy
   counts and may separate from all others on a principal component
   of its own, whatever its biology.
4. A principal component analysis of the samples on the normalized
   log2 counts over the atlas, with sample labels. The size of a
   point is the sample's fragments in peaks. Replicates of a group
   should lie together; samples that separate from their group on PC1
   or PC2 are worth a look.
5. The same PCA without labels.
6. A depth check of the PCA: beside the PCA of page 4, the same PCA
   after every sample is reduced, by drawing its fragments at random,
   to about the depth of the smallest. A sample that is apart from its
   group on the left but not on the right is apart because of its
   depth; a sample apart on both differs from its group in its data.
   The reduced counts are used for this page only.

### Depth table (`Proj_<no>__ATACSeqDepth.csv`)

The numbers of page 3, one row per sample: `SampleID`, `Group`,
`FragmentsInPeaks`, `RatioToMedian`, `FracKeptPeaksLT5` (the share of
the peaks used for the PCA with fewer than 5 fragments) and
`LowDepth` (`TRUE` below half the median).

### Insert-size distributions (`Proj_<no>__postInsDistribution.pdf`)

The insert-size (fragment-length) distribution of each sample's
filtered reads, from Picard CollectInsertSizeMetrics:

1. All samples on one plot, density against insert size (0 to
   1000 bp), colored by group with one line type per replicate.
2. The same on a log scale, which shows the nucleosomal peaks.
3. One panel per sample, six per page.

A successful ATAC-seq library shows a large nucleosome-free peak below
about 100 bp and smaller peaks at about 200 bp (one nucleosome) and
400 bp (two nucleosomes).

The histogram covers 0 to 1000 bp. Its smallest insert is 31 bp:
pairs with an insert of 30 bp or less are removed by the read filter
(see the methods document).

The per-sample files behind the plot are `<sample>_postProcess___INS.txt`,
the Picard metrics file (summary statistics such as the median and
mode insert size, followed by the histogram as `insert_size` and
count columns), and `<sample>_postProcess___INS.pdf`, Picard's own
histogram plot.

### Mitochondrial fraction (`<sample>_postProcess___MITO.txt`)

A header line and one line for the sample, with the columns `SAMPLE`,
`MITO_CONTIG`, `READS`
(all reads of the filtered BAM), `MITO_READS` (those on the
mitochondrial contig) and `FRAC_MITO`. Mitochondrial reads are left
out of the peaks, the atlas and the signal tracks. The fraction is
computed after the read filter, which removes duplicates, so it is
lower than a fraction computed on the raw alignments (as ENCODE
reports it) and should not be read against that standard.

### TSS enrichment (`<sample>_postProcess.tss_enrich.csv`)

The ENCODE TSS enrichment score, the standard measure of
signal-to-noise for ATAC-seq: how much higher the signal is at
transcription start sites than in the flanking regions 2 kb away.
Columns:

| Column | Contents |
| --- | --- |
| `sample_id` | sample |
| `tss_enrich` | the TSS enrichment score |
| `peak_pos_bp` | position of the maximum of the profile, relative to the TSS |
| `avg_noise` | mean signal in the flanking regions, the denominator of the score |
| `n_tss_used`, `n_tss_total`, `pct_tss_used` | TSSs scored, TSSs in the reference set, and their ratio; `pct_tss_used` should be close to 1 |
| `n_mapped_reads` | reads of the filtered BAM |
| `read_len` | read length used for the computation |
| `bam_file` | the filtered BAM |

ENCODE's reading of the score depends on the genome and on the TSS
annotation. The scores here use ENCODE's GENCODE TSS sets (v24 for
hg38, v19 for b37, vM7 for mm10); read the cutoffs as a guide:

| Genome | Concerning | Acceptable | Ideal |
| --- | --- | --- | --- |
| hg38 | < 5 | 5 to 7 | > 7 |
| b37 | < 6 | 6 to 10 | > 10 |
| mm10 | < 10 | 10 to 15 | > 15 |

Two plots accompany the score: `<sample>_postProcess.tss_enrich.pdf`,
the mean signal profile across the 4 kb window centered on the TSS,
and `<sample>_postProcess.large_tss_enrich.png`, a heatmap of the
signal at every TSS, sorted by signal, with the mean profile above it.

## Differential Accessibility

For projects with replicates, pairwise comparisons between groups are
delivered in the `diff` folder as two files, named with the project
number and, when one was used, a run tag:

- `Proj_<no>__DiffPeaksEdgeR_V4.xlsx`: one `Summary` sheet with the
  number of significant regions per comparison, then one sheet per
  comparison, named `<group2>_v_<group1>`, listing the regions
  significant at FDR < 0.05. Columns: `PeakNo` (the atlas region),
  `logFC` (log2 fold change, group 2 over group 1: positive means more
  accessible in group 2), `logCPM` (mean log2 counts per million),
  `PValue`, `FDR` (Benjamini-Hochberg adjusted p-value), the region's
  `Chr`, `Start`, `End` (1-based, from the atlas), `Strand` and
  `Length`, and the ChIPseeker annotation: `SYMBOL` and `GENENAME` of
  the nearest gene, `annotation` (promoter, exon, intron, distal
  intergenic, ...) and `distanceToTSS`.
- `Proj_<no>__DiffPeaks_V4.pdf`: the PCA of the samples (with and
  without labels; the point size is the sample's fragments in peaks,
  as in the QC report), then for each comparison an MA plot (log fold
  change against mean log CPM) and a volcano plot (log fold change
  against p-value, on a log scale). Regions significant at FDR < 0.05
  are colored by the group they are more accessible in, with the
  count of each in the legend; the dashed line on the volcano plot is
  the p-value at which the FDR cut falls.

The names above are those without a run tag (hence the double
underscore); with a tag `<tag>` they are
`Proj_<no>_<tag>_DiffPeaksEdgeR_V4.xlsx` and
`Proj_<no>_<tag>_DiffPeaks_V4.pdf`.

Samples whose group starts with `EXC` were excluded from the
differential analysis.

In our deliveries made before October 2026 the `FDR` column was a
Holm correction, which is much stricter; the same data now give more
regions at FDR < 0.05.
