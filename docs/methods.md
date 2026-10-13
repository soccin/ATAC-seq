# ATAC-seq Pipeline

### Version 2.0.0

## Overview

Samples were processed with the BIC ATAC-seq pipeline, which takes
aligned, duplicate-marked BAM files and produces filtered alignments,
normalized signal tracks (bigWig), MACS2 peak calls for each sample, a
merged peak atlas with a raw fragment count matrix, QC metrics and,
when the design has replicates, a differential accessibility analysis.
The code is at
[https://github.com/soccin/ATAC-seq/releases/tag/v2.0.0](https://github.com/soccin/ATAC-seq/releases/tag/v2.0.0).

Read filtering, the Tn5 shift and peak calling follow the ENCODE
ATAC-seq pipeline (ENCODE-DCC/atac-seq-pipeline v2.2.3). The signal
tracks, the peak atlas and the count normalization follow the methods
of R. Koche (MSKCC Center for Epigenetics Research).

Reads were aligned to one of three reference genomes: GRCh37 (b37),
GRCh38 (the GDC `GRCh38.d1.vd1` build, hg38) or GRCm38 (mm10). The
build of a project is recorded in the BAM headers. Peaks, signal tracks
and the atlas use the autosomes, X and Y only.

## Methods

**Alignment.** Paired-end reads were aligned with BWA-MEM 0.7.19 (Li
2013) using the PEMapper pipeline
([https://github.com/soccin/PEMapper](https://github.com/soccin/PEMapper)),
which trims adapters with cutadapt, aligns with `bwa mem -M`, merges
lanes, sorts, and marks duplicates with Picard MarkDuplicates 2.25.5.
Duplicates are marked in the BAM, not removed. The PEMapper version is
in the `@PG` header line of each BAM.

**Read filtering.** With samtools 1.20 (Danecek 2021) each BAM was
reduced to properly paired reads with mapping quality 10 or more,
without unmapped reads, reads whose mate is unmapped, secondary
alignments, reads failing vendor QC or duplicates (`-q 10 -f 3 -F
1804`). The reads were then name-sorted and passed through `samtools
fixmate -r`, and the same flag filter was applied again, which removes
a read whose mate did not pass the first filter. This is the filter of
the ENCODE ATAC-seq pipeline. There is no filter on insert size. The
result was coordinate-sorted with Picard SortSam 3.4.0; this filtered
BAM is the input to every later step.

**Tn5 shift.** Each read was converted to a BED interval with bedtools
2.31.1 (Quinlan 2010) and shifted to place its 5' end at the Tn5
insertion site: +4 bp on the forward strand and -5 bp on the reverse
strand (Buenrostro 2013). Reads on the mitochondrial contig and on
unplaced, random, decoy and viral contigs were dropped at this step, as
in the ENCODE pipeline; the fraction of mitochondrial reads is reported
as a QC metric. This set of shifted reads is the input to the signal
tracks and to peak calling.

**Signal tracks.** For each sample the coverage of the shifted reads
was computed with bedtools genomecov, each read as its own interval
with no extension, scaled to 10 million reads (the scale factor is 10
million divided by the number of shifted reads of the sample), and
written as a bigWig with the UCSC wigToBigWig tool (Kent 2010). Tracks
are therefore comparable between samples by sequencing depth; the
DESeq2 scale factors below are not applied to them.

**Peak calling.** Peaks were called for each sample with MACS2 2.2.9.1
(Zhang 2008) on the shifted reads, with the settings of the ENCODE
ATAC-seq pipeline: `--nomodel --shift -75 --extsize 150 --keep-dup all
--call-summits -p 0.01` and the genome size `hs` or `mm`. Each read
contributes a 150 bp window centered on its insertion site, every read
is kept (duplicates were already removed by the filter), and peaks are
called at a p-value of 0.01, a permissive cutoff chosen by ENCODE so
that peaks can be ranked downstream. MACS2's own output (narrowPeak,
summits, xls) is delivered unchanged. The peaks were then
post-processed as ENCODE does: sorted by p-value, the 300,000 most
significant kept, clipped to the contig ends, peaks overlapping the
ENCODE blacklist of the genome removed (Amemiya 2019; Boyle lab
blacklist v2 for hg38 and mm10, the DAC consensus excludable list for
b37) and the score column capped at 1000. These `blacklistFiltered`
peaks and summits are delivered beside MACS2's files and are the peaks
used for the atlas.

**Peak atlas.** The blacklist-filtered peaks of all samples were pooled
and merged with bedtools merge, joining peaks within 500 bp of each
other into one region (R. Koche's method; ENCODE does not merge
peaks). Any blacklisted region bridged by a merge was subtracted, so no
atlas region overlaps the blacklist. The result is the peak atlas, one
set of regions shared by every sample of the project, written in SAF
format (1-based, inclusive coordinates) with the regions named
`Peak_1`, `Peak_2`, ... in genome order. Because peaks are merged and
pooled over samples, the number of atlas regions grows with the number
of samples and is not comparable to ENCODE's per-sample peak-count
standards; the per-sample `narrowPeak` files are.

**Counts.** Fragments were counted over the atlas with featureCounts
2.1.1 (Liao 2014) from the filtered BAMs, counting read pairs as single
fragments (`-p --countReadPairs`), with a fragment that overlaps two
atlas regions counted in both (`-O`), and with the same mapping quality
cutoff as the read filter (`-Q 10`). The counts were verified by two
independent checks before use: each sample's fragment total must be
half its mapped reads, within 5%, and the number of fragments in the atlas,
counted a second time with samtools alone, must equal featureCounts'
assigned total exactly. A run whose counts fail either check is
stopped.

**Normalization.** Size factors were computed from the raw count
matrix with DESeq2 (Love 2014; median-of-ratios method, Anders 2010).
They are reported as scale factors, the reciprocal of the DESeq2 size
factor, so that a normalized count is the raw count multiplied by the
sample's scale factor.

**QC.** The insert-size distribution of the filtered BAM was computed
with Picard CollectInsertSizeMetrics and plotted for all samples, to
show the nucleosome-free and nucleosomal fragment populations. The
mitochondrial fraction is the number of reads on the mitochondrial
contig divided by all reads of the filtered BAM. TSS enrichment was
computed with an R port of ENCODE's `encode_task_tss_enrich.py`: the
coverage of the shifted reads in a window of 2 kb on each side of
every unique protein-coding TSS (ENCODE's TSS sets: GENCODE v24 for
hg38, v19 for b37, vM7 for mm10), averaged over TSSs in 400 bins and
divided by the mean coverage of the outer 100 bp on each side; the
score is the maximum of this profile. ENCODE's reading of the score is
given with the output. The fraction of fragments in the atlas and a
principal component analysis of the samples (edgeR TMM-normalized
log2 CPM over the atlas, after filterByExpr) were plotted from the
count matrix.

**Differential accessibility.** When the design has replicates,
pairwise comparisons between groups were tested with edgeR (Robinson
2010; McCarthy 2012). For each comparison, the samples of its two
groups were taken from the count matrix, atlas regions with too few
counts removed with filterByExpr, library sizes normalized with TMM,
dispersions estimated, and a negative binomial generalized linear
model fitted and tested with the likelihood ratio test. P-values were
adjusted with the Benjamini-Hochberg method (Benjamini 1995) and
regions with FDR < 0.05 were called significant. The sign of the log
fold change is group 2 minus group 1 of the comparison: a positive
value means more accessible in group 2. Significant regions were
annotated with the nearest gene and genomic feature using ChIPseeker
(Yu 2015) with the UCSC knownGene annotation of the build and a TSS
region of 5 kb on each side.

## Tool versions

| Tool | Version |
| --- | --- |
| BWA-MEM (PEMapper) | 0.7.19 |
| Picard MarkDuplicates (PEMapper) | 2.25.5 |
| samtools | 1.20 |
| Picard SortSam, CollectInsertSizeMetrics | 3.4.0 |
| bedtools | 2.31.1 |
| UCSC wigToBigWig, bedClip | v369 |
| MACS2 | 2.2.9.1 |
| featureCounts (subread) | 2.1.1 |
| R | 4.5.1 |
| ENCODE ATAC-seq pipeline (settings followed) | 2.2.3 |

<!--
DRAFT NOTES (remove before delivery)

- The Picard, bedtools and wigToBigWig versions are the ones on IRIS
  on 2026-09/10 (docs/SLURM_PORT.md, PEMapper VERSIONS.md). The
  pipeline does not yet record them per run (00.FIXME.md #24). The
  wigToBigWig build is assumed to be v369 like bedClip; check.
- DESeq2, edgeR and ChIPseeker versions are not recorded; add them
  from sessionInfo() of the run, or drop the row.
- The "fraction of fragments in the atlas" is not ENCODE's FRiP
  (00.FIXME.md #6); worded accordingly.
- The mitochondrial fraction is after filtering, unlike ENCODE's
  (00.FIXME.md #7); worded accordingly.
- The TSS sets are ENCODE's v1 GENCODE files, not the v3 files the
  current ENCODE pipeline uses (00.FIXME.md #23); the annotation
  versions are given so the thresholds can be read as a guide.
-->

## References

Amemiya HM, Kundaje A, Boyle AP. The ENCODE Blacklist: identification
of problematic regions of the genome. Sci Rep 2019;9:9354.
doi:10.1038/s41598-019-45839-z

Anders S, Huber W. Differential expression analysis for sequence count
data. Genome Biol 2010;11:R106. doi:10.1186/gb-2010-11-10-r106

Benjamini Y, Hochberg Y. Controlling the false discovery rate: a
practical and powerful approach to multiple testing. J R Stat Soc B
1995;57:289-300.

Buenrostro JD, Giresi PG, Zaba LC, Chang HY, Greenleaf WJ. Transposition
of native chromatin for fast and sensitive epigenomic profiling of open
chromatin, DNA-binding proteins and nucleosome position. Nat Methods
2013;10:1213-1218. doi:10.1038/nmeth.2688

Danecek P, et al. Twelve years of SAMtools and BCFtools. GigaScience
2021;10:giab008. doi:10.1093/gigascience/giab008

ENCODE ATAC-seq pipeline. ENCODE Data Coordination Center.
[https://github.com/ENCODE-DCC/atac-seq-pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline).
Data standards: [https://www.encodeproject.org/atac-seq/](https://www.encodeproject.org/atac-seq/)

Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. BigWig and
BigBed: enabling browsing of large distributed datasets. Bioinformatics
2010;26:2204-2207. doi:10.1093/bioinformatics/btq351

Li H. Aligning sequence reads, clone sequences and assembly contigs with
BWA-MEM. arXiv 2013;1303.3997.

Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose
program for assigning sequence reads to genomic features.
Bioinformatics 2014;30:923-930. doi:10.1093/bioinformatics/btt656

Love MI, Huber W, Anders S. Moderated estimation of fold change and
dispersion for RNA-seq data with DESeq2. Genome Biol 2014;15:550.
doi:10.1186/s13059-014-0550-8

McCarthy DJ, Chen Y, Smyth GK. Differential expression analysis of
multifactor RNA-Seq experiments with respect to biological variation.
Nucleic Acids Res 2012;40:4288-4297. doi:10.1093/nar/gks042

Picard Toolkit. Broad Institute.
[http://broadinstitute.github.io/picard/](http://broadinstitute.github.io/picard/)

Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for
comparing genomic features. Bioinformatics 2010;26:841-842.
doi:10.1093/bioinformatics/btq033

Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for
differential expression analysis of digital gene expression data.
Bioinformatics 2010;26:139-140. doi:10.1093/bioinformatics/btp616

Yu G, Wang LG, He QY. ChIPseeker: an R/Bioconductor package for ChIP
peak annotation, comparison and visualization. Bioinformatics
2015;31:2382-2383. doi:10.1093/bioinformatics/btv145

Zhang Y, et al. Model-based analysis of ChIP-Seq (MACS). Genome Biol
2008;9:R137. doi:10.1186/gb-2008-9-9-r137
