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
ATAC-seq pipeline (ENCODE-DCC/atac-seq-pipeline v2.2.3) with four
differences, described below: the mapping quality cutoff, where
duplicates are marked (by PEMapper, before the read filter), and two
additional read filters. The signal tracks, the peak atlas and the
count normalization follow the methods of R. Koche (MSKCC Center for
Epigenetics Research).

Reads were aligned to one of three reference genomes: GRCh37 (b37),
GRCh38 (the GDC `GRCh38.d1.vd1` build, hg38) or GRCm38 (mm10). The
build of a project is recorded in the BAM headers. Peaks, signal tracks
and the atlas use the autosomes, X and Y only.

## Methods

**Alignment.** Paired-end reads were aligned with BWA-MEM 0.7.19 (Li
2013) using the PEMapper pipeline
([https://github.com/soccin/PEMapper](https://github.com/soccin/PEMapper)),
which trims the TruSeq adapter (`AGATCGGAAGAGC`) with cutadapt, aligns
with `bwa mem -M`, merges lanes, sorts, and marks duplicates with
Picard MarkDuplicates 2.25.5. The Nextera adapter of ATAC-seq
libraries is not trimmed; BWA-MEM soft-clips it at the 3' end of the
read. Duplicates are marked in the BAM, not removed. The PEMapper
version is in the `@PG` header line of each BAM.

**Read filtering.** With samtools 1.20 (Danecek 2021) each BAM was
reduced to properly paired reads with mapping quality 10 or more,
without unmapped reads, reads whose mate is unmapped, secondary
alignments, reads failing vendor QC or duplicates (`-q 10 -f 3 -F
1804`). The cutoff of 10 is R. Koche's; the ENCODE pipeline by
default keeps reads with up to four alignments, or uses a cutoff of 30
when multimapping reads are not kept. The duplicates removed here were
marked by PEMapper on all alignments, before this filter; the ENCODE
pipeline marks duplicates after its own filter. Two further filters,
which the ENCODE pipeline does not have, correct for an artifact of
BWA-MEM local alignment. A read whose sequence is not from the
reference, or only partly from it, can still be reported as a short
alignment of 19 to 35 bp with the rest of the read soft-clipped. Such
reads carry the same non-reference sequence, so hundreds or thousands
of them pile up on single bases, mostly inside peaks; in some
libraries they are about half of all read pairs. Bowtie2 in end-to-end
mode, the ENCODE aligner, cannot produce these alignments, so ENCODE
needs no such filter. The first filter removes every read whose alignment is
soft-clipped at its 5' end, the Tn5 insertion site: a forward-strand
read whose CIGAR string starts with `S` or a reverse-strand read whose
CIGAR string ends with `S`. The reads were then name-sorted and passed
through `samtools fixmate -r`, and the flag filter was applied again,
as in the ENCODE pipeline, which removes a read whose mate did not pass
the filters above. The second filter was applied at the same step: it
removes every pair with an insert of 30 bp or less (|TLEN| <= 30). The
result was coordinate-sorted with Picard SortSam 3.4.0; this filtered
BAM is the input to every later step.

**Tn5 shift.** Each read was converted to a BED interval with bedtools
2.31.1 (Quinlan 2010) and shifted to place its 5' end at the Tn5
insertion site: +4 bp on the forward strand and -5 bp on the reverse
strand (Buenrostro 2013). Reads on the mitochondrial contig were
dropped at this step, as in the ENCODE pipeline, and so were reads on
unplaced, random, decoy and viral contigs, which ENCODE removes only
from its final peaks; the fraction of mitochondrial reads is reported
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
blacklist v2 for hg38 and mm10; for b37 the hg19 DAC consensus
excludable list that ENCODE uses, with the `chr` prefix removed) and
the score column capped at 1000. These `blacklistFiltered`
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
stopped. The result of both checks for each sample is delivered with
the counts.

**Normalization.** Size factors were computed from the raw count
matrix with DESeq2 (Love 2014; median-of-ratios method, Anders 2010)
and written as scale factors, the reciprocal of the DESeq2 size
factor, so that a normalized count is the raw count multiplied by the
sample's scale factor. They are delivered with the counts but are not
applied to any delivered file.

**QC.** The insert-size distribution of the filtered BAM was computed
with Picard CollectInsertSizeMetrics over a fixed range of 0 to 1000
bp, as in the ENCODE pipeline, and plotted for all samples, to show
the nucleosome-free and nucleosomal fragment populations. Its smallest
insert is 31 bp, because of the insert filter above. The
mitochondrial fraction is the number of reads on the mitochondrial
contig divided by all reads of the filtered BAM, so it is computed
after duplicate removal; ENCODE computes its fraction on all
alignments before filtering, and because mitochondrial reads are
highly duplicated the fraction here is the lower of the two. TSS
enrichment was computed with an R port of ENCODE's
`encode_task_tss_enrich.py` on the filtered BAM: the coverage of the
reads, each moved by half the read length so that it is centered on
its 5' end, in a window of 2 kb on each side of every unique
protein-coding TSS, averaged over TSSs in 400 bins and divided by the
mean coverage of the outer 100 bp on each side; the score is the
maximum of this profile. The TSS sets are the GENCODE v24 (hg38), v19
(b37) and vM7 (mm10) files of ENCODE's first genome data release, not
the newer files the current ENCODE pipeline uses, and ENCODE states
its cutoffs for GRCh38 and mm10 on RefSeq TSSs, so the cutoffs given
with the output are a guide. The
fraction of fragments in peaks is the number of fragments
featureCounts assigned to the atlas divided by all fragments of the
filtered BAM, taken from the featureCounts summary. It measures the
same thing as ENCODE's FRiP, signal against background, but is a
different quantity: the regions are the merged atlas of all samples
rather than the sample's own peaks, and the denominator includes
mitochondrial fragments, so it is not comparable to ENCODE's FRiP
cutoffs. It and a principal component analysis of the samples (edgeR
TMM-normalized log2 CPM over the atlas, after filterByExpr) were
plotted from the count matrix, with each sample's fragments in the
atlas (its depth). To tell an effect of depth from one of the data,
the PCA was repeated after the counts of every sample were thinned
at random (binomial, fixed seed) to about the depth of the smallest.

**Differential accessibility.** When the design has replicates,
pairwise comparisons between groups were tested with edgeR (Robinson
2010; McCarthy 2012). For each comparison, the samples of its two
groups were taken from the count matrix, atlas regions with too few
counts removed with filterByExpr (over all samples, then over the
samples of the comparison), library sizes normalized with TMM,
dispersions estimated, and a negative binomial generalized linear
model fitted and tested with the likelihood ratio test. P-values were
adjusted with the Benjamini-Hochberg method (Benjamini 1995) and
regions with FDR < 0.05 were called significant. The sign of the log
fold change is group 2 minus group 1 of the comparison: a positive
value means more accessible in group 2. Significant regions were
annotated with the nearest gene and genomic feature using ChIPseeker
(Yu 2015) with the UCSC knownGene annotation of the build and a TSS
region of 5 kb on each side. For each comparison, every tested region
was plotted as an MA plot (log fold change against mean log2 CPM) and
a volcano plot (log fold change against p-value, on a log scale),
with the significant regions colored by the group they are more
accessible in.

**Pathway enrichment.** When requested, each comparison was tested
for pathway enrichment by preranked gene set enrichment analysis
(Subramanian 2005) with fgsea (Korotkevich 2021). Every tested region
was assigned to the gene with the nearest TSS with ChIPseeker and the
UCSC knownGene annotation of the build, and regions farther from that
TSS than a set distance (by default 1 kb on either side, the
promoter) were dropped. Each region was scored as the sign of its log
fold change times -log10 of its p-value, and each gene took the score
of its region with the largest absolute score. The genes ranked by
this score were tested against MSigDB (Liberzon 2015) gene sets of 15
to 500 ranked genes, from msigdbr, with fgseaMultilevel (no lower
bound on the p-value, fixed seed); p-values were adjusted with the
Benjamini-Hochberg method within each collection of each comparison.
Human data used the human MSigDB collections, mouse data the mouse
MSigDB (Castanza 2023). A positive normalized enrichment score means
the set's genes are more accessible in group 2.

## Tool versions

| Tool | Version |
| --- | --- |
| BWA-MEM (PEMapper) | 0.7.19 |
| Picard MarkDuplicates (PEMapper) | 2.25.5 |
| samtools | 1.20 |
| Picard SortSam, CollectInsertSizeMetrics | 3.4.0 |
| bedtools | 2.31.1 |
| UCSC bedClip | v369 |
| UCSC wigToBigWig | kent source v361 |
| MACS2 | 2.2.9.1 |
| featureCounts (subread) | 2.1.1 |
| R | 4.5.1 |
| DESeq2 | 1.50.2 |
| edgeR | 4.8.2 |
| ChIPseeker | 1.46.1 |
| TxDb.Hsapiens.UCSC.hg38.knownGene (hg38) | 3.22.0 |
| TxDb.Hsapiens.UCSC.hg19.knownGene (b37) | 3.22.1 |
| TxDb.Mmusculus.UCSC.mm10.knownGene (mm10) | 3.10.0 |
| fgsea | 1.36.2 |
| msigdbr | 26.1.1 |
| MSigDB | 2026.1 (Hs, Mm) |
| ENCODE ATAC-seq pipeline (settings followed) | 2.2.3 |

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

Castanza AS, Recla JM, Eby D, Thorvaldsdottir H, Bult CJ, Mesirov JP.
Extending support for mouse data in the Molecular Signatures Database
(MSigDB). Nat Methods 2023;20:1619-1620. doi:10.1038/s41592-023-02014-7

Danecek P, et al. Twelve years of SAMtools and BCFtools. GigaScience
2021;10:giab008. doi:10.1093/gigascience/giab008

ENCODE ATAC-seq pipeline. ENCODE Data Coordination Center.
[https://github.com/ENCODE-DCC/atac-seq-pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline).
Data standards: [https://www.encodeproject.org/atac-seq/](https://www.encodeproject.org/atac-seq/)

Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. BigWig and
BigBed: enabling browsing of large distributed datasets. Bioinformatics
2010;26:2204-2207. doi:10.1093/bioinformatics/btq351

Korotkevich G, Sukhov V, Budin N, Shpak B, Artyomov MN, Sergushichev A.
Fast gene set enrichment analysis. bioRxiv 2021;060012.
doi:10.1101/060012

Li H. Aligning sequence reads, clone sequences and assembly contigs with
BWA-MEM. arXiv 2013;1303.3997.

Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose
program for assigning sequence reads to genomic features.
Bioinformatics 2014;30:923-930. doi:10.1093/bioinformatics/btt656

Liberzon A, Birger C, Thorvaldsdottir H, Ghandi M, Mesirov JP, Tamayo
P. The Molecular Signatures Database (MSigDB) hallmark gene set
collection. Cell Syst 2015;1:417-425. doi:10.1016/j.cels.2015.12.004

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

Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment
analysis: a knowledge-based approach for interpreting genome-wide
expression profiles. Proc Natl Acad Sci USA 2005;102:15545-15550.
doi:10.1073/pnas.0506580102

Yu G, Wang LG, He QY. ChIPseeker: an R/Bioconductor package for ChIP
peak annotation, comparison and visualization. Bioinformatics
2015;31:2382-2383. doi:10.1093/bioinformatics/btv145

Zhang Y, et al. Model-based analysis of ChIP-Seq (MACS). Genome Biol
2008;9:R137. doi:10.1186/gb-2008-9-9-r137
