# ATAC-Seq pipeline

## Version 2.0.0 - 2026-10-12

Most outputs differ from those of v1.x; do not mix the two. See
`CHANGELOG.md` for what changed and why.

Input: paired-end BAMs aligned with PEMapper, duplicates marked. Peak
calling and bigWigs use each read of a pair as a single read; the
count matrix counts fragments. BigWigs and the DESeq2 scale factors
follow R. Koche's methods (`NOTES.md`); read filtering and peak calling
follow the ENCODE ATAC-seq pipeline, except for the MAPQ cutoff, where
duplicates are marked, and two read filters that ENCODE does not have
(see "Post-alignment filtering" below and `CHANGELOG.md`). Open work
is in `00.FIXME.md`.

MACS2 must be installed in a local venv first. See below.

Runs on IRIS/Slurm. From the analysis directory, with `~/bin/sbatch` on
PATH (it exports `SBATCH_SCRIPT_DIR`, which `pipe.sh` needs to find its
own checkout; the partition, qos, time and memory are `#SBATCH` directives
in `pipe.sh`):

```
sbatch /path/to/ATAC-seq/pipe.sh [-q MAPQ] BAM1 [BAM2 ...]
```

To find out whether the run worked, at any time, from the same directory:

```
/path/to/ATAC-seq/bin/checkRun.sh
```

It prints `ATAC STATUS: OK`, `FAILED` or `RUNNING` with the stage, the
failed jobs and their logs, and exits 0, 1 or 2. If `sacct` itself fails
it prints `UNKNOWN` and exits 3; run it again later. The raw status is in
`00.RUNSTATUS.txt`.

One run per analysis directory. To rerun, use a new directory: a run
submitted where `00.RUNSTATUS.txt`, `SLURM.CTRL/jobs.tsv`, `out/`,
`callpeaks/` or `atacSeq/` already exists stops at once and changes
nothing; the reason is in its log, `SLURM.CTRL/<jobid>.out`.

The genome is detected from the BAM headers and recorded as `GENOME=` in
`00.RUNSTATUS.txt`. All BAMs must be on the same build, and the build
must be one of `b37` (GRCh37, contigs `1`, `2`, ...), `hg38` (GRCh38,
`chr1`, ...; called `b38` before v2.0.0) or `mm10`; otherwise the run
stops before submitting any job, with the reason in `00.RUNSTATUS.txt`.
See `docs/GENOMES.md`. hg38 and b37 have been tested on this version;
no mm10 run has been made yet.

Peaks, bigWigs and the atlas use reads on the autosomes, X and Y only.
Mitochondrial reads are left out for every build, as in the ENCODE
ATAC-seq pipeline (before v2.0.0 they were kept for b37 and mm10); each
sample's mitochondrial fraction is in `atacSeq/metrics/*___MITO.txt`.

- Post-alignment filtering:

    - Proper pairs only; reads marked as duplicates, secondary and
      QC-fail reads removed. The pipeline does not mark duplicates: the
      input BAMs must already be marked (PEMapper runs Picard
      MarkDuplicates; BWA-MEM marks nothing)
    - MAPQ 10+ (`-q MAPQ`; the count matrix uses the same cutoff)
    - A read whose mate was removed is removed too (`samtools fixmate
      -r`, as in the ENCODE ATAC-seq pipeline; before v2.0.0 a read
      whose mate failed the MAPQ cutoff was kept)
    - Reads whose 5' end is soft-clipped removed, with their mates. A
      Tn5 read starts at the cut site, so its 5' end must align; BWA-MEM
      aligns a short seed of a read that is not from the reference and
      clips the rest, and such reads pile up on single bases. ENCODE
      has no such filter: its end-to-end aligner (bowtie2) leaves these
      reads unaligned
    - Pairs with an insert of 30 bp or less removed. ENCODE has no such
      filter; in some libraries these fragments are half of all
      fragments and almost all lie in peaks, so they would dominate the
      counts

- Count matrix: featureCounts v2.1.1 counts fragments over the merged
  peak atlas. Two checks must pass before the run continues (the
  fragment total against the mapped reads, and an independent samtools
  count of the fragments in peaks); the results are in
  `peaks_raw_fcCounts.checks.txt`.

- BigWig: coverage of the Tn5-shifted reads, with no read extension,
  scaled to 10 million reads (R. Koche's method, `NOTES.md`). The
  scale uses the reads of the shifted BED, so the mitochondrial reads
  are not in it.

- Peak calling: MACS2 on the Tn5-shifted reads with the ENCODE ATAC-seq
  settings (`--nomodel --shift -75 --extsize 150 --keep-dup all
  --call-summits -p 0.01`). Each sample's directory in `atacSeq/macs`
  has MACS2's own files, unchanged (`*_peaks.narrowPeak`,
  `*_summits.bed`, `*_peaks.xls`), and two files made from them by
  ENCODE's steps: the 300,000 peaks with the smallest p-values, clipped
  to the contig ends, without the peaks that overlap the ENCODE
  blacklist of the genome (`lib/genomes/<TAG>.blacklist.bed`):
  `*_peaks.blacklistFiltered.narrowPeak` (sorted by p-value, score
  capped at 1000) and `*_summits.blacklistFiltered.bed`. The atlas and
  the counts are built from the `blacklistFiltered` peaks, and the
  blacklist is subtracted from the merged atlas too, so no count comes
  from a blacklisted region.

- Atlas: the `blacklistFiltered` peaks of all samples merged within
  500 bp (R. Koche's method; ENCODE does not merge), so ENCODE's
  peak-count standards apply to each sample's peaks, not to the atlas.

- QC: insert-size distribution (Picard), mitochondrial fraction, TSS
  enrichment (the ENCODE algorithm, `R/TSSEnrich/README.md`), and a
  QC PDF with each sample's depth in peaks and the PCA of the
  samples, also after thinning every sample to the smallest depth.

## Output

Delivered in `atacSeq/` in the analysis directory:

| Directory | Contents |
| --- | --- |
| `atlas/` | `macsPeaksMerged.saf` (the atlas; SAF, 1-based), `peaks_raw_fcCounts.txt` (raw counts) and its `.summary` |
| `bigwig/` | `<SID>_postProcess.shifted.10mNorm.bw` |
| `macs/<SID>_postProcess.shifted/` | MACS2's files and the two `blacklistFiltered` files |
| `metrics/` | per sample `___INS.txt`/`.pdf`, `___MITO.txt`, `tss_enrich.csv`/`.pdf`, `large_tss_enrich.png`; for the run `*__postInsDistribution.pdf`, `*__ATACSeqQC.pdf` and `*__ATACSeqDepth.csv` |

The count-check report `peaks_raw_fcCounts.checks.txt`, the DESeq2
scale factors `scaleFactorsDESeq2.csv` (inverted, `NOTES.md`) and
`sampleManifest.csv` stay in the analysis directory.

## Differential analysis

- If there are replicates you can do a differential peak analysis using: `R/diffAnalysisPairwise.R`. Run it from the analysis directory; the genome is read from `00.RUNSTATUS.txt`. Usage:
```
Rscript /path/to/ATAC-seq/R/diffAnalysisPairwise.R SampleManifest.csv Comparisons.csv [RUNTAG]
```

Example inputs:

```
  SampleManifest.csv

        SampleID,Group,MapID
        aNSC_loxp15_1,aNSC_loxp15,s_aNSC_loxp15_1
        aNSC_loxp15_2,aNSC_loxp15,s_aNSC_loxp15_2
        aNSC_p53_1,aNSC_p53,s_aNSC_p53_1

   Comparisons.csv

        aNSC_loxp15,aNSC_p53

   Sign convention X2-X1; e.g., aNSC_p53-aNSC_loxp15
```

Each comparison is an edgeR likelihood ratio test on the samples of its
two groups. A peak is significant at FDR < 0.05 (`p.adjust(method =
"fdr")`, i.e. Benjamini-Hochberg). Before v2.0.0 the `FDR` column was
Holm's correction, which is much stricter. Samples whose `Group` starts
with `EXC` are left out. Every count matrix column must have a manifest
row (`MapID` is the BAM's SM tag), and every row a column.

## MACS2, IDR Installation

MACS2 2.2.9.1 and IDR 2.0.3 both build under the python 3.10 on the IRIS
PATH, so the venv uses whatever `python3` is current.

In root of ATAC-seq repo

```{base}
. 00.SETUP.sh
```

