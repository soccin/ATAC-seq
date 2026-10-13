# ATAC-Seq pipeline

## Version 1.5.1 - 2026-10-10

Single end version which uses both reads from PE-runs. Using methods from R.K. for bigWig generation.

Need to install MACS2 locally in venv. See below.

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

    - Proper pairs only; reads marked as duplicates (the marking is done
      upstream, by the mapper), secondary and QC-fail reads removed
    - MAPQ 10+ (`-q MAPQ`; the count matrix uses the same cutoff)
    - A read whose mate was removed is removed too (`samtools fixmate
      -r`, as in the ENCODE ATAC-seq pipeline); there is no filter on
      insert size (before v2.0.0 pairs with an insert of 30 bp or less
      were removed, and a read whose mate failed the MAPQ cutoff was
      kept)

- Count matrix: featureCounts v2.1.1 counts fragments over the merged
  peak atlas. Two checks must pass before the run continues (the
  fragment total against the mapped reads, and an independent samtools
  count of the fragments in peaks); the results are in
  `peaks_raw_fcCounts.checks.txt`.

- BigWig file generation.

	- convert bam to bed
	- extend reads to library size in the 3’ direction for ChIP, no read extension for ATAC
	- compute density for bigwig formation
		- normalizing to 10 million mapped reads

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


## Differential analysis

- If there are replicates you can do a differential peak analysis using: `R/diffAnalysisPairwise.R`. Run it from the analysis directory; the genome is read from `00.RUNSTATUS.txt`. Usage:
```
Rscript R/diffAnalysisPairwise.R SampleManifest.csv Comparisons.csv [RUNTAG]
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

