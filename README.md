# ATAC-Seq pipeline

## Version 2.0.0 - 2026-10-13

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
failed jobs and their logs, and exits 0, 1 or 2 (`NONE` and exit 1 if
the directory holds no run). If `sacct` itself fails
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
sample's mitochondrial fraction, of the reads that pass the filters
below, is in `atacSeq/metrics/*___MITO.txt`.

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
    - Pairs with an insert of 30 bp or less removed. ENCODE has no
      such filter. They are almost all BWA-MEM seed alignments of the
      kind the previous filter removes, which takes more than 99% of
      them first. Without either filter they were half of all
      fragments in some libraries, and almost all lay in peaks

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

`docs/output.md` describes every delivered file for the client and
`docs/methods.md` the methods; both are rendered to HTML and PDF by
`docs/MAKE_DOCS.sh`, and both forms are delivered in `atacSeq/docs`.

## Differential analysis

With replicates, `R/diffAnalysisPairwise.R` tests pairs of groups
for differential peaks. Run it from the analysis directory; the
genome is read from `00.RUNSTATUS.txt`:

```
Rscript /path/to/ATAC-seq/R/diffAnalysisPairwise.R \
    sampleManifest.csv Comparisons.csv [RUNTAG]
```

`sampleManifest.csv` has a header and the columns `MapID`,
`SampleID` and `Group`. `pipe.sh` writes one if none exists, with
`SampleID` and `Group` guessed from the SM tag; check and edit it
first. `Comparisons.csv` has no header; each row is one comparison,
`Group1,Group2`:

```
sampleManifest.csv

    MapID,SampleID,Group
    s_aNSC_loxp15_1,aNSC_loxp15_1,aNSC_loxp15
    s_aNSC_loxp15_2,aNSC_loxp15_2,aNSC_loxp15
    s_aNSC_p53_1,aNSC_p53_1,aNSC_p53

Comparisons.csv

    aNSC_loxp15,aNSC_p53
```

The sign convention is Group2 - Group1: a positive logFC means
higher in the second group (here aNSC_p53).

Each comparison is an edgeR likelihood ratio test on the samples of its
two groups. A peak is significant at FDR < 0.05 (`p.adjust(method =
"fdr")`, i.e. Benjamini-Hochberg). Before v2.0.0 the `FDR` column was
Holm's correction, which is much stricter. Samples whose `Group` starts
with `EXC` are left out. Every count matrix column must have a manifest
row (`MapID` is the BAM's SM tag), and every row a column.

The results are `<PROJ>_<RUNTAG>_DiffPeaksEdgeR_V5.xlsx` (a summary
sheet and, per comparison, the significant peaks with their nearest
gene) and `<PROJ>_<RUNTAG>_DiffPeaks_V5.pdf` (PCA, MA and volcano
plots). `<PROJ>` is the last `Proj_*` or `B-<digits>` directory in
the path of the analysis directory, else `ATAC`.

## Delivery

From the analysis directory, after `bin/checkRun.sh` exits 0, and
after the differential analysis if there is one:

```bash
/path/to/ATAC-seq/deliverResults.sh -n -d     # dry run: print, write nothing
/path/to/ATAC-seq/deliverResults.sh -d
```

`-d` delivers to `/data1/core002/res/bic/results/<pi>/<invest>/Proj_<no>`,
taken from the path `.../Users/<Pi>/<Invest>/Proj_<no>/...`; any other
delivery folder can be given instead. The run folder is the next
`r_NNN` there, unless the path ends in one. The script stops unless
the run `COMPLETED` and you can read the delivery folder, then:

- writes `project.yaml` (`bin/readme2yaml.R`) from the run's
  `GENOME` and the PI and investigator emails of the first
  `README.txt` found going up from the analysis directory to
  `Proj_<no>`;
- copies, as `svc_core001_bot01`, `atacSeq/` to `r_NNN/atacSeq`,
  `docs/output` and `docs/methods` of this checkout (PDF and HTML)
  to `atacSeq/docs`, and the project's `*_DiffPeaksEdgeR_V5.xlsx` and
  `*_DiffPeaks_V5.pdf` (any RUNTAG) to `atacSeq/diff`;
- writes the email, from `tmpldeliveryEmail_*.txt`, to
  `deliveryEmail_<no>_r_NNN_<yymmdd>.txt`;
- registers the run with the bicdelivery permissions database
  (`~/Code/BIC/Delivery/Version2j`). If this fails, the files are
  already delivered; the script prints the command to rerun this
  step alone.

The email gives the release in `VERSION.md` of the checkout the
script is run from (see "Releases") and links to the docs at that tag.
Diff files with a RUNTAG (the email names the untagged files) give a
`WARNING`: edit the email before sending it. A delivery folder other
than `<pi>/<invest>/Proj_<no>` for the PI and investigator of
`README.txt`, the users the database is given, also gives a
`WARNING`.

## Releases

The client-facing version is the release in `VERSION.md`
(`- **Version:** vX.Y.Z`). It is set by hand when a release is made,
never derived from git, and a git tag of the same name marks the
release, so the doc links in the email resolve. `pipe.sh` records it
as `RELEASE=` in `00.RUNSTATUS.txt`, beside `VERSION=`, the
`git describe` of the exact code, which is for us, not the client.

## Installing MACS2

Once per checkout, from the root of the repository; the script must
be sourced, not run:

```bash
cd /path/to/ATAC-seq
. 00.SETUP.sh
```

It builds `venv/` with the `python3` on PATH (3.10 on IRIS) and
installs MACS2 2.2.9.1 and IDR 2.0.3. The pipeline runs only MACS2;
IDR is installed but never run (`00.FIXME.md` #29). `pipe.sh` stops
if `venv/` is missing.

