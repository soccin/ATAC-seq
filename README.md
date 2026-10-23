# ATAC-Seq pipeline

## Quick Start

Run every command from the analysis directory. Details are in the
sections below.

```
# Run the pipeline (Slurm control job; MACS2 venv must exist)
sbatch /path/to/ATAC-seq/pipe.sh [-q MAPQ] BAM1 [BAM2 ...]

# After editing sampleManifest.csv, redo the QC reports
# (atacSeq/metrics); pipe.sh runs both once
Rscript /path/to/ATAC-seq/plotINSStats.R
Rscript /path/to/ATAC-seq/R/analyzeATAC.R sampleManifest.csv [RUNTAG]

# Differential peaks (needs replicates; writes atacSeq/diff)
Rscript /path/to/ATAC-seq/R/diffAnalysisPairwise.R \
    sampleManifest.csv Comparisons.csv [RUNTAG]

# Pathway GSEA of the diff results (same RUNTAG; writes atacSeq/diff)
Rscript /path/to/ATAC-seq/R/pathwayGSEA.R [-c COLLECTIONS] [-m MODE] [RUNTAG]

# Motif enrichment of the diff results (same RUNTAG; writes
# atacSeq/diff); run as a job, the motif scan uses -n cores, and
# loadTools.sh keeps R's temp files out of /tmp
bsub -n 8 -M 128G -W 6:00:00 -J motifs -o SLURM.MOTIF/ bash -c \
    ". /path/to/ATAC-seq/bin/loadTools.sh; Rscript /path/to/ATAC-seq/R/motifEnrich.R [-w WIDTH] [-m MINPEAKS] [RUNTAG]"
```

## Version 2.1.0 - 2026-10-22

v2.1.0 adds pathway enrichment (GSEA) of the differential analysis
and writes every result to `atacSeq/` directly; the results v2.0.0
made do not change. Most outputs differ from those of v1.x; do not
mix the two. See `CHANGELOG.md` for what changed and why.

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
  `atacSeq/atlas/peaks_raw_fcCounts.checks.txt`.

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

Every result of the run, and of the differential analysis, is written
to `atacSeq/` in the analysis directory, which is what is delivered:

| Directory | Contents |
| --- | --- |
| `atlas/` | `macsPeaksMerged.saf` (the atlas; SAF, 1-based), `peaks_raw_fcCounts.txt` (raw counts) and its `.summary`, the count checks `peaks_raw_fcCounts.checks.txt`, the DESeq2 scale factors `scaleFactorsDESeq2.csv` (inverted, `NOTES.md`) |
| `bigwig/` | `<SID>_postProcess.shifted.10mNorm.bw` |
| `macs/<SID>_postProcess.shifted/` | MACS2's files and the two `blacklistFiltered` files |
| `metrics/` | per sample `___INS.txt`/`.pdf`, `___MITO.txt`, `tss_enrich.csv`/`.pdf`, `large_tss_enrich.png`; for the run `*__postInsDistribution.pdf`, `*__ATACSeqQC.pdf` and `*__ATACSeqDepth.csv` |
| `diff/` | the differential analysis and pathway enrichment, if run (below); the `.rds` files are not delivered |

The rest of the analysis directory is work and run records:
`out/<SID>/` and `callpeaks/` (per-sample work; the delivered files
are copied from them into `atacSeq/macs` and `atacSeq/metrics`),
`SLURM.*/` (logs), `00.RUNSTATUS.txt`, and `sampleManifest.csv`, which
you edit (below).

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

It reads the counts from `atacSeq/atlas` and writes to `atacSeq/diff`.
The results are `<PROJ>_<RUNTAG>_DiffPeaksEdgeR_V5.xlsx` (a summary
sheet and, per comparison, the significant peaks with their nearest
gene) and `<PROJ>_<RUNTAG>_DiffPeaks_V5.pdf` (PCA, MA and volcano
plots). `<PROJ>` is the last `Proj_*` or `B-<digits>` directory in
the path of the analysis directory, else `ATAC`. It also writes
`<PROJ>_<RUNTAG>_DiffPeaksEdgeR_V5.rds`, every tested peak of every
comparison with the logCPM of its samples, the input of the pathway
enrichment; it is not delivered.

### Pathway enrichment (GSEA)

`R/pathwayGSEA.R` runs preranked GSEA (fgsea) of each comparison
against MSigDB gene sets (msigdbr). Run it from the analysis
directory after `diffAnalysisPairwise.R`, with the same RUNTAG:

```
Rscript /path/to/ATAC-seq/R/pathwayGSEA.R [options] [RUNTAG]

  -c, --collections   default MSIGDB_SETS of GENOMES.tsv:
                      H,C2:CP,-C2:CP:KEGG_LEGACY,C5:GO:BP (human),
                      MH,M2:CP,M5:GO:BP (mouse)
  -m, --mode          promoter (default) or distance
  -w, --tss-window    promoter: every gene with a TSS within w (1000)
  -d, --max-distance  distance: nearest-TSS gene(s) within d (50000)
      --min-size, --max-size   gene set size (15, 500)
  -o, --out-tag       added to the output names
```

Each atlas peak is assigned to genes by their TSSs (every transcript
start of the knownGene TxDb, with GenomicRanges): in promoter mode to
every gene with a TSS within `w`, so a promoter shared by two genes
counts for both, in distance mode to the gene with the nearest TSS
(all of them on a tie) within `d`. Each gene takes the
`sign(logFC) * -log10(PValue)` of its tested peak with the largest
absolute value. (Until this, ChIPseeker gave each peak one nearest
gene and broke ties depending on the other peaks it annotated, so a
comparison's results changed with the other comparisons run with
it.) A collection is `COLL` or `COLL:SUBCOLL` as msigdbr
names them; `C2:CP` matches `CP` and every `CP:*`, and an item with
a leading `-` removes those sets from the others. A name msigdbr
does not have stops the run and lists the valid ones. Each item is
its own fgsea run, so `padj` is BH within it. A positive NES means
more accessible in Group2.

Human builds use the human MSigDB and mouse the mouse MSigDB
(`MSIGDB_DB`), so the collection names differ (`H` vs `MH`, `C5` vs
`M5`). The significant sets (padj < 0.05) of each comparison and
collection are collapsed to independent ones (`collapsePathways`;
`CollapsedInto` in the xlsx), and every plot shows only those. The
output, in `atacSeq/diff` beside its input, is, all starting
`<PROJ>_<RUNTAG>_PathwayGSEA`:

| File | Contents |
| --- | --- |
| `.xlsx` | summary, a sheet per comparison with every tested set (not filtered by padj), parameters |
| `.pdf` | NES bar chart, at most 15 sets per side; no page for a comparison and collection with no significant set |
| `_Dot.pdf` | the same sets as a dot plot (size: genes in set; color: padj) |
| `_Heatmap.pdf` | NES across comparisons, per collection, at most 50 sets; only with 2+ comparisons |
| `_Enrichment.pdf` | running-score plots of the top 3 sets per side, and a leading-edge heatmap (logCPM z-score per sample, at most 50 genes) for each |
| `.rds` | the gene ranks and full results; not delivered |

The leading-edge heatmaps use the logCPM that
`diffAnalysisPairwise.R` saves in its `.rds`; an `.rds` written
before they were added has none, and they are left out with a
warning (rerun `diffAnalysisPairwise.R`). Two runs
with different options need different `--out-tag`s, or the second
overwrites the first.

msigdbr (CRAN, 26.1.1 here) downloads MSigDB 2026.1 (md5 checked)
on first use into `tools::R_user_dir("msigdbr", "cache")`
(`~/.cache/R/msigdbr`); it downloads to `tempdir()` first.

### Motif enrichment

`R/motifEnrich.R` tests the differential peaks of each comparison for
enriched transcription factor motifs with monaLisa, against the
JASPAR 2024 CORE vertebrate non-redundant profiles in `lib/motifs`
(`lib/motifs/getJaspar.sh` downloads them again and checks the md5).
Run it from the analysis directory after `diffAnalysisPairwise.R`,
with the same RUNTAG:

```
Rscript /path/to/ATAC-seq/R/motifEnrich.R [options] [RUNTAG]

  -w, --width      region around each peak's strongest summit, odd (501)
  -m, --min-peaks  smallest set tested (50)
  -c, --cores      cores for the motif scan (SLURM_CPUS_PER_TASK, else 1)
```

Each tested peak is the `--width` bp centered on its strongest summit
(the largest MACS2 -log10 p of the `*_summits.blacklistFiltered.bed`
of all samples in `atacSeq/macs`), read from the `FASTA` of the build
in `GENOMES.tsv`. For each comparison the up (FDR < 0.05, logFC > 0)
and down (FDR < 0.05, logFC < 0) sets are tested against the
comparison's other tested peaks, which monaLisa weights to each set's
G+C and k-mer composition; `padj` is BH within the set. The motif scan
takes most of the time, about 26 sequences per second per core on a
cmobic_cpu node, and every tested peak of every comparison is scanned
(the background too), so give the job cores. Proj_18143 (7
comparisons, 1.14 million sequences scanned) took 80 min with 8 cores
and reached the 64G it was given (sacct MaxRSS 64.0G, which counts
the pages the forked workers share once per worker); ask for 128G
(the nodes have about 1 TB). The output, in `atacSeq/diff`, all
starting `<PROJ>_<RUNTAG>_MotifEnrichment`:

| File | Contents |
| --- | --- |
| `.xlsx` | summary, a sheet per comparison with every motif tested for each set (not filtered by padj), parameters |
| `.pdf` | heatmap of the log2 enrichment of the top 20 motifs per set of each comparison, then of the top 50 across comparisons; no page for a comparison with no enriched motif |
| `.rds` | the regions and full results; not delivered |

monaLisa needs TFBSTools, whose dependency DirichletMultinomial
compiles against GSL. IRIS has none; a static GSL 2.8 is in
`/usersoftware/core001/common/RHEL_8/gsl/2.8`, needed only to
install it:
`PATH=/usersoftware/core001/common/RHEL_8/gsl/2.8/bin:$PATH Rscript -e 'BiocManager::install(c("monaLisa", "TFBSTools"))'`.

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
  without the `.rds` files of `atacSeq/diff`, and `docs/output` and
  `docs/methods` of this checkout (PDF and HTML) to `atacSeq/docs`.
  Diff or GSEA files left in the analysis directory by an older
  checkout are not delivered; the script warns and names them;
- writes the email, from `tmpldeliveryEmail_*.txt`, to
  `deliveryEmail_<no>_r_NNN_<yymmdd>.txt`;
- registers the run with the bicdelivery permissions database
  (`~/Code/BIC/Delivery/Version2j`). If this fails, the files are
  already delivered; the script prints the command to rerun this
  step alone.

The email gives the release in `VERSION.md` of the checkout the
script is run from (see "Releases") and links to the docs at that tag.
Diff or GSEA files with a RUNTAG (the email names the untagged
files) give a `WARNING`: edit the email before sending it. A delivery folder other
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

