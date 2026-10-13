# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

To be released as v2.0.0: the counting engine changes, and counts from
any v1.x run may be wrong for one or more samples.

### Breaking changes

- **featureCounts v1.6.0 replaced with v2.1.1** (`-p --countReadPairs`).
  v1.6.0, run with several threads on several BAMs, sometimes did not
  pair the reads of one BAM and counted every read as a fragment, about
  doubling that sample's counts, its DESeq2 scale factor and its place in
  the PCA, with exit 0 and no error. It was found in 5 projects run since
  2022. On the same tests v2.1.1 failed 0 of 250 times, and for samples
  that v1.6.0 handled correctly its counts are identical, cell for cell.
- **POST removes orphan reads, as the ENCODE ATAC-seq pipeline does.**
  After the MAPQ and flag filter, `postMapBamProcessing_ATACSeq.sh` now
  runs `samtools fixmate -r` on the name-sorted reads and the flag
  filter again, which removes a read whose mate failed the MAPQ cutoff;
  such orphans were kept before. The postProcess BAM, the shifted BED,
  and so the bigWigs, peaks, atlas and counts change slightly; Check A
  now gives exactly 1.000. On the downsampled tests each sample lost
  116-371 orphan reads.

  Two filters that ENCODE does not have remove reads that are not from
  the reference. BWA-MEM (PEMapper) is a local aligner: it aligns the
  best 19-35 bp seed of such a read and soft-clips the rest. ENCODE
  aligns end to end with bowtie2, which leaves these reads unaligned.
  The mapper is not changed (C-21 of the ENCODE review). In Proj_18143
  these reads were a tandem repeat, the same in every sample.

  - **New: a read whose 5' end is soft-clipped is removed, and its
    mate with it.** A Tn5 read starts at the cut site, so its 5' end
    must align. Such reads pile up by the hundreds on single bases, at
    the same loci in every sample, and MarkDuplicates cannot remove
    them because the clip length varies. In ED1-1 of the earlier
    Proj_18143 run (insert filter on), the 41 bp
    `chr3:167846660-167846700` held 2,864 reads, 85 of them with the
    5' end aligned, and its atlas row had the 5th highest count of the
    sample. The filter removes 0.5-0.9% of the reads with an insert
    over 30 bp (Proj_18143 BD1-2 and ED1-1), and their mates. A 3'
    clip is kept: it is mostly the Nextera adapter, which PEMapper
    does not trim.
  - The filter that removes pairs with an insert of 30 bp or less
    stays. It was taken out for ENCODE concordance and put back after
    the full-size Proj_18143 run: there 19-30 bp fragments were 14% of
    all fragments and 45-55% in four samples, and 98% of them lay in
    peaks, so they made up 38-75% of the counts of ten samples. They
    are the same seeds; the 5' filter removes 99.6-99.8% of them and
    this filter the rest.
- **MACS2 now runs with ENCODE's settings**, `--shift -75 --extsize
  150 --keep-dup all`. Since the first version of
  `callPeaks_ATACSeq.sh` it ran with `--shift 75` and no `--keep-dup`:
  - The sign of the shift was lost when the settings were copied from
    the Kundaje lab pipeline (ENCODE's predecessor). Each read
    contributed the region 75 to 225 bp downstream of its cut site,
    instead of the 150 bp centered on it, so the signal of an
    accessible region was split between the strands and moved off
    center.
  - The MACS2 default `--keep-dup 1` keeps one read per cut site and
    strand. POST has already removed PCR duplicates, so the reads it
    dropped were distinct fragments; on the full-size Proj_18143 run
    that was 5 to 17% of each sample's reads, more in deeper samples.

  Every narrowPeak, the atlas, the count matrix, the scale factors, the
  QC PDF and any differential analysis of an earlier run were made with
  the old settings. BAMs, shifted BEDs, bigWigs, insert sizes, the
  mitochondrial fraction and TSS enrichment do not use MACS2 and do not
  change. On the downsampled test sets the fraction of a sample's reads
  in its own peaks went from 0.70-0.78 to 0.90-0.95 (hg38) and from
  0.59-0.69 to 0.86-0.89 (b37), with up to 7% (hg38) and 8-17% (b37)
  fewer peaks (`docs/SLURM_PORT.md`, "Validation").
- **Peaks in ENCODE blacklisted regions are left out of the atlas and
  the counts**, as in the ENCODE ATAC-seq pipeline. There was no
  blacklist filter before, so artifact regions (satellite repeats,
  centromeres, high-signal regions) were in the atlas and the counts.
  On the full-size Proj_18143 run they were 1,376 of 654,370 atlas
  peaks with 0.4% of all counts, but 3 of the 10 and 7 of the 100
  peaks with the most counts, which weigh on the scale factors, the PCA
  and the differential analysis. MACS2's own files in
  `callpeaks/<sample>/` (and `atacSeq/macs/<sample>/`) are unchanged.
  Beside them, `callPeaks_ATACSeq.sh` writes two new files:
  `*_peaks.blacklistFiltered.narrowPeak`, MACS2's peaks without every
  peak that overlaps `lib/genomes/<TAG>.blacklist.bed` and with the
  score capped at 1000 (ENCODE's `encode_lib_blacklist_filter.py`),
  and `*_summits.blacklistFiltered.bed`, the summits of those peaks.
  Both are delivered, and the deliverables check requires them and
  MACS2's three files for every sample.
  The atlas is built from the `blacklistFiltered` peaks only.
  `mergePeaksToSAF.sh GENOME MACS_OUTDIR` (it now takes the genome)
  also subtracts the blacklist from the merged atlas, because merging
  within 500 bp can bridge a blacklisted region between two kept
  peaks; such a peak is split in two. Blacklists: hg38 `ENCFF356LFX`,
  mm10 `ENCFF547MET`, b37 the hg19 list `ENCFF001TDO` without `chr`
  and `chrM`, built by `lib/genomes/getBlacklists.sh`
  (`docs/GENOMES.md`). `pipe.sh` stops if the genome's blacklist is
  missing.
- **Each sample keeps at most its 300,000 strongest peaks**, as in the
  ENCODE ATAC-seq pipeline (`encode_task_macs2_atac.py`). Before the
  blacklist filter, `callPeaks_ATACSeq.sh` now sorts MACS2's peaks by
  `-log10 p`, keeps the top 300,000, clamps negative coordinates to 0
  and clips peak ends to the contig ends (`bin/bedClip -truncate`,
  the UCSC v369 build ENCODE uses). The blacklist filter is applied
  to that set, so `*_peaks.blacklistFiltered.narrowPeak` and its
  summits now hold at most 300,000 peaks, sorted by p-value. MACS2's
  own files are unchanged, and the peaks keep MACS2's names (ENCODE
  renames them `Peak_N`). On the full-size Proj_18143 run (old MACS2
  settings) 9 of 15 samples had more than 300,000 peaks, up to
  396,453. On the downsampled tests, no hg38 sample reached the cap;
  3 of 4 b37 samples did, and the b37 atlas went from 758,830 to
  727,125 peaks.
- **The atlas `Start` is 1-based**, as SAF requires.
  `mergePeaksToSAF.sh` wrote the 0-based start of `bedtools merge` into
  `macsPeaksMerged.saf`, and featureCounts reads SAF as 1-based, so
  each atlas peak was counted one base wider on the left. The `Start`
  of `macsPeaksMerged.saf`, `peaks_raw_fcCounts.txt` and the
  `DiffPeaksEdgeR` workbooks is now 1 higher than before, `End` is
  unchanged, and the peaks given to ChIPseeker for gene annotation are
  one base shorter on the left, as they should be.
- The `.summary` file has the 2.x status rows: `Unassigned_Read_Type` and
  `Unassigned_Singleton` are new, `Unassigned_Nonjunction` is now
  `Unassigned_NonSplit`, and the row order changed. The R reports keep
  only `Assigned` and are not affected. The featureCounts log says
  "alignments" where it said "fragments".
- **GRCh38 is tagged `hg38`, not `b38`** (the name the Broad and UCSC
  use), in `GENOME=` of `00.RUNSTATUS.txt`, file names and every
  argument; `b38` is now an unsupported genome. Genome files are
  `lib/genomes/{b37,hg38,mm10}.genome`; the TSS files are
  `R/TSSEnrich/lib/hg38_tss.bed` and `hg38.chrom.sizes`.
- **Mitochondrial reads are left out for every build**, as in the ENCODE
  ATAC-seq pipeline (`docs/GENOMES.md`). Before, the contig filter kept
  `MT` for b37 and `chrM` for mm10 but not `chrM` for b38, so peaks,
  bigWigs (including the 10M-read scaling), the atlas and the counts of
  **every earlier b37 and mm10 run** included them (5.2% of the reads in
  one b37 sample checked). Earlier b38 runs are not affected. The
  postProcess BAM still holds every read.
- **The `FDR` of `R/diffAnalysisPairwise.R` is now the false discovery
  rate** (`p.adjust(method = "fdr")`). Since the first diff script in
  2022 it was `p.adjust()` with no method, which is Holm's
  family-wise correction, so the `FDR` column and the FDR < 0.05 cut
  that picks the significant peaks were Holm. Every earlier
  `DiffPeaksEdgeR` workbook lists far fewer peaks than an FDR of 0.05
  gives; nothing in it is a false call at that level. On the Proj_18741
  count matrix the five comparisons went from 752 to 18,872 significant
  peaks to 14,032 to 100,697 (3.6 to 33 times as many).
- **The PCA in `*__ATACSeqQC.pdf` and `*DiffPeaks_V4.pdf` is the usual
  PCA of the samples** (log-CPM, each peak centered across samples), with
  the variance explained on the axes. Before, the peaks were the
  observations and the plot showed the rotation: its PC1 was the profile
  that all samples share (81% of the variance on Proj_18143, every sample
  at -0.24 to -0.29) and carried no information about them, its PC2 was
  the real PC1, and the real PC2 was not shown.
- **`R/diffAnalysisPairwise.R` takes no genome argument**:
  `Rscript R/diffAnalysisPairwise.R MANIFEST COMPARISONS [RUNTAG]`, run
  from the analysis directory, whose `00.RUNSTATUS.txt` gives the genome.
  It exits 1 on a usage error or a missing or unknown genome (it exited
  0). The `chr` prefix is added to peak contigs for b37 only, by rule,
  not by looking at the first peak.
- `bin/computeTSSEnrich.sh GENOME BAM`: the genome is passed by
  `pipe.sh`, not detected again from the BAM.
- `sCer+sMik_IFO1815` and its genome files moved to `attic/`; it was in
  two stage scripts but could never run. `R/TSSEnrich/v1/` (superseded,
  not called) moved to `R/attic/TSSEnrich_v1/`.

### Features

- **Count checks.** A run now fails, before DESEQ, unless both pass:
  - A: each sample's fragment total is 0.95 to 1.05 times its mapped
    reads / 2 (`samtools idxstats`). An unpaired sample gives 2.0.
  - B: each sample's featureCounts `Assigned` equals, exactly, an
    independent count of the fragments over the atlas made with samtools
    alone (one new `CountCheck` job per sample, run beside Count).

  The first failing line goes into `MESSAGE`; every sample's numbers are
  in the control log and in `peaks_raw_fcCounts.checks.txt`.
  `bin/checkCounts.sh` can be run by hand to audit an earlier run.
- The Count stage is now `bin/countPeaks.sh`. It waits for the BAM
  indexes, writes its temp files to the job's `/localscratch` directory
  (`--tmpDir`), and asks for 8 cores instead of 10.
- The deliverables check also requires
  `atacSeq/atlas/peaks_raw_fcCounts.txt.summary` and
  `peaks_raw_fcCounts.checks.txt`.
- **`MESSAGE` says why a stage failed.** It was empty when a stage job
  failed; the reason was only in that job's log. Now, for each failed
  job (at most 3), it gives the sample, job id, `sacct` state and the
  first line of the log with `FATAL` (else `error` or `fail`, any case),
  e.g. `Count failed: all job 18905720 FAILED: FATAL ERROR: Check A
  failed on peaks_raw_fcCounts.txt.summary`. `checkRun.sh` prints it.
- **TSS enrichment for mm10.** `R/TSSEnrich/lib/mm10_tss.bed` (ENCODE's
  GENCODE vM7 file, rebuilt byte for byte by `lib/getMM10.sh` because the
  ENCODE bucket now returns 404) and `mm10.chrom.sizes`. TSSE now runs
  for every build and is a required deliverable; the mm10 skip and its
  warning are gone. A supported genome without its TSS files stops the
  run before any job is submitted.
- **Mitochondrial fraction per sample**: POST writes
  `<SID>_postProcess___MITO.txt` (reads, mitochondrial reads, fraction,
  counted after the MAPQ/flag filter and before the contig filter),
  staged to `atacSeq/metrics` and required for `COMPLETED`.
- **One genome table**, `lib/genomes/GENOMES.tsv`, read by the new
  `bin/genomeInfo.sh TAG FIELD` and `R/genomeInfo.R`, replaces the
  `case` statements in the stage scripts, `SUPPORTED_GENOMES` in
  `pipe.sh` and the `if` chain in `diffAnalysisPairwise.R`. An unknown
  tag stops with the list of supported ones. The same table, genome
  files and tags as ChIP-seq v1.0.0.
- `bin/getGenomeBuildBAM.sh` is now ChIP-seq's copy, and the two must
  stay identical. It names the PEMapper `human_hg38+mm39` hybrid
  (`hg38+mm39`) and the `mm10+b37` mouse hybrid, both unsupported, so
  such a BAM stops with a named error.

### Fixes

- The insert-size histogram (`___INS.txt`, `___INS.pdf` and so
  `*__postInsDistribution.pdf`) covers 0 to 1000 bp: Picard runs with
  `W=1000`, as ENCODE's does. Without it Picard cut the histogram, and
  the mean and SD, at median + 10 x MAD, which in a sample with many
  short fragments falls inside the nucleosome peaks.
- A BAM whose header cannot be read (missing file, not a BAM) stops the
  run before any job with "cannot identify the genome of [BAM]" in
  `MESSAGE`. Before, `getGenomeBuildBAM.sh` printed the md5 of an empty
  header as an "unknown" genome.
- `bin/loadTools.sh` no longer fails under `set -u` when `TMPDIR` is
  unset (a script with `set -u` run by hand on a login node).
- `R/diffAnalysisPairwise.R` with a manifest whose `SampleID`s differ
  from the `MapID`s (as the generated manifest does for `s_` names), or
  with an `EXC` group, stopped with "Count matrix validation failed".
  Both now work. A count matrix column with no manifest row, a manifest
  row with no column, or a duplicated `MapID` or `SampleID` stops it
  with the names; before, an unmatched PEmap column was left out
  silently.
- A group whose name contains `group` (e.g. `groupA`) gave an invalid
  contrast in `R/diffAnalysisPairwise.R`.
- Output names of `plotINSStats.R`, `R/analyzeATAC.R` and
  `R/diffAnalysisPairwise.R` use the last `Proj_*` or `B-<digits>`
  directory of the path, else `ATAC`. Before, a path with none gave
  names starting with `_`, and a path with two put both in the name.
- `plotINSStats.R` and `R/analyzeATAC.R` color up to 34 groups; with
  `scale_color_uchicago()` alone, groups past the ninth were grey.
- `plotINSStats.R` passes its text to `read_tsv()` with `I()`, as readr
  2.2.0 requires for literal data (it warned).
- `R/analyzeATAC.R` with no arguments exits 1 (it exited 0).
- CountCheck's memory cap is 12G. Its process memory is 4.0G, set by
  `sort -S 4G`, but sacct MaxRSS includes file cache and reached 8G on
  a 46M-fragment sample.
- `RATE_COUNT` lowered from 1.2 to 0.1 min/GB, measured with
  featureCounts v2.1.1 in a job (0.05 min/GB), so Count runs on the
  short partition for up to 600 GB of BAMs. `RATE_CCHECK` stays 1.5
  (measured 0.79 to 1.02).
- `docs/SLURM_PORT.md` has a new "Memory" section with what each stage
  really allocates, measured without file cache.
- `mergePeaksToSAF.sh` sorts every sample's peaks in one `sort`, which
  reads the file names from `find -print0` (`--files0-from`). It used
  `find | xargs sort`: past about 128 KB of file names (roughly 1,000
  samples) xargs ran several sorts, `bedtools merge` stopped at the
  first out-of-order peak, and the run failed; a space in a sample
  name split the file name and failed the run too.

### Not tested

- mm10 was checked only with header-only BAMs (detection, the genome
  table and the file checks); no mm10 run has been made on this version
  (`00.ISSUES.md` #14). hg38 and b37 were run (`docs/SLURM_PORT.md`,
  "Validation").

## [v1.5.1] — 2026-10-10

### Features

- **Run status** (2026-10-09): there is now a way to tell whether a run
  worked. `pipe.sh` writes `00.RUNSTATUS.txt` (`STATUS=RUNNING`, then
  `COMPLETED` or `FAILED` with the stage it failed in) and a manifest of
  every stage job, `SLURM.CTRL/jobs.tsv`. Every stage job runs through
  the new `bin/runStage.sh`, which ends its log with `#ATAC_EXIT=<rc>`.
  The new `bin/checkRun.sh` reads these with `sacct` and exits 0
  (worked), 1 (failed), 2 (running) or 3 (unknown: `sacct` failed, so a
  live run is not reported as failed); it also reports a control job
  that was killed before it could update the status file. `COMPLETED`
  now requires every sample's bigWig, peak file, insert-size metrics and
  TSS enrichment (when that stage runs), not just
  `macsPeaksMerged.saf`.
- `00.POST_RUN.txt` is gone: its notes now follow the status block in
  `00.RUNSTATUS.txt` once the R reports have run, with `$SDIR` expanded
  to the real path.
- **Walltime classes** (2026-10-09): each stage job is `SHORT`
  (`cmobic_short,cpushort`, 1h59m, no qos) or `LONG` (`cmobic_cpu`, 12 h,
  qos priority). POST, BW, CALLP, TSSE and Count choose per job from
  their input size, with rates measured on a full-size run, so large
  samples no longer risk the two-hour limit and small ones stay on the
  short partitions. The estimate is logged in the control log;
  `ATAC_SHORT_MAX_MIN` (default 60) sets the cutoff.

### Fixes

- `pipe.sh` checks the genome of every BAM, not just the first, before
  it submits anything. A run stops at once, with the reason in `MESSAGE`,
  if the BAMs are on different builds or the build is not b37, b38 or
  mm10. Before, a detected build such as hg19, b37_dmp or GRCh37-lite
  failed in every POST job, and mixed builds were not noticed.
- mm10 runs no longer fail after stage 6: there are no mm10 TSS files in
  `R/TSSEnrich/lib/`, so the TSS enrichment stage is skipped with a
  warning (control log and `MESSAGE`), and its output is not required for
  `COMPLETED`. Adding `mm10_tss.bed` and `mm10.chrom.sizes` turns it on.
- The count matrix uses the `-q MAPQ` given to `pipe.sh`; featureCounts
  had a fixed `-Q 10`, so with `-q` below 10 reads under MAPQ 10 were
  left out of the counts.
- A run in a directory that already holds one now stops at once, before
  writing anything, if `00.RUNSTATUS.txt`, `SLURM.CTRL/jobs.tsv`,
  `out/`, `callpeaks/` or `atacSeq/` exists. Rerun in a new directory.
  Before, a second run overwrote the first run's status file and job
  manifest, even while the first was still running, picked up the first
  run's samples from `out/` and `callpeaks/`, and after every stage had
  run stopped before staging into `atacSeq/` (`mkdir` without `-p`). The
  unused `out/postBams`, `out/metrics` and `out/bed` are no longer
  created.
- A transient `squeue` or `sacct` failure no longer aborts a run. Before,
  `bSync.sh` read a failed `squeue` as "stage done", `bCheck` then
  counted the still-running jobs as failed, and the run was cancelled; a
  failed `sacct` call ended `pipe.sh` outright under `set -e`. `bCheck`
  now waits again while any job is active, retries `sacct` for up to 30
  minutes (`ATAC_SACCT_WAIT`), and waits until it has a record for every
  job the stage submitted.
- `mkdir -p SLURM.CTRL` is no longer needed before `sbatch`: Slurm creates
  the log directory.
- No temp files go to the node's `/tmp` (137G, shared by every job on
  the node). `TMPDIR` is now `/localscratch/$USER` (2.8T, node-local;
  override with `ATAC_LOCAL_TMP`), and each stage job gets its own
  directory under it, removed when the job ends. The 16 to 24G sorts in
  `makeBigWigFromBEDZ.sh`, `mergePeaksToSAF.sh` and `mergeSamples.sh`
  pass `-T "$TMPDIR"`. Before this, their spill and the R and Python
  temp files went to `/tmp`.
- An `SBATCH_QOS` (or any other `SBATCH_*` variable) exported in the
  shell that submits `pipe.sh` no longer reaches the stage jobs. With
  `SBATCH_QOS=priority`, as the docs used to advise, every short stage
  job was rejected with `Invalid qos specification`.
- `callPeaks_ATACSeq.sh` runs under `pipefail` and stops if its
  chromosome filter fails; `makeBigWigFromBEDZ.sh` checks its read count
  and bigWig build. Both could exit 0 after a failed step.
- `pipe.sh` `usage` exits 1, so a submission without BAMs no longer shows
  `COMPLETED`.
- `pipe.sh` stops before submitting anything if a BAM has no `@RG SM` tag
  or two BAMs share one.
- `scancel` or a time limit on the control job now cancels the run's
  stage jobs. The `EXIT` trap never ran before: Slurm signals only the
  batch shell, which did not run the trap while `bSync.sh` was in the
  foreground. `bSync` now waits on it in the background.

### Documentation

- Add `00.ISSUES.md` (open issues, numbered in order of work) and
  `attic/zzISSUES.md` (closed issues, with the reason).
- `docs/SLURM_PORT.md`: add the walltime classes with the rates measured
  on the first full-size run, the run-status records and `checkRun.sh`,
  the `/localscratch` temp-file setup, and the 2026-10-10 validation
  runs.
- `README.md`: `checkRun.sh`, one run per directory, the genome check.
- `docs/UPDATE_TO_B38.md` and `docs/CHECKLIST_B38.md` are marked
  historical; b38 shipped in v1.1.0.

### Known issues

- `deliverResults.sh` still points at the JUNO paths
  `/ifs/res/seq/pi/invest` and `~/Code/BIC/Delivery`. Deliver by hand
  (`00.ISSUES.md` #1). The other open issues are in `00.ISSUES.md`.

## [v1.5.0] — 2026-10-09

### Breaking changes

- The pipeline runs only on IRIS/Slurm. LSF support is removed; there is
  no backward compatibility with JUNO. Submit with
  `sbatch /path/to/ATAC-seq/pipe.sh` through the `~/bin/sbatch` wrapper
  after `mkdir -p SLURM.CTRL`. `pipe.sh` stops with an error if
  `SBATCH_SCRIPT_DIR` is not set. Stage logs move from `LSF.*/` to
  `SLURM.*/`.

### Features

- **IRIS/Slurm port** (2026-10-08): JUNO/LSF is gone and the pipeline now
  runs only on Slurm. `pipe.sh` keeps the control-job model: stages are
  submitted with `~/bin/bsub` (the openlava shim), waited on with
  `~/bin/bSync.sh`, and checked with a new `bCheck` in `bin/slurmTools.sh`
  that reads `sacct` State. Every stage runs on the short partitions
  (`cmobic_short,cpushort`); the control job carries `#SBATCH` directives
  for `cmobic_cpu` with `--qos=priority` and is submitted through
  `~/bin/sbatch`, whose `SBATCH_SCRIPT_DIR` export gives `SDIR`.
  `bin/loadTools.sh` loads `samtools/1.20` and
  verifies `bedtools`. Scratch moves to `/scratch/core001/bic/$USER/ATACSeq`;
  `picardV2` is replaced by the `~/bin/picard` wrapper. The
  `MergePeaks -> Count -> DESEQ` chain is sequenced with `bSync` instead of
  `-w post_done()` (Slurm dependencies need numeric ids). `bin/lsfTools.sh`
  moved to `attic/`. See `docs/SLURM_PORT.md`. Validated end to end
  2026-10-08 on 11 b38 BAMs downsampled 100x (control job 18237528,
  70 stage jobs, all COMPLETED, 5m38s wall), and again with the control
  job submitted through `~/bin/sbatch` (18243589, COMPLETED, 5m26s).
- `00.SETUP.sh`: build the venv from the `python3` on PATH; MACS2 2.2.9.1
  and IDR 2.0.3 install under python 3.10, so the 3.9 requirement is gone.
- **b37 TSS enrichment references** (2026-08-23): Stage 7
  (`bin/computeTSSEnrich.sh` / `tss_enrich.R`) can now score b37 alignments.
  Add `R/TSSEnrich/lib/b37_tss.bed` (ENCODE hg19 unique GENCODE TSS with the
  `chr` prefix stripped to match b37 contig names `1`, `2`, …, `MT`) and
  `R/TSSEnrich/lib/b37.chrom.sizes` (contig sizes taken from a b37 BAM
  `@SQ` header, not from ENCODE or UCSC). Add `R/TSSEnrich/lib/getB37.sh` to
  regenerate both files from any b37 BAM.

### Fixes

- `postMapBamProcessing_ATACSeq.sh`: `wait` for the backgrounded
  `CollectInsertSizeMetrics` and fail on any pipeline error. Under Slurm
  the job ends when the script exits, which would have killed picard
  mid-write; and the job status now reflects a failed stage.
- `callPeaks_ATACSeq.sh`: remove the scratch directory after MACS2.
- `R/analyzeATAC.R`: drop the unused `ChIPseeker` import.

### Documentation

- Add `R/TSSEnrich/README.md` documenting the TSS enrichment module: where
  `lib/b38_tss.bed` and `lib/b38.chrom.sizes` actually came from (ENCODE's
  ATAQC reference bucket and the GDC GRCh38.d1.vd1 sequence set respectively,
  with verifying md5s dated 2026-08-23), the juno porting/validation history
  behind `tss_enrich.R`, and download recipes for adding hg19, mm9, and mm10
  lib files or building a TSS BED from a GENCODE GTF for unsupported builds.
- Add `docs/SLURM_PORT.md` (flag translation, resource choices, validation
  runs) and update `README.md` and `CLAUDE.md` for IRIS/Slurm.

### Known issues

- `deliverResults.sh` still points at the JUNO paths
  `/ifs/res/seq/pi/invest` and `~/Code/BIC/Delivery`. Deliver by hand.
- Stage walltimes and memory were measured only on downsampled BAMs. All
  stages run under the two-hour short partition limit.

## [v1.1.0] — 2026-03-12

### Features

- **GRCh38 (b38) support**: The pipeline can now process human BAM files
  aligned to the GRCh38/GDC reference (chr-prefixed chromosome names).
  - Add MD5 checksum for GRCh38 GDC reference to `bin/getGenomeBuildBAM.sh`
    so the genome is auto-detected from the BAM header.
  - Add chromosome sizes files `lib/genomes/human_b38.genome` and
    `lib/genomes/human_b38.genome.bed` (chr1-chr22, chrX, chrY; chrM excluded).
  - Add `b38` case to `postMapBamProcessing_ATACSeq.sh`,
    `callPeaks_ATACSeq.sh`, and `makeBigWigFromBEDZ.sh`.

- **Genome-aware chromosome filtering**: Replace hardcoded `egrep -v` patterns
  (`chrUn|hs37d5|...`) with `bedtools intersect -nonamecheck -b $GENOME_BED`
  in `postMapBamProcessing_ATACSeq.sh`, `callPeaks_ATACSeq.sh`, and
  `makeBigWigFromBEDZ.sh`. Filters are now driven by the allowlist genome BED
  file rather than a denylist regex, eliminating spurious bedtools naming
  warnings from mixed-prefix genomes (e.g. GRCh38_GDC with CMV contig).

- **Genome parameter in postMapBamProcessing**: `postMapBamProcessing_ATACSeq.sh`
  now accepts a required `GENOME` positional argument (`b37 | b38 | mm10`).
  `pipe.sh` passes `$GENOME` automatically.

- **LSF job failure checking** (`bin/lsfTools.sh`): New sourceable utility
  library defining `bCheck JOBNAME`. Call immediately after `bSync` to abort
  the pipeline if any jobs in the group exited with a non-zero status.
  `bCheck` is now called after every `bSync` in `pipe.sh`
  (POST2, BW2, CALLP2, MergePeaks, Count, DESEQ).

- **hg38 implementation docs**: Add `docs/UPDATE_TO_B38.md` (full analysis of
  all files requiring changes) and `docs/CHECKLIST_B38.md` (atomic task
  checklist) to guide remaining b38 work.

### Fixes

- `bin/lsfTools.sh`: Fix `bCheck` false-exit under `set -e` — `grep -c`
  returns exit code 1 on zero matches, killing the pipeline even when all
  jobs succeeded. Changed to `grep -w "EXIT" | wc -l` which always exits 0.

### Refactoring

- Move `getGenomeBuildBAM.sh` to `bin/` alongside other utilities; update
  call site in `pipe.sh`.
- Move `lsfTools.sh` to `bin/`.
- Rename `CMDS.INSTALL.MACS` to `00.SETUP.sh`; update reference in `README.md`.
- `postMapBamProcessing_ATACSeq.sh`: deduplicate repeated usage message into
  a `usage()` function; fix script name in usage string.
- `R/analyzeATAC.R`: remove unimplemented differential analysis stub
  (edgeR contrasts, `doQLFStats`, hardcoded hg19 TxDb) that was guarded
  by `quit()` and never executed. Original archived to `R/attic/analyzeATAC.R`.
- Add BAM validation (`samtools quickcheck`) on the `INPUT_BAM` argument in
  `postMapBamProcessing_ATACSeq.sh` to catch argument-order mistakes early.
- `postMapBamProcessing_ATACSeq.sh`: derive sample ID (SID) and output
  paths from BAM `@RG SM:` header tag instead of stripping `.bam` suffix.
- `R/analyzeATAC.R`, `R/getDESeqScaleFactors.R`, `plotINSStats.R`: wrap
  all `library()`/`require()` calls in `suppressPackageStartupMessages`
  to suppress Bioconductor startup noise in pipeline logs.
- All pipeline R scripts: add `cat("## START/END: <script>\n")` messages
  for log tracing.

---

## [v1.0.1] — previous release (master)

See git log for earlier history.
