# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Features

- **Run status** (2026-10-09): there is now a way to tell whether a run
  worked. `pipe.sh` writes `00.RUNSTATUS.txt` (`STATUS=RUNNING`, then
  `COMPLETED` or `FAILED` with the stage it failed in) and a manifest of
  every stage job, `SLURM.CTRL/jobs.tsv`. Every stage job runs through
  the new `bin/runStage.sh`, which ends its log with `#ATAC_EXIT=<rc>`.
  The new `bin/checkRun.sh` reads these with `sacct` and exits 0
  (worked), 1 (failed) or 2 (running); it also reports a control job
  that was killed before it could update the status file. `COMPLETED`
  now requires every sample's bigWig, peak file, insert-size metrics and
  TSS enrichment, not just `macsPeaksMerged.saf`.
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
- `checkRun.sh` reports `UNKNOWN` (exit 3) when `sacct` fails, instead of
  reporting a live run as `FAILED`.
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
