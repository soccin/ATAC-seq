# Porting the ATAC-seq pipeline from JUNO/LSF to IRIS/Slurm

Done 2026-10-08 on `feat/slurm`; released in v1.5.0. JUNO is gone, so
there is no LSF path left in the code; `attic/lsfTools.sh` is reference
only. The PEMapper port (`Proj_18143_B/PEMapper/docs/LSF_SLURM_PORT.md`)
was the guide; the cluster
facts measured there (hard memory caps, no epilogue in the job log, no
name-glob dependencies, `kill_invalid_depend` off) all apply here.

## Strategy: keep the control job, swap the three primitives

PEMapper rebuilt its pipeline as a dependency graph. This pipeline keeps the
simpler control-job model it already had, because the wrappers for it exist
in `~/bin`:

| LSF | Slurm | Where |
| --- | --- | --- |
| `bsub -q control ... ./pipe.sh` | `sbatch pipe.sh` | `~/bin/sbatch`, which exports `SBATCH_SCRIPT_DIR` |
| `bsub` | `~/bin/bsub` | SchedMD openlava shim with local patches |
| `bSync NAME` | `bSync.sh "^NAME$"` | `~/bin/bSync.sh`, wrapped in `bin/slurmTools.sh` |
| `bCheck NAME` (`bjobs` EXIT count) | `bCheck NAME` (`sacct` State) | `bin/slurmTools.sh` |
| `parseLSF.py` sweep | `bCheckAll REGEX` (`sacct`) | `bin/slurmTools.sh` |

The control job runs on `cmobic_cpu` with `--qos=priority` because a full
run can exceed two hours; those settings are `#SBATCH` directives in
`pipe.sh`. Slurm runs a copy of a batch script from its spool directory, so
`$0` cannot resolve `SDIR`. `~/bin/sbatch` solves this for every script:
it rewrites a temporary copy with `export SBATCH_SCRIPT_DIR=<dir>` placed
after the `#SBATCH` block, and `pipe.sh` sets
`SDIR=${SBATCH_SCRIPT_DIR:-$(dirname $0)}`. The exports are baked into the
script rather than passed with `--export` because any non-`ALL` entry in
`--export` triggers a user-environment probe that is broken on IRIS.

## Flag translation in pipe.sh

| LSF | Slurm (through the shim) | Note |
| --- | --- | --- |
| `-R "rusage[mem=24]"` with `-n 1` | `-M 24G` | The shim has no `-R`. LSF `rusage` was per slot and advisory; `--mem` is the job total and a hard cgroup cap. |
| `-n 3 -R "rusage[mem=6]"` | `-n 3 -M 18G` | Same total, written out. |
| `-n 16 -R "rusage[mem=8]"` (TSS) | `-n 2 -M 32G` | `tss_enrich.R` is single-threaded and reads the BAM in chunks; 128G was never used. |
| `-W 359` / `-W 59` | class `SHORT` (`-W 1:59:00`) or `LONG` (`-W 12:00:00` and qos priority) | See "Walltime classes" below. The shim sends anything under two hours to `cmobic_short,cpushort`, the rest to `cmobic_cpu`. `-W` is sbatch format, so `-W 59` means 59 minutes in both, but `-W 5:00` would be five minutes. |
| `-o LSF.01.POST/` | `-o SLURM.01.POST/` | Trailing slash: the shim creates the directory and writes `%j.out`. |
| `-w "post_done(NAME)"` | removed | The shim maps it to `afterok:NAME`, but Slurm needs numeric ids, and an `afterok` whose parent failed pends forever (`kill_invalid_depend` is off). `MergePeaks -> Count -> DESEQ` now runs in sequence with `bSync`. |
| `-q control` | `#SBATCH -p cmobic_cpu` and `#SBATCH --qos=priority` | Control job only; directives in `pipe.sh`. |

Memory per stage: POST 32G (the `~/bin/picard` wrapper runs `-Xmx24g`, so
the cap has to sit above that heap), BW 24G (`sort -S20g`), CALLP 3 cores
18G, MergePeaks 3 cores 24G (`sort -S 16g`), Count 10 cores 24G
(`featureCounts -T 10`), Index 4G, DESEQ 24G, TSS 2 cores 32G.

## Walltime classes

Every stage job is passed a class by `atacSub`:

| Class | `-W` | Partition (chosen by the shim) | qos |
| --- | --- | --- | --- |
| `SHORT` | `1:59:00` | `cmobic_short,cpushort` | none (`normal`) |
| `LONG` | `12:00:00` | `cmobic_cpu` | `priority`, as `SBATCH_QOS=priority` on that one `bsub` call |

The qos cannot be set for the whole run. `cpushort` allows only
`qos=normal` and `EnforcePartLimits=ALL` is set, so a job sent to
`cmobic_short,cpushort` with priority qos is rejected (`sbatch
--test-only`, 2026-10-09: `Invalid qos specification`). For the same
reason `pipe.sh` unsets every `SBATCH_*` variable it inherits from the
shell that submitted it: bsub submits with `--export=ALL`, so an exported
`SBATCH_QOS` would reach every stage job.

`SHORT` whenever possible: on 2026-10-08, 249 jobs on `cmobic_cpu` waited
41 min on average and up to 7.4 h to start, while `cmobic_short` jobs
start within minutes. So the stages whose run time grows with their input
choose a class per job: `runClass MIN_PER_GB FILE...` multiplies the
input size by a measured rate and picks `SHORT` if the estimate is at
most 60 min (half the `SHORT` walltime; `ATAC_SHORT_MAX_MIN` overrides
it), `LONG` otherwise. The estimate is logged in the control log.

Rates and results from the first full-size run (`Proj_18143_B`, 12 b38
samples, input BAMs 4.6 to 14.1 GB, control job 18335688, 2026-10-09,
every job `SHORT` at the time). No job reached `OUT_OF_MEMORY` or
`TIMEOUT`.

| Stage | Input for the rate | Elapsed (max) | Measured min/GB | `RATE_*` | Largest input that stays `SHORT` |
| --- | --- | --- | --- | --- | --- |
| POST | input BAM | 1:04:05 | 3.5 to 5.3 | 6 | 10 GB |
| BW | `shifted.bed.gz` | 0:27:25 | 24.5 to 26.0 | 30 | 2 GB |
| CALLP | `shifted.bed.gz` | 0:08:34 | 8.1 to 9.5 | 10 | 6 GB |
| TSSE | postProcess BAM | 0:10:13 | 1.6 to 2.9 | 3.5 | 17 GB |
| Count | all postProcess BAMs, summed | 0:35:04 | 0.94 | 1.2 | 50 GB |
| Index | | 0:01:31 | | always `SHORT` | |
| MergePeaks | | 0:00:30 | | always `SHORT` | |
| DESEQ | | 0:00:15 | | always `SHORT` | |
| control job | | 2:22:55 | | `#SBATCH`, 3 days | |

With these rates only the 14.1 GB sample's POST job of that run would
have been `LONG` (estimate 85 min; it took 64). POST is single-threaded
and I/O bound (TotalCPU close to Elapsed; it read 120 GB and wrote 102 GB
for the largest sample), so its time depends on shared filesystem load
as much as on the node.

Memory was not changed. `sacct MaxRSS` counts page cache here: POST
reports 32.0G against its 32G cap and Count 24.0G against 24G, both
without being killed, so the figures say only that the jobs read a lot
of data. An `OUT_OF_MEMORY` state is the signal to raise a cap; none
occurred. The control job peaked at 1.1G of its 16G.

CPU requests are larger than the use. TotalCPU/Elapsed: CALLP about 1.2
of 3 cores, TSSE 1.0 of 2, MergePeaks about 2 of 3, Count 0.75 of 10
(I/O bound across 12 BAMs). They were left as they are (00.ISSUES.md).

## Checking a run

Slurm writes nothing into a job's log when it fails. `bCheck` asks `sacct`
for every job with the stage name since `ATAC_START` and aborts on any
State other than `COMPLETED`. Match State, never ExitCode: an OOM kill
reports `OUT_OF_MEMORY` with `ExitCode=0:125`. `pipe.sh` ends with
`bCheckAll "^qATAC-Seq_.*_<pid>$"`, the replacement for the old
`parseLSF.py | fgrep -v Successfully` sweep.

A failing Slurm client call must not abort a run that is working.
`~/bin/bSync.sh` runs `squeue ... 2>/dev/null` and reads empty output as
"all jobs done", so a single `squeue` failure (controller timeout) makes
it return while the stage is still running. `bCheck` used to retry
`sacct` for 60 s, then count the still-running jobs as failures and abort
the run, and the `EXIT` trap cancelled every job. A failed `sacct` call
ended `pipe.sh` at once, since `rows=$(...)` under `set -e` exits on a
nonzero status. `bCheck NAME NJOBS` now:

- goes back into `bSync` while `sacct` shows any job of the stage active;
- treats a failed `sacct` call, no rows, or fewer rows than `NJOBS` (the
  number of jobs `waitStage` counts in `jobs.tsv`) as "not yet", and
  retries every 30 s for up to `ATAC_SACCT_WAIT` seconds (1800);
- fails the stage only for a job in a final state other than
  `COMPLETED`.

A `bsub` that is rejected still aborts the run; there is no retry,
because a submission whose reply timed out may have created the job, and
a retry would run it twice.

That is enough for the control job to stop on a failure, but not for a
person to find out afterwards whether a run worked. Under LSF every log
ended with `Successfully completed.` or `Exited with exit code N`; under
Slurm a killed control job leaves a log that just stops. `pipe.sh`
therefore keeps three records in the analysis directory, following the
PEMapper port:

| Record | Written by | Content |
| --- | --- | --- |
| `00.RUNSTATUS.txt` | `pipe.sh` | `STATUS=RUNNING`, then `COMPLETED` or `FAILED`; the `STAGE` reached, `RC`, `MESSAGE`, control job id and log, `ATAC_START`, version, genome, samples, BAMs. Once the R reports have run, the post-run notes (formerly `00.POST_RUN.txt`) follow the `KEY=VALUE` block |
| `SLURM.CTRL/jobs.tsv` | `pipe.sh` (`atacSub`) | one line per stage job: job id, stage, sample, log path |
| `#ATAC_EXIT=<rc>` | `bin/runStage.sh` | last line of every stage log; every `bsub` goes through this runner |

`STATUS=COMPLETED` is written only as the last action of `pipe.sh`, after a
check that every sample has its bigWig, peak file, insert-size metrics and
TSS enrichment (when the TSSE stage ran; it is skipped for builds without
files in `R/TSSEnrich/lib/`, such as mm10) and that the atlas, count
matrix, scale factors and QC PDFs exist. Any earlier exit runs the `EXIT` trap, which writes
`STATUS=FAILED` with the stage and cancels the rest of the run's jobs.

`bin/checkRun.sh`, run from the analysis directory (or given its path),
reads the three records with `sacct` and exits 0 if the run worked, 1 if
it failed and 2 if it is still running. If a `sacct` call fails it
prints `UNKNOWN` and exits 3, unless the status file already says
`FAILED` or every stage log trailer shows success; without `sacct` a live
run would otherwise read as failed. It does not take `STATUS=RUNNING`
at its word: if the control job was killed outright (SIGKILL on its
memory cap, node failure) the trap cannot run, the file stays at
`RUNNING`, and `checkRun.sh` reports `FAILED` because `sacct` shows the
control job ended. A stage log without an `#ATAC_EXIT=` line belongs to a
job that never reached the end of its script.

```bash
/path/to/ATAC-seq/bin/checkRun.sh          # in the analysis directory
grep -L "#ATAC_EXIT=0" SLURM.0*/*.out      # stage logs that did not exit 0
sacct -X -u $USER -S <ATAC_START> --format=JobID,JobName%40,State,ExitCode,Elapsed
```

`scancel` and the time limit send SIGTERM to the batch shell only, and
bash runs a trap only once its foreground child exits; Slurm sends
SIGKILL 30 s later (`KillWait`). `bSync` therefore runs `bSync.sh` in the
background and `wait`s on it, so the trap runs at once while the control
job is waiting on a stage, which is nearly all of its run time. During
the R reports at the end of the run the trap may not get to run; then
`checkRun.sh` still reports the run `FAILED` from `sacct`.

The records describe one run, and there is one run per analysis
directory. Before it writes anything, `pipe.sh` stops if
`00.RUNSTATUS.txt`, `SLURM.CTRL/jobs.tsv`, `out/`, `callpeaks/` or
`atacSeq/` exists. `SLURM.CTRL/` itself is not checked: Slurm creates it
for the control log. The refused run exits before the `EXIT` trap is
installed, so it leaves the earlier run's records alone; the reason is
in its own control log. A control job that Slurm requeues after a node
failure (`JobRequeue=1`) is stopped the same way, by the records of its
first attempt, and `checkRun.sh` reports the run `FAILED` from `sacct`.

## Environment changes

| What | Was | Now |
| --- | --- | --- |
| samtools | `module load samtools` | `module load samtools/1.20` in `bin/loadTools.sh` |
| bedtools | `module load bedtools/2.27.1` | `~/bin/bedtools` (2.31.1) on PATH; no module exists |
| picard | `picardV2` from the cluster PATH | `~/bin/picard` (picard 3.4.0, `-Xmx24g`) |
| scratch | `/scratch/socci/_scratch_ATACSeq` | `${ATAC_SCRATCH_ROOT:-/scratch/core001/bic/$USER/ATACSeq}` |
| TMPDIR | node `/tmp` | `${ATAC_LOCAL_TMP:-/localscratch/$USER}`, one directory per stage job |
| venv python | 3.9 required | the `python3` on PATH (3.10); MACS2 2.2.9.1 and IDR 2.0.3 both build |
| R | cluster R | 4.5.1 on PATH; `optparse` must be installed for `tss_enrich.R` |

### Temp files: /localscratch, never /tmp

IRIS compute nodes export `TMPDIR=/tmp`. On the nodes probed on
2026-10-09 (jobs 18394759, 18394762, 18394763: `cpushort`,
`cmobic_short`, `cmobic_cpu`) `/tmp` is a 137G volume shared by every
job on the node, one already 23% full, while `/localscratch` is a 2.8T
node-local disk, world-writable, with no per-user directory made in
advance. `makeBigWigFromBEDZ.sh` (`sort -S20g`) and `mergePeaksToSAF.sh`
(`sort -S 16g`) spilled past their buffers into `/tmp`, and the R
session tempdir and Python tempfiles went there too.

- `bin/loadTools.sh`, sourced by `pipe.sh` and every stage script, sets
  `TMPDIR=${ATAC_LOCAL_TMP:-/localscratch/$USER}` and creates it, unless
  `TMPDIR` is already a directory under that root.
- `bin/runStage.sh` creates `atac.<jobid>.XXXXXX` under that root for
  each stage job, exports it as `TMPDIR`, logs it as `#ATAC_TMPDIR=`,
  and removes it on exit, including after the SIGTERM Slurm sends on
  `scancel` or a time limit.
- The large sorts pass `-T "$TMPDIR"` explicitly.

Not affected: picard (`~/bin/picard` sets `TMP_DIR` and
`java.io.tmpdir` to `/scratch/core001/bic/socci/PICARD/$$`), MACS2
(`callPeaks_ATACSeq.sh` exports `TMPDIR` to its `/scratch` work
directory), and featureCounts (temp files go next to its output).

## Exit-status fixes made along the way

`postMapBamProcessing_ATACSeq.sh` backgrounded `CollectInsertSizeMetrics`
and never waited for it. Under LSF the BED pipeline usually outlasted it;
under Slurm the job ends when the script exits and picard would be killed
mid-write. It now `wait`s, runs under `pipefail`, and exits non-zero on any
failed step so `bCheck` sees it. `makeBigWigFromBEDZ.sh` and
`mergePeaksToSAF.sh` run under `pipefail` for the same reason.

Added after v1.5.0: `callPeaks_ATACSeq.sh` runs under `pipefail` and stops
if its chromosome filter fails, rather than running MACS2 on a truncated
BED. `makeBigWigFromBEDZ.sh` checks the read count that sets its scale
factor (and rejects zero) and the bigWig build itself. `pipe.sh` `usage`
exits 1, and `pipe.sh` stops up front if a BAM has no `@RG SM` tag or two
BAMs share one. It also detects the build of every BAM, not just the
first, and stops up front if they differ or if the build is not one
`postMapBamProcessing_ATACSeq.sh` accepts (b37, b38, mm10); before, a
build such as hg19 failed in every POST job. For mm10 the TSSE stage is
skipped with a warning instead of failing the run after stage 6. The
Count stage passes `-Q $MAPQ` to featureCounts instead of a fixed `-Q
10`.

## Validation

Eleven b38 BAMs downsampled 100x (`PortATACseq/out/picard/DSB/*/*.dn_100.bam`,
about 1.2M reads each) run from `PortATACseq/test01`. Control job
18237528: 70 stage jobs, all `COMPLETED`, 5m38s wall, every deliverable
present under `atacSeq/`. A first attempt (18236503) failed at DESEQ only
because the test harness overrode `R_LIBS_USER`; `bCheck` caught it and
the `EXIT` trap cancelled the queue, which is the failure path working.

The genome check and the TSSE skip were tested on 2026-10-10 with the
same BAMs. A normal run (test15, control job 18460598) completed with 58
of 58 stage jobs, 11 TSS enrichment files and `featureCounts -Q 10`. A
run from a copy of the checkout without `R/TSSEnrich/lib/b38_tss.bed`, standing in
for mm10, with `-q 20` (test16, 18460600) completed with 47 of 47 jobs,
no TSSE jobs, the warning in `MESSAGE`, and `featureCounts -Q 20`.
`pipe.sh` run by hand on a header-only BAM with an unrecognized `@SQ`
set, alone and next to a b38 BAM, stopped with rc 1 and the reason in
`MESSAGE`, before submitting anything.

## Not done

- `deliverResults.sh` still points at `/ifs/res/seq/pi/invest` and
  `~/Code/BIC/Delivery`; not a scheduler change and not touched. It is
  `00.ISSUES.md` #1.
- The `RATE_*` values come from one project on one genome. Check the
  `runClass` lines in the control log against `sacct` Elapsed on new
  projects, and raise a rate if a job gets close to its walltime.
