# Porting the ATAC-seq pipeline from JUNO/LSF to IRIS/Slurm

Done 2026-10-08 on `feat/slurm`; released in v1.5.0. The walltime
classes, run-status records, `/localscratch` temp files and the fixes
under "Added in v1.5.1" came in v1.5.1. JUNO is gone, so
there is no LSF path left in the code; `attic/lsfTools.sh` is reference
only. The PEMapper port (`Proj_18143_B/PEMapper/docs/LSF_SLURM_PORT.md`)
was the guide; the cluster
facts measured there (hard memory caps, no epilogue in the job log, no
name-glob dependencies, `kill_invalid_depend` off) all apply here.

The GRCh38 build was tagged `b38` until v2.0.0 and is `hg38` now
(`docs/GENOMES.md`). The records of runs made before that keep the old
tag and file names (`b38_tss.bed`, `GENOME=b38`), as they were.

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
18G, MergePeaks 3 cores 24G (`sort -S 16g`), Count 8 cores 24G
(`featureCounts -T` from the job's CPUs; 10 cores before v2.0.0),
CountCheck 2 cores 12G (`sort -S 4G`; 8G until 2026-10-12), Index 4G,
DESEQ 24G, TSS 2 cores 32G. See "Memory" below for what each stage
really uses.

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
| Count | all postProcess BAMs, summed | 0:01:44 | 0.05 | 0.1 | 600 GB |
| CountCheck | postProcess BAM | 0:05:20 | 0.79 to 1.02 | 1.5 | 40 GB |
| Index | | 0:01:31 | | always `SHORT` | |
| MergePeaks | | 0:00:30 | | always `SHORT` | |
| DESEQ | | 0:00:15 | | always `SHORT` | |
| control job | | 2:22:55 | | `#SBATCH`, 3 days | |

The Count and CountCheck rows come from a later run, Proj_18143 on
v2.x-pre2 (15 hg38 samples, postProcess BAMs 0.53 to 5.2 GB, 35.0 GB in
all, control job 18983501, 2026-10-12). Count is featureCounts v2.1.1 at
`-T 8`: 104 s for the 35.0 GB, against 35 min for v1.6.0 at `-T 10` on
Proj_18143_B. `RATE_COUNT` is twice the measured rate. CountCheck is
Check B, `bin/assignedCheck.sh`, one job per sample.

With these rates only the 14.1 GB sample's POST job of that run would
have been `LONG` (estimate 85 min; it took 64). POST is single-threaded
and I/O bound (TotalCPU close to Elapsed; it read 120 GB and wrote 102 GB
for the largest sample), so its time depends on shared filesystem load
as much as on the node.

CPU requests are larger than the use. TotalCPU/Elapsed: CALLP about 1.2
of 3 cores, TSSE 1.0 of 2, MergePeaks about 2 of 3 (00.ISSUES.md). Count
with featureCounts v2.1.1 used 6.7 of its 8 cores on Proj_18143, so
`-n 8` stays.

## Memory

`sacct` MaxRSS does not show how much memory a job's processes use.
IRIS gathers it with `jobacct_gather/cgroup` on cgroup v1, so it
includes the file cache of everything the job read and wrote. The kernel
frees that cache when the job reaches its cap, so a job reports MaxRSS at
its cap and still completes. On Proj_18143, 9 of 15 POST jobs reported
32.0G against 32G, the Index jobs on one node 4.0G against 4G
(`samtools index`), and the ED6-2 CountCheck 8.0G against 8G. None was
killed. A job really out of memory ends in State `OUT_OF_MEMORY`, which
`bCheck` and `bin/checkRun.sh` report.

The memory the processes allocate (`total_rss` in the job's
`memory.stat`, without cache) was measured on 2026-10-12 by sampling
the job cgroup every 5 s. Inputs were the largest Proj_18143 sample,
ED6-1: 49M fragments, 5.2 GB postProcess BAM. The probe script is
`/scratch/core001/bic/soccin/AFTER_ACTION_261011/p18143_resources/memProbe.sh`.

| Stage | Cap | Peak without cache | Peak cache | Bounded by |
| --- | --- | --- | --- | --- |
| POST | 32G | 24.1G (SortSam) | 26.3G | the picard heap, `-Xmx24g` in `~/bin/picard`; CollectInsertSizeMetrics used 0.8G |
| CountCheck | 12G (was 8G) | 4.0G | 5.8G | `sort -S 4G` |

In both, the process memory is set by a fixed limit, not by the input:
sort spills to disk once its buffer is full, and the JVM never exceeds
its heap. A larger sample therefore needs no larger cap. CountCheck was
raised to 12G anyway so that its sacct MaxRSS no longer sits at the cap.
For BW (`sort -S20g`, 24G) and MergePeaks (`sort -S 16g`, 24G) the same
bound applies. Count (featureCounts), TSSE, Index, DESEQ and the control
job reported at most 16.2G, 4.8G, 4.0G, 1.1G and 1.4G with cache, well
under their caps. If `~/bin/picard` ever gets a larger `-Xmx`, raise the
POST cap with it.

Checked with test24 (control job 19151733, 2026-10-12, the 11
downsampled BAMs at `-q 20`): `COMPLETED`, 69 of 69 jobs, every
CountCheck job ran with 12G, Count was `SHORT` at 0.1 min/GB, and the
count matrix and count-check report are identical to test23.

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
When the exit came from `bCheck`, the trap also puts the reason in
`MESSAGE`: `bCheck` leaves the `sacct` rows of the failed jobs in
`BCHECK_FAILED`, and `failedJobsMessage` (`bin/slurmTools.sh`) looks up
each job's log in `jobs.tsv` and takes its first line with `FATAL`, else
its first with `error` or `fail` in any case. That also catches the
`error: *** JOB ... DUE TO TIME LIMIT ***` line Slurm writes into the
log of a job it kills (probe job 18919113, `TIMEOUT`); an out-of-memory
kill was not tried. A log with no such line gives the `sacct` state
alone. At most 3 jobs are listed. None of the successful test17 stage
logs has a line that matches.

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
directory). featureCounts writes its temp files next to its output
unless told otherwise; since v2.0.0 `bin/countPeaks.sh` passes
`--tmpDir "$TMPDIR"`.

## Exit-status fixes made along the way

`postMapBamProcessing_ATACSeq.sh` backgrounded `CollectInsertSizeMetrics`
and never waited for it. Under LSF the BED pipeline usually outlasted it;
under Slurm the job ends when the script exits and picard would be killed
mid-write. It now `wait`s, runs under `pipefail`, and exits non-zero on any
failed step so `bCheck` sees it. `makeBigWigFromBEDZ.sh` and
`mergePeaksToSAF.sh` run under `pipefail` for the same reason.

Added in v1.5.1: `callPeaks_ATACSeq.sh` runs under `pipefail` and stops
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
run from a copy of the checkout without
`R/TSSEnrich/lib/b38_tss.bed`, standing in for mm10, with `-q 20`
(test16, 18460600) completed with 47 of 47 jobs,
no TSSE jobs, the warning in `MESSAGE`, and `featureCounts -Q 20`.
`pipe.sh` run by hand on a header-only BAM with an unrecognized `@SQ`
set, alone and next to a b38 BAM, stopped with rc 1 and the reason in
`MESSAGE`, before submitting anything.

### featureCounts v2.1.1 and the count checks (v2.0.0)

Done 2026-10-11 on `fix/featurecounts` with the same 11 downsampled BAMs
(`-q 20`). Test scripts and outputs are in
`/scratch/core001/bic/soccin/AFTER_ACTION_261011/v200_tests/`; the
reference matrix is the v1.6.0 `-T 1` count of the test16 postProcess
BAMs (`repro_T1_cwd/rep1.txt`).

| Test | Result |
| --- | --- |
| test17, normal run (control job 18904482) | `COMPLETED`, 69 of 69 jobs; Checks A and B pass for all 11 samples (ratio 1.0003 to 1.0005, Check B diff 0); count matrix identical to the reference, 0 differing cells. Count 15 s, CountCheck 4 to 6 s, all on `cmobic_short` |
| test18, `--countReadPairs` removed from a copy of the checkout (18904486) | `STATUS=FAILED`, `STAGE=Count`, rc 1; Check A ratio 2.0000 for every sample; the Count job's sacct state `FAILED` |
| test19, `b38_tss.bed` removed from a copy, standing in for mm10 (18904488) | `COMPLETED`, no TSSE jobs, the TSS warning in `MESSAGE`; matrix identical to the reference |
| `checkCounts.sh` on the broken v1.6.0 test16 summary | rc 1; Check A EV_3 ratio 2.0000, Check B EV_3 Assigned 400,560 against 200,961 |
| `checkCounts.sh` on a correct v2.1.1 summary | rc 0 |
| `checkCounts.sh` with one sample's Check B result left out | rc 1, sample named |
| Repeat test, v2.1.1 with the Count options, 20 runs each at `-T 2, 4, 8, 16, 32`, temp dir on `/localscratch` and on `/scratch` | 0 of 200 bad: every run exited 0, passed Check A, and gave a matrix identical to the reference (v1.6.0 at `-T 10` failed 6 of 20 with the shared and 19 of 20 with the local temp dir) |

Check B was calibrated before the code was written: `samtools view -L`
distinct read names over the atlas equalled featureCounts `Assigned`
exactly on all 11 downsampled and all 12 full-size Proj_18143_B BAMs,
as did `bedtools intersect -u -split`, which took about 2.2 times as
long.

In test18 `MESSAGE` was empty; the Check A lines were only in the Count
job log (`00.ISSUES.md` #7, fixed on `fix/message`, below).

### `MESSAGE` on a failed stage (#7)

Done 2026-10-11 on `fix/message`, same BAMs and `-q 20`. Scripts in
`/scratch/core001/bic/soccin/AFTER_ACTION_261011/msg_tests/`.

| Test | Result |
| --- | --- |
| test20, `--countReadPairs` removed from a copy (control job 18918854) | `FAILED`, `STAGE=Count`, rc 1; `MESSAGE=Count failed: all job 18920464 FAILED: FATAL ERROR: Check A failed on peaks_raw_fcCounts.txt.summary`, shown by `checkRun.sh` |
| test21, normal run (18918899) | `COMPLETED`, 69 of 69 jobs, `MESSAGE` empty |
| test22, as test20 and without `b38_tss.bed` (18918900) | `FAILED` at Count; the same message, then `; WARNING: TSS enrichment skipped; no R/TSSEnrich/lib files for b38` |
| probe job 18919113, `-W 0:01:00` running `sleep 600` | `TIMEOUT`; `logReason` returns Slurm's `error: *** JOB ... DUE TO TIME LIMIT ***` line |
| `failedJobsMessage` by hand: 5 failed jobs, missing log, job not in `jobs.tsv`, missing accounting records | 3 jobs listed then `2 more`; state only; `?` for the sample; the `bCheck` reason |

### Genome tags, mitochondrial reads, mm10 TSS (v2.0.0)

Done 2026-10-11 on `fix/genome-tags` (`docs/GENOMES.md`), same BAMs and
`-q 20`. Scripts in
`/scratch/core001/bic/soccin/AFTER_ACTION_261011/genome_probe/`.

| Test | Result |
| --- | --- |
| test23, normal run (control job 18945280) | `COMPLETED`, 69 of 69 jobs, `GENOME=hg38`. Count matrix, atlas, count checks, scale factors, and every sample's shifted BED, narrowPeak, bigWig and TSS CSV identical to test21 (b38; the hg38 genome file is the b38 one renamed). New `___MITO.txt` for all 11: chrM 6.5 to 11.2% of reads. Shifted BEDs hold chr1-22, X, Y only |
| `diffAnalysisPairwise.R sampleManifest.csv comparisons.csv` in test23 | rc 0, `GENOME = hg38` from the status file; no peak at FDR 0.05 on these BAMs, so annotation was tested separately: hg38 and b37 TxDb loaded from `GENOMES.tsv`, 200 of 200 atlas peaks annotated with a symbol (b37 with `chr` added) |
| `diffAnalysisPairwise.R` with no arguments; with no status file; with `GENOME=b38` | rc 1 each, the last two with the reason |
| `pipe.sh` by hand, header-only BAMs from the IRIS dicts | b37 and mm10: detected, pass the genome and TSS file checks, stop at the `@RG SM` check (the header has no `@RG`), 0 jobs. hg19 and hg38+mm39: stop with `unsupported genome [...]; supported: b37 hg38 mm10`. hg38 BAM with a b37 BAM: `BAMs on different genomes` |
| Same, after `bin/getGenomeBuildBAM.sh` was replaced by ChIP-seq's copy (`headerTest2.sh`), plus a missing BAM and a text file named `.bam`, each after an hg38 BAM | detection unchanged for every dict and for a real hg38 BAM; the five cases above stop as before; the two unreadable BAMs stop with `cannot identify the genome of [...]`, 0 jobs |
| `bin/genomeInfo.sh` every tag and field; `b38`, empty and `unknown <md5>` tags; `HOMER` field | values as in the table; rc 1 with the reason for each bad input |
| Stage scripts by hand on one hg38 BAM (EV_1) | POST, BW, CALLP and TSSE rc 0; `___MITO.txt` 39,952 of 488,265 reads, the same as `samtools idxstats`; TSSE `pct_tss_used` 1 |

### Full-size hg38, b37 and the R reports (v2.0.0)

Done 2026-10-12. Scripts in
`/scratch/core001/bic/soccin/AFTER_ACTION_261011/b37_probe/`.

| Test | Result |
| --- | --- |
| Proj_18143 on v2.x-pre2, full size (15 hg38 samples, `Proj_18143/ATACseq`, control job 18983501) | `COMPLETED`, 93 of 93 jobs; Check A 1.0004 to 1.0008, Check B equal for all 15. Against the 2026-03-14 delivery (featureCounts 1.6.0): atlas identical (654,370 peaks), 13 count columns identical, BD6-3 and BRAF-1 about 1.9 times lower (the bug) |
| test25, b37: 4 Proj_18741 BAMs (MCF7-DMSO-1/2, MCF7-LYFUL-1/2) downsampled 100x, MAPQ 10 (control job 19196648) | `COMPLETED`, 27 of 27 jobs, `GENOME=b37` from real BAM headers. `___MITO.txt` contig `MT`, 4.2 to 5.2%, reads and MT reads equal to `samtools idxstats`. Shifted BEDs, narrowPeaks, atlas and bigWigs hold 1-22, X, Y only; bigWig scale from the filtered BED. Checks A and B pass. TSSE 19.0 to 20.6, `pct_tss_used` 1 (full-size delivery 19.7 to 21.3) |
| `diffAnalysisPairwise.R` (`fix/diff-stats`) in test25 with the generated manifest, an `EXC` sample, and `SampleID`s that differ from the `MapID`s | rc 0 for each (the last two stopped with "Count matrix validation failed" before the fix). A manifest missing a sample: rc 1, the sample named |
| Same, on a scratch copy of the full-size Proj_18741 count matrix (b37, 24 samples, 5 comparisons) | rc 0; 14,032 to 100,697 peaks at FDR < 0.05 against 752 to 18,872 in the delivered (Holm) workbook from the same matrix; 99.99% annotated with a gene symbol (`chr` added) |
| `R/analyzeATAC.R` on the same matrix, one group per sample (24 groups) | PCA with 24 distinct colors and the variance explained on the axes (PC1 22%, PC2 15%) |
| test26, the 11 hg38 BAMs with `-q 20` on `fix/diff-stats` (control job 19202406) | `COMPLETED`, 69 of 69 jobs; reports named `ATAC__*.pdf` (no `Proj_` in the path) and accepted by the deliverables check; count matrix, count checks and scale factors identical to test24. `diffAnalysisPairwise.R` rc 0, `GENOME = hg38`, GD1 vs EV 2 peaks at FDR < 0.05, both annotated |

### Read filter: orphans out, no insert-size filter (v2.0.0)

Done 2026-10-12 on `fix/read-filter` (C-06 and C-07 of
`ENCODE_PIPE_CONCORDANCE.md`), against test39 and test40. Script
`/scratch/core001/bic/soccin/AFTER_ACTION_261011/read_filter/compareRuns.sh`
compares two runs; per sample it counts the postProcess reads, the
reads whose name is not there twice (orphans) and the reads with
`|TLEN| <= 30`.

| Test | Result |
| --- | --- |
| test43, the 11 hg38 BAMs with `-q 20` (control job 19283036) | `COMPLETED`, 69 of 69 jobs; Checks A and B pass, Check A exactly 1.0000 for every sample (test39 1.0003-1.0005). 0 orphans (test39 116-305 per sample); 1,126-2,164 reads with `|TLEN| <= 30` kept per sample; postProcess BAMs 0.2-0.4% larger. Filtered peaks 150,408-266,678; 9 samples gained 136-477 and GD1_1 and GD1_2 lost 2,352 and 1,510. Both losses are in MACS2's own output and among peaks with a pileup of 1 (GD1_1 146,348 to 143,782, GD1_2 103,839 to 102,023), which pass `-p 0.01` only just in 100x downsampled data; peaks with a pileup of 2 or more increased slightly. Atlas 902,445 peaks (test39 901,681) |
| test44, the 4 b37 BAMs, MAPQ 10 (control job 19283038) | `COMPLETED`, 27 of 27 jobs; Checks A and B pass, Check A exactly 1.0000. 0 orphans (test40 302-371); 2,188-5,162 reads with `|TLEN| <= 30` kept; filtered peaks within 7 of test40's for the 3 capped samples, 647 more for the fourth. Atlas 727,558 peaks (test40 727,125) |
| POST run time with the name sort and `fixmate` | 25-41 s (hg38) and 43-53 s (b37), against 24-34 s and 33-41 s; to be measured on a full-size run |

### MACS2 `--shift -75 --keep-dup all` (v2.0.0)

Done 2026-10-12 on `fix/macs2-encode`. Scripts in
`/scratch/core001/bic/soccin/AFTER_ACTION_261011/macs2_probe/`
(`compare.sh`, `shiftcheck/run.sh`). FRiP is ENCODE's: reads of the
shifted BED over the sample's own peaks / all its reads. Values of the
earlier run are in parentheses.

| Test | Result |
| --- | --- |
| MACS2 2.2.9.1 on one `+` read with its cut site at 100000 and one `-` read with its cut site at 200000, `-B` pileup | `--shift -75`: 99925-100075 and 199925-200075, centered on both cut sites. `--shift 75`: 100075-100225 and 199775-199925 |
| test27, the 11 hg38 BAMs with `-q 20` (control job 19224611), against test26 | `COMPLETED`, 69 of 69 jobs; Checks A and B pass. Peaks per sample from 7% fewer to 0.4% more (150,451 to 267,009), median width 203-215 bp (208-245), FRiP 0.90-0.95 (0.70-0.78), atlas 903,877 (919,558). Of each sample's top 10,000 peaks, 21-29% have the summit within 500 bp of a TSS (14-28%) |
| test28, the 4 b37 BAMs, MAPQ 10 (control job 19224612), against test25 | `COMPLETED`, 27 of 27 jobs; Checks A and B pass. Peaks 8-17% fewer (265,241 to 366,622), median width 199-209 bp (150-229), FRiP 0.86-0.89 (0.59-0.69), atlas 760,751 (745,067). Top-10,000 summits within 500 bp of a TSS 25-32% (27-31%) |
| `--keep-dup 1` loss, from the MACS2 `.xls` "tags after filtering" | 0.2-0.4% of reads on the downsampled BAMs; 5.0-16.8% on the full-size Proj_18143 run (control job 18983501), the most in the deepest samples |

The downsampled BAMs have 0.3 to 0.9 million reads after filtering, so
their 150,000 to 370,000 peaks at `-p 0.01` are mostly weak; the TSS
figures are the least reliable of these comparisons. The full-size
Proj_18143 run is the test that remains.

### Blacklist filter (v2.0.0)

Done 2026-10-12 on `fix/blacklist`, against test27 and test28 (same
BAMs and MACS2 settings, no blacklist). Scripts
`/scratch/core001/bic/soccin/AFTER_ACTION_261011/macs2_probe/checkBL.sh`
and `checkBL3.sh`.

test29 to test32 ran a first version that wrote the filtered peaks
under MACS2's own names (`*_peaks.narrowPeak`, `*_summits.bed`) and
did not keep MACS2's files. The version kept, tested in test33 and
test34, leaves MACS2's files as written and adds
`*_peaks.blacklistFiltered.narrowPeak` and
`*_summits.blacklistFiltered.bed`. In the test29 to test32 rows,
"filtered narrowPeak" is the file that is now
`*_peaks.blacklistFiltered.narrowPeak`.

| Test | Result |
| --- | --- |
| test29 (hg38, control job 19231490) and test30 (b37, 19231491), first version, peak filter only | `COMPLETED`; filtered peaks correct (as below), but the atlas still overlapped the blacklist at 5 (hg38) and 1 (b37) merged peaks: a 106-383 bp blacklisted region between two kept peaks less than 500 bp apart. Led to the subtract in `mergePeaksToSAF.sh` |
| test31, the 11 hg38 BAMs with `-q 20` (control job 19234652), first version | `COMPLETED`, 69 of 69 jobs; Checks A and B pass. Per sample 264-475 of 150,451-267,009 peaks removed; each filtered narrowPeak is byte-identical to test27's filtered by hand (`intersect -v`, score cap); 0 peaks on the blacklist; one summit per peak; largest score 498. Atlas 901,681 peaks (test27 903,877, of which 2,353 overlapped the blacklist), 0 overlapping it |
| test32, the 4 b37 BAMs, MAPQ 10 (control job 19234654), first version | `COMPLETED`, 27 of 27 jobs; Checks A and B pass. 652-891 of 265,241-366,622 peaks removed per sample, same checks as test31 all pass. Atlas 758,830 (test28 760,751, 1,946 overlapping), 0 overlapping |
| test33 (hg38, as test31, control job 19248967) and test34 (b37, as test32, 19248968), version kept | `COMPLETED`, 69 and 27 jobs; Checks A and B pass. For all 15 samples MACS2's `*_peaks.narrowPeak`, `*_summits.bed` and `*_peaks.xls` (body) are byte-identical to test27/test28, and the `blacklistFiltered` narrowPeak and summits to test31/test32's filtered files; 0 filtered peaks on the blacklist; all five files in `atacSeq/macs/`. Atlas, count matrix, `.summary` and scale factors identical to test31/test32 |
| test35 (hg38, as test33, control job 19253399) and test36 (b37, as test34, 19253400), `fix/blacklist` at `b0cda24`, run to look over the delivered layout | `COMPLETED`, 69 and 27 jobs; Checks A and B pass; same checks and results as test33/test34. No narrowPeak has a negative start or a missing summit (`-1`) |
| test37 (hg38, as test35, control job 19260883) and test38 (b37, as test36, 19260884), `encode-2026` at `bb486a4`: adds the five-file deliverables check (`507116c`) and the one-`sort` atlas merge (`533308e`) | `COMPLETED`, 69 and 27 jobs; Checks A and B pass; same checks and results as test33/test34 |
| Counts in the atlas rows that overlapped the blacklist | test27 0.28% of all counts, 2 of the 10 and 5 of the 100 rows with the most counts; test28 0.29%, 1 and 4; full-size Proj_18143 (v2.x-pre2) 1,376 of 654,370 rows, 0.43% of counts, 3 of the top 10, 7 of the top 100 |

Not run: an mm10 pipeline (`00.ISSUES.md` #14). The mm10 blacklist was
checked only by `getBlacklists.sh` (md5, contigs in `mm10.genome`).

### Peak cap (v2.0.0)

Done 2026-10-12 on `fix/peak-cap` (C-04 of
`ENCODE_PIPE_CONCORDANCE.md`), against test37 and test38. Script
`/scratch/core001/bic/soccin/AFTER_ACTION_261011/peak_cap/compareWithEncode.sh`
applies ENCODE's own commands (`encode_task_macs2_atac.py`,
`bed_clip`, `encode_lib_blacklist_filter.py`) to each sample's MACS2
narrowPeak and compares the result with the pipeline's
`*_peaks.blacklistFiltered.narrowPeak`, every column but the peak
name.

| Test | Result |
| --- | --- |
| test39, the 11 hg38 BAMs with `-q 20` (control job 19268482) | `COMPLETED`, 69 of 69 jobs; Checks A and B pass. No sample over 300,000 peaks (150,451-267,009); every `blacklistFiltered` file equal to ENCODE's commands; MACS2's files and summits equal to test37's; atlas, `.summary` and scale factors byte-identical to test37's |
| test40, the 4 b37 BAMs, MAPQ 10 (control job 19268483) | `COMPLETED`, 27 of 27 jobs; Checks A and B pass. 3 samples capped at 300,000 (315,877, 366,622, 339,604 from MACS2), 299,207-299,329 kept after the blacklist; every `blacklistFiltered` file equal to ENCODE's commands; MACS2's files equal to test38's. Atlas 727,125 peaks (test38 758,830), 0 on the blacklist |
| `bedClip` | ran for all 15 samples; no peak needed clipping |
| The block itself, on crafted files with the cap set to 3 | sort, cap, p-value tie broken by the whole line, clamped start, summit `-1`, end clipped, blacklist, score cap and summits as expected; a peak wholly past its contig end and a contig not in `<TAG>.genome` each stop the job |

### 1-based atlas start (v2.0.0)

Done 2026-10-12 on `fix/saf-start` (C-17 of
`ENCODE_PIPE_CONCORDANCE.md`), against test39 and test40. Script
`/scratch/core001/bic/soccin/AFTER_ACTION_261011/read_filter/compareRuns.sh`
compares two runs: BAMs, shifted BEDs, filtered peaks, the atlas and
every cell of the count matrix.

| Test | Result |
| --- | --- |
| test41, the 11 hg38 BAMs with `-q 20` (control job 19282732) | `COMPLETED`, 69 of 69 jobs; Checks A and B pass. BAMs, shifted BEDs and filtered peaks identical to test39's; atlas equal to test39's with every start plus 1 (901,681 peaks); 55 of 9,918,491 count cells differ, each by -1: a fragment that overlapped only the base left of a peak |
| test42, the 4 b37 BAMs, MAPQ 10 (control job 19282733) | `COMPLETED`, 27 of 27 jobs; Checks A and B pass. Same as test40 but for the atlas starts (727,125 peaks, each start plus 1); 112 of 2,908,500 count cells differ, each by -1 |
| test45, the 11 hg38 BAMs with `-q 20`, on `encode-2026` at `ba329b6` with this fix and the read filter (control job 19287472) | `COMPLETED`, 69 of 69 jobs; Checks A and B pass, Check A 1.0000. BAMs, shifted BEDs and filtered peaks identical to test43's; atlas equal to test43's with every start plus 1 (902,445 peaks); 55 of 9,926,895 count cells differ, each by -1 |
| test46, the 4 b37 BAMs, MAPQ 10, same code (control job 19287473) | `COMPLETED`, 27 of 27 jobs; Checks A and B pass, Check A 1.0000. Identical to test44 but for the atlas starts (727,558 peaks, each start plus 1); 109 of 2,910,232 count cells differ, each by -1 |

## Not done

- `deliverResults.sh` still points at `/ifs/res/seq/pi/invest` and
  `~/Code/BIC/Delivery`; not a scheduler change and not touched. It is
  `00.ISSUES.md` #1.
- The `RATE_*` values come from one project on one genome. Check the
  `runClass` lines in the control log against `sacct` Elapsed on new
  projects, and raise a rate if a job gets close to its walltime.
