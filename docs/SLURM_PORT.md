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
| `-W 359` / `-W 59` | `-W 1:59:00` | Every stage on the short partitions. The shim sends anything under two hours to `cmobic_short,cpushort`. `-W` is sbatch format, so `-W 59` means 59 minutes in both, but `-W 5:00` would be five minutes. |
| `-o LSF.01.POST/` | `-o SLURM.01.POST/` | Trailing slash: the shim creates the directory and writes `%j.out`. |
| `-w "post_done(NAME)"` | removed | The shim maps it to `afterok:NAME`, but Slurm needs numeric ids, and an `afterok` whose parent failed pends forever (`kill_invalid_depend` is off). `MergePeaks -> Count -> DESEQ` now runs in sequence with `bSync`. |
| `-q control` | `#SBATCH -p cmobic_cpu` and `#SBATCH --qos=priority` | Control job only; directives in `pipe.sh`. |

Memory per stage: POST 32G (the `~/bin/picard` wrapper runs `-Xmx24g`, so
the cap has to sit above that heap), BW 24G (`sort -S20g`), CALLP 3 cores
18G, MergePeaks 3 cores 24G (`sort -S 16g`), Count 10 cores 24G
(`featureCounts -T 10`), Index 4G, DESEQ 24G, TSS 2 cores 32G.

## Checking a run

Slurm writes nothing into a job's log when it fails. `bCheck` asks `sacct`
for every job with the stage name since `ATAC_START` and aborts on any
State other than `COMPLETED`. Match State, never ExitCode: an OOM kill
reports `OUT_OF_MEMORY` with `ExitCode=0:125`. `pipe.sh` ends with
`bCheckAll "^qATAC-Seq_.*_<pid>$"`, the replacement for the old
`parseLSF.py | fgrep -v Successfully` sweep.

That is enough for the control job to stop on a failure, but not for a
person to find out afterwards whether a run worked. Under LSF every log
ended with `Successfully completed.` or `Exited with exit code N`; under
Slurm a killed control job leaves a log that just stops. `pipe.sh`
therefore keeps three records in the analysis directory, following the
PEMapper port:

| Record | Written by | Content |
| --- | --- | --- |
| `00.RUNSTATUS.txt` | `pipe.sh` | `STATUS=RUNNING`, then `COMPLETED` or `FAILED`; the `STAGE` reached, `RC`, `MESSAGE`, control job id and log, `ATAC_START`, version, genome, samples, BAMs |
| `SLURM.CTRL/jobs.tsv` | `pipe.sh` (`atacSub`) | one line per stage job: job id, stage, sample, log path |
| `#ATAC_EXIT=<rc>` | `bin/runStage.sh` | last line of every stage log; every `bsub` goes through this runner |

`STATUS=COMPLETED` is written only as the last action of `pipe.sh`, after a
check that every sample has its bigWig, peak file, insert-size metrics and
TSS enrichment and that the atlas, count matrix, scale factors and QC
PDFs exist. Any earlier exit runs the `EXIT` trap, which writes
`STATUS=FAILED` with the stage and cancels the rest of the run's jobs.

`bin/checkRun.sh`, run from the analysis directory (or given its path),
reads the three records with `sacct` and exits 0 if the run worked, 1 if
it failed and 2 if it is still running. It does not take `STATUS=RUNNING`
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

## Environment changes

| What | Was | Now |
| --- | --- | --- |
| samtools | `module load samtools` | `module load samtools/1.20` in `bin/loadTools.sh` |
| bedtools | `module load bedtools/2.27.1` | `~/bin/bedtools` (2.31.1) on PATH; no module exists |
| picard | `picardV2` from the cluster PATH | `~/bin/picard` (picard 3.4.0, `-Xmx24g`) |
| scratch | `/scratch/socci/_scratch_ATACSeq` | `${ATAC_SCRATCH_ROOT:-/scratch/core001/bic/$USER/ATACSeq}` |
| venv python | 3.9 required | the `python3` on PATH (3.10); MACS2 2.2.9.1 and IDR 2.0.3 both build |
| R | cluster R | 4.5.1 on PATH; `optparse` must be installed for `tss_enrich.R` |

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
BAMs share one.

## Validation

Eleven b38 BAMs downsampled 100x (`PortATACseq/out/picard/DSB/*/*.dn_100.bam`,
about 1.2M reads each) run from `PortATACseq/test01`. Control job
18237528: 70 stage jobs, all `COMPLETED`, 5m38s wall, every deliverable
present under `atacSeq/`. A first attempt (18236503) failed at DESEQ only
because the test harness overrode `R_LIBS_USER`; `bCheck` caught it and
the `EXIT` trap cancelled the queue, which is the failure path working.

## Not done

- `deliverResults.sh` still points at `/ifs/res/seq/pi/invest` and
  `~/Code/BIC/Delivery`; not a scheduler change and not touched.
- Stage walltimes were not measured on full-size BAMs. If a stage exceeds
  two hours, raise its `-W` past 2:00:00 (the shim then picks `cmobic_cpu`)
  and export `SBATCH_QOS=priority` for that submission.
