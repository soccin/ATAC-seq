# ATAC-seq — Closed issues

Record of issues removed from `00.ISSUES.md` once addressed. Not a
tracking document; open issues live only in `00.ISSUES.md`.

## Closed — fixed on `feat/runstatus` (part of former #15)

**#15 (partial)** `00.POST_RUN.txt` contained a literal `$SDIR`. Its
notes are now part of `00.RUNSTATUS.txt`, with `$SDIR` expanded.

## Closed — fixed on `fix/resources` (not yet merged)

Measured on the first full-size run, `Proj_18143_B` (12 b38 samples,
input BAMs 4.6 to 14.1 GB, control job 18335688, 2026-10-09). Figures and
reasoning in `docs/SLURM_PORT.md`, "Walltime classes".

**#2 Partition and qos are not chosen per stage** — `atacSub` takes a
class, `SHORT` (`-W 1:59:00`, no qos) or `LONG` (`-W 12:00:00`,
`SBATCH_QOS=priority` on that `bsub` call only). `pipe.sh` unsets every
`SBATCH_*` variable it inherits. `CLAUDE.md` and `docs/SLURM_PORT.md` no
longer say to export `SBATCH_QOS`. Verified: a test run submitted from a
shell with `SBATCH_QOS=priority` exported ran every stage job on
`cmobic_short` at qos `normal`; a `LONG` probe job (18373842) ran on
`cmobic_cpu` at qos `priority`.

**#3 Stage walltimes dropped from 6 h to 1h59m** — Largest Elapsed on the
full-size run: POST 1:04, Count 0:35, BW 0:27, TSSE 0:10, CALLP 0:09; no
`TIMEOUT`. POST, BW, CALLP, TSSE and Count now pick their class per job
from the input size (`runClass` and the `RATE_*` values in `pipe.sh`): a
job is `LONG` if the estimate exceeds 60 min. On that run only the
14.1 GB sample's POST job would have been `LONG`.

**#4 Memory caps are estimates** — No job reached `OUT_OF_MEMORY`. POST
and Count report `MaxRSS` equal to their caps (32G, 24G), which is page
cache, not a sign of pressure (see the original note: do not size from
`MaxRSS`). Caps left unchanged; raise one if a job ends `OUT_OF_MEMORY`.

## Closed — fixed on `feat/runstatus` (merged into `devs`)

**#17 Cancelling the control job leaves its stage jobs running** — The
cause was not what this entry said. On `scancel` Slurm sends SIGTERM to
the batch shell only, and bash runs a trap only after its foreground
child exits, so with `bSync.sh` in the foreground the trap never ran
before the SIGKILL 30 s later (test job 18357259). `bSync` now runs
`bSync.sh` in the background and `wait`s, which a trapped signal
interrupts; `onExit` cancels the run's jobs unless `RUN_DONE=1`, whatever
`$?` is. Verified with job 18357771: trap ran at once, `STATUS=FAILED
rc=143`, both stage jobs `CANCELLED`. During the R reports at the end of
`pipe.sh` the trap still may not run (#26).

## Closed as by-design — do not "fix"

**#18 Dependence on `~/bin` and `~/.Rprofile`** — `sbatch`, `bsub`,
`bSync.sh`, `picard` and `bedtools` come from `~/bin`, and the R scripts
call `len()`, `cc()` and `DATE()` from `~/.Rprofile`. This was the user's
decision for the port. Consequence: only an account with these files can
run the pipeline.

**#19 Control job, not a dependency graph** — The user's decision. The
shim cannot express name-based `afterok` (Slurm needs numeric ids), so
`bSync`/`bCheck` sequencing replaces `-w post_done()`.

**#20 MD5 handshake with a 300 s retry** — Stays, per `CLAUDE.md`. The
consuming stages now start only after `bSync`, so the retry should rarely
trigger; a real mismatch costs five minutes.

**#21 One run per analysis directory** — Job names embed the control
job's pid, so runs in separate directories do not interfere. Two runs in
one directory share `out/` and collide.

**#22 Preemption** — Not a concern. `cmobic_cpu`, `cmobic_short` and
`cpushort` all have `PreemptMode=OFF` (checked 2026-10-09), so a job is
not requeued and restarted from the top.

## Closed — port work

Control job on `cmobic_cpu` with `--qos=priority`, stages on the short
partitions through `~/bin/bsub` (`37a66ee`) · `bin/slurmTools.sh`:
`bSync`, `bCheck` on `sacct` State, `bCheckAll`, `bKill` · `-w post_done()`
replaced by `bSync` sequencing · `bin/loadTools.sh` (`samtools/1.20`,
`~/bin/bedtools`) · `picardV2` → `~/bin/picard` · scratch moved to
`/scratch/core001/bic/$USER/ATACSeq` · `postMapBamProcessing` waits for
`CollectInsertSizeMetrics` and runs under `pipefail` · venv from the
`python3` on PATH, unused `ChIPseeker` dropped (`8ed8a67`) · docs and
`docs/SLURM_PORT.md` (`cf43bf7`) · v1.5.0 release docs (`91fa54e`) ·
outside the repo: `~/bin/sbatch` no longer `exec`s, so its temp copy is
removed.
