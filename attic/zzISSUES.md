# ATAC-seq — Closed issues

Record of issues removed from the open list once addressed. Not a
tracking document; open issues live only in `00.FIXME.md`. Until
2026-10-12 the open list was `00.ISSUES.md` (now
`attic/00.ISSUES.md`), with the ENCODE review's findings (C-NN) in
`attic/ENCODE_PIPE_CONCORDANCE.md`.

Each entry keeps the number the issue had when it was closed. The open
issues were renumbered on 2026-10-10 (old to new: #8 to #1, #9 to #2,
#11 to #3, #12 to #4, #13 to #5, #15 to #6, #25 to #7, #26 to #8, #27 to
#9, #28 to #10, #29 to #11, #30 to #12, #31 to #13), so a number below
refers to the old list unless it says otherwise. Entries from
`00.FIXME.md` say so.

## Closed — fixed on `feat/outdirs-01`

**#11 of `00.FIXME.md`, `makeBigWigFromBEDZ.sh` scale-factor branch
reads the wrong argument** (`bug` `pre-port`) — The script now takes
`GENOME BEDZ OUTDIR [SCALEFACTOR]` and the branch reads `$4`. Still
unused: `pipe.sh` passes three arguments.

**#12 of `00.FIXME.md`, first part: the manifest step leaves `mapid`,
`sid` and `gid`** (`cleanup` `pre-port`) — `pipe.sh` now builds
`sampleManifest.csv` in one pipeline, with no intermediate files; the
rows are the same (checked on test IDs). The `--dirty` part stays open
as #12.

**#16 of `00.FIXME.md`, part: the count checks and the scale factors
are not staged** — Both are now written to `atacSeq/atlas` and
delivered, with sections in `docs/output.md`. The rest of #16 stays
open.

## Closed — fixed on `rel/2.0.0`

**#28 of `00.FIXME.md`, atlas merges peaks within 500 bp** (was
C-18; `docs` `encode`) — The merge itself stays: `mergePeaksToSAF.sh`
uses `bedtools merge -d 500`, R. Koche's method, where ENCODE never
merges, so the atlas grows with the number of samples. The one action
was to note in deliveries that ENCODE's per-sample peak-count
standards do not apply to the atlas; `docs/output.md` ("Atlas") and
`docs/methods.md` ("Peak atlas") say so since v2.0.0.

**#2 of `00.FIXME.md`, `deliverResults.sh` not ported**
(`enhancement`) — Rewritten on `feat/deliver-v2` after adagio's
`bin/deliver.sh`: it delivers to `/data1/core002/res/bic/results` as
`svc_core001_bot01`, picks the next `r_NNN`, stops unless
`bin/checkRun.sh` exits 0, delivers diff files with any RUNTAG and
the client docs `output` and `methods` as PDF and HTML (in
`atacSeq/docs`), and has a dry run (`-n`). The email still names the
untagged diff files; the script prints a `WARNING` to edit it when
there are tagged ones.
Tested with `-n` only (`CHANGELOG.md`, "Not tested").

**#1 of `00.FIXME.md`, no full-size run of v2.0.0 with the 5' clip
filter** (`risk`) — Proj_18143 (15 hg38 samples, MAPQ 10) was run as
`Proj_18143/ATACseq/v2.x-pre5` on tag `v2.x-pre5` (2026-10-13,
control job 19429072) and checked against v2.x-pre4. Every criterion
held: `COMPLETED`; Checks A and B pass, Check A 1.0000; 0.68-1.37% of
the postProcess reads removed per sample and 0 reads left with a 5'
clip or `|TLEN| <= 30` in all 15 BAMs; `MIN_INSERT_SIZE` 31 and the
nucleosome peaks in every `___INS.txt`; the only cut sites with 100 or
more tags lie in the blacklisted `chr1:628903-635104`, and the named
hotspots and their atlas rows are gone from ED1-1's top 20 counts;
TSS enrichment within -0.02 to +0.08 of pre4. No POST job came near
the 1:59 limit (the longest 1:18:32), so `RATE_POST` stays 6. Recorded
in `docs/SLURM_PORT.md`, "Insert-size and 5' clip filters".

**`MUST_BE_DONE_FIRST.md`, the v2.0.0 pre-release review** — A line
by line review of the code commits of 2026-10-12 (`678c04a`,
`5b92c57`, `b391277`, `424d252`, `1f21b5b`, `6eef929`, `5bd921d`,
`507116c`, `3bbbcb3`, `533308e`, `e75a8b7`, `e50dec6`, `c8d17d7`),
against ENCODE and against each other, found no bug that affects the
atlas or the bigWigs. The two code commits of 2026-10-13, `2d469cf`
(insert filter, Picard `W=1000`) and `6549328` (`KEEP_5P`), were
reviewed the same way, with the same result. Its two release items
are done: the full-size run (#1 above; the run times and memory of
every stage are in `docs/SLURM_PORT.md`, "Walltime classes"), and
the `docs/GENOMES.md` MITO wording, which became correct again when
the insert filter came back. Its two non-blocking notes are
`00.FIXME.md` #30 and #31.

**#25 of `00.FIXME.md`, documentation out of step with the code**
(C-22) — `NOTES.md` has a dated note that MACS2 runs on the shifted
BED (`-f BED`), not the BAM. `README.md` says the input BAMs must
have duplicates marked by Picard MarkDuplicates in PEMapper, not by
the mapper. `R/TSSEnrich/README.md` gives ENCODE's thresholds per
genome instead of one rule with a "< 1 failed" cut (the README part
of #23).

## Closed — fixed on `fix/read-filter`

Never numbered here: findings C-06 and C-07 of the ENCODE review.

**Orphan reads kept; pairs with an insert of 30 bp or less dropped**
(`pre-port`) — POST now name-sorts the filtered reads and runs
`samtools fixmate -r` and the flag filter again, as ENCODE does, so a
read whose mate failed the MAPQ cutoff is removed; the `|TLEN| <= 30`
filter, which ENCODE does not have, is gone. Check A now gives exactly
1.0000. Tested in test43 (hg38), test44 (b37), test45 and test46
(`docs/SLURM_PORT.md`, "Validation").

## Closed — fixed on `fix/saf-start`

Never numbered here: finding C-17 of the ENCODE review.

**Atlas `Start` 0-based in a 1-based SAF** (`pre-port`) —
`mergePeaksToSAF.sh` now writes the BED start plus 1. Tested in test41
(hg38) and test42 (b37): only the atlas starts and 55 and 112 count
cells (each -1) changed.

## Closed — fixed on `fix/peak-cap`

Never numbered here: finding C-04 of the ENCODE review.

**No peak cap or contig clipping** — `callPeaks_ATACSeq.sh` now sorts
MACS2's peaks by p-value, keeps the top 300,000, clamps negative
starts, clips to the contig ends with `bin/bedClip` (UCSC v369, as
ENCODE) and then applies the blacklist filter, as
`encode_task_macs2_atac.py` does. Each filtered file equals ENCODE's
commands applied to the same MACS2 output, in test39 (hg38) and test40
(b37).

## Closed — fixed on `fix/atlas-sort`

Never numbered here: found in a review of `fix/blacklist`.

**`find | xargs sort` in `mergePeaksToSAF.sh`** — With enough samples
xargs ran several sorts and `bedtools merge` failed on the unsorted
input; a space in a sample name failed the run. Now one `sort
--files0-from` reads the names from `find -print0`. Tested in test37
(hg38) and test38 (b37).

## Closed — fixed on `fix/blacklist`

Never numbered here: finding C-03 of `ENCODE_PIPE_CONCORDANCE.md`.

**No blacklist filtering of peaks** (`pre-port`) — No peak set was
filtered against a blacklist. `callPeaks_ATACSeq.sh` now writes
`*_peaks.blacklistFiltered.narrowPeak` beside MACS2's own (unchanged)
files: the peaks not overlapping ENCODE's blacklist
(`lib/genomes/<TAG>.blacklist.bed`, built by
`lib/genomes/getBlacklists.sh`), score capped at 1000, as
`encode_lib_blacklist_filter.py` does. The atlas is built from those
only, and `mergePeaksToSAF.sh` subtracts the blacklist from it, since
a merge can otherwise bridge a blacklisted region; `pipe.sh` checks
the file exists. A first version wrote the filtered peaks under
MACS2's own file name and lost MACS2's output; never overwrite a
primary tool's files. Tested in test33 (hg38) and test34 (b37) against
test27, test28, test31 and test32
(`docs/SLURM_PORT.md`, "Validation"; `docs/GENOMES.md`, "Blacklist").

## Closed — fixed on `fix/macs2-encode`

Never numbered here: tracked in `BUGFIXES_v2.md` (deleted once fixed)
and the external ENCODE concordance review of 2026-10-12 (C-01, C-02).

**MACS2 `--shift 75` and the default `--keep-dup 1`** (`pre-port`) —
`callPeaks_ATACSeq.sh` copied its settings from the Kundaje lab pipeline
but passed `+75` where it passes `-75`, and left out `--keep-dup all`.
Now `--shift -75 --extsize 150 --keep-dup all`, as in ENCODE's
`encode_task_macs2_atac.py`. Shift direction checked on a two-read
input; tested in test27 (hg38) and test28 (b37) against test26 and
test25 (`docs/SLURM_PORT.md`, "Validation"; `CHANGELOG.md`, v2.0.0).

## Closed — fixed on `fix/diff-stats`

Numbers here are from the list as renumbered on 2026-10-10.

**#19 `diffAnalysisPairwise.R` output names outside a `Proj_`
directory** — `plotINSStats.R`, `R/analyzeATAC.R` and
`R/diffAnalysisPairwise.R` now use the last `Proj_*` or `B-<digits>`
directory of the path, else `ATAC`. Found during the pre-release review
on 2026-10-12, together with defects that were never issues here: the
diff `FDR` was Holm (`p.adjust()` with no method), the QC PCA plotted
the shared profile as PC1, and an `EXC` group or edited `SampleID`s
stopped the diff analysis (`CHANGELOG.md`, v2.0.0). Tested in test25
(b37), on a scratch copy of the Proj_18741 count matrix, and in test26
(hg38, `docs/SLURM_PORT.md`, "Validation").

## Closed — fixed on `fix/resources-02`

Numbers here are from the list as renumbered on 2026-10-10.

**#15 Count and CountCheck rates not measured in jobs** — Measured on
the Proj_18143 run on v2.x-pre2 (15 hg38 samples, 35.0 GB of
postProcess BAMs, control job 18983501, 2026-10-12). Count with
featureCounts v2.1.1 took 104 s, 0.05 min/GB: `RATE_COUNT` 1.2 to 0.1.
CountCheck took 0.79 to 1.02 min/GB: `RATE_CCHECK` stays 1.5. sacct
MaxRSS includes file cache on IRIS, so the process memory was measured
from the job cgroup on ED6-1, the largest sample: CountCheck 4.0G (set
by `sort -S 4G`), POST 24.1G (the picard `-Xmx24g` heap). The
CountCheck cap went from 8G to 12G, because ED6-2 had reported 8.0G
MaxRSS against 8G. Details in `docs/SLURM_PORT.md`, "Memory". The Count
part of #13 also closed: featureCounts used 6.7 of 8 cores, so `-n 8`
stays.

## Closed — fixed on `fix/message`

Numbers here are from the list as renumbered on 2026-10-10.

**#7 `MESSAGE` is empty when a stage fails** (was #25) — When `bCheck`
aborted the run, `00.RUNSTATUS.txt` had `STAGE` and `RC` but no reason;
in test18 the Check A failure was only in the Count job log. `bCheck`
now leaves the failed jobs' `sacct` rows (or, when accounting records
never arrived, a one-line reason) for the `EXIT` trap, which calls
`failedJobsMessage` in `bin/slurmTools.sh`. For each failed job, at
most 3, `MESSAGE` gets the sample, job id, state and the first `FATAL`
line of its log (else the first with `error` or `fail`, any case). An
mm10 TSS warning follows it. The STAGING sweep (`bCheckAll`) also sets
`MESSAGE`. Verified 2026-10-11: test20 (18918854) and test22 (18918900)
failed at Count with the Check A line in `MESSAGE`, test22 followed by
the TSS warning; test21 (18918899) completed. Details in
`docs/SLURM_PORT.md`, "Validation".

## Closed — fixed on `fix/featurecounts`

**featureCounts v1.6.0 silently counts one BAM by read** (never in
`00.ISSUES.md`; found 2026-10-11 in the after-action review of the
v1.5.1 run of Proj_18143_B, brief in
`PortATACseq/ATACseq-BUGFIX-FeatCount.md`) · `bug` `critical` `pre-port`
— Given several coordinate-sorted paired-end BAMs and more than one
thread, the vendored featureCounts v1.6.0 sometimes did not pair the
reads of one input file and counted every read as a fragment. That
sample's column was about 2x too high (EV_2 in Proj_18143_B: 66,852,385
fragments for 33,426,192 pairs), and so were its DESeq2 scale factor
and its place in the QC PCA. featureCounts exited 0 and nothing checked
the totals. On the downsampled test16 BAMs it failed in 6 of 20 runs at
`-T 10` and 19 of 20 with a node-local temp dir. Found in 5 projects run
since 2022 (listed in the brief, which now holds only the open
items); remediating them is separate from this fix.

Fixed by:
- featureCounts v2.1.1 with `-p --countReadPairs --tmpDir "$TMPDIR"`,
  in the new `bin/countPeaks.sh`. 0 failures in 50 runs
  (`AFTER_ACTION_261011/repro_v211_T10_local/`) plus 200
  on the downsampled BAMs at `-T 2/4/8/16/32` with local and shared temp
  dirs; every matrix identical to the v1.6.0 `-T 1` reference.
- Check A in the Count job: fragment total / (mapped reads / 2) within
  0.95 to 1.05.
- Check B: one `CountCheck` job per sample counts the fragments over the
  atlas with samtools alone; it must equal `Assigned` exactly. Exact on
  12 full-size and 11 downsampled BAMs.

Verified 2026-10-11 on the downsampled BAMs (`docs/SLURM_PORT.md`,
"Validation"): the broken test16 summary fails both checks with EV_3
named; a run without `--countReadPairs` fails at `STAGE=Count`; normal
and no-TSSE runs complete with a matrix identical to the reference.
Not yet run: the full-size rerun of Proj_18143_B.

Also partly addresses #11 (the `.summary` is in the deliverables check)
and #13 (Count asks for 8 cores, not 10); both stay open for the rest.

## Closed — fixed on `fix/genome`

Verified 2026-10-10 on the 11 downsampled b38 BAMs; the jobs and results
are in `docs/SLURM_PORT.md`, "Validation". A normal run (18460598) and a
run without the b38 TSS file, standing in for mm10, with `-q 20`
(18460600) both completed; the second had no TSSE jobs, the warning in
`MESSAGE`, and `featureCounts -Q 20`. Mixed and unrecognized genomes
stopped before any job was submitted.

**#7 Unsupported genomes are not rejected up front** — `pipe.sh` now
runs `getGenomeBuildBAM.sh` on every BAM before it submits anything and
stops, with the reason in `MESSAGE`, if the builds differ or the build
is not in `SUPPORTED_GENOMES` (b37, b38, mm10, the builds
`postMapBamProcessing_ATACSeq.sh` accepts). For a build without files in
`R/TSSEnrich/lib/` (mm10) the TSSE stage is skipped with a warning in the
control log and `MESSAGE`, by the user's decision, instead of adding mm10
TSS files; its output is left out of staging and the deliverables check.

**#14 `featureCounts -Q 10` ignores `-q MAPQ`** — The Count stage passes
`-Q $MAPQ`.

## Closed — fixed on `feat/runstatus` (part of former #15)

**#15 (partial)** `00.POST_RUN.txt` contained a literal `$SDIR`. Its
notes are now part of `00.RUNSTATUS.txt`, with `$SDIR` expanded.

## Closed — fixed on `fix/rerun`

**#6 Rerunning in the same analysis directory fails at the end** and
**#23 A second run in the same directory is not refused** — Closed by
not supporting reruns, the user's decision: rerun and resume logic
belongs in a workflow manager, not in bash. `checkPreviousRun` in
`pipe.sh` stops a run before it writes anything if `00.RUNSTATUS.txt`,
`SLURM.CTRL/jobs.tsv`, `out/`, `callpeaks/` or `atacSeq/` exists. A
working rerun was built and tested first (`mkdir -p`, `cp -alf`, since
`cp -al` stops with "File exists" when a stage has replaced a file the
first run hard-linked into `atacSeq/`; records kept as `.prev`) and then
dropped: the stages glob `out/` and `callpeaks/` for their inputs, so a
rerun without one of the earlier BAMs still put that sample into the
atlas and count matrix. The staging `mkdir` calls are now one
`mkdir -p`, and the unused `out/postBams`, `out/metrics` and `out/bed`
are gone.

## Closed — fixed on `fix/sync`

**#5 A transient `squeue` failure aborts the run** — `~/bin/bSync.sh`
reads a failed `squeue` as "all jobs done"; `bCheck` then counted the
still-running jobs as failed after 60 s and the run was cancelled. A
failed `sacct` call ended `pipe.sh` at once under `set -e`. `bCheck NAME
NJOBS` now goes back into `bSync` while any job is active, retries a
failed, empty or short `sacct` answer for up to `ATAC_SACCT_WAIT` (1800
s), and fails only on a final state other than `COMPLETED`. A rejected
`bsub` still aborts the run: no retry, by the user's decision, because a
submission whose reply timed out may have created the job.

**#16 `mkdir -p SLURM.CTRL` is not required** — Slurm 25.11 creates a
missing `-o` directory (probe job 18351604). The step was dropped from
`README.md`, `CLAUDE.md` and the `pipe.sh` header.

**#24 `checkRun.sh` reads a failed `sacct` call as "no record"** — It
now reports `UNKNOWN` and exits 3 when a `sacct` call fails, unless the
status file says `FAILED` or the log trailers alone show success.

## Closed — fixed on `fix/localscratch`

**#10 `sort` spills to the node's `/tmp`** — `makeBigWigFromBEDZ.sh`
(`sort -S20g`) and `mergePeaksToSAF.sh` (`sort -S 16g`) set no `-T` or
`TMPDIR`; IRIS nodes set `TMPDIR=/tmp` (137G, shared by every job on the
node), so sort spill, R and Python temp files went there. `TMPDIR` is
now `/localscratch/$USER`, one directory per stage job, removed at job
end; the sorts pass `-T "$TMPDIR"`. See `docs/SLURM_PORT.md`, "Temp
files: /localscratch, never /tmp".

## Closed — fixed on `fix/resources`

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
`pipe.sh` the trap still may not run (#26, now #8).

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
one directory share `out/` and collide; since #6/#23 a run in a
directory that already holds one is refused.

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
