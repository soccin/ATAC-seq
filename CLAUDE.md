# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repo is

An ATAC-seq analysis pipeline (bash + R) that runs on the MSKCC IRIS/Slurm
cluster (JUNO/LSF is gone; there is no backward compatibility). It takes
**already-aligned BAMs** as input and produces filtered BAMs, Tn5-shifted
BEDs, bigWigs, MACS2 peak calls, a merged peak atlas with a raw count
matrix, QC metrics, and optional differential peak analysis.

There is no build, no package, and no test suite. Scripts are invoked
directly; changes are validated by running the pipeline on real data.

## Execution model

The repo is a **toolkit directory, not a working directory**. Every script
resolves its own location into `SDIR` and is run from a separate per-project
analysis directory, where all outputs land (`out/`, `callpeaks/`, `SLURM.*/`,
`atacSeq/`). Never assume the repo dir and the cwd are the same.

Install MACS2 + IDR into a local venv once per repo checkout:

```bash
. 00.SETUP.sh          # creates ./venv from the python3 on PATH (3.10 on IRIS)
```

Run the whole pipeline from the analysis directory. `pipe.sh` carries
`#SBATCH` directives for the long partition with `--qos=priority`. Submit it
through `~/bin/sbatch`: Slurm runs a copy of a batch script from its spool
directory, so `$0` is useless, and the wrapper exports `SBATCH_SCRIPT_DIR`,
which `pipe.sh` uses for `SDIR`.

```bash
sbatch /path/to/ATAC-seq/pipe.sh [-q MAPQ] BAM1 [BAM2 ...]
```

`pipe.sh` is the control job: it fans each stage out with `~/bin/bsub` (the
SchedMD openlava shim with local patches), blocks on `bSync JOBNAME`, then
calls `bCheck JOBNAME NJOBS` to abort if any job in the group ended in a
state other than `COMPLETED` in `sacct`. `bSync.sh` reads a failed `squeue`
as "all done", so `bCheck` goes back to `bSync` while `sacct` still shows
a job active, and retries a failed or short `sacct` answer for up to
`ATAC_SACCT_WAIT` seconds (1800). Under `set -e` a failing `$(...)`
assignment ends `pipe.sh`; call Slurm tools in an `if`. Job names embed `$$` so concurrent runs don't
collide; every `sacct` query is bounded by `ATAC_START` because pids
recycle. The helpers are in `bin/slurmTools.sh`; `bSync.sh` itself lives in
`~/bin`. Individual stage scripts can also be run standalone for debugging —
each has its own usage block. See `docs/SLURM_PORT.md` for the flag
translation and the reasoning behind the resource requests.

Whether a run worked is answered by `bin/checkRun.sh` (exit 0 worked, 1
failed, 2 running, 3 unknown because `sacct` failed), run from the
analysis directory. It reads three records
that `pipe.sh` keeps, plus `sacct`:

- `00.RUNSTATUS.txt`: `STATUS=RUNNING|COMPLETED|FAILED`, the `STAGE`
  reached, and the control job id. `COMPLETED` is written only as the last
  line of `pipe.sh`, after the per-sample deliverables check; the `EXIT`
  trap writes `FAILED`. A new control-job step must call `setStage NAME`
  and a new deliverable belongs in that check.
- `SLURM.CTRL/jobs.tsv`: job id, stage, sample and log of every stage job.
  Submit stage jobs with
  `atacSub STAGE SAMPLE LOGDIR CLASS "BSUB_OPTS" CMD ...`
  and wait with `waitStage STAGE`, never with a bare `bsub`, or the job is
  missing from the manifest and its log has no exit trailer.
- `#ATAC_EXIT=<rc>`: the last line of every stage log, written to stderr
  by `bin/runStage.sh`, which `atacSub` puts in front of every command.

One run per analysis directory; reruns are deliberately not supported
(that belongs in a workflow manager, not in bash). Before writing
anything, `pipe.sh` (`checkPreviousRun`) stops if `00.RUNSTATUS.txt`,
`jobs.tsv`, `out/`, `callpeaks/` or `atacSeq/` exists. Do not add rerun
or resume logic; the stages glob `out/` and `callpeaks/` for their
inputs, so an earlier run's samples would be mixed in.

Stage scripts must exit nonzero on any failed step (`pipefail` plus an
explicit check), or `sacct`, `bCheck` and the trailer all report a failed
step as `COMPLETED`.

Every stage job has a walltime class, the `CLASS` argument of `atacSub`,
and `-M` as a hard total-memory cap. `SHORT` is `-W 1:59:00` with no qos;
the shim sends it to `cmobic_short,cpushort`. `LONG` is `-W 12:00:00` with
`SBATCH_QOS=priority` set on that one `bsub` call; the shim sends it to
`cmobic_cpu`. Never give a `SHORT` job a qos: `cpushort` allows only
`normal` and rejects the job. Never export `SBATCH_QOS` (`pipe.sh` unsets
every `SBATCH_*` variable it inherits). POST, BW, CALLP, TSSE, Count and
CountCheck choose their class per job with `runClass MIN_PER_GB FILE...`, which
estimates the run time from the input size using the `RATE_*` values
measured on a full-size run; the other stages are always `SHORT`. A new
stage whose run time grows with its input needs a measured `RATE_*`.

`deliverResults.sh` has not been ported: it still points at the JUNO
paths `/ifs/res/seq/pi/invest` and `~/Code/BIC/Delivery`. Deliver the
`atacSeq/` directory by hand until it is updated (`00.ISSUES.md` #1).

## Pipeline stages (pipe.sh)

1. `postMapBamProcessing_ATACSeq.sh` — MAPQ + flag filter (`-f 3 -F 1804`),
   drop inserts <=30bp, Picard SortSam + CollectInsertSizeMetrics, then
   bamToBed with the Tn5 shift (+4 on `+`, -5 on `-`), keeping only the
   contigs of `lib/genomes/<TAG>.genome` → `*.shifted.bed.gz`; also the
   mitochondrial fraction → `*___MITO.txt`
2. `makeBigWigFromBEDZ.sh` — 10M-read-normalized bigWig via bedtools +
   `bin/wigToBigWig`
3. `callPeaks_ATACSeq.sh` — MACS2 on the BED (not BAM) with ENCODE's
   settings, `--nomodel --shift -75 --extsize 150 --keep-dup all
   --call-summits -p 0.01` (until v2.0.0 `--shift 75` and the default
   `--keep-dup 1`). MACS2's files in `callpeaks/<PREFIX>/` are left as
   written. Then ENCODE's steps, in its order: sort by `-log10 p`, keep
   the top 300,000, clamp negative coordinates, clip to the contig ends
   (`bin/bedClip -truncate`), drop the peaks overlapping
   `lib/genomes/<TAG>.blacklist.bed` and cap the score at 1000. The
   result and its summits are written beside MACS2's files as
   `*_peaks.blacklistFiltered.narrowPeak` and
   `*_summits.blacklistFiltered.bed`; peaks keep MACS2's names (ENCODE
   renames them `Peak_N`). Never overwrite or replace a primary tool's
   output; write derived files under a new, readable name
4. `mergePeaksToSAF.sh GENOME callpeaks` — merge all
   `*_peaks.blacklistFiltered.narrowPeak` (never MACS2's
   `*_peaks.narrowPeak`) within 500bp into `macsPeaksMerged.saf` (the
   peak atlas), then subtract the blacklist again, since a merge can
   bridge a blacklisted region between two kept peaks. SAF is 1-based
   and inclusive: its `Start` is the BED start plus 1 (until v2.0.0 the
   BED start was written as is)
5. `bin/countPeaks.sh` → `bin/featureCounts` (v2.1.1) — raw count matrix
   `peaks_raw_fcCounts.txt` over the atlas, `-O -p --countReadPairs`, with
   `-Q` set to the same MAPQ as stage 1, then Check A. Beside it, one
   `CountCheck` job per sample (`bin/assignedCheck.sh`) for Check B; the
   control job compares both with `bin/checkCounts.sh` (see "Count
   checks")
6. `R/getDESeqScaleFactors.R` — DESeq2 size factors (written **inverted**, per
   R. Koche's convention; see NOTES.md)
7. `bin/computeTSSEnrich.sh` → `R/TSSEnrich/tss_enrich.R` — ENCODE TSS
   enrichment, for every build
8. `plotINSStats.R` and `R/analyzeATAC.R` — insert-size and QC/PCA reports
9. Staging into `atacSeq/{atlas,bigwig,macs,metrics}` for delivery

Differential analysis is a separate manual step (needs replicates), run
from the analysis directory, which gives it the genome
(`00.RUNSTATUS.txt`):

```bash
Rscript R/diffAnalysisPairwise.R sampleManifest.csv Comparisons.csv [RUNTAG]
```

## Genome handling — read this before touching any genome code

Full detail and the reasons are in `docs/GENOMES.md`.

A genome is named by its build tag everywhere: `b37`, `hg38`, `mm10`.
GRCh38 was `b38` until v2.0.0; do not reintroduce that name, species
words (`human`, `mouse`) or stand-ins (`hg19` for b37 data). ChIP-seq
uses the same tags, genome files and `bin/getGenomeBuildBAM.sh` (keep
that file byte-identical in both repos).

Genome is **auto-detected**, never passed by the user: `bin/getGenomeBuildBAM.sh`
md5sums the sorted `@SQ` lines of the BAM header and maps that to a build tag.
`pipe.sh` does this for every BAM before it submits anything, stops
(`STATUS=FAILED`, reason in `MESSAGE`) if the BAMs differ or the tag is
not supported, passes the tag to every stage, and records it as
`GENOME=` in `00.RUNSTATUS.txt`. Nothing detects it again: post-run
tools (`R/diffAnalysisPairwise.R`) read it from the status file and take
no genome argument.

Every genome-dependent value is in `lib/genomes/GENOMES.tsv`; a
supported genome is a row of it. Shell reads it with
`bin/genomeInfo.sh TAG FIELD`, R with `genomeInfo()` in
`R/genomeInfo.R`; both fail on an unknown tag or field. There is no
`case` on genome names anywhere. Adding a genome is data only: its md5
in `getGenomeBuildBAM.sh`, a row in `GENOMES.tsv`,
`lib/genomes/<TAG>.genome`, `lib/genomes/<TAG>.blacklist.bed` (ENCODE's
list, built by `lib/genomes/getBlacklists.sh`), and
`R/TSSEnrich/lib/<TAG>_tss.bed` and `<TAG>.chrom.sizes`. `pipe.sh`
stops if any of the four files is missing.

Chromosome filtering is allowlist-driven: `lib/genomes/<TAG>.genome` lists
autosomes, X and Y only (no mitochondria, no unplaced, random, decoy or
viral contigs). It is the chrom sizes, and the scripts write its BED form
to a temp file for `bedtools intersect -nonamecheck`. Do not reintroduce
`egrep -v` denylists, and never give the 2-column file to
`samtools view -L` (it reads column 2 as one position and selects almost
nothing). Mitochondrial reads are therefore dropped from the shifted BED
for every build, as ENCODE's ATAC pipeline does; POST reports the
fraction in `___MITO.txt`.

## Count checks

The counts are never accepted unchecked, whatever the counter.
featureCounts v1.6.0 (used until v2.0.0), run with several threads on
several BAMs, sometimes did not pair the reads of one file and counted
that sample about twice, exiting 0 with no error. Two checks must pass
before DESEQ runs; `bin/checkCounts.sh SUMMARY [ASSIGNED ...]` does both
and can be run by hand from an analysis directory to audit a run.

- **Check A (pairing)**, in the Count job: each sample's fragment total
  (the sum of every `.summary` row, not rows picked by name) divided by
  mapped reads / 2 (`samtools idxstats`, so Count waits for Index) must
  be within 0.95 to 1.05. Correct runs give 1.000 to 1.005; an unpaired
  BAM, or featureCounts 2.x with `-p` but without `--countReadPairs`,
  gives 2.0.
- **Check B (independent count)**: `bin/assignedCheck.sh` counts the
  distinct read names among reads that `samtools view -L` finds over the
  atlas. It must equal featureCounts' `Assigned` **exactly**; it did for
  all 42 BAMs tested (27 full-size hg38, 11 downsampled hg38, 4
  downsampled b37). A difference fails
  the run and must be explained before anything is changed; do not add a
  tolerance to make it pass. It relies on the postProcess BAM being
  already filtered; if that filter changes, apply it there too.

No Bioconductor package (Rsubread, `GenomicAlignments`, ...) in a
pipeline step that produces counts: their APIs change without warning.
The standalone subread binary is not Bioconductor.

## Cross-stage conventions

- **MD5 handshake**: `postMapBamProcessing` writes `*.bed.gz.md5`; the bigWig
  and peak-calling stages verify it, sleep 300s, retry once, then fail. This
  guards against reading a file still being flushed on the shared filesystem.
  Keep this if you add a stage consuming `.bed.gz`.
- **Sample identity** flows from the BAM `@RG SM:` tag, which becomes the
  `out/<SID>/` directory and the `MapID` column. `sampleManifest.csv`
  (`MapID,SampleID,Group`) is auto-generated by `pipe.sh` if absent, by
  stripping a leading `s_` for `SampleID` and a trailing `[-_]\d+` for `Group`.
  It is usually wrong for real projects and is meant to be hand-edited, after
  which `plotINSStats.R` and `R/analyzeATAC.R` are rerun manually.
- Groups prefixed `EXC` are dropped from differential analysis.
- The R scripts map count-matrix column names back to `SampleID` via the
  manifest; `diffAnalysisPairwise.R` supports two naming conventions (BIC
  `*_postProcess*` and PEmap `*___MD*`). Sample-name mismatches surface here
  as an abort, not silently.
- Comparison sign convention is `X2 - X1`: positive logFC means higher in the
  second column.

## R scripts

- All R scripts depend on helpers from the user's `~/.Rprofile` — `len()`,
  `cc()`, `DATE()`, `halt()`, `suppress()`. They are not defined in this repo
  and are used without qualification; read `~/.Rprofile` before editing.
- `R/diffAnalysisPairwise.R` and `R/TSSEnrich/tss_enrich.R` are the modernized
  scripts (native `|>`, roxygen-style docs, snake_case). `R/analyzeATAC.R`,
  `plotINSStats.R` and `R/getDESeqScaleFactors.R` are older (`%>%`, base R,
  camelCase). Match the file you are in.
- Pipeline R scripts print `## START: <script>` / `## END: <script>` for log
  tracing and wrap library loads in `suppressPackageStartupMessages`.
- Always name the method of `p.adjust()`: its default is Holm, which
  made the diff `FDR` column a Holm correction until v2.0.0. The diff
  FDR is `method = "fdr"`.
- A PCA of samples is `prcomp(t(logCPM))` and its scores `$x` (samples
  as observations), never the `$rotation` of `prcomp(logCPM)`; that
  plotted the profile all samples share as PC1 until v2.0.0.
- `R/attic/` and `attic/` hold superseded versions kept for reference. Do not
  edit them; do not treat them as live code.

## External dependencies not in this repo

`sbatch` (the `SBATCH_SCRIPT_DIR` wrapper), `bsub` (the openlava shim),
`bSync.sh`, `picard` and `bedtools` come from `~/bin`; `samtools` is loaded by `bin/loadTools.sh` via `module load
samtools/1.20`. ImageMagick (`convert`, used only by the unused
`mergePNGs`/`merge_pngs` helpers in the R reports) is `module load
imagemagick` (7.1.2-33); nothing loads it yet. `sacct`, `squeue` and `scancel` are the Slurm client tools.
Scratch for intermediates is
`${ATAC_SCRATCH_ROOT:-/scratch/core001/bic/$USER/ATACSeq}`.

**Nothing may write to `/tmp`.** IRIS nodes set `TMPDIR=/tmp`, a small
volume shared by every job on the node. `bin/loadTools.sh` points
`TMPDIR` at `${ATAC_LOCAL_TMP:-/localscratch/$USER}` (node-local, 2.8T),
and `bin/runStage.sh` gives each stage job its own directory under it,
removed when the job ends. Any new `sort` gets `-T "$TMPDIR"`; any new
tool with its own temp-dir option gets `$TMPDIR` or a scratch dir, never
its default. R is the 4.5.1
on PATH; `tss_enrich.R` needs `optparse`.

`bin/featureCounts`, `bin/wigToBigWig`, `bin/bedGraphToBigWig` and
`bin/bedClip` are vendored Linux x86-64 binaries — they will not run on
macOS, so anything invoking them can only be tested on the cluster.
`bin/bedClip` is the
UCSC v369 build that ENCODE installs
(<https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64.v369/bedClip>,
md5 `7b4d41d51a340f8641e6786383136c63`); the current UCSC build needs
glibc 2.34 and IRIS has 2.28. `bin/featureCounts` is
subread v2.1.1, statically linked, md5
`0eff90d8e30311160dbddda542f17383`, from
`subread-2.1.1-Linux-x86_64.tar.gz` (md5
`ed031bb0734650731f0bfe7ea2053d5a`) on
<https://sourceforge.net/projects/subread/files/subread-2.1.1/>. In 2.x
`-p` only declares paired input; `--countReadPairs` is what makes it
count fragments.

## Docs worth reading

- `NOTES.md` — R. Koche's original method spec; the source of truth for why
  bigWigs and size factors are computed the way they are.
- `QC/QCNotes.md` — which ATAC QC metrics matter and why.
- `docs/SLURM_PORT.md` — how the JUNO/LSF to IRIS/Slurm port was done, the
  flag translation, and the memory and partition choices.
- `docs/GENOMES.md` — genome tags, `GENOMES.tsv`, the genome files, and
  what ENCODE does with mitochondrial reads (with links).
- `docs/UPDATE_TO_B38.md`, `docs/CHECKLIST_B38.md` — the plan for the
  GRCh38 rollout, done in v1.1.0 under the tag `b38` (now `hg38`).
  Historical: the checkboxes were never ticked and the line numbers are
  out of date.
- `CHANGELOG.md` — kept current; add entries for user-visible changes.
- `00.ISSUES.md` — open issues only, numbered in order of work;
  `attic/zzISSUES.md` — the closed ones, with the reason.
