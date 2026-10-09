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
mkdir -p SLURM.CTRL
sbatch /path/to/ATAC-seq/pipe.sh [-q MAPQ] BAM1 [BAM2 ...]
```

`pipe.sh` is the control job: it fans each stage out with `~/bin/bsub` (the
SchedMD openlava shim with local patches), blocks on `bSync JOBNAME`, then
calls `bCheck JOBNAME` to abort if any job in the group did not reach
`COMPLETED` in `sacct`. Job names embed `$$` so concurrent runs don't
collide; every `sacct` query is bounded by `ATAC_START` because pids
recycle. The helpers are in `bin/slurmTools.sh`; `bSync.sh` itself lives in
`~/bin`. Individual stage scripts can also be run standalone for debugging —
each has its own usage block. See `docs/SLURM_PORT.md` for the flag
translation and the reasoning behind the resource requests.

Whether a run worked is answered by `bin/checkRun.sh` (exit 0 worked, 1
failed, 2 running), run from the analysis directory. It reads three records
that `pipe.sh` keeps, plus `sacct`:

- `00.RUNSTATUS.txt`: `STATUS=RUNNING|COMPLETED|FAILED`, the `STAGE`
  reached, and the control job id. `COMPLETED` is written only as the last
  line of `pipe.sh`, after the per-sample deliverables check; the `EXIT`
  trap writes `FAILED`. A new control-job step must call `setStage NAME`
  and a new deliverable belongs in that check.
- `SLURM.CTRL/jobs.tsv`: job id, stage, sample and log of every stage job.
  Submit stage jobs with `atacSub STAGE SAMPLE LOGDIR "BSUB_OPTS" CMD ...`
  and wait with `waitStage STAGE`, never with a bare `bsub`, or the job is
  missing from the manifest and its log has no exit trailer.
- `#ATAC_EXIT=<rc>`: the last line of every stage log, written to stderr
  by `bin/runStage.sh`, which `atacSub` puts in front of every command.

Stage scripts must exit nonzero on any failed step (`pipefail` plus an
explicit check), or `sacct`, `bCheck` and the trailer all report a failed
step as `COMPLETED`.

Every stage runs on the short partitions (`cmobic_short,cpushort`, under
two hours) with `-M` as a hard total-memory cap. If a stage needs longer on
real data, give it `-W` over two hours so the shim picks `cmobic_cpu`, and
export `SBATCH_QOS=priority` for that call only.

`deliverResults.sh` has not been ported: it still points at the JUNO
paths `/ifs/res/seq/pi/invest` and `~/Code/BIC/Delivery`. Deliver the
`atacSeq/` directory by hand until it is updated.

## Pipeline stages (pipe.sh)

1. `postMapBamProcessing_ATACSeq.sh` — MAPQ + flag filter (`-f 3 -F 1804`),
   drop inserts <=30bp, Picard SortSam + CollectInsertSizeMetrics, then
   bamToBed with the Tn5 shift (+4 on `+`, -5 on `-`) → `*.shifted.bed.gz`
2. `makeBigWigFromBEDZ.sh` — 10M-read-normalized bigWig via bedtools +
   `bin/wigToBigWig`
3. `callPeaks_ATACSeq.sh` — MACS2 on the BED (not BAM), `--nomodel --shift 75
   --extsize 150 --call-summits -p 0.01`
4. `mergePeaksToSAF.sh` — merge all narrowPeak within 500bp into
   `macsPeaksMerged.saf` (the peak atlas)
5. `bin/featureCounts` — raw count matrix `peaks_raw_fcCounts.txt` over the atlas
6. `R/getDESeqScaleFactors.R` — DESeq2 size factors (written **inverted**, per
   R. Koche's convention; see NOTES.md)
7. `bin/computeTSSEnrich.sh` → `R/TSSEnrich/tss_enrich.R` — ENCODE TSS enrichment
8. `plotINSStats.R` and `R/analyzeATAC.R` — insert-size and QC/PCA reports
9. Staging into `atacSeq/{atlas,bigwig,macs,metrics}` for delivery

Differential analysis is a separate manual step (needs replicates):

```bash
Rscript R/diffAnalysisPairwise.R GENOME sampleManifest.csv Comparisons.csv [RUNTAG]
```

## Genome handling — read this before touching any genome code

Genome is **auto-detected**, never passed by the user: `bin/getGenomeBuildBAM.sh`
md5sums the sorted `@SQ` lines of the BAM header and maps that to a build tag.
Adding a new reference means adding its md5 there first — everything downstream
keys off that string.

The tag vocabulary is **not uniform across the codebase**, and this is a live
source of bugs:

| Consumer | Accepted tags |
| --- | --- |
| `bin/getGenomeBuildBAM.sh` (emits) | `b37`, `b37_dmp`, `hg19`, `hg19-mainOnly`, `GRCh37-lite`, `b38`, `b37+mm10`, `mm10`, `mm10_hBRAF_V600E`, `mm9Full`, `GRC_m38`, `sCer+sMik_IFO1815` |
| `postMapBamProcessing_ATACSeq.sh` | `b37`, `b38`, `mm10` |
| `callPeaks_ATACSeq.sh`, `makeBigWigFromBEDZ.sh` | `b37`, `b38`, `mm10`, `sCer+sMik_IFO1815` |
| `R/diffAnalysisPairwise.R` | `hg19`, `b38`, `mm10` |

So a detected build that the shell stages accept may still be rejected
downstream, and vice versa. When adding genome support, update **every** case
statement, not just the one that failed.

Chromosome filtering is allowlist-driven: `lib/genomes/<build>.genome` (chrom
sizes) and `lib/genomes/<build>.genome.bed` (regions to keep) are intersected
with `bedtools intersect -nonamecheck`. Do not reintroduce `egrep -v` denylists.

`R/TSSEnrich/lib/` ships only `b38` and `b37` files (`<build>_tss.bed`,
`<build>.chrom.sizes`), so stage 7 fails for any other build until the
matching files are added.

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
- `R/attic/` and `attic/` hold superseded versions kept for reference. Do not
  edit them; do not treat them as live code.

## External dependencies not in this repo

`sbatch` (the `SBATCH_SCRIPT_DIR` wrapper), `bsub` (the openlava shim),
`bSync.sh`, `picard` and `bedtools` come from `~/bin`; `samtools` is loaded by `bin/loadTools.sh` via `module load
samtools/1.20`. `sacct`, `squeue` and `scancel` are the Slurm client tools.
Scratch for intermediates is
`${ATAC_SCRATCH_ROOT:-/scratch/core001/bic/$USER/ATACSeq}`. R is the 4.5.1
on PATH; `tss_enrich.R` needs `optparse`.

`bin/featureCounts`, `bin/wigToBigWig`, and `bin/bedGraphToBigWig` are
vendored Linux x86-64 binaries — they will not run on macOS, so anything
invoking them can only be tested on the cluster.

## Docs worth reading

- `NOTES.md` — R. Koche's original method spec; the source of truth for why
  bigWigs and size factors are computed the way they are.
- `QC/QCNotes.md` — which ATAC QC metrics matter and why.
- `docs/SLURM_PORT.md` — how the JUNO/LSF to IRIS/Slurm port was done, the
  flag translation, and the memory and partition choices.
- `docs/UPDATE_TO_B38.md`, `docs/CHECKLIST_B38.md` — per-file analysis of the
  b38 rollout, including remaining gaps.
- `CHANGELOG.md` — kept current; add entries for user-visible changes.
