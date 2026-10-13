# MUST BE DONE FIRST — v2.0.0 release

Review of every code change of 2026-10-12 (`678c04a`, `5b92c57`,
`b391277`, `424d252`, `1f21b5b`, `6eef929`, `5bd921d`, `507116c`,
`3bbbcb3`, `533308e`, `e75a8b7`, `e50dec6`, `c8d17d7`) for bugs that
affect the atlas or the bigWigs.

## Result

No bug that affects the atlas or the bigWigs was found.

Checked line by line against ENCODE and against each other:
`postMapBamProcessing_ATACSeq.sh` (filter, fixmate, Tn5 shift, contig
filter, MITO count), `makeBigWigFromBEDZ.sh`, `callPeaks_ATACSeq.sh`
(MACS2 flags, p-value sort, clamp, 300,000 cap, bedClip, blacklist,
score cap, summit filter), `mergePeaksToSAF.sh` (one sort, merge,
subtract, SAF start), `bin/countPeaks.sh`, `bin/assignedCheck.sh`,
`bin/checkCounts.sh`, `bin/runStage.sh`, `pipe.sh` (submissions,
memory, staging, deliverables check), `lib/genomes/getBlacklists.sh`
and the three blacklist files (no bad intervals, every contig in its
`.genome` file).

## Must be done before the release

1. **`00.FIXME.md` #1: the full-size run.** Every v2.0.0 change was
   tested only on 100x downsampled data. From the run now in progress,
   record in `docs/SLURM_PORT.md`:
   - POST Elapsed from `sacct` against `RATE_POST=6` (`pipe.sh:438`);
     raise the rate if any POST job is near 1:59:00.
   - The `runClass` lines of every stage against `sacct` Elapsed.
   - `peaks of <PREFIX>: MACS2 N, top N (cap 300000), kept N` from
     each `SLURM.03.CALLP` log: how many samples hit the cap.
   - Atlas size (`wc -l macsPeaksMerged.saf`).
   - `CHECK_A ... ratio=1.0000` and `CHECK_B ... diff=0` for every
     sample in the Count and control logs. Check B must be exactly
     equal; do not add a tolerance.
   - Peak memory of POST (`sacct MaxRSS`): the job runs `samtools
     sort -m 4G` and then picard at `-Xmx24g` under `-M 32G`
     (`pipe.sh:563`).

2. **`docs/GENOMES.md:96`** says the MITO fraction is counted on the
   BAM filtered by "MAPQ, flags, insert size". There is no insert-size
   filter since `c8d17d7`. Change to "MAPQ, flags, mate present".

## Not blocking (noted, no action required for v2.0.0)

- `callPeaks_ATACSeq.sh:169`: at the 300,000 cap, peaks with equal
  -log10 p are ordered by `sort`'s tie-breaking, which may differ from
  ENCODE's. Only matters for a sample with more than 300,000 peaks,
  and only for which tied peaks fall at the boundary.
- `postMapBamProcessing_ATACSeq.sh:117`: `samtools fixmate` has no
  `-u`; it compresses its output before the next `samtools view`
  reads it. Performance only.
