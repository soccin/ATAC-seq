# ENCODE ATAC-seq pipeline concordance review

Review date: 2026-10-12. Checkout reviewed: `v2.x-pre2-11-g58a3322`
(branch `devs`, clean tree). Reference: ENCODE-DCC/atac-seq-pipeline
`master` (v2.2.x), <https://github.com/ENCODE-DCC/atac-seq-pipeline>,
as read on 2026-10-12, plus the ENCODE ATAC-seq data standards page.
The exact files read are listed under "Sources".

Updated 2026-10-12 on branch `fix/peak-cap` (from `encode-2026`):
C-01 to C-04 are fixed (section 1, "Status"); the other findings
are open and their text is brought up to date with the fixed code.

The review covers every stage from the input BAM to the delivered
`atacSeq/` directory: `pipe.sh`, `postMapBamProcessing_ATACSeq.sh`,
`makeBigWigFromBEDZ.sh`, `callPeaks_ATACSeq.sh`, `mergePeaksToSAF.sh`,
`bin/countPeaks.sh`, `bin/assignedCheck.sh`, `bin/checkCounts.sh`,
`bin/computeTSSEnrich.sh`, `R/TSSEnrich/tss_enrich.R`,
`R/getDESeqScaleFactors.R`, `plotINSStats.R`, `R/analyzeATAC.R`,
`R/diffAnalysisPairwise.R`, `lib/genomes/*`, `R/TSSEnrich/lib/*`. The
`attic/` directories were not reviewed.

Conventions used in this document:

- Paths such as `callPeaks_ATACSeq.sh` or `bin/countPeaks.sh` are
  relative to the root of this repository, the ATAC-seq toolkit
  directory that holds `pipe.sh`. Line numbers (`file:NN`) refer to
  branch `fix/peak-cap` and will drift with later commits; the
  quoted code identifies the line when the number has moved.
- Paths such as `out/<SID>/`, `callpeaks/`, `atacSeq/` and
  `00.RUNSTATUS.txt` are in the per-project analysis directory from
  which `pipe.sh` is submitted, not in the repository.
- ENCODE paths (`atac.wdl`, `src/encode_task_*.py`, `scripts/*`,
  `docs/input.md`) are relative to the root of the ENCODE repository
  linked above.
- PEMapper is the BIC mapping pipeline that produces the input BAMs:
  BWA-MEM alignment, then Picard MarkDuplicates with duplicates
  marked, not removed. POST, BW, CALLP, Count, CountCheck, DESEQ and
  TSSE are the stage names `pipe.sh` records in `00.RUNSTATUS.txt`.
- Proj_18143 and the `testNN` runs are validation runs recorded in
  `docs/SLURM_PORT.md`: Proj_18143 is a full-size run of 15 hg38
  samples (made before C-01 to C-03 were fixed). The hg38 test runs
  (test26, test27, test33, test35, test37, test39) use 11 hg38 BAMs
  downsampled 100x at MAPQ 20; the b37 runs (test25, test28, test34,
  test36, test38, test40) use 4 b37 BAMs downsampled 100x at MAPQ 10.
- Findings C-01 to C-22 are numbered in stage order, not by severity.
  High means the finding changes peaks or counts materially; Medium
  means it changes a delivered number or hides a QC failure; Low means
  a small or documented difference. "Design" marks a deliberate choice
  of this pipeline; "inherent" marks a difference fixed by the
  upstream aligner that this pipeline cannot remove.

## 1. Summary

The pipeline follows the ENCODE ATAC-seq pipeline in its overall shape:
proper-pair and flag filtering with `-F 1804`, mitochondrial reads
dropped before peak calling, `bedtools bamtobed` plus the +4/-5 Tn5
shift, and MACS2 on the shifted tagAlign with `--nomodel`, `-p 0.01`,
`--extsize 150`, `--call-summits`. The TSS enrichment module is a
faithful port of ENCODE's and uses ENCODE's own TSS files.

At the reviewed checkout it was not concordant with ENCODE in the
peak-calling step itself (C-01, C-02) and filtered no peaks against a
blacklist (C-03). Those three are fixed, and so is C-04, the peak cap
(see "Status" below); each sample's final peak set is now the one
ENCODE's commands produce. ENCODE's FRiP and QC metrics are still
missing (C-11, C-16), and the atlas, counts, size factors, QC PDF and
differential analysis have no ENCODE counterpart (section 5). The
table lists every finding:

| ID | Finding | Severity |
| --- | --- | --- |
| C-01 | MACS2 `--shift` has the wrong sign (`+75`; ENCODE `-75`) | High; fixed |
| C-02 | MACS2 `--keep-dup all` missing (default `1` collapses cut sites) | High; fixed |
| C-03 | No blacklist filtering of peaks, anywhere | High; fixed |
| C-04 | No peak capping, p-value sorting or coordinate clipping after MACS2 | Medium; fixed |
| C-05 | Read filter: MAPQ 10 on BWA output; ENCODE keeps up to 4 multimappers with no MAPQ filter (or `-q 30` when multimapping is off) | Medium |
| C-09 | Duplicates: Picard MarkDuplicates upstream, removed here by flag, as ENCODE does; but marked on the unfiltered BAM, with no check that the input was marked and no duplicate rate or NRF/PBC1/PBC2 reported | Low-Medium |
| C-10 | `frac_mito` is computed after filtering; ENCODE computes it on the raw BAM | Medium |
| C-11 | "Fraction reads in peaks" is not ENCODE FRiP (atlas instead of own peaks; mito reads in the denominator) | Medium |
| C-16 | Most ENCODE QC metrics are absent (library complexity, fragment-length QC pass/fail, annotation enrichment, number of peaks, JSD, GC bias, reproducibility) | Medium |
| C-06 | Orphan reads (mate failed MAPQ) kept; ENCODE removes them with `fixmate -r` | Low |
| C-07 | Extra filter: fragments with insert size <= 30 bp dropped | Low |
| C-08 | Supplementary alignments (flag 0x800) pass the filter; a BWA-specific consequence of `-F 1804` | Low |
| C-12 | Signal track is a 10M-read-scaled read coverage, not MACS2 fold-enrichment and p-value tracks | Low (design) |
| C-13 | No pooled, pseudoreplicate, overlap or IDR peak sets; the merged atlas replaces them | Low (design) |
| C-14 | TSS enrichment: older ENCODE TSS annotation; 1 bp minus-strand offset for odd read lengths; 0.85% residual difference to the Python original | Low |
| C-15 | Fragment-length QC: histogram only; no NFR fraction, NFR/mono-nucleosome ratio or peak-presence checks | Low |
| C-17 | SAF atlas coordinates are off by one (BED start used as 1-based start) | Low |
| C-18 | Atlas merges peaks within 500 bp; ENCODE never merges | Low (design) |
| C-19 | Tn5 shift lacks ENCODE's guard against start >= end | Very low |
| C-20 | Tool versions differ (MACS2 2.2.9.1 vs 2.2.4, samtools 1.20 vs 1.9, bedtools unknown vs 2.29.0) | Low |
| C-21 | Aligner, trimming and reference differ (BWA-MEM on GDC GRCh38.d1.vd1 vs cutadapt + bowtie2 `-X2000` on the no-alt analysis set) | Inherent |
| C-22 | Documentation out of step with the code: `NOTES.md` describes MACS2 on the BAM, `README.md` says duplicates are marked "by the mapper", `R/TSSEnrich/README.md` gives one threshold rule for every genome | Low (docs) |

### Status

| ID | Branch, commits | Tested on |
| --- | --- | --- |
| C-01, C-02 | `fix/macs2-encode` (`424d252`, `60bdde4`), merged into `encode-2026` as `be43089` | test27 (hg38), test28 (b37) |
| C-03 | `fix/blacklist` (`1f21b5b` to `3bbbcb3`), merged into `encode-2026` as `330ced3` (first version) and `cdf362b` | test33, test35, test37 (hg38); test34, test36, test38 (b37) |
| C-04 | `fix/peak-cap` | test39 (hg38), test40 (b37) |

test33 to test36 ran before the last two code commits of
`fix/blacklist`; test37 and test38 ran on `encode-2026` at `bb486a4`,
which has them both and the one-`sort` fix of `mergePeaksToSAF.sh`
(`fix/atlas-sort`, `533308e`): the five-file deliverables check in
`pipe.sh` (`507116c`) passed, and the peak files, atlas, counts and
scale factors match test27/test28 and test31/test32 as test33/test34
did. `lib/genomes/getBlacklists.sh`
(`3bbbcb3`) is not part of a run; it was tested on its own (a normal
run gives the committed files; a wrong md5, a wrong file served and a
failed download each leave them unchanged).

All ten test runs passed every job and both count checks. Each fix is
in `CHANGELOG.md` under Breaking changes. Not yet done: a full-size run
of Proj_18143 on `encode-2026`. Details are under each finding in
section 3.

## 2. Stage map

| ENCODE task | ENCODE action | This pipeline | Concordant |
| --- | --- | --- | --- |
| `trim_adapter` | cutadapt `-e 0.1 -m 5`, adapters given or auto-detected | Not applicable: input is an aligned BAM | n/a (C-21) |
| `align` | bowtie2 `-k 5 -X2000 --mm` (default `multimapping=4`, passed as `-k multimapping+1`), then `samtools flagstat` on the raw BAM and on the raw BAM without `chrM`; read length from the FASTQ | Upstream BWA-MEM (PEMapper) | No (C-21) |
| `frac_mito` | A second alignment of the trimmed reads to a mitochondria-only bowtie2 index (`align_mito`); `frac_mito` = reads mapped there / (those + reads mapped to non-mito contigs in the main alignment), all unfiltered | `___MITO.txt`: `chrM` reads / all reads of the filtered BAM | Partly (C-10) |
| `filter` step 1 | `samtools view -F 524 -f 2`, name sort, `assign_multimappers.py -k 4`, `fixmate -r`, `samtools view -F 1804 -f 2`; or with `multimapping=0`: `-F 1804 -f 2 -q 30`, `fixmate -r`, `-F 1804 -f 2` | `samtools view -q 10 -f 3 -F 1804`, then drop insert size <= 30 | Partly (C-05, C-06, C-07, C-08) |
| `filter` step 2 | Picard MarkDuplicates on the filtered BAM, then `samtools view -F 1804 -f 2` | Picard MarkDuplicates upstream (PEMapper, marked not removed, on the unfiltered BAM); removed by the same `-F 1804` | Yes, with a caveat (C-09) |
| `filter` step 3 | Remove `chrM`, `MT` from the nodup BAM (`samtools view -L` with the other contigs) | Keep only autosomes, X, Y in the shifted BED; the BAM keeps everything | Yes at the peak level (C-10) |
| `filter` QC | `samstat`, `dup.qc`, `lib_complexity.qc` (NRF, PBC1, PBC2) | None | No (C-09, C-16) |
| `bam2ta` | `bedtools bamtobed -bedpe -mate1` to two tagAlign lines per pair; `+4/-5` shift | `bedtools bamtobed` per read; identical shift awk | Yes (C-19) |
| `call_peak` | `macs2 callpeak -f BED --nomodel --shift -75 --extsize 150 -p 0.01 --keep-dup all --call-summits -B --SPMR` | Same without `-B --SPMR` (signal bedGraphs, not used); until C-01 and C-02 were fixed, `--shift +75` and no `--keep-dup all` | Yes (C-12) |
| `post_call_peak` | sort by `-log10 p`, top 300,000, clip to chromosome, rename `Peak_N`, blacklist filter, `chr[\dXY]+` filter, FRiP, num_peak, peak_region_size, bigBed | Sort, top 300,000, clamp, clip (`bedClip`), blacklist filter and score cap as ENCODE's, written beside MACS2's files as `*.blacklistFiltered.*`, peaks not renamed; `chr[\dXY]+` by the contig allowlist; no FRiP or peak count; the filtered peaks of all samples merged within 500 bp into one SAF, blacklist subtracted | Partly (C-11, C-18) |
| `macs2_signal_track` | `-B --SPMR`, `bdgcmp -m FE` and `-m ppois`, bigWig | `bedtools genomecov -bg -scale 1e7/reads` of the shifted read spans, `wigToBigWig` | No (C-12) |
| `spr`, `pool_ta`, `overlap`, `idr`, `reproducibility` | Pseudoreplicates, pooled peaks, naive overlap, IDR 0.05, rescue and self-consistency ratios | None; IDR is installed in the venv but never run | No (C-13) |
| `tss_enrich` | metaseq on the nodup BAM, +/-2 kb, 400 bins, shift `-read_len/2`, Greenleaf normalization | R port with the same algorithm and ENCODE's TSS files | Yes (C-14) |
| `fraglen_stat_pe` | Picard `CollectInsertSizeMetrics W=1000 STOP_AFTER=5000000` on the nodup BAM; NFR >= 0.4, NFR/mono >= 2.5, peak checks | Picard `CollectInsertSizeMetrics` on the filtered BAM; plots only | Partly (C-15) |
| `jsd`, `gc_bias`, `annot_enrich`, `xcor`, `preseq` | deepTools fingerprint, Picard GC bias, fraction of reads in DHS/promoter/enhancer, cross-correlation (off by default), preseq (off) | None | No (C-16) |
| none | | featureCounts over the atlas, DESeq2 size factors, insert-size PDF, QC PDF with PCA, edgeR differential analysis | Outside ENCODE's scope (section 5) |

## 3. Findings

Each finding gives the ENCODE behavior, this pipeline's behavior, the
effect on the deliverables, and the fix.

### C-01 MACS2 `--shift` has the wrong sign (High; fixed)

ENCODE (`encode_task_macs2_atac.py`):

```python
shiftsize = -int(round(float(smooth_win)/2.0))   # -75
'macs2 callpeak -t {ta} -f BED -n {prefix} -g {gensz} -p {pval_thresh} '
'--shift {shiftsize} --extsize {extsize} '
'--nomodel -B --SPMR --keep-dup all --call-summits'
```

With `--nomodel`, MACS2 moves each read's 5' end by `--shift` in the
read's 5' to 3' direction and extends `--extsize` bases from there. With
`-75 / 150` every cut site becomes a 150 bp window centered on itself.
The pipeline computed `shiftsize=$((smooth_window / 2))`, `+75`, which
puts the window 75 to 225 bp downstream of the cut site, on the read's
strand: reads on the two strands were pushed away from the accessible
region in opposite directions, so pileups were split, off center and
lower. The Kundaje `callpeak_macs2_atac.bds` module it was copied from
passes `-$shiftsize`; the minus sign was lost in the port.

Effect on earlier runs: every narrowPeak, summit, the atlas, the count
matrix, the size factors, the QC PDF and every differential result.

Fixed in `424d252`: `callPeaks_ATACSeq.sh:96`,
`shiftsize=$(( -smooth_window / 2 ))`. `-B --SPMR` are left out on
purpose: they only write signal bedGraphs, which this pipeline does not
use (C-12).

### C-02 `--keep-dup all` is missing (High; fixed)

ENCODE passes `--keep-dup all`. The pipeline passed nothing, so MACS2
used its default `--keep-dup 1`: at most one read per 5' position and
strand. The input already has PCR duplicates removed (flag 1024 is in
`-F 1804`), so what MACS2 removed were distinct fragments whose cut
sites coincide, most at the strongest peaks; on the full-size
Proj_18143 run that was 5 to 17% of each sample's reads.

Fixed in `424d252`: `callPeaks_ATACSeq.sh:114`, `--keep-dup all`.

Test of C-01 and C-02 together (`docs/SLURM_PORT.md`, test27 against
test26 for hg38, test28 against test25 for b37; earlier values in
parentheses):

| | hg38 (11 BAMs) | b37 (4 BAMs) |
| --- | --- | --- |
| Peaks per sample | 7% fewer to 0.4% more | 8-17% fewer |
| Median peak width | 203-215 bp (208-245) | 199-209 bp (150-229) |
| Fraction of reads in own peaks | 0.90-0.95 (0.70-0.78) | 0.86-0.89 (0.59-0.69) |
| Top 10,000 summits within 500 bp of a TSS | 21-29% (14-28%) | 25-32% (27-31%) |

### C-03 No blacklist filtering (High; fixed)

ENCODE (`encode_lib_blacklist_filter.py`, called by
`encode_task_post_call_peak_atac.py` for every peak set):

```bash
bedtools intersect -nonamecheck -v -a peaks -b blacklist \
  | awk 'BEGIN{OFS="\t"} {if ($5>1000) $5=1000; print $0}' \
  | grep -P 'chr[\dXY]+\b'
```

The pipeline had no blacklist file and no filter: not on the
per-sample narrowPeaks, the atlas or the count matrix. Blacklisted
regions (satellite repeats, centromeric and telomeric regions,
high-signal artifacts) passed `-p 0.01` in every sample and were among
the highest-count rows: in the full-size Proj_18143 run 1,376 of
654,370 atlas rows overlapped the blacklist, with 0.43% of all counts,
3 of the 10 and 7 of the 100 rows with the most counts.

Fixed on `fix/blacklist`:

- **Blacklist files**, `lib/genomes/<TAG>.blacklist.bed`, the files
  ENCODE's `download_genome_data.sh` (`VER=v3`) uses. Each genome has
  one under this fixed name (no `GENOMES.tsv` column); `pipe.sh`
  stops before submitting anything if the run's file is missing.

  | Tag | ENCODE file | Regions | Notes |
  | --- | --- | --- | --- |
  | `hg38` | `ENCFF356LFX` | 910 | Boyle lab blacklist v2, as downloaded |
  | `mm10` | `ENCFF547MET` | 164 | Boyle lab blacklist v2, as downloaded |
  | `b37` | `ENCFF001TDO` | 410 | the hg19 list (`wgEncodeDacMapabilityConsensusExcludable`); `chr` stripped, `chrM` (1 region) dropped |

  `lib/genomes/getBlacklists.sh` rebuilds them: it checks each
  download's md5 against the ENCODE portal and that every contig is in
  `<TAG>.genome`, and replaces the committed files only if every check
  passes for all three.
- **Per-sample filter**, `callPeaks_ATACSeq.sh:203-210`: ENCODE's
  `intersect -v` and score cap. The `chr[\dXY]+` step is already done
  by the contig allowlist (`lib/genomes/<TAG>.genome`). Unlike the fix
  first proposed here, MACS2's own files are not changed: the filtered
  peaks are written beside them in `callpeaks/<PREFIX>/` as
  `<PREFIX>_peaks.blacklistFiltered.narrowPeak`, with the summits of
  the kept peaks in `<PREFIX>_summits.blacklistFiltered.bed`. Both sets
  are delivered in `atacSeq/macs/<PREFIX>/`, and the `pipe.sh`
  deliverables check (`pipe.sh:721-734`) requires all five files.
- **Atlas**, `mergePeaksToSAF.sh GENOME MACS_OUTDIR`: merges only the
  `*_peaks.blacklistFiltered.narrowPeak` files, then subtracts the
  blacklist (`mergePeaksToSAF.sh:44-48`), because merging within
  500 bp can bridge a small blacklisted region between two kept peaks
  (5 atlas peaks in the hg38 test, 1 in b37); the subtraction splits
  such a peak in two. No atlas peak, and so no count, covers a
  blacklisted base.

Tests (`docs/SLURM_PORT.md`, "Blacklist filter"): per sample 264-475
hg38 peaks and 652-891 b37 peaks are removed; 0 filtered peaks and 0
atlas peaks overlap the blacklist; MACS2's narrowPeak and summits are
byte-identical to the unfiltered runs (test27, test28), and so is the
body of its `.xls` (its `#` header differs only in the per-job temp
path of the input file); the hg38 atlas has 901,681 peaks
(903,877 before) and the b37 atlas 758,830 (760,751).

### C-04 No capping, sorting or clipping after MACS2 (Medium; fixed)

ENCODE (`encode_task_macs2_atac.py`, then `bed_clip` in
`encode_lib_genomic.py`), on MACS2's narrowPeak:

```bash
LC_COLLATE=C sort -k 8gr,8gr peaks.narrowPeak \
  | awk 'BEGIN{OFS="\t"}{$4="Peak_"NR; if ($2<0) $2=0; if ($3<0) $3=0;
         if ($10==-1) $10=$2+int(($3-$2+1)/2.0); print $0}' > tmp
head -n 300000 tmp > tmp2
bedClip tmp2 chrom.sizes out -truncate -verbose=2
```

`encode_task_post_call_peak_atac.py` then blacklist-filters that
result (C-03). The cap is `cap_num_peak`, 300,000 in `atac.wdl`
(`docs/input.md` says 500,000; the WDL value is the one that runs).

The pipeline had no cap, sort or clip: every blacklist-filtered peak
went into the atlas. With `-p 0.01` on full-size samples that is
several hundred thousand peaks per sample; the low-significance tail is
mostly noise and adds atlas rows that `filterByExpr` removes later. On
the full-size Proj_18143 run (old MACS2 settings) 9 of 15 samples had
more than 300,000 peaks, up to 396,453.

Fixed on `fix/peak-cap`, `callPeaks_ATACSeq.sh:169-201`: ENCODE's
steps in ENCODE's order, before the blacklist filter, which now reads
their result. `bin/bedClip` is the UCSC v369 build ENCODE installs
(`scripts/install_conda_env.sh`); the current UCSC build needs glibc
2.34 and IRIS has 2.28. The output keeps its names,
`*_peaks.blacklistFiltered.narrowPeak` and
`*_summits.blacklistFiltered.bed`, now at most 300,000 peaks sorted
by p-value; MACS2's files are unchanged. Two differences from ENCODE,
neither of which changes a peak, a coordinate or a count:

- Peaks keep MACS2's names; ENCODE renames them `Peak_N`. MACS2's
  names tie each peak to its summit and its `.xls` row.
- For a missing summit ENCODE writes the peak center as a contig
  position, but column 10 is an offset from the start; it is written
  as an offset here. MACS2 with `--call-summits` does not write `-1`,
  so neither form was used in any test.

A peak `bedClip -truncate` would turn into start > end (one wholly
past its contig end, which MACS2 cannot call) stops the job, as does
any `bedClip` error.

Tests (`docs/SLURM_PORT.md`, "Peak cap"): for every sample of test39
(hg38) and test40 (b37), the `blacklistFiltered` narrowPeak equals the
output of ENCODE's commands above plus its blacklist filter, every
column but the name. No hg38 sample reached the cap, and the hg38
atlas, counts and scale factors are byte-identical to test37's. Three
of four b37 samples were capped at 300,000; the b37 atlas went from
758,830 to 727,125 peaks. `bedClip` found no peak to clip in any
sample.

### C-05 Read filter: MAPQ threshold and multimappers (Medium)

ENCODE `rm_unmapped_lowq_reads_pe` has two paths. The default
(`multimapping = 4`) is:

```bash
samtools view -F 524 -f 2 -u raw.bam | samtools sort -n
samtools view -h | assign_multimappers.py -k 4 --paired-end | samtools fixmate -r
samtools view -F 1804 -f 2 -u | samtools sort
```

`assign_multimappers.py` drops a pair only if bowtie2 reported more
than 4 alignments for it (bowtie2 runs with `-k 5`); otherwise it
writes all its records and the final `-F 1804` keeps the primary one.
No MAPQ filter is applied in this path; `mapq_thresh = 30` is used
only when `multimapping = 0`:

```bash
samtools view -F 1804 -f 2 -q 30 -u raw.bam | samtools sort -n
samtools fixmate -r
samtools view -F 1804 -f 2 -u | samtools sort
```

This pipeline (`postMapBamProcessing_ATACSeq.sh:110`) applies
`-q 10 -f 3 -F 1804` once, per read, on BWA-MEM output. The default
MAPQ is 10 (`pipe.sh`); the recent test runs used 10 and 20.

The two cannot be made identical: bowtie2 and BWA score MAPQ
differently, and BWA reports no alternative alignments to count.
ENCODE's default path keeps reads that BWA would give MAPQ 0 to (equal
best hits), so this pipeline is stricter at repeats; its MAPQ 10 is
looser than ENCODE's explicit `-q 30` mode for reads with a close
suboptimal hit. Either way the retained read sets differ at repeats,
which is where blacklist regions sit (C-03).

Fix: none that makes it equal. Document the choice in `README.md`, and
consider `MAPQ=30` as the default, which is the threshold ENCODE uses
whenever it filters on MAPQ at all (ATAC with `multimapping = 0`, and
the ENCODE ChIP-seq pipeline, which does use BWA).

### C-06 Orphan reads are kept (Low)

ENCODE runs `samtools fixmate -r` after the MAPQ filter and then
`-f 2` again, so a read whose mate was removed is removed too. This
pipeline filters each read on its own, so a read whose mate fell below
MAPQ stays, with the proper-pair flag still set. `bin/checkCounts.sh`
documents the effect: Check A ratios of 1.000 to 1.004 are expected,
the excess being singletons whose mate failed the MAPQ filter. The
full-size Proj_18143 run gave 1.0004 to 1.0008 (`docs/SLURM_PORT.md`),
that is under 0.1% of fragments. featureCounts counts them as
fragments; the tagAlign has one line for them instead of two.

Fix (optional): `samtools sort -n | samtools fixmate -r | samtools
view -f 2` between the MAPQ filter and SortSam. Check A's expected
ratio then becomes exactly 1.000.

### C-07 Insert size <= 30 bp filter (Low)

`postMapBamProcessing_ATACSeq.sh:111` drops pairs with `|TLEN| <= 30`.
ENCODE has no such filter. Fragments that short are adapter dimers or
artifacts and are rare among BWA proper pairs, so the effect is small,
but it is an undocumented deviation and it slightly truncates the left
edge of the insert-size histogram that ENCODE reports in full.

Fix: document it in `README.md`, or remove it for concordance.

### C-08 Supplementary alignments pass the filter (Low)

`-F 1804` excludes unmapped, mate unmapped, secondary, QC-fail and
duplicate reads but not supplementary alignments (0x800). ENCODE uses
the same mask, but bowtie2 never emits supplementary records, so for
ENCODE the question does not arise. BWA-MEM does emit them for
chimeric reads, and they carry the proper-pair flag of the primary.
Each one adds an extra tagAlign line and an extra read span to the
bigWig coverage, and is counted as a read in the 10M scale factor
(`bamtobed` appends `/1` and `/2` to names, so `sort | uniq` counts
reads, and a supplementary record shares its primary's name and
suffix).

Fix: measure first (`samtools view -c -f 2048` on a postProcess BAM).
If the fraction is not negligible, use `-F 3852` (1804 + 2048).
`assignedCheck.sh` reads the same BAM, so Check B stays consistent.

### C-09 Duplicates: same tool, different stage, no check and no metrics (Low-Medium)

The input BAMs come from PEMapper with duplicates marked, not removed,
by Picard MarkDuplicates (`REMOVE_DUPLICATES=FALSE`, flag 0x400 set).
This is what ENCODE does too: `mark_dup_picard` runs MarkDuplicates
with `REMOVE_DUPLICATES=FALSE` and `rm_dup_pe` then drops the flagged
reads with `samtools view -F 1804 -f 2`, which is this pipeline's
filter at `postMapBamProcessing_ATACSeq.sh:110`. The tool, the
pairing-aware duplicate definition and the removal are concordant.

The one difference in the marking: ENCODE marks duplicates on the BAM
after the MAPQ, proper-pair and orphan filters; PEMapper marks them
on the unfiltered BAM. MarkDuplicates keeps the highest-scoring pair of
each duplicate set, so when the kept pair later fails the MAPQ filter
and a flagged copy would have passed, ENCODE keeps one fragment and
this pipeline keeps none. The reverse cannot happen. The effect is a
small loss of fragments at low-MAPQ positions, which are mostly
repeats; it does not inflate anything.

What is missing is the record. ENCODE reports the Picard duplicate
metrics (`dup.qc`) and the library complexity metrics from
`pbc_qc_pe`:

```bash
bedtools bamtobed -bedpe -i nmsrt.bam \
  | awk '{print $1,$2,$4,$6,$9,$10}' | grep -v "^chrM\s" | sort | uniq -c \
  | awk '... mt, m0, m1, m2, m0/mt (NRF), m1/m0 (PBC1), m1/m2 (PBC2)'
```

The ENCODE standards are NRF > 0.9, PBC1 > 0.9, PBC2 > 3.

This pipeline trusts that the input was marked and keeps no record of
the duplication. The requirement is not stated anywhere: `README.md:49-50`
says the marking is done "upstream, by the mapper", but BWA-MEM marks
nothing; it is a separate Picard MarkDuplicates step in PEMapper, and a
BAM from any other source has to have had that step too. Nothing
verifies the marking: a BAM that was never run through MarkDuplicates
passes `-F 1804` unchanged, `tss_enrich.R`'s guard sees zero flagged reads and passes,
and every downstream result is computed on duplicated data with no
warning. The Picard metrics file PEMapper writes stays in the mapping
project; the duplicate rate, NRF, PBC1 and PBC2 are not in
`atacSeq/metrics`, so the ENCODE depth standard (25 million
non-duplicate, non-mitochondrial fragments) and the complexity
standards cannot be checked from the delivery.

Fix: state in `README.md` that the input BAMs must be Picard
MarkDuplicates-marked (not removed), and correct "by the mapper". In
POST, count `-f 1024` reads in the input BAM (`samtools flagstat`) and
fail the stage when it is zero unless an explicit option accepts an
unmarked BAM; write the duplicate rate (flagged over
mapped) to a metrics file. Compute NRF, PBC1 and PBC2 with the ENCODE
awk above on the input BAM filtered with `-F 780 -f 2 -q MAPQ` (1804
minus the duplicate bit), so the figures match ENCODE's definition
from a duplicate-marked, not deduplicated, BAM.

### C-10 Mitochondrial reads: same effect on peaks, different metric (Medium)

Concordant part: since v2.0.0 mitochondrial reads are out of the
shifted BED for every build, so out of peaks, bigWigs and the atlas,
as in ENCODE (`filter_chrs = ['chrM', 'MT']`).

Differences:

- ENCODE's `frac_mito` is computed before any filtering. The `align`
  task runs `samtools flagstat` on the raw BAM and on the raw BAM with
  `chrM` removed (`encode_task_post_align.py`); a second `align` call,
  `align_mito`, maps the same trimmed reads to a mitochondria-only
  bowtie2 index. `encode_task_frac_mito.py` then reports
  `frac_mito = Rm / (Rn + Rm)`, with `Rm` the reads mapped to the
  mitochondrial index and `Rn` the reads mapped to non-mitochondrial
  contigs in the main alignment, before any MAPQ, proper-pair or
  duplicate filter. This pipeline's `___MITO.txt` counts reads of the
  filtered BAM (`postMapBamProcessing_ATACSeq.sh:140`). Mitochondrial
  reads are heavily duplicated, so the pipeline's fraction is lower
  than ENCODE's for the same library and the two are not comparable.
- ENCODE removes `chrM` from the nodup BAM itself, so the BAM-based
  QC (insert sizes, TSS, samstats) excludes it. Here the postProcess
  BAM keeps `chrM` and every non-allowlisted contig. `CollectInsertSizeMetrics`
  therefore includes mitochondrial fragments, whose size distribution
  differs from nuclear chromatin, and the featureCounts totals that
  feed `R/analyzeATAC.R` include them (C-11). The `n_mapped_reads`
  column of the TSS CSV (`tss_enrich.R:169-178`, from `samtools
  flagstat` of the whole BAM) includes them too; the score does not.

Fix: compute `___MITO.txt` from the input BAM's mapped reads
(`samtools idxstats` on the input BAM in POST, indexing it first if it
has no `.bai`; `bin/checkCounts.sh` already uses `idxstats` this way on
the postProcess BAM). A mitochondria-only realignment is not possible
without FASTQs, so the raw-BAM fraction is the closest equivalent.
Either write a second, contig-filtered BAM for the metrics or document
that the metrics include `chrM`. Keep the current post-filter number
too if it is wanted; name the two differently.

### C-11 Fraction of reads in peaks is not ENCODE FRiP (Medium)

ENCODE (`encode_lib_frip.py`):

```bash
bedtools intersect -nonamecheck -a ta -b bfilt_peaks -wa -u | wc -l   # over
wc -l ta                                                              # all
```

Reads (tagAlign lines, chrM-free, both mates) of a sample over that
sample's own blacklist-filtered peaks, divided by all its reads. The
standard is FRiP > 0.3, > 0.2 acceptable. ENCODE also reports it for
the pooled, overlap and IDR sets.

`R/analyzeATAC.R:68-81` plots `Assigned / sum(all status rows)` from
the featureCounts summary: fragments over the merged all-sample atlas
(blacklist-filtered since C-03) divided by all fragments of the
postProcess BAM, including `chrM` and non-allowlisted contigs. The
numerator is inflated by other samples' peaks and the 500 bp merging;
the denominator is inflated by mitochondrial fragments. The figure is
useful for comparing samples within a run but is not ENCODE FRiP and
cannot be judged against 0.3.

Fix: add a FRiP step per sample with the ENCODE command on the shifted
BED and the sample's own `*_peaks.blacklistFiltered.narrowPeak`
(capped and filtered, as ENCODE's `bfilt` peaks are), written to
`out/<SID>/` and staged to `atacSeq/metrics`. Relabel the existing
plot as "fragments in atlas".

### C-12 Signal tracks differ by design (Low, design)

ENCODE produces two bigWigs per sample from the same MACS2 run:
fold-enrichment over the local lambda and `-log10 p`, both from
`-B --SPMR` pileups of the shifted-and-extended 150 bp windows,
normalized per million reads (`bdgcmp -m FE`, `-m ppois -S reads/1e6`).

This pipeline (`makeBigWigFromBEDZ.sh:43-71`) builds a coverage track
of the shifted read spans with `bedtools genomecov -bg`, scaled to 10
million reads. This is R. Koche's method (`NOTES.md`) and is a
deliberate choice; it is not an ENCODE output and should not be
compared to ENCODE's signal values. Two details worth knowing:

- The Tn5 shift moves only the 5' end, so each `+` read span is 4 bp
  shorter and each `-` read span 5 bp shorter than the alignment. For a
  coverage track this is harmless but it is not a "cut-site" track
  either. ENCODE's optional `count_signal_track` (`genomecov -5 -bg`,
  per strand) is the cut-site analog.
- The scale factor counts distinct `bamtobed` names, which are reads
  (`/1`, `/2` suffixes), so "10M" means 10 million reads after all
  filters, not fragments. `README.md` says this correctly.

Fix (optional): since MACS2 already runs, adding `-B --SPMR` and the
two `bdgcmp` calls plus the vendored `bedGraphToBigWig` would give
ENCODE-comparable tracks at little cost. The signal tracks come from
the same `callpeak` invocation; the peak cap (C-04) applies to the
peaks only.

### C-13 No reproducibility analysis (Low, design)

ENCODE calls peaks on each replicate, on the pooled tagAlign, on two
pseudoreplicates per replicate and two pooled pseudoreplicates, then
computes the naive overlap set (50% reciprocal overlap) and the IDR
set (threshold 0.05), and reports rescue and self-consistency ratios
(standard: both below 2), with the overlap set as the default ATAC
output. Peak-count standards (replicated > 150,000, IDR > 70,000) refer
to those sets.

This pipeline merges every sample's blacklist-filtered peaks within
500 bp into one atlas and leaves replication to the count-based edgeR
step. IDR 2.0.3
is installed by `00.SETUP.sh` but nothing calls it. None of the
reproducibility figures can be produced.

Fix: out of scope for a count-matrix pipeline. The cheap partial step
is a `num_peak` metric per sample (ENCODE reports it) so that a failed
sample is visible before the atlas is built. Related: the deliverables
check in `pipe.sh:726-734` accepts empty peak files (`-e`, not `-s`;
`00.ISSUES.md` #11), so a sample with zero peaks still gives
`COMPLETED`. ENCODE's `call_peak` task fails on an empty peak file
(`assert_file_not_empty`).

### C-14 TSS enrichment (Low)

Concordant: `tss_enrich.R` reproduces `encode_task_tss_enrich.py`
(+/-2000 bp, 400 bins, shift `-read_len/2` by strand, Greenleaf
normalization on 10 edge bins each side, score = max of the mean
profile), and `hg38_tss.bed`, `b37_tss.bed` (hg19 with `chr`
stripped) and `mm10_tss.bed` are byte-identical to ENCODE's v1 files.

Variances:

- Annotation version. ENCODE's genome data default is now `VER=v3`,
  whose TSS files are `ENCFF493CCB` (hg38) and `ENCFF498BEJ` (mm10);
  the vendored files are the v1 GENCODE v24 / v19 / vM7 unique TSS.
  The portal thresholds (GRCh38 RefSeq: < 5 concerning, 5-7
  acceptable, > 7 ideal; mm10 RefSeq: < 10, 10-15, > 15) are stated
  per annotation. Scores against the v1 files are internally
  consistent with earlier deliveries but are not the numbers the
  current ENCODE pipeline would report. `R/TSSEnrich/README.md:72-73`
  states one rule for every genome, ">= 7 ideal, >= 5 acceptable,
  < 1 failed". The first two are the GRCh38 numbers only; mm10's are
  10 and 15, hg19's 6 and 10; and "< 1 failed" is not an ENCODE cut
  (ENCODE's lowest band is "concerning", below 5 for GRCh38).
- Odd read lengths. The task runs in ENCODE's Python 2 environment
  (`runtime_environment_python2` in `atac.wdl`), so `-read_len/2` is
  integer floor division: `-51` for 101 bp reads, applied as an
  integer to both strands by metaseq (`interval.start += shift` for
  `+`, `-= shift` for `-`). `tss_enrich.R:98` takes `-50.5` and
  truncates toward zero per strand, giving `-51` on `+` and `+50` on
  `-`. The minus-strand shift and the padded fetch window are 1 bp off
  for odd read lengths; even lengths are exact. The effect on a 10 bp
  bin is small and is probably part of the documented 0.85% residual
  (`R/TSSEnrich/README.md:84`: ENCODE's own `encode_task_tss_enrich.py`
  from the v2.2.3 container gave 29.5923 and the R port 29.8447 on the
  same BAM).
- Read extent. metaseq shifts the whole aligned span
  (`interval.start`, `interval.stop`); the R port uses CIGAR M/=/X
  blocks (`parse_cigars`). They differ only for reads with deletions
  or skips, which are rare in ATAC.
- Read length. ENCODE takes the maximum read length over the first
  1,000,000 FASTQ reads; `bin/computeTSSEnrich.sh:35-47` takes the modal
  `^[0-9]+M$` CIGAR length over the first 100,000 such reads. An
  unclipped read's M length is its sequenced length, so the mode over
  unclipped reads is the modal sequenced length however many other
  reads BWA soft-clips; the two definitions differ only when the FASTQ
  had reads of mixed length (then ENCODE takes the longest). A wrong
  read length shifts the whole profile, and the `peak_pos_bp` column of
  the CSV is the symptom: it should sit within a bin or two of 0.
- Input BAM. ENCODE's is chrM-free; here it keeps `chrM`. No effect on
  the score, since only TSS windows are read.

Fix: use `SHIFT <- -(read_len %/% 2)` with integer arithmetic for both
strands and drop the truncation comments; note the annotation version
and the mm10 thresholds in `R/TSSEnrich/README.md` and in deliveries.

### C-15 Fragment-length QC (Low)

ENCODE (`encode_task_fraglen_stat_pe.py`) runs Picard
`CollectInsertSizeMetrics W=1000 STOP_AFTER=5000000` on the nodup BAM
and then checks: fraction of fragments under 150 bp >= 0.4; NFR /
mononucleosome (150-300 bp) >= 2.5; a wavelet peak in 20-90 bp, in
120-250 bp and in 300-500 bp. The results go to `nucleosomal.qc` as
pass/fail lines. The portal standard requires a visible NFR and a
mononucleosome peak (147-294 bp).

This pipeline runs Picard on the full filtered BAM (no `STOP_AFTER`,
`chrM` included, fragments <= 30 bp removed) and plots the
distributions (`plotINSStats.R`). No fraction or pass/fail is
computed, so a sample without nucleosomal periodicity is only visible
to someone reading the PDF.

Fix: compute the NFR fraction, the NFR/mono ratio and the three peak
checks from the `___INS.txt` histogram in `plotINSStats.R` (or a small
R step) and write them next to the TSS CSV in `atacSeq/metrics`.

### C-16 ENCODE QC metrics not produced (Medium as a group)

Enabled by default in `atac.wdl` and absent here:

| ENCODE metric | Standard | Here |
| --- | --- | --- |
| Alignment stats, raw and filtered (`samstats`) | > 95% mapped | none; upstream |
| Duplicate metrics (`dup.qc`) | | written by PEMapper, not carried into the delivery (C-09) |
| NRF, PBC1, PBC2 (`lib_complexity.qc`) | > 0.9, > 0.9, > 3 | none (C-09) |
| `frac_mito` from the unfiltered alignments | | post-filter version (C-10) |
| Fragment-length pass/fail (`nucleosomal.qc`) | NFR and mono peaks present | histogram only (C-15) |
| TSS enrichment | > 7 (hg38), > 15 (mm10) | yes (C-14) |
| FRiP per peak set | > 0.3 | atlas fraction only (C-11) |
| Number of peaks, peak region size | > 150,000 replicated | atlas size only |
| JSD / fingerprint (deepTools) | | none |
| GC bias (Picard) | | none |
| Annotation enrichment (fraction of reads in DHS, promoters, enhancers) | | none |
| IDR / overlap reproducibility | rescue, self-consistency < 2 | none (C-13) |

Off by default in ENCODE and also absent: cross-correlation (`xcor`,
needs phantompeakqualtools), preseq.

### C-17 SAF coordinates are off by one (Low)

`mergePeaksToSAF.sh:44-48` writes the `bedtools merge` and `bedtools
subtract` output (BED, 0-based half-open start) straight into SAF
columns 3 and 4. This includes the second half of a merged peak split
by the blacklist subtraction (C-03), whose start is a blacklist end.
featureCounts reads SAF `Start` as 1-based inclusive, so every atlas
peak is tested as one base wider on the left than the merged region.
`assignedCheck.sh` compensates (`s = $3 - 1`) so Check B tests the
same off-by-one interval, and `R/diffAnalysisPairwise.R` writes
`Start` back out as a BED start, which is coordinate-correct only
because the stored value is in fact 0-based. The effect on counts is
negligible; the inconsistency matters if anyone treats
`peaks_raw_fcCounts.txt` Start as 1-based, which featureCounts
documentation says it is.

Fix: `awk '{print "Peak_"++s, $1, $2+1, $3, "+"}'` as the last step of
`mergePeaksToSAF.sh`, then `s = $3 - 1` in `assignedCheck.sh` is
already right, and `diffAnalysisPairwise.R` must subtract 1 when it
builds the BED for ChIPseeker. Change the three together and rerun
Check B.

### C-18 Atlas merging within 500 bp (Low, design)

ENCODE never merges peaks; its reproducible sets keep MACS2's peak
boundaries. `bedtools merge -d 500` fuses any peaks closer than 500 bp,
including distinct neighboring regulatory elements, and the result
grows with the number of samples. This is R. Koche's method
(`NOTES.md`, "any overlap, or even anything within 500 bp") and is a
deliberate design. Note it in deliveries; ENCODE peak-count standards
do not apply to the atlas.

### C-19 Tn5 shift guard (Very low)

ENCODE's `tn5_shift_ta` adds, after the +4/-5 shift:

```awk
if ($2 >= $3) { if ($6 == "+") {$2 = $3 - 1} else {$3 = $2 + 1} }
```

`postMapBamProcessing_ATACSeq.sh:147` lacks it. A read with an aligned
span of 5 bp or less would produce an inverted interval. BWA-MEM's
minimum seed makes this impossible in practice, but the guard costs
nothing.

### C-20 Tool versions (Low)

ENCODE pins (`scripts/requirements.txt`, `requirements.macs2.txt`):
MACS2 2.2.4, samtools 1.9, bedtools 2.29.0, Picard 2.20.7, bowtie2
2.3.4.3, cutadapt 2.5, IDR 2.0.4.2. This pipeline: MACS2 2.2.9.1,
samtools 1.20, bedtools and Picard from `~/bin` (versions not
recorded), IDR 2.0.3 (unused). No known behavioral change between
MACS2 2.2.4 and 2.2.9.1 for `callpeak` with these options, but the
versions should be recorded in `00.RUNSTATUS.txt` for reproducibility.

### C-21 Aligner, trimming and reference (Inherent)

ENCODE trims adapters with cutadapt (`-e 0.1 -m 5`) and aligns with
bowtie2 end-to-end, `-X2000`, `-k 5`, to the GRCh38 no-alt analysis
set (chrEBV included in v3) or `mm10_no_alt_analysis_set_ENCODE`.
This pipeline takes BWA-MEM BAMs from PEMapper on the GDC
`GRCh38.d1.vd1` (no-alt set plus decoys and viral contigs), UCSC mm10,
or b37 (which ENCODE does not support; the pipeline uses ENCODE's hg19
TSS and blacklist files with `chr` stripped, which apply by
coordinate).

Consequences that cannot be removed by this pipeline: BWA soft-clips
adapters instead of trimming them; proper-pair status is decided by
BWA's insert-size model rather than a fixed 2000 bp cap; reads that
map to the GDC decoys are lost to the primary contigs; MAPQ
distributions differ (C-05). hg38 main contigs are identical in name
and length to ENCODE's, so peak coordinates are directly comparable.

### C-22 Documentation out of step with the code (Low, docs)

- `NOTES.md:5` (R. Koche, 2018): "MACS2 calls on bam files containing
  both reads, but not setting the 'bampe' flag." The script calls
  MACS2 on the shifted BED with `-f BED`, which is the ENCODE choice
  and the right one; the note describes an older method. Add a dated
  remark, as the featureCounts paragraph of `NOTES.md` already has.
- `README.md:49-50`: duplicates marked "by the mapper". It is Picard
  MarkDuplicates in PEMapper; see C-09.
- `R/TSSEnrich/README.md:72-73`: one threshold rule for every genome,
  and a "< 1 failed" cut that ENCODE does not define; see C-14.

## 4. Concordant with ENCODE

- Flag filter `-F 1804` with proper pairs (`-f 3` is `-f 2` plus the
  paired bit, which proper pairs imply).
- Duplicates marked with Picard MarkDuplicates (upstream, in PEMapper)
  and removed by the 0x400 bit of `-F 1804`, the same tool and removal
  step as ENCODE's `filter` task.
- Mitochondrial reads out of peaks, signal and the atlas for every
  build; peak contigs limited to autosomes, X and Y (the ENCODE
  `chr[\dXY]+` filter, applied earlier).
- `bedtools bamtobed` to one line per read and the `+4 / -5` shift,
  same awk as `tn5_shift_ta`.
- MACS2 `-f BED --nomodel --shift -75 --extsize 150 --keep-dup all
  -p 0.01 --call-summits`, `-g hs` for hg38 and b37, `-g mm` for mm10,
  one call per sample on the full tagAlign with no subsampling
  (`subsample_reads = 0`). Since C-01 and C-02.
- Blacklist filter of each sample's peaks with ENCODE's blacklists,
  `intersect -v` and score cap at 1000. Since C-03.
- Before it, sort by `-log10 p`, top 300,000, clamp and `bedClip
  -truncate`, in ENCODE's order; the final per-sample peaks equal
  ENCODE's but for the peak names. Since C-04.
- TSS enrichment algorithm and the v1 ENCODE TSS files; `chrom.sizes`
  from the BAM header as ENCODE requires.
- Picard `CollectInsertSizeMetrics` as the fragment-length tool.

## 5. Outside ENCODE's scope

The merged peak atlas (`macsPeaksMerged.saf`), the featureCounts
matrix and its two checks, the DESeq2 size factors (inverted, per
`NOTES.md`), the insert-size and QC PDFs with PCA, and the edgeR
pairwise differential analysis have no counterpart in the ENCODE
pipeline, which stops at peaks, signal tracks and a QC report. They
were reviewed for internal correctness only; the one coordinate issue
found is C-17. The count checks (A and B) are sound and have no ENCODE
equivalent. One further bug outside ENCODE's scope is already tracked
as `00.ISSUES.md` #5: the three-argument (size-factor) mode of
`makeBigWigFromBEDZ.sh:39` reads the BED path as the scale factor
(`scaleFactor=$2`; it should be `$3`). `pipe.sh:573` passes two
arguments, so the delivered 10M-normalized bigWigs are not affected.

## 6. Recommended order of work

C-01 to C-04 are done (section 1, "Status"). The rest:

1. C-09 duplicate-marking check and NRF/PBC1/PBC2; C-10 raw-BAM
   `frac_mito`; C-11 ENCODE FRiP; C-15 fragment-length pass/fail.
   These are new metrics files, not changes to existing outputs.
2. C-17 SAF off-by-one (three files, one Check B rerun).
3. C-05, C-06, C-07, C-08: decide and document the read filter; a
   change here invalidates comparisons with earlier runs and should be
   a separate, named release.
4. C-14 integer shift in `tss_enrich.R`; README notes on annotation
   version and mm10 thresholds.
5. Record tool versions in `00.RUNSTATUS.txt` (C-20).
6. C-22 documentation fixes.

Add every change that alters peaks or counts to `CHANGELOG.md` under
Breaking changes, as was done for C-01 to C-04.

## 7. Objections considered and rejected

A second, independent review of the same code (2026-10-12) raised the
seven points below, which this report does not adopt. Each is stated
in full with the evidence against it, so no other document is needed
to follow the argument. Everything else that review found is already
part of the findings above.

| Objection | This report | Evidence |
| --- | --- | --- |
| The bigWig divisor counts fragments, because both mates share a read name, so the scale is "a factor of two" off from read ends | The divisor counts reads. `bedtools bamtobed` appends `/1` and `/2` to the names of paired reads, so the two mates are distinct values in column 4 and "10 million mapped reads" in `README.md` is correct | bedtools2 `src/bamToBed/bamToBed.cpp` lines 409-410 and 541-542 (`name += "/1"`, `name += "/2"`); check below |
| ENCODE's peak cap is 500,000 by pipeline input, 300,000 by task default | 300,000 is what runs: `atac.wdl` declares `Int cap_num_peak = 300000` and passes it to the task. `docs/input.md` says 500,000 and is out of date | `atac.wdl` line 192; `encode_task_macs2_atac.py` default 300000 |
| The TSS read length is unreliable because BWA soft-clips often and the unclipped reads may be a minority | The mode over unclipped reads is the sequenced length regardless of how many reads are clipped; an unclipped read's M length equals its sequenced length. It differs from ENCODE's only for mixed-length FASTQs. `peak_pos_bp` is still the right check | `bin/computeTSSEnrich.sh:35-47`; `encode_lib_genomic.py` `get_read_length` (max over 1M FASTQ reads) |
| Cross-correlation (NSC, RSC) is ENCODE QC the pipeline lacks | It is also absent from a default ENCODE ATAC run: `enable_xcor = false` in `atac.wdl` (v1 of the pipeline ran it). Listed in C-16 as off-by-default | `atac.wdl` line 144 |
| A minus-strand alignment in the first bases of a contig becomes a negative or empty interval after the Tn5 shift | Position on the contig does not matter. `+4` on a start and `-5` on an end can only invert an interval whose aligned span is 5 bp or less, which BWA-MEM's seed length rules out. The guard is still worth adding (C-19) | `postMapBamProcessing_ATACSeq.sh:147` |
| Fragments of 30 bp or less are "the short end of the nucleosome-free distribution", so the insert-size filter removes signal | Tn5 cannot produce fragments much under about 38 bp (the dimer footprint), so pairs with TLEN <= 30 are adapter dimers, mis-pairs or alignment artifacts, not accessible-chromatin signal. The filter is still a deviation from ENCODE and should be documented (C-07) | Adey et al. 2010, Genome Biology 11:R119; Buenrostro et al. 2013, Nature Methods 10:1213 |
| The absence of IDR, pooled and pseudoreplicate peak sets is a High-severity finding | Same facts, different weight. The pipeline's product is a count matrix, and ENCODE's reproducible peak set has no role in it; rated design-level (C-13). It becomes High only if a delivery is described as an ENCODE peak set | |

Check for the first row, run in any analysis directory; every name
ends in `/1` or `/2`, one per mate:

```bash
zcat out/*/*.shifted.bed.gz | head -2 | cut -f4
```

## 8. Sources

ENCODE-DCC/atac-seq-pipeline, `master`
(<https://github.com/ENCODE-DCC/atac-seq-pipeline>), read 2026-10-12;
the `atac.wdl` line numbers quoted (144, 192) refer to that revision:

- `atac.wdl` (defaults: `multimapping = 4`, `mapq_thresh = 30`,
  `filter_chrs = ['chrM','MT']`, `subsample_reads = 0`,
  `cap_num_peak = 300000`, `pval_thresh = 0.01`, `smooth_win = 150`,
  `idr_thresh = 0.05`, `enable_tss_enrich`, `enable_fraglen_stat`,
  `enable_jsd`, `enable_gc_bias`, `enable_annot_enrich`,
  `enable_idr` true; `enable_xcor`, `enable_preseq` false)
- `src/encode_task_filter.py`, `src/encode_lib_genomic.py`
  (`remove_chrs_from_bam`, `subsample_ta_pe`, `get_read_length`,
  `bed_clip`)
- `src/assign_multimappers.py`
- `src/encode_task_bam2ta.py`
- `src/encode_task_macs2_atac.py`
- `src/encode_task_post_call_peak_atac.py`,
  `src/encode_lib_blacklist_filter.py`, `src/encode_lib_frip.py`
- `src/encode_task_macs2_signal_track_atac.py`
- `src/encode_task_tss_enrich.py` and metaseq `array_helpers.py`
  (`_local_coverage`)
- `src/encode_task_fraglen_stat_pe.py`
- `src/encode_task_post_align.py` (`samtools flagstat` on the raw and
  the `chrM`-free BAM, read length from the FASTQ),
  `src/encode_task_frac_mito.py` (`Rm / (Rn + Rm)`)
- `src/encode_task_bowtie2.py` (`-k multimapping+1`, `-X2000`),
  `src/encode_task_trim_adapter.py` (`-e 0.1 -m 5`)
- `scripts/download_genome_data.sh` (blacklist, TSS and chrom.sizes
  URLs per genome and version), `scripts/requirements.txt`,
  `scripts/install_conda_env.sh` (UCSC tools v369, `bedClip`)
- `docs/input.md`
- <https://www.encodeproject.org/atac-seq/> and
  <https://www.encodeproject.org/data-standards/atac-seq/atac-encode4/>
  (data standards)
- bedtools2 `src/bamToBed/bamToBed.cpp` (read-name `/1`, `/2`
  suffixes)

Local: `NOTES.md`, `QC/QCNotes.md`, `docs/GENOMES.md` (blacklist
sources), `docs/SLURM_PORT.md` (validation runs Proj_18143 and
test25 to test40), `R/TSSEnrich/README.md` (validation of the R port
against the ENCODE v2.2.3 container), `bin/checkCounts.sh` (Check A
and B definitions), `00.ISSUES.md`, `CHANGELOG.md`.
