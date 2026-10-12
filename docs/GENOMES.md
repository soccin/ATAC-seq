# Genome builds

Written 2026-10-11 for v2.0.0 (`fix/genome-tags`), from the ChIP-seq
v1.0.0 spec (sections 3.6 and 3.7), which makes the same change. The
two repos use the same tags, the same genome files and the same
`bin/getGenomeBuildBAM.sh` (ATAC copied ChIP's).

## The rule

A genome is named by its build tag everywhere: the run record, file
names and every argument. Anything else is an error that lists the
supported tags and exits non-zero.

| Tag | Assembly | Contig names | IRIS reference (PEMapper config) |
| --- | --- | --- | --- |
| `b37` | GRCh37 | `1`..`22`, `X`, `Y`, `MT` | `H.sapiens/b37` (`human_b37`) |
| `hg38` | GRCh38 | `chr1`..`chr22`, `chrX`, `chrY`, `chrM` | `H.sapiens/GRCh38.d1.vd1` (`human_GRCh38`) |
| `mm10` | GRCm38 | `chr1`..`chr19`, `chrX`, `chrY`, `chrM` | `M.musculus/mm10` (`mouse_mm10`) |

GRCh38 was tagged `b38` until v2.0.0. Nobody else uses that name: the
Broad and UCSC both call build 38 with `chr` contigs `hg38`. The GDC
`GRCh38.d1.vd1` main contigs are identical in name and length to UCSC
hg38's (checked 2026-10-11 against the IRIS dicts).

Nobody types the genome of a run's own data. `pipe.sh` detects it from
the BAM headers (`bin/getGenomeBuildBAM.sh`: md5 of the sorted `@SQ`
`SN`/`LN` lines), stops if the BAMs differ or the build is not
supported, passes it to every stage, and records it as `GENOME=` in
`00.RUNSTATUS.txt`. `R/diffAnalysisPairwise.R` reads it from there.

## One table

`lib/genomes/GENOMES.tsv` holds every genome-dependent value; a
supported genome is a row of it. Shell scripts read it through
`bin/genomeInfo.sh TAG FIELD`, R scripts through `genomeInfo()` in
`R/genomeInfo.R`. Both stop with an error for an unknown tag or field.

| Column | Use |
| --- | --- |
| `TAG` | the build tag; `genomeInfo.sh TAG TAG` is the supported-genome check |
| `SPECIES` | reports only; no code branches on it |
| `MACS_GSIZE` | MACS2 `-g` |
| `TXDB`, `ORGDB` | annotation packages for `diffAnalysisPairwise.R` |
| `ADD_CHR` | `yes` if peak contigs need `chr` added to match the UCSC TxDb (b37 only) |
| `MITO` | the mitochondrial contig, for the `___MITO.txt` metric |

The column names match ChIP-seq's where they overlap (ChIP also has
`HOMER`; ATAC has `MITO`).

Per tag there are also:

- `lib/genomes/<TAG>.genome`: the contigs analyzed, autosomes, X and Y,
  with their lengths (b37 24 lines, hg38 24, mm10 21). It is both the
  region filter (POST, CALLP, BW) and the chromosome sizes (BW). The
  scripts write the BED form (`contig 0 length`) to a temp file. Never
  give the 2-column file itself to `samtools view -L`: samtools reads
  the second column as a single position, and on a 488,265-read test BAM
  it selected 0 reads (the 3-column BED selected 447,141).
- `R/TSSEnrich/lib/<TAG>_tss.bed` and `<TAG>.chrom.sizes`: the TSS
  enrichment references. The chrom sizes cover the whole BAM header,
  not just the analyzed contigs (see `R/TSSEnrich/README.md`).

`pipe.sh` checks that all three files exist for the run's genome before
it submits anything.

To add a genome: its header md5 in `bin/getGenomeBuildBAM.sh` (keep the
file identical to ChIP-seq's), a row in `GENOMES.tsv`, `<TAG>.genome`,
and the two TSS files. No code change.

`getGenomeBuildBAM.sh` names more builds than are supported (`hg19`,
`b37_dmp`, `hg38+mm39`, ...), so that an unsupported BAM stops with a
named error rather than "unknown". It exits 1 when it cannot read the
header (a missing or non-BAM file); `pipe.sh` then stops with
"cannot identify the genome of [BAM]" in `MESSAGE`. UCSC hg38 (`H.sapiens/hg38`, md5
`492b7d52...`) is not mapped to `hg38`: no BAMs exist (no bwa index on
IRIS), and `hg38.chrom.sizes` lists the 2,779 GDC contigs.

## Mitochondrial reads

Before v2.0.0 the shifted BED kept the contigs of
`lib/genomes/<build>.genome.bed`, which listed `MT` for b37 and `chrM`
for mm10 but no `chrM` for b38. So b37 and mm10 runs kept mitochondrial
reads in peak calling, the bigWig normalization and the atlas, and b38
runs dropped them. Confirmed on a b37 run (Proj_18741, MCF7-DMSO-1): its
shifted BED has 4,082,846 reads on `MT`, 5.2% of the BAM.

Since v2.0.0 every build drops them, as the ENCODE ATAC-seq pipeline
does. The postProcess BAM still has every read; only the shifted BED,
and so peaks, bigWigs and the atlas, leave them out. POST writes each
sample's fraction to `out/<SID>/<SID>_postProcess___MITO.txt`
(`SAMPLE`, `MITO_CONTIG`, `READS`, `MITO_READS`, `FRAC_MITO`), counted on
the filtered BAM (MAPQ, flags, insert size) before the contig filter. It
is staged to `atacSeq/metrics` and is a required deliverable.

### What ENCODE does

Read 2026-10-11 from the `master` branches (ATAC pipeline v2.2.3, ChIP
pipeline v2.2.2):

- **ATAC-seq pipeline: removes mitochondrial reads before peak calling.**
  `atac.wdl` has `Array[String] filter_chrs = ['chrM', 'MT']`, with the
  help text "mitochondrial reads will be filtered out while filtering".
  `src/encode_task_filter.py` removes those contigs from the final
  (nodup) BAM, and that BAM becomes the TAG-ALIGN given to MACS2.
  `docs/input.md`: "Mitochondrial chromosomes are filtered out by
  default."
- **ChIP-seq pipeline: keeps them.** `chip.wdl` has
  `Array[String] filter_chrs = []` ("no filtering out of specific
  chromosomes").
- **Both** drop non-`chr[\dXY]+` contigs, so `chrM`, from the final
  `.bfilt` peak files (`regex_bfilt_peak_chr_name`, applied in
  `encode_lib_blacklist_filter.py`).
- **QC:** the ATAC pipeline reports `frac_mito_reads = mito / (mito +
  non-mito)` (`src/encode_task_frac_mito.py`), computed in the `align`
  task, before filtering.
- **Data standards page:** mentions mitochondrial reads only in the depth
  standard ("25 million non-duplicate, non-mitochondrial aligned reads");
  no threshold on the fraction.

Sources:

- <https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/atac.wdl>
- <https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_task_filter.py>
- <https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_task_frac_mito.py>
- <https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_lib_blacklist_filter.py>
- <https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/docs/input.md>
- <https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/chip.wdl>
- <https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/master/docs/input.md>
- <https://www.encodeproject.org/atac-seq/>

## Testing

hg38 is tested with full runs on the 11 downsampled GRCh38 BAMs
(`docs/SLURM_PORT.md`). b37 and mm10 are checked only with header-only
BAMs: detection, the table and the file check in `pipe.sh`. No b37 or
mm10 run has been made on this version (`00.ISSUES.md` #14).
