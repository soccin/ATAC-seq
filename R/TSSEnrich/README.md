# TSSEnrich

ENCODE-equivalent TSS enrichment scoring for ATAC-seq, in R.

`tss_enrich.R` is a line-by-line R port of the ENCODE ATAC-seq pipeline's
`encode_task_tss_enrich.py`. It exists so the pipeline can produce the
standard ENCODE TSS enrichment QC number without installing the ENCODE
WDL/Caper stack or running their 2 GB Singularity image on every sample.

`lib/` holds the reference data the script needs. The TSS BED there was **not**
generated here — it is ENCODE's own published file, downloaded once and
vendored verbatim. The chrom.sizes has a different origin and matters just as
much. See [Where lib/ came from](#where-lib-came-from) and
[Getting lib files for other genomes](#getting-lib-files-for-other-genomes).

## Contents

| Path | What it is |
| --- | --- |
| `tss_enrich.R` | The live script. R port of `encode_task_tss_enrich.py`. |
| `lib/hg38_tss.bed` | ENCODE hg38 TSS BED (GENCODE v24, unique protein-coding TSS). Vendored verbatim. |
| `lib/hg38.chrom.sizes` | Chromosome sizes for the **GDC GRCh38.d1.vd1** reference, i.e. what BIC hg38 BAMs are aligned to. **Not** ENCODE's or UCSC's `hg38.chrom.sizes`, despite the name. |
| `lib/b37_tss.bed`, `lib/b37.chrom.sizes` | ENCODE hg19 TSS BED with `chr` stripped; sizes from a b37 BAM header. Made by `lib/getB37.sh`. |
| `lib/mm10_tss.bed`, `lib/mm10.chrom.sizes` | ENCODE mm10 TSS BED (GENCODE vM7), rebuilt byte for byte; sizes from the IRIS mm10 dict. Made by `lib/getMM10.sh`. |

The superseded `v1/` scripts (`computeTSSEnrichment.R`, `encodeTSS.R`,
`TSSEnrichMethods.md`) are in `R/attic/TSSEnrich_v1/`. `TSSEnrichMethods.md`
explains why ATACseqQC scores come out inflated (~310).

## How it is invoked

The pipeline never calls `tss_enrich.R` directly. `bin/computeTSSEnrich.sh
GENOME BAM` is the wrapper, and `pipe.sh` fans it out per sample with the
run's genome tag (`b37`, `hg38`, `mm10`). The wrapper:

1. Checks the tag against `lib/genomes/GENOMES.tsv` (`bin/genomeInfo.sh`).
2. Derives read length from the modal `^[0-9]+M$` CIGAR over the first 100k reads.
3. Calls `tss_enrich.R` with `lib/${GENOME}.chrom.sizes` and `lib/${GENOME}_tss.bed`.

**The lib filenames are built by string interpolation from the build
tag.** Every supported genome (a row of `GENOMES.tsv`) must have both
files; `pipe.sh` stops before submitting anything if either is missing,
and so does the wrapper.

Standalone invocation for debugging:

```bash
Rscript R/TSSEnrich/tss_enrich.R \
    --nodup-bam sample.bam \
    --read-len 101 \
    --chrsz R/TSSEnrich/lib/hg38.chrom.sizes \
    --tss   R/TSSEnrich/lib/hg38_tss.bed \
    --out-dir out/tssEnrich/sample
```

Outputs, prefixed with the BAM basename: `*.tss_enrich.csv` (the QC values),
`*.tss_enrich.pdf` (mean profile), `*.large_tss_enrich.png` (per-TSS heatmap
sorted by signal).

### Duplicates must be REMOVED, not marked

`--nodup-bam` means what it says. A BAM where duplicates are only flagged
(`0x400`) gives inflated scores, because every PCR duplicate is counted as
signal. The script aborts if it sees duplicate-flagged reads. Produce the
input with `samtools view -F 1024 -b marked.bam > nodup.bam`.

### Algorithm parameters

+/-2000 bp window (4001 bp), 400 bins by linear interpolation, per-read strand
shift of `-read_len/2` with C-style truncation, Greenleaf normalization against
the 10 outermost bins on each side (100 bp), score = max of the normalized
mean profile.

### Thresholds

ENCODE's thresholds depend on the genome and the annotation:

| Genome | Concerning | Acceptable | Ideal |
| --- | --- | --- | --- |
| hg38 | < 5 | 5-7 | > 7 |
| b37 (hg19 TSS) | < 6 | 6-10 | > 10 |
| mm10 | < 10 | 10-15 | > 15 |

The ENCODE portal states them per annotation (RefSeq TSS for GRCh38
and mm10). The files in `lib/` are ENCODE's v1 GENCODE files (v24,
v19, vM7; see below), not the v3 files the current ENCODE pipeline
uses, so read the cut-offs as a guide. Scores are consistent with
earlier deliveries, but are not the numbers the current ENCODE
pipeline would report. ENCODE has no "failed" band.

## Where lib/ came from

The porting work was done outside this repo, on juno at
`/juno/bic/work/socci/Work/Users/ChandarS/LiuY16/TSSEnrich/`:

- `TSSTest/` — pulled `atac-seq-pipeline_v2.2.3.sif` from ENCODE and ran the
  real `encode_task_tss_enrich.py` on `ED1-1___MD_postProcess.bam` to get a
  reference score. `RUN_ENCODE_TSS.md` there records the exact download URLs.
- `Port/` — wrote and validated `tss_enrich.R` against that reference.
  Python 29.5923 versus R 29.8447 (0.85 % difference) on the same sample.
  `Port/docs/ENCODE_TSSEScore_PORTING_INFO.md` documents every non-obvious
  implementation decision and the debugging trail. Read that before touching
  the binning, shift, or padding code.

Commit `8a8e953` copied `tss_enrich.R` plus that project's `lib/` into this
repo. Commit `a10bb65` renamed `hg38_*` to `b38_*` to match the pipeline's
genome tag vocabulary of the time; v2.0.0 reverses that, back to `hg38_*`,
when the tag `b38` became `hg38`. Commit `2f0fed6` wired the stage into
`pipe.sh`.

Provenance of the two hg38 files, verified 2026-08-23:

`lib/hg38_tss.bed` is byte-identical to ENCODE's published file:

```bash
curl -sO https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_gencode_tss_unique.bed.gz
gunzip hg38_gencode_tss_unique.bed.gz
# 19,553 lines, md5 4f2803777a4e77dc592105f295d9886c  == lib/hg38_tss.bed
```

19,553 one-base BED6 records, one per unique protein-coding gene TSS, on the
24 primary `chr*` contigs, named by versioned Ensembl gene ID. Built by
ENCODE from GENCODE v24. The upstream file is dated 2019-07-30.

`lib/hg38.chrom.sizes` has a **different** origin, and this is the part that is
easy to get wrong later. Despite its name it is **not** ENCODE's or UCSC's
`hg38.chrom.sizes`, and must not be replaced by either: ENCODE's is the
195-sequence no-alt analysis set; UCSC's is 455 sequences. Neither matches BIC BAMs. The
vendored file is 2,779 sequences — 25 primary (`chr1`-`chr22`, `chrX`, `chrY`,
`chrM`), 42 `*_random`, 2,512 `chrUn_*` including the JTFH decoys, `chrEBV`,
and 199 named viral contigs (HPV, HIV, HCV, CMV, KSHV, HTLV) — sorted by
decreasing length. That is the **GDC `GRCh38.d1.vd1`** sequence set, which is
what `getGenomeBuildBAM.sh` calls `hg38`. The first attempt used the UCSC file
and it was replaced for this reason; the UCSC copy survives on juno as
`TSSTest/hg38.chrom.sizes.orig`.

**Rule: the TSS BED comes from ENCODE, the chrom.sizes must match the BAM
header.** The script inner-joins TSS contigs against the chrom.sizes contigs,
so a chrom.sizes that does not cover the BAM's naming scheme silently drops
TSS windows. Watch `n_tss_used` / `pct_tss_used` in the output CSV: if it is 0
or far below ~0.95, the contig names do not line up.

## Getting lib files for other genomes

All ENCODE ATAQC TSS BEDs live in one bucket with a uniform path. Source of
truth is `scripts/download_genome_data.sh` in
[ENCODE-DCC/atac-seq-pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline)
(the `TSS=` variable per genome). All URLs below were confirmed live
2026-08-23; on 2026-10-11 the whole bucket returned 404 (the hg38 and mm10
files and `mm10/mm10.chrom.sizes` alike), though the script still lists
these URLs. The md5s below still identify the files; the recipe at the end
of this page rebuilds them (see `lib/getMM10.sh`).

| Build | TSS BED URL | Lines | md5 (uncompressed) |
| --- | --- | --- | --- |
| hg19 | `.../hg19/ataqc/hg19_gencode_tss_unique.bed.gz` | 19,741 | `f8edff2c31743548f19ef3a47b5a7137` |
| hg38 | `.../hg38/ataqc/hg38_gencode_tss_unique.bed.gz` | 19,553 | `4f2803777a4e77dc592105f295d9886c` |
| mm9 | `.../mm9/ataqc/mm9_gencode_tss_unique.bed.gz` | 13,983 | `de5725e1b8964e1d7afb211ce49d6295` |
| mm10 | `.../mm10/ataqc/mm10_gencode_tss_unique.bed.gz` | 16,956 | `b0a24d5c8165f8782c5c6c6c98ab511f` |

Prefix is `https://storage.googleapis.com/encode-pipeline-genome-data`.
Annotation versions per ENCODE: hg19 = GENCODE v19, hg38 = v24, mm9 = vM1,
mm10 = vM7. All use UCSC `chr*` naming.

Newer ENCODE genome-data releases point at different TSS files for human.
Only use these if you deliberately want to match a newer ENCODE pipeline
version, not for consistency with existing deliveries:

- v2 hg38: `.../hg38/ataqc/tss.pc.gencode.v29.bed.gz`
- v3 hg38: `https://www.encodeproject.org/files/ENCFF493CCB/@@download/ENCFF493CCB.bed.gz`
- v3 mm10: `https://www.encodeproject.org/files/ENCFF498BEJ/@@download/ENCFF498BEJ.bed.gz`

### mm10 (added in v2.0.0)

Made 2026-10-11 by

```bash
R/TSSEnrich/lib/getMM10.sh /data1/core001/rsrc/genomic/bic/assemblies/M.musculus/mm10/mm10.dict
```

ENCODE's file could not be downloaded (404), so the script rebuilds it
from the GENCODE vM7 GTF (EBI FTP) with the recipe below followed by
`sort -k1,1 -k2,2n -u`. The result is byte-identical to ENCODE's
`mm10_gencode_tss_unique.bed` (16,956 lines, md5
`b0a24d5c8165f8782c5c6c6c98ab511f`, checked by the script). The recipe
alone gives 16,965 lines: 9 TSS are shared by two genes, and ENCODE kept
one, the first when the recipe's output is sorted by the whole line. A
single `sort -u` keeps a different gene for 5 of them, so the script
sorts twice.

`mm10.chrom.sizes` is the 66 contigs of the IRIS mm10 dict, the header
`getGenomeBuildBAM.sh` calls `mm10`, sorted by decreasing length. All
16,956 TSS are on its contigs, and its sizes for chr1-19, X and Y equal
`lib/genomes/mm10.genome`. No mm10 sample has been scored with it yet:
check `pct_tss_used` in the CSV of the first mm10 run.

### Other builds

A new build needs `<TAG>_tss.bed` and `<TAG>.chrom.sizes` named after its
tag, plus a row in `lib/genomes/GENOMES.tsv` (see `docs/GENOMES.md`); the
wrapper interpolates the tag. The chrom.sizes MUST come from a real BAM
header (or the reference `.dict`) of the build, not from ENCODE or UCSC,
or contig sets will not match.

Before dropping files in, check the BAM header naming. ENCODE TSS BEDs are all
`chr`-prefixed. If the BAM uses Ensembl-style names (`1`, `2`, ...), as b37
does, rewrite the BED to match:

```bash
sed 's/^chr//' hg19_tss.bed > b37_tss.bed
```

Mitochondrial naming differs as well (`chrM` versus `MT`), but that contig
carries no TSS in these files so it does not matter. `chrY` does carry TSS (64
in the hg38 file), so do not drop it. A whole-genome prefix mismatch, on the
other hand, zeroes out the score.

### Recipe: build a TSS BED from scratch for an unsupported genome

If ENCODE never published one (yeast, custom references, patched builds),
generate it from a GENCODE/Ensembl GTF using the same command the Kundaje lab
used to make the ENCODE files ([kundajelab/ataqc](https://github.com/kundajelab/ataqc)):

```bash
zcat $GTF |
grep -P '\tgene\t' |
grep 'protein_coding' |
grep -v 'level 3' |
awk -F '[\t|"]' '{ print $1"\t"$4"\t"$5"\t"$10"\t0\t"$7 }' |
awk -F '\t' 'BEGIN{ OFS="\t" } { if ($6=="+") { $3=$2-1; $2=$2-2 } else { $2=$3; $3=$3+1 } print }' |
sort -k1,1 -k2,2n > ${BUILD}_tss.bed
```

Scores from a self-built TSS set are not strictly comparable to ENCODE
thresholds, since the gene set and annotation version differ. Note in the
delivery which annotation was used.

## Related

- `bin/computeTSSEnrich.sh` — the wrapper that chooses the lib files
- `bin/getGenomeBuildBAM.sh` — where build tags are defined
- `docs/GENOMES.md` — genome tags and `lib/genomes/GENOMES.tsv`
- `QC/QCNotes.md` — which ATAC QC metrics matter and why
