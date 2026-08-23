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
| `lib/b38_tss.bed` | ENCODE hg38 TSS BED (GENCODE v24, unique protein-coding TSS). Vendored verbatim. |
| `lib/b38.chrom.sizes` | Chromosome sizes for the **GDC GRCh38.d1.vd1** reference, i.e. what BIC b38 BAMs are aligned to. Not ENCODE's hg38 chrom.sizes. |
| `v1/computeTSSEnrichment.R` | Superseded. Custom method: strand-aware windows, 2000 bp / 10 bp bins, score from +/-200 bp of the TSS. Reference only, not called by the pipeline. |
| `v1/encodeTSS.R` | Superseded. Earlier, looser translation of the ENCODE Python (401 bins, 100 bp background). Reference only. |
| `v1/TSSEnrichMethods.md` | Notes comparing the v1 methods and ATACseqQC. Explains why ATACseqQC scores come out inflated (~310). |

## How it is invoked

The pipeline never calls `tss_enrich.R` directly. `bin/computeTSSEnrich.sh`
is the wrapper, and `pipe.sh` fans it out per sample. The wrapper:

1. Runs `bin/getGenomeBuildBAM.sh` on the BAM to get a build tag (`b38`, `mm10`, ...).
2. Derives read length from the modal `^[0-9]+M$` CIGAR over the first 100k reads.
3. Calls `tss_enrich.R` with `lib/${GENOME}.chrom.sizes` and `lib/${GENOME}_tss.bed`.

**The lib filenames are built by string interpolation from the detected build
tag.** That is the whole genome-support mechanism. To support a new build you
drop two correctly named files into `lib/` and nothing else needs to change.
Currently only the `b38` pair exists, so any other detected build fails at
this stage with a missing-file error.

Standalone invocation for debugging:

```bash
Rscript R/TSSEnrich/tss_enrich.R \
    --nodup-bam sample.bam \
    --read-len 101 \
    --chrsz R/TSSEnrich/lib/b38.chrom.sizes \
    --tss   R/TSSEnrich/lib/b38_tss.bed \
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
mean profile. ENCODE thresholds: >= 7 ideal, >= 5 acceptable, < 1 failed
library.

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
genome tag vocabulary. Commit `2f0fed6` wired the stage into `pipe.sh`.

Provenance of the two vendored files, verified 2026-08-23:

`lib/b38_tss.bed` is byte-identical to ENCODE's published file:

```bash
curl -sO https://storage.googleapis.com/encode-pipeline-genome-data/hg38/ataqc/hg38_gencode_tss_unique.bed.gz
gunzip hg38_gencode_tss_unique.bed.gz
# 19,553 lines, md5 4f2803777a4e77dc592105f295d9886c  == lib/b38_tss.bed
```

19,553 one-base BED6 records, one per unique protein-coding gene TSS, on the
24 primary `chr*` contigs, named by versioned Ensembl gene ID. Built by
ENCODE from GENCODE v24. The upstream file is dated 2019-07-30.

`lib/b38.chrom.sizes` has a **different** origin, and this is the part that is
easy to get wrong later. ENCODE's own `hg38.chrom.sizes` is the 195-sequence
no-alt analysis set; UCSC's is 455 sequences. Neither matches BIC BAMs. The
vendored file is 2,779 sequences — 25 primary (`chr1`-`chr22`, `chrX`, `chrY`,
`chrM`), 42 `*_random`, 2,512 `chrUn_*` including the JTFH decoys, `chrEBV`,
and 199 named viral contigs (HPV, HIV, HCV, CMV, KSHV, HTLV) — sorted by
decreasing length. That is the **GDC `GRCh38.d1.vd1`** sequence set, which is
what `getGenomeBuildBAM.sh` calls `b38`. The first attempt used the UCSC file
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
2026-08-23.

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

### Recipe: add mm10 support

```bash
cd R/TSSEnrich/lib

curl -sL -o mm10_tss.bed.gz \
  https://storage.googleapis.com/encode-pipeline-genome-data/mm10/ataqc/mm10_gencode_tss_unique.bed.gz
gunzip mm10_tss.bed.gz

# chrom.sizes MUST come from a real BAM of the build in question,
# not from ENCODE or UCSC, or contig sets will not match.
samtools view -H /path/to/an/mm10/sample.bam \
  | awk -F'\t' '$1=="@SQ" {
        sn=""; ln="";
        for (i=2; i<=NF; i++) {
            if ($i ~ /^SN:/) sn=substr($i,4);
            if ($i ~ /^LN:/) ln=substr($i,4);
        }
        if (sn != "" && ln != "") print sn"\t"ln
    }' \
  | sort -k2,2nr > mm10.chrom.sizes
```

That is all. `bin/computeTSSEnrich.sh` picks the files up by name once
`getGenomeBuildBAM.sh` reports `mm10`. Verify on one sample and check
`pct_tss_used` in the CSV before trusting the score.

The same recipe covers `mm9Full`, `GRC_m38`, and `mm10_hBRAF_V600E`, except
those tags need files named after the tag itself (`GRC_m38_tss.bed`, and so
on) — the wrapper interpolates the tag, not the species.

Before dropping files in, check the BAM header naming. ENCODE TSS BEDs are all
`chr`-prefixed. If the BAM uses Ensembl-style names (`1`, `2`, ...), as b37 and
b37_dmp do and as `GRC_m38` may, rewrite the BED to match:

```bash
sed 's/^chr//' hg19_tss.bed > b37_tss.bed
```

Mitochondrial naming differs as well (`chrM` versus `MT`), but that contig
carries no TSS in these files so it does not matter. `chrY` does carry TSS (64
in the b38 file), so do not drop it. A whole-genome prefix mismatch, on the
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
- `QC/QCNotes.md` — which ATAC QC metrics matter and why
- `CLAUDE.md` — genome tag vocabulary across the codebase
