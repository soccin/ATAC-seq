# ATAC-seq Pipeline: GRCh38 (hg38/B38) Support - Analysis Report

## Context

The current ATAC-seq pipeline supports human genomes b37/hg19 (GRCh37) and mouse mm10. This report identifies all changes needed to add GRCh38 (hg38/B38) support. A key difference between b37 and hg38 is chromosome naming: b37 uses bare numbers (`1, 2, ... X, Y, MT`) while hg38 uses the "chr" prefix (`chr1, chr2, ... chrX, chrY, chrM`).

Additional note, filter out chrM going forward. Calling events on it are likely to be broken/biased  incorrect without additional processing

---

## Files Requiring Changes

### 1. `getGenomeBuildBAM.sh` — Genome Detection via MD5

**Issue:** Genome is detected by computing an MD5 of the BAM file header and matching it against a lookup table. There is no entry for any hg38 reference.

**Change needed:**
- Determine the MD5 of the `.dict` file for the hg38 reference in use at your institution.
  Example: `md5sum /path/to/hg38.dict` (or extract the MD5 field from the BAM header dict entry).
- Add a new `case` block (after line ~70) mapping that MD5 to `hg38`:

```bash
"<MD5_OF_HG38_DICT>")
    echo "hg38"
    ;;
```

**Action required from user:** Provide the path to the institutional hg38 reference dict, or run:
```bash
samtools view -H <sample.bam> | grep "^@SQ" | md5sum
```
on a known hg38 BAM to compute the checksum.

---

### 2. `callPeaks_ATACSeq.sh` — MACS2 Peak Calling

**Issues:**
1. No `hg38` case in the genome-to-MACS2-code mapping (lines 58-79).
2. **REVIEW CAREFULLY [nds]**: Chromosome exclusion filter at line 53 is b37-specific:
   `chrUn|chrM|_random|_unplaced|GL|NC_|hs37d5`
   — This works for hg38 too (hg38 uses similar non-standard chr naming), but `hs37d5` is b37-
   specific and harmless for hg38. The pattern should be made genome-aware.

**Changes needed:**

Add to the genome case statement:
```bash
hg38)
    genome=hs
    ;;
```

Update chromosome exclusion filter to be genome-aware (or extend the current pattern to cover hg38 alternative sequences like `chrUn|_random|_alt|EBV`):
```bash
# hg38 non-standard contigs
chrUn|_random|_alt|EBV|NC_|GL
```

---

### 3. `makeBigWigFromBEDZ.sh` — BigWig Generation

**Issues:**
1. No `hg38` case in the genome-to-chrom-sizes mapping (lines 38-59).
2. Chromosome filter at lines 74-76 is b37-specific (same issue as above).
3. No chromosome sizes file for hg38 exists in `lib/genomes/`.

**Changes needed:**

Add to the genome case statement:
```bash
hg38)
    GENOME=$SDIR/lib/genomes/human_hg38.genome
    ;;
```

Update chromosome filter to handle hg38.

---

### 4. `postMapBamProcessing_ATACSeq.sh` — BAM Post-Processing / Tn5 Shift

**Issue:** Chromosome exclusion filter at lines 73-78 is hardcoded:
`chrUn|_random|GL|NC_|hs37d5|_unplaced`

The term `hs37d5` is b37-specific. For hg38, alternative sequences use `_alt` suffix. The filter needs to be updated to be genome-aware, or the pattern should be extended to cover both builds.

**Change needed:** Make the filter genome-aware by passing the genome build into the script and selecting the appropriate exclusion pattern:

```bash
# b37/hg19 filter:
chrUn|chrM|_random|GL|NC_|hs37d5|_unplaced

# hg38 filter:
chrUn|chrM|_random|_alt|_decoy|EBV|NC_|GL
```

This script currently does NOT take a genome build parameter — it would need one added.

---

### 5. `lib/genomes/human_hg38.genome` — NEW FILE NEEDED

**Issue:** No chromosome sizes file exists for hg38.

**Change needed:** Create `lib/genomes/human_hg38.genome` with standard hg38 chromosome sizes (chr1-chr22, chrX, chrY, chrM with "chr" prefix). Also create the `.bed` variant.

Get rid of `chrM`

Standard hg38 main chromosome sizes (GRCh38.p14):
```
chr1    248956422
chr2    242193529
chr3    198295559
chr4    190214555
chr5    181538259
chr6    170805979
chr7    159345973
chr8    145138636
chr9    138394717
chr10   133797422
chr11   135086622
chr12   133275309
chr13   114364328
chr14   107043718
chr15   101991189
chr16   90338345
chr17   83257441
chr18   80373285
chr19   58617616
chr20   64444167
chr21   46709983
chr22   50818468
chrX    156040895
chrY    57227415
```

---

### 6. `R/analyzeATAC.R` — Peak Annotation (HARDCODED hg19)

**Issue:** Line 128 hardcodes `TxDb.Hsapiens.UCSC.hg19.knownGene`. There is no genome parameter. This script cannot annotate hg38 peaks.

**Changes needed:**
- Add a command-line genome argument (e.g., `--genome hg38`)
- Add conditional logic to select annotation database:
  - `hg19`/`b37`: `TxDb.Hsapiens.UCSC.hg19.knownGene`
  - `hg38`: `TxDb.Hsapiens.UCSC.hg38.knownGene`
  - `mm10`: `TxDb.Mmusculus.UCSC.mm10.knownGene`
- Ensure `TxDb.Hsapiens.UCSC.hg38.knownGene` R package is installed.
- **Line 170:** `filter(Chr!="MT")` does not exclude the hg38 mitochondrial chromosome. In
  hg38 the chromosome is named `chrM`, so this filter does nothing and chrM peaks pass through
  unchecked. The filter must be made genome-aware (e.g., exclude both `MT` and `chrM`, or use
  a regex that matches either).
- **Line 174:** `mutate(V5=-10*log10(V5), V1=paste0("chr",V1))` unconditionally prepends
  `"chr"` to chromosome names. For b37/hg19 (bare names like `1`, `2`) this is correct, but
  for hg38 (names already prefixed `chr1`, `chr2`) it produces `"chrchr1"`, `"chrchr2"`, etc.,
  causing ChIPseeker to fail or return empty annotations. The fix must mirror the conditional
  check already in `diffAnalysisPairwise.R` (lines 386-393), which tests
  `substr(Chr[1],1,3) != "chr"` before adding the prefix.

---

### 7. `R/diffAnalysisPairwise.R` — Differential Analysis Annotation

**Issue:** The `"human"` genome option (lines 190-192) is hardcoded to hg19: `TxDb.Hsapiens.UCSC.hg19.knownGene`

**Changes needed:**
- Rename genome parameter from `"human"/"mouse"` to explicit build names: `"hg19"`, `"hg38"`, `"mm10"`
  (or add a separate `"hg38"` case)
- Add `hg38` case:
  ```r
  } else if (GENOME == "hg38") {
      txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
      annoDb <- "org.Hs.eg.db"
  }
  ```
- Update chromosome handling: hg38 has "chr" prefix natively (no stripping needed).

---

### 8. `attic/HOMER/doHomer.sh` — HOMER Motif Analysis (Archived)

**Issue:** Line 14 hardcodes `hg19` in the HOMER call:
`findMotifsGenome.pl $HOMERFILE hg19 outHomer/...`

**Change needed** (low priority — archived script):
- Accept genome as a parameter and pass to `findMotifsGenome.pl`
- HOMER natively supports `hg38` as a genome identifier

---

### 9. `pipe.sh` — Main Pipeline Orchestration

**Issues:**
1. **Line 91:** Does not pass genome build to `postMapBamProcessing_ATACSeq.sh`. Once that
   script gains a `-g`/`--genome` parameter (see Section 4), `pipe.sh` must pass it:
   ```bash
   $SDIR/postMapBamProcessing_ATACSeq.sh -q $MAPQ -g $GENOME
   ```
2. **Line 143:** Does not pass genome to `analyzeATAC.R`. Once that script gains a genome
   parameter (see Section 6), `pipe.sh` must pass it:
   ```bash
   Rscript $SDIR/R/analyzeATAC.R sampleManifest.csv $GENOME
   ```
   (or `--genome $GENOME` depending on chosen argument parsing approach)

---

## Summary of Required Changes

| File | Change Type | Priority |
|------|-------------|----------|
| `getGenomeBuildBAM.sh` | Add hg38 MD5 checksum entry | CRITICAL — blocks all hg38 detection |
| `lib/genomes/human_hg38.genome` | Create new file | CRITICAL — needed for BigWig |
| `lib/genomes/human_hg38.genome.bed` | Create new file | CRITICAL — needed for BigWig |
| `makeBigWigFromBEDZ.sh` | Add hg38 case + filter update | HIGH |
| `callPeaks_ATACSeq.sh` | Add hg38 case + filter update | HIGH |
| `postMapBamProcessing_ATACSeq.sh` | Add genome param + hg38 filter | HIGH |
| `pipe.sh` | Pass `$GENOME` to postMap and analyzeATAC.R calls | HIGH |
| `R/analyzeATAC.R` | Add genome param + hg38 annotation + chrM filter + chr prefix fix | HIGH |
| `R/diffAnalysisPairwise.R` | Add hg38 annotation option | HIGH |
| `attic/HOMER/doHomer.sh` | Parameterize genome | LOW (archived) |


---

## Blocker

The single most critical piece of information needed before implementation can begin:

**The MD5 checksum of the hg38 reference dictionary used at your institution.**

Run this on a known hg38 BAM file to get it:
```bash
samtools view -H <sample_hg38.bam> | grep "^@SQ" | md5sum
```

Or on the dict file directly:
```bash
md5sum /path/to/hg38.dict
```

Without this, `getGenomeBuildBAM.sh` cannot auto-detect hg38 BAM files and the whole pipeline will fail to recognize hg38 inputs.
