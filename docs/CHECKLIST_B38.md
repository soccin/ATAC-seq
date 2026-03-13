# hg38 Implementation Checklist

Derived from `UPDATE_TO_B38.md`. Each item is one atomic change.

---

## Prerequisite (Blocker)

- [ ] Obtain MD5 of institutional hg38 reference dict:
      `samtools view -H <hg38.bam> | grep "^@SQ" | md5sum`

---

## CRITICAL

### `getGenomeBuildBAM.sh`

- [ ] Add hg38 MD5 `case` block mapping checksum to `"hg38"` (after line ~70)

### `lib/genomes/`

- [ ] Create `lib/genomes/human_hg38.genome` (chr1-chr22, chrX, chrY — no chrM)
- [ ] Create `lib/genomes/human_hg38.genome.bed`

---

## HIGH

### `callPeaks_ATACSeq.sh`

- [ ] Add `hg38) genome=hs ;;` to genome case statement (lines 58-79)
- [ ] Update chromosome exclusion filter (line 53) to be genome-aware

### `makeBigWigFromBEDZ.sh`

- [ ] Add `hg38) GENOME=$SDIR/lib/genomes/human_hg38.genome ;;` to genome case (lines 38-59)
- [ ] Update chromosome filter for hg38 (lines 74-76)

### `postMapBamProcessing_ATACSeq.sh`

- [ ] Add `-g`/`--genome` parameter to script
- [ ] Update chromosome exclusion filter (lines 73-78): replace `hs37d5` with `_alt|_decoy` for hg38

### `pipe.sh`

- [ ] Line 91: Add `-g $GENOME` to `postMapBamProcessing_ATACSeq.sh` call
- [ ] Line 143: Pass `$GENOME` to `Rscript R/analyzeATAC.R` call

### `R/analyzeATAC.R`

- [ ] Add `--genome` command-line argument
- [ ] Line 128: Replace hardcoded `TxDb.Hsapiens.UCSC.hg19.knownGene` with genome-conditional TxDb selection
- [ ] Install `TxDb.Hsapiens.UCSC.hg38.knownGene` R package
- [ ] Line 170: Fix `filter(Chr!="MT")` to also exclude `chrM` (hg38 mitochondrial name)
- [ ] Line 174: Fix `paste0("chr",V1)` to add prefix conditionally only when not already present

### `R/diffAnalysisPairwise.R`

- [ ] Lines 190-192: Add `hg38` TxDb case (currently `"human"` maps only to hg19)
- [ ] Update chromosome handling for hg38 (chr prefix present natively, no stripping needed)

---

## LOW (archived)

### `attic/HOMER/doHomer.sh`

- [ ] Line 14: Replace hardcoded `hg19` with genome parameter passed to `findMotifsGenome.pl`
