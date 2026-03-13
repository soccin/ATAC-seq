# Changelog

All notable changes to this project will be documented in this file.

## [v1.1.0] — 2026-03-12

### Features

- **GRCh38 (b38) support**: The pipeline can now process human BAM files
  aligned to the GRCh38/GDC reference (chr-prefixed chromosome names).
  - Add MD5 checksum for GRCh38 GDC reference to `bin/getGenomeBuildBAM.sh`
    so the genome is auto-detected from the BAM header.
  - Add chromosome sizes files `lib/genomes/human_b38.genome` and
    `lib/genomes/human_b38.genome.bed` (chr1-chr22, chrX, chrY; chrM excluded).
  - Add `b38` case to `postMapBamProcessing_ATACSeq.sh`,
    `callPeaks_ATACSeq.sh`, and `makeBigWigFromBEDZ.sh`.

- **Genome-aware chromosome filtering**: Replace hardcoded `egrep -v` patterns
  (`chrUn|hs37d5|...`) with `bedtools intersect -nonamecheck -b $GENOME_BED`
  in `postMapBamProcessing_ATACSeq.sh`, `callPeaks_ATACSeq.sh`, and
  `makeBigWigFromBEDZ.sh`. Filters are now driven by the allowlist genome BED
  file rather than a denylist regex, eliminating spurious bedtools naming
  warnings from mixed-prefix genomes (e.g. GRCh38_GDC with CMV contig).

- **Genome parameter in postMapBamProcessing**: `postMapBamProcessing_ATACSeq.sh`
  now accepts a required `GENOME` positional argument (`b37 | b38 | mm10`).
  `pipe.sh` passes `$GENOME` automatically.

- **LSF job failure checking** (`bin/lsfTools.sh`): New sourceable utility
  library defining `bCheck JOBNAME`. Call immediately after `bSync` to abort
  the pipeline if any jobs in the group exited with a non-zero status.
  `bCheck` is now called after every `bSync` in `pipe.sh`
  (POST2, BW2, CALLP2, MergePeaks, Count, DESEQ).

- **hg38 implementation docs**: Add `docs/UPDATE_TO_B38.md` (full analysis of
  all files requiring changes) and `docs/CHECKLIST_B38.md` (atomic task
  checklist) to guide remaining b38 work.

### Fixes

- `bin/lsfTools.sh`: Fix `bCheck` false-exit under `set -e` — `grep -c`
  returns exit code 1 on zero matches, killing the pipeline even when all
  jobs succeeded. Changed to `grep -w "EXIT" | wc -l` which always exits 0.

### Refactoring

- Move `getGenomeBuildBAM.sh` to `bin/` alongside other utilities; update
  call site in `pipe.sh`.
- Move `lsfTools.sh` to `bin/`.
- Rename `CMDS.INSTALL.MACS` to `00.SETUP.sh`; update reference in `README.md`.
- `postMapBamProcessing_ATACSeq.sh`: deduplicate repeated usage message into
  a `usage()` function; fix script name in usage string.
- `R/analyzeATAC.R`: remove unimplemented differential analysis stub
  (edgeR contrasts, `doQLFStats`, hardcoded hg19 TxDb) that was guarded
  by `quit()` and never executed. Original archived to `R/attic/analyzeATAC.R`.
- Add BAM validation (`samtools quickcheck`) on the `INPUT_BAM` argument in
  `postMapBamProcessing_ATACSeq.sh` to catch argument-order mistakes early.
- `postMapBamProcessing_ATACSeq.sh`: derive sample ID (SID) and output
  paths from BAM `@RG SM:` header tag instead of stripping `.bam` suffix.
- `R/analyzeATAC.R`, `R/getDESeqScaleFactors.R`, `plotINSStats.R`: wrap
  all `library()`/`require()` calls in `suppressPackageStartupMessages`
  to suppress Bioconductor startup noise in pipeline logs.
- All pipeline R scripts: add `cat("## START/END: <script>\n")` messages
  for log tracing.

---

## [v1.0.1] — previous release (master)

See git log for earlier history.
