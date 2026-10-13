# Versions to confirm for the client docs

`docs/methods.md` has a tool-version table. Some rows are taken from
notes, not from the pipeline (the pipeline does not record versions
per run; `00.FIXME.md` #24). Confirm these on IRIS from a session with
`bin/loadTools.sh` sourced and the pipeline venv active, then correct
the table and re-render with `docs/MAKE_DOCS.sh`.

## Confirm

| Row | Stated | How to check |
| --- | --- | --- |
| Picard SortSam, CollectInsertSizeMetrics | 3.4.0 | `picard SortSam --version` (the `~/bin/picard` wrapper); or the `# CollectInsertSizeMetrics` header of a `*___INS.txt` does not carry it, so check the jar the wrapper runs |
| bedtools | 2.31.1 | `bedtools --version` |
| UCSC wigToBigWig | v369 (assumed, same as bedClip) | `strings bin/wigToBigWig \| grep -i version`; or `bin/wigToBigWig` with no args prints the build in some versions |
| samtools | 1.20 | `samtools --version` after `module load samtools/1.20` |
| MACS2 | 2.2.9.1 | `venv/bin/macs2 --version` |
| featureCounts | 2.1.1 | `bin/featureCounts -v` |
| R | 4.5.1 | `R --version` |

## Add

The R package versions are not in the table. From the R on PATH:

```r
packageVersion("DESeq2")
packageVersion("edgeR")
packageVersion("ChIPseeker")
```

Add three rows (DESeq2, edgeR, ChIPseeker) to the table in
`docs/methods.md`. If `R/diffAnalysisPairwise.R` loads a build-specific
`TxDb` package for ChIPseeker, record that too.

## Then

- Remove the draft note in `docs/methods.md` that points here.
- Delete this file.
