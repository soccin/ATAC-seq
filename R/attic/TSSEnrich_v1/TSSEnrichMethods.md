# TSS Enrichment Methods Comparison

## ENCODE Exact Definition

Source: `encode_task_tss_enrich.py`, function `make_tss_plot()`, in the ENCODE DCC
ATAC-seq pipeline (`github.com/ENCODE-DCC/atac-seq-pipeline`).

### Parameters (hardcoded defaults)
```python
bins=400, bp_edge=2000   # 400 bins of 10 bp over +-2000 bp
greenleaf_norm=True
```

### Read shift
```python
shift_width = -read_len / 2   # e.g. -25 for 50 bp reads
```
All reads shifted by the same amount regardless of strand. Centers each read on its
midpoint (approximates fragment center). NOT the Tn5 +4/-5 bp strand-specific shift.

### Strand handling
`stranded=True` passed to metaseq. Negative-strand windows are reverse-complemented
before aggregation so that position 0 always = TSS.

### Normalization (Greenleaf method)
```python
num_edge_bins = int(100 / (2*bp_edge/bins))   # = int(100/10) = 10 bins = 100 bp
bin_means     = bam_array.mean(axis=0)         # mean across all TSSs per bin
avg_noise     = (sum(bin_means[:10]) + sum(bin_means[-10:])) / (2*10)
bam_array    /= avg_noise
```
Background = mean of the **outermost 100 bp on each side** (10 bins x 10 bp).
The entire 2D array (one row per TSS) is divided by this single scalar.

### Score
```python
tss_point_val = max(bam_array.mean(axis=0))
```
Score = **maximum of the per-bin mean profile** (after normalization) over the
**entire +-2000 bp window**. No smoothing. No restriction to +-200 bp.

---

## Method Comparison

| | ENCODE exact | Our BED method | ATACseqQC default |
|---|---|---|---|
| Input | shifted BED/BAM | shifted BED | BAM |
| Window | +-2000 bp | +-2000 bp | +-1000 bp |
| Bin size | 10 bp | 10 bp | 100 bp |
| Strand-aware | yes | yes | no |
| Read shift | -read_len/2 (all reads) | Tn5 +4/-5 bp (strand-specific) | internal |
| Background | outermost **100 bp** each end (10 bins) | outermost **200 bp** each end (20 bins) | 100 bp outside +-1000 bp |
| Score region | max over **all 400 bins** | max within **+-200 bp** only | max of loess-smoothed ratios |

Differences between ENCODE and our BED method are minor. Both give ~15 for BD1-1.
A score of 310 from ATACseqQC with matched 10 bp bins is a normalization artifact:
individual edge bins near +-2000 bp have near-zero coverage, collapsing the
per-bin background denominator toward zero.

## Notes

**Shift sensitivity:** ATACseqQC uses 100 bp bins. The Tn5 shift is only 4-5 bp, which
is completely invisible at that resolution -- reads land in the same bin whether shifted
or not. Our 10 bp bins are sensitive to this shift. This explains why shifting has no
effect on ATACseqQC scores.

## Scores for BD1-1 (reference sample)

| Method | Score |
|---|---|
| Our BED method (`computeTSSEnrichment.R`) | 15.03 |
| ATACseqQC default (unshifted BAM) | 14.15 |
| ATACseqQC default (shifted BAM) | 14.29 |
| ATACseqQC matched params (10 bp bins) | 310.19 (artifact) |
| ENCODE-exact R (`encodeTSS.R`) | ~15 (expected) |

## Scripts

- `R/computeTSSEnrichment.R` -- our BED method (background 200 bp, score +-200 bp)
- `R/encodeTSS.R` -- ENCODE-exact reimplementation (background 100 bp, score all bins)
- `R/validateTSS.R` -- Tests 1, 2, 4: profile plots and parameter-matched ATACseqQC run
- `R/syntheticTSS.R` -- Test 3: synthetic data ground truth test
- `R/testTSS.R` -- ATACseqQC method wrapper
