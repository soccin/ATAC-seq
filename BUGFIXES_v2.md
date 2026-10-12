# Bugs to fix before v2.0.0

Open items only. Each one is deleted from here once it is fixed and
tested; the record goes to `CHANGELOG.md` and `attic/zzISSUES.md`.

## MACS2 peak calling does not match ENCODE

`callPeaks_ATACSeq.sh:86-96` takes its settings from the Kundaje lab
ATAC pipeline (cited at line 72), which became the ENCODE ATAC-seq
pipeline, but it got two of them wrong. Both are `pre-port`: they go
back to the first version of the script.

### 1. `--shift` has the wrong sign

The script computes `shiftsize=$((smooth_window / 2))` (75) and passes
`--shift $shiftsize`, i.e. `--shift 75`. The Kundaje and ENCODE
pipelines pass the negative value, `--shift -75`.

The input is the Tn5-shifted BED, one line per read, and MACS2 (`-f BED`,
`--nomodel`) takes each read's 5' end, which is the Tn5 cut site. It
moves that point by `--shift` toward the read's 3' end and then extends
`--extsize` bases. So:

| Setting | Region each read contributes |
| --- | --- |
| `--shift -75 --extsize 150` (ENCODE) | cut site -75 to +75: a 150 bp window centered on the cut site |
| `--shift 75 --extsize 150` (ours) | cut site +75 to +225: entirely on one side of the cut site |

With our setting, reads on the two strands are pushed away from the cut
site in opposite directions, so the signal around an accessible region
is spread out and off center. That affects the peaks, the summits, the
merged atlas and therefore every count.

### 2. `--keep-dup` is missing

The script does not set `--keep-dup`, so MACS2 uses its default, `1`: it
keeps one read per position and strand. ENCODE uses `--keep-dup all`.
Duplicates are already removed before this step (POST drops flag 1024,
set by MarkDuplicates on the pair). With `1`, MACS2 also drops distinct
fragments that happen to share a cut site, which happens often at
strong ATAC peaks, so it caps the signal where it is highest.

### Correct convention

ENCODE, `encode_task_macs2_signal_track_atac.py` (chip-seq-pipeline2,
used by atac-seq-pipeline):

```python
shiftsize = -int(round(float(smooth_win)/2.0))
...
'macs2 callpeak '
'-t {ta} -f BED -n {prefix} -g {gensz} -p {pval_thresh} '
'--shift {shiftsize} --extsize {extsize} '
'--nomodel -B --SPMR '
'--keep-dup all --call-summits '
```

Kundaje `modules/callpeak_macs2_atac.bds` (the source cited in our
script):

```
shiftsize := round( smooth_window.parseReal()/2.0 )
... --nomodel --shift -$shiftsize --extsize $smooth_window -B --SPMR --keep-dup all --call-summits
```

The fix in `callPeaks_ATACSeq.sh`. `-B --SPMR` only write the signal
bedGraphs, which this pipeline does not use (its bigWigs come from
`makeBigWigFromBEDZ.sh`), so they stay out:

```bash
smooth_window=150

#
# Negative: move each cut site back half the window, so --extsize
# centers a smooth_window-wide region on it (ENCODE:
# shiftsize = -int(round(float(smooth_win)/2.0))).
#
shiftsize=$(( -smooth_window / 2 ))

pval_thres=0.01

$MACS callpeak \
    -t $TDIR/cleanBED.bed.gz \
    -f BED \
    -n $PREFIX \
    -g $genome \
    -p $pval_thres \
    --nomodel \
    --shift $shiftsize \
    --extsize $smooth_window \
    --keep-dup all \
    --call-summits \
    --outdir $ODIR
```

### What the fix changes

Every MACS2 peak file, the merged atlas, the count matrix, the scale
factors, the QC PDF and any differential analysis. The bigWigs, BAMs,
BEDs, insert sizes, mitochondrial fraction and TSS enrichment do not
use MACS2 and do not change.

To do:

1. Make the change on its own branch and run the hg38 and b37 test sets;
   compare peak counts, peak widths, summit positions relative to TSSs,
   the atlas size and FRiP (fraction of fragments in peaks) with test26
   and test25.
2. Decide whether Proj_18143 is redelivered on v2.0.0 with this fix (a
   new run, and the difference report then covers it too) or on the
   v2.x-pre2 run already made.
3. Add both to `CHANGELOG.md` under Breaking changes: peaks from every
   earlier run were called with these settings.

Sources:

- <https://raw.githubusercontent.com/ENCODE-DCC/chip-seq-pipeline2/master/src/encode_task_macs2_signal_track_atac.py>
- <https://raw.githubusercontent.com/kundajelab/atac_dnase_pipelines/803e3d92435f52ef31b02fa264fecb5706d550e3/modules/callpeak_macs2_atac.bds>
