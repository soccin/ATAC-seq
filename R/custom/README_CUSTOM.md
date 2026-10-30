# Custom analyses

Scripts here are not part of the pipeline. Each was written for one
client request and is kept as a template for the next one like it.
None is run by `pipe.sh` or `deliverResults.sh`. Copy the script
into the project, edit its settings block, and run it there.

## plotESR1Locus.R: accessibility at one gene

Written for Proj_18143. The request was: "Would it be possible to
specifically examine chromatin accessibility at the ESR1 promoter,
gene body, and nearby regulatory regions across all six groups? If
available, normalized accessibility tracks or peak-level
visualizations would be especially helpful."

It replaces stepping through the locus in IGV with three files:

| File | Contents |
|---|---|
| `<PROJ>__<GENE>_Tracks.pdf` | Page 1: the gene +/- `FLANK`, mean coverage per group. Page 2: the same region, every sample. Then one page per cluster of the gene's TSSs (TSSs within `TSS_GAP` of each other), +/- `ZOOM_PAD`. Each page has coverage tracks on one y scale, the atlas peaks, the peaks with FDR < `FDR_CUT` in each comparison (colored by sign) and the gene's transcripts. |
| `<PROJ>__<GENE>_Peaks.pdf` | One panel per peak, log2 CPM per sample by group: the gene's promoter peaks first, then the peaks in the gene with FDR < `FDR_CUT` in any comparison, each set ordered by its smallest FDR. 6 panels per page. |
| `<PROJ>__<GENE>_Peaks.xlsx` | Every atlas peak in the gene +/- `FLANK`: location (promoter, exon, intron of any gene, or intergenic), distance to the nearest TSS of the gene, mean log2 CPM per group, logFC and FDR per comparison, log2 CPM per sample. The `Info` sheet defines the columns. |

### What it reads

All from the analysis directory of a completed run on which
`R/diffAnalysisPairwise.R` has also been run:

- `sampleManifest.csv`
- `atacSeq/bigwig/<MapID>_postProcess.shifted.10mNorm.bw`
- `atacSeq/atlas/peaks_raw_fcCounts.txt`
- `atacSeq/diff/<PROJ>__DiffPeaksEdgeR_V5.rds`. `<PROJ>` is the run
  tag given to `diffAnalysisPairwise.R`.
- `<ATAC_SDIR>/R/palettes.R`, from this repo

### Methods

- **Coverage:** the 10mNorm bigWigs (reads per 10M reads) as made by
  the pipeline, averaged in bins: at most `N_BINS` bins per panel,
  each at least `MIN_BIN` bp. Group tracks are the mean of the
  replicates.
- **Normalized counts (peak plot and xlsx):** edgeR log2 CPM (prior
  count 2), with TMM factors computed over all samples after
  `filterByExpr`, as `diffAnalysisPairwise.R` does before its
  per-comparison step. The comparisons' logFC and FDR use
  per-comparison TMM factors, so the two can differ slightly.
- **Gene models:** TxDb.Hsapiens.UCSC.hg38.knownGene. The track draws
  only the gene's transcripts, one row per distinct exon structure.
  The xlsx `Location` names any gene.
- **Promoter:** a peak within `PROMOTER` bp of a TSS.
- **Peak plot window:** the gene, from `PROMOTER` bp upstream of its
  first TSS to its 3' end. Peaks in the flanks are left out on
  purpose. A peak 50 kb away inside another gene is not an element of
  this gene without other evidence (e.g. chromatin interaction data),
  and a biologist reading a figure titled with the gene will not
  expect it there. The flanks remain in the tracks and the xlsx.

### Using it for another gene or project

1. Copy the script into the project, e.g.
   `atacSeq/custom/<GENE>/plot<GENE>Locus.R`. Putting it under
   `atacSeq/` means `deliverResults.sh` delivers it with the outputs;
   put it elsewhere if that is not wanted.
2. Edit the settings block at the top:
   - `PROJ`: the run tag of the diff results file
   - `GENE`: HGNC symbol
   - `GROUP_ORDER`: every `Group` of `sampleManifest.csv`, in plot
     order. The script stops if a group is missing.
   - `OUT_DIR`: the output folder
   - `ATAC_SDIR`: path from the analysis directory to this repo's
     checkout
   - `FLANK`, `TSS_GAP`, `ZOOM_PAD`, `PROMOTER`, `FDR_CUT` as needed.
     A gene with many alternative promoters gives many zoom pages;
     raise `TSS_GAP` to merge them.
   - Update the header comment (gene, file names).
3. Run it from the analysis directory as one Slurm task. Pass
   `-n 1`, or `srun` may start several copies that write the same
   files:

   ```bash
   srun -n 1 -N 1 -p cpushort,cmobic_short -c 2 --mem=16G -t 30 \
       Rscript atacSeq/custom/<GENE>/plot<GENE>Locus.R
   ```

   It takes about a minute. It prints the region, the number of atlas
   peaks in it and in the peak plot, and the significant peaks per
   comparison.
4. Look at every page before sending. Check that the peak plot holds
   only peaks the client asked about.

### Limits

- **hg38 only:** the TxDb and org.Hs.eg.db packages are written into
  the script. For another build, swap them for the `TXDB` and `ORGDB`
  of the build in `lib/genomes/GENOMES.tsv`.
- **BIC BAM names only:** the bigWig and count-matrix column names
  assume `<MapID>_postProcess`, not PEmap's `___MD`.
- **No enhancer-gene links:** distal elements are not assigned to the
  gene. To include them, add a public enhancer-gene map (e.g.
  GeneHancer, ENCODE cCREs) and select flank peaks with it.
- **One gene per run.**
