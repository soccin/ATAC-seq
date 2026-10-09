# ATAC-Seq pipeline

## Version 1.1.0 - 2026-03-12

Single end version which uses both reads from PE-runs. Using methods from R.K. for bigWig generation.

Need to install MACS2 locally in venv. See below.

Runs on IRIS/Slurm. From the analysis directory, with `~/bin/sbatch` on
PATH (it exports `SBATCH_SCRIPT_DIR`, which `pipe.sh` needs to find its
own checkout; the partition, qos, time and memory are `#SBATCH` directives
in `pipe.sh`):

```
mkdir -p SLURM.CTRL
sbatch /path/to/ATAC-seq/pipe.sh [-q MAPQ] BAM1 [BAM2 ...]
```

- Post-alignment filtering:

    - Mark Duplicates
    - MAPQ 10+

- BigWig file generation.

	- convert bam to bed
	- extend reads to library size in the 3’ direction for ChIP, no read extension for ATAC
	- compute density for bigwig formation
		- normalizing to 10 million mapped reads


## Differential analysis

- If there are replicates you can do a differential peak analysis using: `R/diffAnalysisPairwise.R`. Usage:
```
Rscript R/diffAnalysisPairwise.R GENOME SampleManifest.csv Comparisons.csv [RUNTAG]
```

Example inputs:

```
  SampleManifest.csv

        SampleID,Group,MapID
        aNSC_loxp15_1,aNSC_loxp15,s_aNSC_loxp15_1
        aNSC_loxp15_2,aNSC_loxp15,s_aNSC_loxp15_2
        aNSC_p53_1,aNSC_p53,s_aNSC_p53_1

   Comparisons.csv

        aNSC_loxp15,aNSC_p53

   Sign convention X2-X1; e.g., aNSC_p53-aNSC_loxp15
```

## MACS2, IDR Installation

MACS2 2.2.9.1 and IDR 2.0.3 both build under the python 3.10 on the IRIS
PATH, so the venv uses whatever `python3` is current.

In root of ATAC-seq repo

```{base}
. 00.SETUP.sh
```

