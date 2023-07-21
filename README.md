# ATAC-Seq pipeline

## Version 4 - DEVS 2022-06-20

Single end version which uses both reads from PE-runs. Using methods from R.K. for bigWig generation.

Need to install MACS2 locally in venv. See below

- Post-alignment filtering:

    - Mark Duplicates
    - MAPQ 10+

- BigWig file generation.

	- convert bam to bed
	- extend reads to library size in the 3’ direction for ChIP, no read extension for ATAC
	- compute density for bigwig formation
		- normalizing to 10 million mapped reads


## Differential analysis

- If there are replicates you can do a differential peak analysis using: `R/diffAnalysisATACPairwise.R`. Usage:
```
analyzeATAC.R GENOME SampleManifest.csv Comparisons.csv [RUNTAG]
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

Note we need python 3.8+ for idr. Python 3.8 on JUNO is broken (SSL junk) so use
private version. 3.10 does not work because MACS2 version checking is broken so
us 3.9.

In root of ATAC-seq repo

```{base}
. CMDS.INSTALL.MACS
```

