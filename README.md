# ATAC-Seq pipeline

## Version 4 / JUNO (2020-09-26)

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


## MACS2 Installation

In root of ATAC-seq repo

```{base}
python3 -m venv venv
. venv/bin/activate
pip install --upgrade pip
pip install numpy
pip install MACS2
deactivate
```
