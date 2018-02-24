# ATAC-Seq pipeline

## Version 3.3 (2018-02-24)

### Notes from R.K.

- Pre-alignment filtering:

    - Trim raw read for quality (?) & adapter sequence

- Alignment:

    - bowtie2 or bwa mem, depending on what we’re trying to do

- Post-alignment filtering:

    - Mark Duplicates
    - MAPQ 10+

- BigWig file generation.

	- convert bam to bed
	- extend reads to library size in the 3’ direction for ChIP, no read extension for ATAC
	- compute density for bigwig formation
		- normalizing to 10 million mapped reads

