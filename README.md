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


## MACS2, IDR Installation

Note we need python 3.8+ for idr. Python 3.8 on JUNO is broken (SSL junk) so use
private version. 3.10 does not work because MACS2 version checking is broken so
us 3.9.

In root of ATAC-seq repo

```{base}
/juno/work/bic/socci/opt/common/CentOS_7/python/python-3.9.7/bin/python3 -m venv venv
. venv/bin/activate
pip install --upgrade pip
pip install numpy
pip install MACS2
pip install matplotlib
cd code/idr
python3 setup.py install
cd ../..
deactivate
```
