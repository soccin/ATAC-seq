# Pipeline notes from R. Koche

## 2018-02-27 Email

1) MACS2 calls on bam files containing both reads, but not setting the ‘bampe’ flag.

2) Merge all peaks (any overlap, or even anything within 500 bp, as we do) and calculating the raw counts matrix
Note: if using featureCounts and PE reads, set the ‘-p’ flag so fragments are counted instead of individual reads.
Also, in featureCounts, multi-mapped reads are filtered with the NH tag, which I don’t believe is set by bwa mem, but this shouldn’t be an issue because we’re filtering by MAPQ, etc anyway.

An example of how we create the raw counts matrix, based on a gtf or SAF (below) formatted merged peaks file:

subread-1.6.0-Linux-x86_64/bin/featureCounts -O -Q 10 -p -T 10 -F SAF -a  <PeakCoordinates_inSAFformat>  -o Peaks_raw_ftCounts.txt  <tab-separated list of all bam files to be processed>

3) Get the DESeq2 sizeFactors for the above matrix, and as you’ve already noticed, they’re the inverse of what we’re used to dealing with.
(and if curious, can make scatter plot of [inverse] size factor vs the sequencing depth-based normalization for comparison)

4) Bigwig file generation from PE read bams

A) This is the way we’ve traditionally done it, works very well, and even if it’s a bit tedious, it gives us control at every step. All tools are from bedtools except for wigToBigWig, which you can find here:

http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/

Here is a snippet of bash code that we use, with the bigwig created in a single line below:

sortedbam=<bamfile> [filtered for quality etc]
genome=<hg19, mm9, etc>
chromsizes=~/Annotation/Reference_Sizes/chromsizes_${genome}.tab
scalefactor=<factor from seq depth>
=> or sizeFactor scaling [once again, the inverse of the sizeFactor output from DESeq2]
bwout=${sortedbam/%.bam/.10mNorm.bw}
Or
bwout=${sortedbam/%.bam/.sizeFactorNorm.bw}

### ATAC (no read extension)

bamToBed -split -i $sortedbam |  slopBed -i stdin -g  $chromsizes -s -l 0 -r 0 | grep -vi "rand\|Un\|gl"  |  genomeCoverageBed -split -g  $chromsizes -i stdin -bg -scale $scalefactor |  wigToBigWig stdin $chromsizes $bwout

### Or for ChIP, with read extension (-r) accounting for avg frag length: 200 bp (or whatever is empirically determined)

bamToBed -split -i $sortedbam |  slopBed -i stdin -g  $chromsizes -s -l 0 -r 200 | grep -vi "rand\|Un"  |  genomeCoverageBed -split -g  $chromsizes -i stdin -bg -scale $scalefactor |  wigToBigWig stdin $chromsizes $bigwigout


B) the bamCoverage option in the deepTools package also converts bams to bigwigs. We’ve only played around with this thus far and have noted some minor discrepancies, so we’re sticking with the above for the time being, though they largely produced similar results. We were only testing, so I don’t have a set of standard options we’re using, so I’d just follow the examples here if you want to use this:

http://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html
