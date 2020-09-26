Koche, Richard P./Sloan Kettering Institute <kocher@mskcc.org> wrote:

We haven't found one single 'one size fits all' metric for ATAC-seq success, but the combination of the following tends to work well, though I think you're already familiar with all of these based on our past discussions:


1. ratio of TSS height to flanking region (+/-5kb) for housekeeping genes. We have a list of mouse and human housekeeping genes, the majority of which tend to be highly accessible in nearly all cell and tissue types, so the mean composite signal is often helpful. I'm happy to give you the lists we use.


2. strand cross-correlation coefficient (though use only one read of PE alignment, otherwise will get artificial boost in signal). Shows that signal is not randomly distributed and looks 'peaky.'


3. simple fragment length distribution from bam (different from the above because looking at both read pairs). Expect to see some nucleosomal periodicity, which indicates the experiment worked.


[I like the above because none of them depend on peak calling, whereas the following do.]


4. the number of peaks. It's very simple, of course, but this is often very helpful, especially for samples with more noise.


5. SPoT: signal portion of tags, the number of reads in peaks


6. the inter and intra-group separation of samples in a PCA, based on the peak atlas of all samples in a data set. If one or more samples didn't work, or if the type of normalization used isn't fitting, it often shows here.


For code/packages:


1. custom

2. phantom peak qual tools

3. picard

4. macs2 or caller of choice, vs input (not computational background)

5. custom

6. R


The above are in addition to the basics for all seq data, like alignment stats. And note that any one metric above 'failing' doesn't mean the experiment didn't work, but if most fail then it usually does.


Still, there a few scenarios where the signal is borderline and it's unclear if it actually worked.
