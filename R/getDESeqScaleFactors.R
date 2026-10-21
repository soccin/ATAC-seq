cat("## START: getDESeqScaleFactors.R\n")

suppressPackageStartupMessages(require(DESeq2))
dd=read.delim("atacSeq/atlas/peaks_raw_fcCounts.txt",comment="#")
ds=as.matrix(dd[,7:ncol(dd)])
rownames(ds)=dd$Geneid
cond=factor(colnames(ds))
dds=DESeqDataSetFromMatrix(ds,colData=data.frame(Sample=cond),design=~Sample)
dds=estimateSizeFactors(dds)
sfm=data.frame(scaleFactor=1/sizeFactors(dds))
write.csv(sfm,"atacSeq/atlas/scaleFactorsDESeq2.csv")

cat("## END: getDESeqScaleFactors.R\n")
