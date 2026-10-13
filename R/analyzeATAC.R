cat("## START: analyzeATAC.R\n")

args=commandArgs(trailing=T)
if(len(args)<1) {
    cat("\n   Usage: analyzeATAC.R SampleManifest.csv [RUNTAG]\n\n")
    quit(status=1)
}

fixSampleNames<-function(ss) {
    sampRename[str_remove(ss,"_postProcess.*") %>% basename] %>% unname
}

reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv,
              log_breaks(base = base), domain = c(1e-100, Inf))
}

pngCairo<-function(filename,width=14,height=8.5,pointsize=12,res=150) {

    png(filename,type="cairo",units="in",
        width=width,height=height,pointsize=pointsize,res=res)

}


mergePNGs<-function(fileSpec) {
    fileRe=gsub("_%\\d+d",".*",fileSpec)
    pdfFile=gsub("_%\\d+d.*",".pdf",fileSpec)
    system2("convert",c(sort(dir_ls(regex=fileRe)),pdfFile),stderr=cc("stderr","mergePNGs","convert",DATE()))
}

#halt("INCLUDE")

suppressPackageStartupMessages(library(AnnotationDbi))

suppressPackageStartupMessages(require(patchwork))
suppressPackageStartupMessages(require(scales))
suppressPackageStartupMessages(require(edgeR))
suppressPackageStartupMessages(require(ggrepel))
suppressPackageStartupMessages(require(ggsci))

suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(fs))
suppressPackageStartupMessages(require(openxlsx))

ds=read_tsv("peaks_raw_fcCounts.txt.summary",show_col_types = FALSE)

MANIFEST_FILE=args[1]
if(len(args)==2) {
    RUNTAG=args[2]
} else {
    RUNTAG=""
}

manifest=read_csv(MANIFEST_FILE,show_col_types = FALSE) %>% arrange(SampleID)

sampRename=manifest$SampleID
names(sampRename)=manifest$MapID

ds=ds %>%
    gather(Sample,Count,-Status) %>%
    mutate(Sample=fixSampleNames(Sample)) %>%
    mutate(Status=gsub("_.*","",Status)) %>%
    group_by(Sample,Status) %>%
    summarize(Counts=sum(Count),.groups="drop") %>%
    mutate(Status=ifelse(Status=="Assigned","InPeaks","Outside")) %>%
    mutate(Status=factor(Status,levels=c("Outside","InPeaks")))

pg0=ggplot(ds,aes(Sample,Counts,fill=Status)) +
    theme_light(base_size=16) +
    scale_fill_brewer(palette="Paired") +
    coord_flip()

pg1=pg0 + ggtitle("Mapped Reads in MACS Peaks") +
    geom_bar(stat="identity") +
    scale_y_continuous(labels = label_number(suffix = " M", scale = 1e-6))

pg2=pg0 + geom_bar(stat="identity",position="fill") +
    ggtitle("Fraction Reads in MACS Peaks") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ylab("Percentage")

dd=read_tsv("peaks_raw_fcCounts.txt",comment="#",show_col_types = FALSE)

peak.annote=dd %>% select(PeakNo=Geneid,Chr,Start,End,Strand,Length)

#
# Remove excluded points
#

manifest=manifest %>% filter(!grepl("^EXC",Group))

d=dd %>%
    select(PeakNo=Geneid,matches("out/")) %>%
    data.frame(check.names=F) %>%
    column_to_rownames("PeakNo")

colnames(d)=fixSampleNames(colnames(d))
d=d[,manifest$SampleID]
group=factor(manifest$Group)
y <- DGEList(counts=d,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

#
# Depth: each sample's fragments in peaks, the column sum of the count
# matrix (the counts the PCA uses). A sample much shallower than the
# others has noisier counts and separates from them on a PC of its own,
# whatever its biology: in Proj_18143 ED1-1, at 0.28 x the median, was
# alone on PC2, and any deep sample thinned to its depth did the same.
#
LOW_DEPTH=0.5
LOW_COUNT=5

fragsInPeaks=colSums(d)
depth=tibble(
        SampleID=colnames(d),
        Group=manifest$Group,
        FragmentsInPeaks=fragsInPeaks,
        RatioToMedian=round(fragsInPeaks/median(fragsInPeaks),3),
        FracKeptPeaksLT5=round(colMeans(y$counts<LOW_COUNT),3)
    ) %>%
    mutate(LowDepth=RatioToMedian<LOW_DEPTH)

lowDepth=depth %>% filter(LowDepth) %>% pull(SampleID)
if(len(lowDepth)>0) {
    cat("WARNING: low depth (fragments in peaks <",LOW_DEPTH,"x the median):",
        lowDepth,"\n")
}

#
# 34 group colors; scale_color_uchicago() alone has 9 and draws any
# further group in grey
#
groupColors=c(pal_uchicago("default")(9),pals::cols25())

pd1=depth %>%
    mutate(SampleID=fct_reorder(SampleID,FragmentsInPeaks)) %>%
    mutate(Label=ifelse(LowDepth,
        str_glue("{percent(FracKeptPeaksLT5,accuracy=1)}  LOW DEPTH"),
        percent(FracKeptPeaksLT5,accuracy=1))) %>%
    ggplot(aes(FragmentsInPeaks,SampleID,fill=Group)) +
    theme_light(base_size=12) +
    geom_col() +
    geom_vline(xintercept=LOW_DEPTH*median(fragsInPeaks),linetype=2) +
    geom_text(aes(label=Label),hjust=-0.1,size=3) +
    scale_fill_manual(values=groupColors) +
    scale_x_continuous(labels=label_number(suffix=" M",scale=1e-6),
        expand=expansion(mult=c(0,.25))) +
    labs(x="Fragments in peaks",y=NULL,
        subtitle=str_glue(
            "Label: share of the {format(nrow(y),big.mark=',')} peaks ",
            "used for the PCA with < {LOW_COUNT} fragments\n",
            "Dashed line: {LOW_DEPTH} x the median; a sample left of it ",
            "may separate on a PC of its own (see the PCA depth check)"))

depthTable=depth %>%
    arrange(FragmentsInPeaks) %>%
    transmute(
        Sample=SampleID,
        Group,
        `Fragments in peaks`=format(FragmentsInPeaks,big.mark=","),
        `x median`=sprintf("%.2f",RatioToMedian),
        "Peaks < {LOW_COUNT}":=percent(FracKeptPeaksLT5,accuracy=1),
        Depth=ifelse(LowDepth,"LOW",""))

pgDepth=(pd1 | wrap_elements(gridExtra::tableGrob(depthTable,rows=NULL,
            theme=gridExtra::ttheme_minimal(base_size=9)))) +
    plot_annotation(title="Sequencing depth in peaks",
        theme=theme(plot.title=element_text(size=16)))

#
# PCA of the samples: the samples are the observations, so each peak is
# centered across samples and PC1 is the largest difference between them.
# (Before v2.0.0 the peaks were the observations and the plot showed the
# rotation; its PC1 was the profile all samples share.)
#
pr=prcomp(t(cpm(y,log=T)),scale=F)
pct=round(100*summary(pr)$importance[2,1:2])
pcaAxes=labs(x=str_glue("PC1 ({pct[1]}%)"),y=str_glue("PC2 ({pct[2]}%)"))
dp=pr$x %>%
  data.frame %>%
  rownames_to_column("SampleID") %>%
  left_join(manifest,by = join_by(SampleID)) %>%
  left_join(depth %>% select(SampleID,FragmentsInPeaks),by = join_by(SampleID))

#
# Point size is the sample's fragments in peaks, so a sample apart
# because of its depth can be seen as such
#
depthSize=scale_size(name="Fragments\nin peaks",range=c(1.5,7),
    labels=label_number(suffix=" M",scale=1e-6))

pcaLabels=geom_label_repel(
        color="black",
        max.overlaps=Inf,
        min.segment.length = 0,
        size=3,
        force=16,
        max.time=10,
        max.iter=100000)

pp1=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=16) + geom_point(aes(size=FragmentsInPeaks),alpha=.6) + scale_color_manual(values=groupColors) + depthSize + pcaAxes
pp2=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=16) + geom_point(aes(size=FragmentsInPeaks),alpha=.6) + scale_color_manual(values=groupColors) + depthSize + pcaAxes +
    pcaLabels

#
# Depth check: the same PCA after every sample is thinned at random
# (fixed seed) to about the smallest sample's fragments in peaks. A
# sample apart from its group in the PCA as run but not here is apart
# because of its depth; a sample apart in both differs in its counts. A
# check only: nothing else uses the thinned counts.
#
# The thinning is binomial, each fragment kept with probability
# smallest / total, which is what sequencing less gives. Not
# edgeR::thinCounts(target.size=): it removes a multinomial draw, with
# replacement, and clips at 0, so a sample thinned 10-fold gets about 10
# times the variance; on Proj_18143 that put the deepest samples apart.
#
set.seed(20261013)
nThin=min(fragsInPeaks)
dThin=as.matrix(d)
dThin[]=rbinom(length(dThin),size=dThin,prob=rep(nThin/fragsInPeaks,each=nrow(dThin)))
yThin=DGEList(counts=dThin,group=group)
yThin=yThin[filterByExpr(yThin),,keep.lib.sizes=FALSE]
yThin=calcNormFactors(yThin)
prThin=prcomp(t(cpm(yThin,log=T)),scale=F)
pctThin=round(100*summary(prThin)$importance[2,1:2])
dpThin=prThin$x %>%
  data.frame %>%
  rownames_to_column("SampleID") %>%
  left_join(manifest,by = join_by(SampleID))

pp3a=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=12) + geom_point(size=2) + scale_color_manual(values=groupColors) + pcaAxes +
    pcaLabels + ggtitle("As run")
pp3b=ggplot(dpThin,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=12) + geom_point(size=2) + scale_color_manual(values=groupColors) +
    labs(x=str_glue("PC1 ({pctThin[1]}%)"),y=str_glue("PC2 ({pctThin[2]}%)")) +
    pcaLabels +
    ggtitle(str_glue("Every sample thinned to about {format(nThin,big.mark=',')} fragments in peaks"))
pp3=(pp3a | pp3b) + plot_layout(guides="collect") +
    plot_annotation(
        title="PCA depth check",
        subtitle=str_c(
            "Right: the PCA after every sample is thinned at random to about the depth of the smallest.\n",
            "A sample apart from its group on the left but not on the right is apart because of its depth;\n",
            "a sample apart on both differs in its counts. The thinned counts are used for this page only."),
        theme=theme(plot.title=element_text(size=16)))

#
# Project for output naming: the last directory of the working directory
# path named Proj_* or B-<digits>, else ATAC
#
pp=strsplit(getwd(),"/")[[1]]
projNo=tail(grep("^Proj_|^B-\\d+",pp,value=T),1)
if(len(projNo)==0) {
    projNo="ATAC"
}

pfile=cc(projNo,RUNTAG,"ATACSeqQC.pdf")
pdf(file=pfile,width=11,height=8.5)
print(pg2)
print(pg1)
print(pgDepth)
print(pp2)
print(pp1)
print(pp3)
dev.off()

write_csv(depth,cc(projNo,RUNTAG,"ATACSeqDepth.csv"))

cat("## END: analyzeATAC.R\n")

