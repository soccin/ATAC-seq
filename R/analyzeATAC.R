cat("## START: analyzeATAC.R\n")

args=commandArgs(trailing=T)
if(len(args)<1) {
    cat("\n   Usage: analyzeATAC.R SampleManifest.csv [RUNTAG]\n\n")
    quit()
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

suppressPackageStartupMessages(library(ChIPseeker))
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

pr=prcomp(cpm(y,log=T),scale=F)
dp=pr$rotation %>%
  data.frame %>%
  rownames_to_column("SampleID") %>%
  left_join(manifest,by = join_by(SampleID))

pp1=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=16) + geom_point(size=4,alpha=.6) + scale_color_uchicago()
pp2=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=16) + geom_point(size=2) + scale_color_uchicago() +
    geom_label_repel(
        color="black",
        max.overlaps=Inf,
        min.segment.length = 0,
        size=3,
        force=16,
        max.time=10,
        max.iter=100000)

pp=strsplit(getwd(),"/")[[1]]
projNo=grep("^Proj_|^B-\\d+",pp,value=T)

pfile=cc(projNo,RUNTAG,"ATACSeqQC.pdf")
pdf(file=pfile,width=11,height=8.5)
print(pg2)
print(pg1)
print(pp2)
print(pp1)
dev.off()

cat("## END: analyzeATAC.R\n")

