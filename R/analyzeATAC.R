args=commandArgs(trailing=T)
if(len(args)<1) {
    cat("\n   Usage: analyzeATAC.R SampleManifest.csv [RUNTAG]\n\n")
    quit()
}

fixSampleNames<-function(ss) {
    sampRename[gsub("_postProcess.*","",ss) %>% gsub(".*_s_","s_",.)] %>%
        unname
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

library(ChIPseeker)
library(AnnotationDbi)

require(patchwork)
require(scales)
require(edgeR)
require(ggrepel)

require(tidyverse)
require(fs)
require(openxlsx)

ds=read_tsv("peaks_raw_fcCounts.txt.summary")

MANIFEST_FILE=args[1]
if(len(args)==2) {
    RUNTAG=args[2]
} else {
    RUNTAG=""
}

manifest=read_csv(MANIFEST_FILE) %>% arrange(SampleID)

sampRename=manifest$SampleID
names(sampRename)=manifest$MapID

ds=ds %>%
    gather(Sample,Count,-Status) %>%
    mutate(Sample=fixSampleNames(Sample)) %>%
    mutate(Status=gsub("_.*","",Status)) %>%
    group_by(Sample,Status) %>%
    summarize(Counts=sum(Count)) %>%
    mutate(Status=ifelse(Status=="Assigned","InPeaks","Outside"))

pg0=ggplot(ds,aes(Sample,Counts,fill=Status)) +
    theme_light(base_size=16) +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(guide=guide_axis(n.dodge=2))

pg1=pg0 + ggtitle("Mapped Reads in MACS Peaks") +
    geom_bar(stat="identity") +
    scale_y_continuous(labels = label_number(suffix = " M", scale = 1e-6))

pg2=pg0 + geom_bar(stat="identity",position="fill") +
    ggtitle("Fraction Reads in MACS Peaks") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ylab("Percentage")

dd=read_tsv("peaks_raw_fcCounts.txt",comment="#")

peak.annote=dd %>% select(PeakNo=Geneid,Chr,Start,End,Strand,Length)

#
# Remove excluded points
#

manifest=manifest %>% filter(!grepl("^EXC",Group))

d=dd %>%
    select(PeakNo=Geneid,matches("Proj.*_s_")) %>%
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
dp=pr$rotation %>% data.frame %>% rownames_to_column("SampleID") %>% left_join(manifest)

pp1=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=16) + geom_point(size=4) + scale_color_brewer(palette="Paired")
pp2=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=16) + geom_point(size=2) + scale_color_brewer(palette="Paired") + geom_label_repel(color="black")

design <- model.matrix(~0+group)
y <- estimateDisp(y,design)

# fit <- glmQLFit(y,design)
# qlf <- glmQLFTest(fit,coef=2)
# topTags(qlf)

txdb=TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene

doQLFStats<-function(y,design,contrast,fdrCut=0.05) {

    #fit <- glmQLFit(y,design)
    #qlf <- glmQLFTest(fit,contrast=contrast)
    #model=qlf

    fit <- glmFit(y,design)
    lrt <- glmLRT(fit,contrast=contrast)
    model=lrt

    cp=model$comparison
    compTag=paste0(sort(strsplit(cp," ")[[1]],decreasing=T),collapse="") %>%
        gsub("-1\\*","-",.) %>%
        gsub("^1\\*","",.) %>%
        gsub("group","",.)

    tt=model$table %>%
        data.frame %>%
        rownames_to_column("PeakNo") %>%
        tibble %>%
        mutate(FDR=p.adjust(PValue)) %>%
        arrange(FDR)

    max.logFC=max(abs(tt$logFC))

    p.ma=ggplot(arrange(tt,desc(PValue)),aes(logCPM,logFC,color=FDR<fdrCut)) +
        theme_light(base_size=16) +
        geom_point() +
        scale_color_manual(values=c("#7f7f7f33","#e31a1c")) +
        ggtitle(compTag) +
        scale_y_continuous(limits=c(-1,1)*max.logFC)

    p.vc=ggplot(arrange(tt,desc(PValue)),aes(logFC,PValue,color=FDR<fdrCut)) +
        theme_light(base_size=16) +
        geom_point() +
        scale_y_continuous(trans=reverselog_trans(10)) +
        scale_x_continuous(limits=c(-1,1)*max.logFC) +
        scale_color_manual(values=c("#7f7f7f33","#e31a1c")) +
        ggtitle(compTag)

    ans=tt %>% filter(FDR<fdrCut) %>% left_join(peak.annote) %>% select(-matches("^F$|^LR")) %>% filter(Chr!="MT")

    sig.peaks=ans %>%
        select(V1=Chr,V2=Start,V3=End,V4=PeakNo,V5=PValue) %>%
        mutate(V5=-10*log10(V5),V1=paste0("chr",V1)) %>%
        data.frame

    sig.peaks=ChIPseeker:::peakDF2GRanges(sig.peaks)

    aa=annotatePeak(sig.peaks,TxDb=txdb,tssRegion=c(-5000,5000),annoDb="org.Hs.eg.db")

    peak.annote=as.data.frame(aa) %>% tibble %>% dplyr::select(PeakNo=V4,SYMBOL,GENENAME,annotation,distanceToTSS)

    tbl=ans %>% left_join(peak.annote)

    list(tbl=tbl,comparison=compTag,p.ma=p.ma,p.vc=p.vc,model=model)

}

pp=strsplit(getwd(),"/")[[1]]
projNo=grep("^Proj_|^B-\\d+",pp,value=T)

pfile=cc(projNo,RUNTAG,"ATACSeqQC.pdf")
pdf(pfile,width=11,height=8.5)
print(pg1)
print(pg2)
print(pp1)
print(pp2)
dev.off()

cat("\n\nNOT IMPLEMENTED: Set contrasts\n\n")
quit()

# Need to write code to allow contrasts to be specified from command line

contrast=c(-1,1,0,0)

res=list()
res[[1]]=doQLFStats(y,design,c(-1,1,0,0))
res[[2]]=doQLFStats(y,design,c(0,0,-1,1))

tbls=map(res,"tbl")
names(tbls)=map(res,"comparison") %>% unlist
stats=map(tbls,nrow) %>% bind_rows %>% gather(Comparison,NumSig)

write.xlsx(c(list(Summary=stats),tbls),cc(projNo,RUNTAG,"DiffPeaksEdgeR.xlsx"))

pfile=cc(projNo,RUNTAG,"DiffPeaks_%02d.png")
pngCairo(pfile,width=11,height=8.5)
print(pg1)
print(pg2)
print(pp1)
print(pp2)
for(ii in seq(len(res))) {
    print(res[[ii]]$p.ma)
    print(res[[ii]]$p.vc)
}
dev.off()
mergePNGs(pfile)


