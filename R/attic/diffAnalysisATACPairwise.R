args=commandArgs(trailing=T)
if(len(args)<3) {
    cat("
    Usage: analyzeATAC.R GENOME SampleManifest.csv Comparisons.csv [RUNTAG]

        Comparisons.csv (NO Column names)

            aNSC_loxp15,aNSC_p53

            Sign convention X2-X1; e.g., aNSC_p53-aNSC_loxp15


")
    quit()
}

#
# Example inputs:
#   SampleManifest.csv
#
#        SampleID,Group,MapID
#        aNSC_loxp15_1,aNSC_loxp15,s_aNSC_loxp15_1
#        aNSC_loxp15_2,aNSC_loxp15,s_aNSC_loxp15_2
#        aNSC_p53_1,aNSC_p53,s_aNSC_p53_1
#
#   Comparisons.csv
#
#        aNSC_loxp15,aNSC_p53
#
#   Sign convention X2-X1; e.g., aNSC_p53-aNSC_loxp15
#

fixSampleNames<-function(ss) {
    if(grepl("___MD",ss[1])) {
        fixSampleNamesPEmap(ss)
    } else {
        fixSampleNamesBIC(ss)
    }
}

fixSampleNamesBIC<-function(ss) {
    sampRename[gsub("_postProcess.*","",ss) %>% gsub(".*_s_","s_",.)] %>%
        unname
}

fixSampleNamesPEmap<-function(ss) {
    sampRename[basename(ss) %>% gsub("___.*","",.)] %>% unname
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
require(ggsci)

require(tidyverse)
require(fs)
require(openxlsx)

require(ggrastr)

ds=read_tsv("peaks_raw_fcCounts.txt.summary")

GENOME=args[1]
MANIFEST_FILE=args[2]
COMPARISON_FILE=args[3]
if(len(args)==4) {
    RUNTAG=args[4]
} else {
    RUNTAG=""
}

if(GENOME=="human") {
    txdb=TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    annoDb="org.Hs.eg.db"
} else if(GENOME=="mouse") {
    txdb=TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
    annoDb="org.Mm.eg.db"
} else {
    cat("\n\tUnknown GENOME",GENOME,"\n")
    cat("\tValid genomes: human, mouse\n\n")
    quit()
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

dd=read_tsv("peaks_raw_fcCounts.txt",comment="#")

peak.annote=dd %>% select(PeakNo=Geneid,Chr,Start,End,Strand,Length)

#
# Remove excluded points
#

manifest=manifest %>% filter(!grepl("^EXC",Group))

d=dd %>%
    select(PeakNo=Geneid,matches(".bam$")) %>%
    data.frame(check.names=F) %>%
    column_to_rownames("PeakNo")

colnames(d)=fixSampleNames(colnames(d))

if(!all(colnames(d) %in% manifest$MapID)) {
    cat("\nERROR in creation of count matrix 'd'\n")
    cat("LINE-141\n\n")
    rlang::abort("ERROR")
}

d=d[,manifest$SampleID]
group=factor(manifest$Group)
y <- DGEList(counts=d,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

pr=prcomp(cpm(y,log=T),scale=F)
dp=pr$rotation %>% data.frame %>% rownames_to_column("SampleID") %>% left_join(manifest)

colors1=c(pal_uchicago("default")(9),pals::cols25())

pp1=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=16) + geom_point(size=4,alpha=.6) + scale_color_manual(values=colors1)
pp2=ggplot(dp,aes(PC1,PC2,color=Group,label=SampleID)) + theme_light(base_size=16) + geom_point(size=2) + scale_color_manual(values=colors1) +
    geom_label_repel(
        color="black",
        max.overlaps=Inf,
        min.segment.length = 0,
        size=3,
        force=16,
        max.time=10,
        max.iter=100000)

design <- model.matrix(~0+group)
y <- estimateDisp(y,design)

# fit <- glmQLFit(y,design)
# qlf <- glmQLFTest(fit,coef=2)
# topTags(qlf)

doQLFStats<-function(y,design,contrast,fdrCut=0.05) {

    #fit <- glmQLFit(y,design)
    #qlf <- glmQLFTest(fit,contrast=contrast)
    #model=qlf

    cat("Num peaks =",nrow(y),"\n")

    cat("\nRefilter peaks to just those in comp\n\n")
    #browser()

    compSamps=which(rowSums(design[,contrast!=0])>0) %>% unname

    gComp=factor(group[compSamps])
    yComp=getCounts(y)[,compSamps]
    designComp <- model.matrix(~0+gComp)

    yComp <- DGEList(counts=yComp,group=gComp)
    keep <- filterByExpr(yComp,design=designComp)
    yComp <- yComp[keep,,keep.lib.sizes=FALSE]
    yComp <- calcNormFactors(yComp)

    yComp <- estimateDisp(yComp,designComp)

    contrastComp=contrast[contrast!=0]

    fit <- glmFit(yComp,designComp)
    lrt <- glmLRT(fit,contrast=contrastComp)
    model=lrt

    cp=model$comparison
    compTag=paste0(sort(strsplit(cp," ")[[1]],decreasing=T),collapse="") %>%
        gsub("-1\\*","-",.) %>%
        gsub("^1\\*","",.) %>%
        gsub("gComp","",.)

    tt=model$table %>%
        data.frame %>%
        rownames_to_column("PeakNo") %>%
        tibble %>%
        mutate(FDR=p.adjust(PValue)) %>%
        arrange(FDR,PValue) %>%
        mutate(PValue.mod=ifelse(PValue<.Machine$double.eps^2,.Machine$double.eps^2,PValue))

    require(IHW)
    #ihwRes=ihw(PValue ~ logCPM,dat=tt,alpha=.25)
    ihwRes=ihw(PValue ~ Length,dat=tt%>%left_join(peak.annote),alpha=.25)
    #browser()
    cat("ihw rejections =",rejections(ihwRes),"\n")

    #
    # PValue.mod is a clipped PValue to make the volcano plot look reasonable
    #
    max.logFC=max(abs(tt$logFC))

    p.ma=ggplot(arrange(tt,desc(PValue)),aes(logCPM,logFC,color=FDR<fdrCut)) +
        theme_light(base_size=16) +
        rasterize(geom_point()) +
        scale_color_manual(values=c("#7f7f7f33","#e31a1c")) +
        ggtitle(compTag) +
        scale_y_continuous(limits=c(-1,1)*max.logFC)

    p.vc=ggplot(arrange(tt,desc(PValue)),aes(logFC,PValue.mod,color=FDR<fdrCut)) +
        theme_light(base_size=16) +
        rasterize(geom_point()) +
        scale_y_continuous(trans=reverselog_trans(10)) +
        scale_x_continuous(limits=c(-1,1)*max.logFC) +
        scale_color_manual(values=c("#7f7f7f33","#e31a1c")) +
        ggtitle(compTag)

    ans=tt %>% filter(FDR<fdrCut) %>% left_join(peak.annote) %>% select(-matches("^F$|^LR")) %>% filter(Chr!="MT")

    sig.peaks=ans %>%
        select(V1=Chr,V2=Start,V3=End,V4=PeakNo,V5=PValue) %>%
        mutate(V5=-10*log10(V5))


    if(nrow(ans)>0) {
        if(substr(ans$Chr[1],1,3)!="chr") {
            sig.peaks=sig.peaks %>% mutate(V1=paste0("chr",V1)) %>% data.frame
        } else {
            sig.peaks=sig.peaks %>% data.frame
        }

        sig.peaks=ChIPseeker:::peakDF2GRanges(sig.peaks)

        aa=annotatePeak(sig.peaks,TxDb=txdb,tssRegion=c(-5000,5000),annoDb=annoDb)

        peak.annote=as.data.frame(aa) %>% tibble %>% dplyr::select(PeakNo=V4,SYMBOL,GENENAME,annotation,distanceToTSS)

        tbl=ans %>% left_join(peak.annote)
    } else {
        cat("\n   No significant peaks",compTag,"at FDR",fdrCut,"\n\n")
        tbl=ans
    }

    list(tbl=tbl,comparison=compTag,p.ma=p.ma,p.vc=p.vc,model=model)

}

cwd=strsplit(getwd(),"/")[[1]]
projNo=grep("^Proj_|^B-\\d+",cwd,value=T)

gNames=gsub("group","",colnames(design))
comps=read_csv(COMPARISON_FILE,col_names=F)

res=list()

for(ci in transpose(comps)) {
    contrast=as.numeric(gNames==ci$X2)-as.numeric(gNames==ci$X1)
    cat("========================================================\n")
    cat(str(ci))
    if(!(sum(contrast!=0)>0 & sum(contrast)==0)) {
        cat("\n\n")
        cat("   Invalid contrasts")
        cat("\n\n")
        stop("FATAL")
    }
    res[[len(res)+1]]=doQLFStats(y,design,contrast)
}
cat("========================================================\n")
cat("========================================================\n\n")

tbls=map(res,"tbl")
names(tbls)=map(res,"comparison") %>% unlist %>% substr(.,1,31)
stats=map(tbls,nrow) %>% bind_rows %>% gather(Comparison,NumSig)

write.xlsx(c(list(Summary=stats),tbls),cc(projNo,RUNTAG,"DiffPeaksEdgeRv2.xlsx"))

# pfile=cc(projNo,RUNTAG,"DiffPeaksV2_%02d.png")
# pngCairo(pfile,width=11,height=8.5)

pfile=cc(projNo,RUNTAG,"DiffPeaksV2.pdf")
pdf(pfile,width=11,height=8.5)

print(pp1)
print(pp2)
for(ii in seq(len(res))) {
    print(res[[ii]]$p.ma)
    print(res[[ii]]$p.vc)
}
dev.off()

# mergePNGs(pfile)


