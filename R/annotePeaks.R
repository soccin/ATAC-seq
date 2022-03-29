source("analyzeATAC.R")
contrast
    fit <- glmQLFit(y,design)
    qlf <- glmQLFTest(fit,contrast=contrast)
    tt=qlf$table %>%
        data.frame %>%
        rownames_to_column("PeakNo") %>%
        tibble %>%
        mutate(FDR=p.adjust(PValue)) %>%
        arrange(FDR)
    p.ma=ggplot(arrange(tt,desc(PValue)),aes(logCPM,logFC,color=FDR<fdrCut)) +
        theme_light(base_size=16) +
        geom_point() +
        scale_color_manual(values=c("#7f7f7f33","#e31a1c")) +
        ggtitle(qlf$comparison)
    p.vc=ggplot(arrange(tt,desc(PValue)),aes(logFC,PValue,color=FDR<fdrCut)) +
        theme_light(base_size=16) +
        geom_point() +
        scale_y_continuous(trans=reverselog_trans(10)) +
        scale_x_continuous(limits=c(-1,1)*8) +
        scale_color_manual(values=c("#7f7f7f33","#e31a1c")) +
        ggtitle(qlf$comparison)
p.vc
fdrCut=0.05
p.vc
tt
max(abs(tt$logFC))
    fit <- glmQLFit(y,design)
    qlf <- glmQLFTest(fit,contrast=contrast)
    tt=qlf$table %>%
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
        ggtitle(qlf$comparison) +
        scale_y_continuous(limits=c(-1,1)*max.logFC)
    p.vc=ggplot(arrange(tt,desc(PValue)),aes(logFC,PValue,color=FDR<fdrCut)) +
        theme_light(base_size=16) +
        geom_point() +
        scale_y_continuous(trans=reverselog_trans(10)) +
        scale_x_continuous(limits=c(-1,1)*max.logFC) +
        scale_color_manual(values=c("#7f7f7f33","#e31a1c")) +
        ggtitle(qlf$comparison)
p.ma
tt
tt %>% filter(FDR<fdrCut)
tt %>% filter(FDR<fdrCut) %>% left_join(peak.annote)
ans=tt %>% filter(FDR<fdrCut) %>% left_join(peak.annote) %>% select(-F)
ans
library(ChipSeeker)
library(ChIPseeker)
ChIPseeker::readPeakFile
files <- getSampleFiles()
ChIPseeker::peak2DF
ChIPseeker:::peak2DF
ChIPseeker:::peak2DF(files[[4]])
ChIPseeker:::peak2DF(files[[4]]) %>% tibble
ans
ans %>% select(V1=Chr,V2=Start,V3=End,V4=PeakNo,V5=-10*log10(PValue))
ans %>% select(V1=Chr,V2=Start,V3=End,V4=PeakNo) %>% mutate(V5=-10*log10(PValue))
ans %>% select(V1=Chr,V2=Start,V3=End,V4=PeakNo) %>% mutate(V5=-10*log10(PValue))
ans
ans %>% select(V1=Chr,V2=Start,V3=End,V4=PeakNo) 
ans %>% select(V1=Chr,V2=Start,V3=End,V4=PeakNo,PValue)  %>% mutate(V5=FDR)
ans %>% select(V1=Chr,V2=Start,V3=End,V4=PeakNo,V5=PValue) %>% mutate(V5=-10*log10(V5))
ans %>% select(V1=Chr,V2=Start,V3=End,V4=PeakNo,V5=PValue) %>% mutate(V5=-10*log10(V5),V1=paste0("chr",V1))
sig.peaks=ans %>% select(V1=Chr,V2=Start,V3=End,V4=PeakNo,V5=PValue) %>% mutate(V5=-10*log10(V5),V1=paste0("chr",V1))
readPeakFile
ChIPseeker:::peakDF2GRanges(sig.peaks)
ChIPseeker:::peak2DF(files[[4]])
pdf=ChIPseeker:::peak2DF(files[[4]])
head(pdf)
head(sig.peaks)
ChIPseeker:::peakDF2GRanges(data.frame(sig.peaks))
sig.peaks=ans %>% select(V1=Chr,V2=Start,V3=End,V4=PeakNo,V5=PValue) %>% mutate(V5=-10*log10(V5),V1=paste0("chr",V1)) %>% data.frame
ChIPseeker:::peakDF2GRanges(sig.peaks)
sig.peaks=ChIPseeker:::peakDF2GRanges(sig.peaks)
sig.peaks
annotePeaks
annotatePeak
annotatePeak(sig.peaks)
sig.peaks
ans
ans %>% count(Chr)
ans %>% count(Chr) %>% data.frame
    ans=tt %>% filter(FDR<fdrCut) %>% left_join(peak.annote) %>% select(-F) %>% filter(Chr!="MT")
    sig.peaks=ans %>%
        select(V1=Chr,V2=Start,V3=End,V4=PeakNo,V5=PValue) %>%
        mutate(V5=-10*log10(V5),V1=paste0("chr",V1)) %>%
        data.frame
    sig.peaks=ChIPseeker:::peakDF2GRanges(sig.peaks)
annotatePeak(sig.peaks)
aa=annotatePeak(sig.peaks)
head(aa)
aa
plotAnnoPie(aa)
plotAnnoPar(aa)
plotAnnoBar(aa)
vennpie(aa)
upsetplot(aa)
upsetplot(aa,vennpie=T)
install.packages("ggimage")
plotDistToTSS(aa)

install.packages("ReactomePA")
library(ReactomePA)
install.packages("ReactomePA")
BiocManager::install("ReactomePA")
as.data.frame(aa)
as.data.frame(aa) %>% tibble
as.data.frame(aa)$geneId
ReactomePA
library(ReactomePA)
ReactomePA::enrichPathway(as.data.frame(aa)$geneId)
BiocManager::install("org.Hs.eg.db")
ReactomePA::enrichPathway(as.data.frame(aa)$geneId)
pathway1=ReactomePA::enrichPathway(as.data.frame(aa)$geneId)
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
seqlevelsStyle(edb) <- "UCSC"
BiocManager::install("EnsDb.Hsapiens.v75")
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
seqlevelsStyle(edb) <- "UCSC"
aa=annotatePeak(sig.peaks,tssRegion=c(-5000,5000),TxDb=edb,annoDb="org.Hs.eg.db")
aa
sig.peaks
aa
as.data.frame(aa)
as.data.frame(aa) %>% tibble
as.data.frame(aa) %>% tibble %>% select(GENENAME)
as.data.frame(aa) %>% as_tibble %>% select(GENENAME)
as.data.frame(aa) %>% as_tibble
aa
data.frame(aa)
tibble(data.frame(aa))
tibble(data.frame(aa)) %>% select(GENENAME)
tibble(data.frame(aa)) %>% dplyr::select(GENENAME)
tibble(data.frame(aa)) %>% dplyr::select(SYMBOL,GENENAME)
tibble(data.frame(aa))
tibble(data.frame(aa)) %>% dplyr::select(SYMBOL,GENENAME,annotation)
tibble(data.frame(aa)) %>% dplyr::select(SYMBOL,GENENAME,annotation,distanceToTSS) 
tibble(data.frame(aa)) %>% dplyr::select(SYMBOL,GENENAME,annotation,distanceToTSS)  %>% data.frame
plotAnnoPie(aa)
tibble(data.frame(aa)) %>% dplyr::select(SYMBOL,GENENAME,annotation,distanceToTSS)  %>% data.frame
tibble(data.frame(aa)) %>% dplyr::select(SYMBOL,GENENAME,annotation,distanceToTSS)
tibble(data.frame(aa)) %>% dplyr::select(PeakNo=V5,SYMBOL,GENENAME,annotation,distanceToTSS)
tibble(data.frame(aa)) %>% dplyr::select(PeakNo=V4,SYMBOL,GENENAME,annotation,distanceToTSS)
aa=annotatePeak(sig.peaks,tssRegion=c(-5000,5000),TxDb=edb,annoDb="org.Hs.eg.db")
