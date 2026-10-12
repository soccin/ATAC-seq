cat("## START: plotINSStats.R\n")

read_insdat<-function(ff) {
    rr=readLines(ff)
    start=grep("^insert_size",rr)
    read_tsv(I(paste0(rr[start:len(rr)],collapse="\n")),show_col_types=FALSE) %>%
        gather(Sample,Count,-insert_size) %>%
        mutate(Sample=gsub("\\.fr_count","",Sample)) %>%
        rename(MapID=Sample)
}

suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(fs))
suppressPackageStartupMessages(require(ggsci))
suppressPackageStartupMessages(require(ggforce))


#
# 34 group colors; scale_color_uchicago() alone has 9 and draws any
# further group in grey
#
groupColors=c(pal_uchicago("default")(9),pals::cols25())

insFiles=dir_ls("out",recurs=T,regex="___INS.txt")
dd=map(insFiles,read_insdat) %>% bind_rows %>% group_by(MapID) %>% mutate(Density=Count/sum(Count))

manifest=read_csv("sampleManifest.csv") %>% group_by(Group) %>% mutate(Rep=as.factor(row_number())) %>% ungroup
dd=left_join(dd,manifest)

pg1=ggplot(dd,aes(insert_size,1000*Density,color=Group,group=SampleID)) +
    theme_light(base_size=16) +
    geom_line(aes(linetype=Rep)) +
    scale_color_manual(values=groupColors) +
    guides(color=guide_legend(override.aes=list(size=2.4))) +
    scale_x_continuous(breaks=c(0,rep(1:10)*100),limits=c(0,1000))

minDensity=10^ceiling(log10(max(dd$Density)/1e3))

pg2=ggplot(dd,aes(insert_size,Density,color=Group,group=SampleID)) +
    theme_light(base_size=16) +
    geom_line(aes(linetype=Rep)) +
    scale_color_manual(values=groupColors) +
    guides(color=guide_legend(override.aes=list(size=2.4))) +
    scale_x_continuous(breaks=c(0,rep(1:10)*100),limits=c(0,1000)) +
    scale_y_log10(limits=c(minDensity,NA))

#
# Project for output naming: the last directory of the working directory
# path named Proj_* or B-<digits>, else ATAC
#
pp=strsplit(getwd(),"/")[[1]]
projNo=tail(grep("^Proj_|^B-\\d+",pp,value=T),1)
if(len(projNo)==0) {
    projNo="ATAC"
}

pp2=ggplot(dd,aes(insert_size,1000*Density,color=Group,group=SampleID)) + theme_light(base_size=12) + geom_line() + scale_color_manual(values=groupColors) + guides(color=guide_legend(override.aes=list(size=2.4))) + scale_x_continuous(breaks=c(0,rep(1:5)*200),limits=c(0,1000))
pw=pp2 + facet_wrap_paginate(~SampleID,nrow=2,ncol=3) + theme(legend.position="none")
nPages=n_pages(pw)

pdf(file=cc(projNo,"_postInsDistribution.pdf"),width=11,height=8.5)
print(pg1)
print(pg2)
for(ii in seq(nPages)) {
    print(pp2 + facet_wrap_paginate(~SampleID,nrow=2,ncol=3,page=ii))
}
dev.off()

cat("## END: plotINSStats.R\n")

