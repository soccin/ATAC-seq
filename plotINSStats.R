read_insdat<-function(ff) {
    rr=readLines(ff)
    start=grep("^insert_size",rr)
    read_tsv(paste0(rr[start:len(rr)],collapse="\n")) %>%
        gather(Sample,Count,-insert_size) %>%
        mutate(Sample=gsub("\\.fr_count","",Sample)) %>%
        rename(MapID=Sample)
}

suppressPackageStartupMessages(require(tidyverse))


insFiles=dir("out/metrics",pattern="___INS.txt",full.names=T)
dd=map(insFiles,read_insdat) %>% bind_rows %>% group_by(MapID) %>% mutate(Density=Count/sum(Count))

manifest=read_csv("sampleManifest.csv") %>% group_by(Group) %>% mutate(Rep=as.factor(row_number())) %>% ungroup
dd=left_join(dd,manifest)

pg1=ggplot(dd,aes(insert_size,1000*Density,color=Group,group=SampleID)) +
    theme_light(base_size=16) +
    geom_line(aes(linetype=Rep)) +
    scale_color_brewer(palette="Set1") +
    scale_x_continuous(breaks=c(0,rep(1:10)*100))

minDensity=10^ceiling(log10(max(dd$Density)/1e3))

pg2=ggplot(dd,aes(insert_size,Density,color=Group,group=SampleID)) +
    theme_light(base_size=16) +
    geom_line(aes(linetype=Rep)) +
    scale_color_brewer(palette="Set1") +
    scale_x_continuous(breaks=c(0,rep(1:10)*100)) +
    scale_y_log10(limits=c(minDensity,NA))


pdf(file="postInsDistribution.pdf",width=11,height=8.5)
print(pg1)
print(pg2)
dev.off()
