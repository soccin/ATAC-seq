read_insdat<-function(ff) {
    rr=readLines(ff)
    start=grep("^insert_size",rr)
    read_tsv(paste0(rr[start:len(rr)],collapse="\n")) %>%
        gather(Sample,Count,-insert_size) %>%
        mutate(Sample=gsub("\\.fr_count","",Sample))
}

suppressPackageStartupMessages(require(tidyverse))


insFiles=dir(pattern="___INS.txt")
dd=map(insFiles,read_insdat) %>% bind_rows %>% group_by(Sample) %>% mutate(Density=Count/sum(Count))

pg=ggplot(dd,aes(insert_size,Count,color=Sample)) +
    theme_light(base_size=16) +
    geom_line() +
    scale_x_continuous(breaks=c(0,rep(1:10)*100)) +
    scale_color_brewer(palette="Paired")

pg2=ggplot(dd,aes(insert_size,Count,color=Sample)) +
    theme_light(base_size=16) +
    geom_line() +
    scale_color_brewer(palette="Paired") +
    scale_x_continuous(breaks=c(0,rep(1:10)*100))


pdf(file="postInsDistribution.pdf",width=11,height=8.5)
print(pg)
print(pg2)
dev.off()
