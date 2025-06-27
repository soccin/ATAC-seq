suppressPackageStartupMessages({
    require(tidyverse)
    require(readxl)
})

args=commandArgs(trailing=T)
if(len(args)!=2) {
    cat("\n\n\tgetHomerBED.R DIFF_PEAKS_FILE SHEET_NO\n\n")
    quit()
}

diffPeaksFile=args[1]
sheetNo=as.numeric(args[2])

compName=excel_sheets(diffPeaksFile)[sheetNo]
dp=read_xlsx(diffPeaksFile,sheet=sheetNo) %>% mutate(Chr=paste0("chr",Chr))

homerFile=cc(gsub(".xlsx","",diffPeaksFile),"_",compName,".tsv")

dp %>%
    select(Chr,Start,End,PeakNo,SYMBOL,Strand) %>%
    write_tsv(homerFile,col_names=F)

cat(homerFile,"\n")

#dp %>% select(Chr,Start,End,PeakNo,distanceToTSS,Strand) %>% filter(distanceToTSS==0) %>% write_tsv("B-101-131_Pass1_DiffPeaksEdgeR___DmRp_Rm___TSSeq0.tsv",col_names=F)

