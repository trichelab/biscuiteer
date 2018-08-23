library(rtracklayer) 

# hg19
CGI.hg19 <- import("HMM_CpG_islands.hg19.bed.gz", genome="hg19")
resorts.hg19 <- trim(resize(CGI.hg19, width(CGI.hg19) + 8000, fix="center"))
openSeas.hg19 <- subset(gaps(resorts.hg19), strand == "*")
save(openSeas.hg19, file="../../data/openSeas.hg19.rda")

# hg38
CGI.hg38 <- import("HMM_CpG_islands.hg38.bed.gz", genome="hg38")
resorts.hg38 <- trim(resize(CGI.hg38, width(CGI.hg38) + 8000, fix="center"))
openSeas.hg38 <- subset(gaps(resorts.hg38), strand == "*")
save(openSeas.hg38, file="../../data/openSeas.hg38.rda")
