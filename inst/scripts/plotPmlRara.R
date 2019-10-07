### Example scritp for plotting PML RARA circle plot
### Written by: Tim Triche, Jr.
### Last edit: October 7, 2019
library(Homo.sapiens)
library(GenomicRanges)
library(StructuralVariantAnnotation)

file <- system.file("extdata", "PML_RARA_TCGA_AB_2998.hg19.tsv",
                    package="biscuiteer")
gr <- as(read.table(file, header=TRUE), "GRanges")
seqinfo(gr) <- seqinfo(Homo.sapiens)[seqlevels(gr)]

paired <- .makePairs(split(gr, substr(gr$name, nchar(gr$name), nchar(gr$name))))
mcols(paired)$name <- sapply(paired,
                             function(x) substr(mcols(first(x))$name, 1, 6))
mcols(paired)$score <- sapply(paired, function(x) mcols(first(x))$QUAL)

library(circlize)
circos.initializeWithIdeogram()
circos.genomicLink(as.data.frame(S4Vectors::first(paired)),
                   as.data.frame(S4Vectors::second(paired)),
                   col="red", lwd=3, border="darkred")

bp <- subset(genes(Homo.sapiens, columns="SYMBOL"), SYMBOL %in% c("PML",
                                                                  "RARA"))
names(bp) <- sapply(bp$SYMBOL, `[[`, 1)
circos.genomicLabels(bp, labels=names(bp), side="inside")
dev.copy2pdf(file="PML_RARA.pdf")


# Helper function
.makePairs <- function(...) {
  Pairs(as.list(...)[[1]], as.list(...)[[2]])
}
