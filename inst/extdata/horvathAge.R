library(rtracklayer)

# Horvath coefficients 
maps <- read.csv("HorvathAdditionalFile3.csv", row.names=1) 
CpGmaps <- subset(maps, grepl("^cg", rownames(maps)))[, c(1,2)]
CpGs <- rownames(CpGmaps) 

# reminder: the epigenetic clock was designed against HumanMethylation27 probes 
horvath <- list() 
for (i in c("hg19", "hg38")) {
  manifest <- paste("hm27", i, "manifest", "rds", sep=".")
  mappings <- granges(subset(readRDS(manifest), names %in% CpGs))
  bsgenome <- paste0("BSgenome.Hsapiens.UCSC.", i)
  library(bsgenome, character.only=TRUE) 
  genome(mappings) <- i
  seqinfo(mappings) <- seqinfo(Hsapiens)[seqlevels(mappings)] 
  mappings$name <- names(mappings)
  mappings$score <- maps[names(mappings), 1]
  mappings$shrunken <- maps[names(mappings), 2]
  mappings$ENSR <- NA
  mappings$ENSRstart <- NA 
  mappings$ENSRend <- NA 
  mappings$ENSRtype <- NA
  mappings$sign <- as.factor(ifelse(mappings$score > 0, "+", "-")) 
  horvath[[i]] <- unstrand(trim(mappings)) 
}

# since ENSEMBL uses GRCh37/GRCh38...
harmonizeAndOverlap <- function(x, y) {
  seqlevelsStyle(x) <- seqlevelsStyle(y)
  res <- unstrand(suppressWarnings(trim(sort(subsetByOverlaps(x, y)))))
  names(res) <- res$ID
  return(res)
}

# ENSEMBL build 93
regBuildFiles <- c(
hg19="homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20161117.gff.gz",
hg38="homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20161111.gff.gz"
)

# regulatory build from ENSEMBL 
regBuilds <- list(
  hg19=harmonizeAndOverlap(import(regBuildFiles["hg19"]), horvath$hg19),
  hg38=harmonizeAndOverlap(import(regBuildFiles["hg38"]), horvath$hg38)
)

# find overlaps with regulatory build regions where possible
for (i in c("hg19", "hg38")) {
  ha <- horvath[[i]]
  rb <- regBuilds[[i]]
  ol <- suppressWarnings(findOverlaps(trim(ha), trim(rb)))
  mcols(ha)[queryHits(ol), "ENSR"] <- mcols(rb)[subjectHits(ol), "ID"]
  mcols(ha)[queryHits(ol), "ENSRstart"] <- 
    as.integer(mcols(rb)[subjectHits(ol), "bound_start"])
  mcols(ha)[queryHits(ol), "ENSRend"] <- 
    as.integer(mcols(rb)[subjectHits(ol), "bound_end"])
  mcols(ha)[queryHits(ol), "ENSRtype"] <- 
    mcols(rb)[subjectHits(ol), "feature_type"]
  mcols(ha)$ENSRtype <- 
    as.factor(sub("TF binding site", "TFBS", 
                  sub("CTCF Binding Site", "CTCF", 
                      sub("Promoter Flanking Region", "Flanking", 
                          sub("Open chromatin", "Accessible", 
                              mcols(ha)$ENSRtype)))))
  horvath[[i]] <- ha
}

# shrunken version:
horvathAgeShrunken <- lapply(horvath, function(x) { 
  x <- subset(x, !is.na(x$shrunken))
  x$score <- x$shrunken
  x$shrunken <- NULL 
  return(x)
})
horvathAgeShrunken$intercept <- maps[1,2]
save(horvathAgeShrunken, file="horvathAgeShrunken.rda")

# unshrunken version:
horvathAge <- lapply(horvath, function(x) { 
  x$shrunken <- NULL
  return(x) 
})
horvathAge$intercept <- maps[1,1]
save(horvathAge, file="horvathAge.rda")

# clean up 
rm(horvath)
