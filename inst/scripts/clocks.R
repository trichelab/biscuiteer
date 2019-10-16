library(rtracklayer)

# old and new coefficients 
horvath <- read.csv("clocks_input/HorvathClock.csv", row.names=1)[, 1:2]
coefs <- list(
  horvath=subset(horvath, !is.na(CoefficientTraining))[, 1, drop=FALSE],
  horvathshrunk=subset(horvath, !is.na(CoefficientTrainingShrunk))[, 2, FALSE],
  skinandblood=read.csv("clocks_input/SkinAndBloodClock.csv", row.names=1)[, 1, drop=FALSE],
  hannum=read.csv("clocks_input/HannumClock.csv", row.names=1)[, 5, drop=FALSE]
)

# merge into one big data.frame of coefficients
reRowName <- function(x) { rownames(x) <- x[,1]; return(x[,-1]) }
for (i in names(coefs)) colnames(coefs[[i]]) <- i
clocks <- reRowName(
  merge(with(coefs, merge(horvath,horvathshrunk,by="row.names",all=TRUE)),
        with(coefs, merge(hannum, skinandblood, by="row.names",all=TRUE)),
        by="Row.names", all=TRUE)
)
clocks["(Intercept)", "hannum"] <- 0 # none provided in Molecular Cell paper

# regulatory builds from ENSEMBL; processed previously, here for reference 
regBuildFiles <- c(
hg19="homo_sapiens.GRCh37.Regulatory_Build.regulatory_features.20161117.gff.gz",
hg38="homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20161111.gff.gz"
)

# since ENSEMBL uses GRCh37/GRCh38...
harmonizedSubset <- function(ENSGFF, UCSCGR) {
  seqlevelsStyle(UCSCGR) <- "NCBI"
  res <- unstrand(suppressWarnings(trim(sort(import(ENSGFF, which=UCSCGR)))))
  seqlevelsStyle(res) <- seqlevelsStyle(UCSCGR) <- "UCSC"
  res <- keepSeqlevels(res, seqlevelsInUse(UCSCGR), pruning.mode="coarse")
  seqinfo(res) <- seqinfo(UCSCGR)[seqlevels(res)] 
  genome(res) <- unique(genome(UCSCGR)) 
  names(res) <- res$ID
  return(sort(trim(res)))
}

# all manifests (remapped) are from http://zwdzwd.github.io/InfiniumAnnotation
reload <- FALSE 
CpGs <- rownames(clocks)[-1]
chips <- c(hm27="hm27", hm450="hm450")
assemblies <- c(hg19="hg19", hg38="hg38") 
getMfest <- function(chip, a) {
  res <- unstrand(granges(readRDS(paste(chip, a,"manifest","rds", sep="."))))
  res <- subset(sort(res[intersect(CpGs, names(res))]), seqnames != "*") 
  seqinf <- paste0("seqinfo.", a)
  if (!exists(seqinf)) data(list=seqinf, character.only=TRUE)
  seqinfo(res) <- get(seqinf)[seqlevels(res)]
  genome(res) <- a
  return(keepSeqlevels(res, seqlevelsInUse(res), pruning.mode="coarse"))
}
if (reload) {
  mfests <- lapply(assemblies, function(a) lapply(chips, getMfest, a=a))
  regBuilds <- list(
    hg19=harmonizedSubset(regBuildFiles["hg19"], Reduce(union, mfests$hg19)),
    hg38=harmonizedSubset(regBuildFiles["hg38"], Reduce(union, mfests$hg38))
  )
} else { 
  mfests <- readRDS("clocks_input/mfests.rds") # subsetted (much much smaller)
  regBuilds <- readRDS("clocks_input/regBuilds.rds") 
}

# hg19 mappings
clocks$hg19chrom <- NA
clocks$hg19start <- NA
clocks$hg19end <- NA
clocks$hg19HMMI <- NA
clocks$hg19ENSR <- NA

# hg38 mappings
clocks$hg38chrom <- NA
clocks$hg38start <- NA
clocks$hg38end <- NA
clocks$hg38HMMI <- NA
clocks$hg38ENSR <- NA

# hidden markov model (HMM) based CpG islands (systematic approach)
if (!exists("HMM_CpG_islands.hg19")) data(HMM_CpG_islands.hg19)
if (!exists("HMM_CpG_islands.hg38")) data(HMM_CpG_islands.hg38)

# add this information
for (i in assemblies) { 
  for (j in chips) {
    mfest <- subset(mfests[[i]][[j]], names %in% CpGs)
    clocks[names(mfest), paste0(i, "chrom")] <- seqnames(mfest)
    clocks[names(mfest), paste0(i, "start")] <- start(mfest)
    clocks[names(mfest), paste0(i, "end")] <- end(mfest)

    # are features within CpG islands? Which ones?
    HMMIs <- get(paste("HMM_CpG_islands", i, sep="."))
    HMMI_ol <- suppressWarnings(findOverlaps(HMMIs, mfest))
    HMMI_names <- names(HMMIs)[queryHits(HMMI_ol)]
    clocks[names(mfest)[subjectHits(HMMI_ol)], paste0(i, "HMMI")] <- HMMI_names
    
    # are features within ENSEMBL regulatory build regions? Which ones? 
    ENSRs <- regBuilds[[i]]
    ENSR_ol <- suppressWarnings(findOverlaps(ENSRs, mfest))
    ENSR_names <- names(ENSRs)[queryHits(ENSR_ol)]
    ENSR_types <- names(ENSRs)[queryHits(ENSR_ol)]
    clocks[names(mfest)[subjectHits(ENSR_ol)], paste0(i, "ENSR")] <- ENSR_names
  }
}

ENSR_subset.hg19 <- regBuilds$hg19
save(ENSR_subset.hg19, file="../../data/ENSR_subset.hg19.rda")
ENSR_subset.hg38 <- regBuilds$hg38
save(ENSR_subset.hg38, file="../../data/ENSR_subset.hg38.rda")
save(clocks, file="../../data/clocks.rda")
