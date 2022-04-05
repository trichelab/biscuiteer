---
title: "Biscuiteer User Guide"
date: "5 April 2022"
package: "biscuiteer 1.9.6"
output:
  BiocStyle::html_document:
    highlight: pygments
    toc_float: true
    fig_width: 8
    fig_height: 6
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{Biscuiteer User Guide}
  %\VignetteEncoding[utf8]{inputenc}
---

knitr::opts_chunk$set(
    warning=FALSE,
    message=FALSE
  )

# Biscuiteer

`biscuiteer` is package to process output from
[biscuit](https://github.com/huishenlab/biscuit) into
[bsseq](https://bioconductor.org/packages/bsseq) objects. It includes a number
of features, such as VCF header parsing, shrunken M-value calculations (which
can be used for compartment inference), and age inference. However, the task of
locus- and region-level differential methylation inference is delegated to
other packages (such as `dmrseq`).

# Quick Start

## Installing

From Bioconductor,

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("biscuiteer")
```

A development version is available on GitHub and can be installed via:

```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("trichelab/biscuiteerData")
BiocManager::install("trichelab/biscuiteer")
```

## Loading Methylation Data

`biscuiteer` can load either headered or header-free BED files produced from
`biscuit vcf2bed` or `biscuit mergecg`. In either case, a VCF file is needed
when loading `biscuit` output. For practical purposes, only the VCF header is
for `biscuiteer`. However, it is encouraged that the user keep the entire VCF,
as `biscuit` can be used to call SNVs and allows for structural variant
detection in a similar manner to typical whole-genome sequencing tools.
Furthermore, `biscuit` records the version of the software and the calling
arguments used during processing the output VCF, which allows for better
reproducibility.

NOTE: Both the input BED and VCF files must be tabix'ed before being input to
`biscuiteer`. This can be done by running `bgzip biscuit_output.xxx` followed by
`tabix -p xxx biscuit_output.xxx.gz`, where `xxx` is either `bed` or `vcf`.

Data can be loaded using the `readBiscuit` function in `biscuiteer`:

```r
library(biscuiteer)
```

```
## Loading required package: biscuiteerData
```

```
## Loading required package: ExperimentHub
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which.max, which.min
```

```
## Loading required package: AnnotationHub
```

```
## Loading required package: BiocFileCache
```

```
## Loading required package: dbplyr
```

```
## Loading biscuiteerData.
```

```
## Loading required package: bsseq
```

```
## Loading required package: GenomicRanges
```

```
## Loading required package: stats4
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
```

```
## Loading required package: IRanges
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: SummarizedExperiment
```

```
## Loading required package: MatrixGenerics
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'MatrixGenerics'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## 
## Attaching package: 'Biobase'
```

```
## The following object is masked from 'package:MatrixGenerics':
## 
##     rowMedians
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     anyMissing, rowMedians
```

```
## The following object is masked from 'package:ExperimentHub':
## 
##     cache
```

```
## The following object is masked from 'package:AnnotationHub':
## 
##     cache
```

```
## 
```

```
## 
```

```
## 
```

```
## Warning: replacing previous import 'BiocParallel::bpstart' by 'QDNAseq::bpstart'
## when loading 'biscuiteer'
```

```
## Warning: replacing previous import 'AnnotationHub::hubUrl' by
## 'rtracklayer::hubUrl' when loading 'annotatr'
```

```r
orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
                        package="biscuiteer")
orig_vcf <- system.file("extdata", "MCF7_Cunha_header_only.vcf.gz",
                        package="biscuiteer")
bisc <- readBiscuit(BEDfile = orig_bed, VCFfile = orig_vcf,
                    merged = FALSE)
```

```
## Checking /Library/Frameworks/R.framework/Versions/4.1/Resources/library/biscuiteer/extdata/MCF7_Cunha_chr11p15.bed.gz for import...
```

```
## Extracting sample names from /Library/Frameworks/R.framework/Versions/4.1/Resources/library/biscuiteer/extdata/MCF7_Cunha_header_only.vcf.gz...
```

```
## /Library/Frameworks/R.framework/Versions/4.1/Resources/library/biscuiteer/extdata/MCF7_Cunha_chr11p15.bed.gz does not have a header. Using VCF file header information to help set column names.
```

```
## Assuming unmerged data. Checking now... ...The file might be alright. Double check if you're worried.
## /Library/Frameworks/R.framework/Versions/4.1/Resources/library/biscuiteer/extdata/MCF7_Cunha_chr11p15.bed.gz has 254147 indexed loci.
## /Library/Frameworks/R.framework/Versions/4.1/Resources/library/biscuiteer/extdata/MCF7_Cunha_chr11p15.bed.gz looks valid for import.
## Reading unmerged input from /Library/Frameworks/R.framework/Versions/4.1/Resources/library/biscuiteer/extdata/MCF7_Cunha_chr11p15.bed.gz...
## Excluding CpG sites with uniformly zero coverage...
## Loaded /Library/Frameworks/R.framework/Versions/4.1/Resources/library/biscuiteer/extdata/MCF7_Cunha_chr11p15.bed.gz. Creating bsseq object......Done!
```

Metadata from the `biscuit` output can be viewed via:

```r
biscuitMetadata(bisc)
```

```
## CharacterList of length 3
## [["Reference genome"]] hg19.fa
## [["Biscuit version"]] 0.1.3.20160324
## [["Invocation"]] biscuit pileup -r /primary/vari/genomicdata/genomes/hg19/hg1...
```

If further information about the VCF header is desired,

```r
metadata(bisc)$vcfHeader
```

```
## class: VCFHeader 
## samples(1): MCF7_Cunha
## meta(5): fileformat reference source contig program
## fixed(1): FILTER
## info(3): NS CX N5
## geno(7): GT DP ... GL GQ
```

## Combining Methylation Results

In the instance where you have two separate BED files that you would like to
analyze in a single bsseq object, you can combine the files using `unionize`,
which is a wrapper around the BiocGenerics function, `combine`.

```r
shuf_bed <- system.file("extdata", "MCF7_Cunha_chr11p15_shuffled.bed.gz",
                        package="biscuiteer")
shuf_vcf <- system.file("extdata",
                        "MCF7_Cunha_shuffled_header_only.vcf.gz",
                        package="biscuiteer")
bisc2 <- readBiscuit(BEDfile = shuf_bed, VCFfile = shuf_vcf,
                     merged = FALSE)
```

```
## Checking /Library/Frameworks/R.framework/Versions/4.1/Resources/library/biscuiteer/extdata/MCF7_Cunha_chr11p15_shuffled.bed.gz for import...
```

```
## Extracting sample names from /Library/Frameworks/R.framework/Versions/4.1/Resources/library/biscuiteer/extdata/MCF7_Cunha_shuffled_header_only.vcf.gz...
```

```
## /Library/Frameworks/R.framework/Versions/4.1/Resources/library/biscuiteer/extdata/MCF7_Cunha_chr11p15_shuffled.bed.gz does not have a header. Using VCF file header information to help set column names.
```

```
## Assuming unmerged data. Checking now... ...The file might be alright. Double check if you're worried.
## /Library/Frameworks/R.framework/Versions/4.1/Resources/library/biscuiteer/extdata/MCF7_Cunha_chr11p15_shuffled.bed.gz has 254147 indexed loci.
## /Library/Frameworks/R.framework/Versions/4.1/Resources/library/biscuiteer/extdata/MCF7_Cunha_chr11p15_shuffled.bed.gz looks valid for import.
## Reading unmerged input from /Library/Frameworks/R.framework/Versions/4.1/Resources/library/biscuiteer/extdata/MCF7_Cunha_chr11p15_shuffled.bed.gz...
## Excluding CpG sites with uniformly zero coverage...
## Loaded /Library/Frameworks/R.framework/Versions/4.1/Resources/library/biscuiteer/extdata/MCF7_Cunha_chr11p15_shuffled.bed.gz. Creating bsseq object......Done!
```

```r
comb <- unionize(bisc, bisc2)
```

## Loading epiBED files

The epiBED file format provides an easy way to analyze read- or fragment-level
methylation and genetic information at the same time. `readEpibed` provides
functionality for parsing the RLE strings found in the epiBED file into a
GRanges object for analysis in R.

NOTE: The input file must be bgzip'ed and tabix'ed.

```r
epibed.nome <- system.file("extdata", "hct116.nome.epiread.gz",
                           package="biscuiteer")
epibed.bsseq <- system.file("extdata", "hct116.bsseq.epiread.gz",
                            package="biscuiteer")
epibed.nome.gr <- readEpibed(epibed = epibed.nome, is.nome = TRUE,
                             genome = "hg19", chr = "chr1")
```

```
## Decoding RLE and converting to GRanges
```

```r
epibed.bsseq.gr <- readEpibed(epibed = epibed.bsseq,
                              genome = "hg19", chr = "chr1")
```

```
## Decoding RLE and converting to GRanges
```

# Analysis Functionality

A handful of analysis paths are available in `biscuiteer`, including A/B
comparment inference, age estimation from WGBS data, hypermethylation of
Polycomb Repressor Complex (PRC) binding sites, and hypomethylation of
CpG-poor "partially methylated domains" (PMDs).

## Inputs for A/B Compartment Inference

When performing A/B compartment inference, the goal is to have something that
has roughly gaussian error. `getLogitFracMeth` uses Dirichlet smoothing to turn
raw measurements into lightly moderated, logit-transformed methylated-fraction
estimates, which can the be used as inputs to
[compartmap](https://bioconductor.org/packages/release/bioc/html/compartmap.md)


```r
reg <- GRanges(seqnames = rep("chr11",5),
               strand = rep("*",5),
               ranges = IRanges(start = c(0,2.8e6,1.17e7,1.38e7,1.69e7),
                                end= c(2.8e6,1.17e7,1.38e7,1.69e7,2.2e7))
              )

frac <- getLogitFracMeth(bisc, minSamp = 1, r = reg)
frac
```

```
## GRanges object with 5 ranges and 1 metadata column:
##       seqnames            ranges strand | MCF7_Cunha
##          <Rle>         <IRanges>  <Rle> |  <numeric>
##   [1]    chr11         0-2800000      * |   1.340682
##   [2]    chr11  2800000-11700000      * |   0.575875
##   [3]    chr11 11700000-13800000      * |   1.162989
##   [4]    chr11 13800000-16900000      * |   0.581874
##   [5]    chr11 16900000-22000000      * |   0.442985
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

## Age Estimation

`biscuiteer` has the functionalitity to guess the age of the sample(s) provided
using the Horvath-style "clock" models (see
[Horvath, 2013](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-10-r115)
for more information).

NOTE: The prediction accuracy of this function is entirely dependent on the
parameters set by the user. As such, the defaults (as shown in the example
below) should only be used as a starting point for exploration by the user.

NOTE: Please cite the appropriate papers for the epigenetic "clock" chosen:
* For `horvath` or `horvathshrunk`
    * Horvath, Genome Biology, 2013
* For `hannum`
    * Hannum et al., Molecular Cell, 2013
* For `skinandblood`
    * Horvath et al., Aging, 2018


```r
ages <- WGBSage(comb, "horvath")
```

```
## Assessing coverage across age-associated regions...
```

```
## All regions in all samples appear to be sufficiently covered.
```

```r
ages
```

```
## $call
## WGBSage(comb, "horvath")
## 
## $droppedSamples
## NULL
## 
## $droppedRegions
## NULL
## 
## $intercept
## [1] 0.6955073
## 
## $methcoefs
## GRanges object with 2 ranges and 3 metadata columns:
##                           seqnames            ranges strand | MCF7_Cunha
##                              <Rle>         <IRanges>  <Rle> |  <numeric>
##     chr11:6678129-6678158    chr11   6678129-6678158      * |   0.800000
##   chr11:12030629-12030658    chr11 12030629-12030658      * |   0.833333
##                           MCF7_Cunha_shuffled        coefs
##                                     <numeric>    <numeric>
##     chr11:6678129-6678158            0.250000  0.000792206
##   chr11:12030629-12030658            0.247732 -0.138857398
##   -------
##   seqinfo: 22 sequences from hg19 genome
## 
## $age
##                         [,1]
## MCF7_Cunha          33.18896
## MCF7_Cunha_shuffled 34.88742
```
