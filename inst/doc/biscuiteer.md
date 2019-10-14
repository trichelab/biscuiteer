---
title: "Biscuiteer User Guide"
date: "14 October 2019"
package: "biscuiteer 0.99.1"
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

# Biscuiteer

`biscuiteer` is package to process output from
[biscuit](https://github.com/zwdzwd/biscuit) into
[bsseq](https://bioconductor.org/packages/bsseq) objects. It includes a number
of features, such as VCF header parsing, shrunken M-value calculations (which
can be used for compartment inference), and age inference are included. However,
the task of locus- and region-level differential methylation inference is
delegated to other packages (such as `dmrseq`).

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

## Loading Data

`biscuiteer` can load either headered of header-free BED files produced from
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

Data can be loaded using the `read.biscuit` function in `biscuiteer`:

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
##     anyDuplicated, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter,
##     Find, get, grep, grepl, intersect, is.unsorted, lapply, Map,
##     mapply, match, mget, order, paste, pmax, pmax.int, pmin,
##     pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
##     setdiff, sort, table, tapply, union, unique, unsplit, which,
##     which.max, which.min
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
## The following object is masked from 'package:base':
## 
##     expand.grid
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
## Loading required package: DelayedArray
```

```
## Loading required package: matrixStats
```

```
## 
## Attaching package: 'matrixStats'
```

```
## The following objects are masked from 'package:Biobase':
## 
##     anyMissing, rowMedians
```

```
## Loading required package: BiocParallel
```

```
## 
## Attaching package: 'DelayedArray'
```

```
## The following objects are masked from 'package:matrixStats':
## 
##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges
```

```
## The following objects are masked from 'package:base':
## 
##     aperm, apply, rowsum
```

```
## Registered S3 method overwritten by 'R.oo':
##   method        from       
##   throw.default R.methodsS3
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

```r
orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
                        package="biscuiteer")
orig_vcf <- system.file("extdata", "MCF7_Cunha_header_only.vcf.gz",
                        package="biscuiteer")
bisc <- read.biscuit(BEDfile = orig_bed, VCFfile = orig_vcf,
                     merged = FALSE)
```

```
## Checking /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/MCF7_Cunha_chr11p15.bed.gz for import...
```

```
## Extracting sample names from /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/MCF7_Cunha_header_only.vcf.gz...
```

```
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/MCF7_Cunha_chr11p15.bed.gz does not have a header. Using VCF file header information to help set column names.
```

```
## Assuming unmerged data. Checking now... ...The file might be alright. Double check if you're worried.
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/MCF7_Cunha_chr11p15.bed.gz has 254147 indexed loci.
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/MCF7_Cunha_chr11p15.bed.gz looks valid for import.
## Reading unmerged input from /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/MCF7_Cunha_chr11p15.bed.gz...
## Excluding CpG sites with uniformly zero coverage...
## Loaded /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/MCF7_Cunha_chr11p15.bed.gz. Creating bsseq object...
```

Metadata from the `biscuit` output can be viewed via:

```r
biscuitMetadata(bisc)
```

```
## CharacterList of length 3
## [["Reference genome"]] hg19.fa
## [["Biscuit version"]] 0.1.3.20160324
## [["Invocation"]] biscuit pileup -r /primary/vari/genomicdata/genomes/hg1...
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

## Combining Results

In the instance where you have two separate BED files that you would like to
analyze in a single bsseq object, you can combine the files using `unionize`,
which is a wrapper around the bsseq function, `combine`.

```r
shuf_bed <- system.file("extdata", "MCF7_Cunha_chr11p15_shuffled.bed.gz",
                        package="biscuiteer")
shuf_vcf <- system.file("extdata",
                        "MCF7_Cunha_shuffled_header_only.vcf.gz",
                        package="biscuiteer")
bisc2 <- read.biscuit(BEDfile = shuf_bed, VCFfile = shuf_vcf,
                      merged = FALSE)
```

```
## Checking /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/MCF7_Cunha_chr11p15_shuffled.bed.gz for import...
```

```
## Extracting sample names from /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/MCF7_Cunha_shuffled_header_only.vcf.gz...
```

```
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/MCF7_Cunha_chr11p15_shuffled.bed.gz does not have a header. Using VCF file header information to help set column names.
```

```
## Assuming unmerged data. Checking now... ...The file might be alright. Double check if you're worried.
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/MCF7_Cunha_chr11p15_shuffled.bed.gz has 254147 indexed loci.
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/MCF7_Cunha_chr11p15_shuffled.bed.gz looks valid for import.
## Reading unmerged input from /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/MCF7_Cunha_chr11p15_shuffled.bed.gz...
## Excluding CpG sites with uniformly zero coverage...
## Loaded /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/MCF7_Cunha_chr11p15_shuffled.bed.gz. Creating bsseq object...
```

```r
comb <- unionize(bisc, bisc2)
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
[compartmap](https://bioconductor.org/packages/release/bioc/html/compartmap.html)


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
##       seqnames            ranges strand |        MCF7_Cunha
##          <Rle>         <IRanges>  <Rle> |         <numeric>
##   [1]    chr11         0-2800000      * |  1.34068175424637
##   [2]    chr11  2800000-11700000      * | 0.575874918312675
##   [3]    chr11 11700000-13800000      * |   1.1629889150866
##   [4]    chr11 13800000-16900000      * | 0.581873982746806
##   [5]    chr11 16900000-22000000      * | 0.442984923927319
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
## $meth
##                         MCF7_Cunha MCF7_Cunha_shuffled
## chr11:6678129-6678158    0.8000000           0.2500000
## chr11:12030629-12030658  0.8333333           0.2477324
## 
## $coefs
##   chr11:6678129-6678158 chr11:12030629-12030658 
##             0.000792206            -0.138857398 
## 
## $age
##                         [,1]
## MCF7_Cunha          33.18896
## MCF7_Cunha_shuffled 34.88742
```

## Hypermethylation of PRCs and Hypomethylation of PMDs

To generate a shorthand summary of the hypermethylation of PRCs and
hypomethylation of PMDs, the `CpGindex` function in `biscuiteer` can be used.


```r
bisc.CpGindex <- CpGindex(bisc)
```

```
## Computing hypermethylation indices...
```

```
## Loading HMM_CpG_islands.hg19...
```

```
## Loading H9state23unmeth.hg19...
```

```
## Computing hypomethylation indices...
```

```
## Loading PMDs.hg19.rda from biscuiteerData...
```

```
## Loading Zhou_solo_WCGW_inCommonPMDs.hg19.rda from biscuiteerData...
```

```
## Computing indices...
```

```r
show(bisc.CpGindex)
```

```
## CpGindex with 1 row and 3 columns
##     hyper.MCF7_Cunha   hypo.MCF7_Cunha ratio.MCF7_Cunha
##            <numeric>         <numeric>        <numeric>
## 1 0.0690734126984127 0.199261516805161 0.34664702851757
##   -------
## This object is just a DataFrame that has an idea of where it came from:
## Hypermethylation was tallied across 120 regions (see bisc.CpGindex@hyperMethRegions). 
## Hypomethylation was tallied across 13127 regions (see bisc.CpGindex@hypoMethRegions).
```

```r
bisc.CpGindex@hyperMethRegions
```

```
## GRanges object with 120 ranges and 1 metadata column:
##       seqnames            ranges strand |              score
##          <Rle>         <IRanges>  <Rle> |          <numeric>
##     1     chr1 32230201-32230224      * | 0.0399999991059303
##     2     chr1 43638401-43638449      * | 0.0399999991059303
##     3     chr1 44884001-44884005      * | 0.0399999991059303
##     4     chr1 46860401-46860406      * | 0.0599999986588955
##     5     chr1 51435801-51436075      * | 0.0499999998137355
##   ...      ...               ...    ... .                ...
##   116    chr20   8112392-8112400      * | 0.0299999993294477
##   117    chr20 17207801-17208191      * | 0.0599999986588955
##   118    chr22 20004801-20004802      * | 0.0350000001490116
##   119    chr22 37252601-37252731      * | 0.0599999986588955
##   120    chr22 43781850-43781952      * | 0.0449999999254942
##   -------
##   seqinfo: 21 sequences from hg19 genome
```
