---
title: "Biscuiteer User Guide"
date: "13 September 2019"
package: "biscuiteer 0.9.93"
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
## Loading required package: bsseq
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
tcga_bed <- system.file("extdata", "TCGA_BLCA_A13J_chr11p15_merged.bed.gz",
                        package = "biscuiteer")
tcga_vcf <- system.file("extdata", "TCGA_BLCA_A13J_header_only.vcf.gz",
                        package = "biscuiteer")
bisc <- read.biscuit(BEDfile = tcga_bed, VCFfile = tcga_vcf,
                     merged = TRUE, genome = "hg38")
```

```
## Checking /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_merged.bed.gz for import...
```

```
## Extracting sample names from /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_header_only.vcf.gz...
```

```
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_merged.bed.gz does not have a header. Using VCF file header information to help set column names.
```

```
## Assuming merged data. Checking now... ...The file might be alright. Double check if you're worried.
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_merged.bed.gz has 211270 indexed loci.
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_merged.bed.gz looks valid for import.
## Reading merged input from /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_merged.bed.gz...
## Excluding CpG sites with uniformly zero coverage...
## Loaded /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_merged.bed.gz. Creating bsseq object...
```

Metadata from the `biscuit` output can be viewed via:

```r
biscuitMetadata(bisc)
```

```
## CharacterList of length 3
## [["Reference genome"]] hg38.fa
## [["Biscuit version"]] 0.3.5.20180313
## [["Invocation"]] biscuit pileup -q 5 /primary/vari/genomicdata/genomes/h...
```

If further information about the VCF header is desired,

```r
metadata(bisc)$vcfHeader
```

```
## class: VCFHeader 
## samples(1): TCGA_BLCA_A13J_markdup
## meta(5): fileformat reference source contig program
## fixed(1): FILTER
## info(4): NS CX N5 AB
## geno(9): GT DP ... GL1 GQ
```

## Combining Results

In the instance where you have two separate BED files that you would like to
analyze in a single bsseq object, you can combine the files using `unionize`,
which is a wrapper around the bsseq function, `combine`.

```r
tcga_mrg <- system.file("extdata",
                        "TCGA_BLCA_A13J_chr11p15_unmerged.bed.gz",
                        package = "biscuiteer")
tcga_shf <- system.file("extdata",
                        "TCGA_BLCA_A13J_chr11p15_shuffled_unmerged.bed.gz",
                         package = "biscuiteer")
tcga_mvcf <- system.file("extdata",
                         "TCGA_BLCA_A13J_header_only.vcf.gz",
                         package = "biscuiteer")
tcga_svcf <- system.file("extdata",
                         "TCGA_BLCA_A13J_shuffled_header_only.vcf.gz",
                         package = "biscuiteer")
bisc1 <- read.biscuit(BEDfile = tcga_mrg, VCFfile = tcga_mvcf,
                      merged = FALSE, genome = "hg38")
```

```
## Checking /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_unmerged.bed.gz for import...
```

```
## Extracting sample names from /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_header_only.vcf.gz...
```

```
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_unmerged.bed.gz does not have a header. Using VCF file header information to help set column names.
```

```
## Assuming unmerged data. Checking now... ...The file might be alright. Double check if you're worried.
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_unmerged.bed.gz has 211350 indexed loci.
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_unmerged.bed.gz looks valid for import.
## Reading unmerged input from /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_unmerged.bed.gz...
## Excluding CpG sites with uniformly zero coverage...
## Loaded /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_unmerged.bed.gz. Creating bsseq object...
```

```r
bisc2 <- read.biscuit(BEDfile = tcga_shf, VCFfile = tcga_svcf,
                      merged = FALSE, genome = "hg38")
```

```
## Checking /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_shuffled_unmerged.bed.gz for import...
## Extracting sample names from /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_shuffled_header_only.vcf.gz...
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_shuffled_unmerged.bed.gz does not have a header. Using VCF file header information to help set column names.
## Assuming unmerged data. Checking now... ...The file might be alright. Double check if you're worried.
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_shuffled_unmerged.bed.gz has 211350 indexed loci.
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_shuffled_unmerged.bed.gz looks valid for import.
## Reading unmerged input from /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_shuffled_unmerged.bed.gz...
## Excluding CpG sites with uniformly zero coverage...
## Loaded /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_shuffled_unmerged.bed.gz. Creating bsseq object...
```

```r
comb <- unionize(bisc1, bisc2)
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
tcga_bed <- system.file("extdata", "TCGA_BLCA_A13J_chr11p15_merged.bed.gz",
                        package = "biscuiteer")
tcga_vcf <- system.file("extdata", "TCGA_BLCA_A13J_header_only.vcf.gz",
                        package = "biscuiteer")
bisc <- read.biscuit(BEDfile = tcga_bed, VCFfile = tcga_vcf,
                     merged = TRUE, genome = "hg38")
```

```
## Checking /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_merged.bed.gz for import...
```

```
## Extracting sample names from /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_header_only.vcf.gz...
```

```
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_merged.bed.gz does not have a header. Using VCF file header information to help set column names.
```

```
## Assuming merged data. Checking now... ...The file might be alright. Double check if you're worried.
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_merged.bed.gz has 211270 indexed loci.
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_merged.bed.gz looks valid for import.
## Reading merged input from /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_merged.bed.gz...
## Excluding CpG sites with uniformly zero coverage...
## Loaded /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_merged.bed.gz. Creating bsseq object...
```

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
##                         TCGA_BLCA_A13J_markdup
## chr11:0-2800000                      0.5518970
## chr11:2800000-11700000               0.7735474
## chr11:11700000-13800000              0.9581771
## chr11:13800000-16900000              0.6452412
## chr11:16900000-22000000              0.6409875
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
tcga_mrg <- system.file("extdata",
                        "TCGA_BLCA_A13J_chr11p15_unmerged.bed.gz",
                        package = "biscuiteer")
tcga_shf <- system.file("extdata",
                        "TCGA_BLCA_A13J_chr11p15_shuffled_unmerged.bed.gz",
                         package = "biscuiteer")
tcga_mvcf <- system.file("extdata",
                         "TCGA_BLCA_A13J_header_only.vcf.gz",
                         package = "biscuiteer")
tcga_svcf <- system.file("extdata",
                         "TCGA_BLCA_A13J_shuffled_header_only.vcf.gz",
                         package = "biscuiteer")
bisc1 <- read.biscuit(BEDfile = tcga_mrg, VCFfile = tcga_mvcf,
                      merged = FALSE, genome = "hg38")
```

```
## Checking /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_unmerged.bed.gz for import...
```

```
## Extracting sample names from /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_header_only.vcf.gz...
```

```
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_unmerged.bed.gz does not have a header. Using VCF file header information to help set column names.
```

```
## Assuming unmerged data. Checking now... ...The file might be alright. Double check if you're worried.
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_unmerged.bed.gz has 211350 indexed loci.
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_unmerged.bed.gz looks valid for import.
## Reading unmerged input from /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_unmerged.bed.gz...
## Excluding CpG sites with uniformly zero coverage...
## Loaded /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_unmerged.bed.gz. Creating bsseq object...
```

```r
bisc2 <- read.biscuit(BEDfile = tcga_shf, VCFfile = tcga_svcf,
                      merged = FALSE, genome = "hg38")
```

```
## Checking /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_shuffled_unmerged.bed.gz for import...
## Extracting sample names from /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_shuffled_header_only.vcf.gz...
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_shuffled_unmerged.bed.gz does not have a header. Using VCF file header information to help set column names.
## Assuming unmerged data. Checking now... ...The file might be alright. Double check if you're worried.
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_shuffled_unmerged.bed.gz has 211350 indexed loci.
## /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_shuffled_unmerged.bed.gz looks valid for import.
## Reading unmerged input from /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_shuffled_unmerged.bed.gz...
## Excluding CpG sites with uniformly zero coverage...
## Loaded /secondary/projects/shen/tools/morrison/anaconda3/envs/r_env_3.6/lib/R/library/biscuiteer/extdata/TCGA_BLCA_A13J_chr11p15_shuffled_unmerged.bed.gz. Creating bsseq object...
```

```r
comb <- unionize(bisc1, bisc2)
ages <- WGBSage(comb, "horvath")
```

```
## Assessing coverage across age-associated regions...
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
##                         TCGA_BLCA_A13J_markdup TCGA_BLCA_A13J_shuffled
## chr11:6656898-6656927                        0               0.1000000
## chr11:12009082-12009111                      0               0.4607843
## 
## $coefs
##   chr11:6656898-6656927 chr11:12009082-12009111 
##             0.000792206            -0.138857398 
## 
## $age
##                             [,1]
## TCGA_BLCA_A13J_markdup  35.60565
## TCGA_BLCA_A13J_shuffled 34.26367
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
## Loading HMM_CpG_islands.hg38...
```

```
## Loading H9state23unmeth.hg38...
```

```
## Computing hypomethylation indices...
```

```
## Loading PMDs.hg38...
```

```
## Loading Zhou_solo_WCGW_inCommonPMDs.hg38...
```

```
## Computing indices...
```

```r
show(bisc.CpGindex)
```

```
## CpGindex with 1 row and 3 columns
##   hyper.TCGA_BLCA_A13J_markdup hypo.TCGA_BLCA_A13J_markdup
##                      <numeric>                   <numeric>
## 1           0.0654991319444445           0.443562305099473
##   ratio.TCGA_BLCA_A13J_markdup
##                      <numeric>
## 1            0.147666136620324
##   -------
## This object is just a DataFrame that has an idea of where it came from:
## Hypermethylation was tallied across 121 regions (see bisc.CpGindex@hyperMethRegions). 
## Hypomethylation was tallied across 13123 regions (see bisc.CpGindex@hypoMethRegions).
```

```r
bisc.CpGindex@hyperMethRegions
```

```
## GRanges object with 121 ranges and 0 metadata columns:
##       seqnames            ranges strand
##          <Rle>         <IRanges>  <Rle>
##     1     chr1 31764600-31764623      *
##     2     chr1 43172730-43172778      *
##     3     chr1 44418329-44418333      *
##     4     chr1 46394729-46394734      *
##     5     chr1 50970129-50970403      *
##   ...      ...               ...    ...
##   117    chr20   8131745-8131753      *
##   118    chr20 17227156-17227546      *
##   119    chr22 20017278-20017279      *
##   120    chr22 36856559-36856689      *
##   121    chr22 43385844-43385946      *
##   -------
##   seqinfo: 21 sequences from hg38 genome
```
