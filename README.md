# biscuiteer

[![Build Status](https://travis-ci.org/ttriche/biscuiteer.png?branch=master)](https://travis-ci.org/ttriche/biscuiteer)  [![codecov](https://codecov.io/gh/ttriche/biscuiteer/branch/master/graph/badge.svg)](https://codecov.io/gh/ttriche/biscuiteer)

## The original luxury biscuit boutique

Wait, no, that's [these guys](https://www.biscuiteers.com/). ```biscuiteer```, on the other hand, is a package to process [biscuit](https://github.com/zwdzwd/biscuit) output quickly into [bsseq](https://bioconductor.org/packages/bsseq) objects. A number of features such as VCF header parsing, shrunken M-value calculations (for compartment inference), and tumor/normal copy number segmentation are also included, but the task of locus- and region-level differential methylation inference is delegated to other packages (such as ```dmrseq```).

## Installing

```R
    install.packages("BiocManager")
    library(BiocManager)
    install("trichelab/biscuiteer")
```

## Usage

### Loading data 

```biscuiteer``` can load headered or header-free BED-like files as produced from ```biscuit vcf2bed``` or ```biscuit mergecg```, but we encourage users to keep their VCF headers (or just the entire VCF, which you will want to do anyways, as biscuit calls SNVs and allows for structural variant detection in a manner similar to typical whole-genome sequencing tools).  Since ```biscuit``` records the version of the software and the calling arguments used to process a set of files in the output VCF, this allows for much better reproducibility:

```R
# from SRP080893
# To load some runs:
#
library(biscuiteer)

# We just need the BED output and the VCF header to proceed... 
# It is possible to get away without the VCF header, but better to use it. 
#
SvNS <- read.biscuit(BEDfile="SvNS_HSPC_WGBS.merged.hg19.bed.gz",
                     VCFfile="SvNS_HSPC_WGBS.hg19.vcf.gz",
                     simplify=TRUE)

# for reproducibility:
#
biscuitMetadata(SvNS)
#
# CharacterList of length 3
# [["Reference genome"]] hg19.fa
# [["Biscuit version"]] 0.3.8.20180515
# [["Invocation"]] biscuit pileup -o SvNS_HSPC_WGBS.hg19.vcf /home/tim.triche...

# this is all drawn from the VCF header:
#
metadata(SvNS)$vcfHeader
# 
# class: VCFHeader 
# samples(6): 052314-2-N.hg19.sorted.markDups.hg19
#   052314-2-S.hg19.sorted.markDups.hg19 ...
```

### Downstream bits 

A/B compartment inference, age estimation from WGBS, copy number segmentation, and so forth (examples to appear).    
To wit, a shorthand summary of hypo- and hyper-methylation at some regions commonly associated with each:

```R
SvNS.CpGindex <- CpGindex(SvNS) 
# Computing hypermethylation indices...
# Loading HMM_CpG_islands.hg19...
# Loading H9state23unmeth.hg19...
# Computing hypomethylation indices...
# Loading PMDs.hg19...
# Loading Zhou_solo_WCGW_inCommonPMDs.hg19...
# Computing indices...

SvNS.CpGindex
# CpGindex with 6 rows and 3 columns
#                         hyper              hypo              ratio
#                     <numeric>         <numeric>          <numeric>
# 052314-2-N 0.0257280213640668  0.85695245688856 0.0300226939747399
# 052314-2-S 0.0199967524923102 0.778438421291039 0.0256882907438531
# 052314-3-N 0.0259634680610258 0.866410364528792 0.0299667099148177
# 052314-3-S 0.0223123531642581 0.862420672085278 0.0258717745138324
# 060614-2-N 0.0352164982188626 0.858982911495614 0.0409979031568225
# 060614-2-S  0.026870818880881 0.858271082834001 0.0313080790187569
#   -------
# This object is just a DataFrame that has an idea of where it came from:
# Hypermethylation was tallied across 120 regions (see x@hyperMethRegions). 
# Hypomethylation was tallied across 15315 regions (see x@hypoMethRegions). 

SvNS.CpGindex@hyperMethRegions
# GRanges object with 120 ranges and 1 metadata column:
#       seqnames            ranges strand |              score
#          <Rle>         <IRanges>  <Rle> |          <numeric>
#     1     chr1 32230201-32230224      * | 0.0399999991059303
#     2     chr1 43638401-43638449      * | 0.0399999991059303
#     3     chr1 44884001-44884005      * | 0.0399999991059303
#     4     chr1 46860401-46860406      * | 0.0599999986588955
#     5     chr1 51435801-51436075      * | 0.0499999998137355
#   ...      ...               ...    ... .                ...
#   116    chr20   8112392-8112400      * | 0.0299999993294477
#   117    chr20 17207801-17208191      * | 0.0599999986588955
#   118    chr22 20004801-20004802      * | 0.0350000001490116
#   119    chr22 37252601-37252731      * | 0.0599999986588955
#   120    chr22 43781850-43781952      * | 0.0449999999254942
#   -------
#   seqinfo: 21 sequences from hg19 genome

```


### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.
