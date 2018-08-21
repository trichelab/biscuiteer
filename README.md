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

Examples of A/B compartment inference, age estimation from WGBS, copy number segmentation, and so forth will appear here presently.


### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.
