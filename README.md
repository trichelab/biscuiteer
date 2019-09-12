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
# Load the package
#
library(biscuiteer)

# To read in some data, you need:
#   A BED file (with or without a header)
#   A VCF file (only really needs the header)
#   Whether the file is merged CG data or not
#   Any other additional function inputs (genome, region to load, etc)
#
tcga_bed <- system.file("extdata", "TCGA_BLCA_A13J_chr11p15_merged.bed.gz",
                        package = "biscuiteer")
tcga_vcf <- system.file("extdata", "TCGA_BLCA_A13J_header_only.vcf.gz",
                        package = "biscuiteer")
bisc <- read.biscuit(BEDfile = tcga_bed, VCFfile = tcga_vcf,
                     merged = TRUE, genome = "hg38")

# To print metadata information from the loaded file:
#
biscuitMetadata(bisc)
#
# CharacterList of length 3
# [["Reference genome"]] hg38.fa
# [["Biscuit version"]] 0.3.5.20180313
# [["Invocation"]] biscuit pileup -q 5 /primary/vari/genomicdata/genomes/hg38/h...

# This is all drawn from the VCF header:
#
metadata(bisc)$vcfHeader
#
# class: VCFHeader 
# samples(1): TCGA_BLCA_A13J_markdup
# meta(5): fileformat reference source contig program
# fixed(1): FILTER
# info(4): NS CX N5 AB
# geno(9): GT DP ... GL1 GQ
```

### Downstream bits 

A/B compartment inference, age estimation from WGBS, copy number segmentation, and so forth (examples to appear).    
To wit, a shorthand summary of hypo- and hyper-methylation at some regions commonly associated with each:

```R
bisc.CpGindex <- CpGindex(bisc)
#
# Computing hypermethylation indices...
# Loading HMM_CpG_islands.hg38...
# Loading H9state23unmeth.hg38...
# Computing hypomethylation indices...
# Loading PMDs.hg38...
# Loading Zhou_solo_WCGW_inCommonPMDs.hg38...
# Computing indices...

show(bisc.CpGindex)
#
# CpGindex with 1 row and 3 columns
#      hyper.TCGA_BLCA_A13J_markdup hypo.TCGA_BLCA_A13J_markdup
#      `<`numeric`>`                `<`numeric`>`
#   1  0.0654991319444445           0.443562305099473
#      ratio.TCGA_BLCA_A13J_markdup
#      `<`numeric`>`
#   1  0.147666136620324
#   -------
#   This object is just a DataFrame that has an idea of where it came from:
#   Hypermethylation was tallied across 121 regions (see bisc.CpGindex@hyperMethRegions). 
#   Hypomethylation was tallied across 13123 regions (see bisc.CpGindex@hypoMethRegions).

bisc.CpGindex@hyperMethRegions
#
# GRanges object with 121 ranges and 0 metadata columns:
#   seqnames            ranges strand
#  `<`Rle`>`     `<`IRanges`>`  `<`Rle`>`
# 1     chr1 31764600-31764623      *
# 2     chr1 43172730-43172778      *
# 3     chr1 44418329-44418333      *
# 4     chr1 46394729-46394734      *
# 5     chr1 50970129-50970403      *
# ...    ...               ...    ...
# 117  chr20   8131745-8131753      *
# 118  chr20 17227156-17227546      *
# 119  chr22 20017278-20017279      *
# 120  chr22 36856559-36856689      *
# 121  chr22 43385844-43385946      *
# -------
# seqinfo: 21 sequences from hg38 genome
```
### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.
