# biscuiteer

[![Build Status](https://travis-ci.org/ttriche/biscuiteer.png?branch=master)](https://travis-ci.org/ttriche/biscuiteer)  [![codecov](https://codecov.io/gh/ttriche/biscuiteer/branch/master/graph/badge.svg)](https://codecov.io/gh/ttriche/biscuiteer)

## The original luxury biscuit boutique

Wait, no, that's [these guys](https://www.biscuiteers.com/). ```biscuiteer```, on the other hand, is a package to process [biscuit](https://github.com/zwdzwd/biscuit) output quickly into [bsseq](https://bioconductor.org/packages/bsseq) objects. A number of features such as VCF header parsing, shrunken M-value calculations (for compartment inference), and tumor/normal copy number segmentation are also included, but the task of locus- and region-level differential methylation inference is delegated to other packages (such as ```dmrseq```).

## Installing

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("biscuiteer")
```

A development version is available on GitHub and can be installed via:
```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("trichelab/biscuiteerData")
BiocManager::install("trichelab/biscuiteer")
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
orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
                        package="biscuiteer")
orig_vcf <- system.file("extdata", "MCF7_Cunha_header_only.vcf.gz",
                        package="biscuiteer")
bisc <- readBiscuit(BEDfile = orig_bed, VCFfile = orig_vcf,
                    merged = FALSE)

# To print metadata information from the loaded file:
#
biscuitMetadata(bisc)
#
# CharacterList of length 3
# [["Reference genome"]] hg19.fa
# [["Biscuit version"]] 0.1.3.20160324
# [["Invocation"]] biscuit pileup -r /primary/vari/genomicdata/genomes/hg19/hg1...

# This is all drawn from the VCF header:
#
metadata(bisc)$vcfHeader
#
# class: VCFHeader 
# samples(1): MCF7_Cunha
# meta(5): fileformat reference source contig program
# fixed(1): FILTER
# info(3): NS CX N5
# geno(7): GT DP ... GL GQ
```

### Downstream bits 

A/B compartment inference, age estimation from WGBS,  and so forth (examples to appear).    
To wit, a shorthand summary of hypo- and hyper-methylation at some regions commonly associated with each:

```R
bisc.CpGindex <- CpGindex(bisc)
#
# Computing hypermethylation indices...
# Loading HMM_CpG_islands.hg19...
# Loading H9state23unmeth.hg19...
# Computing hypomethylation indices...
# Loading PMDs.hg19.rda from biscuiteerData...
# Loading Zhou_solo_WCGW_inCommonPMDs.hg19.rda from biscuiteerData...
# Computing indices...

show(bisc.CpGindex)
#
# CpGindex with 1 row and 3 columns
#     hyper.MCF7_Cunha   hypo.MCF7_Cunha ratio.MCF7_Cunha
#            <numeric>         <numeric>        <numeric>
# 1 0.0690734126984127 0.199261516805161 0.34664702851757
#   -------
# This object is just a DataFrame that has an idea of where it came from:
# Hypermethylation was tallied across 120 region (see 'object@hyperMethRegions'). 
# Hypomethylation was tallied across 13127 region (see 'object@hypoMethRegions').

bisc.CpGindex@hyperMethRegions
#
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
