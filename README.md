# biscuiteer

[![Build Status](https://travis-ci.org/ttriche/biscuiteer.png?branch=master)](https://travis-ci.org/ttriche/biscuiteer)  [![codecov](https://codecov.io/gh/ttriche/biscuiteer/branch/master/graph/badge.svg)](https://codecov.io/gh/ttriche/biscuiteer)

## The original luxury biscuit boutique

Wait, no, that's [these guys](https://www.biscuiteers.com/).  ```biscuiteer``` is a package to process [biscuit](https://github.com/zwdzwd/biscuit) output very quickly into [bsseq](https://bioconductor.org/packages/bsseq) objects for downstream processing.

## Installing

```R
    install.packages("BiocManager")
    library(BiocManager)
    install("trichelab/biscuiteer")
```

### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.
