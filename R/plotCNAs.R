#' Analyze CNAs
#'
#' @param filenames  List of files with tumor and normal Excel sheets
#' @param minWidth   Scaling factor (DEFAULT: 1000)
#'
#' @return           Heatmap of CNAs
#'
#' @importFrom grid unit gpar
#' @import readxl
#' @import GenomeInfoDb
#' @import ComplexHeatmap
#' @import RaggedExperiment
#'
#' @export
#'
plotCNA <- function(filenames,
                    minWidth = 1000) { 
  
  excelCNAs <- .readExcelCNAs(filenames)

  rse <- .makeRse(excelCNAs)
  seqlevelsStyle(rse) <- "NCBI" # more compact
  rowData(rse)$units <- floor(width(rse) / minWidth)
  rse <- subset(rse, rowData(rse)$units > 0)

  toPlot <- assays(rse)$score
  toPlot[is.na(toPlot)] <- 0 # yuck
  toKeep <- which(rowSums(abs(toPlot)) != 0)
  rse <- sortSeqlevels(rse)
  chrom <- seqnames(rse) # keep as factor!
  toPlot <- toPlot[toKeep, ]
  chrom <- chrom[toKeep]

  Heatmap(t(toPlot), name="log2(CN)", 
          cluster_columns=FALSE, column_split=chrom, row_gap=unit(3, "mm"),
          clustering_method_rows="ward.D2", row_split=rse$subject,
          column_title_gp=gpar(fontsize=8), show_column_names=FALSE, 
          row_title_rot=0, na_col="#FFFFFF", column_gap=unit(3, "mm"))
}


# Helper function
.readExcelCNAs <- function(filenames) { 
  res <- GRangesList(lapply(filenames, .readExcelCNA))
  names(res) <- vapply(strsplit(filenames, "_"), `[`, 1,
                       FUN.VALUE=character(1))
  for (nm in names(res)) {
    mcols(res[[nm]])$name <- nm
  }
  sort(sortSeqlevels(unlist(res)))
}


# Helper function
.readExcelCNA <- function(filename) { 
  message("Processing ", filename, "...")
  res1 <- read_excel(filename, sheet=1) 
  res2 <- read_excel(filename, sheet=2)
  res <- rbind(res1, res2) 
  gr <- .parseCoords(res)
  gr$score <- res$log2Ratio
  return(gr)  
}


# Helper function
.parseCoords <- function(res) {
  gr <- as(sapply(strsplit(res$Altered_Region, "\\("), `[`, 1), "GRanges")
  seqlevelsStyle(gr) <- "UCSC"
  return(gr)
}


# Helper function
.makeRse <- function(gr) { 
  disjoinSummarizedExperiment(.makeRagExp(gr), mean, 1)
}


# Helper function
.makeRagExp <- function(gr, div=10) { 
  gr$score <- pmax(-2, pmin(2, round(gr$score*div)/div))
  grl <- split(gr, gr$name)
  re <- RaggedExperiment(grl)
  colData(re)$subject <- substr(colnames(re), 1, 9) 
  colData(re)$portion <- substr(colnames(re), 11, 11) 
  return(re)
}
