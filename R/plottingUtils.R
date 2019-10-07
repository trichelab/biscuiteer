#' Plot variant annotation matrix
#'
#' @param VRL      VRangesList of variants
#' @param cutoff   Minimum value to plot (DEFAULT: 0.25)
#' @param ANNonly  Show annotations only (DEFAULT: FALSE)
#'
#' @return         A Heatmap of the variants
#'
#' @importFrom circlize colorRamp2
#' @import ComplexHeatmap
#' @import VariantAnnotation
#'
#' @export
#'
plotVafMat <- function(VRL,
                       cutoff = 0.25,
                       ANNonly = FALSE) { 

  vafMat <- .getVafMat(VRL, cutoff=cutoff, ANNonly=ANNonly)

  if (ANNonly) {
    rowSym <- vapply(strsplit(rownames(vafMat), "\\|"), `[`, 1,
                     FUN.VALUE=character(1))
    vafMat <- rowsum(vafMat, rowSym)
  }

  cols <- colorRamp2(c(      0,   0.5,      0.75,       1),
                     c('white', 'red', 'darkred', 'black'))
  subj <- substr(colnames(vafMat), 1, 9)
  Heatmap(vafMat, name="VAF", show_row_names=ANNonly, 
          clustering_distance_rows="manhattan",
          clustering_method_rows="ward.D2",
          col=cols, column_split=subj)
}

#' Plot variant annotation matrix and CN
#'
#' @param VRL       VRangesList of variants
#' @param rse       RaggedExperiment of CNs
#' @param cutoff    Minimum value to plot for variants (DEFAULT: 0.25)
#' @param minCN     Minimum CN to plot (DEFAULT: 0.25)
#' @param minWidth  Scaling factor for CNs (DEFAULT: 1000)
#'
#' @return          A Heatmap of the variants and CNs
#'
#' @importFrom grid unit
#' @importFrom stats relevel
#' @importFrom circlize colorRamp2
#' @import ComplexHeatmap
#' @import VariantAnnotation
#'
#' @export
#'
plotVafAndCN <- function(VRL, rse, cutoff=0.25, minCN=0.25, minWidth=1e3) { 

  bothMat <- .getVafAndCN(VRL, rse, cutoff=cutoff,
                         minCN=minCN, minWidth=minWidth)

  chr <- as.factor(.getChr(bothMat)) 
  for (ch in rev(c(paste0("chr", c(seq_len(22), "X", "Y"))))) {
    if (ch %in% levels(chr)) {
      chr <- relevel(chr, ch)
    }
  }

  cols <- colorRamp2(c(        -1,  -0.75,        0,  0.75,         1), 
                     c('darkblue', 'blue', 'gray99', 'red', 'darkred'))

  subj <- substr(rownames(bothMat), 1, 9)
  annoSyms <- .annotateSymbols(bothMat)

  Heatmap(bothMat, name="VAF/CN", 
          border="gray50", col=cols,
          column_title_rot=90, 
          show_column_names=FALSE,
          show_column_dend=FALSE,
          cluster_column_slices=FALSE,
          column_gap=unit(2, "mm"),
          row_gap=unit(2, "mm"),
          show_parent_dend_line=FALSE,
          clustering_distance_columns="manhattan",
          clustering_method_columns="ward.D2", 
          clustering_distance_rows="manhattan",
          clustering_method_rows="ward.D2",
          row_split=subj, column_split=chr,
          bottom_annotation=annoSyms)
}

# Helper function
.getChr <- function(bothMat) {
  toSplit <- grep("\\|", colnames(bothMat))
  coords <- colnames(bothMat) 
  coords[toSplit] <- vapply(strsplit(coords[toSplit], "\\|"), `[`, 2,
                            FUN.VALUE=character(1)) 
  vapply(strsplit(coords, ":"), `[`, 1,
         FUN.VALUE=character(1)) 
}

# Helper function
#' @importFrom grid unit gpar
.annotateSymbols <- function(bothMat, annoFontSize=9) {
  toAnnotate <- grep("\\|", colnames(bothMat))
  geneLabels <- vapply(strsplit(colnames(bothMat)[toAnnotate], "\\|"), `[`, 1,
                       FUN.VALUE=character(1))
  columnAnnotation(link=anno_mark(at=toAnnotate, 
                                  side="bottom",
                                  labels=geneLabels, 
                                  extend=unit(2, "cm"),
                                  link_height=unit(1.5, "cm"),
                                  labels_gp=gpar(fontsize=annoFontSize)),
                   height=unit(3, "cm"))
}

# Helper function
.getVafAndCN <- function(VRL, rse, cutoff=0.25, minCN=0.25, minWidth=1e3) {
  vafMat <- .getVafMat(VRL, cutoff=cutoff, ANNonly=TRUE)
  cnaMat <- .getCnaMat(rse, cutoff=minCN, minWidth=minWidth)/2 
  cnaMat[ which(is.na(cnaMat)) ] <- 0
  stopifnot(identical(colnames(vafMat), colnames(cnaMat)))
  cnaMat <- cnaMat/2 # scale to make sense alongside VAF
  bothMat <- t(rbind(cnaMat, vafMat))
  return(bothMat)
}

# Helper function 
#' @importFrom matrixStats rowMaxs
.getVafMat <- function(VRL, cutoff=0.25, ANNonly=FALSE) { 
  if (ANNonly) VRL <- VRangesList(lapply(VRL, .annOnly))
  VRL <- VRangesList(lapply(VRL, .addLocation))
  vafMat <- .asMatrix(VRL, ANNonly=ANNonly)
  vafMat <- subset(vafMat, rowMaxs(vafMat) >= cutoff)
  return(vafMat)
}

# Helper function 
#' @importFrom matrixStats rowMaxs
.getCnaMat <- function(rse, cutoff=0.25, minWidth=1e3) {
  wideEnough <- width(rse) >= minWidth
  deepEnough <- rowMaxs(abs(assays(rse)$score), na.rm=TRUE) >= cutoff
  assays(subset(rse, deepEnough & wideEnough))[[1]]
}

# Helper function
.asMatrix <- function(VRL, ANNonly=FALSE) {
  VRL <- VRangesList(.addVaf(VRL))
  if (ANNonly) VRL <- VRangesList(lapply(VRL, .addLocationToName))
  variants <- sort(Reduce(union, lapply(VRL, names)))
  bigMat <- matrix(0, nrow=length(variants), ncol=length(VRL))
  rownames(bigMat) <- variants
  colnames(bigMat) <- names(VRL)
  for (i in names(VRL)) bigMat[names(VRL[[i]]), i] <- VRL[[i]]$VAF
  return(bigMat)
}

# Helper function
.addVaf <- function(x) { 
  if (is(x, "VRangesList")) {
    lapply(x, .addVaf)
  } else {
    x$VAF <- altDepth(x) / totalDepth(x)
    return(x)
  }
}

# Helper function
.addLocationToName <- function(vr) { 
  names(vr) <- paste(names(vr), vr$location, sep="|")
  return(vr)
}

# Helper function
.annOnly <- function(vr) {
  vr <- subset(vr, sapply(grepl("(HIGH|MODERATE)", vr$ANN), any))
  names(vr) <- .getSymbol(vr)
  return(vr)
}

# Helper function
.getSymbol <- function(vr) { 
  sapply(vr$ANN, 
         function(x) as.character(sapply(strsplit(x, "\\|"), `[`, 5)[1]))
} 

# Helper function
.addLocation <- function(vr) {
  vr$location <- as.character(vr) # for later placement 
  return(vr)
}
