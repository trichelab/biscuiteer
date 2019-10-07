### Example script for plotting mutation signatures
### Written by: Tim Triche, Jr.
### October 7, 2019
library(YAPSA) 
library(circlize)
library(matrixStats)
library(ComplexHeatmap) 
library(VariantAnnotation)

VCF_files <- list.files(pattern=".vcf$") 
names(VCF_files) <- sapply(strsplit(VCF_files, "_"), `[`, 1)
VCFs <- lapply(VCF_files, readVcf, genome="hg19")
VRL <- VRangesList(lapply(VCFs, as, "VRanges"))
for (nm in names(VRL)) VRL[[nm]]$pid <- nm

# helper fn
toYapsaDf <- function(VR) { 
  columns <- c("seqnames","start","ref","alt","pid")
  res <- as(unname(VR), "data.frame")[, columns]
  names(res) <- c("CHROM","POS","REF","ALT","PID")
  res$POS <- as.integer(res$POS)
  return(res)
}

data("exchange_colour_vector")
DF <- do.call(rbind, lapply(VRL, toYapsaDf))
DF <- DF[order(DF$PID, DF$CHROM, DF$POS),]
DF$change <- attribute_nucleotide_exchanges(DF) 
DF$col <- exchange_colour_vector[DF$change]
DF$SUBJECT <- substr(DF$PID, 1, 9)
DF$PORTION <- substr(DF$PID, 11, 11)

library(BSgenome.Hsapiens.UCSC.hg19)
mutCats <- 
  create_mutation_catalogue_from_df(DF,
                                    this_seqnames.field="CHROM", 
                                    this_start.field="POS",
                                    this_end.field="POS",
                                    this_PID.field="PID",
                                    this_subgroup.field="SUBJECT",
                                    this_refGenome=BSgenome.Hsapiens.UCSC.hg19,
                                    this_wordLength=3)

# use AlexCosmicValid_sig_df and AlexCosmicValid_sigInd_df
data(sigs)
cutoff <- 0.033 # 3.3% representation
mutcat_df <- as.data.frame(mutCats$matrix)
general_cutoff_vector <- rep(cutoff, dim(AlexCosmicValid_sig_df)[2])
cutoff_list <- LCD_complex_cutoff(in_mutation_catalogue_df=mutcat_df,
                                  in_signatures_df=AlexCosmicValid_sig_df,
                                  in_cutoff_vector=general_cutoff_vector,
                                  in_sig_ind_df=AlexCosmicValid_sigInd_df)
subgroups <- make_subgroups_df(DF, cutoff_list$exposures, 
                               in_subgroup.field="SUBJECT")
exposures_barplot(in_title="Mutational signatures", 
                  in_exposures_df=cutoff_list$norm_exposures,
                  in_signatures_ind_df=cutoff_list$out_sig_ind_df,
                  in_subgroups_df=subgroups)
dev.copy2pdf(file="mutationSignatures.NOSH_0003_and_0005.normalized.pdf") 

exposures_barplot(in_title="Mutational signatures", 
                  in_exposures_df=cutoff_list$exposures,
                  in_signatures_ind_df=cutoff_list$out_sig_ind_df,
                  in_subgroups_df=subgroups)
dev.copy2pdf(file="mutationSignatures.NOSH_0003_and_0005.absolute.pdf") 

# AC1 (green) is C->T, clocklike
# AC3 (goldenrod) is DNA repair defects (HRR)
# AC9 (brown) is POLE/SHM related
# AC12 (darkblue) is MMR related
# AC16 (violet) is endogenous.unknown
# AC19 (deep pink) is exogenous.unknown
# AC28 (navy blue) is (completely) unknown 
complex_heatmap_exposures(in_exposures_df=cutoff_list$norm_exposures, 
                          in_method="manhattan",
                          in_subgroups_df=subgroups,
                          in_subgroup_column="subgroup",
                          in_subgroup_colour_column="col",
                          in_signatures_ind_df=cutoff_list$out_sig_ind_df)
dev.copy2pdf(file="mutationSignatures.NOSH_0003_and_0005.clustered.pdf") 

plotExchangeSpectra(mutcat_df, in_show_triplets=TRUE) 
ggsave("mutationSignatures.NOSH_0003_and_0005.spectra.png")
