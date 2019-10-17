# External data (extdata) source and preprocessing

The data found in the extdata directory are derived from the paper,
[The RON Receptor Tyrosine Kinase Promotes Metastasis by Triggering
MBD4-Dependent DNA Methylation Reprogramming](https://www.ncbi.nlm.nih.gov/pubmed/24388747).
The associated data has been uploaded to the Gene Expression Omnibus (GEO) with
the GEO Accession number [GSM127125](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1274125).

## MCF7_Cunha_header_only.vcf.gz(.tbi) and MCF7_Cunha_chr11p15.bed.gz(.tbi)

WGBS data accessed from the GEO link above was aligned to the hg19 genome using
[biscuit](https://github.com/zwdzwd/biscuit) version 1.3.20160324. The output
VCF was converted to a BED file (using biscuit vcf2bed) and then chromosome 
11p15 (chr11p15) was selected as the data is purely for example purposes.
Following the reduction to chr11p15, the BED file was gzipped and tabixed using
bgzip and tabix, which are available through samtools. Due to the large size of
the corresponding VCF file, the header of the file was selected out using the
command line tool `head` (`head -n 125 vcf_file.vcf`). This reduced file was 
then gzipped and tabixed in the same way as the BED file.

## MCF7_Cunha_shuffled_header_only.vcf.gz(.tbi) and MCF7_Cunha_chr11p15_shuffled.bed.gz(.tbi)

The "shuffled" versions of the MCF7 Cunha BED and VCF files are derived from the
corresponding BED and VCF files described previously. This was done as these 
files are merely to produce examples of biscuiteer capabilities and are not to
be used in a normal analysis. With that in mind, the VCF file was modified from
the original by changing out references to "MCF7_Cunha" with
"MCF7_Cunha_shuffled". The BED file was generated with the code found in
`scripts/genShuffledMCF7CunhaBed.py`. Please see the code for documentation on
how the shuffled BED file was created. Once the shuffled VCF and BED files were
generated, they were gzipped and tabixed for use in biscuiteer.
