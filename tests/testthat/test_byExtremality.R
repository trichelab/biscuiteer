test_that("test='byExtremality' gives correct errors", {
    
    orig_bed <- system.file("extdata", "MCF7_Cunha_chr11p15.bed.gz",
                            package="biscuiteer")
    orig_vcf <- system.file("extdata",
                            "MCF7_Cunha_header_only.vcf.gz",
                            package="biscuiteer")
    bisc2 <- readBiscuit(BEDfile = orig_bed, VCFfile = orig_vcf,
                         merged = FALSE)
    expect_error(byExtremality(comb))
})
