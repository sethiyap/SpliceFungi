# Splice Fungi

## Splicing Analysis on fungal RNASeq data

RNASeq data was generated for Cn cells treated with drug. SpliceWiz,
R-package was implemented to process BAM files and further analyse the
splicing results.

### Install SpliceFungi

``` r
if(require("devtools")){
        options(repos = BiocManager::repositories())
        devtools::install_github("sethiyap/SpliceFungi",build = FALSE)
} else{
        options(repos = BiocManager::repositories())
        install.packages("devtools")
        devtools::install_github("sethiyap/SpliceFungi", build = FALSE)
}
```

#### Process BAM files

Fastq files were aligned to ensembl genome using STAR aligner and are
processed with SpliceWiz (as below).

``` r
ref_path <- "/Users/pooja/Documents/CDK7_project/SpliceFungi/Reference"


bam_path <- SpliceWiz::findBAMS(sample_path = "/Users/pooja/Documents/CDK7_project/RNASeq/RNASeq_5h/Mev_5h_ensembl/", level = 0)

pb_path=file.path("/Users/pooja/Documents/CDK7_project/SpliceFungi", "spliceWiz_Mev_5h")

SpliceWiz::processBAM(
    bamfiles = bam_path$path,
    sample_names = bam_path$sample,
    reference_path = ref_path,
    n_threads = 1,
    output_path = pb_path,
    run_featureCounts = T
)

# Load gene counts
gene_counts <- readRDS(file.path(pb_path, "main.FC.Rds"))

# Access gene counts:
gene_counts$counts
```

## Sample analysis

### WT+SY-1365 (1h)



    |                                                 |condition |batch |
    |:------------------------------------------------|:---------|:-----|
    |WT_DMSO_Set1_star_alignAligned.sortedByCoord.out |WT_1h     |set1  |
    |WT_DMSO_Set2_star_alignAligned.sortedByCoord.out |WT_1h     |set2  |
    |WT_DMSO_Set3_star_alignAligned.sortedByCoord.out |WT_1h     |set3  |
    |WT_Mevo_Set1_star_alignAligned.sortedByCoord.out |WT_Mev_1h |set1  |
    |WT_Mevo_Set2_star_alignAligned.sortedByCoord.out |WT_Mev_1h |set2  |
    |WT_Mevo_Set3_star_alignAligned.sortedByCoord.out |WT_Mev_1h |set3  |

    [1] "Only PSI plot"

![](README_files/figure-commonmark/mev_1h_splicing-1.png)

![](README_files/figure-commonmark/mev_1h_splicing-2.png)

![](README_files/figure-commonmark/mev_1h_splicing-3.png)

### ASE coverage plots

#### Intron Exclusion

``` r
event_list <- c("CNAG_03333/AFR96554_Intron1/clean", "CNAG_04112/AFR96843_Intron8/clean")

SpliceFungi::coverage_plot_for_list(event_list = event_list, res_edgeR = res_edgeR_Mev_1h, se_object = se_Mev_1h, control = "WT_1h",treatment = "WT_Mev_1h")
```

    [1] "CNAG_03333/AFR96554_Intron1/clean" "CNAG_04112/AFR96843_Intron8/clean"

![](README_files/figure-commonmark/ase_coverage_plot_IR%20exclude-1.png)

#### Intron retention

``` r
event_list <- c("CNAG_00399/AFR92532_Intron6/clean", "CNAG_06032/AFR98271_Intron4/clean")

SpliceFungi::coverage_plot_for_list(event_list = event_list, res_edgeR = res_edgeR_Mev_1h, se_object = se_Mev_1h, control = "WT_1h",treatment = "WT_Mev_1h")
```

    [1] "CNAG_00399/AFR92532_Intron6/clean" "CNAG_06032/AFR98271_Intron4/clean"

![](README_files/figure-commonmark/ase_coverage_plot_IR-1.png)

#### A3Ss

``` r
event_list <- c("A3SS:CNAG_00075-novelTr003-exon2;AFR92212-exon2", "A3SS:AFR95560-exon6;CNAG_07164-novelTr002-exon2")

SpliceFungi::coverage_plot_for_list(event_list = event_list, res_edgeR = res_edgeR_Mev_1h, se_object = se_Mev_1h, control = "WT_1h",treatment = "WT_Mev_1h")
```

    [1] "A3SS:CNAG_00075-novelTr003-exon2;AFR92212-exon2"
    [2] "A3SS:AFR95560-exon6;CNAG_07164-novelTr002-exon2"

![](README_files/figure-commonmark/ase_coverage_plot_A3S-1.png)
