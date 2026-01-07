# Regional Microglial Transcriptomics and Epigenomics

Alexander V. Margetts1-3, Samara J. Vilca1,2, Florence Bourgain-Guglielmetti1-3 & Luis M. Tuesta1-3*

1 Department of Psychiatry & Behavioral Sciences
2 Center for Therapeutic Innovation
3 Sylvester Comprehensive Cancer Center
  University of Miami Miller School of Medicine, Miami, FL 33136
* Corresponding author (ltuesta@miami.edu)
* Site Maintainer (alexmargetts@med.miami.edu)

Microglia, the innate immune cells in the central nervous system, exhibit distinct transcriptional profiles
across brain regions that are important for facilitating their specialized function. There has been recent
interest in identifying the epigenetic modifications associated with these distinct transcriptional profiles,
as these may improve our understanding of the underlying mechanisms governing the functional
specialization of microglia. One obstacle to achieving this goal is the large number of microglia required
to obtain a genome-wide profile for a single histone modification. Given the cellular and regional
heterogeneity of the brain, this would require pooling many samples which would impede biological
applications that are limited by numbers of available animals. To overcome this obstacle, we have
adapted a method of chromatin profiling known as Cleavage Under Targets and Tagmentation
(CUT&Tag-Direct) to profile histone modifications associated with regional differences in gene
expression throughout the brain reward system. Consistent with previous studies, we find that
transcriptional profiles of microglia vary by brain region. However, here we report that these regional
differences also exhibit transcriptional network signatures specific to each region. Additionally, we find
that these region-dependent network signatures are associated with differential deposition of H3K27ac
and H3K7me3, and while the H3K27me3 landscape is remarkably stable across brain regions, the
H3K27ac landscape is most consistent with the anatomical location of microglia which explain their
distinct transcriptional profiles. Altogether, these findings underscore the established role of H3K27me3
in cell fate determination and support the active role of H3K27ac in the dynamic regulation of microglial
gene expression. In this study, we report a molecular and computational framework that can be applied
to improve our understanding of the role of epigenetic regulation in microglia in both health and disease,
using as few as 2,500 cells per histone mark.

[Manuscript]


# Navigating our Repository:

Overall data files can be found in: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277756
Necessary RNA Sequencing data files can be found in : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277755
Necessary CUT&Tag Sequencing data files can be found in: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277753


## RNA Sequencing Data:

### Differential Gene Expression Analysis

analysis/index.Rmd

Necessary input files: 

1. data/RNASeq_DESeq_Data/sample_info_mic_br.csv
2. data/RNASeq_DESeq_Data/gene_count_matrix.csv

Outputs:
1. Figure 1

### Weighted Gene Correlation Network Analysis

analysis/index2.Rmd

Necessary input files:

1. data/RNASeq_WCGNA/BRMIC_gene_count_matrix.csv
2. data/RNASeq_WCGNA/BRMIC_Metadata.tsv

Outputs:
1. Figure 2

### DEPG Analyses:

#### Generating Peak Count Matrices using Significant Peak lists and Associated BAM files -- Not necessary to run

analysis/index3.Rmd 

Necessary Input Files:
1. GEO Peak Files (CATD)
2. GEO BAM files (CATD)

#### H3K27me3 Differential Peak Deposition Analysis

analysis/index4.Rmd

Necessary Input Files:
1. data/CATD_DESeq_Data/sample_info_K27me3.csv
2. data/CATD_DESeq_Data/peak_count_matrix_CATD_Final_K27me3.csv

Outputs:
1. Figure 3b
2. Figure 3c
3. Figure 3f
4. Figure 3g

#### H3K27ac Differential Peak Deposition Analysis

analysis/index5.Rmd

Necessary Input Files:
1. data/CATD_DESeq_Data/sample_info_K27ac.csv
2. data/CATD_DESeq_Data/peak_count_matrix_CATD_Final_K27ac.csv

Outputs:
1. Figure 3d
2. Figure 3e
3. Figure 3f
4. Figure 3h

#### Combined H3K27me3 and H3K27ac Differential Peak Deposition Analysis

analysis/index6.Rmd

Necessary Input Files:
1. data/CATD_DESeq_Data/sample_info_combined.csv
2. data/CATD_DESeq_Data/peak_count_matrix_CATD_Final.csv

Outputs:

1. Figure 3a

#### Region specific epigenetic regulation of microglial genes -- overlapping significant results

analysis/index7.Rmd

Necessary Input Files:
1. output/RNASeq_DEGs/Results_Sig_CSVs/***All comparisons****.csv (From index.Rmd)
2. output/CATD_DEPGs/H3K27me3/***All comparisons****.csv (From index4.Rmd)
3. output/CATD_DEPGs/H3K27ac/***All comparisons****.csv (From index5.Rmd)

Outputs:

1. Figure 4a
2. Figure 4b


#### Region specific epigenetic regulation of microglial genes -- Generating IGV Views

analysis/index8.Rmd

Necessary Input Files:
1. https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz (Download and save to folder)
2. https://github.com/fl-yu/CUT-RUNTools-2.0/tree/master/assemblies/chrom.mm10 (Download and save to folder)
3. Need "begGraphtoBigWig" and "bigWigMerge" programs from UCSC Tools
4. All RNA-Seq BAM Files from GEO or proccess FASTQ to BAM files using our script (code/BRMic_RNA_Seq_PreProcessingCodeFinal)
5. All CATD BAM Files from GEO or proccess FASTQ to BAM files using our script (code/BRMic_CATD_Seq_PreProcessingCodeFinal)

Outputs:

1. Figure 4c

#### Region specific epigenetic regulation of microglial genes -- correlating regional gene expression

analysis/index9.Rmd

Necessary Input Files:
1. All FPKM Files from RNA-Seq, which are generated during the StringTie step of the RNA-Seq processing pipeline.  (data/RNASeq_FPKMs)
2. Significantly called peak files (output from SEACR) (Download from GEO)
3. HOMER Program (Specifcally: annotatePeaks.pl)

Ouputs:

1. Supplementary Figure 3



[Manuscript]: https://www.biorxiv.org/content/10.1101/2024.08.08.607229v1



[workflowr]: https://github.com/workflowr/workflowr
