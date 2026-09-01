# Regional Microglial Transcriptomics and Epigenomics

Alexander V. Margetts¹⁻³, Samara J. Vilca¹˒², Lauren L. Bystrom¹˒², Florence Bourgain-Guglielmetti¹˒², and Luis M. Tuesta¹⁻³*

¹ Department of Psychiatry & Behavioral Sciences
² Center for Therapeutic Innovation
³ Sylvester Comprehensive Cancer Center
University of Miami Miller School of Medicine, Miami, FL 33136

* Corresponding author: ltuesta@miami.edu
* Site maintainer: alexmargetts@med.miami.edu

Microglia, innate CNS immune cells in the brain, exhibit distinct transcriptional profiles across brain regions that are important for facilitating their specialized function. Identification of the epigenetic modifications associated with these distinct transcriptional profiles may improve our understanding of the mechanisms governing the functional specialization of microglia. However, the large number of microglia required to obtain a genome-wide profile for histone modifications has been an obstacle. To this end, we combined a method of chromatin profiling known as Cleavage Under Targets and Tagmentation with RNA sequencing to profile histone modifications associated with regional differences in gene expression in microglia. We find that region-dependent microglial transcriptomics are associated with differential deposition of H3K27ac and H3K27me3. Although the epigenetic landscape is remarkably stable across brain regions, H3K27ac may help explain their distinct transcriptional profiles. In this study, we report a molecular and computational framework to improve our understanding of epigenetic regulation in microglia using as few as 2,500 cells per histone mark.

[Manuscript]

# Navigating Our Repository

## Requirements

The analyses were developed using:

- R 4.2.1
- Bioconductor 3.16
- Ubuntu 22.04 on an x86_64 Linux system
- Git
- RStudio is recommended for opening the project and running the analysis files interactively

Exact R package versions, recursive dependencies, and GitHub package commit hashes are recorded in `renv.lock`. Package-version records for the individual analysis files are available in `code/RPackageVersions/`.

## Install Git and Clone the Repository

Git installers and platform-specific installation instructions are available from the [official Git website](https://git-scm.com/downloads).

On Ubuntu or Debian Linux, Git can be installed with:

```bash
sudo apt update
sudo apt install git
```

Confirm that Git is available:

```bash
git --version
```

Clone the repository and enter its root directory:

```bash
git clone https://github.com/avm27/RegionalMicroglialTranscriptomicsandEpigenomics.git
cd RegionalMicroglialTranscriptomicsandEpigenomics
```

## Obtain the Input Data

Processed input files needed for the principal R analyses are stored in `data/` when their size permits. Additional raw and processed files are available through the Gene Expression Omnibus:

- Complete dataset: [GSE277756](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277756)
- RNA-sequencing dataset: [GSE277755](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277755)
- CUT&Tag dataset: [GSE277753](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277753)

The CUT&Tag GEO dataset contains:

- Raw FASTQ files
- Peak-count matrices
- SEACR peak files
- Normalized BigWig files ending in `_norm.bw`

The RNA-sequencing GEO dataset contains:

- Raw FASTQ files
- The gene-count matrix
- Normalized-count files

BAM files are not provided for either CUT&Tag or RNA sequencing. Any analysis that requires BAM files must first recreate them by aligning the corresponding raw FASTQ files. RNA-sequencing BigWig files are also not provided and must be generated from recreated RNA-sequencing BAM files when genomic coverage tracks are required.

The raw and large processed sequencing files are not duplicated in this GitHub repository. Analyses that require these files are identified below.

## Restore the R Package Environment

Install R 4.2.1 before restoring the project library. Open `RegionalMicroglialTranscriptomicsandEpigenomics.Rproj` in RStudio, or start R from the root of the cloned repository.

Run:

```r
source("renv/activate.R")
renv::restore(prompt = FALSE)
```

If R cannot find `renv`, install it and then repeat the commands above:

```r
install.packages("renv", repos = "https://cloud.r-project.org")
```

`renv::restore()` creates an isolated project library containing the package versions recorded in `renv.lock`. It does not rerun any analyses and does not replace packages in the user's global R library.

Confirm the environment after restoration:

```r
R.version.string
BiocManager::version()
renv::status()
```

The expected versions are R 4.2.1 and Bioconductor 3.16. Run `renv::restore()` again after pulling a commit that changes `renv.lock`. Do not run `renv::update()` when the goal is to reproduce the original computational environment.

Some R packages require operating-system libraries. If `renv::restore()` reports missing system packages, install the packages it identifies and rerun `renv::restore()`. For example, on Ubuntu the restoration process may request:

```bash
sudo apt install gsfonts libmagick++-dev
```

## Run the Analysis Files

The files in `analysis/` are records of the analyses performed for this project. They are intended to be run interactively and selectively from RStudio. Rebuilding the workflowr website or rendering every R Markdown file is not required.

For each analysis file:

1. Open `RegionalMicroglialTranscriptomicsandEpigenomics.Rproj` in RStudio.
2. Confirm that `renv` is active and that `renv::status()` does not report unexpected inconsistencies.
3. Open the required `.Rmd` file from `analysis/`.
4. Restart the R session before starting a different analysis file.
5. Run the setup chunk first.
6. Run the required R chunks in their original order.

Do not use the **Knit** button or `workflowr::wflow_build()` unless intentionally rebuilding an entire analysis page. Some analyses require substantial data, memory, external software, and computation time.

The analysis files identify the repository root from the opened RStudio project and use repository-relative paths for files stored in `data/`, `output/`, and `code/`. These paths should not need to be edited after cloning the repository.

## Configure External Data

Analyses that require files too large to store in GitHub use the `BRMIC_EXTERNAL_DATA` environment variable. Set this variable to the absolute path of the directory containing the downloaded external data.

From R:

```r
Sys.setenv(BRMIC_EXTERNAL_DATA = "/absolute/path/to/BRMIC_external_data")
```

From a Linux or macOS terminal:

```bash
export BRMIC_EXTERNAL_DATA="/absolute/path/to/BRMIC_external_data"
```

Do not commit a computer-specific external-data path to the repository. The analyses that use this variable check for the required directories and files before proceeding.

The relevant external-data structure is:

```text
BRMIC_external_data/
├── fastq/
├── BAMs/
├── FASTQ/
├── bigwig_withnorm/
├── peakCalling/
│   └── SEACR/
│       ├── normalized/
│       └── HOMER/
└── reference/
    ├── mm10-blacklist.v2.bed
    └── mm10.chrom.sizes
```

Only the directories relevant to the analysis being run are required. The `BAMs/` directories shown above contain locally recreated BAM files; BAM files are not available from the GEO records.

## Recommended Analysis Order

The analysis files do not all need to be run consecutively. When using the processed inputs supplied with the repository or deposited on GEO, `index3.Rmd` and the preprocessing chunks in `index8.Rmd` and `index9.Rmd` can be skipped.

1. Restore the R environment with `renv::restore()`.
2. Run `index.Rmd` and `index2.Rmd` independently for the RNA-sequencing analyses.
3. Run `index3.Rmd` only when regenerating peak-count matrices from the deposited SEACR peak files and BAM files recreated from the raw CUT&Tag FASTQ files.
4. Run `index4.Rmd`, `index5.Rmd`, and `index6.Rmd` from their corresponding processed peak-count matrices.
5. Run `index7.Rmd` after the significant-result files have been generated by `index.Rmd`, `index4.Rmd`, and `index5.Rmd`.
6. Run the required plotting chunks in `index8.Rmd` after the necessary BigWig files are available.
7. Run `index9.Rmd` after the RNA-sequencing FPKM files and HOMER-annotated CUT&Tag peak files are available.

## Analysis Files, Inputs, and Starting Points

### `analysis/index.Rmd` — Differential Gene-Expression Analysis

**Required inputs:**

- `data/RNASeq_DESeq_Data/sample_info_mic_br.csv`
- `data/RNASeq_DESeq_Data/gene_count_matrix.csv`

**How to start:** Open `analysis/index.Rmd`, restart R, run the setup chunk, and then run the analysis chunks in order.

**Primary output:** Figure 1.

### `analysis/index2.Rmd` — Weighted Gene Co-Expression Network Analysis

**Required inputs:**

- `data/RNASeq_WGCNA/BRMIC_gene_count_matrix.csv`
- `data/RNASeq_WGCNA/BRMIC_Metadata.tsv`

**How to start:** Open `analysis/index2.Rmd`, restart R, run the setup chunk, and then run the analysis chunks in order. This analysis is independent of `index.Rmd`.

**Primary output:** Figure 2.

### `analysis/index3.Rmd` — Peak-Count Matrix Generation

**Required inputs:**

- Significant CUT&Tag SEACR peak files available through GEO
- CUT&Tag BAM files recreated by aligning the raw CUT&Tag FASTQ files available through GEO

**How to start:** Download the SEACR peak files and raw CUT&Tag FASTQ files from GEO. Recreate the required BAM files by running the CUT&Tag alignment and processing workflow, configure their external locations as documented in the analysis file, open `analysis/index3.Rmd`, restart R, and run the required chunks in order.

This analysis does not need to be run when using the deposited processed peak-count matrices.

**Primary outputs:**

- H3K27me3 peak-count matrix used by `analysis/index4.Rmd`
- H3K27ac peak-count matrix used by `analysis/index5.Rmd`
- Combined peak-count matrix used by `analysis/index6.Rmd`

### `analysis/index4.Rmd` — H3K27me3 Differential Peak-Deposition Analysis

**Required inputs:**

- `data/CATD_DESeq_Data/sample_info_K27me3.csv`
- `data/CATD_DESeq_Data/peak_count_matrix_CATD_Final_K27me3.csv`

**How to start:** Open `analysis/index4.Rmd`, restart R, run the setup chunk, and then run the analysis chunks in order. `index3.Rmd` is not required when the processed count matrix is present.

**Primary outputs:** Figures 3b, 3c, 3f, and 3g.

### `analysis/index5.Rmd` — H3K27ac Differential Peak-Deposition Analysis

**Required inputs:**

- `data/CATD_DESeq_Data/sample_info_K27ac.csv`
- `data/CATD_DESeq_Data/peak_count_matrix_CATD_Final_K27ac.csv`

**How to start:** Open `analysis/index5.Rmd`, restart R, run the setup chunk, and then run the analysis chunks in order. `index3.Rmd` is not required when the processed count matrix is present.

**Primary outputs:** Figures 3d, 3e, 3f, and 3h.

### `analysis/index6.Rmd` — Combined H3K27me3 and H3K27ac Analysis

**Required inputs:**

- `data/CATD_DESeq_Data/sample_info_combined.csv`
- `data/CATD_DESeq_Data/peak_count_matrix_CATD_Final.csv`

**How to start:** Open `analysis/index6.Rmd`, restart R, run the setup chunk, and then run the analysis chunks in order. This analysis can start directly from the combined processed count matrix.

**Primary output:** Figure 3a.

### `analysis/index7.Rmd` — Transcriptomic and Epigenomic Overlap Analysis

**Required inputs:**

- Significant RNA-sequencing comparison files in `output/RNASeq_DEGs/Results_Sig_CSVs/`
- Significant H3K27me3 comparison files in `output/CATD_DEPGs/H3K27me3/`
- Significant H3K27ac comparison files in `output/CATD_DEPGs/H3K27ac/`

These files are generated by `analysis/index.Rmd`, `analysis/index4.Rmd`, and `analysis/index5.Rmd`, respectively.

**How to start:** Confirm that the required significant-result files exist, open `analysis/index7.Rmd`, restart R, run the setup chunk, and then run the analysis chunks in order.

**Primary outputs:** Figures 4a and 4b.

### `analysis/index8.Rmd` — Genomic Track Plots

**Required inputs:**

- Normalized CUT&Tag BigWig files ending in `_norm.bw`, available through GEO
- RNA-sequencing BigWig files generated from recreated RNA-sequencing BAM files
- Raw RNA-sequencing FASTQ files from GEO when regenerating the RNA-sequencing BAM and BigWig files
- Raw CUT&Tag FASTQ files from GEO when regenerating the deposited CUT&Tag BigWig files
- [mm10 blacklist](https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz)
- [mm10 chromosome-size file](https://github.com/fl-yu/CUT-RUNTools-2.0/tree/master/assemblies/chrom.mm10)
- The mouse genome-annotation file used by the track-plotting code

**External software required only when regenerating BigWig files:**

- deepTools `bamCoverage`
- UCSC `bigWigMerge`
- UCSC `bedGraphToBigWig`

Neither the CUT&Tag nor RNA-sequencing BAM files are deposited on GEO. The deposited normalized CUT&Tag BigWig files can be used directly.

RNA-sequencing coverage tracks require the raw RNA-sequencing FASTQ files to be aligned to recreate the BAM files, followed by BigWig generation. Recreating the CUT&Tag BigWig files similarly requires alignment of the raw CUT&Tag FASTQ files before BigWig generation.

**How to start:** Download the deposited normalized CUT&Tag BigWig files. If RNA-sequencing coverage tracks are required, recreate the RNA-sequencing BAM and BigWig files from the deposited raw RNA-sequencing FASTQ files. Set `BRMIC_EXTERNAL_DATA`, open `analysis/index8.Rmd`, restart R, and run the setup and plotting chunks in order.

The shell preprocessing chunks have `eval=FALSE` and must be run manually only when regenerating BigWig files. They are not executed during an RStudio preview.

**Primary output:** Figure 4c.

### `analysis/index9.Rmd` — Regional Gene Expression and Histone-Mark Correlation

**Required inputs:**

- RNA-sequencing FPKM files in `data/RNASeq_FPKMs/`
- Significant CUT&Tag peak files generated by SEACR and available through GEO when regenerating the annotations
- HOMER-annotated H3K27ac and H3K27me3 peak files under `BRMIC_EXTERNAL_DATA/peakCalling/SEACR/HOMER/`

**External software and services:**

- HOMER `annotatePeaks.pl` when regenerating the annotation files
- Ensembl BioMart, accessed through the R package `biomaRt`, when performing the live gene-annotation query

**How to start:** Set `BRMIC_EXTERNAL_DATA`, open `analysis/index9.Rmd`, restart R, run the setup chunk, and then run the R chunks in order.

The HOMER shell chunk has `eval=FALSE` and should be run manually only when the annotation files need to be regenerated. Generation of the HOMER annotation files uses the deposited SEACR peak files and does not require BAM files.

An internet connection is required for the BioMart query unless the corresponding annotation results have been saved locally.

**Primary output:** Supplementary Figure 3.

## Session and Version Records

The package-version code in each analysis file writes the following records to `code/RPackageVersions/`:

- Package names and versions used by that analysis
- R and Bioconductor versions
- Complete `sessionInfo()` output

These files supplement `renv.lock`; the lockfile is the authoritative record used to restore the complete package environment.

This project was made with [workflowr]

[Manuscript]: https://www.biorxiv.org/content/10.1101/2024.08.08.607229v1
[workflowr]: https://github.com/workflowr/workflowr
