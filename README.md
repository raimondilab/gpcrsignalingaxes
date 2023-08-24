# GPCR Ligand Signaling Dysregulation in Cancer Transcriptomics

Welcome to the repository for the manuscript titled "Exploring Dysregulation of GPCR Ligand Signaling Systems in Cancer Transcriptomics." This repository houses the data, code, and resources used in the study to investigate potential therapeutic opportunities in oncology through the analysis of G Protein-Coupled Receptor (GPCR) ligand signaling dysregulation in cancer transcriptomics datasets.

![cancer_cell_2.svg](https://github.com/raimondilab/gpcrsignalingaxes/blob/main/cancer_cell_2.svg)

## Repository Contents

data/: This directory contains the relevant datasets used in the study.

scripts/: Here, you will find the scripts and code used for data preprocessing, analysis, and visualization.

Please feel free to explore the repository and utilize the provided resources. If you have any questions or require further clarification, don't hesitate to reach out to the authors.

## Following is a brief explanation of the essential scripts:

*1. deseq2_whole.R* : the script performs comprehensive analysis on transcriptomics data from different cancer types, including data preprocessing, differential expression analysis, annotation of results, and generation of various types of visualizations. Additionally, there's a placeholder section for running GSEA, which sis designed to uncover enriched gene sets associated with the identified differential expression patterns in cancer transcriptomics data.

Usage:

Make sure you have R and the required packages (like DESeq2, ggplot2, org.Hs.eg.db, WebGestaltR) installed.
Place this script in your project directory or repository.
Update the paths, filenames, and any other necessary parameters according to your data and requirements.
Run the script using the command line or your preferred R environment.
Please note that this script assumes the existence of certain files and directories, so ensure that your directory structure matches the assumptions made in the script. 

here are the input files that are needed (present in '/data/'):

cancer_list.csv: A CSV file containing a list of cancer types. This file is read to determine which cancer types to process.

/data/TCGA_GTEX_meta_count_merged.rds: A merged dataset in RDS format that contains metadata and count data from TCGA and GTEX samples (original source: https://xenabrowser.net/datapages/?cohort=TCGA%20TARGET%20GTEx&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443).

probemap.csv: A CSV file containing a mapping between gene IDs and gene symbols.

Output: The results are present in '/data/TCGAvGTEX_CN'

*2. deseq2_subtypes.R* : This script is similar to the script above, except that it performs the same operations for the subtypes of a given cancer type. The subtype information is retrieved via the TCGAbiolinks library. 

Usage of the script:
To use the script, you would typically run it from the command line and provide a specific cancer type as a command-line argument. For example:

Rscript script_name.R Lung

This instructs the script to analyze data for the "Lung" cancer type. The script reads the data, performs DE analysis, generates plots, and conducts GSEA for each subtype associated with the given cancer type.


Input:

args: Command-line argument specifying the cancer type to analyze.

Merged dataset (/data/TCGA_GTEX_meta_count_merged.rds): Contains metadata and count data from TCGA and GTEX samples.

probemap.csv: A mapping file between gene IDs and gene symbols.

PanCancerAtlas_subtypes(): Retrieves subtype information for different cancer types using the TCGAbiolinks package.

Output:

Subdirectories for each cancer type (e.g., /data/subtypes_TvG_all/Lung/).

Various CSV files containing DE analysis results for different subtypes of the specified cancer type.

PNG images of count comparison plots, volcano plots, and PCA plots for each subtype.

GSEA analysis results in an output directory for each subtype.

The script generates output files and images that provide insights into gene expression differences, subtype-specific patterns, and enriched gene sets related to the provided cancer type.

*3. survival_whole.R* : This script performs survival analysis for different cancer types by analyzing the relationships between receptor-ligand pairs (or receptor-enzyme pairs) and survival outcomes. It calculates hazard ratios, p-values, confidence intervals, and other statistics for different scenarios involving the receptor, ligand, and their combined effect. The script iterates over a list of cancer types and processes survival data for each type. The results of the analysis are saved in CSV files.

Usage of the script:

The script is intended to be run in an R environment. It reads receptor-ligand (or receptor-enzyme) pairs from a CSV file, and for each cancer type, it reads survival data from separate CSV files. It performs Cox proportional hazards survival analysis on the provided data to analyze the impact of receptor and ligand expression levels on survival outcomes. The results are stored in an output CSV file for each cancer type.


Input:

Receptor-ligand (or receptor-enzyme) pairs CSV file (e.g., 'RL_pairs.csv' or 'RE_pairs.csv'): This file contains information about the pairs of genes being studied.
Individual CSV files for each cancer type (e.g., 'Bile duct.csv', 'Adrenal gland.csv', etc.): These files contain survival data as well as expression data specific to each cancer type.

Output:

Individual CSV output files for each cancer type (e.g., 'Bile duct.csv', 'Adrenal gland.csv', etc.): These files contain the results of the survival analysis for each receptor-ligand (or receptor-enzyme) pair, including hazard ratios, p-values, confidence intervals, and other statistics. Each output file contains rows representing different pairs and their analysis results.

*survival_subtypes.R*: This script performs the same operations as the earlier script but for the subtypes. The input files for this script are:

Receptor-Ligand or Receptor-Enzyme Pairs CSV File:

The script loads receptor-ligand or receptor-enzyme pairs information from a CSV file. The specific file used is determined by the lines:

#receptor_enzyme file
#rl_list <- read.csv('/data/survival/RE_pairs.csv')

#receptor_ligand file
rl_list <- read.csv('/data/survival/RL_pairs.csv')

You can choose to have either 'RE_pairs.csv' or 'RL_pairs.csv' in the /data/survival/ directory, depending on which pairs you are studying (receptor-enzyme or receptor-ligand).

Cancer Subtype List File:

The script reads a list of cancer subtypes from a CSV file. The specific file used is determined by the lines:

can<-scan('/data/survival/cancer_files.csv', what="", sep="\n")

There's a CSV file named 'cancer_files.csv' in the /data/survival/ directory. This file contains a list of cancer subtype filenames, each on a separate line.

Cancer Subtype Data Files:

The script iterates through each cancer subtype in the list obtained from the 'cancer_files.csv' file and reads individual data CSV files. The data files are named after the cancer subtype and are located in the /data/survival/data_subtypes/ directory. The specific lines for reading data files are:

data1<-read.csv('/data/survival/data_subtypes/'+cancer)


The script processes these input files to perform survival analysis for each cancer subtype and receptor-ligand/enzyme pair, and then generates output CSV files containing the analysis results for each cancer subtype.




