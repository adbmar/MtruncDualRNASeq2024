# MtruncDualRNASeq2024
Data and code repository for dual RNAseq project with Medicago truncatula

This repository holds scripts and data associated with the study found at ____.

The snakemake directory contains the script needed to download data from the Sequence Read Archive and the snakemake script to perform the bioinformatic analyses outlined in the publication.
- The adapter directory contains the adapter sequences needed for Trimmomatic
- The env_ files contain the enrivonments needed to run snakemake
- rawReads files should be downloaded from SRA from the following BioProject
- All other directories will be populated as the snakemake pipeline (snakemake.smk) is run and final data will be put into the counts directory
- The snakemake pipeline does not include making bowtie2 genome indices or HISAT2 genome indices. So indices must be made and put into the snakemake directory first. See or edit snakemake parameters for proper naming of indices.

The DGE directory contains all the scripts and data necessary to take the output of the snakemake pipeline, and perform the various differential expression and Gene Ontology overrepresentation analyses outlined in the publication. Within the DGE directory:
- samples.csv file contains sample information necessary to perform analyses in the PrimaryAnalysis.R script.
- Sudirectories (Medicago, Nematode, Rhizobe21, Rhizobe22) contain the count data (output of the snakemake pipeline and HTSeq-union method) and genomic annotations for each organism (see publication for detailed sources of these annotations).
- PrimaryAnalysis.R script is the first script that should be run and performs a number of differential expression analyses and outputs the results of these analyses into a large table in the out directory. See the publication supporting information for details about these analyses and the naming convention of columns in this file. The PrimaryAnalysis.R script also calls the Quality Check.R script which outputs various quality control metrics into the out directory as well.
- Final Figures Fig1.R, Final Figures Fig2.R, and SupFig1.R scripts must be run after the PrimaryAnalysis.R script has been run and will output figures in the out directory, many of which are found in the publication. These scripts rely on the Parially paired t test.R script.
- SupplementalDatasets.R script can be used to generate the datasets found in the supporting information of the publication and must be run only after PrimarilyAnalysis.R script has run.
- GO Enrichment Analysis.R and GO Enrichment Analysis with Parents.R are two scripts that run Gene Ontology enrichment analyses and output the results to the out directory. These scripts will create a EnrichmentDatabase directory that will house SQL databases needed to perform the analyses.
- Utility Scripts directory contains various utility scripts that may be used for data exploration. Many of these scripts require running the PrimaryAnalysis.R script first.
