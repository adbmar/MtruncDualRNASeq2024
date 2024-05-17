#!/usr/bin/env zsh
source ~/.zshrc

# This script downloads the NCBI Ensifer genomes
conda activate
mamba activate ncbi-datasets
folder_data="/Users/cychang/Desktop/lab/ncbi-ensifer/data"

# Download the Ensifer genomes
datasets download genome accession \
    GCF_002197445.1 GCF_013315775.1 \
    --include genome,gff3,gbff \
    --filename $folder_data/ensifer_ncbi.zip

# EM1021 https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000006965.1/
# USDA1022 https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_004004435.1/

# Clean up the files and names
cd $folder_data
unzip ensifer_ncbi.zip

mkdir -p $folder_data/genomes/
cp $folder_data/ncbi_dataset/data/GCF_002197445.1/*_genomic.fna $folder_data/genomes/USDA1021.fasta
cp $folder_data/ncbi_dataset/data/GCF_013315775.1/*_genomic.fna $folder_data/genomes/WSM1022.fasta

# mkdir -p $folder_data/gffs/
# cp $folder_data/ncbi_dataset/data/GCF_002197445.1/genomic.gff $folder_data/gffs/USDA1021.gff
# cp $folder_data/ncbi_dataset/data/GCF_013315775.1/genomic.gff $folder_data/gffs/WSM1022.gff

rm -rf README.md ensifer_ncbi.zip ncbi_dataset

# USDA1021 https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002197445.1/
# WSM1022 https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_013315775.1/
