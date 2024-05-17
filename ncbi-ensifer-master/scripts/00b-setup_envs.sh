#!/usr/bin/env zsh

# This script creates the mamba envs for various bioconda tools
# In priniciple, one env per tool

cd
mkdir -p ~/bioinformatics

# Install ncbi-datasets v15.27.1
mamba create -y -n ncbi-datasets
mamba activate ncbi-datasets
mamba install --yes -c bioconda ncbi-datasets=15.27.1

# Install prokka v1.14.5
mamba create -n prokka
mamba activate prokka
mamba install --yes -c bioconda prokka=1.14.5

# Install fastANI v1.31
# fastANI is developed for fast alignment-free computation of whole-genome Average Nucleotide Identity (ANI)
mamba create -n fastani
mamba activate fastani
mamba install --yes -c bioconda fastani=1.31

# Install panaroo v1.3.4
# panaroo is A Bacterial Pangenome Analysis Pipeline that can call large structural variants
mamba create -n panaroo
mamba activate panaroo
mamba install --yes -c bioconda panaroo=1.3.4
