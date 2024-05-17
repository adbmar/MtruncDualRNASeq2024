#!/usr/bin/env zsh
source ~/.zshrc

# This script implements pangenome analysis
folder_data="/Users/cychang/Desktop/lab/ncbi-ensifer/data"

# Move the annotated genomes to one folder
mkdir -p $folder_data/gffs
cp $folder_data/prokka/USDA1021/annotated.gff $folder_data/gffs/USDA1021.gff
cp $folder_data/prokka/WSM1022/annotated.gff $folder_data/gffs/WSM1022.gff

# Create a list of gffs
mkdir -p $folder_data/pangenome

for i in USDA1021 WSM1022
do
    echo -e $folder_data/gffs/$i.gff
done >| $folder_data/pangenome/list_gffs.txt

conda activate
mamba activate panaroo

panaroo \
    -i $folder_data/pangenome/list_gffs.txt \
    -o $folder_data/pangenome \
    --core_threshold 0.95 \
    -t 10 \
    --clean-mode strict \
    --remove-invalid-genes

