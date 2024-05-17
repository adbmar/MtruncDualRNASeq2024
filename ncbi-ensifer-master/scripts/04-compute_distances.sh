#!/usr/bin/env zsh
source ~/.zshrc

# This script implements pangenome analysis
folder_data="/Users/cychang/Desktop/lab/ncbi-ensifer/data"

# 1. ANI
mkdir -p $folder_data/ani

## Create a list of genome fasta files
for i in USDA1021 WSM1022; do;
    echo -e $folder_data/genomes/$i.fasta
done >| $folder_data/ani/list_genomes.txt

## Compute ani
mamba activate fastani
fastANI -t 10 \
    --ql $folder_data/ani/list_genomes.txt \
    --rl $folder_data/ani/list_genomes.txt \
    -o $folder_data/ani/ani_genomes.txt
# `--ql` list of names of query sequences in fasta
# `--rl` list of names of reference sequences in fasta


# 2. kmers
mkdir -p $folder_data/kmer

## Create kmer signatures
mamba activate sourmash
for i in USDA1021 WSM1022; do;
    sourmash sketch dna -p k=31,scaled=1000 \
        -o $folder_data/kmer/$i.sig \
        $folder_data/genomes/$i.fasta
done

## Create a list of kmer signatures
for i in USDA1021 WSM1022; do;
    echo $folder_data/kmer/$i.sig
done |> $folder_data/kmer/list_sigs.txt

## Compare signatures
mamba activate sourmash
sourmash compare \
    --from-file $folder_data/kmer/list_sigs.txt \
    --ksize 31 \
    --distance-matrix \
    --csv $folder_data/kmer/kmer.txt




