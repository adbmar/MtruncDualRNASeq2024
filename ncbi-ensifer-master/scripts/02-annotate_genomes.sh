# #!/usr/bin/env zsh
source ~/.zshrc

# This annotates the genome
conda activate
mamba activate prokka
folder_data="/Users/cychang/Desktop/lab/ncbi-ensifer/data"

# Annotate genoms
mkdir $folder_data/prokka

for i in USDA1021 WSM1022
do
    prokka --force --outdir $folder_data/prokka/$i --kingdom Bacteria --locustag $i --prefix annotated --gcode 11 $folder_data/genomes/$i.fasta
    # `--force` force overwriting existing output folder
    # `--outdir` output folder
    # `--kingdom`
    # `--locustag` locus tax prefix
    # `--prefix` filname output prefix
    # `--gcode` genetic code / translation table (set if --kingdom is set)
done

