#!/usr/bin/bash

eval "$(conda shell.bash hook)"
conda activate /work/Home/myenv # activat enviroment

# Specify directory ---------------------------------------
mkdir -p /work/Home/new_alignment/genome

dir_genome="/work/Home/islets_multiome_hfd/data/genome"

# Genome assembly files from cellranger-arc
dir_fasta="/work/Home/islets_multiome/data/assembly/refdata-cellranger-arc-mm10-2020-A-2.0.0/fasta/genome.fa"
dir_genes="/work/Home/islets_multiome/data/assembly/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf"
# Generate Genome -----------------------------------------
STAR --runMode genomeGenerate --runThreadN 60 --genomeDir $dir_genome --genomeFastaFiles $dir_fasta --genomeSAindexNbases 14 --genomeChrBinNbits 18 --genomeSAsparseD 3 --limitGenomeGenerateRAM 17179869184 --sjdbGTFfile $dir_genes
