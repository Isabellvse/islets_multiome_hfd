#!/usr/bin/bash
## Alignment of snMultiome RNA data

eval "$(conda shell.bash hook)"
conda activate /work/Home/myenv # activat enviroment

# STARsolo manual: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
# Define directory
dir_assembly=/work/Home/islets_multiome/data/assembly
dir_genome=/work/Home/new_alignment/genome/
dir_data=/work/Home/new_alignment/RNA
dir_out=/work/Home/new_alignment/star_output

mkdir -p $dir_out

# go to data folder
cd $dir_data

# Load STAR genome into memory
STAR --genomeDir "$dir_genome" --genomeLoad LoadAndExit
rm Log* Aligned.out.sam SJ.out.tab

# Align RNA data Genefull.
ls *.fastq.gz | awk 'BEGIN{FS=OFS="_"}{print $1,$2}'| sort | uniq | while read line
do
echo "***Start alignment on sample $line***"

# Specify reads and lanes. R1 = barcode and UMI reads R2 = cDNA reads
R1_L1=${line}_L001_R1_001.fastq.gz
R1_L2=${line}_L002_R1_001.fastq.gz
R2_L1=${line}_L001_R2_001.fastq.gz
R2_L2=${line}_L002_R2_001.fastq.gz

echo "***Read 1 Lane 1 and 2 = $R1_L1 and $R1_L2 . Read 2 Lane 1 and 2 = $R2_L1 and $R2_L2***"
STAR \
--genomeLoad LoadAndKeep \
--runThreadN 60 \
--genomeDir $dir_genome \
--readFilesCommand zcat \
--readFilesIn $R2_L1,$R2_L2 $R1_L1,$R1_L2 \
--soloUMIlen 12 \
--soloCellFilter None \
--soloType CB_UMI_Simple \
--soloCBmatchWLtype 1MM \
--clipAdapterType CellRanger4 \
--soloCBwhitelist $dir_assembly/737K-arc-v1.txt \
--soloFeatures GeneFull_Ex50pAS Gene

mkdir -p "$dir_out/$line/"
mv Solo.out Log* "$dir_out/$line/"
rm Aligned.out.sam SJ.out.tab
    
echo "***FINISH alignment of ${line}***"
done

STAR --genomeDir "$dir_genome" --genomeLoad Remove
rm Log* Aligned.out.sam SJ.out.tab 
