#!/usr/bin/bash

# Align multiome data using cellranger-arc

# Activate conda
eval "$(conda shell.bash hook)"
conda activate /work/Home/myenv # activate enviroment

# set path to cellranger-arc
export PATH=$PATH:/work/Home/programs/cellranger-arc/cellranger-arc-2.0.1:$PATH

# Define directory
dir_data=/work/Home/islets_multiome_hfd/data-raw/ATAC
dir_library=/work/Home/islets_multiome_hfd/data-raw/ATAC/cellranger_arc_library
dir_out=/work/Home/islets_multiome_hfd/data/cellranger_output
dir_genome=/work/Home/islets_multiome_hfd/data/assembly

cd $dir_library

ls *_library.csv | awk 'BEGIN{FS=OFS="_"}{print $1}'| sort | uniq | while read line
do
echo "***Start alignment on sample $line***"

# define library

lib=${line}_library.csv

echo "**** library is $lib****"

# Go to output directory
cd $dir_out

cellranger-arc count --id=${line} \
                       --reference=$dir_genome/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                       --libraries=$dir_library/$lib \
                       --localcores=16 \
                       --localmem=64

echo "***FINISHED alignment on sample $line***"
done
