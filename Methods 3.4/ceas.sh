#!/usr/bin/bash

# CEAS: Genomic annotation
eval "$(conda shell.bash hook)" #Run this in order for conda to work 
conda activate Python27
mkdir ceas
VARIANTS=("H1.0" "H1.2" "H1.4" "H1.5" "H1X")
for variants in "${VARIANTS[@]}";
do
        echo "------------CEAS: ${variants} --------------"
        ceas -b merged_islands/*${variants}* -g /home/ajvlab/Softwares/hg19_ceas.refGene
        mkdir ${variants}_ceas
        mv *${variants}*.R *${variants}*.xls *${variants}*.pdf ${variants}_ceas/ 
        mv ${variants}_ceas/ ceas/
done  

echo "CEAS is OVER!"
