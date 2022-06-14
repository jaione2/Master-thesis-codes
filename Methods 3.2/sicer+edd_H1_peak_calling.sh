#!/usr/bin/bash

#------------------------------------------------------------------------------------------------------------------------------#
# EXPLANATION OF THE CODE 
# This code will take the BAM files of the alignments of the ChIP-Seq and do the peak calling.There will be 2 arguments: The
# path where the BED files are located (-p | --path) and the suffix that must be removed from the name of the file (e.g. ".bam")
# in order to obtain the "core" name (-s|--suffix). Peak calling will be done with SICER and EDD and then they will be merged.
#------------------------------------------------------------------------------------------------------------------------------#

# CHECK ARGUMENTS
INSTRUCTIONS="INSTRUCTIONS= $0 -p path -s suffix_to_remove"
if [ "$#" -ne "4" ]; then
	echo "Incorrect number of argument: $INSTRUCTIONS"
	exit 1
fi
case $1 in
	-p|--path)
	PATH_TO_FILE="$2"
	-s|--suffix)
	SUFFIX_TO_REMOVE="$2"
esac

case $3 in
	-p|--path)
	PATH_TO_FILE="$4"
	-s|--suffix)
	SUFFIX_TO_REMOVE="$4"
esac

#SICER ISLAND CALLING
cd $PATH_TO_FILE
for files in ChIP*.bam
do
	echo "---------- Peak calling: $files ----------"
	folder_name=${files//$SUFFIX_TO_REMOVE/}
	mkdir ${folder_name}
	cd $folder_name 
	sicer -t $PATH_TO_FILE/$files -c $PATH_TO_FILE/Input*bed -s hg19 -w 200 -rt 1 -f 150 -egf 0.80 -fdr 0.01 -g 600
	mv *FDR* ${folder_name}_sicer_island.bed
	cp *sicer*
	cd ..	
done
mkdir sicer_islands 
mv *sicer* sicer_islands

#EDD ISLANDS CALLING 
eval "$(conda shell.bash hook)" #Run this in order for conda to work 
conda activate Python27
 
for i in ChIP*.bam
do
        echo "EDD peak calling: $i"
        time edd /home/ajvlab/Softwares/hg19.chrom.sizes ../empty.bed $i Input*.bam ${i}_edd_output
done    

mkdir EDD_islands 
mv *edd* EDD_islands

#Merge SICER and EDD islands
mkdir merged_islands
VARIANTS=("H1.0" "H1.2" "H1.4" "H1.5" "H1X")
for variants in "${VARIANTS[@]}";
do
        sort -k1,1 -k2,2n sicer_islands/*${variants}* | sed 's/[ \t]\+$//' > tmp_file_1.bed
        sort -k1,1 -k2,2n EDD_islands/*${variants}* > tmp_file_2.bed
        cat tmp_file_1.bed tmp_file_2.bed | sort -k1,1 -k2,2n | bedtools merge > ${variants}_merged_island.bed
        mv *merged_island.bed merged_islands
done 

#------------------------------------------------------------------------------------------------------------------------------#
# The OUTPUT of the peak calling of a folder called "merged_islands" that have 5 files, 1 for each H1 variant. Moreover it will
# also generate "sicer_islands" and "EDD_islands" folder with the separated files. 
#------------------------------------------------------------------------------------------------------------------------------#
