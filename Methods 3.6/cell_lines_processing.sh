#!usr/bin/bash
#------------------------------------------------------------------------------------------------------------------------------#
# This script will take ChIPSeq data for H1 variants of different cell lines and process in order to obtain the mapping of the # H1 signal to the UCSC LADs, merged islands (SICER + EDD), general summarized information of the islands (mean and median
# length, number of bases, number of islands) and CEAS genomic annotation. The input will be 2 folders: BDG_files and
# BAM_files. There will be 2 arguments: The path in which the "BDG_files" and "BAM_files" are located and the cell line we are
# analyzing.
#------------------------------------------------------------------------------------------------------------------------------#
# CHECK ARGUMENTS
INSTRUCTIONS="INSTRUCTIONS= $0 -p path -c cell_line"
if [ "$#" -ne "4" ]; then
	echo "Incorrect number of argument: $INSTRUCTIONS"
	exit 1
fi
case $1 in
	-p|--path)
	PATH_TO_FILE="$2"
	-c|--cell_line)
	CELL_LINE="$2"
esac

case $3 in
	-p|--path)
	PATH_TO_FILE="$4"
	-c|--cell_line)
	CELL_LINE="$4"
esac

#------------------------------------------------------------------------------------------------------------------------------#

# 1) Peak calling with SICER and EDD 
### SICER ISLANDS 
cd $PATH_TO_FILE
cd BAM_files
mkdir H1_islands
SUFFIX_TO_REMOVE=.bam 
for files in *H1*.bam
do
        echo "---------- SICER Peak calling: $files ----------"
        folder_name=${files//$SUFFIX_TO_REMOVE/}
        mkdir ${folder_name}
        cd $folder_name 
        pwd
        sicer -t ../$files -c ../*Input*.bam -s hg19 -w 200 -rt 1 -f 150 -egf 0.80 -fdr 0.01 -g 600
        mv *FDR0.01-island.bed ${folder_name}_island.bed
        cp *island.bed ..
        cd ..
done

mv *island.bed H1_islands/

#Change the name of the SICER files  
cd H1_islands 
mv *H1.0* H1.0_${CELL_LINE}_sicer_island.bed
mv *H1.2* H1.2_${CELL_LINE}_sicer_island.bed
mv *H1.4* H1.4_${CELL_LINE}_sicer_island.bed
mv *H1.5* H1.5_${CELL_LINE}_sicer_island.bed
mv *H1X* H1X_${CELL_LINE}_sicer_island.bed

### EDD ISLANDS
eval "$(conda shell.bash hook)" #Run this in order for conda to work 
conda activate Python27
mkdir EDD_islands
for i in *H1*.bam
do
        echo "EDD peak calling: $i"
        time edd /home/ajvlab/Softwares/hg19.chrom.sizes ../empty.bed $i Input*.bam ${i}_output
	cd ${i}_output
	mv *.bed ../EDD_islands
	cd ..
done

#Change the names of the EDD files 
cd EDD_islands 
mv *H1.0* H1.0_${CELL_LINE}_edd_island.bed
mv *H1.2* H1.2_${CELL_LINE}_edd_island.bed
mv *H1.4* H1.4_${CELL_LINE}_edd_island.bed
mv *H1.5* H1.5_${CELL_LINE}_edd_island.bed
mv *H1X* H1X_${CELL_LINE}_edd_island.bed

echo "Peak calling is OVER!"

#------------------------------------------------------------------------------------------------------------------------------#

# 2) Merge SICER and EDD files
mkdir merged_islands
VARIANTS=("H1.0" "H1.2" "H1.4" "H1.5" "H1X")
for variants in "${VARIANTS[@]}";
do
        sort -k1,1 -k2,2n H1_islands/*${variants}* | sed 's/[ \t]\+$//' > tmp_file_1.bed #Order the file and remove TABs
        sort -k1,1 -k2,2n EDD_islands/*${variants}* > tmp_file_2.bed #Order the files 
        cat tmp_file_1.bed tmp_file_2.bed | sort -k1,1 -k2,2n | bedtools merge > ${variants}_merged_island.bed #Merge files 
        mv *merged_island.bed merged_islands #Move to the corresponding file
        rm tmp_file*
done 

echo "Peak merging is OVER!"

#------------------------------------------------------------------------------------------------------------------------------#

# 3) Intersect H1 islands with UCSC LADs 
cd merged_islands/

WORD_TO_REMOVE=_merged_island

mkdir UCSC_intersection
for i in H1*
do 
        j=${i//$WORD_TO_REMOVE/}
        sed 's/[ \t]\+$//' $i > tmp_file #Remove the tabs at the end of the line, else it does not work
        bedtools intersect -a tmp_file -b $PATH_TO_UCSC_FILE > intersect_ucsc_${j}
        mv intersect_ucsc_${j} UCSC_intersection
	rm tmp_file
done
cd ..
echo "Intersection is OVER!"

#------------------------------------------------------------------------------------------------------------------------------#

# 4) CEAS: Genomic annotation
eval "$(conda shell.bash hook)" #Run this in order for conda to work 
conda activate Python27
mkdir ceas 
VARIANTS=("H1.0" "H1.2" "H1.4" "H1.5" "H1X")
for variants in "${VARIANTS[@]}";
do
        echo "------------CEAS: ${variants} --------------"
        ceas -b  merged_islands/*${variants}* -g /home/ajvlab/Softwares/hg19_ceas.refGene
        mkdir ${variants}_ceas
        mv *${variants}*.R *${variants}*.xls *${variants}*.pdf ${variants}_ceas/ 
        mv ${variants}_ceas/ ceas/
done  

echo "CEAS is OVER!"

#------------------------------------------------------------------------------------------------------------------------------#
# As output of this code we will obtain: "signal_H1_UCSC.bed" with the mapped signal, SICER islands (H1_islands), EDD islands
# (EDD_islands), EDD + SICER merged islands (merged islands), the intersection of the merged islands with UCS LADs
# (UCSC_intersection) and the genomic annotation of the H1 islands (ceas)
