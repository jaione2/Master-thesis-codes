#!/bin/bash

#------------------------------------------------------------------------------------------------------------------------------#
# EXPLANATION OF THE CODE
# The input of this code will be the BAM file obtained from the alignment of ChIPSeq experiment data. There will be 2 arguments
# given in the command line: The path to the working directory where the bam files are (-p,--path) and the suffix to remove 
# (-s --suffix). The suffix to remove will be a single (e.g. .bam) or a group of extensions (.fq.gz.bam) apart from the "core"
# name of the file. This suffix will be used to obtain the core name in a variable, and generate the file name in a more
# compacted way.
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

#------------------------------------------------------------------------------------------------------------------------------#

#Enter into the working directory and iterate over the file with .bam extension
cd $PATH_TO_FILE
for i in *.bam 
do
 	echo "---------- Processing file: $i ----------"
        i=${i//$SUFFIX_TO_REMOVE/} #Obtain the core name of the file
# 1) BAM to BedGraph conversion
	echo "---------- BAM to BedGraph conversion ----------"
        mkdir bedgraph_files
        N=$(samtools view -c $i.fq.gz.bam) #N the number of reads in the file after filtering them 
        echo "$N"
        RPM=$(echo "1 / $N * 1000000" | bc -l) 
        echo "$RPM" 
        bedtools genomecov -ibam $i.fq.gz.bam -scale $RPM -bga > $i.bdg
        mv $i.bdg bedgraph_files
done

#------------------------------------------------------------------------------------------------------------------------------#

# 2) Input subtraction 
eval "$(conda shell.bash hook)" #MACs need to work with conda because it has the previous version of python (2.7).
conda activate Python27

SUFFIX_TO_REMOVE=.bdg
 
mkdir subtracted_files
cd bedgraph_files
for j in *ChIP*.bdg
do	
	echo "---------- Input subtraction: $j ----------"
        j=${j//$SUFFIX_TO_REMOVE/}
        macs2 bdgcmp -m subtract -t $j.bdg -c Input*.bdg -o $j.subtract.bdg
        mv $j.subtract.bdg $PATH_TO_FILE/subtracted_files
done 
cd $PATH_TO_FILE

#------------------------------------------------------------------------------------------------------------------------------#

# 3) Map ChIP-Seq signal against genome windows of 100kb. The genome is divided in windows of 100kb in /home/ajvlab/Nuria/Files/Coords100k.bed

SUFFIX_TO_REMOVE=.bdg
cd $PATH_TO_FILE
mkdir $PATH_TO_FILE/H1_100kb_bed_files
cd subtracted_files
for k in *.bdg
do 
	echo "---------- BedGraph mapping: $j ----------"
        k=${k//$SUFFIX_TO_REMOVE/}
        bedtools map -a /home/ajvlab/Nuria/Files/Coords100k.bed -b <(sort -k1,1 -k2,2n ${k}.bdg) -c 4 -o mean > ${k}_100kb.bed
        mv ${k}_100kb.bed $PATH_TO_FILE/H1_100kb_bed_files
done

cd $PATH_TO_FILE

#------------------------------------------------------------------------------------------------------------------------------#

# 4) Join in a single file the signal for each of the files

# In the file were we have stored BED files, first we are going to take the first file and copy the first 3 columns into a new # file (abundances_H1_Lamins.bed). Then, the loop will go for all the BED files. This loop will take as base the previously
# generated file, and it will print in the output what the file already had + the forth column of the $i file. This will be
# printed in a new temporal file (abundances_tmp.bed). #We will remove the "abundances_H1_Lamins.bed" file and rename the tmp
# file with the same new. This will allow us to continue iterating and adding the abundance column. 

cd $PATH_TO_FILE/H1_100kb_bed_files
files=(*)

awk '{print $1,$2,$3}' ${files[0]} > abundances_H1_Lamins.bed

for i in *100kb.bed
do
	echo "---------- Adding $i file ----------"
	awk '{getline f1 <"abundances_H1_Lamins.bed" ; print f1,$4}' < $i > abundances_tmp.bed
	rm abundances_H1_Lamins.bed
	mv abundances_tmp.bed abundances_H1_Lamins.bed
done 

#At last, the generated file is separated by space so translate the spaces to tabs. 
cat abundances_H1_Lamins.bed | tr " " "\t" > abundances_H1_Lamins_tab.bed

# THE OUTPUT FILE (abundances_H1_Lamins_tab.bed) will contain 3 columns (chr, init, end) + the number of files from which the
# signal has been obtained. This file should be introduced to Rstudio H1_Lamins_cor.r code to visualize the results.  

#------------------------------------------------------------------------------------------------------------------------------#
