 #!/usr/bin/bash

#------------------------------------------------------------------------------------------------------------------------------#
# EXPLANATION OF THE CODE 
# This code will take the filtered and sorted bam files of the ChIPSeq signal of the lamins as input. Enirched Domain Detector (EDD) is used for the peak calling because Lamin associated domains (LADs) are longer in length and SICER do not capture them as a island but rather as noise. Due to this, we need to use EDD for bigger domain detection. The code will have 2 arguments: The path where the BED files are located (-p | --path) and the suffix it must be removed from the name of the file in order to obtain the "core" name (-s|--suffix).  
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
eval "$(conda shell.bash hook)" #EDD needs to work with conda because it has the previous version of python (2.7).
conda activate Python27

cd $PATH_TO_FILE
for files in ChIP*.bam
do
	echo "---------- Peak calling: $files ----------"
	folder_name=${files//$SUFFIX_TO_REMOVE/}
	mkdir ${folder_name}
	edd hg19.chrom.sizes empty.bed $files Input*bam ${folder_name} --bin-size 11 -g 5
done
#------------------------------------------------------------------------------------------------------------------------------#
# The OUTPUT of this code is a BED file containing the coordinates for the LADs identified for each of the ChIPSeq signal. 
#------------------------------------------------------------------------------------------------------------------------------#
