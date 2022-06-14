#!usr/bin/bash
#------------------------------------------------------------------------------------------------------------------------------#
# This script will take ChIP-Seq data for H1 variants and it maps against 3 different coordinates: UCSC LADs, Inter LADs (the
# regions that are not LADs) and random LADs (randomized UCSC LADs' coordinates). Moreover in this 3 cases, the blacklist
# regions were removed. The input files are BDG files of each of the H1 variants, and as output 3 bed files, each of them
# having the mapped signal of all 5 H1 variants against ILADs,LADs or random LADs. Those 3 file must be downloaded to later be
# introduced as input in H1_mapping_UCSC_LADs.r. 

# The parameters to be introduced are the path were the original UCSC LADs are stored ($PATH_TO_FILE,-p) and the other will be
# the working path were all the BDG files are stored and where new folders will be created ($WORKING_PATH, -w).
#------------------------------------------------------------------------------------------------------------------------------#
# CHECK ARGUMENTS
INSTRUCTIONS="INSTRUCTIONS= $0 -p path_ucsc_lad -w working_path"
if [ "$#" -ne "4" ]; then
	echo "Incorrect number of argument: $INSTRUCTIONS"
	exit 1
fi
case $1 in
	-p|--path_ucsc_lad)
	PATH_TO_FILE="$2"
	-w|--working_path)
	WORKING_PATH="$2"
esac

case $3 in
	-p|--path_ucsc_lad)
	PATH_TO_FILE="$4"
	-w|--working_path)
	WORKING_PATH="$2"
esac

# 1) MAP signal to UCSC LADs, inter LADs and randoms LADs 
SUFFIX_TO_REMOVE=.bdg
PATH_TO_UCSC_FILE=/data/AJV_TFM_Jaione/LADs_hg19.bed
cd $WORKING_PATH
bedtools shuffle -i $PATH_TO_UCSC_FILE -g /home/ajvlab/Softwares/hg19.chrom.sizes > random_LAD.bed
bedtools subtract -a random_LAD.bed -b /data/AJV_TFM_Jaione/Signal_mapping_No_BL/hg19-blacklist.v2.bed > random_LAD_BL.bed
mkdir map_UCSC
for k in H1*.bdg
do 
        echo "---------- BedGraph mapping: $k ----------"
        kn=${k//$SUFFIX_TO_REMOVE/} #New name for the output
        #Map in UCSC 
        echo "Map UCSC"
        bedtools map -a <(sort -k1,1 -k2,2n ../LADs_hg19_BL.bed)  -b <(sort -k1,1 -k2,2n $k) -c 4 -o mean > ${kn}_ucsc.bed
        #Map in Random LADs
        echo "Map RANDOM UCSC"
        bedtools map -a <(sort -k1,1 -k2,2n random_LAD_BL.bed) -b <(sort -k1,1 -k2,2n $k) -c 4 -o mean > ${kn}_zrandom.bed
        #Map in inter LADs
        echo "Map iLADs"
        bedtools map -a <(sort -k1,1 -k2,2n ../Inter_LADs_hg19_BL.bed) -b <(sort -k1,1 -k2,2n $k) -c 4 -o mean > ${kn}_iLAD.bed
        mv ${kn}_*.bed map_UCSC
done


cd map_UCSC
mkdir ilad
mv *iLAD.bed ilad
cd ilad

LAD=ILAD

files=(*)

# This piece of code will take the files generated in the previous step and merge all of the file into a single file with each 
# column storing the information of 1 file. 
awk '{print $1,$2,$3}' ${files[0]} > signal_H1_${LAD}.bed #Open an initial file with the coordinates of the LADs (3 columns)
file=signal_H1_${LAD}.bed
head signal_H1_${LAD}.bed 
for i in *H1*.bed
do
        echo "---------- Adding the signal of $i to file ----------"
        awk '{getline f1 <"signal_H1_ILAD.bed" ; print f1,$4}' < $i > signal_tmp.bed #New temporal file with 3 columns + 1 signal 
        rm signal_H1_${LAD}.bed #We remove the file to in the next step change the name to the signal_tmp file and go adding cols.
        mv signal_tmp.bed signal_H1_${LAD}.bed
        head signal_H1_${LAD}.bed #BORRAR
done 

#Add tabs to the file
cat signal_H1_${LAD}.bed | tr " " "\t" > signals_H1_${LAD}.bed
rm signal_H1_${LAD}.bed
cd ..
cd ..
echo "\n Signal mapping is in iLAD OVER! \n"

cd ..

mkdir lad
mv *ucsc.bed lad
cd lad 

LAD=UCSC
files=(*)

# This piece of code will take the files generated in the previous step and merge all of the file into a single file with each 
# column storing the information of 1 file. 
awk '{print $1,$2,$3}' ${files[0]} > signal_H1_${LAD}.bed #Open an initial file with the coordinates of the LADs (3 columns)
file=signal_H1_${LAD}.bed
head signal_H1_${LAD}.bed 
for i in *H1*.bed
do
        echo "---------- Adding the signal of $i to file ----------"
        awk '{getline f1 <"signal_H1_UCSC.bed" ; print f1,$4}' < $i > signal_tmp.bed #New temporal file with 3 columns + 1 signal 
        rm signal_H1_${LAD}.bed #We remove the file to in the next step change the name to the signal_tmp file and go adding cols.
        mv signal_tmp.bed signal_H1_${LAD}.bed
        head signal_H1_${LAD}.bed #BORRAR
done 

#Add tabs to the file
cat signal_H1_${LAD}.bed | tr " " "\t" > signals_H1_${LAD}.bed
rm signal_H1_${LAD}.bed
cd ..
cd ..
echo "\n Signal mapping in LADis OVER! \n"

cd ..

mkdir random
mv *zrandom.bed random 

LAD=random
files=(*)

# This piece of code will take the files generated in the previous step and merge all of the file into a single file with each 
# column storing the information of 1 file. 
awk '{print $1,$2,$3}' ${files[0]} > signal_H1_${LAD}.bed #Open an initial file with the coordinates of the LADs (3 columns)
file=signal_H1_${LAD}.bed
head signal_H1_${LAD}.bed 
for i in *H1*.bed
do
        echo "---------- Adding the signal of $i to file ----------"
        awk '{getline f1 <"signal_H1_random.bed" ; print f1,$4}' < $i > signal_tmp.bed #New temporal file with 3 columns + 1 signal 
        rm signal_H1_${LAD}.bed #We remove the file to in the next step change the name to the signal_tmp file and go adding cols.
        mv signal_tmp.bed signal_H1_${LAD}.bed
        head signal_H1_${LAD}.bed #BORRAR
done 

#Add tabs to the file
cat signal_H1_${LAD}.bed | tr " " "\t" > signals_H1_${LAD}.bed
rm signal_H1_${LAD}.bed
cd ..
cd ..
echo "\n Signal mapping in random is OVER! \n"

cd ..
