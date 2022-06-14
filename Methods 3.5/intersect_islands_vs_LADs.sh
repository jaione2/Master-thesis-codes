# 4) Intersect H1 islands with UCSC LADs 
#--------------------------------------------------------------------------------------------------------------------------------#
# This piece of code is run in the file were the BED files of the islands are located. In this case they are intersected against
# UCSC LADs, but the feature against which we want to intersect can be changed. The instersection output file for each of the H1	
# variants will be stored in the folder created (UCSC_intersection)
#--------------------------------------------------------------------------------------------------------------------------------#
PATH_TO_UCSC_FILE=/data/AJV_TFM_Jaione/LADs_hg19.bed
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
