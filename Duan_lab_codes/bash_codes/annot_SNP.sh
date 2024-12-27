#!/bin/bash

# Siwei 11 Feb 2022
# !! Ensure that a functional Perl > 5.26 is active in your env !!

annovar="/home/zhangs3/Data/Tools/annovar/annotate_variation.pl"

mkdir -p "annot_output"


for eachfile in *.avinput
do
	echo $eachfile
	$annovar \
		--geneanno \
		-out annot_output/${eachfile/%.avinput/.geneanno} \
		-buildver hg38 \
		-dbtype knownGene \
		$eachfile \
		/home/zhangs3/Data/Tools/annovar/humandb

#        $annovar \
#		--regionanno \
#	        -out annot_output/${eachfile/%.avinput/.regionanno} \
#	        -buildver hg38 \
#                -dbtype gwasCatalog \
#		$eachfile \
#	        /home/zhangs3/Data/Tools/annovar/humandb
done


