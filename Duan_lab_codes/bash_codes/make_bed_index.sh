#!/bin/bash

# Siwei 12 Nov 2022

# make idx files for bed lists

gatk4="/home/zhangs3/Data/Tools/gatk-4.2.6.1/gatk"

for eachfile in *.bed
do
#	$gatk4 \
#		IndexFeatureFile \
#		-I $eachfile
	cat $eachfile \
		| bgzip \
		> $eachfile.bgzip

	tabix \
		-p bed \
		-f $eachfile.bgzip
done
