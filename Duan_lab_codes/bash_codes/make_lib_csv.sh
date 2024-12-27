#!/bin/bash
# Siwei 11 Mar 2022

readarray -t sample_names < sample_names.txt

ATAC_path=$PWD'/ATAC'
GEX_path=$PWD'/GEX'

ATAC_suffix='_ATAC_FL,Chromatin Accessibility'
GEX_suffix='_GEX_FL,Gene Expression'


for each_sample in ${sample_names[*]}
do
	echo $each_sample
	echo 'fastqs,sample,library_type' >> $each_sample"_input_lib.csv"
	echo "$ATAC_path,$each_sample$ATAC_suffix" >> $each_sample"_input_lib.csv"
	echo "$GEX_path,$each_sample$GEX_suffix" >> $each_sample"_input_lib.csv"

done

