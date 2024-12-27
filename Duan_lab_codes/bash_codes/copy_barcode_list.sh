#!/bin/bash

# Siwei 12 Nov 2022

time_prefix="6hr"
path_prefix="../../barcodes_by_time_4_WASP_rmdup/"

for eachfile in *.txt
do
#	echo $eachfile
	echo "$path_prefix${eachfile/0hr/$time_prefix}"
	cp "$path_prefix${eachfile/0hr/$time_prefix}" .
done

mkdir -p processed

mv 0hr* processed/
