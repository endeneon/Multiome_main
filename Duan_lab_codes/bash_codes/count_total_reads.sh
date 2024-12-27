#!/bin/bash

# Siwei 05 Jul 2022

for eachfile in *.bam
do
	samtools flagstat -@ 10 \
		$eachfile \
		| grep '0 mapped' \
		| cut -f 1 \
		>> total_counts.txt
done

cat total_counts.txt \
	| cut -d ' ' -f 1 \
	| paste -s -d '+' \
	| bc
