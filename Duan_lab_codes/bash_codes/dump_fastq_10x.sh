#!/bin/bash

# Siwei 14 Nov 2022

bam2fastq="/home/zhangs3/Data/Tools/bam2fastq_10x/bamtofastq_linux"

mkdir -p dumped_fastqs

for eachfile in *.bam
do
	echo $eachfile
	$bam2fastq \
		--nthreads=16 \
		--relaxed \
		$eachfile \
		"dumped_fastqs"/${eachfile/%hg38only.bam/}
done
