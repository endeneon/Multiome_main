#!/bin/bash

echo "file_name\tread_count_(2x_pairs)" > count_results.txt

for eachfile in *.bam
do
	echo $eachfile
	samtools index -@ 16 $eachfile
	printf "$eachfile\t" >> count_results.txt

	samtools idxstats $eachfile \
		| cut -f 3 \
		| paste -s -d+ - \
		| bc \
		>> count_results.txt
done


