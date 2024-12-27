#! /bin/bash
# Siwei 25 Apr 2025

input_bed=$1



findMotifsGenome.pl \
	$input_bed \
	hg38 \
	${input_bed/%.bed/} \
	-size 50 \
	-nomotif \
	-p 20 \
	-bits

