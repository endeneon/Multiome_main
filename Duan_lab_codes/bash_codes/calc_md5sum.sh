#!/bin/bash

# Siwei 09 Mar 2022

rm md5_cal_ed.txt

for eachfile in *.gz
do
	echo $eachfile
	md5sum $eachfile >> md5_calc_ed.txt
done
