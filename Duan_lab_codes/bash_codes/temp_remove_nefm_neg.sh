#!/bin/bash
# Siwei 06 Dec 2022

# temporary remove nefm_neg from the barcode files
# since this cell category has been subsetted

for eachdir in group*
do
	echo $eachdir

	mkdir $eachdir/temp_not_to_process
	mv $eachdir/*NEFM_neg_glut_barcodes.txt $eachdir/temp_not_to_process

done
