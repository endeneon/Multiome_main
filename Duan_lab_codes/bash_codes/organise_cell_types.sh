#!/bin/bash

# Siwei 28 Jul 2022

mkdir -p hr_0_GABA
mkdir -p hr_0_NEFM_pos
mkdir -p hr_0_NEFM_neg

rm *.idx

mv *_GABA_*.vcf hr_0_GABA/
mv *_NEFM_pos_*.vcf hr_0_NEFM_pos/
mv *_NEFM_neg_*.vcf hr_0_NEFM_neg/
