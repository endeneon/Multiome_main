#! /bin/bash

# find and pack these files:
# atac_fragments.*
# cloupe.cloupe
# filtered_feature_bc_matrix.h5
# gex_molecule_info.h5
# Use py_312 conda env

depo_dir="../GEO_deposit"
mkdir -p $depo_dir


shopt -s globstar;

for eachdir in Duan_024_*
do
	echo $eachdir

	tar \
		--use-compress-program="pigz -9 -p 32 --recursive | pv" \
		-cvf \
		$depo_dir"/"$eachdir".tar.gz" \
		"$eachdir/"**/atac_fragments.* \
		"$eachdir/"**/cloupe.cloupe \
		"$eachdir/"**/filtered_feature_bc_matrix.h5 \
		"$eachdir/"**/gex_molecule_info.h5
done
