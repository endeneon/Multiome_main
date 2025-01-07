while IFS=' ' read -r col1 col2 col3; do
    ./caQTL_mapping.sh "$col1" "$col2" "$col3"
done < caQTL_mapping_inputs.txt

