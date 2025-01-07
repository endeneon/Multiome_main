while IFS=' ' read -r arg1 arg2 arg3; do
    ./eval_covar.sh "$arg1" "$arg2" "$arg3"
done < eval_covar_inputs.txt

