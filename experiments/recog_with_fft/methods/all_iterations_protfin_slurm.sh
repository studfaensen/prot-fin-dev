# bash protfin_slurm.sh ../../../materials/protein.fa ../../../materials/mapmanreferencebins.results.txt

function sbatch_script() {
    name=protfin_WINSIZE_${window_size}_NPEAKS_${peaks}_OVERLAP_${overlap}_SKIPK_${skip_k_freqs}_SELMETH_${selection_method}ALPHA${significance}
    result_name=../results/${exp}/${name}
    echo "#!/bin/bash -l

# name
#SBATCH --job-name=${name}

# cpu
#SBATCH --ntasks=1

#SBATCH --mem-per-cpu=5GB

#SBATCH --output=../results/${exp}/_logs/%x%j_slurm.out
#SBATCH --error=../results/${exp}/_logs/%x%j_slurm.err

echo createdb >&2
python3 protfin.py create-db -s ${significance} -w ${window_size} -o ${overlap} -n ${peaks} -m ${selection_method} -k ${skip_k_freqs} -p ${result_name}.pickle $1

echo findmatches >&2
python3 protfin.py find-matches -d ${result_name}.pickle ../results/${exp}/_test_selection.fa > ${result_name}.matches

echo eval >&2
python3 evaluation.py eval ${result_name}.matches $2 > ${result_name}.summary.csv
"
}

exp=_v0.4-exp-selection_method
mkdir ../results/${exp}
mkdir ../results/${exp}/_logs
python3 evaluation.py select-samples $2 $1 -s 7 > ../results/${exp}/_test_selection.fa


window_size=50
overlap=25
peaks=3
significance=0.001
skip_k_freqs=2
selection_method=none
cat <(sbatch_script $*)
#sbatch --chdir . <(sbatch_script $*)
