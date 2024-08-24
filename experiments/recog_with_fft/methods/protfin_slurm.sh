# bash protfin_slurm.sh ../../../materials/protein.fa ../../../materials/mapmanreferencebins.results.txt

function sbatch_script() {
    name=protfin_WINSIZE_${window_size}_NPEAKS_${peaks}_OVERLAP_${overlap}_SKIPK_${skip_k_freqs}_SELMETH_${selection_method}_ALPHA_${significance}
    result_name=../results/${exp}/${name}
    echo "#!/bin/bash -l

# name
#SBATCH --job-name=${name}

# cpu
#SBATCH --ntasks=1

#SBATCH --mem-per-cpu=5GB

#SBATCH --output=../results/${exp}/_logs/%x_%j_slurm.out
#SBATCH --error=../results/${exp}/_logs/%x_%j_slurm.err

echo createdb >&2
python3 protfin.py create-db -s ${significance} -w ${window_size} -o ${overlap} -n ${peaks} -m ${selection_method} -k ${skip_k_freqs} -p ${result_name}.pickle $1

echo findmatches >&2
python3 protfin.py find-matches -d ${result_name}.pickle ../results/${exp}/_test_selection.fa > ${result_name}.matches

echo eval >&2
python3 evaluation.py eval ${result_name}.matches $2 > ${result_name}.summary.csv

echo matchfamily >&2
python3 protfin.py match-family <(awk '"'BEGIN{FS="\t"; print "Family_ID,Protein_ID"} $5=="T"{x=$1","toupper($3);gsub(/'"'\"'\"'"'/, "", x);print x}'"' "$2") -d ${result_name}.pickle > ${result_name}.match_fam.csv
"
}

exp=_v0.4-exp-selection_method
mkdir ../results/${exp}
mkdir ../results/${exp}/_logs
python3 evaluation.py select-samples $2 $1 -s 7 > ../results/${exp}/_test_selection.fa

for selection_method in deviation absolute none; do
    for (( window_size=50; window_size >= 30; window_size-=10 )); do
        declare -a overlaps=()
        declare -i overlap
        for i in 2 4; do
            overlap=${window_size}-${window_size}/${i}
            overlaps+=(${overlap})
        done
        overlap=${window_size}-3*${window_size}/4
        overlaps+=(${overlap})

        for overlap in "${overlaps[@]}"; do
            for peaks in 5 3 0; do
                for significance in 0.01 0.001 0.1 5; do
                    for skip_k_freqs in 0 1 2 3; do
                        sbatch --chdir . <(sbatch_script $*)
                    done
                done
            done
        done
    done
done