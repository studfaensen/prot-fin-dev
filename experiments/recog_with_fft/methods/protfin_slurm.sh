#!/bin/bash -l

# name
#SBATCH --job-name=frequencies

# cpu
#SBATCH --ntasks=80

#SBATCH --mem-per-cpu=1GB

#SBATCH --output=../results/_v0.4-exp-frequency_selection/_logs/%x_%j_slurm.out
#SBATCH --error=../results/_v0.4-exp-frequency_selection/_logs/%x_%j_slurm.err

exp=_frequencies
mkdir ../results/${exp}
mkdir ../results/${exp}/_logs

csv=../results/${exp}/mean_frequencies.csv
echo Window_Size,N_Peaks,Overlap,Significance,Mean_Frequencies_Per_Window,Mean_Frequencies_Per_Window_Without_First_npeaks_Freqs,Mean_Quantile_Without_First_npeaks_Freqs > ${csv}

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
        for peaks in 3 2 1; do
            for significance in 5 .1 .01 .001; do
                name=frequencies_WINSIZE_${window_size}_NPEAKS_${peaks}_OVERLAP_${overlap}_SIGNIFICANCE_${significance}
                result_name=../results/${exp}/${name}

                echo -en ${window_size},${peaks},${overlap},${significance}, >> ${csv}

                python3 evaluation.py plot-frequencies -n ${peaks} -w ${window_size} -o ${overlap} -s ${significance} ../../../materials/protein.fa ${result_name}.png -c 80 >> ${csv}
            done
        done
    done
done
