import pandas as pd
from sys import argv as args

data = pd.read_csv(args[1], sep=",")

freqs_per_win = data["Mean_Frequencies_Per_Window_Without_First_npeaks_Freqs"]
window_size = data["Window_Size"]
n_peaks = data["N_Peaks"]
overlap = data["Overlap"]
target_zone = 8
kidera_factors = 10
median_sequence_length = 300

window_count = median_sequence_length / (window_size - overlap)
freqs_absolute = window_count * freqs_per_win
hashes_per_freq = target_zone * freqs_per_win

hash_count = hashes_per_freq * freqs_absolute

db_size_per_prot = hash_count * kidera_factors * 4  # Byte
proteins = 41669
data["est_DB_Size_MB"] = round(db_size_per_prot * proteins / (1024 ** 2), 2)
print(
    data[data["est_DB_Size_MB"] < 20]
    # data
    .sort_values("Mean_Frequencies_Per_Window_Without_First_npeaks_Freqs", ascending=False)
    .to_csv(index=False)
)
