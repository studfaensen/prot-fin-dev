from tools import *
from .algorithm import create_constellation, get_aa_vector
import pickle
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap as Cmap, Normalize
from matplotlib.cm import ScalarMappable


def plot_frequencies(prot_file: str, out_file: str):
    selected_freqs = {}
    fasta = Fasta(prot_file)
    protein_count = len(fasta)

    for _, _, seq in fasta:
        constellation_map = create_constellation(get_aa_vector(seq))

        for _, freq in constellation_map:
            freq_count, prot_count = selected_freqs.get(freq, (0, 0))
            selected_freqs[freq] = (freq_count + 1, prot_count)

        for freq in set(f for _, f in constellation_map):
            freq_count, prot_count = selected_freqs[freq]
            selected_freqs[freq] = (freq_count, prot_count + 1)

    cmap = Cmap.from_list(
        'custom_cmap',
        [(0, 'orange'), (1/3, 'red'), (1, 'black')]
    )
    norm = Normalize(vmin=100, vmax=protein_count)

    freqs, counts = zip(*sorted(selected_freqs.items()))
    freq_counts, prot_counts = zip(*counts)

    fig = plt.figure(figsize=(10, 5))
    plt.scatter(freqs, freq_counts, marker='o', c=list(map(lambda c: cmap(norm(c)), prot_counts)))
    for outlier_freq, (outlier_freq_count, _) in filter(lambda i: i[1][0] > max(freq_counts) / 2, selected_freqs.items()):
        plt.text(
            outlier_freq,
            outlier_freq_count - max(freq_counts) / 100 * 1.5,
            round(outlier_freq, 2),
            ha="center",
            va="top"
        )
    plt.colorbar(ScalarMappable(cmap=cmap, norm=norm)).set_label("Sequences", rotation=-90, va="bottom")
    plt.xlabel("Frequencies")
    plt.ylabel("Absolute count in all sequences")
    plt.title("Occurences of selected STFT frequences")
    plt.savefig(out_file)
