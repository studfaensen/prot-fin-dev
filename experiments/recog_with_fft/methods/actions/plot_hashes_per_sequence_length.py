import pickle
import matplotlib.pyplot as plt
from .algorithm.hash_gen import FREQUENCY_BITS, DIFFERENCE_BITS


def plot_hashes_per_sequence_length(database: str, out_file: str):
    with open(database, "rb") as f:
        db = pickle.load(f)

    hashes_per_seqlen = {}
    for seqlen, hash_count in db.lookup.values():
        win_count = seqlen // (db.config.window_size - db.config.overlap) + 1
        hashes_per_seqlen[win_count] = hashes_per_seqlen.get(win_count, []) + [hash_count]

    plt.figure(figsize=(20, 10))
    plt.boxplot(hashes_per_seqlen.values(), positions=list(hashes_per_seqlen.keys()))

    xmin, xmax = plt.xlim()
    plt.xticks(range(int(xmin), int(xmax), 20), range(int(xmin), int(xmax), 20))
    plt.xlim(xmin - .5, xmax + .5)
    plt.xlabel("Sequence length in number of STFT windows")
    plt.ylabel("Hash counts")
    plt.title("Hash counts per sequence length")

    plt.savefig(out_file, bbox_inches='tight')
