import pickle
import matplotlib.pyplot as plt
from .algorithm.hash_gen import FREQUENCY_BITS, DIFFERENCE_BITS


def plot_prots_per_windist(database: str, out_file: str):
    with open(database, "rb") as f:
        database = pickle.load(f).db

    prots_per_windist = {}
    for hash_, prots in database.items():
        windist = (hash_ >> 2*FREQUENCY_BITS) % 2**DIFFERENCE_BITS
        prots_per_windist[windist] = prots_per_windist.get(windist, []) + [len(prots)]

    plt.figure(figsize=(20, 10))
    plt.boxplot(prots_per_windist.values(), positions=list(prots_per_windist.keys()))

    xmin, xmax = plt.xlim()
    plt.xticks(range(int(xmin), int(xmax), 20), range(int(xmin), int(xmax), 20))
    plt.xlim(xmin - .5, xmax + .5)
    plt.xlabel("Window distance")
    plt.ylabel("Proteins")
    plt.title("Protein counts per hashes' window distances")

    plt.savefig(out_file, bbox_inches='tight')
