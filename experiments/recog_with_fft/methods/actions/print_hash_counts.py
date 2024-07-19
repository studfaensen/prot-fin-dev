import pickle


def print_hash_counts(database: str):
    with open(database, "rb") as f:
        lookup = pickle.load(f).lookup

    _, hash_counts = zip(*lookup.values())
    print(*hash_counts, sep=",", end="")
