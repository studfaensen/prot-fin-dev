import pickle


def print_prots_per_hash(database: str):
    with open(database, "rb") as f:
        database, _ = pickle.load(f)

    print(*(len(prots) for prots in database.values()), sep=",", end="")
