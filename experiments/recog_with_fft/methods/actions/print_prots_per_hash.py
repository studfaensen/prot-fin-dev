import pickle


def print_prots_per_hash(database: str):
    with open(database, "rb") as f:
        database = pickle.load(f).db

    print(*(len(prots) for prots in database.values()), sep=",", end="")
