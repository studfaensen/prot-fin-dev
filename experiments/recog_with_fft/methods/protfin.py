#! python3
from typing import List, Dict, Tuple, TextIO
import pickle
import argparse
from tools import *
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)  # fixes weird python error, look: https://newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html

# experimental
KIDERA_FACTOR = ["Helix/bend preference", "Side-chain size", "Extended structure preference", "Hydrophobicity", "Double-bend preference", "Partial specific volume", "Flat extended preference", "Occurrence in alpha region", "pK-C", "Surrounding hydrophobicity"]\
    .index("Hydrophobicity")
WINDOW_SIZE = 30
WINDOW_TYPE = ["boxcar", "triang", "blackman", "hamming", "hann", "bartlett", "flattop", "parzen", "bohman", "blackmanharris", "nuttall", "barthann", "cosine", "exponential", "tukey", "taylor", "lanczos"]\
              [0]
OVERLAP = 15
N_PEAKS = 0  # 0 means all
MAX_ANKER_LENGTH = 20 * WINDOW_SIZE # in shazam algorithm it is the twentyfold of windowlength


def main():
    parser = argparse.ArgumentParser(
        prog='prot-fin',
        description='protein recognition by transforming their sequences into music',
    )
    sub_commands = parser.add_subparsers(required=True)

    create_db_parser = sub_commands.add_parser("create-db", help="Create Database")
    create_db_parser.add_argument("fasta-file")
    create_db_parser.set_defaults(func=lambda a: create_db(getattr(a, "fasta-file")))

    find_match_parser = sub_commands.add_parser("find-match", help="Find Matches for Proteins")
    find_match_parser.add_argument("fasta-file")
    find_match_parser.set_defaults(func=lambda a: find_match(getattr(a, "fasta-file")))

    args = parser.parse_args()
    args.func(args)


def find_match(fasta_file: str):
    with open(fasta_file) as f:
        for (prot_id, description, seq) in iter_fasta(f):
            hashes = hashes_from_seq(seq)

            database = pickle.load(open('database.pickle', 'rb'))
            prot_index_lookup = pickle.load(open("prot_index_lookup.pickle", "rb"))
            scores = score_prots(hashes, database, prot_index_lookup)

            for prot_index, score in scores:
                match_prot_description, _ = prot_index_lookup[prot_index]
                print(f"{prot_index} - {match_prot_description}: Jaccard Index of {score[2]} : Score of {score[1]}")

            print("\nSeems to be: %s - %s" % (scores[0][0], prot_index_lookup[scores[0][0]][0]))
            print("\nInput:       %s - %s\n             %s" % (prot_id, description, seq))
            print("\nFound hashes: %d" % len(hashes))


def score_prots(hashes, database, prot_index_lookup):
    matches_per_prot = {}
    for hash_, (sample_index, _) in hashes.items():
        if hash_ in database:
            matching_occurences = database[hash_]
            for source_index, match_prot_index in matching_occurences:
                if match_prot_index not in matches_per_prot:
                    matches_per_prot[match_prot_index] = []
                matches_per_prot[match_prot_index].append((sample_index, source_index))

    scores = {}

    for match_prot_index, matches in matches_per_prot.items():
        # jaccard similarity index
        search_size = len(hashes)
        _, match_size = prot_index_lookup[match_prot_index]
        intersection_size = len(matches)
        union_size = search_size + match_size - intersection_size
        jacc_sim_index = intersection_size / union_size

        # original scoring
        prot_scores_by_offset = {}
        for sample_index, source_index in matches:
            delta = source_index - sample_index
            if delta not in prot_scores_by_offset:
                prot_scores_by_offset[delta] = 0
            prot_scores_by_offset[delta] += 1

        max_ = (0, 0, jacc_sim_index)
        for offset, score in prot_scores_by_offset.items():
            if score > max_[1]:
                max_ = (offset, score, jacc_sim_index)
        scores[match_prot_index] = max_
        # scores[match_prot_index] = (0, jacc_sim_index)

    # Sort the scores for the user
    scores = list(sorted(scores.items(), key=lambda x: x[1][2]*x[1][1], reverse=True))

    return scores


def create_db(prot_file: str):
    import numpy as np
    prot_index_lookup = {}
    database: Dict[str, List[Tuple[int, int]]] = {}

    with open(prot_file) as f:
        for prot_id, description, seq in iter_fasta(f):
            if (proc := len(prot_index_lookup)) == 1000:
                break
            if proc % 100 == 0:
                print("%d%%" % int(proc / 10))
            hashes = hashes_from_seq(seq, prot_id)
            prot_index_lookup[prot_id] = (description, len(hashes))
            for hash_, index_prot_id_pair in hashes.items():
                if hash_ not in database:
                    database[hash_] = []
                database[hash_].append(index_prot_id_pair)
    # %%
    with open("database.pickle", 'wb') as db:
        pickle.dump(database, db, pickle.HIGHEST_PROTOCOL)
    with open("prot_index_lookup.pickle", 'wb') as prots:
        pickle.dump(prot_index_lookup, prots, pickle.HIGHEST_PROTOCOL)


def hashes_from_seq(seq: str, prot_id=None):
    aa_vec = get_aa_vector(seq, KIDERA_FACTOR)
    constellation = create_constellation(aa_vec, WINDOW_SIZE, N_PEAKS, overlap=OVERLAP)
    return create_hashes(constellation, prot_id)


def iter_fasta(fasta_file: TextIO) -> Tuple[str, str]:
    while True:
        prot_desc = fasta_file.readline()

        if "\n" not in prot_desc:  # -> EOF
            break
        if prot_desc[0] == ">":
            prot_id, _, description = prot_desc.split(" ", 2)
            seq = fasta_file.readline()
            if seq[-1] == "\n":
                seq = seq[:-1]

            yield prot_id[1:], description[:-1], seq


if __name__ == '__main__':
    main()
