#! python3
from typing import List, Dict, Tuple, TextIO
import numpy as np
import pandas as pd
import pickle
import argparse


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


def create_constellation(seq_vec):
    constellation_map = []

    for idx, aa_vec in enumerate(seq_vec):
        for peak in aa_vec:
            constellation_map.append([idx, int(peak * 100)])

    return constellation_map


def find_match(fasta_file: str):
    aa_vec_map = get_aa_vectors()
    with open(fasta_file) as f:
        for (prot_id, description, seq) in iter_fasta(f):
            hashes = hashes_from_seq(seq, aa_vec_map)

            database = pickle.load(open('database.pickle', 'rb'))
            scores = score_prots(hashes, database)

            prot_index_lookup = pickle.load(open("prot_index.pickle", "rb"))
            for prot_index, score in scores:
                match_prot_description = prot_index_lookup[prot_index]
                print(f"{prot_index} - {match_prot_description}: Score of {score[1]} at {score[0]}")

            print("\nSeems to be: %s - %s" % (scores[0][0], prot_index_lookup[scores[0][0]]))
            print("\nInput:       %s - %s\n             %s" % (prot_id, description, seq))


def score_prots(hashes, database):
    matches_per_prot = {}
    for hash_, (sample_index, _) in hashes.items():
        if hash_ in database:
            matching_occurences = database[hash_]
            for source_index, prot_index in matching_occurences:
                if prot_index not in matches_per_prot:
                    matches_per_prot[prot_index] = []
                matches_per_prot[prot_index].append((sample_index, source_index))

    scores = {}
    for prot_index, matches in matches_per_prot.items():
        prot_scores_by_offset = {}
        for sample_index, source_index in matches:
            delta = source_index - sample_index
            if delta not in prot_scores_by_offset:
                prot_scores_by_offset[delta] = 0
            prot_scores_by_offset[delta] += 1

        max_ = (0, 0)
        for offset, score in prot_scores_by_offset.items():
            if score > max_[1]:
                max_ = (offset, score)
        scores[prot_index] = max_

    # Sort the scores for the user
    scores = list(sorted(scores.items(), key=lambda x: x[1][1], reverse=True))

    return scores


def create_db(prot_file: str):
    prot_index = {}
    database: Dict[str, List[Tuple[int, int]]] = {}
    aa_vec_map = get_aa_vectors()
    with open(prot_file) as f:
        for prot_id, description, seq in iter_fasta(f):
            if (proc := len(prot_index)) == 1000:
                break
            if proc % 100 == 0:
                print("%d%%" % int(proc/2))
            prot_index[prot_id] = description
            hashes = hashes_from_seq(seq, aa_vec_map, prot_id)
            for hash_, index_prot_id_pair in hashes.items():
                if hash_ not in database:
                    database[hash_] = []
                database[hash_].append(index_prot_id_pair)
    # %%
    with open("database.pickle", 'wb') as db:
        pickle.dump(database, db, pickle.HIGHEST_PROTOCOL)
    with open("prot_index.pickle", 'wb') as prots:
        pickle.dump(prot_index, prots, pickle.HIGHEST_PROTOCOL)


def create_hashes(constellation_map, prot_id=None):
    hashes = {}
    # assume pre-sorted
    # Iterate the constellation
    for idx, (index, freq) in enumerate(constellation_map):
        # Iterate the next 100 pairs to produce the combinatorial hashes
        for other_index, other_freq in constellation_map[idx : idx + 100]:
            diff = other_index - index
            # If the index difference between the pairs is too small or large
            # ignore this set of pairs
            if diff <= 1 or diff > 10:
                continue

            # Place the frequencies (in Hz) into a 1024 bins

            # Produce a 32 bit hash
            hash = int(freq) | (int(other_freq) << 10) | (int(diff) << 20)
            hashes[hash] = (index, prot_id)
    return hashes


def hashes_from_seq(seq: str, aa_vec_map, prot_id=None):
    audio = seq_to_vectors(seq, aa_vec_map)
    constellation = create_constellation(audio)
    return create_hashes(constellation, prot_id)


def seq_to_vectors(seq: str, aa_vec_map) -> list:
    return [aa_vec_map[aa] for aa in seq]


def get_aa_vectors(file="../../../materials/Amino_Acid_Kidera_Factors.csv") -> Dict[str, np.ndarray]:
    kidera = pd.read_csv(file).loc[:, "A":]

    special_aa = {
        "X": list(kidera),  # any aminoacid
        "B": ["D", "N"],
        "Z": ["E", "Q"],
        "J": ["I", "L"],
        "Ψ": ["I", "L", "M", "V"],
        "Ω": ["F", "W", "Y", "H"],
        "Φ": ["I", "L", "M", "V", "F", "W", "Y"],
        "ζ": ["D", "E", "H", "K", "N", "Q", "R", "S", "T"],
        "Π": ["A", "G", "P", "S"],
        "+": ["K", "R", "H"],
        "-": ["D", "E"]
    }
    for s in special_aa:
        kidera[s] = kidera[special_aa[s]].to_numpy().mean(axis=1)

    return kidera


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
