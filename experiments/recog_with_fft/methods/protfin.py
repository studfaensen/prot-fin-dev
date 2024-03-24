#! python3
"""
A tool to find a protein and its relatives based on chemical-physical
features of their amino acid sequences
"""

import pickle
import argparse
from tools import *
from os import environ as env
from typing import List, Dict, Tuple
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)  # fixes weird python error, look: https://newbebweb.blogspot.com/2012/02/python-head-ioerror-errno-32-broken.html

KIDERA_FACTOR = ["Helix/bend preference", "Side-chain size", "Extended structure preference", "Hydrophobicity", "Double-bend preference", "Partial specific volume", "Flat extended preference", "Occurrence in alpha region", "pK-C", "Surrounding hydrophobicity"]\
    .index("Hydrophobicity")

# parameters for the STFT
WINDOW_SIZE = int(env.get("WINDOW_SIZE", 30))
WINDOW_TYPE = env.get("WINDOW_TYPE", ["boxcar", "triang", "blackman", "hamming", "hann", "bartlett", "flattop", "parzen", "bohman", "blackmanharris", "nuttall", "barthann", "cosine", "exponential", "tukey", "taylor", "lanczos"]\
                                     [0])
OVERLAP = int(env.get("OVERLAP", 15))
N_PEAKS = int(env.get("N_PEAKS", 0))  # 0 means all

DB_DEFAULT = "database.pickle"
LOOKUP_DEFAULT = "protein_lookup.pickle"


def main():
    """
    This is where the command line interface is defined to interpret
    the passed arguments to execute the respective functions.
    """

    parser = argparse.ArgumentParser(
        prog='prot-fin',
        description=__doc__,
    )
    sub_commands = parser.add_subparsers(required=True)

    # protfin.py create-db <fasta-file>
    create_db_parser = sub_commands.add_parser("create-db", help="Create Database")
    create_db_parser.add_argument("fasta-file")
    create_db_parser.add_argument("-p", "--path", default=DB_DEFAULT)
    create_db_parser.set_defaults(func=lambda args: create_db(getattr(args, "fasta-file"), db_out=args.path))

    # protfin.py find-match <fasta-file>
    find_match_parser = sub_commands.add_parser("find-match", help="Find Matches for Proteins")
    find_match_parser.add_argument("fasta-file")
    find_match_parser.add_argument("-d", "--database", default=DB_DEFAULT)
    find_match_parser.set_defaults(func=lambda args: find_match(getattr(args, "fasta-file"), db_in=args.database))

    args = parser.parse_args()
    args.func(args)


def create_db(
        prot_file: str,
        db_out=DB_DEFAULT,
        lookup_out=LOOKUP_DEFAULT
        ):
    """
    Creates a database of combinatorial hashes for all protein sequences in the
    given FASTA file and also a lookup database of protein identifiers pointing
    to their description and number of hashes

    ...

    Parameters
    ----------
    prot_file : str
        The path to the FASTA formatted file containing all train sequences of
        known protein sequences
    db_out : str
        Name of the file to write the database to
    lookup_out : str
        Name of the file to write the protein lookup to
    """

    protein_lookup: ProteinLookup = {}
    database: Database = {}

    # train the database with hashes pointing to their matching proteins
    for prot_id, description, seq in Fasta(prot_file):

        # calculate the combinatorial hashes for the sequence
        hashes = hashes_from_seq(seq, prot_id)

        # save the protein related data
        protein_lookup[prot_id] = (description, len(hashes))

        # sort the hashes into the database by using them as key pointing to
        # their matching proteins
        for hash_, index_prot_id_pair in hashes.items():
            if hash_ not in database:
                database[hash_] = []
            database[hash_].append(index_prot_id_pair)

    # write the databases into files
    with open(db_out, 'wb') as db:
        pickle.dump(database, db, pickle.HIGHEST_PROTOCOL)
    with open(lookup_out, 'wb') as prots:
        pickle.dump(protein_lookup, prots, pickle.HIGHEST_PROTOCOL)


def find_match(
        fasta_file: str,
        db_in=DB_DEFAULT,
        lookup_in=LOOKUP_DEFAULT
        ):
    """
    Find matches for the proteins defined in the FASTA file
    and print them to stdout

    ...

    Parameters
    ----------
    fasta_file : str
        The path to the FASTA formatted file containing all protein sequences
        of interest
    db_in : str
        Name of the file storing the trained database
    lookup_in : str
        Name of the file storing the protein lookup
    """

    for (prot_id, description, seq) in Fasta(fasta_file):

        # create the combinatorial hashes for the sequence
        hashes: Hashes = hashes_from_seq(seq, prot_id)

        # load databases
        with open(db_in, 'rb') as f:
            database: Database = pickle.load(f)
        with open(lookup_in, 'rb') as f:
            protein_lookup: ProteinLookup = pickle.load(f)

        # calculate the scores for proteins in the database
        scores: Scores = \
            score_prots(hashes, database, protein_lookup)

        if scores:
            top_scored = filter(lambda score: scores[0][1][1] == score[1][1] and scores[0][1][2] == score[1][2], scores)
            scores = list(top_scored)

        # print the matches with description and score
        for prot_index, (_, score, jacc_sim_index) in scores:
            match_prot_description, _ = protein_lookup[prot_index]
            print(f"{prot_index} - {match_prot_description}: Jaccard Index of {jacc_sim_index} : Score of {score}")

        # print first match as most likely match
        if len(scores):
            first_match_id: str = scores[0][0]
            first_match_description: str = protein_lookup[first_match_id][0]
            print("\nSeems to be: %s - %s" % (first_match_id, first_match_description))
        else:
            print("\nNo matches found")

        print("\nInput:       %s - %s" % (prot_id, description))
        print(seq)
        _, score, jsi = sorted(scores, key=lambda match: match[0] == prot_id)[-1][1]
        print("Input-JSI: %s  -  Input-Score: %s" % (jsi, score))
        print("\nFound hashes: %d" % len(hashes))


def hashes_from_seq(seq: str, prot_id: ProteinID) -> Hashes:
    """
    Generate the combinatorial hashes from an amino acid sequence

    ...

    Parameters
    ----------
    seq : str
        A sequence of amino acid in their one letter codes
    prot_id : ProteinID
        The identifier for the protein of the passed sequence. If unknown,
        just pass a custom
    """

    # start the pipeline
    aa_vec = get_aa_vector(seq, KIDERA_FACTOR)
    constellation = create_constellation(aa_vec, WINDOW_SIZE, N_PEAKS, overlap=OVERLAP)
    return create_hashes(constellation, prot_id)


def score_prots(
        hashes: Hashes,
        database: Database,
        protein_lookup: ProteinLookup
        ) -> Scores:
    """
    Scores the proteins of the database by their suitability to the given hashes

    ...

    Parameters
    ----------
    hashes : Hashes
        A dictionary of combinatorial hashes generated from a
        sample protein sequence pointing to the index of their occurence and
        the id of the protein
    database: Database
        A dictionary of combinatorial hashes generated from known protein
        sequences pointing to all proteins they occur in, including the position
    protein_lookup: ProteinLookup
        A dictionary of protein identifiers pointing to their description and
        number of combinatorial hashes they have, necessary to calclulate the
        Jaccard Similarity Score

    Returns
    -------
    A list of match protein identifiern and their respective scores.
        The float value is the Jaccard Similarity Index and describes the ratio
        of how many hashes of all possibles for a sequence are found,
        the second value is an additional score to describe how good the found
        hashes fit to the respective protein,
        the first value describes the offset of the matching hashes to their
        original position in the protein sequence
    """

    # stores all found hashes that exist for a known protein
    matches_per_prot: Dict[ProteinID, Tuple[int, int]] = {}

    # find the matches per protein
    for hash_, (sample_index, _) in hashes.items():
        if hash_ in database:
            matching_occurences = database[hash_]

            # for each known protein for the hash, keep the indices of its
            # occurrence in the sample and the known sequence
            for source_index, match_prot_id in matching_occurences:
                if match_prot_id not in matches_per_prot:
                    matches_per_prot[match_prot_id] = []
                matches_per_prot[match_prot_id].append((sample_index, source_index))

    # stores all identifiers of proteins that have matching hashes, pointing
    # to the distance to the match region in their sequence, the score and the
    # Jaccard Similarity Index
    scores_map: Dict[ProteinID, Tuple[int, int, float]] = {}

    # calculate the scores for each protein based on the matching hashes
    for match_prot_id, matches in matches_per_prot.items():

        # calculate the jaccard similarity index as one scoring value
        sample_cardinality = len(hashes)
        _, match_cardinality = protein_lookup[match_prot_id]

        intersection_cardinality = len(matches)
        union_cardinality = sample_cardinality + match_cardinality - intersection_cardinality

        jacc_sim_index = intersection_cardinality / union_cardinality

        # stores all offsets of hashes in the sample to their original position
        # in the match sequence pointing to the number of hashes having this offset
        prot_scores_by_offset: Dict[int, int] = {}

        for sample_index, source_index in matches:
            # calculate the offset
            delta = source_index - sample_index

            # initialize with zero if necessary and count
            if delta not in prot_scores_by_offset:
                prot_scores_by_offset[delta] = 0
            prot_scores_by_offset[delta] += 1

        # find the most related offset, as it describes the biggest matching
        # constellation of found hashes for a protein
        max_ = (0, 0, jacc_sim_index)
        for offset, frequency in prot_scores_by_offset.items():
            if frequency > max_[1]:
                max_ = (offset, frequency, jacc_sim_index)

        # take the frequency of the most related offset as score
        scores_map[match_prot_id] = max_

    # Sort the scores for the user descending by the Jaccard Similarity Index
    # at first level and the calculated score on second level
    scores: Scores = list(sorted(
        scores_map.items(),
        key=lambda x: (x[1][2], x[1][1]),
        reverse=True
    ))

    return scores


if __name__ == '__main__':
    main()
