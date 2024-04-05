from typing import Dict, List, Tuple
from tools import *
from .algorithm import hashes_from_seq
import pickle

Matches = List[Tuple[WindowIndex, WindowIndex]]
MatchesPerProt = Dict[ProteinID, Matches]
ScoresByOffset = Dict[WindowIndex, Score]
ScoresMap = Dict[ProteinID, Tuple[WindowIndex, Score, JSI]]


def find_matches(
        fasta_file: str,
        db_in: str
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
    """

    for (prot_id, description, seq) in Fasta(fasta_file):

        # create the combinatorial hashes for the sequence
        hashes: Hashes = hashes_from_seq(seq, prot_id)

        # load databases
        with open(db_in, 'rb') as f:
            database, protein_lookup = pickle.load(f)

        # calculate the scores for proteins in the database
        scores: Scores = score_prots(hashes, database, protein_lookup)

        top_matches: Scores = get_top_matches(scores)

        print_result(top_matches, protein_lookup)

        print_input_info(prot_id, description, seq, scores, hashes)


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

    matches_per_prot: MatchesPerProt = get_matches_per_prot(hashes, database)

    # stores all identifiers of proteins that have matching hashes, pointing
    # to the distance to the match region in their sequence, the score and the
    # Jaccard Similarity Index
    scores_map: ScoresMap = {}

    # calculate the scores for each protein based on the matching hashes
    for match_prot_id, matches in matches_per_prot.items():

        # calculate the jaccard similarity index as one scoring value
        sample_cardinality: int = len(hashes)
        _, match_cardinality = protein_lookup[match_prot_id]

        intersection_cardinality: int = len(matches)
        union_cardinality: int = sample_cardinality + match_cardinality - intersection_cardinality

        # union cardinality can't be zero, because:
        # sample_cardinality = 0 -> no found hashes -> no matches per prot
        # match_cardinality = 0 -> no found hashes -> not in database -> not in matches per prot
        # -> union_cardinality >= 2
        jacc_sim_index: JSI = intersection_cardinality / union_cardinality

        # count all offsets of hashes in the sample to their original position
        # in the match sequence pointing to the number of hashes having this offset
        prot_scores_by_offset: ScoresByOffset = count_offsets(matches)

        # find the most related offset, as it describes the biggest matching
        # constellation of found hashes for a protein
        max_offset, max_frequency = get_max_offset(prot_scores_by_offset)

        # take the frequency of the most related offset as score
        scores_map[match_prot_id] = (max_offset, max_frequency, jacc_sim_index)

    return sort_scores_descending(scores_map)


def sort_scores_descending(scores_map: ScoresMap) -> Scores:
    # Sort the scores for the user descending by the Jaccard Similarity Index
    # at first level and the calculated score on second level
    scores: Scores = list(sorted(
        scores_map.items(),
        key=lambda x: (x[1][2], x[1][1]),
        reverse=True
    ))

    return scores


def get_max_offset(prot_scores_by_offset: ScoresByOffset) -> Tuple[WindowIndex, Score]:
    max_ = (0, 0)
    for offset, frequency in prot_scores_by_offset.items():
        if frequency > max_[1]:
            max_ = (offset, frequency)

    return max_


def count_offsets(matches: Matches) -> ScoresByOffset:
    prot_scores_by_offset: ScoresByOffset = {}

    for sample_index, source_index in matches:
        # calculate the offset
        delta: WindowIndex = source_index - sample_index

        # initialize with zero if necessary and count
        if delta not in prot_scores_by_offset:
            prot_scores_by_offset[delta] = 0
        prot_scores_by_offset[delta] += 1

    return prot_scores_by_offset


def get_top_matches(scores: Scores) -> Scores:
    # filter to the best
    if scores:
        top_scored = filter(lambda score: scores[0][1][1] == score[1][1] and scores[0][1][2] == score[1][2], scores)
        scores: Scores = list(top_scored)

    return scores


def print_result(matches: Scores, protein_lookup: ProteinLookup):
    # print the matches with description and score
    for prot_index, (_, score, jacc_sim_index) in matches:
        match_prot_description, _ = protein_lookup[prot_index]
        print(f"{prot_index} - {match_prot_description}: Jaccard Index of {jacc_sim_index} : Score of {score}")

    # print first match as most likely match
    if len(matches):
        first_match_id: str = matches[0][0]
        first_match_description: str = protein_lookup[first_match_id][0]
        print("\nSeems to be: %s - %s" % (first_match_id, first_match_description))
    else:
        print("\nNo matches found")


def print_input_info(prot_id: ProteinID, description: str, seq: str, scores: Scores, hashes: Hashes):
    print("\nInput:       %s - %s" % (prot_id, description))
    print(seq)
    if scores:
        _, score, jsi = sorted(scores, key=lambda match: match[0] == prot_id)[-1][1]
    else:
        score, jsi = Score(0), JSI(0)
    print("Input-JSI: %s  -  Input-Score: %s" % (jsi, score))
    print("\nFound hashes: %d" % len(hashes))


def get_matches_per_prot(hashes: Hashes, database: Database) -> MatchesPerProt:
    # stores all found hashes that exist for a known protein
    matches_per_prot: MatchesPerProt = {}

    # find the matches per protein
    for hash_, (sample_index, _) in hashes.items():
        if hash_ in database:
            matching_occurences: List[HashOccurence] = database[hash_]

            # for each known protein for the hash, keep the indices of its
            # occurrence in the sample and the known sequence
            for source_index, match_prot_id in matching_occurences:
                if match_prot_id not in matches_per_prot:
                    matches_per_prot[match_prot_id] = []
                matches_per_prot[match_prot_id].append((sample_index, source_index))

    return matches_per_prot
