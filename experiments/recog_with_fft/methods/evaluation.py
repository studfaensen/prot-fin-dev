"""
This module is to do some statistic analysis on the protfin results
"""

from protfin import find_match, DB_DEFAULT
from tqdm import tqdm
from tools import Fasta
import pandas as pd
from sys import stdin
from numpy import random
import argparse
import re

random.seed(1)


def main():
    """
    This is where the command line interface is defined to interpret
    the passed arguments to execute the respective functions.
    """

    parser = argparse.ArgumentParser(
        prog='evaluation',
        description=__doc__,
    )
    sub_commands = parser.add_subparsers(required=True)

    # protfin.py ... | evaluation.py eval
    eval_parser = sub_commands.add_parser("eval", help="Evaluate protfin find-match output from stdin")
    eval_parser.add_argument("protfin-out-file")
    eval_parser.set_defaults(func=lambda args: evaluate_protfin(getattr(args, "protfin-out-file")))

    # evaluation.py select-samples
    eval_parser = sub_commands.add_parser("select-samples", help="Select samples from reference")
    eval_parser.add_argument("mapman-file")
    eval_parser.add_argument("protein-file")
    eval_parser.add_argument("-s", "--samples-per-family", default=1, type=int)
    eval_parser.set_defaults(func=lambda args:
                             select_samples(getattr(args, "mapman-file"),
                                            getattr(args, "protein-file"),
                                            args.samples_per_family
                                            )
                             )

    args = parser.parse_args()
    args.func(args)


def evaluate_protfin(protfin_out_file: str):
    """
    Summarizes the output of prot-fin's find-match output.
    For each sample protein, it collects the following information into csv:
      - sample protein id                         -> Sample_ID
      - number of found matches                   -> First_Match_Count
      - occurs the sample in the found matches    -> Sample_In_First_Matches
      - sequence length of sample                 -> Sequence_Length
      - created hashes for sample                 -> Sample_Hashes
      - Jaccard Similarity Score of sample        -> Sample_JSI
      - Score of sample                           -> Sample_Score
      - Jaccard Similarity Score of found matches -> Sample_JSI
      - Score of found matches                    -> Sample_Score

    The data is printed to stdout

    ...

    Parameters
    ----------
    protfin_out_file : str
        Path to the file the output of protfin was written to
    """
    with open(protfin_out_file) as f:

        # just counting number proteins used as protfin input
        result_count = 0
        while (buffer := f.read(1024 ** 2)):
            result_count += len(re.findall("Seems to be", buffer))
        f.seek(0)

        # data will be collected into a dataframe
        evaluation = pd.DataFrame(columns=["Sample_ID", "First_Match_Count", "Sample_In_First_Matches", "Sequence_Length", "Sample_Hashes", "Sample_JSI", "Sample_Score", "Top_JSI", "Top_Score"])
        for _ in tqdm(range(result_count)):
            matches = {}

            # iterate through the lines of matches and extract their information
            while (match := f.readline()) != "\n":
                prot_id, jacc_index, score = re.findall("(.*) - .*Jaccard Index of (\d+\.\d+) : Score of (\d+)", match)[0]
                matches[prot_id] = (jacc_index, score)

            # now crawl the few lines below a protfin result to collect data
            # about the input protein
            f.readline()  # Seems to be...
            f.readline()  #
            input_sample = re.findall(".*Input: +(.*) - .*", f.readline())[0]  # Input:  ...
            seq = f.readline()[:-1]  # sequence line
            jsi, score = re.findall("JSI: (\d+\.\d).*Score: (\d+)", f.readline())[0]  # Input-JSI ...
            f.readline()  #
            hashes = re.findall("Found hashes: (\d+)", f.readline())[0]  # Found hashes: ...

            # insert the data into the dataframe
            evaluation.loc[len(evaluation.index)] = (
                input_sample,  # Sample_ID
                len(matches),  # First_Match_Count
                input_sample in matches,  # Sample_In_First_Matches
                len(seq),  # Sequence_Length
                hashes,  # Sample_Hashes
                jsi,  # Sample_JSI
                score,  # Sample_Score
                list(matches.values())[0][0],  # Top_JSI
                list(matches.values())[0][1]  # Top_Score
            )

        # write to stdout
        print(evaluation.to_csv(index=False), end="")


def select_samples(mapman: str, protein_file: str, samples_per_family: int):
    """
    Analyzes the mapman bins to select random functional different
    proteins as train data.

    The data is printed to stdout

    ...

    Parameters
    ----------
    mapman : str
        Path to the mapman result for the proteins.
    protein_file : str
        Path to the FASTA formatted file storing the amino acid sequences for
        the proteins used for mapman.
    samples_per_family : int
        Number of samples that will be selected per each family.
    """

    # read the mapman reference
    map_data = pd.read_csv(mapman, sep="\t").loc[:, ["BINCODE", "IDENTIFIER"]]

    # now erase unnecessary rows
    mask = map_data['BINCODE'].str.contains('\\.') & (map_data['IDENTIFIER'] == "''")  # all entries defining a subfamily
    map_data = map_data[~mask]  # keep only rows with Identifiern or rows introducing a new family -> BINCODE has no dot
    map_data.reset_index(drop=True, inplace=True)

    # The following groups the rows by family:
    # mask: every index that introduces a new family has value of 0 (=False), the others have 1 (True) -> [0, 1, 1, 1, ...]
    # ~mask: invert the values -> [1, 0, 0, 0, ...]
    # cumsum: The cumulative sum enumerates the families, as the family starts with 1 and the members have value of zero -> [1, 1, 1, 1, ...]
    mask = map_data["BINCODE"].str.contains("\\.")
    groups = (~mask).cumsum()

    # erase the rows that introduce one family, as their indices in mask have the value of 0 (=False) -> [2, 2, 2, ...]
    groups = groups[mask]

    # get the families as lists of the indices that are part of the family
    families = [v.index.tolist() for _, v in map_data[mask].groupby(groups)]

    selected_samples = {}

    for fam in families:
        chosen_idx = random.choice(fam, samples_per_family, replace=False)
        for fam_id, prot_id in map_data.iloc[chosen_idx].to_numpy():
            fam_id = fam_id[1:-1]  # because it has surrounding apostrophes
            prot_id = prot_id[1:-1].upper()
            selected_samples[prot_id] = fam_id

    for prot_id, _, seq in Fasta(protein_file):
        if prot_id in selected_samples:
            print(">%s  %s\n%s" % (prot_id, selected_samples[prot_id], seq))


if __name__ == '__main__':
    main()
