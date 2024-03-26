from typing import List, Dict
from tools import Fasta, ProteinID
import pandas as pd
from numpy import random

Families = List[List[int]]
Samples = Dict[ProteinID, str]

random.seed(1)


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
    map_data: pd.DataFrame = pd.read_csv(mapman, sep="\t").loc[:, ["BINCODE", "IDENTIFIER"]]

    # get families' members' indices
    families: Families = get_groups(map_data)

    samples: Samples = get_samples(families, samples_per_family, map_data)

    print_samples(samples, protein_file)


def print_samples(samples, protein_file):
    for prot_id, _, seq in Fasta(protein_file):
        if prot_id in samples:
            print(">%s  %s\n%s" % (prot_id, samples[prot_id], seq))


def get_samples(families: Families, samples_per_family: int, map_data: pd.DataFrame) -> Samples:
    selected_samples: Samples = {}
    for fam in families:
        chosen_idx = random.choice(fam, samples_per_family, replace=False)
        for fam_id, prot_id in map_data.iloc[chosen_idx].to_numpy():

            fam_id: str = fam_id[1:-1]  # because it has surrounding apostrophes
            prot_id: ProteinID = prot_id[1:-1].upper()

            selected_samples[prot_id] = fam_id

    return selected_samples


def get_groups(map_data: pd.DataFrame) -> Families:
    # The following groups the rows by family:
    # mask: every index that introduces a new family has value of 1 (=True), the others have 0 (False) -> [1, 0, 0, 0, ..., 1, 0, ...]
    # cumsum: The cumulative sum enumerates the families, as the family starts with 1 and the members have value of zero -> [1, 1, ..., 2, 2, ...]
    mask = ~map_data["BINCODE"].str.contains("\\.")
    groups = mask.cumsum()

    # erase the rows that introduce a family
    mask = map_data["IDENTIFIER"] != "''"
    groups = groups[mask]
    map_data = map_data[mask]

    # get the families as lists of the indices that are part of the family
    families: Families = [v.index.tolist() for _, v in map_data.groupby(groups)]

    return families
