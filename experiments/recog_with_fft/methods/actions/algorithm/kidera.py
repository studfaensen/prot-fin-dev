from tools import *
import pandas as pd
import numpy as np
from os import environ as env

KIDERA_FACTOR = ["Helix/bend preference", "Side-chain size", "Extended structure preference", "Hydrophobicity", "Double-bend preference", "Partial specific volume", "Flat extended preference", "Occurrence in alpha region", "pK-C", "Surrounding hydrophobicity"]\
    .index("Hydrophobicity")


def get_aa_vector(
        seq: str,
        factor=KIDERA_FACTOR,
        normalize=True,
        file="../../../materials/Amino_Acid_Kidera_Factors.csv"
        ) -> np.ndarray:
    """
    Transform an amino acid sequence into a vector of floats from the
    Kidera factor table for the selected factor

    ...

    Parameters
    ----------
    seq : str
        The amino acid sequence to be transformed
    factor : int
        The index of the Kidera factor

    Returns
    -------
    A numpy array of 32-bit floats
    """

    # read the table and select the amino acidic data only
    kidera: pd.DataFrame = pd.read_csv(file).loc[:, "A":]

    # normalizing to non-negatives by adding the absolute of the global minimum
    if normalize:
        kidera: pd.DataFrame = normalize_values(kidera)

    # select the specified Kidera factor from the table
    sel_factor: pd.Series = select_factor(kidera, factor)

    # define symbols representing multiple amino acids
    extend_selected_factor(sel_factor)

    # transform the sequence
    transformed_seq = transform_seq(seq, sel_factor)

    return np.array(transformed_seq)


def normalize_values(kidera: pd.DataFrame) -> pd.DataFrame:
    min_val = kidera.to_numpy().min()
    return kidera.transform(lambda aa_vec: aa_vec + abs(min_val))


def select_factor(kidera: pd.DataFrame, factor: int) -> pd.Series:
    return kidera.iloc[factor].astype(np.float32)


def extend_selected_factor(sel_factor: pd.Series):
    special_aa = {
        "X": sel_factor.keys(),  # any aminoacid
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

    # extend the factor data with the multi-representing symbols
    for s in special_aa:
        sel_factor[s] = sel_factor[special_aa[s]].to_numpy().mean()


def transform_seq(seq, sel_factor) -> List[float]:
    transformed_seq: List[float] = []
    for aa in seq:
        value = sel_factor.get(aa)

        # currently 'O' and 'U' are unknown
        if value is None:
            value = 0
            warn(f"No known values for Kidera factors for {aa} -> treating as zero")
        transformed_seq.append(value)

    return transformed_seq
