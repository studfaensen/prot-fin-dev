from tools import *
from .kidera import get_aa_vector
from .constellation import create_constellation
from .hash_gen import create_hashes
from os import environ as env

KIDERA_FACTOR = ["Helix/bend preference", "Side-chain size", "Extended structure preference", "Hydrophobicity", "Double-bend preference", "Partial specific volume", "Flat extended preference", "Occurrence in alpha region", "pK-C", "Surrounding hydrophobicity"]\
    .index("Hydrophobicity")

# parameters for the STFT
WINDOW_SIZE = int(env.get("WINDOW_SIZE", 30))
WINDOW_TYPE = env.get("WINDOW_TYPE", ["boxcar", "triang", "blackman", "hamming", "hann", "bartlett", "flattop", "parzen", "bohman", "blackmanharris", "nuttall", "barthann", "cosine", "exponential", "tukey", "taylor", "lanczos"]\
                                     [0])
OVERLAP = int(env.get("OVERLAP", 15))
N_PEAKS = int(env.get("N_PEAKS", 0))  # 0 means all


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
    constellation = create_constellation(aa_vec, WINDOW_SIZE, N_PEAKS, WINDOW_TYPE, overlap=OVERLAP)
    return create_hashes(constellation, prot_id)
