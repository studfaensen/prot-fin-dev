from tools import *
from .kidera import get_aa_vector
from .constellation import create_constellation
from .hash_gen import create_hashes


def hashes_from_seq(seq: str) -> Hashes:
    """
    Generate the combinatorial hashes from an amino acid sequence

    ...

    Parameters
    ----------
    seq : str
        A sequence of amino acid in their one letter codes
    """

    # start the pipeline
    aa_vec = get_aa_vector(seq)
    constellation = create_constellation(aa_vec)
    return create_hashes(constellation)
