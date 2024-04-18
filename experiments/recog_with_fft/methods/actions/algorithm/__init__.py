from tools import *
from .kidera import get_aa_vector
from .constellation import create_constellation
from .hash_gen import create_hashes


def hashes_from_seq(seq: str, prot_id: str) -> Hashes:
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
    aa_vec = get_aa_vector(seq)
    constellation = create_constellation(aa_vec)
    return create_hashes(constellation, prot_id)
