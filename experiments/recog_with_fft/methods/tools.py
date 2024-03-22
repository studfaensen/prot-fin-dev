"""
Some essential and useful functions for the algorithm behind prot-fin
"""

from typing import List, Dict, Tuple, Generator
from scipy import signal
from sys import stderr
from tqdm import tqdm
import pandas as pd
import numpy as np
import re

# type aliases
Hash = int
ProteinID = str
Hashes = Dict[Hash, Tuple[int, ProteinID]]
Scores = List[Tuple[ProteinID, Tuple[int, int, float]]]
Database = Dict[Hash, List[Tuple[int, ProteinID]]]
ProteinLookup = Dict[ProteinID, Tuple[str, int]]
ConstellationMap = List[Tuple[int, float]]


class Fasta():
    """
    A class used for convenient iteration over a FASTA file's contents.

    ...

    Attributes
    ----------
    file_name : str
        The name of the FASTA formatted file
    protein_count : int
        The number of sequences stored in the FASTA file
    """
    def __init__(self, file_name: str):
        """
        Parameters
        ----------
        file_name : str
            The name of the FASTA formatted file
        """
        with open(file_name) as f:
            self.file_name = file_name
            self.protein_count = 0

            # count sequences
            while (buffer := f.read(1024 ** 2)):
                self.protein_count += len(re.findall('^>', buffer, re.MULTILINE))
            f.seek(0)

            # validate ... TODO

    def __len__(self):
        return self.protein_count

    def __iter__(self) -> Generator[Tuple[ProteinID, str, str], None, None]:
        with open(self.file_name) as f:

            # create a progress bar and iterate over the FASTA file
            for _ in tqdm(range(len(self))):

                # find line of sequence description
                while True:
                    prot_desc = f.readline()
                    if prot_desc[0] == ">":
                        break

                # extract information from describing line
                prot_id, _, description = prot_desc.split(" ", 2)
                seq = f.readline()
                if seq[-1] == "\n":
                    seq = seq[:-1]

                # yield the extracted values, remove '>' from identifier and
                # '\n' from description
                yield prot_id[1:], description[:-1], seq


def get_aa_vector(
        seq: str,
        factor: int,
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
    kidera = pd.read_csv(file).loc[:, "A":]

    # normalizing to non-negatives by adding the absolute of the global minimum
    if normalize:
        min_val = kidera.to_numpy().min()
        kidera = kidera.transform(lambda aa_vec: aa_vec + abs(min_val))

    # select the specified Kidera factor from the table
    assert len(kidera) > factor or factor < 0, f"get_aa_vector: Kidera factor index has to be between 0-9, got: {factor}"
    sel_factor = kidera.loc[factor, "A":].astype(np.float32)

    # define symbols representing multiple amino acids
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

    # transform the sequence
    transformed_seq = []
    for aa in seq:
        value = sel_factor.get(aa)

        # currently 'O' and 'U' are unknown
        if value is None:
            value = 0
            warn(f"No known values for Kidera factors for {aa} -> treating as zero")
        transformed_seq.append(value)

    return np.array(transformed_seq)


def create_constellation(
        aa_vec: np.ndarray,
        window_size: int,
        n_peaks=0,
        window="boxcar",
        **kwargs
        ) -> ConstellationMap:
    """
    The function carries out a windowed fast fourier transformation,
    often called short time FFT (STFT), on the given vector and creates a list
    of frequencies, found by the STFT, and the index of the window they are found

    ...

    Parameter
    ---------
    aa_vec : np.ndarray
        An array of floats representing an amino acid sequence
    window_size : int
        size of a window, related to amino acids -> size 2 means two amino acids
    n_peaks : int, optional
        the number of frequency peaks that are selected from the STFT result,
        if 0 then all are selected (defaults to 0)
    window : str
        The window type for the STFT, read https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.get_window.html#scipy.signal.get_window
        for supported ones (defaults to "boxcar")
    **kwargs:
        overlap : int
            the overlap between two windows during doing the STFT (defaults to half the window size)

    Returns
    -------
     A ConstellationMap, a list of coordinates, meaning the pairs of
     window index and prominent frequency peak of the window
    """

    # extract overlap from keyword arguments
    overlap = kwargs.get("overlap", window_size // 2)

    assert type(n_peaks) is int, "create_constellation: n_peaks must be integer, got: %s" % type(n_peaks)
    assert type(overlap) is int, "create_constellation: Overlap must be integer, got: %s" % type(overlap)
    assert type(window_size) is int, "create_constellation: Window Size must be integer, got: %s" % type(window_size)
    assert n_peaks > -1, f"create_constellation: n_peaks must be positive, got: {n_peaks}"
    assert overlap > -1, f"create_constellation: Overlap must be positive, got: {overlap}"
    assert window_size > 0, f"create_constellation: Window Size must be true positive, got: {window_size}"

    # adjust window size and overlap if invalid
    if len(aa_vec) < window_size:
        window_size = len(aa_vec)
    if overlap >= window_size:
        overlap = window_size - 1

    # executing the STFT
    frequencies, window_indexes, stft = signal.stft(
        aa_vec,
        nperseg=window_size,
        noverlap=overlap,
        window=window
    )

    constellation_map: ConstellationMap = []

    # find and collect the most prominent frequencies from STFT per window
    for window_idx, amplitudes in zip(window_indexes, stft.T):

        # get rid of complex values to make them comparable
        spectrum = abs(amplitudes)

        # find peaks
        # prominence=0 includes all peaks, but weights their prominence as well
        peaks, props = signal.find_peaks(spectrum, prominence=0)

        # Only want the most prominent peaks
        peaks = sorted(zip(props["prominences"], peaks), reverse=True)
        if n_peaks:
            peaks = peaks[:n_peaks]

        for _, peak in peaks:
            frequency = frequencies[peak]
            constellation_map.append([window_idx, frequency])

    return constellation_map


def create_hashes(
        constellation_map: ConstellationMap,
        prot_id: ProteinID
        ) -> Hashes:
    """
    Creates combinatorial Hashes from a constellation map for efficient
    database searches

    ...

    Parameters
    ----------
    constellation_map : ConstellationMap
        A list of coordinates, index-frequency pairs, the result of the STFT
        in create_constellation. It is assumed as pre-sorted by index
    prot_id : ProteinID
        The identifier for the protein the hashes are generated for. If it is
        unknown, just pass a custom

    Returns
    -------
    A dictionary of hashes pointing to the index of their occurence and the
    protein they are created for
    """

    hashes: Hashes = {}

    # Iterate through the constellation map
    for idx, (index, freq) in enumerate(constellation_map):
        FREQUENCY_BITS = 10
        DIFFERENCE_BITS = 12

        # Iterate through the next pairs to produce combinatorial hashes
        for other_index, other_freq in constellation_map[idx:]:
            diff = other_index - index

            # If the index difference between the pairs is too small,
            # don't create the hash from them
            if diff <= 1 or diff > 2 ** DIFFERENCE_BITS:
                continue

            assert freq < 2 ** FREQUENCY_BITS, \
                "Frequency %s bigger than %s" % (freq, 2 ** FREQUENCY_BITS)
            assert other_freq < 2 ** FREQUENCY_BITS, \
                "Frequency %s bigger than %s" % (freq, 2 ** FREQUENCY_BITS)
            assert diff < 2 ** DIFFERENCE_BITS, \
                "Distance %s bigger than %s" % (diff, 2 ** DIFFERENCE_BITS)

            # Produce a 32 bit hash
            hash_ = int(freq) | (int(other_freq) << FREQUENCY_BITS) | (int(diff) << 2 * FREQUENCY_BITS)
            hashes[hash_] = (index, prot_id)
    return hashes


def eprint(*args, **kwargs):
    print(*args, file=stderr, **kwargs)


def warn(*args, **kwargs):
    eprint("WARNING:", *args, **kwargs)
