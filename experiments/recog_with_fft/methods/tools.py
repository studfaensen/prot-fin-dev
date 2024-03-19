import pandas as pd
import numpy as np
from scipy import signal
from protfin import MAX_ANKER_LENGTH


def get_aa_vector(seq: str, factor: int, normalize=True, file="../../../materials/Amino_Acid_Kidera_Factors.csv") -> np.ndarray:
    # read
    kidera = pd.read_csv(file).loc[:, "A":]

    # normalize
    if normalize:
        min_val = kidera.to_numpy().min()
        kidera.transform(lambda aa_vec: aa_vec + abs(min_val))

    assert len(kidera) > factor or factor < 0, f"get_aa_vector: Kidera factor index has to be between 0-9, got: {factor}"
    sel_factor = kidera.loc[factor, "A":].astype(np.float32)

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
    for s in special_aa:
        sel_factor[s] = sel_factor[special_aa[s]].to_numpy().mean()

    return np.array([sel_factor[i] for i in seq])


def create_constellation(aa_vec: np.ndarray, window_size: int, n_peaks=0, window="boxcar", **kwargs):
    """
    The function carries out a windowed fast fourier transformation,
    often called short time FFT (STFT).

    window_size: size of a window, related to amino acids -> size 2 means two amino acids
    n_peaks: if 0 then all selected

    """
    overlap = kwargs.get("overlap", window_size // 2)

    assert type(n_peaks) is int, "create_constellation: n_peaks must be integer, got: %s" % type(window_size)
    assert type(overlap) is int, "create_constellation: Overlap must be integer, got: %s" % type(overlap)
    assert type(window_size) is int, "create_constellation: Window Size must be integer, got: %s" % type(window_size)
    assert n_peaks > -1, f"create_constellation: n_peaks must be positive, got: {n_peaks}"
    assert overlap > -1, f"create_constellation: Overlap must be positive, got: {overlap}"
    assert window_size > 0, f"create_constellation: Window Size must be true positive, got: {window_size}"

    freq, aa_indexes, stft = signal.stft(
        aa_vec,
        nperseg=window_size,
        noverlap=overlap,
        window=window
    )

    constellation_map = []

    for aa_idx, window in zip(aa_indexes, stft.T):
        # get rid of complex nums
        spectrum = abs(window)
        # Find peaks - these correspond to interesting features
        # prominence=0 includes all peaks, but gives weights weights their prominence
        peaks, props = signal.find_peaks(spectrum, prominence=0)

        # Only want the most prominent peaks
        peaks = sorted(zip(props["prominences"], peaks), reverse=True)
        if n_peaks:
            peaks = peaks[:n_peaks]

        for _, peak in peaks:
            frequency = freq[peak]
            constellation_map.append([aa_idx, frequency])

    return constellation_map


def create_hashes(constellation_map, prot_id=None):
    assert len(constellation_map) > 0
    hashes = {}
    # assume pre-sorted
    # Iterate the constellation
    for idx, (index, freq) in enumerate(constellation_map):
        # Iterate the next 100 pairs to produce the combinatorial hashes
        for other_index, other_freq in constellation_map[idx:]:
            diff = other_index - index
            # If the index difference between the pairs is too small or large
            # ignore this set of pairs
            if diff <= 1 or diff > MAX_ANKER_LENGTH:
                continue

            assert freq < 2 ** 10, f"Frequency {freq} bigger than 1024"
            assert other_freq < 2 ** 10, f"Frequency {other_freq} bigger than 1024"
            assert diff < 2 ** 12, f"Frequency {diff} bigger than 4096"
            # assert int(freq) * int(other_freq) > 0, "frequency 0"

            # Produce a 32 bit hash
            hash = int(freq) | (int(other_freq) << 10) | (int(diff) << 20)
            hashes[hash] = (index, prot_id)
    return hashes
