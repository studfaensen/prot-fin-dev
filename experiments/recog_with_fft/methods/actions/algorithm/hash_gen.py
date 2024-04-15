from tools import *


def create_hashes(
        constellation_map: ConstellationMap,
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

    Returns
    -------
    A dictionary of hashes pointing to the index of their occurence
    """

    hashes: Hashes = {}
    FREQUENCY_BITS = 10
    DIFFERENCE_BITS = 12

    # Iterate through the constellation map
    for idx, (index, freq) in enumerate(constellation_map):
        # Iterate through the next pairs to produce combinatorial hashes
        for other_index, other_freq in constellation_map[idx:]:
            diff = other_index - index

            # If the index difference between the pairs is too small,
            # don't create the hash from them
            if diff <= 1 or diff >= 2 ** DIFFERENCE_BITS:
                continue

            # Produce a 32 bit hash
            hash_: Hash = create_hash((diff, DIFFERENCE_BITS), (other_freq, FREQUENCY_BITS), (freq, FREQUENCY_BITS))
            hashes[hash_] = index
    return hashes


def create_hash(*args) -> Hash:
    hash_: Hash = 0
    bits = 0
    for val, shift in args:
        assert int(val) < 2 ** shift, "%s too big for %s bit" % (val, shift)

        # move bits to the left to make space for the next value
        hash_ <<= shift

        # insert the value's bits
        hash_ |= int(val)

        bits += shift
    assert bits == 32, "Hash exceeds 32 bit"

    return hash_
