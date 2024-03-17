#! python3
from typing import List, Dict, Tuple, TextIO
import numpy as np
import pandas as pd
import pickle
import argparse

# %%
from create_constellations import create_constellation

SAMPLERATE = 44100


def main():
    parser = argparse.ArgumentParser(
        prog='prot-fin',
        description='protein recognition by transforming their sequences into music',
    )
    sub_commands = parser.add_subparsers(required=True)

    create_db_parser = sub_commands.add_parser("create-db", help="Create Database")
    create_db_parser.add_argument("fasta-file")
    create_db_parser.set_defaults(func=lambda a: create_db(getattr(a, "fasta-file")))

    find_match_parser = sub_commands.add_parser("find-match", help="Find Matches for Proteins")
    find_match_parser.add_argument("fasta-file")
    find_match_parser.set_defaults(func=lambda a: find_match(getattr(a, "fasta-file")))

    args = parser.parse_args()
    args.func(args)


def find_match(fasta_file: str):
    chord_map = get_aa_chords()
    with open(fasta_file) as f:
        for (prot_id, description, seq) in iter_fasta(f):
            hashes = hashes_from_seq(seq, chord_map)

            database = pickle.load(open('database.pickle', 'rb'))
            scores = score_songs(hashes, database)

            song_index_lookup = pickle.load(open("song_index.pickle", "rb"))
            for song_index, score in scores:
                match_prot_description = song_index_lookup[song_index]
                print(f"{song_index} - {match_prot_description}: Score of {score[1]} at {score[0]}")

            print("\nSeems to be: %s - %s" % (scores[0][0], song_index_lookup[scores[0][0]]))
            print("\nInput:       %s - %s\n             %s" % (prot_id, description, seq))


def score_songs(hashes, database):
    matches_per_song = {}
    for hash_, (sample_time, _) in hashes.items():
        if hash_ in database:
            matching_occurences = database[hash_]
            for source_time, song_index in matching_occurences:
                if song_index not in matches_per_song:
                    matches_per_song[song_index] = []
                matches_per_song[song_index].append((sample_time, source_time))

    scores = {}
    for song_index, matches in matches_per_song.items():
        song_scores_by_offset = {}
        for sample_time, source_time in matches:
            delta = source_time - sample_time
            if delta not in song_scores_by_offset:
                song_scores_by_offset[delta] = 0
            song_scores_by_offset[delta] += 1

        max_ = (0, 0)
        for offset, score in song_scores_by_offset.items():
            if score > max_[1]:
                max_ = (offset, score)
        scores[song_index] = max_

    # Sort the scores for the user
    scores = list(sorted(scores.items(), key=lambda x: x[1][1], reverse=True))

    return scores


def create_db(prot_file: str):
    song_index = {}
    database: Dict[str, List[Tuple[int, int]]] = {}
    chord_map = get_aa_chords()
    with open(prot_file) as f:
        for prot_id, description, seq in iter_fasta(f):
            if (proc := len(song_index)) == 200:
                break
            if proc % 20 == 0:
                print("%d%%" % int(proc/2))
            song_index[prot_id] = description
            hashes = hashes_from_seq(seq, chord_map, prot_id)
            for hash_, time_index_pair in hashes.items():
                if hash_ not in database:
                    database[hash_] = []
                database[hash_].append(time_index_pair)
    # %%
    with open("database.pickle", 'wb') as db:
        pickle.dump(database, db, pickle.HIGHEST_PROTOCOL)
    with open("song_index.pickle", 'wb') as songs:
        pickle.dump(song_index, songs, pickle.HIGHEST_PROTOCOL)


def create_hashes(constellation_map, song_id=None):
    upper_frequency = 23_000
    frequency_bits = 10
    hashes = {}
    # assume pre-sorted
    # Iterate the constellation
    for idx, (time, freq) in enumerate(constellation_map):
        # Iterate the next 100 pairs to produce the combinatorial hashes
        for other_time, other_freq in constellation_map[idx : idx + 100]:
            diff = other_time - time
            # If the time difference between the pairs is too small or large
            # ignore this set of pairs
            if diff <= 1 or diff > 10:
                continue

            # Place the frequencies (in Hz) into a 1024 bins
            freq_binned = freq / upper_frequency * (2 ** frequency_bits)
            other_freq_binned = other_freq / upper_frequency * (2 ** frequency_bits)

            # Produce a 32 bit hash
            hash = int(freq_binned) | (int(other_freq_binned) << 10) | (int(diff) << 20)
            hashes[hash] = (time, song_id)
    return hashes


def hashes_from_seq(seq: str, chord_map, song_id=None):
    audio = seq_to_audio(seq, chord_map)
    constellation = create_constellation(audio, SAMPLERATE)
    return create_hashes(constellation, song_id)


def get_wave(freq, duration=0.05):
    amplitude = 4096
    t = np.linspace(0, duration, int(SAMPLERATE * duration))
    wave = amplitude * np.sin(2 * np.pi * freq * t)

    return wave


def seq_to_audio(seq: str, chord_map) -> np.ndarray:
    return np.concatenate([chord_map[aa] for aa in seq])


def get_aa_chords(file="../../../materials/Amino_Acid_Kidera_Factors.csv") -> Dict[str, np.ndarray]:
    kidera = pd.read_csv(file).loc[:, "A":]

    base_freq = 261.63  # Frequency of Note C4
    base_chord = pd.Series(get_wave(base_freq * pow(2, (i / 12))) for i in range(len(kidera.index)))

    special_aa = {
        "X": list(kidera),  # any aminoacid
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
        kidera[s] = kidera[special_aa[s]].to_numpy().mean(axis=1)

    chord_map = {}
    for aa in kidera:
        chord = sum(base_chord * (5 + kidera[aa] * 2))
        chord *= 16300 / np.max(chord)
        chord_map[aa] = chord.astype(np.int16)

    return chord_map


def iter_fasta(fasta_file: TextIO) -> Tuple[str, str]:
    while True:
        prot_desc = fasta_file.readline()

        if "\n" not in prot_desc:  # -> EOF
            break
        if prot_desc[0] == ">":
            prot_id, _, description = prot_desc.split(" ", 2)
            seq = fasta_file.readline()
            if seq[-1] == "\n":
                seq = seq[:-1]

            yield prot_id[1:], description[:-1], seq


if __name__ == '__main__':
    main()
