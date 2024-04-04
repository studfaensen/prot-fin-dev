"""
Some essential and useful functions for the algorithm behind prot-fin
"""

from typing import List, Dict, Tuple, Generator, TextIO
from sys import stderr
from tqdm import tqdm
import re

# type aliases
Hash = int
ProteinID = str
Score = int
JSI = float
Hashes = Dict[Hash, Tuple[int, ProteinID]]
Scores = List[Tuple[ProteinID, Tuple[int, Score, JSI]]]
HashOccurence = Tuple[float, ProteinID]
Database = Dict[Hash, List[HashOccurence]]
ProteinLookup = Dict[ProteinID, Tuple[str, int]]
ConstellationMap = List[Tuple[int, float]]


class Fasta:
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
            self.protein_count = count_appearances_in_file("^>", f)

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


def count_appearances_in_file(pattern, file: TextIO):
    count = 0
    file.seek(0)
    while (buffer := file.read(1024 ** 2)):
        count += len(re.findall(pattern, buffer, re.MULTILINE))
    file.seek(0)

    return count


def eprint(*args, **kwargs):
    print(*args, file=stderr, **kwargs)


def warn(*args, **kwargs):
    eprint("WARNING:", *args, **kwargs)
