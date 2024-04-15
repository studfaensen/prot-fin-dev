from tools import *
from .algorithm import hashes_from_seq
import pickle


def create_db(
        prot_file: str,
        db_out: str
        ):
    """
    Creates a database of combinatorial hashes for all protein sequences in the
    given FASTA file and also a lookup database of protein identifiers pointing
    to their description and number of hashes

    ...

    Parameters
    ----------
    prot_file : str
        The path to the FASTA formatted file containing all train sequences of
        known protein sequences
    db_out : str
        Name of the file to write the databases to
    """

    protein_lookup: ProteinLookup = {}
    database: Database = {}

    # train the database with hashes pointing to their matching proteins
    for prot_id, description, seq in Fasta(prot_file):

        # calculate the combinatorial hashes for the sequence
        hashes: Hashes = hashes_from_seq(seq)

        # save the protein related data
        protein_lookup[prot_id] = (description, len(hashes))

        # sort the hashes into the database by using them as key pointing to
        # their matching proteins
        for hash_, index in hashes.items():
            if hash_ not in database:
                database[hash_] = []
            database[hash_].append((index, prot_id))

    # write the databases into files
    with open(db_out, 'wb') as db:
        pickle.dump((database, protein_lookup), db, pickle.HIGHEST_PROTOCOL)
