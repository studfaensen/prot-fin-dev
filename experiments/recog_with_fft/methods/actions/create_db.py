from tools import *
from .algorithm import hashes_from_seq
from multiprocessing import Pool
from sys import getsizeof as objsize
from os.path import getsize as filesize
import pickle


def create_db(
        prot_file: str,
        db_out: str,
        cpu_count=1,
        **kwargs
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
    db_config = DBConfig(**kwargs)

    fasta = Fasta(prot_file)

    if cpu_count > 1:
        with Pool(cpu_count - 1) as p:
            subprocesses = p.map_async(_process, ((fasta, slice(i, None, cpu_count, db_config)) for i in range(1, cpu_count)))
            database, protein_lookup = _process((fasta, slice(0, None, cpu_count, db_config)))

            for sub_db, sub_lookup in subprocesses.get():
                protein_lookup = {**protein_lookup, **sub_lookup}

                for hash_, occs in sub_db.items():
                    if hash_ not in database:
                        database[hash_] = occs
                    else:
                        database[hash_].extend(occs)

    else:
        database, protein_lookup = _process((fasta, slice(None), db_config))

    # write the databases into files
    with open(db_out, 'wb') as db:
        pickle.dump(DB(database, protein_lookup, db_config), db, pickle.HIGHEST_PROTOCOL)


def _process(args) -> Tuple[Database, ProteinLookup]:
    fasta, slc, db_config = args
    cores = slc.indices(1)[-1]

    database: Database = {}
    protein_lookup: ProteinLookup = {}

    fasta_size = filesize(fasta.file_name)

    # train the database with hashes pointing to their matching proteins
    for i, (prot_id, _, seq) in enumerate(fasta[slc], 1):

        # calculate the combinatorial hashes for the sequence
        hashes: Hashes = hashes_from_seq(seq, prot_id, db_config)

        # save the protein related data
        protein_lookup[prot_id] = (len(seq), len(hashes))

        # sort the hashes into the database by using them as key pointing to
        # their matching proteins
        for hash_, occ in hashes.items():
            if hash_ not in database:
                database[hash_] = []
            database[hash_].append(occ)

        if i % 100 == 0:
            if objsize(pickle.dumps(DB(database, protein_lookup, db_config), pickle.HIGHEST_PROTOCOL)) > 6 * fasta_size / cores:
                raise MemoryError("database too big")

    return database, protein_lookup
