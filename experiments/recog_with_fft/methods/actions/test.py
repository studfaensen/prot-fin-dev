from tools import *
import unittest as ut
import pickle
import sys
import os

from .create_db import *
from .find_matches import *


class TestCreateDB(ut.TestCase):
    protein_file = "test/create_db.fa"
    db_out = "test/create_db.pickle"

    @classmethod
    def setUpClass(cls):
        create_db(cls.protein_file, cls.db_out)

    @classmethod
    def tearDownClass(cls):
        os.remove(cls.db_out)

    def test_create_db(self):
        with open(self.db_out, "rb") as f:
            db = pickle.load(f)
            self.assertTrue(verify_type(db, Tuple[Database, ProteinLookup]), "The loaded pickle-file stores a wrong type")

            db, lookup = db

        self.assertTrue(verify_type(db, Database), "Database has a wrong type")
        self.assertTrue(verify_type(lookup, ProteinLookup), "Database has a wrong type")
        self.assertEqual(len(lookup), 3, "Protein lookup is not complete")

        for hash_, idx_prot_pairs in db.items():

            for pair in idx_prot_pairs:
                self.assertEqual(len(pair), 2, "Protein index pair in DB has different length")
                _, prot = pair
                self.assertIn(prot, lookup, f"'{prot}' not in protein lookup")

        for prot_id, prot_info in lookup.items():
            self.assertEqual(len(prot_info), 2, "Lookup value is differs in length")

            description, hash_count = prot_info
            self.assertFalse(hash_count < 0, "Hash count in protein lookup below zero")


class TestFindMatches(ut.TestCase):
    protein_file = "test/create_db.fa"
    db_in = "test/find_matches.pickle"
    stdout_pipe = "test/find_matches.matches.tmp"

    @classmethod
    def setUpClass(cls):
        create_db(cls.protein_file, cls.db_in)

    @classmethod
    def tearDownClass(cls):
        os.remove(cls.db_in)

    def tearDown(self):
        sys.stdout = sys.__stdout__
        if os.path.exists(self.stdout_pipe):
            os.remove(self.stdout_pipe)

    def test_find_matches(self):
        # Just testing if the function runs without errors

        with open(self.stdout_pipe, "w") as f:
            sys.stdout = f
            self.assertIsNone(find_matches(self.protein_file, self.db_in))

    def test_score_prots(self):
        with open(self.db_in, "rb") as f:
            db, lookup = pickle.load(f)

        hashes = {hash_: (i, "") for i, hash_ in enumerate(db.keys())}
        self.assertTrue(verify_type(hashes, Hashes), "Created test hashes of wrong type")
        scores = score_prots(hashes, db, lookup)
        self.assertTrue(verify_type(scores, Scores), "Returned scores of wrong type")

    def test_get_matches_per_prot(self):
        pass

    def test_count_offsets(self):
        pass

    def get_max_offset(self):
        pass

    def test_sort_scores_descending(self):
        pass

    def test_get_top_matches(self):
        pass

    def test_print_result(self):
        pass

    def test_input_info(self):
        pass


class TestEvaluateProtfin(ut.TestCase):
    ...


class TestSelectSamples(ut.TestCase):
    ...
