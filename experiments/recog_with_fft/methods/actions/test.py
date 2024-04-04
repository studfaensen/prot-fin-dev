import unittest as ut
import pickle
import os

from .create_db import *


class TestCreateDB(ut.TestCase):
    protein_file = "test/create_db.fa"
    db_out = "test/create_db.pickle"

    def tearDown(self):
        os.remove(self.db_out)

    def test_create_db(self):
        create_db(self.protein_file, self.db_out)
        with open(self.db_out, "rb") as f:
            db = pickle.load(f)
            self.assertIsInstance(db, tuple, "The loaded pickle-file doesn't store a tuple")
            self.assertEqual(len(db), 2, "The loaded tuple from pickle-file has a different length")

            db, lookup = db

        self.assertIsInstance(db, dict, "Database is not a dictionary")
        self.assertIsInstance(lookup, dict, "Protein Lookup is not a dictionary")
        self.assertEqual(len(lookup), 3, "Protein lookup is not complete")

        for hash_, idx_prot_pairs in db.items():
            self.assertIsInstance(hash_, int, "Hash in DB not integer")
            self.assertIsInstance(idx_prot_pairs, list, "DB value not list")

            for pair in idx_prot_pairs:
                self.assertIsInstance(pair, tuple, "DB value's list item not tuple")
                self.assertEqual(len(pair), 2, "Protein index pair in DB has different length")

                idx, prot = pair
                self.assertIsInstance(idx, float, "Index in protein index pair not float for index")
                self.assertIsInstance(prot, str, "Index in protein index pair not str for protein identifier")

                self.assertIn(prot, lookup, f"'{prot}' not in protein lookup")

        for prot_id, prot_info in lookup.items():
            self.assertIsInstance(prot_id, str, "Protein identifier not str in protein lookup")
            self.assertIsInstance(prot_info, tuple, "Lookup value not tuple")
            self.assertEqual(len(prot_info), 2, "Lookup value is differs in length")

            description, hash_count = prot_info
            self.assertIsInstance(description, str, "Protein description not str in protein lookup")
            self.assertIsInstance(hash_count, int, "Hash count not integer in protein lookup")
            self.assertFalse(hash_count < 0, "Hash count in protein lookup below zero")


class TestFindMatches(ut.TestCase):
    ...


class TestEvaluateProtfin(ut.TestCase):
    ...


class TestSelectSamples(ut.TestCase):
    ...
