from typing import Type, Any
from functools import reduce
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
        cls.assertIsNone(cls, create_db(cls.protein_file, cls.db_out))

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

    def create_valid(self, ty: Type, value: Any) -> Type:
        type_name = str(locals()["ty"])
        self.assertTrue(verify_type(value, ty), "Value of wrong type, expected: " + type_name)
        return value

    @classmethod
    def setUpClass(cls):
        cls.assertIsNone(cls, create_db(cls.protein_file, cls.db_in))
        with open(cls.db_in, "rb") as f:
            cls.db, cls.lookup = pickle.load(f)

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

        hashes: Hashes = self.create_valid(
            Hashes,
            {hash_: (WindowIndex(i), ProteinID("")) for i, hash_ in enumerate(db.keys())}
        )
        scores: Scores = self.create_valid(
            Scores,
            score_prots(hashes, db, lookup)
        )
        self.assertEqual(db, self.db, "Database has changed")
        self.assertEqual(lookup, self.lookup, "Database has changed")

    def test_get_matches_per_prot(self):
        hashes: Hashes = self.create_valid(
            Hashes,
            {Hash(123): (WindowIndex(1), ProteinID(""))}
        )
        proteins = [(WindowIndex(0), ProteinID(i)) for i in range(10)]
        db = self.create_valid(
            Database,
            {Hash(123): proteins}
        )
        matches: MatchesPerProt = self.create_valid(
            MatchesPerProt,
            get_matches_per_prot(hashes, db)
        )
        self.assertEqual(len(matches), len(proteins), "Returned matches not complete")
        for idx, prot_id in proteins:
            self.assertIn(prot_id, matches, f"Protein identifier '{prot_id}' missing")
            original_matches = [(list(hashes.values())[0][0], idx)]
            self.assertEqual(
                matches[prot_id],
                original_matches,
                f"Found matches for protein id '{prot_id}' wrong"
            )

    def test_count_offsets(self):
        matches_per_offset = {
            WindowIndex(-1): [
                (WindowIndex(10), WindowIndex(11)),
                (WindowIndex(100), WindowIndex(101))
            ],
            WindowIndex(0): [
                (WindowIndex(11), WindowIndex(11)),
            ],
            WindowIndex(1): [
                (WindowIndex(12), WindowIndex(11)),
                (WindowIndex(2), WindowIndex(1)),
                (WindowIndex(102), WindowIndex(101))
            ]
        }
        matches: Matches = self.create_valid(
            Matches,
            reduce(lambda a, b: a + b, matches_per_offset.values(), [])
        )
        scored_offsets: ScoresByOffset = self.create_valid(
            ScoresByOffset,
            count_offsets(matches)
        )
        self.assertEqual(len(scored_offsets), len(matches_per_offset), "Returned scored offsets not complete")
        for offset, score in scored_offsets.items():
            self.assertIn(offset, matches_per_offset, f"Offset '{offset}' missing")
            self.assertEqual(
                scored_offsets[offset],
                score,
                f"Returned score for offset '{offset}' wrong"
            )

    def test_get_max_offset(self):
        scores = list(range(5)) + list(range(3, 0, -1))
        scored_offsets: ScoresByOffset = self.create_valid(
            ScoresByOffset,
            {WindowIndex(i): score for i, score in enumerate(scores)}
        )
        max_offset: (WindowIndex, Score) = self.create_valid(
            Tuple[WindowIndex, Score],
            get_max_offset(scored_offsets)
        )
        max_score = sorted(scores)[-1]
        self.assertEqual(max_offset, (WindowIndex(scores.index(max_score)), max_score))

    def test_sort_scores_descending(self):
        scores_map: ScoresMap = self.create_valid(
            ScoresMap,
            {ProteinID(i): (WindowIndex(i), Score(1 + i//2), JSI(i//3 * .1)) for i in range(30)}
        )
        descending: Scores = self.create_valid(
            Scores,
            sort_scores_descending(scores_map)
        )
        prev_score, prev_jsi = float("+inf"), float("+inf")
        for _, (_, score, jsi) in descending:
            self.assertTrue(score <= prev_score, f"Score {score} is greater than its previous {prev_score}")
            self.assertTrue(jsi <= prev_jsi, f"JSI {jsi} is greater than its previous {prev_jsi}")
            prev_score, prev_jsi = score, jsi

    def test_get_top_matches(self):
        scores: Scores = self.create_valid(
            Scores,
            [(ProteinID(i), (WindowIndex(i), Score(1 + i//2), JSI(i//3 * .1))) for i in range(29, 0, -1)]
        )
        top: Scores = self.create_valid(
            Scores,
            get_top_matches(scores)
        )
        self.assertNotEqual(len(top), len(scores), "Length of input scores has changed")
        self.assertEqual(top, scores[:2], "Top scores not as expected")

    def test_print_result(self):
        proteins = [ProteinID(i) for i in range(10)]
        scores_map: ScoresMap = self.create_valid(
            ScoresMap,
            {p: (WindowIndex(i), Score(1 + i//2), JSI(i//3 * .1)) for i, p in enumerate(proteins)}
        )
        scores: Scores = self.create_valid(
            Scores,
            list(scores_map.items())
        )
        lookup: ProteinLookup = self.create_valid(
            ProteinLookup,
            {p: (f"description of {p}", 1) for p in proteins}
        )

        with open(self.stdout_pipe, "w") as f:
            sys.stdout = f
            print_result(scores, lookup)
        with open(self.stdout_pipe, "r") as f:
            for line, prot_id in zip(f, proteins):
                self.assertEqual(
                    line,
                    f"{prot_id} - description of {prot_id}: Jaccard Index of {scores_map[prot_id][2]} : Score of {scores_map[prot_id][1]}\n",
                    "Printed match output different"
                )

    def test_input_info(self):
        prot_id = ProteinID("prot_id")
        description = "very cool"
        seq = "MAAGCSTTTYRAGGG"
        scores: Scores = self.create_valid(
            Scores,
            [(prot_id, (WindowIndex(1), Score(10), JSI(.5)))]
        )
        hashes: Hashes = self.create_valid(
            Hashes,
            {Hash(i): (WindowIndex(i), ProteinID(0)) for i in range(19)}
        )
        with open(self.stdout_pipe, "w") as f:
            sys.stdout = f
            print_input_info(prot_id, description, seq, scores, hashes)
        with open(self.stdout_pipe, "r") as f:
            content = f.read()
            expected = f"""
Input:       {prot_id} - {description}
{seq}
Input-JSI: {scores[0][1][2]}  -  Input-Score: {scores[0][1][1]}

Found hashes: {len(hashes)}
"""
            self.assertEqual(
                content,
                expected,
                "Printed input info output different"
            )


class TestEvaluateProtfin(ut.TestCase):
    ...


class TestSelectSamples(ut.TestCase):
    ...
