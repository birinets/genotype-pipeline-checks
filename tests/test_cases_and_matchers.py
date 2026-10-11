# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
import copy
import unittest

from genotype_checks.cli import EXPECTED_NAIVE_FAILURES
from genotype_checks.matchers import KNOWN_BAD, reference_matches
from genotype_checks.runner import failed_ids, load_cases, run_cases


class CaseFileTest(unittest.TestCase):
    """The case file is read by other languages too, so its shape is a contract."""

    def setUp(self):
        self.cases = load_cases()

    def test_ids_are_unique(self):
        ids = [c["id"] for c in self.cases]
        self.assertEqual(len(ids), len(set(ids)))

    def test_expectations_are_known_values(self):
        for c in self.cases:
            self.assertIn(c["expect"], ("match", "no_match"), c["id"])

    def test_positions_are_distinct(self):
        keys = [(c["record"]["chrom"], c["record"]["pos"]) for c in self.cases]
        self.assertEqual(len(keys), len(set(keys)))

    def test_imputed_rows_carry_their_own_ref_and_alt(self):
        for c in self.cases:
            g = c["genotype"]
            if g["source"] == "imputed":
                self.assertTrue(g.get("row_ref") and g.get("row_alt"), c["id"])
                self.assertTrue(set(g["alleles"]) <= {g["row_ref"], g["row_alt"]}, c["id"])

    def test_array_genotypes_are_two_single_letters_without_row_alleles(self):
        for c in self.cases:
            g = c["genotype"]
            if g["source"] == "array":
                self.assertEqual([len(a) for a in g["alleles"]], [1, 1], c["id"])
                self.assertNotIn("row_ref", g, c["id"])
                self.assertNotIn("row_alt", g, c["id"])

    def test_every_case_explains_itself(self):
        for c in self.cases:
            self.assertGreater(len(c["why"]), 20, c["id"])

    def test_there_are_both_kinds_of_case(self):
        self.assertEqual({c["expect"] for c in self.cases}, {"match", "no_match"})


class MatcherTest(unittest.TestCase):
    def test_reference_rule_passes_every_case(self):
        self.assertEqual(failed_ids(run_cases(reference_matches)), set())

    def test_each_known_bad_rule_fails_exactly_its_cases(self):
        self.assertEqual(set(KNOWN_BAD), set(EXPECTED_NAIVE_FAILURES))
        for name, rule in KNOWN_BAD.items():
            with self.subTest(rule=name):
                self.assertEqual(failed_ids(run_cases(rule)), EXPECTED_NAIVE_FAILURES[name])

    def test_every_no_match_case_catches_at_least_one_known_bad_rule(self):
        caught = set().union(*EXPECTED_NAIVE_FAILURES.values())
        for c in load_cases():
            if c["expect"] == "no_match":
                self.assertIn(c["id"], caught, "a trap that no known mistake falls into proves nothing")

    def test_a_matcher_cannot_change_the_cases(self):
        def vandal(record, genotype):
            record["alt"] = "Z"
            genotype["alleles"] = []
            return True

        cases = load_cases()
        before = copy.deepcopy(cases)
        run_cases(vandal, cases)
        self.assertEqual(cases, before)

    def test_identical_ref_and_alt_never_match(self):
        self.assertFalse(reference_matches({"ref": "A", "alt": "A"}, {"source": "array", "alleles": ["A", "A"]}))
        self.assertFalse(reference_matches(
            {"ref": "A", "alt": "A"},
            {"source": "imputed", "alleles": ["C", "C"], "row_ref": "C", "row_alt": "C"}))

    def test_a_row_with_only_one_of_ref_and_alt_never_matches(self):
        record = {"ref": "A", "alt": "C"}
        self.assertFalse(reference_matches(record, {"source": "imputed", "alleles": ["A", "C"], "row_ref": "A"}))
        self.assertFalse(reference_matches(record, {"source": "imputed", "alleles": ["A", "C"], "row_alt": "C"}))

    def test_missing_alleles_never_match(self):
        record = {"chrom": "1", "pos": 1, "ref": "A", "alt": "C"}
        self.assertFalse(reference_matches(record, {"source": "array", "alleles": []}))
        self.assertFalse(reference_matches(record, {"source": "array", "alleles": ["A"]}))

    def test_lower_case_input_is_accepted(self):
        record = {"chrom": "1", "pos": 1, "ref": "a", "alt": "c"}
        self.assertTrue(reference_matches(record, {"source": "array", "alleles": ["a", "c"]}))

    def test_documented_limit_letters_on_the_other_strand_are_not_resolved(self):
        # The README states this limit. If strand handling is ever added, this test and that text change together.
        record = {"chrom": "1", "pos": 1, "ref": "A", "alt": "G"}
        self.assertFalse(reference_matches(record, {"source": "array", "alleles": ["T", "C"]}))


if __name__ == "__main__":
    unittest.main()
