# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
import unittest

from genotype_checks.cli import EXPECTED_NAIVE_FAILURES
from genotype_checks.matchers import naive_position_only, naive_string_equal, reference_is_carrier
from genotype_checks.runner import failed_ids, load_cases, run_cases


class CaseFileTest(unittest.TestCase):
    def setUp(self):
        self.cases = load_cases()

    def test_ids_are_unique(self):
        ids = [c["id"] for c in self.cases]
        self.assertEqual(len(ids), len(set(ids)))

    def test_expectations_are_known_values(self):
        for c in self.cases:
            self.assertIn(c["expect"], ("carrier", "not_carrier"), c["id"])

    def test_positions_are_distinct(self):
        keys = [(c["record"]["chrom"], c["record"]["pos"]) for c in self.cases]
        self.assertEqual(len(keys), len(set(keys)))

    def test_imputed_rows_carry_their_own_ref_and_alt(self):
        for c in self.cases:
            g = c["genotype"]
            if g["source"] == "imputed":
                self.assertTrue(g.get("row_ref") and g.get("row_alt"), c["id"])
                self.assertTrue(set(g["alleles"]) <= {g["row_ref"], g["row_alt"]}, c["id"])
            else:
                self.assertEqual(g["source"], "array", c["id"])
                self.assertNotIn("row_ref", g, c["id"])

    def test_every_case_explains_itself(self):
        for c in self.cases:
            self.assertGreater(len(c["why"]), 20, c["id"])

    def test_there_are_both_kinds_of_case(self):
        kinds = {c["expect"] for c in self.cases}
        self.assertEqual(kinds, {"carrier", "not_carrier"})


class MatcherTest(unittest.TestCase):
    def test_reference_rule_passes_every_case(self):
        self.assertEqual(failed_ids(run_cases(reference_is_carrier)), set())

    def test_position_only_matcher_is_caught(self):
        self.assertEqual(failed_ids(run_cases(naive_position_only)),
                         EXPECTED_NAIVE_FAILURES["position-only"])

    def test_string_equal_matcher_fails_only_the_different_anchor_case(self):
        self.assertEqual(failed_ids(run_cases(naive_string_equal)),
                         EXPECTED_NAIVE_FAILURES["string-equal"])

    def test_a_matcher_cannot_change_the_cases(self):
        before = load_cases()

        def vandal(record, genotype):
            record["alt"] = "Z"
            genotype["alleles"] = []
            return True

        run_cases(vandal)
        self.assertEqual(load_cases(), before)

    def test_missing_alleles_are_never_a_carrier(self):
        record = {"chrom": "1", "pos": 1, "ref": "A", "alt": "T"}
        self.assertFalse(reference_is_carrier(record, {"source": "array", "alleles": []}))
        self.assertFalse(reference_is_carrier(record, {"source": "array", "alleles": ["A"]}))

    def test_lower_case_input_is_accepted(self):
        record = {"chrom": "1", "pos": 1, "ref": "a", "alt": "t"}
        self.assertTrue(reference_is_carrier(record, {"source": "array", "alleles": ["a", "t"]}))


if __name__ == "__main__":
    unittest.main()
