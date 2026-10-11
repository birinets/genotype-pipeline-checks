# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
import unittest

from genotype_checks.alleles import same_allele, trim_alleles


class TrimAllelesTest(unittest.TestCase):
    def test_snv_is_unchanged(self):
        self.assertEqual(trim_alleles(100, "A", "T"), (100, "A", "T"))

    def test_deletion_with_one_letter_anchor(self):
        self.assertEqual(trim_alleles(100, "CA", "C"), (101, "A", ""))

    def test_deletion_with_longer_anchor_trims_to_the_same_thing(self):
        self.assertEqual(trim_alleles(100, "CAA", "CA"), (101, "A", ""))

    def test_input_is_upper_cased(self):
        self.assertEqual(trim_alleles(5, "ca", "c"), (6, "A", ""))

    def test_missing_input_becomes_empty(self):
        self.assertEqual(trim_alleles(5, None, None), (5, "", ""))


class SameAlleleTest(unittest.TestCase):
    def test_different_anchor_is_the_same_variant(self):
        self.assertTrue(same_allele("CA", "C", "CAA", "CA"))

    def test_different_deleted_base_is_not(self):
        self.assertFalse(same_allele("CA", "C", "CT", "C"))

    def test_snvs_compare_exactly(self):
        self.assertTrue(same_allele("A", "T", "A", "T"))
        self.assertFalse(same_allele("A", "T", "A", "G"))

    def test_indel_and_snv_are_different(self):
        self.assertFalse(same_allele("AA", "C", "A", "C"))

    def test_missing_input_is_never_a_match(self):
        self.assertFalse(same_allele(None, "C", "A", "C"))
        self.assertFalse(same_allele("A", None, "A", "C"))


if __name__ == "__main__":
    unittest.main()
