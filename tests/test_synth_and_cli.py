# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from genotype_checks.runner import failed_ids, load_cases
from genotype_checks.synth import compare_reported, read_reported, write_all

ROOT = Path(__file__).resolve().parent.parent
ENV = {**os.environ, "PYTHONPATH": str(ROOT / "src")}


def run_cli(*args):
    return subprocess.run([sys.executable, "-m", "genotype_checks", *args],
                          capture_output=True, text=True, env=ENV, cwd=ROOT)


class SynthTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = Path(self.tmp.name)
        write_all(self.out)
        self.by_pos = {int(c["record"]["pos"]): c for c in load_cases()}

    def tearDown(self):
        self.tmp.cleanup()

    def test_four_files_are_written(self):
        self.assertEqual(sorted(p.name for p in self.out.iterdir()),
                         ["array.txt", "expected.tsv", "imputed.vcf", "records.tsv"])

    def test_vcf_rows_carry_each_imputed_case_as_written(self):
        lines = (self.out / "imputed.vcf").read_text().splitlines()
        self.assertEqual(lines[0], "##fileformat=VCFv4.2")
        header = next(l for l in lines if l.startswith("#CHROM"))
        self.assertEqual(header.split("\t")[-1], "SYNTH")
        rows = [l.split("\t") for l in lines if not l.startswith("#")]
        imputed = [c for c in load_cases() if c["genotype"]["source"] == "imputed"]
        self.assertEqual(len(rows), len(imputed))
        for row in rows:
            self.assertEqual(len(row), 10)
            g = self.by_pos[int(row[1])]["genotype"]
            self.assertEqual((row[3], row[4]), (g["row_ref"], g["row_alt"]))
            alt_copies = sum(a == g["row_alt"] for a in g["alleles"])
            self.assertEqual(row[9].count("1"), alt_copies)

    def test_array_rows_carry_each_array_case_as_written(self):
        rows = [l.split("\t") for l in (self.out / "array.txt").read_text().splitlines()
                if not l.startswith("#")]
        array = [c for c in load_cases() if c["genotype"]["source"] == "array"]
        self.assertEqual(len(rows), len(array))
        for row in rows:
            self.assertEqual(list(row[3]), self.by_pos[int(row[2])]["genotype"]["alleles"])

    def test_records_and_expected_cover_every_case(self):
        ids = {c["id"] for c in load_cases()}
        records = {l.split("\t")[0] for l in (self.out / "records.tsv").read_text().splitlines()[1:]}
        expected = {l.split("\t")[0] for l in (self.out / "expected.tsv").read_text().splitlines()[1:]}
        self.assertEqual(records, ids)
        self.assertEqual(expected, ids)

    def test_a_perfect_pipeline_scores_every_case(self):
        perfect = {c["id"] for c in load_cases() if c["expect"] == "match"}
        self.assertEqual(failed_ids(compare_reported(perfect)), set())

    def test_reported_list_can_use_coordinates(self):
        reported = set()
        for c in load_cases():
            if c["expect"] == "match":
                r = c["record"]
                reported.add((str(r["chrom"]), int(r["pos"]), r["ref"], r["alt"]))
        self.assertEqual(failed_ids(compare_reported(reported)), set())

    def test_a_pipeline_that_reports_everything_fails_every_no_match_case(self):
        everything = {c["id"] for c in load_cases()}
        traps = {c["id"] for c in load_cases() if c["expect"] == "no_match"}
        self.assertEqual(failed_ids(compare_reported(everything)), traps)

    def test_read_reported_accepts_ids_coordinates_and_records_files(self):
        path = self.out / "reported.tsv"
        path.write_text("# comment\ncase_id\tchrom\tpos\tref\talt\n"
                        "array-snv-exact\t1\t6000\tA\tC\n"
                        "chr1\t7000\tA\tG\n"
                        "imputed-snv-exact\n", encoding="utf-8")
        self.assertEqual(read_reported(path),
                         {("1", 6000, "A", "C"), ("1", 7000, "A", "G"), "imputed-snv-exact"})


class RoundTripTest(unittest.TestCase):
    """Write the files, read them back as a toy pipeline would, and score the result."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = Path(self.tmp.name)
        write_all(self.out)
        sys.path.insert(0, str(ROOT / "examples"))

    def tearDown(self):
        self.tmp.cleanup()
        sys.path.remove(str(ROOT / "examples"))

    def test_reference_rule_gets_every_case_right_through_the_files(self):
        import pipeline_demo
        from genotype_checks.matchers import reference_matches
        reported = {row[0] for row in pipeline_demo.report_matches(self.out, reference_matches)}
        self.assertEqual(failed_ids(compare_reported(reported)), set())

    def test_demo_script_and_compare_command_agree(self):
        demo = subprocess.run([sys.executable, str(ROOT / "examples" / "pipeline_demo.py"),
                               str(self.out), "reference"], capture_output=True, text=True, env=ENV)
        self.assertEqual(demo.returncode, 0, demo.stderr)
        reported = self.out / "reported.tsv"
        reported.write_text(demo.stdout, encoding="utf-8")
        compare = run_cli("compare", "--reported", str(reported))
        self.assertEqual(compare.returncode, 0, compare.stdout)
        n = len(load_cases())
        self.assertIn(f"{n}/{n} cases passed.", compare.stdout)


class CliTest(unittest.TestCase):
    def test_selftest_passes(self):
        result = run_cli("selftest")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("selftest: ok", result.stdout)

    def test_list_shows_every_case(self):
        result = run_cli("list")
        self.assertEqual(result.returncode, 0)
        for c in load_cases():
            self.assertIn(c["id"], result.stdout)

    def test_check_with_the_reference_rule_passes(self):
        result = run_cli("check", "--adapter", "genotype_checks.matchers:reference_matches")
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_check_with_a_bad_rule_fails_with_exit_code_one(self):
        result = run_cli("check", "--adapter", "genotype_checks.matchers:naive_position_only")
        self.assertEqual(result.returncode, 1)
        self.assertIn("FAIL", result.stdout)

    def test_a_malformed_adapter_is_refused(self):
        result = run_cli("check", "--adapter", "nonsense")
        self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
