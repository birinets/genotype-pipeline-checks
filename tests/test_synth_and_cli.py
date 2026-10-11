# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from genotype_checks.runner import failed_ids, load_cases
from genotype_checks.synth import compare_carriers, read_carriers, write_all

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

    def tearDown(self):
        self.tmp.cleanup()

    def test_four_files_are_written(self):
        self.assertEqual(sorted(p.name for p in self.out.iterdir()),
                         ["array.txt", "expected.tsv", "imputed.vcf", "records.tsv"])

    def test_vcf_is_well_formed(self):
        lines = (self.out / "imputed.vcf").read_text().splitlines()
        self.assertEqual(lines[0], "##fileformat=VCFv4.2")
        header = next(l for l in lines if l.startswith("#CHROM"))
        self.assertEqual(header.split("\t")[-1], "SYNTH")
        rows = [l.split("\t") for l in lines if not l.startswith("#")]
        self.assertEqual(len(rows), sum(c["genotype"]["source"] == "imputed" for c in load_cases()))
        for row in rows:
            self.assertEqual(len(row), 10)
            self.assertEqual(row[9], "0/1")  # every imputed case is a heterozygote

    def test_array_file_has_one_row_per_array_case(self):
        rows = [l for l in (self.out / "array.txt").read_text().splitlines() if not l.startswith("#")]
        self.assertEqual(len(rows), sum(c["genotype"]["source"] == "array" for c in load_cases()))
        self.assertTrue(any(r.endswith("\tGG") for r in rows))  # the homozygous case

    def test_records_and_expected_cover_every_case(self):
        ids = {c["id"] for c in load_cases()}
        records = {l.split("\t")[0] for l in (self.out / "records.tsv").read_text().splitlines()[1:]}
        expected = {l.split("\t")[0] for l in (self.out / "expected.tsv").read_text().splitlines()[1:]}
        self.assertEqual(records, ids)
        self.assertEqual(expected, ids)

    def test_a_perfect_pipeline_scores_every_case(self):
        perfect = {c["id"] for c in load_cases() if c["expect"] == "carrier"}
        self.assertEqual(failed_ids(compare_carriers(perfect)), set())

    def test_carrier_list_can_use_coordinates(self):
        called = set()
        for c in load_cases():
            if c["expect"] == "carrier":
                r = c["record"]
                called.add((str(r["chrom"]), int(r["pos"]), r["ref"], r["alt"]))
        self.assertEqual(failed_ids(compare_carriers(called)), set())

    def test_a_pipeline_that_calls_everything_fails_the_traps(self):
        everything = {c["id"] for c in load_cases()}
        traps = {c["id"] for c in load_cases() if c["expect"] == "not_carrier"}
        self.assertEqual(failed_ids(compare_carriers(everything)), traps)

    def test_read_carriers_accepts_ids_coordinates_and_records_files(self):
        path = self.out / "carriers.tsv"
        path.write_text("# comment\ncase_id\tchrom\tpos\tref\talt\n"
                        "array-snv-exact\t1\t6000\tA\tT\n"
                        "chr1\t7000\tA\tG\n"
                        "imputed-snv-exact\n", encoding="utf-8")
        called = read_carriers(path)
        self.assertEqual(called, {("1", 6000, "A", "T"), ("1", 7000, "A", "G"), "imputed-snv-exact"})


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

    def _failures(self, matcher):
        import pipeline_demo
        called = {row[0] for row in pipeline_demo.call_carriers(self.out, matcher)}
        return failed_ids(compare_carriers(called))

    def test_reference_rule_gets_every_case_right_through_the_files(self):
        from genotype_checks.matchers import reference_is_carrier
        self.assertEqual(self._failures(reference_is_carrier), set())

    def test_known_bad_rules_fail_through_the_files_exactly_as_in_memory(self):
        from genotype_checks.cli import EXPECTED_NAIVE_FAILURES
        from genotype_checks.matchers import naive_position_only, naive_string_equal
        self.assertEqual(self._failures(naive_position_only), EXPECTED_NAIVE_FAILURES["position-only"])
        self.assertEqual(self._failures(naive_string_equal), EXPECTED_NAIVE_FAILURES["string-equal"])

    def test_demo_script_and_compare_command_agree(self):
        demo = subprocess.run([sys.executable, str(ROOT / "examples" / "pipeline_demo.py"),
                               str(self.out), "reference"], capture_output=True, text=True, env=ENV)
        self.assertEqual(demo.returncode, 0, demo.stderr)
        carriers = self.out / "carriers.tsv"
        carriers.write_text(demo.stdout, encoding="utf-8")
        compare = run_cli("compare", "--carriers", str(carriers))
        self.assertEqual(compare.returncode, 0, compare.stdout)
        self.assertIn("9/9 cases passed.", compare.stdout)


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
        result = run_cli("check", "--adapter", "genotype_checks.matchers:reference_is_carrier")
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
