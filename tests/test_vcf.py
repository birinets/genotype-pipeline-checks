# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
import gzip
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from genotype_checks.matchers import MATCH, NO_MATCH, REF_ALT_SWAPPED
from genotype_checks.runner import load_cases
from genotype_checks.synth import write_all
from genotype_checks.vcf import NO_GENOTYPE, TSV_COLUMNS, check_vcf

ROOT = Path(__file__).resolve().parent.parent
ENV = {**os.environ, "PYTHONPATH": str(ROOT / "src")}


def run_cli(*args):
    return subprocess.run([sys.executable, "-m", "genotype_checks", *args],
                          capture_output=True, text=True, env=ENV, cwd=ROOT)


def vcf_text(header_columns, rows):
    lines = ["##fileformat=VCFv4.2", "#" + "\t".join(header_columns)]
    return "\n".join(lines + ["\t".join(r) for r in rows]) + "\n"


SITE_COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]


class SyntheticFilesTest(unittest.TestCase):
    def test_the_synthetic_files_give_every_case_its_outcome(self):
        """Records for every case against the imputed rows, through the command line.

        Cases with a data row get their own outcome; array cases have no row in the VCF, so
        no_genotype. The data is read plain and as two concatenated gzip members, which is how
        BGZF files are built.
        """
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            write_all(out)
            raw = (out / "imputed.vcf").read_bytes()
            middle = len(raw) // 2
            (out / "imputed.vcf.gz").write_bytes(gzip.compress(raw[:middle]) + gzip.compress(raw[middle:]))
            expected = {c["id"]: c["expect"] if c["genotype"]["source"] == "imputed" else NO_GENOTYPE
                        for c in load_cases()}
            for data in ("imputed.vcf", "imputed.vcf.gz"):
                with self.subTest(data=data):
                    result = run_cli("vcf", "--records", str(out / "records.vcf"), "--data", str(out / data))
                    self.assertEqual(result.returncode, 0, result.stderr)
                    header, *rows = [line.split("\t") for line in result.stdout.splitlines()]
                    self.assertEqual(header, list(TSV_COLUMNS))
                    self.assertTrue(all(len(row) == len(TSV_COLUMNS) for row in rows))
                    self.assertEqual({row[2]: row[5] for row in rows}, expected)


class DataRowsTest(unittest.TestCase):
    """Row shapes the case file cannot express, under each way of naming one chromosome."""

    # (records, data): a leading "chr" on either side, and M and MT for the mitochondrial genome.
    CHROMOSOME_NAMES = [("chr1", "1"), ("1", "chr1"), ("chrM", "MT"), ("MT", "chrM")]

    RECORDS = [
        ("100", "multi-allelic-holds-alt", "A", "T"),
        ("200", "multi-allelic-gt-index", "A", "G"),
        ("300", "no-data-row", "A", "G"),
        ("400", "missing-call", "A", "G"),
        ("500", "half-call", "A", "G"),
        ("600", "one-allele-call", "A", "G"),
        ("700", "two-alt-record", "A", "G,C"),
        ("800", "symbolic-record", "A", "<DEL>"),
        ("900", "split-rows", "A", "T"),
        ("1000", "flipped-and-swapped", "A", "G"),
        ("1100", "row-without-alt", "A", "G"),
        ("1200", "split-rows-other-order", "A", "T"),
        ("1300", "rows-disagree", "A", "G"),
        ("1400", "record-with-n", "A", "N"),  # no data row: the record itself confirms nothing
        ("1500", "gt-with-trailing-space", "A", "G"),
    ]
    DATA = [
        ("100", "A", "G,T", "GT", "1/2"),
        ("200", "A", "T,G", "GT", "0/1"),
        ("400", "A", "G", "GT", "./."),
        ("500", "A", "G", "GT", "./1"),
        ("600", "A", "G", "GT", "1"),
        ("700", "A", "G", "GT:DP", "0/1:12"),
        ("800", "A", "G", "GT", "0/1"),
        ("900", "A", "T", "GT", "0/1"),  # a split 1/2 call: the record's row first, then the other
        ("900", "A", "G", "GT", "1/0"),
        ("1000", "C", "T", "GT", "0/1"),
        ("1100", "A", ".", "GT", "0/0"),
        ("1200", "A", "G", "GT", "1/0"),  # the same pair the other way round, so keeping one row fails
        ("1200", "A", "T", "GT", "0/1"),
        ("1300", "A", "G", "GT", "0/1"),  # two rows for the record's variant that contradict each other
        ("1300", "A", "G", "GT", "0/0"),
        ("1500", "A", "G", "GT", "0/1 "),
    ]
    EXPECTED = {  # (id, alt): (outcome, genotype)
        ("multi-allelic-holds-alt", "T"): (MATCH, "G/T"),
        ("multi-allelic-gt-index", "G"): (NO_MATCH, "A/T"),  # GT 1 is T here, not the record's G
        ("no-data-row", "G"): (NO_GENOTYPE, "."),
        ("missing-call", "G"): (NO_GENOTYPE, "."),
        ("half-call", "G"): (NO_GENOTYPE, "."),
        ("one-allele-call", "G"): (NO_GENOTYPE, "."),
        ("two-alt-record", "G"): (MATCH, "A/G"),
        ("two-alt-record", "C"): (NO_MATCH, "A/G"),
        ("symbolic-record", "<DEL>"): (NO_MATCH, "."),
        ("split-rows", "T"): (MATCH, "A/T"),
        ("flipped-and-swapped", "G"): (REF_ALT_SWAPPED, "C/T"),
        ("row-without-alt", "G"): (NO_MATCH, "A/A"),
        ("split-rows-other-order", "T"): (MATCH, "A/T"),
        ("rows-disagree", "G"): (NO_MATCH, "."),
        ("record-with-n", "N"): (NO_MATCH, "."),
        ("gt-with-trailing-space", "G"): (MATCH, "A/G"),
    }
    # Where several rows meet one record, the reason must name what the other rows said.
    REASON_MENTIONS = {
        ("split-rows", "T"): ["no_match (G/A)"],
        ("rows-disagree", "G"): ["describe this variant and disagree", "match (A/G)", "no_match (A/A)"],
    }

    def test_each_row_shape_gets_its_outcome_and_genotype(self):
        for records_chrom, data_chrom in self.CHROMOSOME_NAMES:
            with self.subTest(records=records_chrom, data=data_chrom), tempfile.TemporaryDirectory() as tmp:
                records, data = Path(tmp) / "records.vcf", Path(tmp) / "data.vcf"
                records.write_text(vcf_text(SITE_COLUMNS, [[records_chrom, pos, vid, ref, alt, ".", ".", "."]
                                                           for pos, vid, ref, alt in self.RECORDS]))
                data.write_text(vcf_text(SITE_COLUMNS + ["FORMAT", "S1"],
                                         [[data_chrom, pos, ".", ref, alt, ".", "PASS", ".", fmt, sample]
                                          for pos, ref, alt, fmt, sample in self.DATA]))
                name, checked = check_vcf(records, data)
                self.assertEqual(name, "S1")
                self.assertEqual({(c.id, c.alt): (c.outcome, c.genotype) for c in checked}, self.EXPECTED)
                self.assertEqual(len(checked), len(self.EXPECTED))
                reasons = {(c.id, c.alt): c.reason for c in checked}
                for key, mentions in self.REASON_MENTIONS.items():
                    for text in mentions:
                        self.assertIn(text, reasons[key], key)


class UnusableInputTest(unittest.TestCase):
    """Input the check cannot use stops it with exit code 2 and says why, never with a traceback."""

    RECORDS = vcf_text(SITE_COLUMNS, [["1", "100", "r1", "A", "G", ".", ".", "."]])
    DATA = vcf_text(SITE_COLUMNS + ["FORMAT", "S1"], [["1", "100", ".", "A", "G", ".", "PASS", ".", "GT", "0/1"]])
    BCF = gzip.compress(b"BCF\x02\x02" + bytes([0xC8, 0, 0, 0]) + b"##fileformat=VCFv4.2\n\x00")
    # name: (records file, data file, extra arguments, text the message must contain).
    # None means the file is not there; {out} stands for a path inside a directory that does not exist.
    CASES = {
        "missing file": (None, DATA, [], "cannot be read"),
        "truncated gzip": (RECORDS, gzip.compress(DATA.encode())[:-12], [], "cannot be read"),
        "not UTF-8": (RECORDS, DATA.replace("VCFv4.2", "VCFv4.2 \xe9").encode("latin-1"), [], "not UTF-8"),
        "BCF": (RECORDS, BCF, [], "must first be converted to VCF"),
        "short records row": ("1\t100\tr1\tA\n", DATA, [], "at least 5 tab-separated columns"),
        "position not a number": ("1\tx\tr1\tA\tG\n", DATA, [], "is not a whole number"),
        "no header line": (RECORDS, "##fileformat=VCFv4.2\n", [], "no #CHROM header line"),
        "row before the header": (RECORDS, "1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\n", [], "comes before the #CHROM"),
        "row without tabs": (RECORDS, DATA + "1 100 . A G\n", [], "not a tab-separated VCF row"),
        "no sample column": (RECORDS, vcf_text(SITE_COLUMNS, []), [], "no sample column"),
        "row without the sample": (RECORDS, DATA + "1\t100\t.\tA\tG\t.\t.\t.\tGT\n", [], "no column for sample"),
        "several samples": (RECORDS, vcf_text(SITE_COLUMNS + ["FORMAT", "S1", "S2"], []), [], "--sample"),
        "unknown sample": (RECORDS, DATA, ["--sample", "S9"], "no sample named 'S9'"),
        "unwritable output": (RECORDS, DATA, ["--out", "{out}"], "cannot write"),
    }

    def test_each_unusable_input_exits_with_two_and_a_reason(self):
        for name, (records_text, data_text, extra, message) in self.CASES.items():
            with self.subTest(name), tempfile.TemporaryDirectory() as tmp:
                paths = []
                for stem, content in (("records", records_text), ("data", data_text)):
                    path = Path(tmp) / f"{stem}.vcf"
                    if isinstance(content, bytes):
                        path.write_bytes(content)
                    elif content is not None:
                        path.write_text(content)
                    paths.append(str(path))
                out = str(Path(tmp) / "no-such-directory" / "outcomes.tsv")
                result = run_cli("vcf", "--records", paths[0], "--data", paths[1],
                                 *[arg.replace("{out}", out) for arg in extra])
                self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
                self.assertIn(message, result.stderr)


class SampleChoiceTest(unittest.TestCase):
    def test_a_named_sample_is_the_one_read(self):
        with tempfile.TemporaryDirectory() as tmp:
            records, data = Path(tmp) / "records.vcf", Path(tmp) / "data.vcf"
            records.write_text(vcf_text(SITE_COLUMNS, [["1", "100", "r1", "A", "G", ".", ".", "."]]))
            data.write_text(vcf_text(SITE_COLUMNS + ["FORMAT", "S1", "S2"],
                                     [["1", "100", ".", "A", "G", ".", "PASS", ".", "GT", "0/0", "1/1"]]))
            for sample, row in (("S1", [NO_MATCH, "A/A"]), ("S2", [MATCH, "G/G"])):
                with self.subTest(sample=sample):
                    out = Path(tmp) / f"{sample}.tsv"
                    result = run_cli("vcf", "--records", str(records), "--data", str(data),
                                     "--sample", sample, "--out", str(out))
                    self.assertEqual(result.returncode, 0, result.stderr)
                    self.assertEqual(out.read_text().splitlines()[1].split("\t")[5:7], row)


if __name__ == "__main__":
    unittest.main()
