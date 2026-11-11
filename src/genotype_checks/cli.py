# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
"""Command line: list, selftest, check, synth, compare, vcf."""

from __future__ import annotations

import argparse
import importlib
import sys

from . import __version__
from .matchers import KNOWN_BAD, reference_matches, reference_outcome
from .runner import UnknownOutcome, failed_ids, format_report, load_cases, run_cases
from .synth import compare_reported, read_reported, write_all
from .vcf import VcfError, check_vcf, summary, write_tsv

# What each known-bad matcher must get wrong. If these sets change, the cases changed.
EXPECTED_NAIVE_FAILURES = {
    "position-only": {
        "delins-vs-array-snp",
        "deletion-anchor-equals-alt",
        "array-genotype-outside-record-alleles",
        "imputed-indel-different-deletion",
        "imputed-row-reference-disagrees",
        "imputed-indel-different-anchor",
        "array-palindromic-homozygous",
        "imputed-ref-alt-swapped",
        "imputed-palindromic-reference-disagrees",
    },
    "any-non-reference": {
        "delins-vs-array-snp",
        "deletion-anchor-equals-alt",
        "array-alt-mismatch",
        "imputed-alt-mismatch",
        "array-genotype-outside-record-alleles",
        "imputed-indel-different-deletion",
        "imputed-row-reference-disagrees",
        "imputed-reference-disagrees-in-trimmed-base",
        "array-strand-flip",
        "imputed-strand-flip",
        "array-palindromic-homozygous",
        "imputed-ref-alt-swapped",
        "imputed-palindromic-reference-disagrees",
    },
    "string-equal": {"imputed-indel-different-anchor"},
    "trim-only": {"imputed-reference-disagrees-in-trimmed-base"},
    "as-written": {"array-palindromic-homozygous"},
    "auto-flip": {
        "array-strand-flip",
        "imputed-strand-flip",
        "array-palindromic-homozygous",
        "array-palindromic-homozygous-reference",
    },
    "allele-set": {"imputed-ref-alt-swapped", "imputed-palindromic-reference-disagrees"},
}


def _cmd_list(_args) -> int:
    for case in load_cases():
        print(f"{case['id']}\n    expect: {case['expect']}\n    why:    {case['why']}")
    return 0


def _cmd_selftest(_args) -> int:
    ok = True
    ref_results = run_cases(reference_outcome)
    print("Reference rule, outcomes:")
    print(format_report(ref_results))
    ok &= not failed_ids(ref_results)

    yes_no = failed_ids(run_cases(reference_matches))
    ok &= not yes_no
    print(f"\nReference rule, yes/no form: {'ok' if not yes_no else 'FAILED ' + str(sorted(yes_no))}")

    for name, fn in KNOWN_BAD.items():
        got = failed_ids(run_cases(fn))
        want = EXPECTED_NAIVE_FAILURES[name]
        same = got == want
        ok &= same
        print(f"\nKnown-bad matcher '{name}' fails {len(got)} case(s): "
              f"{'as expected' if same else 'UNEXPECTED'}")
        for case_id in sorted(got):
            print(f"    {case_id}")
        if not same:
            print(f"    expected: {sorted(want)}")
    print("\nselftest:", "ok" if ok else "FAILED")
    return 0 if ok else 1


def _load_adapter(spec: str):
    module_name, _, attr = spec.partition(":")
    if not module_name or not attr:
        raise SystemExit("--adapter must look like package.module:function")
    sys.path.insert(0, ".")
    return getattr(importlib.import_module(module_name), attr)


def _cmd_check(args) -> int:
    matcher = _load_adapter(args.adapter)
    try:
        results = run_cases(matcher)
    except UnknownOutcome as exc:
        print(f"genotype-checks check: {exc}", file=sys.stderr)
        return 2
    print(format_report(results))
    return 1 if failed_ids(results) else 0


def _cmd_synth(args) -> int:
    for path in write_all(args.outdir):
        print(path)
    return 0


def _cmd_compare(args) -> int:
    results = compare_reported(read_reported(args.reported))
    print(format_report(results))
    return 1 if failed_ids(results) else 0


def _cmd_vcf(args) -> int:
    try:
        sample, checked = check_vcf(args.records, args.data, args.sample)
    except VcfError as exc:
        print(f"genotype-checks vcf: {exc}", file=sys.stderr)
        return 2
    if args.out:
        try:
            with open(args.out, "w", encoding="utf-8", newline="") as fh:
                write_tsv(checked, fh)
        except OSError as exc:
            print(f"genotype-checks vcf: cannot write {args.out} ({exc.strerror or exc})", file=sys.stderr)
            return 2
    else:
        write_tsv(checked, sys.stdout)
    print(summary(sample, checked), file=sys.stderr)
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="genotype-checks",
        description="Allele-identity checks for genotype pipelines. Synthetic cases only.")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("list", help="show every case and why it is there").set_defaults(fn=_cmd_list)
    sub.add_parser("selftest", help="check the kit itself").set_defaults(fn=_cmd_selftest)

    p = sub.add_parser("check", help="run your matcher function over the cases")
    p.add_argument("--adapter", required=True,
                   help="package.module:function taking (record, genotype) and returning True/False "
                        "or an outcome name")
    p.set_defaults(fn=_cmd_check)

    p = sub.add_parser("synth", help="write the cases as made-up subject files")
    p.add_argument("outdir")
    p.set_defaults(fn=_cmd_synth)

    p = sub.add_parser("compare", help="score the records your pipeline reported on the synthetic files")
    p.add_argument("--reported", required=True,
                   help="file with one case id per line, or tab-separated chrom, pos, ref, alt")
    p.set_defaults(fn=_cmd_compare)

    p = sub.add_parser("vcf", help="check a records VCF against a data VCF: one outcome per record")
    p.add_argument("--records", required=True,
                   help="VCF of the variants a pipeline may report (plain or gzip/BGZF; samples ignored)")
    p.add_argument("--data", required=True, help="VCF with the person's genotypes in GT (plain or gzip/BGZF)")
    p.add_argument("--sample", help="sample column to read; required when the data VCF has several")
    p.add_argument("--out", help="write the TSV here instead of to standard output")
    p.set_defaults(fn=_cmd_vcf)

    args = parser.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
