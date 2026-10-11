# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
"""Command line: list, selftest, check, synth, compare."""

from __future__ import annotations

import argparse
import importlib
import sys

from . import __version__
from .matchers import naive_position_only, naive_string_equal, reference_is_carrier
from .runner import failed_ids, format_report, load_cases, run_cases
from .synth import compare_carriers, read_carriers, write_all

# What the two known-bad matchers must get wrong. If these sets change, the cases changed.
EXPECTED_NAIVE_FAILURES = {
    "position-only": {
        "indel-vs-array-snp",
        "array-genotype-outside-record-alleles",
        "imputed-indel-different-deletion",
        "imputed-indel-different-anchor",
    },
    "string-equal": {"imputed-indel-different-anchor"},
}


def _cmd_list(_args) -> int:
    for case in load_cases():
        print(f"{case['id']}\n    expect: {case['expect']}\n    why:    {case['why']}")
    return 0


def _cmd_selftest(_args) -> int:
    ok = True
    ref_results = run_cases(reference_is_carrier)
    print("Reference rule:")
    print(format_report(ref_results))
    ok &= not failed_ids(ref_results)

    for name, fn in (("position-only", naive_position_only), ("string-equal", naive_string_equal)):
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
    results = run_cases(_load_adapter(args.adapter))
    print(format_report(results))
    return 1 if failed_ids(results) else 0


def _cmd_synth(args) -> int:
    for path in write_all(args.outdir):
        print(path)
    return 0


def _cmd_compare(args) -> int:
    results = compare_carriers(read_carriers(args.carriers))
    print(format_report(results))
    return 1 if failed_ids(results) else 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="genotype-checks",
        description="Allele-identity checks for genotype pipelines. Synthetic cases only.")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("list", help="show every case and why it is there").set_defaults(fn=_cmd_list)
    sub.add_parser("selftest", help="check the kit itself").set_defaults(fn=_cmd_selftest)

    p = sub.add_parser("check", help="run your matcher function over the cases")
    p.add_argument("--adapter", required=True, help="package.module:function taking (record, genotype)")
    p.set_defaults(fn=_cmd_check)

    p = sub.add_parser("synth", help="write the cases as made-up subject files")
    p.add_argument("outdir")
    p.set_defaults(fn=_cmd_synth)

    p = sub.add_parser("compare", help="score the carriers your pipeline reported on the synthetic files")
    p.add_argument("--carriers", required=True,
                   help="file with one case id per line, or tab-separated chrom, pos, ref, alt")
    p.set_defaults(fn=_cmd_compare)

    args = parser.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
