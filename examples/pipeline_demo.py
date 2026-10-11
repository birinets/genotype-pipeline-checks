# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
"""A toy pipeline that reads the synthetic files and reports matches, to show the whole loop.

    python -m genotype_checks synth out/
    python examples/pipeline_demo.py out/ reference > reported.tsv
    python -m genotype_checks compare --reported reported.tsv

The second argument picks the matcher: reference, position-only, any-non-reference or
string-equal. A real pipeline replaces this file: read ``array.txt`` and ``imputed.vcf``, look
up each variant in ``records.tsv``, and print the records it would report for the subject, in
the same columns.
"""

from __future__ import annotations

import sys
from pathlib import Path

from genotype_checks.matchers import KNOWN_BAD, reference_matches

MATCHERS = {"reference": reference_matches, **KNOWN_BAD}


def read_subject(indir: Path) -> dict:
    """Return {(chrom, pos): genotype} from array.txt and imputed.vcf."""
    subject = {}
    for line in (indir / "array.txt").read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        _rsid, chrom, pos, letters = line.split("\t")
        subject[(chrom, int(pos))] = {"source": "array", "alleles": list(letters)}
    for line in (indir / "imputed.vcf").read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        chrom, pos, _id, ref, alt, _q, _f, _info, _fmt, sample = line.split("\t")
        options = {"0": ref, "1": alt}
        subject[(chrom, int(pos))] = {
            "source": "imputed",
            "alleles": [options[i] for i in sample.split("/")],
            "row_ref": ref,
            "row_alt": alt,
        }
    return subject


def report_matches(indir: Path, matcher) -> list[tuple[str, str, int, str, str]]:
    """Return (case_id, chrom, pos, ref, alt) for every record the matcher reports for the subject."""
    subject = read_subject(indir)
    reported = []
    for line in (indir / "records.tsv").read_text(encoding="utf-8").splitlines()[1:]:
        case_id, chrom, pos, ref, alt = line.split("\t")
        record = {"chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt}
        genotype = subject.get((chrom, int(pos)))
        if genotype and matcher(record, genotype):
            reported.append((case_id, chrom, int(pos), ref, alt))
    return reported


def main(argv=None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if len(args) != 2 or args[1] not in MATCHERS:
        print(f"usage: pipeline_demo.py DIR {{{'|'.join(MATCHERS)}}}", file=sys.stderr)
        return 2
    print("case_id\tchrom\tpos\tref\talt")
    for row in report_matches(Path(args[0]), MATCHERS[args[1]]):
        print("\t".join(str(x) for x in row))
    return 0


if __name__ == "__main__":
    sys.exit(main())
