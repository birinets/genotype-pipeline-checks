# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
"""Write the cases out as made-up subject files, so a whole pipeline can be tried end to end.

Output, in one directory:

    array.txt     array-style raw genotype file (tab-separated: rsid, chromosome, position, genotype)
    imputed.vcf   one-sample VCF for the cases whose data rows carry their own ref and alt
    records.tsv   the variant records (case id, chrom, pos, ref, alt), as a variant database would list them
    expected.tsv  the answer for each case id: carrier or not_carrier

Everything is synthetic. Nothing here describes a real person or a real variant.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .runner import Result, load_cases

SAMPLE = "SYNTH"


def _rsid(index: int) -> str:
    return f"rsSYN{index + 1:03d}"


def _vcf_gt(alleles: list[str], ref: str, alt: str) -> str:
    code = {ref: "0", alt: "1"}
    try:
        return "/".join(sorted(code[a] for a in alleles))
    except KeyError as exc:
        raise ValueError(f"allele {exc} is neither row ref {ref!r} nor row alt {alt!r}") from None


def array_text(cases: list[dict[str, Any]]) -> str:
    lines = [
        "# genotype-pipeline-checks: synthetic array file. Not a real person.",
        "# rsid\tchromosome\tposition\tgenotype",
    ]
    for i, case in enumerate(cases):
        g = case["genotype"]
        if g["source"] != "array":
            continue
        if any(len(a) != 1 for a in g["alleles"]):
            raise ValueError(f"{case['id']}: array alleles must be single letters")
        r = case["record"]
        lines.append(f"{_rsid(i)}\t{r['chrom']}\t{r['pos']}\t{''.join(g['alleles'])}")
    return "\n".join(lines) + "\n"


def imputed_vcf(cases: list[dict[str, Any]]) -> str:
    rows = []
    for i, case in enumerate(cases):
        g = case["genotype"]
        if g["source"] != "imputed":
            continue
        r = case["record"]
        gt = _vcf_gt(g["alleles"], g["row_ref"], g["row_alt"])
        rows.append((r["chrom"], int(r["pos"]),
                     f"{r['chrom']}\t{r['pos']}\t{_rsid(i)}\t{g['row_ref']}\t{g['row_alt']}\t.\tPASS\t.\tGT\t{gt}"))
    rows.sort(key=lambda t: (t[0], t[1]))
    contigs = sorted({t[0] for t in rows})
    head = [
        "##fileformat=VCFv4.2",
        "##source=genotype-pipeline-checks (synthetic, not a real person)",
        *[f"##contig=<ID={c}>" for c in contigs],
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}",
    ]
    return "\n".join(head + [t[2] for t in rows]) + "\n"


def records_tsv(cases: list[dict[str, Any]]) -> str:
    lines = ["case_id\tchrom\tpos\tref\talt"]
    for case in cases:
        r = case["record"]
        lines.append(f"{case['id']}\t{r['chrom']}\t{r['pos']}\t{r['ref']}\t{r['alt']}")
    return "\n".join(lines) + "\n"


def expected_tsv(cases: list[dict[str, Any]]) -> str:
    lines = ["case_id\texpect"] + [f"{c['id']}\t{c['expect']}" for c in cases]
    return "\n".join(lines) + "\n"


def write_all(outdir: str | Path, cases: list[dict[str, Any]] | None = None) -> list[Path]:
    """Write the four files into ``outdir`` (created if needed) and return their paths."""
    cases = cases if cases is not None else load_cases()
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    files = {
        "array.txt": array_text(cases),
        "imputed.vcf": imputed_vcf(cases),
        "records.tsv": records_tsv(cases),
        "expected.tsv": expected_tsv(cases),
    }
    paths = []
    for name, text in files.items():
        path = out / name
        path.write_text(text, encoding="utf-8")
        paths.append(path)
    return paths


def compare_carriers(called: set, cases: list[dict[str, Any]] | None = None):
    """Score a list of carrier calls from a pipeline run on the synthetic files.

    ``called`` holds case ids and/or ``(chrom, pos, ref, alt)`` tuples for every record the
    pipeline reported the subject as carrying. Returns a list of ``runner.Result``.
    """
    results = []
    for case in cases if cases is not None else load_cases():
        r = case["record"]
        key = (str(r["chrom"]), int(r["pos"]), r["ref"].upper(), r["alt"].upper())
        got = case["id"] in called or key in called
        results.append(Result(case["id"], case["expect"] == "carrier", got, case["why"]))
    return results


def read_carriers(path: str | Path) -> set:
    """Read a carrier list: one case id per line, or tab-separated chrom, pos, ref, alt."""
    called: set = set()
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if parts[0] in ("case_id", "chrom"):  # header row
            continue
        if len(parts) >= 5 and parts[2].isdigit():  # records.tsv style: case id first
            called.add((parts[1].removeprefix("chr"), int(parts[2]), parts[3].upper(), parts[4].upper()))
        elif len(parts) >= 4 and parts[1].isdigit():  # chrom, pos, ref, alt
            called.add((parts[0].removeprefix("chr"), int(parts[1]), parts[2].upper(), parts[3].upper()))
        else:
            called.add(parts[0])
    return called
