# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
"""Check a records VCF against a data VCF and give one outcome per record.

The records VCF lists the variants a pipeline is about to score: CHROM, POS, ID, REF and ALT are
read and any sample columns are ignored. Each alternate allele of a row is checked as its own record.
The data VCF holds one person's genotypes in its GT field. Either file may be plain text or
gzip-compressed, including BGZF, which is a series of gzip members.

A record is compared with the data rows that start at the same position on the same chromosome (a
leading ``chr`` is ignored on both sides), using the reference rule in ``matchers``. Multi-allelic data
rows are read as written: each GT index picks its allele. If several data rows start at one position,
each is compared. When two or more of them describe the record's variant (its two alleles, as written
or complemented) and their outcomes differ, the record is ``no_match``: the kit does not choose between
rows about one variant, and a row without a usable call counts as disagreeing. Rows that describe
other variants at the position, such as split multi-allelic rows or an indel beside an SNV, do not
count against it: the most confirming outcome is kept, and the reason names the others.

Besides the reference rule's five outcomes there is ``no_genotype``: the data has no usable two-allele
genotype at the record's position (no row, a missing call, a half call or a one-allele call). A missing
row is not a reference call.

Not done here: a genome-build check, left-alignment, data rows that start before the record's position
(such as a deletion spanning it), and index-based reading: the whole data file is read.
"""

from __future__ import annotations

import gzip
import re
import zlib
from dataclasses import dataclass
from typing import IO, Iterator, Optional

from .matchers import NO_MATCH, OUTCOMES, describes_variant, reference_outcome

NO_GENOTYPE = "no_genotype"

#: Every outcome this mode writes, from the most confirmed to the least.
VCF_OUTCOMES = OUTCOMES + (NO_GENOTYPE,)

TSV_COLUMNS = ("chrom", "pos", "id", "ref", "alt", "outcome", "genotype", "reason")

_LETTERS = re.compile(r"[ACGTN]+")
_GT_SEPARATOR = re.compile(r"[/|]")


class VcfError(Exception):
    """A file that cannot be read as the VCF this mode needs."""


@dataclass(frozen=True)
class CheckedRecord:
    chrom: str
    pos: int
    id: str
    ref: str
    alt: str
    outcome: str
    genotype: str
    reason: str


@dataclass(frozen=True)
class _DataRow:
    ref: str
    alts: tuple
    letters: Optional[tuple]  # None when the row has no usable two-allele call
    problem: str = ""         # why, when letters is None


def _open_text(path) -> IO[str]:
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def _lines(path) -> Iterator[tuple[int, str]]:
    try:
        with _open_text(path) as fh:
            for number, line in enumerate(fh, 1):
                yield number, line.rstrip("\r\n")
    except UnicodeDecodeError:
        raise VcfError(f"{path}: not a text VCF (a BCF file must first be converted, "
                       "for example with bcftools view)") from None
    except (OSError, EOFError, zlib.error) as exc:
        raise VcfError(f"{path}: cannot be read ({exc})") from None


def _position(text: str, path, number: int) -> int:
    try:
        return int(text)
    except ValueError:
        raise VcfError(f"{path}, line {number}: POS {text!r} is not a whole number") from None


def _chrom_key(chrom: str) -> str:
    return chrom[3:] if chrom[:3].lower() == "chr" else chrom


def read_records(path) -> list[tuple[str, int, str, str, str]]:
    """Return ``(chrom, pos, id, ref, alt)`` for every alternate allele of every row, in file order."""
    records = []
    for number, line in _lines(path):
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 5:
            raise VcfError(f"{path}, line {number}: a VCF row needs at least 5 tab-separated columns")
        chrom, pos, vid, ref, alts = cols[:5]
        position = _position(pos, path, number)
        for alt in alts.split(","):
            records.append((chrom, position, vid, ref, alt))
    return records


def _pick_sample(header: list, sample: Optional[str], path) -> tuple[int, str]:
    samples = header[9:]
    if not samples:
        raise VcfError(f"{path}: the data VCF has no sample column")
    if sample is None:
        if len(samples) > 1:
            raise VcfError(f"{path}: the data VCF has {len(samples)} samples; choose one with --sample")
        return 9, samples[0]
    if sample not in samples:
        raise VcfError(f"{path}: the data VCF has no sample named {sample!r}")
    return 9 + samples.index(sample), sample


def _count(n: int, noun: str) -> str:
    return f"one {noun}" if n == 1 else f"{n} {noun}s"


def _data_row(cols: list, column: int) -> _DataRow:
    ref = cols[3].upper()
    alts = () if cols[4] == "." else tuple(a.upper() for a in cols[4].split(","))
    keys = cols[8].split(":")
    if "GT" not in keys:
        return _DataRow(ref, alts, None, "The data row has no GT field.")
    values = cols[column].split(":")
    index = keys.index("GT")
    gt = values[index] if index < len(values) else "."
    calls = _GT_SEPARATOR.split(gt) if gt else ["."]
    if all(c == "." for c in calls):
        return _DataRow(ref, alts, None, f"The genotype is missing ({gt or 'empty'}).")
    if "." in calls:
        return _DataRow(ref, alts, None, f"The genotype {gt} has a missing allele; the kit compares two "
                                         "called alleles only.")
    if len(calls) != 2:
        return _DataRow(ref, alts, None, f"The genotype {gt} has {_count(len(calls), 'allele')}; the kit "
                                         "compares two-allele genotypes only.")
    alleles = (ref,) + alts
    letters = []
    for call in calls:
        if not call.isdigit() or int(call) >= len(alleles):
            return _DataRow(ref, alts, None, f"The genotype {gt} refers to an allele the row does not list.")
        letters.append(alleles[int(call)])
    return _DataRow(ref, alts, tuple(letters))


def read_data(path, wanted: set, sample: Optional[str] = None) -> tuple[str, dict]:
    """Return the sample's name and ``{(chrom key, pos): [rows]}`` for the positions in ``wanted``."""
    column, name = None, ""
    rows: dict = {}
    for number, line in _lines(path):
        if not line or line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            column, name = _pick_sample(line.split("\t"), sample, path)
            continue
        if line.startswith("#"):
            continue
        if column is None:
            raise VcfError(f"{path}, line {number}: a data row comes before the #CHROM header line")
        first = line.find("\t")
        second = line.find("\t", first + 1)
        if first < 0 or second < 0:
            raise VcfError(f"{path}, line {number}: not a tab-separated VCF row")
        key = (_chrom_key(line[:first]), _position(line[first + 1:second], path, number))
        if key not in wanted:
            continue
        cols = line.split("\t", column + 1)
        if len(cols) <= column:
            raise VcfError(f"{path}, line {number}: the row has no column for sample {name}")
        rows.setdefault(key, []).append(_data_row(cols, column))
    if column is None:
        raise VcfError(f"{path}: no #CHROM header line, so the sample columns are unknown")
    return name, rows


def _record_problem(ref: str, alt: str) -> Optional[str]:
    if alt == ".":
        return "The record lists no alternate allele (ALT is .), so there is nothing to confirm."
    if alt.startswith("<") or alt == "*" or "[" in alt or "]" in alt:
        return f"The record's alternate allele {alt} is symbolic or a breakend; the kit compares letters only."
    if not _LETTERS.fullmatch(ref.upper()) or not _LETTERS.fullmatch(alt.upper()):
        return f"The record {ref}>{alt} uses characters other than A, C, G, T and N."
    return None


def _compare(record: dict, row: _DataRow) -> tuple[str, str, str]:
    if row.letters is None:
        return NO_GENOTYPE, ".", row.problem
    genotype = {"source": "vcf", "alleles": list(row.letters), "row_ref": row.ref, "row_alt": list(row.alts)}
    found = reference_outcome(record, genotype)
    return found.outcome, "/".join(row.letters), found.reason


def _labelled(result: tuple[str, str, str]) -> str:
    outcome, letters, _ = result
    return outcome if letters == "." else f"{outcome} ({letters})"


def _check_one(ref: str, alt: str, rows: list) -> tuple[str, str, str]:
    problem = _record_problem(ref, alt)
    if problem:
        return NO_MATCH, ".", problem
    if not rows:
        return NO_GENOTYPE, ".", ("The data has no row at this position. A missing row is not a reference "
                                  "call: it says nothing about this record.")
    record = {"ref": ref, "alt": alt}
    results = [_compare(record, row) for row in rows]
    if len(results) == 1:
        return results[0]
    # Rows about this variant must agree; a row without a usable call disagrees with one that has
    # it. Rows about other variants at the position (split multi-allelic rows, an SNV and an
    # indel) say nothing against it, so the most confirming row counts.
    about = [result for row, result in zip(rows, results) if describes_variant(record, row.ref, row.alts)]
    if len({outcome for outcome, _, _ in about}) > 1:
        return NO_MATCH, ".", (f"{len(about)} rows describe this variant and disagree: "
                               f"{', '.join(_labelled(r) for r in about)}. The kit does not choose between "
                               "rows about one variant.")
    i = min(range(len(results)), key=lambda k: VCF_OUTCOMES.index(results[k][0]))
    (outcome, letters, reason), others = results[i], results[:i] + results[i + 1:]
    gives = "the other gives" if len(others) == 1 else "the others give"
    return outcome, letters, (f"{reason} The data has {len(rows)} rows at this position; this one gives the "
                              f"most confirming outcome, and {gives} {', '.join(_labelled(r) for r in others)}.")


def check_vcf(records_path, data_path, sample: Optional[str] = None) -> tuple[str, list[CheckedRecord]]:
    """Check every record against the data. Return the sample's name and one result per record."""
    records = read_records(records_path)
    wanted = {(_chrom_key(chrom), pos) for chrom, pos, _, _, _ in records}
    name, rows_at = read_data(data_path, wanted, sample)
    checked = []
    for chrom, pos, vid, ref, alt in records:
        outcome, letters, reason = _check_one(ref, alt, rows_at.get((_chrom_key(chrom), pos), []))
        checked.append(CheckedRecord(chrom, pos, vid, ref, alt, outcome, letters, reason))
    return name, checked


def _field(value) -> str:
    return str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def write_tsv(checked: list[CheckedRecord], fh: IO[str]) -> None:
    """Write the header and one tab-separated row per checked record."""
    fh.write("\t".join(TSV_COLUMNS) + "\n")
    for c in checked:
        values = (c.chrom, c.pos, c.id, c.ref, c.alt, c.outcome, c.genotype, c.reason)
        fh.write("\t".join(_field(v) for v in values) + "\n")


def summary(sample: str, checked: list[CheckedRecord]) -> str:
    """One line counting the outcomes."""
    counts = {outcome: 0 for outcome in VCF_OUTCOMES}
    for c in checked:
        counts[c.outcome] += 1
    parts = [f"{counts[o]} {o}" for o in VCF_OUTCOMES if counts[o]]
    return f"Checked {_count(len(checked), 'record')} against sample {sample}: {', '.join(parts) or 'none'}."
