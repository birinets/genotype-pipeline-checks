# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
"""One reference rule for "does this data match this recorded variant?", and four known-bad ones.

A matcher takes a *record* and a *genotype* and returns True if the data matches the record, so
that a pipeline may report the record for this subject.

    record    {"chrom": "1", "pos": 100, "ref": "A", "alt": "C"}    as a variant database writes it
    genotype  {"source": "array" | "imputed",
               "alleles": ["A", "C"],                               the two alleles the subject holds
               "row_ref": "A", "row_alt": "C"}                      imputed or sequenced rows only

An array reports only the letters it observed, so ``row_ref`` and ``row_alt`` are absent.

False means "this data does not confirm this record". It is not a negative result for the person.

The reference rule is one reasonable reading, not the only one. Its main choices:

* an array genotype can confirm only a single-letter record;
* array letters are taken as written on the forward strand of the record's genome build, so
  strand is not resolved here;
* a row with its own ``ref`` and ``alt`` must agree with the record about the reference bases
  they share, and is then compared by minimal representation;
* malformed input (identical ref and alt, a row with only one of ``row_ref`` and ``row_alt``)
  never matches.

If your pipeline disagrees on purpose, the case tells you where.
"""

from __future__ import annotations

from typing import Any

from .alleles import same_allele

Record = dict[str, Any]
Genotype = dict[str, Any]


def _clean(record: Record, genotype: Genotype):
    ref = (record.get("ref") or "").upper()
    alt = (record.get("alt") or "").upper()
    alleles = [str(a).upper() for a in (genotype.get("alleles") or []) if a]
    row_ref = (genotype.get("row_ref") or "").upper() or None
    row_alt = (genotype.get("row_alt") or "").upper() or None
    return ref, alt, alleles, row_ref, row_alt


def _usable(ref: str, alt: str, alleles: list) -> bool:
    return bool(ref and alt and ref != alt and len(alleles) == 2)


def _reference_bases_agree(ref: str, row_ref: str) -> bool:
    """Both start at the same position, so the shorter reference must be a prefix of the longer."""
    n = min(len(ref), len(row_ref))
    return ref[:n] == row_ref[:n]


def _array_rule(ref: str, alt: str, observed: set) -> bool:
    if len(ref) == 1 and len(alt) == 1:
        return alt in observed and observed <= {ref, alt}
    return False


def reference_matches(record: Record, genotype: Genotype) -> bool:
    """The reference rule. See the module docstring for its main choices."""
    ref, alt, alleles, row_ref, row_alt = _clean(record, genotype)
    if not _usable(ref, alt, alleles):
        return False
    observed = set(alleles)
    if row_ref and row_alt:
        if not _reference_bases_agree(ref, row_ref):
            return False
        return same_allele(row_ref, row_alt, ref, alt) and row_alt in observed
    if row_ref or row_alt:
        return False
    return _array_rule(ref, alt, observed)


def naive_position_only(record: Record, genotype: Genotype) -> bool:
    """Known-bad: matches whenever the record's alternate allele appears in the genotype.

    This is the position-only join that lets an indel inherit a SNP genotype. Included so the
    cases can be shown to catch it.
    """
    _, alt, alleles, _, _ = _clean(record, genotype)
    return bool(alt) and alt in set(alleles)


def naive_any_non_reference(record: Record, genotype: Genotype) -> bool:
    """Known-bad: matches whenever the genotype shows anything other than the record's reference.

    "There is a variant at a recorded position" is not "this is the recorded variant".
    """
    ref, _, alleles, _, _ = _clean(record, genotype)
    return bool(ref) and any(a != ref for a in alleles)


def naive_string_equal(record: Record, genotype: Genotype) -> bool:
    """Known-bad: like the reference rule, but compares ref and alt letter for letter.

    It rejects every wrong match, so only the case where the same variant is written with a
    different anchor length shows it up.
    """
    ref, alt, alleles, row_ref, row_alt = _clean(record, genotype)
    if not _usable(ref, alt, alleles):
        return False
    observed = set(alleles)
    if row_ref and row_alt:
        return (row_ref, row_alt) == (ref, alt) and row_alt in observed
    if row_ref or row_alt:
        return False
    return _array_rule(ref, alt, observed)


def naive_trim_only(record: Record, genotype: Genotype) -> bool:
    """Known-bad: compares minimal representations without checking the reference bases agree.

    Trimming throws away the bases two encodings share, including a base on which they
    disagree. This was this kit's own reference rule in version 0.1.0.
    """
    ref, alt, alleles, row_ref, row_alt = _clean(record, genotype)
    if not _usable(ref, alt, alleles):
        return False
    observed = set(alleles)
    if row_ref and row_alt:
        return same_allele(row_ref, row_alt, ref, alt) and row_alt in observed
    if row_ref or row_alt:
        return False
    return _array_rule(ref, alt, observed)


KNOWN_BAD = {
    "position-only": naive_position_only,
    "any-non-reference": naive_any_non_reference,
    "string-equal": naive_string_equal,
    "trim-only": naive_trim_only,
}
