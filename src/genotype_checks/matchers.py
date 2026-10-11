# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
"""One reference rule for "does this subject carry this recorded variant?", and two known-bad ones.

A matcher takes a *record* and a *genotype* and returns True if the subject should be counted as
a carrier.

    record    {"chrom": "1", "pos": 100, "ref": "A", "alt": "T"}    as a variant database writes it
    genotype  {"source": "array" | "imputed",
               "alleles": ["A", "T"],                               the two alleles the subject holds
               "row_ref": "A", "row_alt": "T"}                      imputed or sequenced rows only

An array reports only the letters it observed, so ``row_ref`` and ``row_alt`` are absent.

The reference rule is one reasonable reading, not the only one. Its two main choices are that
an array genotype can confirm only a single-letter record, and that rows with their own ``ref``
and ``alt`` are compared by minimal representation. If your pipeline disagrees on purpose, the
case tells you where.
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


def reference_is_carrier(record: Record, genotype: Genotype) -> bool:
    """The reference rule. See the module docstring for its two main choices."""
    ref, alt, alleles, row_ref, row_alt = _clean(record, genotype)
    if not (ref and alt and len(alleles) == 2):
        return False
    observed = set(alleles)
    if row_ref and row_alt:
        return same_allele(row_ref, row_alt, ref, alt) and row_alt in observed
    if len(ref) == 1 and len(alt) == 1:
        return alt in observed and observed <= {ref, alt}
    return False


def naive_position_only(record: Record, genotype: Genotype) -> bool:
    """Known-bad: counts a carrier whenever the record's alternate allele appears in the genotype.

    This is the position-only join that let an indel inherit a SNP genotype. Included so the
    cases can be shown to catch it.
    """
    _, alt, alleles, _, _ = _clean(record, genotype)
    return bool(alt) and alt in set(alleles)


def naive_string_equal(record: Record, genotype: Genotype) -> bool:
    """Known-bad: like the reference rule, but compares ref and alt letter for letter.

    It rejects every wrong match, so only the case where a true carrier is written with a
    different anchor length shows it up.
    """
    ref, alt, alleles, row_ref, row_alt = _clean(record, genotype)
    if not (ref and alt and len(alleles) == 2):
        return False
    observed = set(alleles)
    if row_ref and row_alt:
        return (row_ref, row_alt) == (ref, alt) and row_alt in observed
    if len(ref) == 1 and len(alt) == 1:
        return alt in observed and observed <= {ref, alt}
    return False
