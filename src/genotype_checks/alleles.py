# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
"""Allele comparison by minimal representation.

Two valid encodings of the same variant can differ letter for letter: a different anchor length,
for example, writes one deletion as ``CA>C`` or as ``CAA>CA``. Comparing the raw strings misses
the match. Trimming the shared suffix and then the shared prefix (the "minimal representation"
of vt and Tan et al. 2015) makes the comparison independent of how the variant was written.

This is trimming only. It does not left-align repeats against a reference genome, so two
encodings of a deletion inside a repeat that sit at different positions still compare as
different. A single-letter change (SNV) is returned unchanged, so this never loosens SNV matching.
"""

from __future__ import annotations


def trim_alleles(pos: int, ref: str | None, alt: str | None) -> tuple[int, str, str]:
    """Return ``(pos, ref, alt)`` in minimal representation, upper-cased.

    The shared suffix is trimmed first, then the shared prefix; each trimmed prefix letter moves
    ``pos`` one to the right. A pure deletion or insertion ends with an empty ``alt`` or ``ref``.
    """
    ref = (ref or "").upper()
    alt = (alt or "").upper()
    pos = int(pos)
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref, alt, pos = ref[1:], alt[1:], pos + 1
    return pos, ref, alt


def same_allele(ref1: str | None, alt1: str | None, ref2: str | None, alt2: str | None) -> bool:
    """True if two ``ref>alt`` encodings written at the same locus are the same variant.

    The caller has already matched the position, so only the trim offset and the trimmed alleles
    are compared. Missing input is never a match.
    """
    if ref1 is None or alt1 is None or ref2 is None or alt2 is None:
        return False
    return trim_alleles(0, ref1, alt1) == trim_alleles(0, ref2, alt2)
