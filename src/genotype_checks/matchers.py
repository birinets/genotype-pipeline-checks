# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 Ruslan Plakhuta
"""One reference rule for "does this data match this recorded variant?", and known-bad ones.

A matcher takes a *record* and a *genotype*:

    record    {"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}    as a variant database writes it
    genotype  {"source": "array" | "imputed",
               "alleles": ["A", "G"],                               the two alleles the subject holds
               "row_ref": "A", "row_alt": "G"}                      imputed or sequenced rows only

An array reports only the letters it observed, so ``row_ref`` and ``row_alt`` are absent. A data
row with more than one alternate allele (read from a VCF) may give ``row_alt`` as a list.

``reference_outcome`` returns one of five outcomes, with a one-line reason:

    match              the data confirms the record as written
    match_after_flip   the letters fit the record only when complemented (A with T, C with G)
    ambiguous_strand   an A/T or C/G site where the two strands give different answers
    ref_alt_swapped    a data row gives the record's reference and alternate the other way round
    no_match           the data does not confirm the record

Only ``match`` lets a pipeline report the record as it stands. ``no_match`` is not a negative result
for the person. ``reference_matches`` is the yes/no form of the same rule: True only for ``match``.

The reference rule is one reasonable reading, not the only one. Its main choices:

* an array genotype can confirm only a single-letter record;
* array letters may be on either strand: at an A/T or C/G site a heterozygous call holds one copy
  of the alternate allele either way, while a homozygous call depends on the strand;
* a data row whose reference agrees with the record is taken as written, even at an A/T or C/G site,
  and is compared by minimal representation; a disagreeing reference is never a match;
* strand flips and swaps are recognised for single-letter records only;
* malformed input never matches.

If your pipeline disagrees on purpose, the case tells you where.
"""

from __future__ import annotations

from typing import Any, NamedTuple, Optional

from .alleles import same_allele

Record = dict[str, Any]
Genotype = dict[str, Any]

MATCH = "match"
MATCH_AFTER_FLIP = "match_after_flip"
AMBIGUOUS_STRAND = "ambiguous_strand"
REF_ALT_SWAPPED = "ref_alt_swapped"
NO_MATCH = "no_match"

#: Every outcome of the reference rule, from the most confirmed to the least.
OUTCOMES = (MATCH, MATCH_AFTER_FLIP, AMBIGUOUS_STRAND, REF_ALT_SWAPPED, NO_MATCH)

_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


class Outcome(NamedTuple):
    outcome: str
    reason: str


# ---------------------------------------------------------------------------------------------
# The reference rule
# ---------------------------------------------------------------------------------------------

def _alts(value: Any) -> tuple:
    """``row_alt`` as a tuple of upper-case alleles. A list stands for a multi-allelic row."""
    if value is None or value == "":
        return ()
    if isinstance(value, str):
        return (value.upper(),)
    return tuple(str(a).upper() for a in value if a)


def _parse(record: Record, genotype: Genotype):
    ref = (record.get("ref") or "").upper()
    alt = (record.get("alt") or "").upper()
    letters = [str(a).upper() for a in (genotype.get("alleles") or []) if a]
    row_ref = (genotype.get("row_ref") or "").upper()
    has_row_alt = genotype.get("row_alt") not in (None, "")
    return ref, alt, letters, row_ref, has_row_alt, _alts(genotype.get("row_alt"))


def _complement(allele: str) -> Optional[str]:
    if allele and all(base in _COMPLEMENT for base in allele):
        return "".join(_COMPLEMENT[base] for base in allele)
    return None


_COPY_WORDS = {0: "no copy", 1: "one copy", 2: "two copies"}


def _copies(n: int, allele: str) -> str:
    words = _COPY_WORDS.get(n, f"{n} copies")
    return f"{words} of {allele}"


def _reference_bases_agree(ref: str, row_ref: str) -> bool:
    """Both start at the same position, so the shorter reference must be a prefix of the longer."""
    n = min(len(ref), len(row_ref))
    return ref[:n] == row_ref[:n]


def reference_outcome(record: Record, genotype: Genotype) -> Outcome:
    """The reference rule: one of ``OUTCOMES`` and a one-line reason. See the module docstring."""
    ref, alt, letters, row_ref, has_row_alt, alts = _parse(record, genotype)
    if not ref or not alt or ref == alt:
        return Outcome(NO_MATCH, "The record needs a reference allele and a different alternate allele.")
    if len(letters) != 2:
        return Outcome(NO_MATCH, f"The kit compares two-allele genotypes, and this one has {len(letters)}.")
    if row_ref and has_row_alt:
        return _row_outcome(ref, alt, letters, row_ref, alts)
    if row_ref or has_row_alt:
        return Outcome(NO_MATCH, "The data row gives only one of its reference and alternate alleles.")
    return _array_outcome(ref, alt, letters)


def reference_matches(record: Record, genotype: Genotype) -> bool:
    """The reference rule as a yes/no answer: True only when the outcome is ``match``."""
    return reference_outcome(record, genotype).outcome == MATCH


def _array_outcome(ref: str, alt: str, letters: list) -> Outcome:
    gt, rec = "/".join(letters), f"{ref}>{alt}"
    if len(ref) != 1 or len(alt) != 1:
        return Outcome(NO_MATCH, f"An array reports single letters, which cannot confirm {rec}: "
                                 "it is not a single-letter change.")
    if any(len(a) != 1 for a in letters):
        return Outcome(NO_MATCH, f"An array genotype is two single letters, and {gt} is not.")
    if ref not in _COMPLEMENT or alt not in _COMPLEMENT:
        return Outcome(NO_MATCH, f"The record {rec} uses a letter other than A, C, G and T.")
    pair = {ref, alt}
    flipped = [_COMPLEMENT.get(a) for a in letters]
    fits_as_written = set(letters) <= pair
    fits_flipped = None not in flipped and set(flipped) <= pair

    if _COMPLEMENT[ref] == alt:  # A/T or C/G: both strands give the same pair of letters
        if not fits_as_written:
            return Outcome(NO_MATCH, f"{gt} holds a letter that is not an allele of {rec} on either strand.")
        n = letters.count(alt)
        if n == 1:
            return Outcome(MATCH, f"{rec} reads the same on both strands, so the letters cannot show the "
                                  f"data's strand. It does not matter here: heterozygous {gt} holds one "
                                  f"copy of {alt} on either strand.")
        return Outcome(AMBIGUOUS_STRAND, f"{rec} reads the same on both strands, and here the strand decides: "
                                         f"{gt} is {_copies(n, alt)} on the record's strand and "
                                         f"{_copies(2 - n, alt)} on the other. The letters cannot tell which.")
    if fits_as_written:
        n = letters.count(alt)
        if n:
            return Outcome(MATCH, f"{gt} holds {_copies(n, alt)}, the record's alternate allele.")
        return Outcome(NO_MATCH, f"{gt} holds only the record's reference allele {ref}.")
    if fits_flipped:
        fgt, n = "/".join(flipped), flipped.count(alt)
        if n:
            return Outcome(MATCH_AFTER_FLIP, f"{gt} fits {rec} only on the other strand: complemented it "
                                             f"reads {fgt}, {_copies(n, alt)}.")
        return Outcome(NO_MATCH, f"{gt} fits {rec} only on the other strand, and complemented it reads "
                                 f"{fgt}, the reference allele only.")
    return Outcome(NO_MATCH, f"{gt} fits {rec} on neither strand: as written and complemented, it holds a "
                             "letter that is not one of the record's alleles.")


def _row_outcome(ref: str, alt: str, letters: list, row_ref: str, alts: tuple) -> Outcome:
    gt, rec = "/".join(letters), f"{ref}>{alt}"
    row = f"{row_ref}>{','.join(alts) if alts else '.'}"

    if _reference_bases_agree(ref, row_ref):
        if not alts:
            return Outcome(NO_MATCH, f"The data row {row} lists no alternate allele, so it cannot confirm {rec}.")
        verb = "describes" if len(alts) == 1 else "includes"
        for row_alt in alts:
            if same_allele(row_ref, row_alt, ref, alt):
                n = letters.count(row_alt)
                if n:
                    return Outcome(MATCH, f"The data row {row} {verb} the record's variant {rec}, and the "
                                          f"genotype {gt} holds {_copies(n, row_alt)}.")
                return Outcome(NO_MATCH, f"The data row {row} {verb} the record's variant {rec}, but the "
                                         f"genotype {gt} holds no copy of it.")
        if len(alts) > 1:
            return Outcome(NO_MATCH, f"None of the alternate alleles in the data row {row} is the record's "
                                     f"variant {rec}.")
        return Outcome(NO_MATCH, f"The data row {row} describes a different variant at this position from "
                                 f"the record {rec}.")

    if len(ref) == 1 and len(alt) == 1 and len(row_ref) == 1 and {ref, alt, row_ref} <= set(_COMPLEMENT):
        swapped = row_ref == alt and ref in alts
        flipped = _COMPLEMENT[row_ref] == ref and _COMPLEMENT[alt] in alts
        if swapped and flipped:  # only possible at an A/T or C/G site
            return Outcome(AMBIGUOUS_STRAND, f"The record gives {ref} as the reference base and the data row "
                                             f"{row} gives {row_ref}. {rec} reads the same on both strands, so "
                                             "the row may be written on the other strand or list the alleles "
                                             "the other way round, and the letters cannot tell which.")
        if flipped:
            complemented = [_complement(a) for a in letters]
            if None not in complemented:
                cgt, n = "/".join(complemented), complemented.count(alt)
                if n:
                    return Outcome(MATCH_AFTER_FLIP, f"Complemented, the data row {row} reads {rec}: the record "
                                                     "and the row are written on opposite strands, or one is "
                                                     f"wrong. Complemented, the genotype {gt} reads {cgt}, "
                                                     f"{_copies(n, alt)}.")
                return Outcome(NO_MATCH, f"Complemented, the data row {row} reads {rec}, but the complemented "
                                         f"genotype {cgt} holds no copy of {alt}.")
        if swapped:
            return Outcome(REF_ALT_SWAPPED, f"The data row {row} gives the record's alternate {alt} as its "
                                            f"reference and the record's reference {ref} as an alternate. Two "
                                            "correct descriptions on one genome build cannot disagree about "
                                            "the reference base: one of them is on another build, lists its "
                                            "alleles in another order, or is wrong.")
        if _COMPLEMENT[row_ref] == alt and _COMPLEMENT[ref] in alts:
            return Outcome(REF_ALT_SWAPPED, f"Complemented, the data row {row} gives the record's alternate "
                                            f"{alt} as its reference and the record's reference {ref} as an "
                                            "alternate: it is on the other strand and lists the alleles the "
                                            "other way round, or it is wrong.")
    return Outcome(NO_MATCH, f"The record {rec} and the data row {row} give different reference bases at the "
                             "same position, so one of them is on another build or coordinate system, or is "
                             "wrong.")


def describes_variant(record: Record, row_ref: str, row_alt: Any) -> bool:
    """True if a data row is about the record's variant rather than another variant at its position.

    Such a row holds the record's two alleles, as written or complemented, in either order. An indel
    counts only as written, compared by minimal representation, with reference bases that agree.
    The VCF check uses this to tell rows that disagree about one variant from rows that describe
    other variants at the same position, such as split multi-allelic rows.
    """
    ref = (record.get("ref") or "").upper()
    alt = (record.get("alt") or "").upper()
    row_ref, alts = (row_ref or "").upper(), _alts(row_alt)
    if not ref or not alt or not row_ref:
        return False
    if _reference_bases_agree(ref, row_ref) and any(same_allele(row_ref, a, ref, alt) for a in alts):
        return True
    if len(ref) != 1 or len(alt) != 1 or len(row_ref) != 1 or not {ref, alt, row_ref} <= set(_COMPLEMENT):
        return False
    c = _COMPLEMENT
    return ((row_ref == alt and ref in alts)             # the record's alleles the other way round
            or (c[row_ref] == ref and c[alt] in alts)     # complemented
            or (c[row_ref] == alt and c[ref] in alts))    # complemented and the other way round


# ---------------------------------------------------------------------------------------------
# Known-bad rules. Each answers yes/no. They exist so the cases can be shown to catch them.
# ---------------------------------------------------------------------------------------------

def _clean(record: Record, genotype: Genotype):
    """Input as version 0.2.0 read it: ``row_alt`` is a single string."""
    ref = (record.get("ref") or "").upper()
    alt = (record.get("alt") or "").upper()
    alleles = [str(a).upper() for a in (genotype.get("alleles") or []) if a]
    row_ref = (genotype.get("row_ref") or "").upper() or None
    row_alt = (genotype.get("row_alt") or "").upper() or None
    return ref, alt, alleles, row_ref, row_alt


def _usable(ref: str, alt: str, alleles: list) -> bool:
    return bool(ref and alt and ref != alt and len(alleles) == 2)


def _array_rule_as_written(ref: str, alt: str, observed: set) -> bool:
    if len(ref) == 1 and len(alt) == 1:
        return alt in observed and observed <= {ref, alt}
    return False


def naive_position_only(record: Record, genotype: Genotype) -> bool:
    """Known-bad: matches whenever the record's alternate allele appears in the genotype.

    This is the position-only join that lets an indel inherit a SNP genotype.
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
    """Known-bad: like the reference rule, but compares a data row's ref and alt letter for letter.

    It rejects every wrong match, so only the case where the same variant is written with a
    different anchor length shows it up.
    """
    ref, alt, letters, row_ref, has_row_alt, alts = _parse(record, genotype)
    if row_ref and has_row_alt and _usable(ref, alt, letters) and _reference_bases_agree(ref, row_ref):
        return any((row_ref, a) == (ref, alt) and a in letters for a in alts)
    return reference_matches(record, genotype)


def naive_trim_only(record: Record, genotype: Genotype) -> bool:
    """Known-bad: like the reference rule, but compares a data row by minimal representation
    without first checking that the reference bases agree.

    Trimming throws away the bases two encodings share, including a base on which they
    disagree. This was this kit's own reference rule in version 0.1.0.
    """
    ref, alt, letters, row_ref, has_row_alt, alts = _parse(record, genotype)
    if row_ref and has_row_alt:
        return _usable(ref, alt, letters) and any(
            same_allele(row_ref, a, ref, alt) and a in letters for a in alts)
    return reference_matches(record, genotype)


def naive_as_written(record: Record, genotype: Genotype) -> bool:
    """Known-bad: takes array letters as written on the record's strand.

    At an A/T or C/G site it reports a homozygous call that, on the other strand, would be the
    opposite homozygote. This was this kit's own reference rule in version 0.2.0.
    """
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
    return _array_rule_as_written(ref, alt, observed)


def _complement_every_letter(allele: Any) -> str:
    return "".join(_COMPLEMENT.get(base, base) for base in str(allele).upper())


def naive_auto_flip(record: Record, genotype: Genotype) -> bool:
    """Known-bad: when the letters do not fit as written, complements every letter and tries again.

    It reports a strand flip as a plain match, and at an A/T or C/G site it turns either
    homozygote into a match.
    """
    if naive_as_written(record, genotype):
        return True
    flipped = dict(genotype)
    flipped["alleles"] = [_complement_every_letter(a) for a in genotype.get("alleles") or []]
    for key in ("row_ref", "row_alt"):
        if genotype.get(key):
            flipped[key] = _complement_every_letter(genotype[key])
    return naive_as_written(record, flipped)


def naive_allele_set(record: Record, genotype: Genotype) -> bool:
    """Known-bad: like the reference rule, but treats a data row's ref and alt as an unordered pair.

    A row that gives the record's alleles the other way round then matches, although the two
    disagree about the reference base. (Records that list an effect and an other allele, rather
    than a reference and an alternate, have to be oriented against the reference genome first.)
    """
    ref, alt, letters, row_ref, has_row_alt, alts = _parse(record, genotype)
    if row_ref and has_row_alt and _usable(ref, alt, letters):
        if alt in letters and any({row_ref, a} == {ref, alt} for a in alts):
            return True
    return reference_matches(record, genotype)


KNOWN_BAD = {
    "position-only": naive_position_only,
    "any-non-reference": naive_any_non_reference,
    "string-equal": naive_string_equal,
    "trim-only": naive_trim_only,
    "as-written": naive_as_written,
    "auto-flip": naive_auto_flip,
    "allele-set": naive_allele_set,
}
