# Changelog

## 0.2.0 (2026-10-11)

This release corrects wording and cases that were wrong or misleading in 0.1.0.

**Corrected**

- The answer is now `match` or `no_match`. Version 0.1.0 said `carrier` or `not_carrier`, which overstated what the data shows: no match means the data does not confirm the record, not that the person is free of the variant.
- `imputed-indel-different-deletion` in 0.1.0 paired `CA>C` with `CT>C` and called them two different deletions. Two correct records at one position cannot have different reference bases, so that pair is a reference disagreement. It is now its own case, `imputed-row-reference-disagrees`, and the "different deletion" case uses `CA>C` against `CAA>C`.
- `indel-vs-array-snp` is now `delins-vs-array-snp`. `AA>C` is a deletion-insertion, and its explanation now says the record is not ruled out.
- `array-genotype-outside-record-alleles` now explains both strand readings instead of saying the genotype is "not evidence".
- `snv-reference-mismatch` is now `array-alt-mismatch`: the reference matched, the alternate allele differed.
- `array-snv-exact` now uses `A>C` and `imputed-snv-exact` uses `G>A`. The earlier `A>T` and `C>G` read the same on both strands.
- The reference rule in 0.1.0 matched a record and a data row that disagree about a reference base when trimming removed that base (`CAA>CA` against `CAT>CT`). It now checks that the reference bases they share agree. The old behaviour is kept as the known-bad rule `trim-only`.
- Identical reference and alternate alleles, and a row with only one of `row_ref` and `row_alt`, no longer match.
- The trimming is no longer described as the "minimal representation of vt and Tan et al. 2015". It is only the trimming half of that normalisation, without left-alignment.

**Added**

- Cases `deletion-anchor-equals-alt` and `imputed-reference-disagrees-in-trimmed-base`.
- Known-bad rules `any-non-reference` and `trim-only`, so that every "no match" case is failed by at least one known-bad rule.
- README sections "Not covered" and "Related tools", and the statements that no match is not a negative result and that strand is not tested.

**Changed**

- `reference_is_carrier` is now `reference_matches`; `compare --carriers` is now `compare --reported`.
- Removed the "Planned, not promised" section. No scorer or build-conversion check is planned here; `rtg vcfeval` and `hap.py` already compare call sets.
- The repository's home is now `indaga-ai/genotype-pipeline-checks`.

## 0.1.0 (2026-10-11)

First public version: nine cases, a reference rule, two known-bad rules, a synthetic file writer and a toy pipeline.
