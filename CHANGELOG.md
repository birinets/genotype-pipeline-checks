# Changelog

## Unreleased (0.3.0)

Strand flips, A/T and C/G sites, swapped reference and alternate, and a VCF check. Not released yet; the code says `0.3.0.dev0`.

**Corrected**

- A homozygous call at an A/T or C/G site was matched as written: `T/T` against `A>T` was a match, although on the other strand it is the opposite homozygote. It is now `ambiguous_strand`. The 0.2.0 README listed this as a limit. The old behaviour is kept as the known-bad rule `as-written`. No existing case changes its expectation.
- A single-letter record with a letter other than A, C, G or T (for example `A>N` against `A/N`) is no longer a match. 0.2.0 matched it as written.
- The same holds for data rows and the VCF check: a record or a held allele with a letter other than A, C, G or T, such as N, is `no_match`, because such letters confirm nothing. 0.2.0 matched the record `A>N` against a data row `A>N` holding `N`.

**Added**

- Three outcomes besides `match` and `no_match`: `match_after_flip`, `ambiguous_strand` and `ref_alt_swapped`, each with a one-line reason. `reference_outcome(record, genotype)` returns them; `reference_matches` stays the yes/no form and is True only for `match`.
- Eight cases: `array-strand-flip`, `imputed-strand-flip`, `array-palindromic-homozygous`, `array-palindromic-homozygous-reference`, `array-palindromic-heterozygous`, `imputed-ref-alt-swapped`, `imputed-palindromic-reference-disagrees` and `imputed-palindromic-row-agrees`. Unlike the first twelve, they were written for this kit.
- Known-bad rules `as-written`, `auto-flip` and `allele-set`. Every case whose outcome is not `match` is failed by at least one known-bad rule.
- `genotype-checks vcf`: checks a records VCF against a data VCF (plain or gzip/BGZF, multi-allelic rows included) and writes one TSV row per record with its outcome and reason. It adds the outcome `no_genotype` for positions without a usable two-allele call.
- The VCF check exits with 0 when it ran, whatever the outcomes, and with 2 and a one-line message when it cannot read its input or write `--out`. It refuses a BCF file with a hint to convert it, treats M and MT as the same chromosome, and ignores spaces around a GT value.
- In the VCF check, data rows at one position that describe the record's variant must agree: if their outcomes differ, the record is `no_match` and the reason lists each row. Rows about other variants at that position, such as split multi-allelic rows, do not count against it, and the most confirming outcome is kept.
- `synth` also writes `records.vcf`.
- A data row may give `row_alt` as a list, for a row with several alternate alleles. The case file still uses single strings.
- A draft WDL workflow (`workflows/genotype_checks.wdl`) and `.dockstore.yml`, neither run nor registered, and the design note `docs/v0.3-design.md`.

**Changed**

- `cases.json`: `expect` can be any of five outcomes. Readers that test `expect == "match"` keep working; readers that assume two values must be updated. The file now carries `"schema_version": 2`, separate from the package version; the 0.2.0 file has none, which counts as 1. Its `about` text no longer says array letters are taken as written on the forward strand.
- `check --adapter` and `run_cases` accept True or False, scored as before (yes only for `match`), or an outcome name or an `(outcome, reason)` pair, scored on the exact outcome. A tuple is read by its first element, so `(False, reason)` is a no; 0.2.0 counted any non-empty tuple as yes. A name that is not an outcome stops the run: `run_cases` raises `runner.UnknownOutcome` (a `ValueError`), and `check` exits with 2. `runner.Result` now holds the expected outcome and the answer (`expect`, `got`) instead of two booleans.
- `string-equal` and `trim-only` now build on the new reference rule, so each still fails exactly its one case. `position-only` and `any-non-reference` are unchanged and also fail some of the new cases.
- The README's strand statements now say what is tested and what is not. "Not covered" lost strand flips and swaps for single-letter records, A/T and C/G sites, and rows with several alternate alleles (the VCF check reads them), and gained the limits of the new rules.
- The unit test that pinned "letters on the other strand are not resolved" is replaced by the case `array-strand-flip`.

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
- The README now points to the copy at `indaga-ai/genotype-pipeline-checks`. The source repository stays at `birinets/genotype-pipeline-checks`, and the two are kept identical.

## 0.1.0 (2026-10-11)

First public version: nine cases, a reference rule, two known-bad rules, a synthetic file writer and a toy pipeline.
