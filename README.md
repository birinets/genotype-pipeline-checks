# genotype-pipeline-checks

Small checks for one kind of mistake in genotype pipelines: reporting a database variant for a person whose data does not show it.

It has twelve synthetic cases, a reference rule that passes them, four known-bad rules that the cases catch, and a way to run a whole pipeline against made-up subject files. Python standard library only.

## Read this first

- **It is not a clinical validation.** It does not test variant calling, imputation accuracy, annotation or interpretation. It tests one thing: whether a recorded variant is matched to the right genotype.
- **No match is not a negative result.** It means this data does not confirm this record, so a pipeline must not report it. It does not mean the person is free of the variant.
- **Array genotypes are taken as written on the forward strand of the record's genome build.** Strand is not tested: a true carrier reported on the other strand (record `A>G`, array `T/C`) comes back as no match.
- **Twelve cases are a start, not coverage.** Passing means your pipeline avoids these mistakes. See "Not covered" below.
- **The reference rule is one reasonable reading, not the only one.** An array genotype is allowed to confirm only a single-letter record. Rows that carry their own reference and alternate allele must agree with the record about the reference bases they share, and are then compared by minimal representation, which trims shared letters but does not left-align repeats against a genome. If your pipeline differs on purpose, a failing case shows you where.
- **Everything is synthetic.** Coordinates, ids and letters are made up and describe no real person or real variant. Nothing is downloaded.

## The cases

| Case | Does the data match the record? | What it shows |
|---|---|---|
| `delins-vs-array-snp` | no | The record `AA>C` is a deletion-insertion. An array reports two single letters, and `A/C` fits a plain SNP as well as this record. The record must not be reported, and it is not ruled out either. |
| `deletion-anchor-equals-alt` | no | The record `CA>C` is a deletion whose alternate allele is the single letter `C`. A person without the deletion also shows `C`. A rule that looks for the alternate letter reports everyone. |
| `array-alt-mismatch` | no | The genotype holds `G` and `C`. The record's alternate allele is `A`. |
| `imputed-alt-mismatch` | no | The data row's alternate allele is `A`, the record's is `T`. |
| `array-genotype-outside-record-alleles` | no | Genotype `G/T` against a record `A>G`. Forward strand: the record's `G` and a third allele. Reverse strand: `C/A`, no `G`. The array alone cannot say which. |
| `imputed-indel-different-deletion` | no | Both delete bases after the same anchor: one base in the record, two in the data row. |
| `imputed-row-reference-disagrees` | no | The record and the data row give different reference bases at the same coordinates, so one is on another build or coordinate system, or is simply wrong. Do not match, and flag the row. |
| `imputed-reference-disagrees-in-trimmed-base` | no | `CAA>CA` against `CAT>CT`. After trimming both look like the same deletion, but they disagree about the third reference base, which trimming throws away. |
| `array-snv-exact` | yes | Exact single-letter match, one copy. |
| `array-snv-homozygous-alt` | yes | Exact single-letter match, two copies. |
| `imputed-snv-exact` | yes | The data row's alleles equal the record's. |
| `imputed-indel-different-anchor` | yes | The same one-base deletion written with a longer anchor (`CAA>CA` against `CA>C`). Comparing the letters misses it. |

`genotype-checks list` prints each case with its reason. The data is in `src/genotype_checks/cases.json` and does not depend on Python.

## Not covered

- Strand flips, and SNPs whose two letters are complementary (`A/T`, `C/G`), which read the same on both strands. The reference rule still matches these as written; a homozygous call cannot be told from the opposite homozygote.
- Insertion and deletion codes in array files (`I`, `D`) and no-calls (`--`).
- One-letter calls on X, Y and the mitochondrial genome.
- Rows with more than one alternate allele.
- Records and rows with reference and alternate swapped.
- Deletions inside repeats written at different positions (this needs a reference genome).

## Quick start

```bash
git clone https://github.com/indaga-ai/genotype-pipeline-checks
cd genotype-pipeline-checks
PYTHONPATH=src python3 -m genotype_checks selftest
```

Or install it (`pip install .`) and run `genotype-checks selftest`.

`selftest` checks that the reference rule passes all twelve cases, and that four known-bad rules fail exactly the cases they should:

- `position-only` matches whenever the record's alternate allele appears in the genotype. It fails six cases.
- `any-non-reference` matches whenever the genotype shows anything other than the record's reference. It fails all eight "no" cases.
- `string-equal` compares letters instead of minimal representations. It fails only `imputed-indel-different-anchor`.
- `trim-only` compares minimal representations without checking that the reference bases agree. It fails only `imputed-reference-disagrees-in-trimmed-base`. This was the kit's own reference rule in 0.1.0.

## Three ways to use it

**1. A matcher function.** Write a function `matches(record, genotype) -> bool` and point the checker at it:

```bash
genotype-checks check --adapter mypipeline.matching:matches
```

The inputs are plain dictionaries:

```python
record   = {"chrom": "1", "pos": 9000, "ref": "CA", "alt": "C"}      # as a variant database writes it
genotype = {"source": "imputed",                                      # or "array"
            "alleles": ["CAA", "CA"],                                 # the two alleles the subject holds
            "row_ref": "CAA", "row_alt": "CA"}                        # imputed or sequenced rows only
```

An array reports only the letters it observed, so array genotypes have no `row_ref` or `row_alt`. The command exits with 1 if any case fails.

**2. A whole pipeline.** Write the cases out as made-up subject files, run your pipeline on them, and score what it reports. The commands assume `pip install .`; without installing, run them from the repository root with `PYTHONPATH=src` in front and `python3 -m genotype_checks` in place of `genotype-checks`:

```bash
genotype-checks synth out/
# out/array.txt     array-style raw file (rsid, chromosome, position, genotype)
# out/imputed.vcf   one-sample VCF
# out/records.tsv   the variant records, one per case
# out/expected.tsv  the answers

python3 examples/pipeline_demo.py out/ reference > reported.tsv     # a toy pipeline; replace it with yours
genotype-checks compare --reported reported.tsv
```

`reported.tsv` lists the records your pipeline reports for the subject: either the case id, or tab-separated `chrom`, `pos`, `ref`, `alt`. `examples/pipeline_demo.py` also accepts `position-only`, `any-non-reference`, `string-equal` and `trim-only`, to show what failure looks like.

**3. Just the cases.** Read `cases.json` and write the tests in your own language.

## Related tools

This kit covers a narrow gap: matching array letters or VCF rows to a database record without a reference genome. With a VCF and a reference genome, use `bcftools norm` or `vt normalize` to normalise variants, and `rtg vcfeval` or `hap.py` to compare call sets. They do far more than the trimming in `alleles.py`.

## Where the cases come from

They come from a regression test written after a release-blocking bug in the allele-matching step of a commercial product that makes reports from consumer DNA files. A database record for an indel was reported because an array SNP genotype contained one of its letters. The cases were rewritten here as synthetic, engine-independent data, and `alleles.py` is adapted from that product's code, released here by its author.

The maintainer, Ruslan Plakhuta, runs that product ([Indaga](https://indaga.ai)). The product's code, data and history are not part of this repository.

## Contributing

Issues and pull requests are welcome at the source repository, https://github.com/birinets/genotype-pipeline-checks. The copy at https://github.com/indaga-ai/genotype-pipeline-checks is kept identical. The most useful contribution is a new case: say in one line why it matters, and cite a public source or show how it was derived. A disagreement with the reference rule or with a case's explanation is also useful. Version 0.2.0 exists because three explanations in 0.1.0 were wrong and its reference rule had a hole; see `CHANGELOG.md`.

## Licence

Code and data files: Apache-2.0, see `LICENSE`. Prose documentation (this README): CC BY 4.0, https://creativecommons.org/licenses/by/4.0/.
