# genotype-pipeline-checks

Small checks for one kind of mistake in genotype pipelines: counting someone as a carrier of the wrong variant.

It has nine synthetic cases, a reference rule that passes them, two known-bad rules that the cases catch, and a way to run a whole pipeline against made-up subject files. Python standard library only.

## Read this first

- **It is not a clinical validation.** It does not test variant calling, imputation accuracy, annotation or interpretation. It tests one thing: whether a recorded variant is matched to the right genotype.
- **Nine cases are a start, not coverage.** Passing means your pipeline avoids these nine mistakes. It says nothing about the ones nobody has written down yet.
- **The reference rule is one reasonable reading, not the only one.** An array genotype is allowed to confirm only a single-letter record. Rows that carry their own reference and alternate allele are compared by minimal representation, which trims shared letters but does not left-align repeats against a genome. If your pipeline differs on purpose, a failing case shows you where.
- **Everything is synthetic.** Coordinates, ids and letters are made up and describe no real person or real variant. Nothing is downloaded.

## The cases

| Case | A carrier? | What it shows |
|---|---|---|
| `indel-vs-array-snp` | no | A record `AA>C` is an indel. An array reports the letters `A` and `C`, which also fit a plain SNP. A `C` in the genotype does not make the subject a carrier of the indel. |
| `snv-reference-mismatch` | no | The subject carries `C`. The record's alternate allele is `A`. |
| `imputed-alt-mismatch` | no | The data row's alternate allele is `A`, the record's is `T`. |
| `array-genotype-outside-record-alleles` | no | Genotype `G/T` against a record `A>G`: `T` belongs to neither allele of the record. |
| `imputed-indel-different-deletion` | no | Both rows delete one base at the same position, but different bases (`A` against `T`). |
| `array-snv-exact` | yes | Exact SNV match, one copy. |
| `array-snv-homozygous-alt` | yes | Exact SNV match, two copies. |
| `imputed-snv-exact` | yes | The data row's alleles equal the record's. |
| `imputed-indel-different-anchor` | yes | The same one-base deletion written with a longer anchor (`CAA>CA` against `CA>C`). Comparing the letters misses it. |

`genotype-checks list` prints each case with its reason. The data is in `src/genotype_checks/cases.json` and does not depend on Python.

## Quick start

```bash
git clone https://github.com/birinets/genotype-pipeline-checks
cd genotype-pipeline-checks
PYTHONPATH=src python3 -m genotype_checks selftest
```

Or install it (`pip install .`) and run `genotype-checks selftest`.

`selftest` checks that the reference rule passes all nine cases, and that two known-bad rules fail the cases they should. The first, `position-only`, counts a carrier whenever the record's alternate allele appears in the genotype. It fails four cases. The second, `string-equal`, compares letters instead of minimal representations. It fails only `imputed-indel-different-anchor`.

## Three ways to use it

**1. A matcher function.** Write a function `is_carrier(record, genotype) -> bool` and point the checker at it:

```bash
genotype-checks check --adapter mypipeline.matching:is_carrier
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

python3 examples/pipeline_demo.py out/ reference > carriers.tsv     # a toy pipeline; replace it with yours
genotype-checks compare --carriers carriers.tsv
```

`carriers.tsv` lists the records your pipeline reports as carried: either the case id, or tab-separated `chrom`, `pos`, `ref`, `alt`. `examples/pipeline_demo.py` also accepts `position-only` and `string-equal`, to show what failure looks like.

**3. Just the cases.** Read `cases.json` and write the tests in your own language.

## Where the cases come from

They come from a regression test written after a release-blocking bug in the allele-matching step of a commercial genotype interpretation product. A database record for an indel was counted as carried because an array SNP genotype contained one of its letters. The cases were rewritten here as synthetic, engine-independent data, and `alleles.py` is adapted from that product's code, released here by its author.

The maintainer, Ruslan Plakhuta, runs that product ([Indaga](https://indaga.ai)). The product's code, data and history are not part of this repository.

## Planned, not promised

- A scorer that checks calls against the public Genome in a Bottle truth files, downloading them at run time and verifying checksums.
- A check for converting coordinates between genome builds, with a slot for any converter.

## Contributing

Issues and pull requests are welcome. The most useful contribution is a new case: say in one line why it matters, and cite a public source or show how it was derived. A disagreement with the reference rule is also useful. Please open an issue and explain which reading you think is right.

## Licence

Code and data files: Apache-2.0, see `LICENSE`. Prose documentation (this README): CC BY 4.0, https://creativecommons.org/licenses/by/4.0/.
