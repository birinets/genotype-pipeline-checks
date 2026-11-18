# genotype-pipeline-checks

Small checks for one kind of mistake in genotype pipelines: reporting a database variant for a person whose data does not show it.

It has twenty synthetic cases, a reference rule that passes them, seven known-bad rules that the cases catch, a way to run a whole pipeline against made-up subject files, and a check of the variants a pipeline may report (a records VCF) against a person's VCF. Python standard library only.

## Read this first

- **It is not a clinical validation.** It does not test variant calling, imputation accuracy, annotation or interpretation. It tests one thing: whether a recorded variant is matched to the right genotype.
- **Only `match` lets a pipeline report a record as it stands.** `no_match` means this data does not confirm this record, so a pipeline must not report it. It does not mean the person is free of the variant. The other three outcomes say why the data cannot confirm the record as written. No outcome other than `match` is a negative result: `ambiguous_strand` and `match_after_flip` can hide one or two copies, and `no_genotype` says nothing at all.
- **Strand is tested for single-letter records only.** Array letters may be on either strand: letters that fit the record only when complemented, and then hold its alternate allele, are `match_after_flip`, and a homozygous call at an A/T or C/G site, whose two letters are the same pair on both strands, is `ambiguous_strand`. A data row whose reference agrees with the record is taken as written on the forward strand, as the VCF specification expects (REF is the reference genome's base at POS); a row with the record's alleles the other way round is `ref_alt_swapped`, or `ambiguous_strand` at an A/T or C/G site. Strand for indels, and anything that needs the array manifest, allele frequencies or a reference genome, is not tested.
- **Twenty cases are a start, not coverage.** Passing means your pipeline gives the expected answer on these twenty inputs, not that it avoids these mistakes everywhere. See "Not covered" below.
- **The reference rule is one reasonable reading, not the only one.** An array genotype is allowed to confirm only a single-letter record. A heterozygous call at an A/T or C/G site is a match, because it holds one copy of the alternate allele on either strand; a stricter reading would flag every such site. This is a per-genotype answer: at one A/T or C/G site a one-copy call is reported and a two-copy call is not. Score tools decide per site for everyone; do not use these outcomes to filter genotypes for a score. Rows that carry their own reference and alternate allele must agree with the record about the reference bases they share, and are then compared by minimal representation, which trims shared letters but does not left-align repeats against a genome. A record or held allele with a letter other than A, C, G and T, such as N, never matches. If your pipeline differs on purpose, a failing case shows you where.
- **The cases are synthetic.** Their coordinates, ids and letters are made up and describe no real person or real variant. Nothing is downloaded. The VCF check reads the files you give it, and its output is as sensitive as they are.

## The outcomes

| Outcome | Meaning | What a pipeline should do |
|---|---|---|
| `match` | The data confirms the record as written. | It may report the record. |
| `match_after_flip` | The letters fit the record only when complemented (A with T, C with G), and complemented they hold the record's alternate allele: the data and the record are written on opposite strands, or one of them is wrong. Complemented letters that hold only the reference are `no_match`, and only the reason names the flip. | Count it only if the pipeline corrects strand on purpose, for example from the array manifest. Never report it as a plain match. |
| `ambiguous_strand` | An A/T or C/G site where the letters cannot show the strand and it could change the answer: a homozygous array call, or a data row that gives the record's two letters the other way round, whatever its genotype, because a row whose reference disagrees is never a match. | Do not report it from these letters. The array manifest's strand, or the data provider's stated strand, could resolve it. Allele frequencies cannot for one person: they can only suggest a strand for many people typed on the same array, and tools treat even that as unreliable. This kit uses neither. |
| `ref_alt_swapped` | A data row gives the record's reference as its alternate and the record's alternate as its reference, as written or after complementing. At an A/T or C/G site such a row is `ambiguous_strand`. Two descriptions that both give the genome's reference base at a position cannot disagree about it. | Do not report it. Check the genome build and the allele order of both files. |
| `no_match` | The data does not confirm the record. | Do not report it. It is not a negative result for the person. |

The VCF check adds `no_genotype`: the data has no usable two-allele genotype at the record's position.

## The cases

| Case | Outcome | What it shows |
|---|---|---|
| `delins-vs-array-snp` | no_match | The record `AA>C` is a deletion-insertion. An array reports two single letters, and `A/C` fits a plain SNP as well as this record. The record must not be reported, and it is not ruled out either. |
| `deletion-anchor-equals-alt` | no_match | The record `CA>C` is a deletion whose alternate allele is the single letter `C`. A person without the deletion also shows `C`. A rule that looks for the alternate letter reports everyone. |
| `array-alt-mismatch` | no_match | The genotype holds `G` and `C`. The record's alternate allele is `A`. |
| `imputed-alt-mismatch` | no_match | The data row's alternate allele is `A`, the record's is `T`. |
| `array-genotype-outside-record-alleles` | no_match | Genotype `G/T` against a record `A>G`. As written: the record's `G` and a third allele. Complemented: `C/A`, no `G`. The array alone cannot say which. |
| `imputed-indel-different-deletion` | no_match | Both delete bases after the same anchor: one base in the record, two in the data row. |
| `imputed-row-reference-disagrees` | no_match | The record and the data row give different reference bases at the same coordinates, so one is on another build or coordinate system, or is simply wrong. Do not match, and flag the row. |
| `imputed-reference-disagrees-in-trimmed-base` | no_match | `CAA>CA` against `CAT>CT`. After trimming both look like the same deletion, but they disagree about the third reference base, which trimming throws away. |
| `array-snv-exact` | match | Exact single-letter match, one copy. |
| `array-snv-homozygous-alt` | match | Exact single-letter match, two copies. |
| `imputed-snv-exact` | match | The data row's alleles equal the record's. |
| `imputed-indel-different-anchor` | match | The same one-base deletion written with a longer anchor (`CAA>CA` against `CA>C`). Comparing the letters misses it. |
| `array-strand-flip` | match_after_flip | `T/C` against `A>G`. Complemented, the letters read `A/G`, so they fit only if the data and the record are on opposite strands. |
| `imputed-strand-flip` | match_after_flip | The same flip in a data row: `T>C` complemented reads `A>G`, the record. |
| `array-palindromic-homozygous` | ambiguous_strand | `T/T` against `A>T`: two copies of `T` on the record's strand, none on the other. Version 0.2.0 called this a match. |
| `array-palindromic-homozygous-reference` | ambiguous_strand | `C/C` against `C>G`: none or two copies. A rule that complements whenever the letters as written are not a match turns it into a homozygous match. |
| `array-palindromic-heterozygous` | match | `G/C` against `G>C`: one copy of `C` on either strand, so the strand does not matter. |
| `imputed-ref-alt-swapped` | ref_alt_swapped | The data row `G>A` has the record's alleles the other way round. `G/G` is two copies of the record's alternate by letter, and homozygous reference by the row's own labels. |
| `imputed-palindromic-reference-disagrees` | ambiguous_strand | The data row `T>A` against `A>T` may be on the other strand or swapped, and for `T/T` the two readings give opposite answers. |
| `imputed-palindromic-row-agrees` | match | The data row `C>G` agrees with the record about the reference, so the kit takes it as written on the forward strand: `G/G` is two copies of `G`. A row on the other strand with its alleles also swapped would look the same (see Not covered). |

`genotype-checks list` prints each case with its reason. The data is in `src/genotype_checks/cases.json` and does not depend on Python.

## Not covered

- Strand and swaps for indels and multi-letter changes. An indel on the other strand sits at a different position, and an insertion and a deletion at one position can be two real, different variants.
- Resolving A/T and C/G sites with the array manifest's strand information, or, for a whole dataset typed on one array, with allele frequencies.
- A data row on the other strand that also has its alleles swapped. At an A/T or C/G site it reads like an agreeing row, and the kit takes it as one. In VCFs made from array data without checking REF against a genome this is one slip, not two: PLINK writes its A2 allele as REF (flagged `PR`), and Illumina TOP-strand letters are not always the forward strand. Run `bcftools +fixref` on such files first.
- Letters or a data row that fit the record only when complemented but hold no alternate allele are `no_match`; only the reason shows the strand problem, so outcome counts understate strand problems in a file.
- Sites with more than two alleles. Strand reasoning assumes the site has only the record's two alleles; at a site with more, letters that fit only when complemented could be two other alleles on the record's strand.
- Insertion and deletion codes in array files (`I`, `D`) and no-calls (`--`).
- One-letter calls on X, Y and the mitochondrial genome. The VCF check reports them as `no_genotype`.
- A case for data rows with more than one alternate allele. The VCF check reads such rows; the case file cannot express them yet.
- A genome-build check. The kit compares a record with the data, never either of them with a reference genome. Many reference disagreements usually mean the files are on different builds; many swaps more often mean that one file orders its alleles another way (major and minor, PLINK's A1 and A2, or effect and other allele).
- Call quality in the VCF check: FILTER, GQ, DP, GP, DS and imputation quality are not read.
- Deletions inside repeats written at different positions (this needs a reference genome), and data rows that start before the record's position, such as a deletion that spans it.

## Quick start

```bash
git clone https://github.com/indaga-ai/genotype-pipeline-checks
cd genotype-pipeline-checks
PYTHONPATH=src python3 -m genotype_checks selftest
```

Or install it (`pip install .`) and run `genotype-checks selftest`.

`selftest` checks that the reference rule gives every case its outcome, and that seven known-bad rules fail exactly the cases they should:

- `position-only` matches whenever the record's alternate allele appears in the genotype. It fails nine cases.
- `any-non-reference` matches whenever the genotype shows anything other than the record's reference. It fails every case whose outcome is not `match` except `array-palindromic-homozygous-reference`, whose letters are all the record's reference.
- `string-equal` compares letters instead of minimal representations. It fails only `imputed-indel-different-anchor`.
- `trim-only` compares minimal representations without checking that the reference bases agree. It fails only `imputed-reference-disagrees-in-trimmed-base`. This was the kit's own reference rule in 0.1.0.
- `as-written` takes array letters as written on the record's strand. It fails only `array-palindromic-homozygous`. This was the kit's own reference rule in 0.2.0.
- `auto-flip` complements every letter whenever the letters as written are not a match, and tries again. It fails both flip cases and both homozygous A/T and C/G array cases.
- `allele-set` treats a data row's reference and alternate as an unordered pair. It fails `imputed-ref-alt-swapped` and `imputed-palindromic-reference-disagrees`.

## Four ways to use it

**1. A matcher function.** Write a function `matches(record, genotype)` and point the checker at it:

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

An array reports only the letters it observed, so array genotypes have no `row_ref` or `row_alt`.

The function can answer in two ways. **True or False**, where True means "report this record": it passes a case when it says True exactly for the cases whose outcome is `match`, so every other outcome expects False. Or **an outcome name** such as `"match_after_flip"`, or an `(outcome, reason)` pair: it is then scored on the exact outcome. A tuple is read by its first element, so `(False, reason)` is a no. The kit's own rule is available in both forms, `genotype_checks.matchers:reference_matches` and `genotype_checks.matchers:reference_outcome`. The command exits with 1 if any case fails, and with 2 if the function answers with a name that is not one of the five outcomes.

**2. A whole pipeline.** Write the cases out as made-up subject files, run your pipeline on them, and score what it reports. The commands assume `pip install .`; without installing, run them from the repository root with `PYTHONPATH=src` in front and `python3 -m genotype_checks` in place of `genotype-checks`:

```bash
genotype-checks synth out/
# out/array.txt     array-style raw file (rsid, chromosome, position, genotype)
# out/imputed.vcf   one-sample VCF
# out/records.tsv   the variant records, one per case
# out/records.vcf   the same records as a VCF, for the VCF check below
# out/expected.tsv  the outcomes

python3 examples/pipeline_demo.py out/ reference > reported.tsv     # a toy pipeline; replace it with yours
genotype-checks compare --reported reported.tsv
```

`reported.tsv` lists the records your pipeline reports for the subject: either the case id, or tab-separated `chrom`, `pos`, `ref`, `alt`. `examples/pipeline_demo.py` also accepts the name of any known-bad rule, to show what failure looks like. `compare` scores every reported record as a match, so a pipeline that corrects strand on purpose and reports the corrected records fails the two flip cases here: from a list of reported records, a deliberate correction cannot be told from a silent one.

**3. A VCF check.** Check the variants your pipeline may report against one person's VCF before reporting them:

```bash
genotype-checks vcf --records records.vcf.gz --data person.vcf.gz --out outcomes.tsv
```

The outcomes answer one question: may this record be reported as present for this person? They are not a filter for a polygenic score. A well-called homozygous reference row is `no_match` here but a valid dose of zero in a score, and A/T and C/G sites are judged per genotype here, while score tools drop or resolve them per site for everyone.

Both files may be plain or gzip/BGZF; convert a BCF file with `bcftools view` first. From the records VCF the check reads CHROM, POS, ID, REF and ALT, and checks each alternate allele as its own record. From the data VCF it reads the GT field; if the file has several samples, name one with `--sample`. A record is compared with the data rows that start at the same position, with or without a leading `chr` in the chromosome name (M and MT are the same chromosome). Rows with several alternate alleles are read as written.

The check reads GT only. It ignores FILTER, genotype quality, depth, genotype probabilities and imputation quality, so in this mode `match` means the called letters agree with the record, not that the call is reliable. Filter the data first, for example with `bcftools view -f PASS` and an imputation-quality threshold for imputed files.

If several data rows start at the record's position, each is compared. When two or more of them describe the record's variant (its two alleles, as written or complemented) and give different outcomes, the record is `no_match` and the reason lists what each row said; a row without a usable call counts as disagreeing. Rows that describe other variants at the position, such as split multi-allelic rows or an SNV beside an indel, do not count against it: the most confirming outcome is kept, and the reason names the others.

The output is a TSV with one row per record: `chrom`, `pos`, `id`, `ref`, `alt`, `outcome`, `genotype` (the data's alleles as letters) and `reason`. Besides the five outcomes, `no_genotype` means the data has no usable two-allele genotype there: no row, a missing or half call, or a one-allele call. A missing row is not a reference call. A count of outcomes goes to standard error. The output contains genotypes from the data file, so treat it like that file. The command exits with 0 when it ran, whatever the outcomes, and with 2 when it is called wrongly or cannot read its input or write its output.

The check reads the whole data file, which is slow for large multi-sample files; subset those to the records' positions and one sample first, as the workflow below does with `bcftools view`. Try it on the synthetic files with `genotype-checks synth out/` and then `genotype-checks vcf --records out/records.vcf --data out/imputed.vcf`.

**4. Just the cases.** Read `cases.json` and write the tests in your own language. A case's `expect` is one of the five outcomes; if your tests only ask "report or not", treat `match` as yes and everything else as no. The file's `schema_version` is 2, the version in which `expect` gained the three outcomes besides `match` and `no_match`. The 0.2.0 file has no `schema_version`; read a missing one as 1. The kit's own loader reads either.

## A workflow (draft)

`workflows/genotype_checks.wdl` (WDL 1.0) wraps the VCF check for workflow engines such as Cromwell. It takes one or more bgzipped, indexed data files (for example one per chromosome), cuts each down with bcftools to the records' positions and to one sample, joins the pieces, and runs the check. It also asks for the genome build of the records and the data, which must be the same: the kit matches chromosome and position only. The check runs in an image with the kit installed, built from the `Dockerfile` by `.github/workflows/image.yml` for release tags only, so nothing is installed while the workflow runs. Both images are pinned by digest, the kit's own once it has been built.

`.dockstore.yml` describes the workflow for Dockstore and reads release tags only. Its test inputs, `workflows/test.json`, are the synthetic files in `workflows/test-data/` (written by `genotype-checks synth`, then compressed with `bgzip` and indexed with `tabix`), linked as raw GitHub files under the v0.3.0 tag.

All of this is a draft: the workflow has not been run on a workflow engine, and nothing is registered anywhere. The image is built for each release tag.

## Related tools

This kit covers a narrow gap: matching array letters or VCF rows to a database record without a reference genome. With a VCF and a reference genome, use `bcftools norm` or `vt normalize` to normalise variants, and `rtg vcfeval` or `hap.py` to compare call sets. They do far more than the trimming in `alleles.py`.

For strand and reference alleles, `bcftools norm --check-ref` checks REF alleles against a reference genome, and the `bcftools +fixref` plugin can flip or swap REF and ALT to match it; its `flip` mode leaves A/T and C/G sites alone. This kit uses no reference genome, so it can say that a record and the data disagree, but not which of them is right.

Tools that fix strand against a reference already exist, and most do more than this kit: [Genotype Harmonizer](https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer) and [conform-gt](https://faculty.washington.edu/browning/conform-gt.html) align a dataset to a reference panel and resolve A/T and C/G sites from allele frequency or linkage disequilibrium; Will Rayner's [HRC/1000G check tool](https://www.chg.ox.ac.uk/~wrayner/tools/) compares a PLINK dataset with a panel's sites and writes flip, swap and exclusion lists; [snpflip](https://github.com/biocore-ntnu/snpflip) finds reverse-strand and ambiguous SNPs from a reference FASTA; PLINK's `--flip` and `--ref-from-fa` apply such fixes; the Michigan and TOPMed imputation servers check allele switches and strand flips at upload ([Michigan](https://genepi.github.io/michigan-imputationserver/pipeline/), [TOPMed](https://topmedimpute.readthedocs.io/en/latest/pipeline/)); and [pgsc_calc](https://pgsc-calc.readthedocs.io/) labels every scoring-file match by allele order, strand flip and strand ambiguity before it scores. This kit does none of that fixing. What it adds is a reference-free check of one record against one person's data, and the synthetic cases, which test whether your own matching step draws the same distinctions.

pgsc_calc's matching log ([its output documentation](https://pgsc-calc.readthedocs.io/en/latest/explanation/output.html)) describes how each scoring-file variant matched a target variant, for all samples at once. This kit's outcomes also depend on one person's genotype: a `0/0` row is `no_match` here and a match with a dose of zero there. So the outcomes only roughly correspond to its labels:

| This kit | Roughly corresponds in pgsc_calc's log |
|---|---|
| `match` | A `match_type` of `refalt` or `altref`, without `_flip`. |
| `match_after_flip` | A `match_type` ending in `_flip` (`match_flipped` is true). |
| `ambiguous_strand` | `ambiguous` is true. |
| `ref_alt_swapped` | No label of its own. pgsc_calc accepts both allele orders, because a scoring file names an effect allele and an other allele, not a reference and an alternate. A records VCF names the reference, so here the other order is a disagreement. |
| `no_match`, `no_genotype` | A `match_status` of `unmatched` when no target variant matches. A homozygous reference call or a missing call has no label there, because pgsc_calc matches variants, not genotypes. |

pgsc_calc drops matches to A/T and C/G sites unless its [`--keep_ambiguous`](https://pgsc-calc.readthedocs.io/en/latest/reference/params.html) option is set, and Will Rayner's tool removes A/T and G/C SNPs with a minor allele frequency above 0.4. Both decide per site, for everyone in the data. The TOPMed server's documentation says its allele-switch check ignores A/T and C/G sites. This kit decides per genotype instead.

## Where the cases come from

The first twelve cases come from a regression test written after a release-blocking bug in the allele-matching step of a commercial product that makes reports from consumer DNA files. A database record for an indel was reported because an array SNP genotype contained one of its letters. The cases were rewritten here as synthetic, engine-independent data, and `alleles.py` is adapted from that product's code, released here by its author. The maintainer, Ruslan Plakhuta, runs that product ([Indaga](https://indaga.ai)). The product's code, data and history are not part of this repository.

The eight strand and swap cases added in 0.3 were written for this kit. They follow distinctions that strand-fixing tools such as `bcftools +fixref` also draw: flipped alleles, swapped reference and alternate, and A/T or C/G sites that the letters alone cannot resolve.

## Contributing

Issues and pull requests are welcome at the source repository, https://github.com/birinets/genotype-pipeline-checks. The copy at https://github.com/indaga-ai/genotype-pipeline-checks is kept identical. The most useful contribution is a new case: say in one line why it matters, and cite a public source or show how it was derived. A disagreement with the reference rule or with a case's explanation is also useful. Version 0.2.0 exists because three explanations in 0.1.0 were wrong and its reference rule had a hole; see `CHANGELOG.md`.

## Licence

Code and data files: Apache-2.0, see `LICENSE`. Prose documentation (this README and the files in `docs/`): CC BY 4.0, https://creativecommons.org/licenses/by/4.0/.
