version 1.0

# genotype-pipeline-checks, VCF mode as a workflow: check the variants a pipeline will score (a records
# VCF) against one person's VCF before scoring, and write one outcome per record.
#
# DRAFT. This file has not been run or syntax-checked yet. Use it with public or synthetic data only.
# The default kit_source points at a release tag that does not exist yet.

workflow genotype_checks {
  input {
    File records_vcf
    File data_vcf
    String? sample
    String docker_image = "python:3.12-slim"
    String kit_source = "https://github.com/birinets/genotype-pipeline-checks/archive/refs/tags/v0.3.0.tar.gz"
  }

  call check_vcf {
    input:
      records_vcf = records_vcf,
      data_vcf = data_vcf,
      sample = sample,
      docker_image = docker_image,
      kit_source = kit_source
  }

  output {
    File outcomes = check_vcf.outcomes
    File summary = check_vcf.summary
  }

  meta {
    author: "Ruslan Plakhuta"
    description: "Checks allele identity and strand between a records VCF and a data VCF before scoring. Writes one TSV row per record with an outcome (match, match_after_flip, ambiguous_strand, ref_alt_swapped, no_match or no_genotype) and a short reason. Synthetic test data only in this repository."
  }

  parameter_meta {
    records_vcf: "VCF, plain or gzip/BGZF, of the variants a pipeline will score. CHROM, POS, ID, REF and ALT are read; sample columns are ignored."
    data_vcf: "VCF, plain or gzip/BGZF, with the person's genotypes in GT."
    sample: "Sample column to read. Required when data_vcf has more than one sample."
    docker_image: "Container image with Python 3.9 or later and pip. Default: the official python:3.12-slim image."
    kit_source: "pip-installable source of genotype-pipeline-checks, pinned to a release."
    outcomes: "TSV: chrom, pos, id, ref, alt, outcome, genotype, reason. It contains genotypes from data_vcf."
    summary: "One line counting the outcomes."
  }
}

task check_vcf {
  input {
    File records_vcf
    File data_vcf
    String? sample
    String docker_image
    String kit_source
    Int extra_disk_gb = 10
  }

  Int disk_gb = ceil(size(records_vcf, "GB") + size(data_vcf, "GB")) + extra_disk_gb

  command <<<
    set -euo pipefail
    python3 -m pip install --quiet --no-cache-dir --disable-pip-version-check "~{kit_source}"
    genotype-checks vcf \
      --records "~{records_vcf}" \
      --data "~{data_vcf}" \
      ~{"--sample " + sample} \
      --out outcomes.tsv 2> summary.txt || { cat summary.txt >&2; exit 1; }
    cat summary.txt >&2
  >>>

  output {
    File outcomes = "outcomes.tsv"
    File summary = "summary.txt"
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "2 GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 2
  }
}
