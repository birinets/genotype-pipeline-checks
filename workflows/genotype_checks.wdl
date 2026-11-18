version 1.0

# genotype-pipeline-checks, VCF mode as a workflow: check the variants a pipeline may report (a records
# VCF) against one person's VCF before reporting them, and write one outcome per record.
#
# DRAFT. It has not been run on a workflow engine yet. Use it with public or synthetic data only.
# Its default image is built from the v0.3.0 tag by .github/workflows/image.yml.

workflow genotype_checks {
  input {
    File records_vcf
    Array[File] data_vcfs
    Array[File] data_vcf_indexes
    String genome_build
    String? sample
    # Built from this repository by .github/workflows/image.yml when the v0.3.0 tag is pushed. Once it is
    # built, this tag is replaced by the image's digest (ghcr.io/birinets/genotype-pipeline-checks@sha256:...).
    String docker_image = "ghcr.io/birinets/genotype-pipeline-checks:0.3.0"
    # bcftools 1.23.1 from BioContainers (tag 1.23.1--hb2cee57_0), pinned by digest.
    String bcftools_image = "quay.io/biocontainers/bcftools@sha256:8cf74dd0738e0a86fa3f8dc598596e8160e7e4e71380ab6db943a3b53f993869"
  }

  scatter (files in zip(data_vcfs, data_vcf_indexes)) {
    call subset_vcf {
      input:
        records_vcf = records_vcf,
        data_vcf = files.left,
        data_vcf_index = files.right,
        sample = sample,
        docker_image = bcftools_image
    }
  }

  call concat_vcfs {
    input:
      vcfs = subset_vcf.subset,
      docker_image = bcftools_image
  }

  call check_vcf {
    input:
      records_vcf = records_vcf,
      data_vcf = concat_vcfs.data,
      sample = sample,
      genome_build = genome_build,
      docker_image = docker_image
  }

  output {
    File outcomes = check_vcf.outcomes
    File summary = check_vcf.summary
  }

  meta {
    author: "Ruslan Plakhuta"
    description: "Checks allele identity and strand between a records VCF and one person's VCF before a pipeline reports the records. Cuts the data down to the records' positions with bcftools, then writes one TSV row per record with an outcome (match, match_after_flip, ambiguous_strand, ref_alt_swapped, no_match or no_genotype) and a short reason. Synthetic test data only in this repository."
  }

  parameter_meta {
    records_vcf: "VCF, plain or gzip/BGZF, of the variants a pipeline may report. CHROM, POS, ID, REF and ALT are read; sample columns are ignored. Its positions also choose the data rows to keep."
    data_vcfs: "One or more VCF or BCF files, bgzipped and indexed, with the person's genotypes in GT, for example one per chromosome. Each is cut down with bcftools to the records' positions (and to sample, if given) before the check."
    data_vcf_indexes: "The .tbi or .csi index of each data file, in the same order as data_vcfs."
    genome_build: "The genome build of both the records and the data, such as GRCh38. They must be on the same build: the kit matches chromosome and position only and cannot detect a build mismatch. The value is written into the summary; nothing checks it."
    sample: "Sample to keep and read. Required when the data files have more than one sample. A name containing a single quote is not supported."
    docker_image: "Image with genotype-pipeline-checks installed, built from this repository."
    bcftools_image: "Image with bcftools, for cutting down and joining the data files."
    outcomes: "TSV: chrom, pos, id, ref, alt, outcome, genotype, reason. It contains genotypes from the data files."
    summary: "The genome build as given, and one line counting the outcomes."
  }
}

task subset_vcf {
  input {
    File records_vcf
    File data_vcf
    File data_vcf_index
    String? sample
    String docker_image
    Int extra_disk_gb = 10
  }

  Int disk_gb = ceil(2 * size(data_vcf, "GB") + size(records_vcf, "GB")) + extra_disk_gb

  command <<<
    set -euo pipefail

    # bcftools looks for the index next to the data file, so link both into the working directory.
    name="$(basename "~{data_vcf}")"
    index="~{data_vcf_index}"
    ln -s "~{data_vcf}" "$name"
    case "$index" in
      *.csi) ln -s "$index" "$name.csi" ;;
      *) ln -s "$index" "$name.tbi" ;;
    esac

    # The records' positions under each spelling of their chromosome (with and without "chr", and M or
    # MT), because bcftools matches chromosome names exactly and the kit does not.
    bcftools query --format '%CHROM\t%POS\n' "~{records_vcf}" > positions.tsv
    while IFS=$'\t' read -r chrom pos; do
      key="${chrom#[Cc][Hh][Rr]}"
      case "$key" in
        [Mm] | [Mm][Tt]) keys="M MT" ;;
        *) keys="$key" ;;
      esac
      for k in $keys; do
        printf '%s\t%s\nchr%s\t%s\n' "$k" "$pos" "$k" "$pos"
      done
    done < positions.tsv > regions.tsv

    bcftools view ~{if defined(sample) then "--samples '" + select_first([sample]) + "'" else ""} \
      --regions-file regions.tsv --output-type z --output subset.vcf.gz "$name"
  >>>

  output {
    File subset = "subset.vcf.gz"
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "2 GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}

task concat_vcfs {
  input {
    Array[File] vcfs
    String docker_image
    Int extra_disk_gb = 10
  }

  Int disk_gb = ceil(2 * size(vcfs, "GB")) + extra_disk_gb

  command <<<
    set -euo pipefail
    # One file per chromosome or region, in any order of chromosomes.
    bcftools concat --file-list "~{write_lines(vcfs)}" --output-type z --output data.vcf.gz
  >>>

  output {
    File data = "data.vcf.gz"
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "2 GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}

task check_vcf {
  input {
    File records_vcf
    File data_vcf
    String? sample
    String genome_build
    String docker_image
    Int extra_disk_gb = 10
  }

  Int disk_gb = ceil(size(records_vcf, "GB") + size(data_vcf, "GB")) + extra_disk_gb

  command <<<
    set -euo pipefail
    printf 'Genome build given for the records and the data: %s. The kit does not check it.\n' \
      '~{genome_build}' > summary.txt
    genotype-checks vcf \
      --records "~{records_vcf}" \
      --data "~{data_vcf}" \
      ~{if defined(sample) then "--sample '" + select_first([sample]) + "'" else ""} \
      --out outcomes.tsv 2>> summary.txt || { cat summary.txt >&2; exit 1; }
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
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}
