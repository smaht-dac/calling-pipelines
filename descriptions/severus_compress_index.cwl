#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: Workflow

requirements:
  MultipleInputFeatureRequirement: {}

inputs:
  - id: input_file_cram
    type: File
    secondaryFiles:
      - .crai
    doc: Primary input CRAM (sorted + indexed)

  # Reference file
  - id: tandem_repeats_bed
    type: File
    doc: BED file with tandem repeat regions

  # Optional arguments
  - id: nthreads
    type: int
    default: 16
    doc: Number of threads to use [16]

outputs:
  output_file_vcf_gz:
    type: File
    outputSource: compress_index_vcf/output_file_vcf_gz

steps:
  severus:
    run: severus.cwl
    in:
      input_file_cram:
        source: input_file_cram
      tandem_repeats_bed:
        source: tandem_repeats_bed
      nthreads:
        source: nthreads
    out:
      - output_file_vcf

  compress_index_vcf:
    run: compress_index_vcf.cwl
    in:
      input_file_vcf:
        source: severus/output_file_vcf
    out:
      - output_file_vcf_gz

doc: |
    Long-read SNV and Indel calling using Severus, followed by compression and indexing of the resulting VCF file
