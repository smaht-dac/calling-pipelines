#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: Workflow

requirements:
  MultipleInputFeatureRequirement: {}

inputs:
  - id: input_file_vcf_gz
    type: File
    secondaryFiles:
      - .tbi
    doc: Input VCF (bgzipped) with .tbi index

  - id: germline_input_file_vcf_gz
    type: File
    secondaryFiles:
      - .tbi
    doc: Germline SNV calls input VCF (bgzipped) with .tbi index

  - id: genome_reference_fasta
    type: File
    secondaryFiles:
      - ^.dict
      - .fai
    doc: Reference FASTA with index files

  - id: input_files_sr_cram
    type:
      -
        items: File
        type: array
    secondaryFiles:
      - .crai
    doc: Short-read CRAM files (with .crai)

  - id: input_files_pb_cram
    type:
      -
        items: File
        type: array
    secondaryFiles:
      - .crai
    doc: PacBio CRAM files (with .crai)

  - id: input_files_ont_cram
    type:
      -
        items: File
        type: array
    secondaryFiles:
      - .crai
    doc: ONT CRAM files (with .crai)

  - id: sample_id
    type: string
    doc: Output file prefix

  - id: additional_args
    type: string
    default: "-c -C -Q 20 -q 30 -s 0"
    doc: Additional minipileup args (string)

  - id: group_intervals
    type: int
    default: 100
    doc: Group size for interval batching

  - id: sex
    type: string
    doc: Donor sex (male, female, unknown)

  - id: threads
    type: int
    doc: Number of threads to use

outputs:
  output_file_vcf_gz:
    type: File
    outputSource: phase_mosaic_snvs/output_file_vcf_gz

steps:
  minipileup_parallel:
    run: minipileup-parallel.cwl
    in:
      input_file_vcf_gz:
        source: input_file_vcf_gz
      genome_reference_fasta:
        source: genome_reference_fasta
      input_files_sr_cram:
        source: input_files_sr_cram
      input_files_pb_cram:
        source: input_files_pb_cram
      input_files_ont_cram:
        source: input_files_ont_cram
      additional_args:
        source: additional_args
      group_intervals:
        source: group_intervals
    out:
      - output_file_vcf_gz

  tier_filter_variants_SR_PB_ONT:
    run: tier_filter_variants_SR_PB_ONT.cwl
    in:
      input_file_vcf_gz:
        source: input_file_vcf_gz
      minipileup_vcf_gz:
        source: minipileup_parallel/output_file_vcf_gz
    out:
      - output_file_vcf_gz

  phase_mosaic_snvs:
    run: phase_mosaic_snvs.cwl
    in:
      input_file_vcf_gz:
        source: tier_filter_variants_SR_PB_ONT/output_file_vcf_gz
      germline_input_file_vcf_gz:
        source: germline_input_file_vcf_gz
      input_files_pb_cram:
        source: input_files_pb_cram
      sample_id:
        source: sample_id 
      sex:
        source: sex 
      threads:
        source: threads 
    out:
      - output_file_vcf_gz

doc: |
  Filters an SNV VCF file to retain high-confidence variants. |
  Step-2 filters: run minipileup using short-read and long-read CRAM files to compute read support for each SNV, |
  then filter and tier based on read support, strand balance (Fisher) and germline deviation (binomial), |
  then phase/filter based on nearest germline variant
