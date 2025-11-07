#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: Workflow

requirements:
  MultipleInputFeatureRequirement: {}

inputs:
  # Required primary positional BAM/CRAM
  - id: input_file_cram
    type: File
    secondaryFiles:
      - .crai
    doc: Primary input CRAM (sorted + indexed)

  # Optional extra CRAM files
  - id: additional_input_files_cram
    default: null
    type:
      -
        items: File
        type: array
    secondaryFiles:
      - .crai
    doc: Additional CRAM files (sorted + indexed)

  # Required reference files
  - id: genome_reference_fasta
    type: File
    secondaryFiles:
      - ^.dict
      - .fai
    doc: Reference FASTA with index files

  - id: mei_reference_fasta
    type: File
    doc: Reference FASTA for MEI (Mobile Element Insertion)

  # Required arguments
  - id: sample_name
    type: string
    doc: Sample name in VCF

  # Optional arguments
  - id: platform
    type: string
    default: "--hifi"
    doc: Sequencing platform flag --hifi / --ont [--hifi]

  - id: nthreads
    type: int
    default: 8
    doc: Number of threads to use [8]

outputs:
  output_file_vcf_gz:
    type: File
    outputSource: compress_index_vcf/output_file_vcf_gz

steps:
  longcallD:
    run: longcallD.cwl
    in:
      input_file_cram:
        source: input_file_cram
      additional_input_files_cram:
        source: additional_input_files_cram
      genome_reference_fasta:
        source: genome_reference_fasta
      mei_reference_fasta:
        source: mei_reference_fasta
      sample_name:
        source: sample_name
      platform:
        source: platform
      nthreads:
        source: nthreads
    out:
      - output_file_vcf

  compress_index_vcf:
    run: compress_index_vcf.cwl
    in:
      input_file_vcf:
        source: longcallD/output_file_vcf
    out:
      - output_file_vcf_gz

doc: |
    Long-read SNV and Indel calling using LongcallD, followed by compression and indexing of the resulting VCF file
