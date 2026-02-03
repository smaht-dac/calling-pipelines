#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: ACCOUNT/calling_utils:VERSION

baseCommand: [phase_mosaic_vars.sh]

inputs:
  - id: input_file_vcf_gz
    type: File
    inputBinding:
      prefix: -w
      position: 1
    secondaryFiles:
      - .tbi
    doc: Input VCF (bgzipped) with .tbi index

  - id: germline_input_file_vcf_gz
    type: File
    inputBinding:
      prefix: -v
      position: 2
    secondaryFiles:
      - .tbi
    doc: Germline SNV calls input VCF (bgzipped) with .tbi index

  - id: genome_reference_fasta
    type: File
    inputBinding:
      prefix: -r
      position: 3
    secondaryFiles:
      - ^.dict
      - .fai
    doc: Reference FASTA with index files

  - id: input_files_pb_cram_donor_pooled
    type:
      -
        items: File
        type: array
        inputBinding:
          prefix: --pb-cram
    secondaryFiles:
      - .crai
    inputBinding:
      position: 4
    doc: PacBio CRAM files (with .crai)

  - id: sample_id
    type: string
    default: "sample"
    inputBinding:
      prefix: -i
      position: 5
    doc: Output file prefix

  - id: sex
    type: string
    default: "unknown"
    inputBinding:
      prefix: -s
      position: 6
    doc: Donor sex (male, female, unknown)

outputs:
  - id: output_file_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.sample_id + ".phased.vcf.gz")
    secondaryFiles:
      - .tbi
    doc: Compressed VCF with index

doc: |
  Run phasing on tiered variants using PacBio CRAM files. |
  Outputs a compressed VCF with index
