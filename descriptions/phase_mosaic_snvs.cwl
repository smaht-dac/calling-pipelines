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

  - id: input_files_pb_cram
    type:
      -
        items: File
        type: array
        inputBinding:
          prefix: --pb-cram
    secondaryFiles:
      - .crai
    inputBinding:
      position: 3
    doc: PacBio CRAM files (with .crai)

  - id: sample_id
    type: string
    default: "sample"
    inputBinding:
      prefix: -i
      position: 4
    doc: Output file prefix

  - id: sex
    type: string
    default: "unknown"
    inputBinding:
      prefix: -s
      position: 5
    doc: Donor sex (male, female, unknown)

  - id: threads
    type: int
    default: 1 
    inputBinding:
      prefix: -t
      position: 6
    doc: Number of threads to use

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
