#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/calling_utils:VERSION

baseCommand: [tier_filter_variants_SR_PB_ONT.py]

inputs:
  - id: input_file_vcf_gz
    type: File
    inputBinding:
      prefix: --input_vcf
    secondaryFiles:
      - .tbi
    doc: Input VCF (bgzipped) with corresponding tabix index (.tbi)

  - id: minipileup_vcf_gz
    type: File
    inputBinding:
      prefix: --minipileup_vcf
    secondaryFiles:
      - .tbi
    doc: Minipileup VCF input (bgzipped) with corresponding tabix index (.tbi)

  - id: current_tissue
    type: string
    inputBinding:
      prefix: --current_tissue
    doc: Tissue identifier for current run (e.g. SMHT009-3A)

  - id: output_file_name
    type: string
    default: "tiered.vcf.gz"
    inputBinding:
      prefix: --output_vcf
    doc: Output VCF filename; use .vcf.gz to enable bgzip + tabix indexing

outputs:
  - id: output_file_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_file_name)
    secondaryFiles:
      - .tbi

doc: |
    Tier and filter variants using short-read (SR), PacBio, and ONT minipileup counts. |
    Applies read-support thresholds, Fisher strand-bias, and binomial germline-deviation tests
