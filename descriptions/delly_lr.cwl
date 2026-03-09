#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: ACCOUNT/delly:VERSION

baseCommand: [run_delly_lr.sh]

inputs:
  - id: output_file_prefix
    type: string
    inputBinding:
      prefix: -n
    doc: Output file name prefix

  - id: input_file_long_cram
    type: 
      - 
        type: File
    secondaryFiles:
      - .crai
    inputBinding:
      prefix: -l
    doc: Long-read input CRAM file (with .crai)

  - id: genome_reference_fasta
    type: File
    secondaryFiles:
      - .fai
      - ^.dict
    inputBinding:
      prefix: -r
    doc: Reference FASTA with index files

outputs:
  - id: output_file_vcf_gz
    type: File
    secondaryFiles:
      - .tbi
    outputBinding:
      glob: $(inputs.output_file_prefix + "-Delly-LR.vcf.gz")

doc: |
  Run Delly on long-read data (PacBio and ONT). |
  Produce structural variant calls in VCF format. 