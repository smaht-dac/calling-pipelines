#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: ACCOUNT/severus:VERSION

baseCommand: [run_kanpig_plup.sh]

inputs:
  - id: input_file_cram
    type: File
    inputBinding:
      prefix: -i
    secondaryFiles:
      - .crai
    doc: Primary input CRAM (sorted + indexed)

  - id: reference_fasta
    type: File
    inputBinding:
      prefix: -r
    secondaryFiles:
      - .fai
    doc: Reference FASTA (indexed)

  - id: output_prefix
    type: string
    default: "output"
    inputBinding:
      prefix: -o
    doc: Output file prefix

  - id: nthreads
    type: int
    default: 8
    inputBinding:
      prefix: -t
    doc: Number of threads to use [8]

outputs:
  - id: output_file_plup
    type: File
    secondaryFiles:
      - .tbi
    outputBinding:
      glob: $(inputs.output_prefix + ".plup.gz")
    doc: Plup file is output as output.plup.gz

doc: |
  Kanpig structural variant pileups on long-read sequencing data (PacBio HiFi / ONT).
