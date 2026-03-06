#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: ACCOUNT/severus:VERSION

baseCommand: [severus]
arguments: [
    "--min-sv-size", "50"
]

inputs:
  - id: input_file_cram
    type: File
    inputBinding:
      prefix: --target-bam
    secondaryFiles:
      - .crai
    doc: Primary input CRAM (sorted + indexed)

  - id: output_dir
    type: string
    default: "."
    inputBinding:
      prefix: --out-dir
    doc: Output directory

  - id: tandem_repeats_bed
    type: File
    inputBinding:
      prefix: --vntr-bed
    doc: BED file with tandem repeat regions

  - id: nthreads
    type: int
    default: 16
    inputBinding:
      prefix: -t
    doc: Number of threads to use [16]

outputs:
  - id: output_file_vcf
    type: File
    outputBinding:
      glob: $(inputs.output_dir + "/all_SVs/severus_all.vcf")
    doc: VCF is output as ./all_SVs/severus_all.vcf

doc: |
  Severus structural variant calling on long-read sequencing data (PacBio HiFi / ONT).
