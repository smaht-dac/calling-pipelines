#!/usr/bin/env cwl-runner


# severus
# --target-bam $single_Bam
# --out-dir ${outputDir}/${outputFile}
# --vntr-bed $tandem_repeat_bed
# -t 16
# --min-sv-size 50


cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: ACCOUNT/severus:VERSION

baseCommand: [severus]
arguments: [
    "-t", "16",
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

  - id: tandem_repeats_bed
    type: File
    inputBinding:
      prefix: --vntr-bed
    doc: BED file with tandem repeat regions

  - id: nthreads
    type: int
    default: 16
    inputBinding:
      prefix: --threads
    doc: Number of threads to use [16]

  - id: output_file_name
    type: string
    default: "output.vcf.gz"
    inputBinding:
      prefix: --out-dir
    doc: Output VCF file name (.vcf.gz) [output.vcf.gz]

outputs:
  - id: output_file_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_file_name)
    secondaryFiles:
      - .tbi

doc: |
  Severus structural variant calling on long-read sequencing data (PacBio HiFi / ONT).
