#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: ACCOUNT/sniffles:VERSION

baseCommand: [sniffles]
arguments: [
    "--allow-overwrite", 
    "--minsvlen", "50",
    "--mosaic",
    "--mosaic-include-germline",
    "--mosaic-af-min", "0.001",
    "--cluster-binsize", "50",
    "--cluster-merge-len", "0.10",
    "--minsupport", "1",
    "--no-qc"
]

inputs:
  - id: input_file_cram
    type: File
    inputBinding:
      prefix: --input
    secondaryFiles:
      - .crai
    doc: Primary input CRAM (sorted + indexed)

  - id: genome_reference_fasta
    type: File
    inputBinding:
      prefix: --reference
    secondaryFiles:
      - ^.dict
      - .fai
    doc: Reference FASTA with index files

  - id: tandem_repeats_bed
    type: File
    inputBinding:
      prefix: --tandem-repeats
    doc: BED file with tandem repeat regions

  - id: nthreads
    type: int
    default: 32
    inputBinding:
      prefix: --threads
    doc: Number of threads to use [32]

  - id: output_file_name
    type: string
    default: "output.vcf.gz"
    inputBinding:
      prefix: --vcf
    doc: Output VCF file name (.vcf.gz) [output.vcf.gz]

outputs:
  - id: output_file_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_file_name)
    secondaryFiles:
      - .tbi

doc: |
  Sniffles structural variant calling on long-read sequencing data (PacBio HiFi / ONT).
