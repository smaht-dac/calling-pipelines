#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: ACCOUNT/longcalld:VERSION

baseCommand: [longcallD, call]
arguments: [-s] # mosaic / somatic

inputs:
  - id: sample_name
    type: string
    inputBinding:
      prefix: -n
      position: 1
    doc: Sample name in VCF

  - id: platform
    type: string
    default: "--hifi"
    inputBinding:
      position: 2
    doc: Sequencing platform flag --hifi / --ont [--hifi]

  - id: nthreads
    type: int
    default: 8
    inputBinding:
      prefix: -t
      position: 3
    doc: Number of threads to use [8]

  - id: output_file_name
    type: string
    default: "out.vcf"
    inputBinding:
      prefix: -o
      position: 4
    doc: Output VCF file name [out.vcf]

  - id: mei_reference_fasta
    type: File
    inputBinding:
      prefix: -T
      position: 5
    doc: Reference FASTA for MEI (Mobile Element Insertion)

  # Required positional reference
  - id: genome_reference_fasta
    type: File
    inputBinding:
      position: 6
    secondaryFiles:
      - ^.dict
      - .fai
    doc: Reference FASTA with index files

  # Required primary positional BAM/CRAM
  - id: input_file_cram
    type: File
    inputBinding:
      position: 7
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
        inputBinding:
          prefix: -X
    secondaryFiles:
      - .crai
    inputBinding:
      position: 8
    doc: Additional CRAM files (sorted + indexed)

outputs:
  - id: output_file_vcf
    type: File
    outputBinding:
      glob: $(inputs.output_file_name)

doc: |
  LongcallD somatic (-s) variant caller for long-read sequencing data (PacBio HiFi / ONT).
