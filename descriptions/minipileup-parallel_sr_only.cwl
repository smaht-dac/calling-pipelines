#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: ACCOUNT/calling_utils:VERSION

baseCommand: [minipileup-parallel_sr_only.sh]

inputs:
  - id: input_file_vcf_gz
    type: File
    inputBinding:
      prefix: -i
      position: 1
    secondaryFiles:
      - .tbi
    doc: Input VCF (bgzipped) with .tbi index, from bcftools_regions (only easy allowed)

  - id: genome_reference_fasta
    type: File
    inputBinding:
      prefix: -r
      position: 2
    secondaryFiles:
      - ^.dict
      - .fai
    doc: Reference FASTA with index files

  - id: input_files_sr_cram_donor_pooled
    type:
      -
        items: File
        type: array
        inputBinding:
          prefix: --sr-cram
    secondaryFiles:
      - .crai
    inputBinding:
      position: 3
    doc: Short-read CRAM files (with .crai)

  - id: input_files_tissue_descriptors_sr
    type:
      -
        items: File
        type: array
        inputBinding:
          prefix: --sr-tissue
    inputBinding:
      position: 4
    doc: Tissue identifiers (e.g. SMHT009-3A)

  - id: output_prefix
    type: string
    default: "minipileup"
    inputBinding:
      prefix: -o
      position: 5
    doc: Output file prefix

  - id: additional_args
    type: string
    default: "-c -C -Q 30 -q 30 -s 0"
    inputBinding:
      prefix: --args
      position: 6
    doc: Additional minipileup args (string)

  - id: group_intervals
    type: int
    default: 100
    inputBinding:
      prefix: --group
      position: 7
    doc: Group size for interval batching

outputs:
  - id: output_file_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + ".vcf.gz")
    secondaryFiles:
      - .tbi
    doc: Compressed VCF with index

doc: |
  Run minipileup in parallel on multiple short-read CRAM files.
  Outputs a compressed VCF with index
