#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: ACCOUNT/calling_utils:VERSION

baseCommand: [minipileup-parallel_v2.sh]

inputs:
  - id: input_file_vcf_gz
    type: File
    inputBinding:
      prefix: -i
      position: 1
    secondaryFiles:
      - .tbi
    doc: Input VCF (bgzipped) with .tbi index

  - id: genome_reference_fasta
    type: File
    inputBinding:
      prefix: -r
      position: 2
    secondaryFiles:
      - ^.dict
      - .fai
    doc: Reference FASTA with index files

  - id: input_files_sr_cram_tissue_specific 
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
    doc: Short-read CRAM files for tissue (with .crai)
  
  - id: input_files_core_ids_sr
    type:
      -
        items: string
        type: array
        inputBinding:
          prefix: --sr-core
    inputBinding:
      position: 4
    doc: Core identifiers for short-read tissue-specific CRAM files (1:1 match) 

  - id: input_files_all_long_read_cram_donor_pooled
    type:
      -
        items: File
        type: array
        inputBinding:
          prefix: --lr-cram
    secondaryFiles:
      - .crai
    inputBinding:
      position: 5
    doc: PacBio + ONT CRAM files for donor (with .crai)

  - id: input_files_tissue_descriptors_all_long_read
    type:
      -
        items: string
        type: array
        inputBinding:
          prefix: --lr-tissue
    inputBinding:
      position: 6
    doc: Tissue identifiers for PacBio + ONT (1:1 match) (e.g. SMHT009-3A)

  - id: input_files_types_all_long_read
    type:
      -
        items: string
        type: array
        inputBinding:
          prefix: --lr-type
    inputBinding:
      position: 7
    doc: Sequencing type identifiers for PacBio + ONT (1:1 match) (e.g. PB ONT...)
  
  - id: input_files_core_ids_all_long_read
    type:
      -
        items: string
        type: array
        inputBinding:
          prefix: --lr-core
    inputBinding:
      position: 8
    doc: Core identifiers for PacBio + ONT CRAM files (1:1 match) 

  - id: output_prefix
    type: string
    default: "minipileup2"
    inputBinding:
      prefix: -o
      position: 9
    doc: Output file prefix

  - id: additional_args
    type: string
    default: "-D -c -C -Q 20 -q 30 -s 0"
    inputBinding:
      prefix: --args
      position: 10
    doc: Additional minipileup args (string)

outputs:
  - id: output_file_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + ".vcf.gz")
    secondaryFiles:
      - .tbi
    doc: Compressed VCF with index

doc: |
  Run minipileup2 on multiple CRAM files (short-read, PacBio, ONT). |
  De-duplicate counts for read fragments that overlap at variant position. |
  Outputs a compressed VCF with index
