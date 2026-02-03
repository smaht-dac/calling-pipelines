#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/calling_utils:VERSION

baseCommand: [remove_germline_variants.sh]

inputs:
  - id: input_file_vcf_gz
    type: File
    inputBinding:
      prefix: -i
    secondaryFiles:
      - .tbi
    doc: Input file in VCF format. Compressed with the corresponding index file

  - id: germline_input_file_vcf_gz
    type: File
    inputBinding:
      prefix: -g
      position: 2
    secondaryFiles:
      - .tbi
    doc: Germline SNV calls input VCF (bgzipped) with .tbi index

  - id: output_prefix
    type: string
    default: "output"
    inputBinding:
      prefix: -o
    doc: Output file prefix

outputs:
  - id: output_file_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + ".nogermline.vcf.gz")
    secondaryFiles:
      - .tbi
    doc: Compressed VCF file with index

doc: |
    Run remove_germline_variants.sh to remove any germline variants from somatic calls
