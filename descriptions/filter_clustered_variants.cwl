#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/calling_utils:VERSION

baseCommand: [filter_clustered_variants.py]

inputs:
  - id: input_file_vcf_gz
    type: File
    inputBinding:
      position: 2
    secondaryFiles:
      - .tbi
    doc: Input file in VCF format. Compressed with the corresponding index file

  - id: output_file_name
    type: string
    default: "output.vcf.gz"
    inputBinding:
      position: 3
    doc: Output file in VCF format. Compressed with bgzip (use .vcf.gz suffix)

  - id: window_size
    type: int
    default: 100
    inputBinding:
      position: 1
      prefix: --window
    doc: Window size in bp [100]

outputs:
  - id: output_file_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_file_name)
    secondaryFiles:
      - .tbi

doc: |
    This tool filters clustered variants from a VCF file. |
    Variants are considered clustered if there is another variant within a specified window size (default 100 bp)