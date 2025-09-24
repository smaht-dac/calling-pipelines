#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/calling_utils:VERSION

baseCommand: [compress_index_vcf.sh]

inputs:
  - id: input_file_vcf
    type: File
    inputBinding:
      position: 1
    doc: Input file in VCF format

  - id: output_prefix
    type: string
    default: output
    inputBinding:
      position: 2
    doc: Output file prefix

outputs:
  - id: output_file_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + ".vcf.gz")
    secondaryFiles:
      - .tbi
    doc: Compressed VCF file with index
