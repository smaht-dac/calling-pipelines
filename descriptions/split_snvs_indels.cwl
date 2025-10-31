#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/calling_utils:VERSION

baseCommand: [split_snvs_indels.sh]

inputs:
  - id: input_file_vcf_gz
    type: File
    inputBinding:
      prefix: -i
    secondaryFiles:
      - .tbi
    doc: Input file in VCF format. Compressed with the corresponding index file

  - id: output_prefix
    type: string
    default: "output"
    inputBinding:
      prefix: -o
    doc: Output file prefix

outputs:
  - id: output_snvs_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + "_snvs.vcf.gz")
    secondaryFiles:
      - .tbi
    doc: Compressed VCF file with index. SNVs

  - id: output_indels_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + "_indels.vcf.gz")
    secondaryFiles:
      - .tbi
    doc: Compressed VCF file with index. Indels
