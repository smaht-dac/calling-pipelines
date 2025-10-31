#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/calling_utils:VERSION

baseCommand: [filter_panel_errors.sh]

inputs:
  - id: input_file_vcf_gz
    type: File
    inputBinding:
      prefix: -i
    secondaryFiles:
      - .tbi
    doc: Input file in VCF format. Compressed with the corresponding index file

  - id: panel_of_normals_fasta
    type: File
    inputBinding:
      prefix: -f
    secondaryFiles:
      - .fai
    doc: Panel of normals in FASTA format with the corresponding index file

  - id: output_prefix
    type: string
    default: "pon"
    inputBinding:
      prefix: -o
    doc: Output file prefix

outputs:
  - id: output_file_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + ".vcf.gz")
    secondaryFiles:
      - .tbi
    doc: Compressed VCF file with index

doc: |
    Run filter_panel_errors.sh to filter variants that are present in a panel of normals from BSMN
