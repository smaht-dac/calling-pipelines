#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/calling_utils:VERSION

baseCommand: [granite, filterByTag]

inputs:
  - id: input_file_vcf_gz
    type: File
    inputBinding:
      prefix: -i
    secondaryFiles:
      - .tbi
    doc: Input file in VCF format. Compressed with the corresponding index file

  - id: tag_filters
    type: string[]
    inputBinding:
      prefix: -t
    doc: One or more tag filters. |
         Format "name/value/operator/type/logic[/field=sep][/entry=sep]"

  - id: output_file_name
    type: string
    default: output.vcf
    inputBinding:
      prefix: -o
    doc: Output file in VCF format

  - id: filters_logic
    type: string
    default: null
    inputBinding:
      prefix: -l
    doc: Across-tag logic (combine multiple tag filters). |
         Accept "any" or "all" [any]

  - id: info_separator
    type: string
    default: null
    inputBinding:
      prefix: --separator
    doc: Tag separator within INFO field [;]

outputs:
  - id: output_file_vcf
    type: File
    outputBinding:
      glob: $(inputs.output_file_name)

doc: |
    Run granite filterByTag to filter variants by tag in input VCF file
