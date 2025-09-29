#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/calling_utils:VERSION

baseCommand: [bcftools_regions.sh]

inputs:
  - id: input_file_vcf_gz
    type: File
    inputBinding:
      prefix: -i
    secondaryFiles:
      - .tbi
    doc: Input file in VCF format. Compressed with the corresponding index file

  - id: input_files_bed
    type:
      -
        items: File
        type: array
        inputBinding:
          prefix: -x
    doc: List of BED files with regions to exclude from the VCF

  - id: output_prefix
    type: string
    default: output
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
    Run bcftools_regions.sh to filter VCF by excluding multiple BED region lists
