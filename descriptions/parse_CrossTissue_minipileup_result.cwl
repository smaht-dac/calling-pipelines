#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/calling_utils:VERSION

baseCommand: [parse_CrossTissue_minipileup_result.sh]

inputs:
  - id: input_file_vcf_gz
    type: File
    inputBinding:
      prefix: -i
    secondaryFiles:
      - .tbi
    doc: Input VCF (bgzipped) with corresponding tabix index (.tbi)

  - id: minipileup_vcf_gz
    type: File
    inputBinding:
      prefix: --mp
    secondaryFiles:
      - .tbi
    doc: Minipileup VCF input (bgzipped) with corresponding tabix index (.tbi)

  - id: current_tissue
    type: string
    inputBinding:
      prefix: --tissue
    doc: Current tissue being run (e.g. 3A or SMHT009-3A)

  - id: output_prefix
    type: string
    default: "output"
    inputBinding:
      prefix: -o
    doc: Output prefix

outputs:
  - id: output_file_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + ".final.vcf.gz")
    secondaryFiles:
      - .tbi

doc: |
    Parse results of SR only minipileup to look at CrossTissue patterns |
    return a final vcf file with HighConf, LowConf and . filter designations
