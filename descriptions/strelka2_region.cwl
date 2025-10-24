#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/strelka2:VERSION

baseCommand: [strelka2_region.sh]

inputs:
  - id: input_normal_cram
    type: File
    inputBinding:
      position: 1
      prefix: -n
    secondaryFiles:
      - .crai
    doc: Input CRAM file with index (.crai). Normal sample

  - id: input_tumor_cram
    type: File
    inputBinding:
      position: 2
      prefix: -t
    secondaryFiles:
      - .crai
    doc: Input CRAM file with index (.crai). Tumor sample

  - id: genome_reference_fasta
    type: File
    inputBinding:
      position: 3
      prefix: -f
    secondaryFiles:
      - .fai
      - ^.dict
    doc: Reference FASTA with index files

  - id: regions_bed_gz
    type: File
    inputBinding:
      position: 4
      prefix: -b
    secondaryFiles:
      - .tbi
    doc: Regions BED file, bgzipped (.bed.gz) with .tbi index

  - id: line_index
    type: int
    default: null
    inputBinding:
      position: 5
      prefix: -i
    doc: Optional 1-based line index into regions BED (no header assumed)

  - id: output_prefix
    type: string
    default: output
    inputBinding:
      position: 6
      prefix: -o
    doc: Output prefix -> <prefix>.snvs.vcf.gz / <prefix>.indels.vcf.gz

  - id: threads
    type: int
    default: null
    inputBinding:
      position: 7
      prefix: -j
    doc: Parallel jobs for runWorkflow.py (default nproc inside the script)

outputs:
  - id: output_snvs_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + ".snvs.vcf.gz")
    secondaryFiles:
      - .tbi
    doc: Strelka2 somatic SNVs VCF file (bgzipped) with index

  - id: output_indels_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + ".indels.vcf.gz")
    secondaryFiles:
      - .tbi
    doc: Strelka2 somatic indels VCF file (bgzipped) with index

doc: |
  Wrapper for strelka2_region.sh to run Strelka2 paired somatic variant calling on specific regions |
  from a BED file (optionally a specific line in the BED file). Use CRAM files as input
