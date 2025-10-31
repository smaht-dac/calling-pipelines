#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/vep:VERSION

baseCommand: [vep-parallel.sh]

inputs:
  - id: input_file_vcf_gz
    type: File
    inputBinding:
      prefix: -i
    secondaryFiles:
      - .tbi
    doc: Input file in VCF format. Compressed with the corresponding index file

  - id: regions_list_txt
    type: File
    inputBinding:
      prefix: -l
    doc: Regions list file (one region per line, e.g. chr1:1-1000000)

  - id: vep_database_archive
    type: File
    inputBinding:
      prefix: -v
    doc: Compressed VEP database archive (tar.gz)

  - id: gnomad_vcf_gz
    type: File
    inputBinding:
      prefix: -g
    secondaryFiles:
      - .tbi
    doc: Compressed gnomAD VCF file with corresponding index file

  - id: genome_reference_fasta
    type: File
    inputBinding:
      prefix: -r
    secondaryFiles:
      - ^.dict
      - .fai
    doc: Genome reference in FASTA format with the corresponding index files

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
      glob: $(inputs.output_prefix + ".vep.vcf.gz")
    secondaryFiles:
      - .tbi
    doc: Compressed VCF file with index

doc: |
    Run VEP in parallel mode using multiple regions.
