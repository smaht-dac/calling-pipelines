#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/rufus:VERSION

baseCommand: [run_rufus.sh]

inputs:
  - id: input_file_cram
    type: File
    inputBinding:
      prefix: -s
    secondaryFiles:
      - .crai
    doc: Input CRAM file (with .crai)

  - id: genome_reference_fasta
    type: File
    inputBinding:
      prefix: -r
    secondaryFiles:
      - ^.dict
      - .fai
    doc: Reference FASTA with index files

  - id: genome_reference_bwt
    type: File
    inputBinding:
      prefix: -b
      valueFrom: $(self.path.match(/(.*)\.[^.]+$/)[1])
    secondaryFiles:
      - ^.ann
      - ^.amb
      - ^.pac
      - ^.sa
    doc: Genome reference in BWT format with the corresponding index files

  - id: region_file_tsv
    type: File
    inputBinding:
      prefix: -f
    doc: TSV File with regions grouped by shard |
         shard_index <TAB> region, e.g. chr1:1-1000000

  - id: shard_index
    type: int
    inputBinding:
      prefix: -i
    doc: Index to use to extract the right set of regions for the shard |
         from inputs.region_file_tsv

  - id: rufus_directory
    type: Directory
    default:
      class: Directory
      path: "/data1/input-mounted-smaht-wolf-application-files/RUFUS_HASH"
      # /data1/input-mounted-smaht-production-application-files/RUFUS_HASH
    inputBinding:
      prefix: -d
    doc: Path to main bucket containing Jhash files. |
         See run_rufus.sh help for expected directory structure

  - id: nthreads
    type: int
    default: 8
    inputBinding:
      prefix: -t
    doc: Number of threads to use per region

  - id: kg1_version
    type: string
    default: v3.0
    inputBinding:
      prefix: -g
    doc: KG1 hash version

  - id: control_version
    type: string
    default: v1.0
    inputBinding:
      prefix: -c
    doc: Control hash version

  - id: output_prefix
    type: string
    default: output
    inputBinding:
      prefix: -o
    doc: <OUTPUT_PREFIX>.vcf.gz

outputs:
  - id: output_file_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + ".vcf.gz")
    secondaryFiles:
      - .tbi

doc: |
  Run RUFUS per shard. Each shard can run multiple regions using parallel jobs
