#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: ACCOUNT/manta:VERSION

baseCommand: [run_manta.sh]

inputs:
  - id: input_file_cram
    type: File
    inputBinding:
      position: 1
      prefix: -t
    secondaryFiles:
      - .crai
    doc: Input CRAM file with index (.crai). Treated as tumor sample in Manta.

  - id: genome_reference_fasta
    type: File
    inputBinding:
      position: 2
      prefix: -r
    secondaryFiles:
      - .fai
      - ^.dict
    doc: Reference FASTA with index files

  - id: output_prefix
    type: string
    default: "output"
    inputBinding:
      position: 3
      prefix: -o
    doc: Output prefix -> <prefix>.tumorSV.vcf.gz / <prefix>.candidateSV.vcf.gz

  - id: nthreads
    type: int
    default: null
    inputBinding:
      position: 4
      prefix: -j
    doc: Parallel jobs for runWorkflow.py (default nproc inside the script)

outputs:
  - id: output_file_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + ".tumorSV.vcf.gz")
    secondaryFiles:
      - .tbi
    doc: Manta structural variant VCF (scored)

  - id: output_file_candidate_vcf_gz
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + ".candidateSV.vcf.gz")
    secondaryFiles:
      - .tbi
    doc: Manta unscored candidate structural variant VCF file (bgzipped) with index

doc: |
  Wrapper for run_manta.sh to run Manta in tumor-only mode for whole-genome |
  structural variant calling. Use a CRAM file as input
