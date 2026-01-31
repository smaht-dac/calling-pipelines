#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: Workflow

requirements:
  MultipleInputFeatureRequirement: {}

inputs:
  - id: input_file_vcf_gz
    type: File
    secondaryFiles:
      - .tbi
    doc: Input VCF (bgzipped) with .tbi index

  - id: germline_input_file_vcf_gz
    type: File
    secondaryFiles:
      - .tbi
    doc: Germline SNV calls input VCF (bgzipped) with .tbi index

  - id: genome_reference_fasta
    type: File
    secondaryFiles:
      - ^.dict
      - .fai
    doc: Reference FASTA with index files

  - id: input_files_sr_cram_tissue_specific
    type:
      -
        items: File
        type: array
    secondaryFiles:
      - .crai
    doc: Short-read CRAM files for current tissue (with .crai)

  - id: input_files_sr_cram_donor_pooled
    type:
      -
        items: File
        type: array
    secondaryFiles:
      - .crai
    doc: Short-read CRAM files for donor (with .crai)

  - id: input_files_pb_cram_donor_pooled
    type:
      -
        items: File
        type: array
    secondaryFiles:
      - .crai
    doc: PacBio CRAM files for donor (with .crai)

  - id: input_files_ont_cram_donor_pooled
    default: []
    type:
      -
        items: File
        type: array
    secondaryFiles:
      - .crai
    doc: ONT CRAM files for donor (with .crai)

  - id: input_files_tissue_descriptors_pb
    type:
      -
        items: File
        type: array
    doc: Tissue identifiers for donor pooled PacBio files (e.g. SMHT009-3A)

  - id: input_files_tissue_descriptors_ont
    default: []
    type:
      -
        items: File
        type: array
    doc: Tissue identifiers for donor pooled ONT files (e.g. SMHT009-3A)

  - id: additional_args
    type: string
    default: "-c -C -Q 20 -q 30 -s 0"
    doc: Additional minipileup args (string)

  - id: group_intervals
    type: int
    default: 100
    doc: Group size for interval batching

  - id: sex
    type: string
    doc: Donor sex (male, female, unknown)

  - id: threads
    type: int
    doc: Number of threads to use

  - id: input_files_bed
    type:
      -
        items: File
        type: array
    doc: List of BED files with regions to exclude from the VCF

  - id: input_files_tissue_descriptors_sr
    type:
      -
        items: File
        type: array
    doc: Tissue identifiers for donor specific short read files (e.g. SMHT009-3A)

  - id: current_tissue 
    type: string
    doc: Tissue identifier for current run (e.g. SMHT009-3A)

outputs:
  output_file_vcf_gz:
    type: File
    outputSource: phase_mosaic_snvs/output_file_vcf_gz

steps:
  minipileup_parallel:
    run: minipileup-parallel.cwl
    in:
      input_file_vcf_gz:
        source: input_file_vcf_gz
      genome_reference_fasta:
        source: genome_reference_fasta
      input_files_sr_cram_tissue_specific:
        source: input_files_sr_cram_tissue_specific
      input_files_pb_cram_donor_pooled:
        source: input_files_pb_cram_donor_pooled 
      input_files_ont_cram_donor_pooled:
        source: input_files_ont_cram_donor_pooled 
      input_files_tissue_descriptors_pb:
        source: input_files_tissue_descriptors_pb
      input_files_tissue_descriptors_ont:
        source: input_files_tissue_descriptors_ont
      additional_args:
        source: additional_args
      group_intervals:
        source: group_intervals
    out:
      - output_file_vcf_gz

  tier_filter_variants_SR_PB_ONT:
    run: tier_filter_variants_SR_PB_ONT.cwl
    in:
      input_file_vcf_gz:
        source: input_file_vcf_gz
      minipileup_vcf_gz:
        source: minipileup_parallel/output_file_vcf_gz
      current_tissue:
        source: current_tissue
    out:
      - output_file_vcf_gz

  phase_mosaic_snvs:
    run: phase_mosaic_snvs.cwl
    in:
      input_file_vcf_gz:
        source: tier_filter_variants_SR_PB_ONT/output_file_vcf_gz
      germline_input_file_vcf_gz:
        source: germline_input_file_vcf_gz
      genome_reference_fasta:
        source: genome_reference_fasta
      input_files_pb_cram_donor_pooled:
        source: input_files_pb_cram_donor_pooled 
      sex:
        source: sex 
      threads:
        source: threads 
    out:
      - output_file_vcf_gz

  bcftools_regions:
    run: bcftools_regions.cwl
    in:
      input_file_vcf_gz:
        source: phase_mosaic_snvs/output_file_vcf_gz
      input_files_bed:
        source: input_files_bed
    out:
      - output_file_vcf_gz

  minipileup_parallel_sr_only:
    run: minipileup-parallel_sr_only.cwl
    in:
      input_file_vcf_gz:
        source: bcftools_regions/output_file_vcf_gz
      genome_reference_fasta:
        source: genome_reference_fasta
      input_files_sr_cram_donor_pooled:
        source: input_files_sr_cram_donor_pooled
      input_files_tissue_descriptors_sr:
        source: input_files_tissue_descriptors_sr
      additional_args:
        source: additional_args
      group_intervals:
        source: group_intervals
    out:
      - output_file_vcf_gz

  parse_CrossTissue_minipileup_result:
    run: parse_CrossTissue_minipileup_result.cwl
    in:
      input_file_vcf_gz:
        source: phase_mosaic_snvs/output_file_vcf_gz
      minipileup_vcf_gz:
        source: minipileup_parallel_sr_only/output_file_vcf_gz
      current_tissue:
        source: current_tissue
    out:
      - output_file_vcf_gz


doc: |
  Filters an SNV VCF file to retain high-confidence variants. |
  Step-2 filters: run minipileup using short-read and long-read CRAM files to compute read support for each SNV, |
  then filter and tier based on read support, strand balance (Fisher) and germline deviation (binomial), |
  then phase/filter based on nearest germline variant |
  then check all other SR within the donor for CrossTissue variants |
  then output a final vcf with HighConf, LowConf, and . filters annotated
