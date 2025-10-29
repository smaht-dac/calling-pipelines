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
    doc: Input file in VCF format. Compressed with the corresponding index file

  - id: genome_reference_fasta
    type: File
    secondaryFiles:
      - ^.dict
      - .fai
    doc: Genome reference in FASTA format with the corresponding index files

  - id: input_files_bed
    type:
      -
        items: File
        type: array
    doc: List of BED files with regions to exclude from the VCF

  - id: panel_of_normals_fasta
    type: File
    secondaryFiles:
      - .fai
    doc: Panel of normals in FASTA format with the corresponding index file

  - id: regions_list_txt
    type: File
    doc: Regions list file (one region per line, e.g. chr1:1-1000000)

  - id: vep_database_archive
    type: File
    doc: Compressed VEP database archive (tar.gz)

  - id: gnomad_vcf_gz
    type: File
    secondaryFiles:
      - .tbi
    doc: Compressed gnomAD VCF file with corresponding index file

  - id: tag_filters
    type: string[]
    doc: One or more tag filters. |
         Format "name/value/operator/type/logic[/field=sep][/entry=sep]"

  - id: window_size
    type: int
    default: 100
    doc: Window size in bp [100]

outputs:
  output_snvs_vcf_gz:
    type: File
    outputSource: split_snvs_indels/output_snvs_vcf_gz

  output_indels_vcf_gz:
    type: File
    outputSource: split_snvs_indels/output_indels_vcf_gz

steps:
  bcftools_PASS_norm_dedup:
    run: bcftools_PASS_norm_dedup.cwl
    in:
      input_file_vcf_gz:
        source: input_file_vcf_gz
      genome_reference_fasta:
        source: genome_reference_fasta
    out:
      - output_file_vcf_gz

  bcftools_regions:
    run: bcftools_regions.cwl
    in:
      input_file_vcf_gz:
        source: bcftools_PASS_norm_dedup/output_file_vcf_gz
      input_files_bed:
        source: input_files_bed
    out:
      - output_file_vcf_gz

  filter_panel_errors:
    run: filter_panel_errors.cwl
    in:
      input_file_vcf_gz:
        source: bcftools_regions/output_file_vcf_gz
      panel_of_normals_fasta:
        source: panel_of_normals_fasta
    out:
      - output_file_vcf_gz

  filter_clustered_variants:
    run: filter_clustered_variants.cwl
    in:
      input_file_vcf_gz:
        source: filter_panel_errors/output_file_vcf_gz
      window_size:
        source: window_size
    out:
      - output_file_vcf_gz

  vep_parallel:
    run: vep-parallel.cwl
    in:
      input_file_vcf_gz:
        source: filter_clustered_variants/output_file_vcf_gz
      regions_list_txt:
        source: regions_list_txt
      vep_database_archive:
        source: vep_database_archive
      gnomad_vcf_gz:
        source: gnomad_vcf_gz
      genome_reference_fasta:
        source: genome_reference_fasta
    out:
      - output_file_vcf_gz

  granite_filterByTag:
    run: granite_filterByTag.cwl
    in:
      input_file_vcf_gz:
        source: vep_parallel/output_file_vcf_gz
      tag_filters:
        source: tag_filters
    out:
      - output_file_vcf

  compress_index_vcf:
    run: compress_index_vcf.cwl
    in:
      input_file_vcf:
        source: granite_filterByTag/output_file_vcf
    out:
      -  output_file_vcf_gz

  split_snvs_indels:
    run: split_snvs_indels.cwl
    in:
      input_file_vcf_gz:
        source: compress_index_vcf/output_file_vcf_gz
    out:
      - output_snvs_vcf_gz
      - output_indels_vcf_gz

doc: |
  Filters raw SNVs and Indels VCF files to retain high-confidence variants. |
  Step-1 filters: PASS calls; normalize and atomize variants; remove genomic regions; |
  Brain Somatic Mosaicism Network (BSMN) filter; remove clustered variants; run VEP and filter by allele frequency. |
  Splits the filtered variants into separate SNVs and Indels VCF files
