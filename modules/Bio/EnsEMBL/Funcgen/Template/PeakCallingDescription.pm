package Bio::EnsEMBL::Funcgen::Template::PeakCallingDescription;

use strict;

use base qw( Exporter );
use vars qw( @EXPORT_OK );

our @EXPORT_OK = qw(
  PEAK_CALLING_TXT_TEMPLATE
);


sub PEAK_CALLING_TXT_TEMPLATE {

  my $description_template = <<TEMPLATE

Peak calling on [% peak_calling.display_label %]
[%- FILTER indent('    ') %] 

Summary
=======

[%- FILTER indent('    ') %] 
- The Epigenome:   [%- PROCESS summarise_epigenome epigenome = peak_calling.fetch_Epigenome %]
- was assayed for: [%- PROCESS summarise_feature_type feature_type = peak_calling.fetch_FeatureType %]
- which creates:   
[%- IF peak_calling.fetch_FeatureType.creates_broad_peaks -%]
Broad peaks
[%- ELSE -%]
Narrow peaks
[%- END %].
- The peak caller: [% peak_calling.fetch_Analysis.display_label %] was used
- and yielded:     [% format_number(peak_calling.num_peaks) %] peaks.
[%- END %] 

Fraction of reads in peaks
==========================

[%- FILTER indent('    ') %] 
- Fraction of reads in peaks: [%- round_percent(100 * peak_calling.fetch_Frip.frip) %]
- Total number of reads: [%- format_number(peak_calling.fetch_Frip.total_reads) %]
[%- END %] 


Idr
===

[%- FILTER indent('    ') %] 
  [%- PROCESS summarise_idr idr = peak_calling.fetch_Idr %]
[%- END %] 

Signal alignment
================

[%- FILTER indent('    ') %]
  [%- signal_alignment = peak_calling.fetch_signal_Alignment -%]
  
  Description of merged alignment
  -------------------------------

  [%- FILTER indent('    ') %]
    [%- PROCESS summarise_alignment alignment = signal_alignment %]
      Phantom peak scores
      -------------------
      [%- FILTER indent('    ') %]
        [% PROCESS summarise_phantom_peak phantom_peak = alignment.fetch_PhantomPeak %]
      [%- END %]

      Chance
      ------
      [%- FILTER indent('    ') %] 
        [%- PROCESS summarise_chance chance = signal_alignment.fetch_Chance_by_control_Alignment(peak_calling.fetch_control_Alignment) %]
      [%- END %]
  [%- END %]
  
  [%- replicate_alignments = signal_alignment.fetch_all_deduplicated_replicate_Alignments %]
  [% IF replicate_alignments.length>0 %]
  Description of replicate alignments
  -----------------------------------

  [%- FILTER indent('    ') %]
    [%- FOR replicate_alignment IN replicate_alignments %]
      [% PROCESS summarise_alignment alignment = replicate_alignment -%]
      
      Phantom peak scores
      -------------------
      [%- FILTER indent('    ') %]
        [% PROCESS summarise_phantom_peak phantom_peak = replicate_alignment.fetch_PhantomPeak %]
      [%- END %]

      Chance
      ------
      [%- FILTER indent('    ') %] 
        [%- PROCESS summarise_chance chance = replicate_alignment.fetch_Chance_by_control_Alignment(peak_calling.fetch_control_Alignment) %]
      [%- END %]

    [%- END -%]
  [%- END %]
  [%- END %]
[%- END %]

Control alignment
=================

[%- FILTER indent('    ') %]
  [%  PROCESS summarise_alignment alignment = peak_calling.fetch_control_Alignment -%]
[%- END -%]

[%- END -%]

[% BLOCK summarise_feature_type %]
[%- feature_type.name %] ([% feature_type.description %])
[%- END -%]

[%- BLOCK summarise_epigenome -%]
[% epigenome.gender.ucfirst %] [% epigenome.display_label %]: [% epigenome.description %].
[%- IF epigenome.efo_accession -%]
  [% epigenome.efo_accession %]
[%- ELSE -%]
  No efo accession available.
[%- END %]
[%- END -%]

[% BLOCK summarise_phantom_peak %]
- Estimated fragment length: [% phantom_peak.est_frag_len %]
- NSC: [% phantom_peak.nsc %]
- RSC: [% phantom_peak.rsc %]
- Quality Tag: [% phantom_peak.quality_tag %]
[%- END -%]

[% BLOCK summarise_chance %]
- Percent of genome enriched: [% round_percent(chance.percent_genome_enriched) %]
- z-score: [% default_round(chance.z_score) %]
- Control enrichment stronger than chip at bin: [% format_number(chance.control_enrichment_stronger_than_chip_at_bin) %]
- PCR amplification bias in Input coverage of 1% of genome: [% chance.pcr_amplification_bias_in_Input_coverage_of_1_percent_of_genome %]
- Differential percentage enrichment: [% default_round(chance.differential_percentage_enrichment) %]
- Divergence: [% scientific_notation(chance.divergence) %]
[%- END -%]

[% BLOCK summarise_idr %]
- Experimental configuration: [% idr.fetch_Experiment.summarise_replicate_configurations %]
- Type of IDR: [% idr.type             %]
- Max peaks threshold: [% format_number(idr.max_peaks)        %]
- Failures:            [% idr.failed_idr_pairs %]
[%- END -%]

[%- BLOCK summarise_alignment %]
  - Name: [% alignment.name %]
[% IF alignment.has_bam_DataFile -%]
  - Alignment in in bam format:   [% canonpath(alignment.fetch_bam_DataFile.relative_ftp_site_path) %]
[%- END %]
[% IF alignment.has_bigwig_DataFile -%]
  - Signal file in bigwig format: [% canonpath(alignment.fetch_bigwig_DataFile.relative_ftp_site_path) %]
[%- END %]
  - Aligned to gender: [% alignment.to_gender      %]
  - Is control: [% bool_to_yes_no(alignment.is_control)     %]
  - Uses all reads from experiment: [% bool_to_yes_no(alignment.is_complete)    %]
  - Has duplicates: [% bool_to_yes_no(alignment.has_duplicates) %]
  - Analysis: [% alignment.fetch_Analysis.logic_name %]
  
  [%- IF alignment.has_duplicates -%]
  - Read files used:
    [%- FILTER indent('        ') -%] 
      [%- PROCESS summarise_read_files read_files = alignment.fetch_all_ReadFiles %]
    [%- END %]
  [%- END -%]
  
  [%- IF not(alignment.has_duplicates) -%]

  - Source: Alignment with duplicates:
    [%- FILTER indent('    ') %]
      [%- PROCESS summarise_alignment alignment = alignment.fetch_source_Alignment %]
    [%- END -%]

  [%- END -%]
[%- END -%]

[%- BLOCK summarise_read_files -%]
  [%- FOR read_file IN read_files %]
    - [%- read_file.name %]

        FastQC
        ------
        [%- FILTER indent('        ') -%] 
          [% PROCESS summarise_fastqc fastqc = read_file.fetch_FastQC %]
        [%- END %]

  [%- END -%]
[%- END -%]

[% BLOCK summarise_fastqc %]
  [%- IF not(fastqc) %]
    - FastQC has not been run on this read file.
  [% END -%]
  [% IF fastqc %]
    - Basic statistics:              [% fastqc.basic_statistics             %]
    - Per base sequence quality:     [% fastqc.per_base_sequence_quality    %]
    - Per tile sequence quality:     [% fastqc.per_tile_sequence_quality    %]
    - Per sequence quality scores:   [% fastqc.per_sequence_quality_scores  %]
    - Per base sequence content:     [% fastqc.per_base_sequence_content    %]
    - Per sequence gc content:       [% fastqc.per_sequence_gc_content      %]
    - Per base n content:            [% fastqc.per_base_n_content           %]
    - Sequence length distribution:  [% fastqc.sequence_length_distribution %]
    - Sequence duplication levels:   [% fastqc.sequence_duplication_levels  %]
    - Overrepresented sequences:     [% fastqc.overrepresented_sequences    %]
    - Adapter content:               [% fastqc.adapter_content              %]
    - Kmer content:                  [% fastqc.kmer_content                 %]
  [%- END -%]
[%- END -%]

TEMPLATE
;
  return $description_template;
}


