-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
-- Copyright [2016-2018] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

# patch_55_56_c.sql
#
# title: analysis_description.display_label_unique_key
#
# description:
# Change the analysis_id key to be unique and alter display_label to be not null
# Also update the ProbeAlign analyses 

#ALTER TABLE analysis_description MODIFY `display_label` varchar(255) NOT NULL;

DROP TABLE IF EXISTS tmp_analysis_description; 

CREATE TABLE `tmp_analysis_description` (
  `analysis_id` int(10) unsigned NOT NULL,
  `description` text,
  `display_label` varchar(255) NOT NULL,
  `displayable` tinyint(1) NOT NULL DEFAULT '1',
  `web_data` text,
  UNIQUE KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;


insert ignore into tmp_analysis_description select * from analysis_description;

DROP TABLE analysis_description;

CREATE TABLE `analysis_description` (
  `analysis_id` int(10) unsigned NOT NULL,
  `description` text,
  `display_label` varchar(255) NOT NULL,
  `displayable` tinyint(1) NOT NULL DEFAULT '1',
  `web_data` text,
  UNIQUE KEY `analysis_idx` (`analysis_id`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;

insert into analysis_description select * from tmp_analysis_description;


DROP TABLE tmp_analysis_description;

###

update analysis_description ad, analysis a set ad.description='Microarray probes from manufacturers are aligned to the genome by Ensembl, if the probe sequences are provided. The mapping is a two-step procedure outlined <a href="/info/docs/microarray_probe_set_mapping.html">here</a>.', ad.display_label='Probe2Transcript Annotation' where a.logic_name='probe2transcript' and a.analysis_id=ad.analysis_id;


insert ignore into analysis_description(analysis_id) select analysis_id from analysis where logic_name='AFFY_UTR_ProbeAlign';
update analysis_description ad, analysis a set ad.description='Genomic alignments for AFFY_UTR arrays', display_label='AFFY_UTR_ProbeAlign', displayable=1, web_data="{'type' => '_oligo', 'key' => 'array_chip', 'colourset' => 'feature', 'display' =>'off' }" where a.logic_name='AFFY_UTR_ProbeAlign' and a.analysis_id=ad.analysis_id;

insert ignore into analysis_description(analysis_id) select analysis_id from analysis where logic_name='AFFY_UTR_ProbeTranscriptAlign';
update analysis_description ad, analysis a set ad.description='Transcript alignments for AFFY_UTR arrays', display_label='AFFY_UTR_ProbeTranscriptAlign', displayable=1, web_data="{'type' => '_oligo', 'key' => 'array_chip', 'colourset' => 'feature', 'display' =>'off' }" where a.logic_name='AFFY_UTR_ProbeTranscriptAlign' and a.analysis_id=ad.analysis_id;


insert ignore into analysis_description(analysis_id) select analysis_id from analysis where logic_name='AFFY_ST_ProbeAlign';
update analysis_description ad, analysis a set ad.description='Genomic alignments for AFFY_ST arrays', display_label='AFFY_ST_ProbeAlign', displayable=1, web_data="{'type' => '_oligo', 'key' => 'array_chip', 'colourset' => 'feature', 'display' =>'off' }" where a.logic_name='AFFY_ST_ProbeAlign' and a.analysis_id=ad.analysis_id;

insert ignore into analysis_description(analysis_id) select analysis_id from analysis where logic_name='AFFY_ST_ProbeTranscriptAlign';
update analysis_description ad, analysis a set ad.description='Transcript alignments for AFFY_ST arrays', display_label='AFFY_ST_ProbeTranscriptAlign', displayable=1, web_data="{'type' => '_oligo', 'key' => 'array_chip', 'colourset' => 'feature', 'display' =>'off' }" where a.logic_name='AFFY_ST_ProbeTranscriptAlign' and a.analysis_id=ad.analysis_id;




insert ignore into analysis_description(analysis_id) select analysis_id from analysis where logic_name='CODELINK_ProbeAlign';
update analysis_description ad, analysis a set ad.description='Genomic alignments for CODELINK arrays', display_label='CODELINK_ProbeAlign', displayable=1, web_data="{'type' => '_oligo', 'key' => 'array_chip', 'colourset' => 'feature', 'display' =>'off' }" where a.logic_name='CODELINK_ProbeAlign' and a.analysis_id=ad.analysis_id;

insert ignore into analysis_description(analysis_id) select analysis_id from analysis where logic_name='AFFY_UTR_ProbeTranscriptAlign';
update analysis_description ad, analysis a set ad.description='Transcript alignments for CODELINK arrays', display_label='CODELINK_ProbeTranscriptAlign', displayable=1, web_data="{'type' => '_oligo', 'key' => 'array_chip', 'colourset' => 'feature', 'display' =>'off' }" where a.logic_name='CODELINK_ProbeTranscriptAlign' and a.analysis_id=ad.analysis_id;



insert ignore into analysis_description(analysis_id) select analysis_id from analysis where logic_name='AGILENT_ProbeAlign';
update analysis_description ad, analysis a set ad.description='Genomic alignments for AGILENT arrays', display_label='AGILENT_ProbeAlign', displayable=1, web_data="{'type' => '_oligo', 'key' => 'array_chip', 'colourset' => 'feature', 'display' =>'off' }" where a.logic_name='AGILENT_ProbeAlign' and a.analysis_id=ad.analysis_id;

insert ignore into analysis_description(analysis_id) select analysis_id from analysis where logic_name='AGILENT_ProbeTranscriptAlign';
update analysis_description ad, analysis a set ad.description='Transcript alignments for AGILENT arrays', display_label='AGILENT_ProbeTranscriptAlign', displayable=1, web_data="{'type' => '_oligo', 'key' => 'array_chip', 'colourset' => 'feature', 'display' =>'off' }" where a.logic_name='AGILENT_ProbeTranscriptAlign' and a.analysis_id=ad.analysis_id;


insert ignore into analysis_description(analysis_id) select analysis_id from analysis where logic_name='ILLUMINA_WG_ProbeAlign';
update analysis_description ad, analysis a set ad.description='Genomic alignments for ILLUMINA_WG arrays', display_label='ILLUMINA_WG_ProbeAlign', displayable=1, web_data="{'type' => '_oligo', 'key' => 'array_chip', 'colourset' => 'feature', 'display' =>'off' }" where a.logic_name='ILLUMINA_WG_ProbeAlign' and a.analysis_id=ad.analysis_id;

insert ignore into analysis_description(analysis_id) select analysis_id from analysis where logic_name='ILLUMINA_WG_ProbeTranscriptAlign';
update analysis_description ad, analysis a set ad.description='Transcript alignments for ILLUMINA_WG arrays', display_label='ILLUMINA_WG_ProbeTranscriptAlign', displayable=1, web_data="{'type' => '_oligo', 'key' => 'array_chip', 'colourset' => 'feature', 'display' =>'off' }" where a.logic_name='ILLUMINA_WG_ProbeTranscriptAlign' and a.analysis_id=ad.analysis_id;


insert ignore into analysis_description(analysis_id) select analysis_id from analysis where logic_name='PHALANX_ProbeAlign';
update analysis_description ad, analysis a set ad.description='Genomic alignments for PHALANX arrays', display_label='PHALANX_ProbeAlign', displayable=1, web_data="{'type' => '_oligo', 'key' => 'array_chip', 'colourset' => 'feature', 'display' =>'off' }" where a.logic_name='PHALANX_ProbeAlign' and a.analysis_id=ad.analysis_id;

insert ignore into analysis_description(analysis_id) select analysis_id from analysis where logic_name='PHALANX_ProbeTranscriptAlign';
update analysis_description ad, analysis a set ad.description='Transcript alignments for PHALANX arrays', display_label='PHALANX_ProbeTranscriptAlign', displayable=1, web_data="{'type' => '_oligo', 'key' => 'array_chip', 'colourset' => 'feature', 'display' =>'off' }" where a.logic_name='PHALANX_ProbeTranscriptAlign' and a.analysis_id=ad.analysis_id;








INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_55_56_c.sql|analysis_description.display_label_not_null');

