-- Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

/**
@header patch_91_92_e.sql - Peak_calling table changes
@desc   Peak_calling table changes
*/

alter table peak_calling add CONSTRAINT peak_calling_name_unique UNIQUE (name);
-- alter table peak_calling add    column experiment_id int(15) unsigned default null;
alter table peak_calling change column alignment_id signal_alignment_id int(23) unsigned default NULL;
alter table peak_calling add    column control_alignment_id int(23) unsigned default null;

-- alter table peak_calling drop   column epigenome_id;
-- alter table peak_calling add column epigenome_id int(10) default null;
-- update peak_calling, experiment set peak_calling.epigenome_id=experiment.epigenome_id where peak_calling.experiment_id=experiment.experiment_id;


-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_92_93_e.sql|Peak_calling table changes');