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


/**
@header patch_75_76_b.sql - result_set/experiment.display_url/archive_id
@desc   Move input_subset display_url and archive_id to the experiment table, and 
        add an experiment_id field to result_set. Drop redundant input_subset fields.
*/


ALTER TABLE experiment ADD COLUMN archive_id varchar(60) DEFAULT NULL; 
-- To handle concated archive_ids
ALTER TABLE experiment ADD COLUMN display_url varchar(255) DEFAULT NULL;
ALTER TABLE result_set ADD COLUMN experiment_id int(10) unsigned default NULL;
ALTER TABLE result_set ADD KEY    experiment_idx (experiment_id);
ALTER TABLE feature_set ADD COLUMN experiment_id int(10) unsigned default NULL;
ALTER TABLE feature_set ADD KEY    experiment_idx (experiment_id);

/** Data patch for release 75 data **/ 
 
UPDATE feature_set fs, input_set inp set fs.experiment_id=inp.experiment_id where fs.input_set_id=inp.input_set_id;

UPDATE result_set rs, supporting_set ss, data_set ds, feature_set fs set rs.experiment_id=fs.experiment_id 
   WHERE fs.experiment_id is not NULL and fs.feature_set_id=ds.feature_set_id and ds.data_set_id=ss.data_set_id 
   AND   ss.supporting_set_id=rs.result_set_id and ss.type='result';
  
UPDATE experiment e, input_set inp, input_set_input_subset isiss, input_subset iss set e.archive_id=iss.archive_id, e.display_url=iss.display_url
   WHERE e.experiment_id=inp.experiment_id and inp.input_set_id=isiss.input_set_id and isiss.input_subset_id=iss.input_subset_id 
   AND iss.name not like "%WCE%" and iss.is_control=0;

/** Update for release 76 data
 * This was loaded using 75 schema and patched up mid-release cycle
 * 
 * UPDATE experiment e, input_subset iss, result_set_input rsi, result_set rs 
 * SET e.display_url=iss.display_url, e.archive_id=iss.archive_id, rs.experiment_id=e.experiment_id
 * WHERE e.experiment_id=iss.experiment_id and iss.is_control=0 and iss.input_subset_id=rsi.table_id
 * AND   rsi.table_name='input_subset' and rsi.result_set_id=rs.result_set_id;
 * 
 *  Just result_set.experiment_id update after the first patch  
 * UPDATE input_subset iss, result_set_input rsi, result_set rs 
 * SET rs.experiment_id=iss.experiment_id
 * WHERE iss.is_control=0 and iss.input_subset_id=rsi.table_id
 * AND   rsi.table_name='input_subset' and rsi.result_set_id=rs.result_set_id;  
 * 
 * UPDATE feature_set fs JOIN data_set ds USING(feature_set_id) 
 *   JOIN supporting_set ss USING(data_set_id)
 * JOIN result_set rs ON ss.supporting_set_id=rs.result_set_id
 * SET fs.experiment_id=rs.experiment_id WHERE ss.type='result';
**/  
  
ALTER TABLE feature_set DROP COLUMN input_set_id;
ALTER TABLE input_subset DROP COLUMN archive_id;
ALTER TABLE input_subset DROP COLUMN display_url;
analyze table result_set;
optimize table result_set;
  
-- patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_75_76_b.sql|result/feature_set.experiment_id & experiment/input_subset.display_url/archive_id');


