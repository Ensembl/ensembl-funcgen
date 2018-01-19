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
@header patch_74_75_d.sql - result_set_input.table_name_input_subset
@desc   input_subset was actually added to the enum in patch_73_74_c
        This patch is simply converting the result_set_input records 
        from input_set to input_subset
*/

INSERT INTO 
  result_set_input(
    result_set_id, 
    table_name, 
    table_id
  ) 
  SELECT  
    rsi.result_set_id, 
    'input_subset', 
    isiss.input_subset_id 
  FROM 
    input_set              inp, 
    input_set_input_subset isiss, 
    result_set_input       rsi 
  WHERE 
    rsi.table_id     = inp.input_set_id AND 
    rsi.table_name   ='input_set'       AND 
    inp.input_set_id = isiss.input_set_id
    ;

DELETE FROM 
  result_set_input 
WHERE 
  table_name='input_set';

 
analyze table  result_set_input;
optimize table result_set_input;
 
# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_74_75_d.sql|result_set_input.table_name_input_subset');


