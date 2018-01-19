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
@header patch_88_89_g.sql - Remove probe mappings from the xref tables.
@desc   Remove probe mappings from the xref tables. 
*/

delete 
  xref
from 
  xref, object_xref
where 
  object_xref.xref_id=xref.xref_id
  and object_xref.ensembl_object_type="Probe";

delete from 
  object_xref
where 
  object_xref.ensembl_object_type="Probe";

--  Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_88_89_g.sql|Removed probe mappings from the xref tables');
