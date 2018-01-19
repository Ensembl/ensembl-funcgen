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
@header patch_88_89_f.sql - Remove probe features from object_xref and xref table
@desc   Remove probe features from xref table, they are now in probe_feature_transcript
*/

delete 
  xref
from 
  xref, object_xref
where 
  object_xref.xref_id=xref.xref_id
  and object_xref.ensembl_object_type="ProbeFeature";

delete from 
  object_xref
where 
  object_xref.ensembl_object_type="ProbeFeature";

--  Patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_88_89_f.sql|Removed probe features from object_xref and xref table');
