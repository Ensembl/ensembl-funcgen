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

--This patch is to correct a problem with DBs which do not follow
--the standard naming convention(drosophila, c elegans). Non standard names result in new
--coord systems and seq_regions being generated when they are in fact the 
--same as a previous release.  This can cause data loaded on an old coord system
--not to be included in the new coord system, hence data will not be accessable using 
--the latest dnadb.

--The fix is to simple patch the seq_region_ids in the feature tables affected.
--osr = old_seq_region
--nsr = new_seq_region
 
--This is a self join
--will give large product due to nr records!

--select osr.seq_region_id, nsr.seq_region_id, nsr.name from (select seq_region_id, name from seq_region where schema_build='55_54c') osr, (select seq_region_id, name from seq_region where schema_build='56_513a') nsr where osr.name=nsr.name; 

--Could we nest a self join rather than deriving both tables?
--Or is this just hte same?

--select sr.old_seq_region_id, sr.new_seq_region_id, sr.name from (select sr2.seq_region_id as old_seq_region_id, sr1.seq_region_id as new_seq_region_id, sr2.name from seq_region sr2 join seq_region sr1 on sr1.name=sr2.name and sr2.schema_build='55_54c' and sr1.schema_build='56_513a') sr;

--No, this uses a key and one less derived table!


update external_feature ef, (select sr.old_seq_region_id as osr, sr.new_seq_region_id as nsr from (select sr2.seq_region_id as old_seq_region_id, sr1.seq_region_id as new_seq_region_id, sr2.name from seq_region sr2 join seq_region sr1 on sr1.name=sr2.name and sr2.schema_build='55_54c' and sr1.schema_build='56_513a') sr) sr3 set ef.seq_region_id=sr3.nsr where ef.seq_region_id=sr3.osr; 
