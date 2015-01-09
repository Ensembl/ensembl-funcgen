-- Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

-- pre reg build QC on supporting feature_sets
-- requires regbuild.feature_set_ids to be update, so we need to archive and setup the new set before this is run.
-- basically we need to check we have an old version of this value for the correct schema_build/version?
--This now counts all features from regulatory feature supporting sets over all available seq_region for a given fset
-- Takes about 30-40 mins
--Need to restrict this to a particular coord_system
select fs.name as 'FeatureSet', fs.feature_set_id, sr.name as 'Chr', count(*) as 'FeatureCount' from feature_set fs, annotated_feature af left join (select distinct(seq_region_id), name from seq_region where schema_build='55_37') sr on sr.seq_region_id=af.seq_region_id where af.feature_set_id=fs.feature_set_id and fs.feature_set_id in( select a.feature_set_id from (select substring_index(substring_index(m.meta_value, ',',@r:=@r+1),',',-1) feature_set_id from (select @r:=0) x, (select ss.supporting_set_id from supporting_set ss, data_set ds where ds.name='RegulatoryFeatures' and ds.data_set_id=ss.data_set_id) z, meta m where m.meta_key='regbuild.feature_set_ids') a) group by fs.name, sr.name order by sr.name, fs.name;

-- We don't get counts of 0 for those sr's which are not present, we can resolve this if we wrap it in a perl script


-- compare counts feature counts across chromosome and release
-- sr product is account for by distinct sid
select count(distinct(rf.regulatory_feature_id)), sr.name, fs.name, count(distinct(rf.regulatory_feature_id)) from regulatory_feature rf, feature_set fs, seq_region sr where fs.name in ('RegulatoryFeatures', 'RegulatoryFeatures_v6') and fs.feature_set_id=rf.feature_set_id and rf.seq_region_id=sr.seq_region_id group by rf.feature_set_id, rf.seq_region_id order by sr.name, rf.feature_set_id;

--This takes a while and can probably be done better
--maybe do the nested left join where select distinct sr's as above


-- Total counts across feature_sets
select fs.name, count(*)  from regulatory_feature rf, feature_set fs where fs.feature_set_id in(69,112) and fs.feature_set_id=rf.feature_set_id group by rf.feature_set_id;


-- post analysis feature type QC


select count(distinct(rf.regulatory_feature_id)), ft.name, fs.name from regulatory_feature rf, feature_set fs, feature_type ft where fs.name in ('RegulatoryFeatures', 'RegulatoryFeatures_v6') and fs.feature_set_id=rf.feature_set_id and rf.feature_type_id=ft.feature_type_id group by rf.feature_set_id, rf.feature_type_id order by rf.feature_type_id, rf.feature_set_id;


-- by chr?

-- This counts all inclusions of regulatory attributes, so there may be duplicates
-- From inclusions in two separate features
-- Is there a way around this? distinct(af.annotated_feature_id)?

 select fs1.name, sr.name, fs1.feature_set_id, count(*) from annotated_feature af, regulatory_attribute ra, feature_set fs, feature_set fs1, regulatory_feature rf left join (select distinct(seq_region_id), name from seq_region) sr on sr.seq_region_id=af.seq_region_id where fs.name='RegulatoryFeatures' and fs.feature_set_id=rf.feature_set_id and rf.regulatory_feature_id=ra.regulatory_feature_id and ra.attribute_feature_id=af.annotated_feature_id and ra.attribute_feature_table='annotated' and af.feature_set_id=fs1.feature_set_id group by fs1.name, sr.name order by sr.name, fs1.name;


-- now we want stats(mean, median, mode) on length of rfs and core regions
-- genome wide or by chr? Will with rollup work here?
-- mean, no problem
-- median would require tmp tables
-- mode is not really possible without proper ditribution modelling?
-- Just leave out mode for now.

-- Roll up does do correct averages here!
-- Can we get longest feature for core and full for each seq_region?

-- We need to nest the left join in to the sub select!!!

--Can't be temporary table as we wan tot refer to it more than once inthe following query
DROP TABLE IF EXISTS rf_lengths;

CREATE TABLE rf_lengths select sr.name, rf.regulatory_feature_id, rf.feature_set_id, (rf.seq_region_end - rf.seq_region_start + 1) as core_length,  (rf.bound_seq_region_end - rf.bound_seq_region_start + 1) as full_length,  rf.seq_region_start, rf.seq_region_end from regulatory_feature rf left join (select distinct(seq_region_id), name from seq_region) sr on sr.seq_region_id=rf.seq_region_id join feature_set fs ON fs.name='RegulatoryFeatures' where fs.feature_set_id=rf.feature_set_id;

alter table rf_lengths add key `name_core_length` (name, core_length);
alter table rf_lengths add key `name_full_length` (name, full_length);


--No keys on here so this select will take ages!

DROP TABLE IF EXISTS rf_stats;
create table rf_stats SELECT rfl.name, cmin.min as 'min_core_length', AVG(rfl.core_length) as 'avg_core_length', cmax.max as 'max_core_length', fmin.min as 'min_full_length', AVG(rfl.full_length) as 'avg_full_length',  fmax.max as 'max_full_length'
--rfl_cmin.regulatory_feature_id as 'Min. core id', rfl_cmax.regulatory_feature_id as 'Max. core id', rfl_fmin.regulatory_feature_id as 'Min. full length id', rfl_fmax.regulatory_feature_id as 'Max. full length id'
	from rf_lengths rfl, 
	-- rf_lengths rfl_cmin, rf_lengths rfl_cmax, rf_lengths rfl_fmin, rf_lengths rfl_fmax, 
	--fetching id with MIN etc does not work
	(select name, MIN(core_length) as 'min' from rf_lengths group by name) as cmin, 
	(select name, MAX(core_length) as 'max' from rf_lengths group by name) as cmax,
	(select name, MIN(full_length) as 'min' from rf_lengths group by name) as fmin, 
	(select name, MAX(full_length) as 'max' from rf_lengths group by name) as fmax
	where rfl.name=cmin.name and rfl.name=cmax.name
	AND rfl.name=fmin.name and rfl.name=fmax.name
--	AND rfl.name=rfl_cmin.name AND cmin.min=rfl_cmin.core_length
--	AND rfl.name=rfl_cmax.name AND cmax.max=rfl_cmax.core_length
--	AND rfl.name=rfl_fmin.name AND fmin.min=rfl_fmin.full_length
--	AND rfl.name=rfl_fmax.name AND fmax.max=rfl_fmax.full_length
	group by rfl.name with rollup; -- as 'Genome wide';

-- We get some truncation warnings here, but the data looks good?
-- Should really define table beforehand

alter table rf_stats add column min_core_id  int(10) unsigned NOT NULL;
alter table rf_stats add column max_core_id  int(10) unsigned NOT NULL;
alter table rf_stats add column min_full_id  int(10) unsigned NOT NULL;
alter table rf_stats add column max_full_id  int(10) unsigned NOT NULL;


-- We may have some duplicated lengths here
-- Would be better to show the top 5 lengths, ids, fset_ids?
update rf_stats rfs, rf_lengths rfl set rfs.min_core_id=rfl.regulatory_feature_id where rfl.name=rfs.name and rfs.min_core_length=rfl.core_length;
update rf_stats rfs, rf_lengths rfl set rfs.max_core_id=rfl.regulatory_feature_id where rfl.name=rfs.name and rfs.max_core_length=rfl.core_length;
update rf_stats rfs, rf_lengths rfl set rfs.min_full_id=rfl.regulatory_feature_id where rfl.name=rfs.name and rfs.min_full_length=rfl.full_length;
update rf_stats rfs, rf_lengths rfl set rfs.max_full_id=rfl.regulatory_feature_id where rfl.name=rfs.name and rfs.max_full_length=rfl.full_length;


-- 2027013, 2068586, 2474712, 1970333
--select * from misc_feature mf, misc_attrib ma where mf.misc_feature_id = ma.misc_feature_id and ma.value = 'Encode Duke excluded regions' limit 10;


select rfl.regulatory_feature_id, rfl.name, rfl.seq_region_start, rfl.seq_region_end, bl.seq_region_start as bl_start, bl.seq_region_end as bl_end 
	from rf_lengths rfl,
	(select sr.name, mf.seq_region_start, mf.seq_region_end 
		from homo_sapiens_core_57_37b.misc_feature mf, homo_sapiens_core_57_37b.misc_attrib ma, homo_sapiens_core_57_37b.seq_region sr 
		where mf.misc_feature_id = ma.misc_feature_id 
		and ma.value = 'Encode Duke excluded regions' 
		and mf.seq_region_id=sr.seq_region_id) bl 
	where bl.name=rfl.name 
	and bl.seq_region_end>=rfl.seq_region_start 
	and bl.seq_region_start<=rfl.seq_region_end; 
--	and rfl.regulatory_feature_id in(2027013, 2068586, 2474712, 1970333);


--We may get duplicaed lines if we have more than one feature with the same min or max length
--could change this to group concat the longest/shortest 5 ids/lengths

--Would we be better generating this through a few queries in a tmp table?


--select af.* from annotated_feature af, regulatory_attribute rf where rf.regulatory_feature_id=1970333 and rf.attribute_feature_id=af.annotated_feature_id and af.feature_set_id in(149,20,21,22) order by seq_region_start;
-- For the Y chr proglem feature this show us that the focus features do not over lap. So are they being integrated by a larger containing attribute feature?

-- can't order by with rollup
-- need extended join syntax here for including mutltiple table sin a left join statement
-- otherwise we get weird errors:ERROR 1054 (42S22): Unknown column 'rf.seq_region_id' in 'on clause'

SELECT sr.name, AVG(rf.seq_region_end - rf.seq_region_start + 1) as 'Average core length', AVG(rf.bound_seq_region_end - rf.bound_seq_region_start + 1) as 'Average full length' from regulatory_feature rf left join (select distinct(seq_region_id), name from seq_region) sr on sr.seq_region_id=rf.seq_region_id join feature_set fs on fs.feature_set_id=rf.feature_set_id and fs.name='RegulatoryFeatures_v6' group by sr.name with rollup; 





--#CREATE TEMPORARY TABLE rf_medians SELECT (rf1.seq_region_end - rf1.seq_region_start + 1) as rl1 FROM
--#data x, data y GROUP BY x.val HAVING
--#((COUNT(*)-(SUM(y.val-x.val))))) <=
--#floor((COUNT(*) +1)/2)) and
--#((COUNT(*)-(SUM(SIGN(1+SIGN(y.val-x.val))))) <=
--#floor((COUNT(*) +1)/2));
