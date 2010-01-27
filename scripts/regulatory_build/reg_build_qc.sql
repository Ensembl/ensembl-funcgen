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
select count(distinct(rf.regulatory_feature_id)), sr.name, fs.name, count(distinct(rf.regulatory_feature_id)) from regulatory_feature rf, feature_set fs, seq_region sr where fs.name in ('RegulatoryFeatures', 'RegulatoryFeatures_v5') and fs.feature_set_id=rf.feature_set_id and rf.seq_region_id=sr.seq_region_id group by rf.feature_set_id, rf.seq_region_id order by sr.name, rf.feature_set_id;

--This takes a while and can probably be done better
--maybe do the nested left join where select distinct sr's as above


-- Total counts across feature_sets
select fs.name, count(*)  from regulatory_feature rf, feature_set fs where fs.feature_set_id in(69,112) and fs.feature_set_id=rf.feature_set_id group by rf.feature_set_id;


-- post analysis feature type QC


select count(distinct(rf.regulatory_feature_id)), ft.name, fs.name from regulatory_feature rf, feature_set fs, feature_type ft where fs.name in ('RegulatoryFeatures', 'RegulatoryFeatures_v4') and fs.feature_set_id=rf.feature_set_id and rf.feature_type_id=ft.feature_type_id group by rf.feature_set_id, rf.feature_type_id order by rf.feature_type_id, rf.feature_set_id;


-- by chr?

-- This counts all inclusions of regulatory attributes, so there may be duplicates
-- From inclusions in two separate features
-- Is there a way around this? distinct(af.annotated_feature_id)?

 select fs1.name, sr.name, fs1.feature_set_id, count(*) from annotated_feature af, regulatory_attribute ra, feature_set fs, feature_set fs1, regulatory_feature rf left join (select distinct(seq_region_id), name from seq_region) sr on sr.seq_region_id=af.seq_region_id where fs.name='RegulatoryFeatures' and fs.feature_set_id=rf.feature_set_id and rf.regulatory_feature_id=ra.regulatory_feature_id and ra.attribute_feature_id=af.annotated_feature_id and ra.attribute_feature_table='annotated' and af.feature_set_id=fs1.feature_set_id group by fs1.name, sr.name order by sr.name, fs1.name;



