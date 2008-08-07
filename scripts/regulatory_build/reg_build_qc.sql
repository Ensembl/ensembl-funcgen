 select count(distinct(rf.regulatory_feature_id)), sr.name, fs.name from regulatory_feature rf, feature_set fs, seq_region sr where fs.name in ('RegulatoryFeatures', 'RegulatoryFeatures_v3') and fs.feature_set_id=rf.feature_set_id and rf.seq_region_id=sr.seq_region_id group by rf.feature_set_id, rf.seq_region_id order by rf.seq_region_id, rf.feature_set_id;


--This takes a while and can probably be done better


select fs.name, count(*)  from regulatory_feature rf, feature_set fs where fs.feature_set_id in(69,112) and fs.feature_set_id=rf.feature_set_id group by rf.feature_set_id;


-- post analysis feature type QC




select count(distinct(rf.regulatory_feature_id)), ft.name, fs.name from regulatory_feature rf, feature_set fs, feature_type ft where fs.name in ('RegulatoryFeatures', 'RegulatoryFeatures_v3') and fs.feature_set_id=rf.feature_set_id and rf.feature_type_id=ft.feature_type_id group by rf.feature_set_id, rf.feature_type_id order by rf.feature_type_id, rf.feature_set_id;
