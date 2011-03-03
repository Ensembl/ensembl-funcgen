/** 
@header patch_61_62_c.sql - interdb_stable_id
@desc   Add 'internal' motif_feature/external_feature stable_ids.  
        Required for inter-DB linking without using internal db ids
        i.e variation consequences
*/


ALTER table motif_feature ADD `interdb_stable_id` mediumint(8) unsigned DEFAULT NULL;

set @n = 0;
update motif_feature mf
  join (select motif_feature_id, @n := @n + 1 new_stable_id
          from motif_feature
          order by motif_feature_id) v
    on mf.motif_feature_id = v.motif_feature_id
  set mf.interdb_stable_id = v.new_stable_id;

ALTER table motif_feature ADD UNIQUE KEY `interdb_stable_id_idx` (`interdb_stable_id`);


analyze table motif_feature;
optimize table motif_feature;


ALTER table external_feature ADD `interdb_stable_id` mediumint(8) unsigned DEFAULT NULL;

set @n = 0;
update external_feature mf
  join (select external_feature_id, @n := @n + 1 new_stable_id
          from external_feature
          order by external_feature_id) v
    on mf.external_feature_id = v.external_feature_id
  set mf.interdb_stable_id = v.new_stable_id;

ALTER table external_feature ADD UNIQUE KEY `interdb_stable_id_idx` (`interdb_stable_id`);


analyze table external_feature;
optimize table external_feature;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_61_62_b.sql|interdb_stable_id');
