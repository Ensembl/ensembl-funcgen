/** 
@header Assign_new_interdb_stable_ids
@desc   Creates new interdb_stable_id values for NULL
        values, starting at the next logical value.
*/


SELECT "YOU MUST EDIT THIS SQL MANUALLY TO POINT TO THE RELEVANT TABLE";

exit;

-- Valid tables are current motif_feature or external_feature

set @n =(select max(interdb_stable_id) +1 from motif_feature);
update motif_feature mf
  join (select motif_feature_id, @n := @n + 1 new_stable_id
          from motif_feature where interdb_stable_id is NULL
          order by motif_feature_id) v
    on mf.motif_feature_id = v.motif_feature_id
  set mf.interdb_stable_id = v.new_stable_id;

analyze table motif_feature;
optimize table motif_feature;
