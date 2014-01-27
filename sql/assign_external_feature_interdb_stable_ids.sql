/** 
@header Assign_new_interdb_stable_ids
@desc   Creates new interdb_stable_id values for NULL
        values, starting at the next logical value.
*/


-- No point in doing this as we can't guarantee they won't be reused between releases
-- if we selectively delete some of the last IDs
-- set @n =(select max(interdb_stable_id) +1 from motif_feature);

update external_feature set interdb_stable_id=NULL;
set @n = 0;

update external_feature mf
  join (select external_feature_id, @n := @n + 1 new_stable_id
          from external_feature where interdb_stable_id is NULL
          order by external_feature_id) v
    on mf.external_feature_id = v.external_feature_id
  set mf.interdb_stable_id = v.new_stable_id;

