-- update coord_system_ids and meta_coord and meta table
-- This will disappear when we implement version to assembly mapping 
-- Tidy up analysis table and ids;

-- update meta table
delete from meta where meta_key="schema.version";
insert into meta values('', 'schema_version', '43');

-- update coord_system_table dependent on species
select "Need to manually update the schema_build in coord_system to match the species schema_build";
--update coord_system set schema_build='43_36e' where schema_build='42_36d';
--update coord_system set schema_build='43_36d' where schema_build='42_36c';



-- Tidy up analysis tables/ids
update predicted_feature pf, analysis a set pf.analysis_id=a.analysis_id where a.logic_name='Nessie';
delete a, ad from analysis a, analysis_description ad where a.logic_name='TilingHMM' and a.analysis_id=ad.analysis_id ;
