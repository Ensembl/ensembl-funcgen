update regulatory_feature set activity_as_enum = 
case activity
    when 0 then 'INACTIVE'
    when 1 then 'ACTIVE'
    when 2 then 'POISED'
    when 3 then 'REPRESSED'
    when 4 then 'NA'
    else null
end;

alter table regulatory_feature drop column activity;
alter table regulatory_feature change activity_as_enum activity ENUM('INACTIVE', 'REPRESSED', 'POISED', 'ACTIVE', 'NA');

-- patch identifier
insert into meta (species_id, meta_key, meta_value) values (null, 'patch', 'patch_84_85_j.sql|Make activity an enum.');