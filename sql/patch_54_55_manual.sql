# patch_54_55_c.sql
#
# title: manual patch
#
# description:
# This contains a list of manual patches which are species or data dependant


select "This patch dependant on data or species, please perform this manually";
exit

#coord_system tidy
#This is only for those species which do not have versions for unassembled level
#And have been integrated into eFG with inheriting the assembled levels version

update mus_musculus_funcgen_55_37h.coord_system set version=NULL where name not in ('scaffold', 'chromosome');
update mus_musculus_funcgen_55_37h.coord_system set attrib='default_version' where name='supercontig';




INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_54_55_manual.sql|manual');

