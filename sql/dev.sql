-- List features above 2kb.



--SELECT (t.seq_region_end - t.seq_region_start + 1 ) as max, t.annotated_feature_id, t.feature_set_id FROM annotated_feature t where max >2000 order by max desc limit 100;

-- Can't use alias in where/order

-- Ignore eQTLs, these should really be in a different table if they are likely to be much larger than other features
-- This will speed up the range queries

 SELECT (t.seq_region_end - t.seq_region_start + 1 ) as max, t.annotated_feature_id, t.feature_set_id FROM annotated_feature t where (t.seq_region_end - t.seq_region_start + 1 ) >2000 and t.feature_set_id!=70 order by (t.seq_region_end - t.seq_region_start + 1 ) desc;

-- get feature_sets effected
-- SELECT distinct(fs.name), t.annotated_feature_id, t.feature_set_id FROM annotated_feature t where (t.seq_region_end - t.seq_region_start + 1 ) >2000 order by (t.seq_region_end - t.seq_region_start + 1 ) desc;



-- Cannot do this for probe feature as we have the sanger PCR array, need to do a not in for the chip_cannel_ids.


-----------------------------------

--- associated_feature_types
--- should really change associated_feature_type to use a number for feature_table?
--- This is only true if we are actually matching the field on a query?
--- This will be thes case when we are retrieving using a Feature
--- But not when we are retrieving using a FeatureType
--- Main usage id by Feature, so number would be better
--- Not the case for RegFeat as we use the reg_feat_id as the key, not the table_id/name
--- We need an extra index on feature_type_id.
--- Don't want to obfuscate table_names in DB, maybe we could have them in meta?
--- And get the SetFeature adaptor to load them when instantiated
--- We would lose some speed on converting numbers to text when retrieving reg_attrs
--- Or pulling features back by their group feature_types

--- Is this always going to represent groups, in the sense of a motif group?
--- Or can this be used generically as simply a group of feature_types?
--- This is simply associating multiple feature types with the same feature
--- Maybe feature_type_set, no confusion with other sets.
--- feature_attribute? Does not have feature_type focus and maybe conflict with the idea of a core attribute.
--- Which 

--- groups should not really be added as FeatureTypes as they are a little more abstract.
--- Less of a region on the genome and more of a classification of a group of feature_types.


--- regulatory_attribute should be renamed regulatory_supporting_feature?
--- Then we could have annotated_supporting_feature..if required, or do we just consider all combined
--- features a regulatory_feature, and make the distinction with the feature set?

--- in core terms attrs can be a count of somthing e.g has a value, can be 1 or another metric.

--- So we're really talking about classes here, which may overlap, kinda like a DAG.
--- Then we have ambiguity with feature_type.class

--- Just implement this for now and think about a more generic solution?
--- methods would be get_FeatureTypes? get_attribute_FeatureTypes?
--- get_group_FeatureTypes get_associated_FeatureTypes?




----------- array


--	 #This is used to store Arrays

--	 #Do we need to re-introduce class here?
--	 #EXPRESSION is the usage really. i.e. We can use a tiling array for an expression experiment
--	 #but the format is definitely tiling, but then we can also have targetting tiling or whole genome, cnv etc.
--	 #Could class be overall design focus. i.e. if it was designed with expression in mind 
--	 #then we would use expression.
--	 #So we could have VENDOR CLASS FORMAT
--	 #AFFY EXPRESSION ST
--	 #AFFY EXPRESSION UTR
--	 #NIMBLEGEN TILING WHOLE_GENOME
--	 #NIMBLEGEN TILING TARGETTED
--	 #ILLUMINA ??? ???
--	 #CODELINK ??? ???

--	 #Let's take a look at some other vendors before we decide on this

 -
--	 #These formats are differing between generic and vendor specific designations

--	 #Is type redundant? Only sanger use PCR probes, any other types of probe?

-- Remeber to change ImportArrays config if we change this, as well as Parsers




--- Affy Exon array checks

--|             38 |            NULL |             4 | HuEx_1_0_st_v2   | AFFY |
--|             39 |            NULL |             4 | HuGene_1_0_st_v1 | AFFY |

mysql> select count(distinct(probeset)) from oligo_probe where oligo_array_id in(38,39) group by oligo_array_id ;
+---------------------------+
| count(distinct(probeset)) |
+---------------------------+
|                   1432033 |
|                     33843 |
+---------------------------+

mysql> select * from external_db where db_name like 'AFFY%Hu%st%';
+----------------+-----------------------+------------+--------+------------------------+------------------------+----------+----------------------------------------+-------+-------------------+--------------------+
| external_db_id | db_name               | db_release | status | dbprimary_acc_linkable | display_label_linkable | priority | db_display_name                        | type  | secondary_db_name | secondary_db_table |
+----------------+-----------------------+------------+--------+------------------------+------------------------+----------+----------------------------------------+-------+-------------------+--------------------+
|           3135 | AFFY_HuEx_1_0_st_v2   | 1          | XREF   |                      1 |                      0 |        1 | Affymx Microarray Human Exon 1.0 ST v2 | ARRAY | NULL              | NULL               |
|           3136 | AFFY_HuGene_1_0_st_v1 | 1          | XREF   |                      1 |                      0 |        1 | Affymx Microarray Human Gene 1.0 ST    | ARRAY | NULL              | NULL               |
+----------------+-----------------------+------------+--------+------------------------+------------------------+----------+----------------------------------------+-------+-------------------+--------------------+

mysql> select count(distinct(x.xref_id)) from xref x, object_xref ox where x.external_db_id in(3135,3136) and x.xref_id=ox.xref_id group by x.external_db_id;
+----------------------------+
| count(distinct(x.xref_id)) |
+----------------------------+
|                     112104 |
|                       1557 |
+----------------------------+
2 rows in set (1.99 se


-- this is just the number of mapped probesets.

-- same for ox.xref_id i.e. the number of probeset to transcript mappings. So each is just mapped once




--- Count type of unmapped reasons
---This will do it for all affy.
select ur.summary_description, count(uo.unmapped_reason_id) from unmapped_reason ur, unmapped_object uo where uo.identifier like '%:%' and uo.type='probe2transcript' and uo.unmapped_reason_id=ur.unmapped_reason_id group by ur.unmapped_reason_id;

-- This will restrict to st arrays, but maybe inlfated counts due to product of oligo_probe table, 4 probes per probeset.
select oa.name, ur.summary_description, count(uo.unmapped_reason_id) from oligo_array oa, oligo_probe op, unmapped_reason ur, unmapped_object uo where oa.name like '%st%' and oa.oligo_array_id=op.oligo_array_id and uo.identifier like concat(op.probeset, ':%') and uo.type='probe2transcript' and uo.unmapped_reason_id=ur.unmapped_reason_id group by ur.unmapped_reason_id;

-- We could simply do a direct match instead of a like for the whole probeset.


select oa.name, ur.summary_description, count(uo.unmapped_reason_id) from oligo_array oa, oligo_probe op, unmapped_reason ur, unmapped_object uo where oa.name like '%st%' and oa.oligo_array_id=op.oligo_array_id and uo.identifier=concat(op.probeset, ':', op.oligo_probe_id) and uo.type='probe2transcript' and uo.unmapped_reason_id=ur.unmapped_reason_id group by ur.unmapped_reason_id, oa.name;


-- So these are probe level unmapped_object
-- Need to add probeset level transcript unmapped_object to get probeset level counts.

v50

+------------------+---------------------+------------------------------+
| name             | summary_description | count(uo.unmapped_reason_id) |
+------------------+---------------------+------------------------------+
| HuEx_1_0_st_v2   | Unmapped anti-sense |                      9266349 |
| HuGene_1_0_st_v1 | Unmapped anti-sense |                      2622990 |
| HuEx_1_0_st_v2   | Unmapped intronic   |                      2403697 |
| HuGene_1_0_st_v1 | Unmapped intronic   |                       299899 |
+------------------+---------------------+------------------------------+


v51

+------------------+---------------------+------------------------------+
| name             | summary_description | count(uo.unmapped_reason_id) |
+------------------+---------------------+------------------------------+
| HuEx_1_0_st_v2   | Unmapped anti-sense |                      8884857 |
| HuGene_1_0_st_v1 | Unmapped anti-sense |                      2542473 |
| HuEx_1_0_st_v2   | Unmapped intronic   |                      2252280 |
| HuGene_1_0_st_v1 | Unmapped intronic   |                       271810 |
+------------------+---------------------+------------------------------+
4 rows in set (3 hours 46 min 34.75 sec)




-- Get total probeset unmapped counts


--Note: when using distinct, have to have distinct fields called first, or MySQL gives an error????


 mysql> select t.name, ur.summary_description, count(uo.unmapped_reason_id) from unmapped_reason ur, unmapped_object uo, (select distinct(op.probeset) as probeset, oa.name as name from oligo_array oa, oligo_probe op  where oa.name like '%st%' and oa.oligo_array_id=op.oligo_array_id) as t where uo.identifier=t.probeset and uo.type='probe2transcript' and uo.unmapped_reason_id=ur.unmapped_reason_id group by ur.summary_description, t.name;            
+------------------+------------------------+------------------------------+
| name             | summary_description    | count(uo.unmapped_reason_id) |
+------------------+------------------------+------------------------------+
| HuEx_1_0_st_v2   | Insufficient hits      |                        30044 |
| HuGene_1_0_st_v1 | Insufficient hits      |                         9959 |
| HuEx_1_0_st_v2   | No transcript mappings |                      1361347 |
| HuGene_1_0_st_v1 | No transcript mappings |                        32751 |
+------------------+------------------------+------------------------------+




v50 oligo_probe gene st counts


Seems we have duplicate probes, altho this should not effect mappings, just relative counts?

mysql> select * from oligo_probe where probeset =7893306 limit 10;
+----------------+----------------+----------+--------+-------------+--------+
| oligo_probe_id | oligo_array_id | probeset | name   | description | length |
+----------------+----------------+----------+--------+-------------+--------+
|        6770764 |             39 | 7893306  | 793829 | NULL        |     25 |
|        6774819 |             39 | 7893306  | 793829 | NULL        |     25 |
|        6775012 |             39 | 7893306  | 793829 | NULL        |     25 |
|        6789237 |             39 | 7893306  | 793829 | NULL        |     25 |
|        6789368 |             39 | 7893306  | 793829 | NULL        |     25 |
|        6791621 |             39 | 7893306  | 793829 | NULL        |     25 |
|        6805857 |             39 | 7893306  | 793829 | NULL        |     25 |
|        6809391 |             39 | 7893306  | 793829 | NULL        |     25 |
|        6821455 |             39 | 7893306  | 793829 | NULL        |     25 |
|        6829290 |             39 | 7893306  | 793829 | NULL        |     25 |
+----------------+----------------+----------+--------+-------------+--------+
10 rows in set (6.10 sec)


-- Total probe records
mysql> select count(name) from oligo_probe where oligo_array_id =39;
+-------------+
| count(name) |
+-------------+
|      862560 |
+-------------+

-- Total probesets
mysql> select count(distinct(probeset)) from oligo_probe where oligo_array_id =39;
+---------------------------+
| count(distinct(probeset)) |
+---------------------------+
|                     33843 |
+---------------------------+
1 row in set (6.89 sec)

--should be 4 probes per set?

mysql> select distinct(name) from oligo_probe where probeset =7893306 and oligo_array_id=39;
+--------+
| name   |
+--------+
| 793829 |
+--------+

-- This is a control probe set! Do we want to map these?
>probe:HuGene-1_0-st-v1:793829;536:198; TranscriptClusterID=7893306; Sense; ProbeSetType=control->affx
>probe:HuGene-1_0-st-v1:793829;952:1010; TranscriptClusterID=7893306; Sense; ProbeSetType=control->affx
>probe:HuGene-1_0-st-v1:793829;697:513; TranscriptClusterID=7893306; Sense; ProbeSetType=control->affx


--These are not being imported properly, using the TranscriptClusterID as the probeset value, this is true for normal probes
--The name is not being parsed properly. Do we want the TranscriptClusterID here or 793829.
-- So the duplicates we are getting are most likely all from control sets. Therefore not the end of the world.
-- Still don't know why we aren't mapping many.

mysql> select * from oligo_probe where probeset =7893306 limit 10;
+----------------+----------------+----------+--------+-------------+--------+
| oligo_probe_id | oligo_array_id | probeset | name   | description | length |
+----------------+----------------+----------+--------+-------------+--------+
|        6770764 |             39 | 7893306  | 793829 | NULL        |     25 |
|        6774819 |             39 | 7893306  | 793829 | NULL        |     25 |
|        6775012 |             39 | 7893306  | 793829 | NULL        |     25 |




-- Just the one, but over a thousand of them !!? 

-- let's try with a different one.

mysql> select distinct(name) from oligo_probe where probeset =8145611 and oligo_array_id=39;
+-----------------+
| name            |
+-----------------+
| 1024586-8145611 |
| 268291-8145611  |
| 59368-8145611   |
| 311496-8145611  |
| 197117-8145611  |
| 870225-8145611  |
| 816217-8145611  |
| 700852-8145611  |
| 1050452-8145611 |
| 362470-8145611  |
| 38088-8145611   |
| 656514-8145611  |
| 224779-8145611  |
| 334249-8145611  |
| 624033-8145611  |
| 1010192-8145611 |
| 840073-8145611  |
| 883823-8145611  |
| 683325-8145611  |
| 361210-8145611  |
| 571709-8145611  |
| 588393-8145611  |
| 625134-8145611  |
| 1068184-8145611 |
| 71224-8145611   |
| 17316-8145611   |
| 29014-8145611   |
| 318823-8145611  |
| 785869-8145611  |
| 191191-8145611  |
| 579491-8145611  |
| 1034303-8145611 |
| 775529-8145611  |
| 1035837-8145611 |
| 1081813-8145611 |
| 967305-8145611  |
| 719983-8145611  |
+-----------------+
37 rows in set (14.68 sec)

-- So 36 probes per set for some gene probesets? Some have 23, 24.


-- Now let's try the exon array. All 4.


-- Count probeset sizes for gene st

select t.name_count, count(*) from(select count(name) as name_count, oligo_array_id from oligo_probe where oligo_array_id in(38) group by probeset) as t group by t.name_count;

+------------+----------+
| name_count | count(*) |
+------------+----------+
|          1 |    64509 |
|          2 |    43631 |
|          3 |    42110 |
|          4 |  1281576 |
|          5 |        1 |
|          7 |        1 |
|          8 |        1 |
|          9 |        2 |
|         10 |        1 |
|         11 |        2 |
|         13 |        4 |
|         14 |        3 |
|         15 |        5 |
|         16 |        2 |
|         21 |        1 |
|         22 |       21 |
|         23 |        1 |
|         26 |        1 |
|         40 |        9 |
|         41 |        2 |
|         44 |        1 |
|         52 |        1 |
|         53 |        1 |
|         57 |        1 |
|         60 |        1 |
|         68 |        2 |
|         69 |        1 |
|         70 |        2 |
|         71 |        2 |
|         73 |        1 |
|         74 |        1 |
|         75 |        2 |
|         78 |        1 |
|         79 |        2 |
|         80 |        2 |
|         81 |        1 |
|         82 |        1 |
|         84 |        3 |
|         85 |        1 |
|         86 |        3 |
|         87 |        4 |
|         88 |        5 |
|         89 |        4 |
|         90 |        8 |
|         91 |       18 |
|         92 |        8 |
|         93 |       12 |
|         94 |        8 |
|         95 |        3 |
|         96 |        6 |
|         97 |        2 |
|        132 |        1 |
|        136 |        1 |
|        160 |        1 |
|        162 |        1 |
|        203 |        1 |
|        204 |        1 |
|        285 |        1 |
|        300 |        1 |
|        347 |        1 |
|        348 |        1 |
|        350 |        1 |
|        355 |        1 |
|        406 |        1 |
|        407 |        1 |
|        419 |        1 |
|        430 |        1 |
|        431 |        1 |
|        442 |        1 |
|        454 |        1 |
|        456 |        1 |
|        458 |        2 |
|        463 |        1 |
|        468 |        1 |
|        470 |        1 |
|        472 |        1 |
|        474 |        1 |
|        475 |        2 |
|        477 |        1 |
|        479 |        1 |
|        480 |        3 |
|        481 |        1 |
|        483 |        2 |
|        484 |        3 |
|        485 |        1 |
|        489 |        1 |




Now for gene



+------------+----------+
| name_count | count(*) |
+------------+----------+
|          1 |      334 |
|          2 |      356 |
|          3 |      304 |
|          4 |     4165 |
|          5 |       84 |
|          6 |       90 |
|          7 |       66 |
|          8 |      186 |
|          9 |      193 |
|         10 |      176 |
|         11 |      183 |
|         12 |      232 |
|         13 |      212 |
|         14 |      243 |
|         15 |      273 |
|         16 |      267 |
|         17 |      265 |
|         18 |      338 |
|         19 |      380 |
|         20 |      515 |
|         21 |      649 |
|         22 |      747 |
|         23 |      953 |
|         24 |     2063 |
|         25 |     4725 |
|         26 |     1758 |
|         27 |     1384 |
|         28 |     1526 |
|         29 |      934 |
|         30 |     1277 |
|         31 |      822 |
|         32 |     1068 |
|         33 |      807 |
|         34 |      703 |
|         35 |      485 |
|         36 |      736 |
|         37 |      354 |
|         38 |      422 |
|         39 |      283 |
|         40 |      408 |
|         41 |      265 |
|         42 |      323 |
|         43 |      209 |
|         44 |      292 |
|         45 |      196 |
|         46 |      203 |
|         47 |      151 |
|         48 |      138 |
|         49 |       42 |
|         50 |       66 |
|         51 |       55 |
|         52 |       67 |
|         53 |       48 |
|         54 |       59 |
|         55 |       37 |
|         56 |       43 |
|         57 |       41 |
|         58 |       37 |
|         59 |       27 |
|         60 |       35 |
|         61 |       27 |
|         62 |       34 |
|         63 |       14 |
|         64 |       29 |
|         65 |       31 |
|         66 |       19 |
|         67 |       19 |
|         68 |       20 |
|         69 |       12 |
|         70 |       24 |
|         71 |       17 |
|         72 |       11 |
|         73 |       17 |
|         74 |       11 |
|         75 |       11 |
|         76 |        5 |
|         77 |       15 |
|         78 |        3 |
|         79 |       13 |
|         80 |        6 |
|         81 |        6 |
|         82 |        5 |
|         83 |        6 |
|         84 |        8 |
|         85 |        6 |
|         86 |       11 |
|         87 |        5 |
|         88 |        6 |
|         89 |        5 |
|         90 |        5 |
|         91 |        4 |
|         92 |        2 |
|         93 |        7 |
|         94 |        3 |
|         95 |        3 |
|         96 |        2 |
|         97 |        1 |
|         98 |        5 |
|         99 |        1 |
|        100 |        2 |
|        101 |        2 |
|        102 |        1 |
|        104 |        2 |
|        105 |        2 |
|        106 |        3 |
|        107 |        1 |
|        108 |        2 |
|        109 |        4 |
|        110 |        2 |
|        111 |        4 |
|        112 |        3 |
|        113 |        3 |
|        114 |        1 |
|        115 |        1 |
|        116 |        1 |
|        117 |        4 |
|        118 |        1 |
|        120 |        1 |
|        121 |        2 |
|        123 |        1 |
|        124 |        1 |
|        125 |        3 |
|        132 |        1 |
|        133 |        3 |
|        134 |        1 |
|        136 |        1 |
|        140 |        1 |
|        141 |        1 |
|        144 |        1 |
|        149 |        2 |
|        151 |        1 |
|        153 |        1 |
|        156 |        2 |
|        160 |        1 |
|        162 |        1 |
|        180 |        1 |
|        184 |        1 |
|        187 |        1 |
|        203 |        1 |
|        204 |        1 |
|        210 |        1 |
|        242 |        1 |
|        285 |        1 |
|        300 |        1 |
|        315 |        1 |
|        347 |        1 |
|        348 |        1 |
|        350 |        1 |
|        355 |        1 |
|        380 |        1 |
|        406 |        1 |
|        407 |        1 |
|        417 |        1 |
|        419 |        1 |
|        430 |        1 |
|        431 |        1 |
|        442 |        1 |
|        454 |        1 |
|        456 |        1 |
|        458 |        2 |
|        463 |        1 |
|        468 |        1 |
|        470 |        1 |
|        472 |        1 |
|        474 |        1 |
|        475 |        2 |
|        477 |        1 |
|        479 |        1 |
|        480 |        3 |
|        481 |        1 |
|        483 |        2 |
|        484 |        3 |
|        485 |        1 |
|        489 |        1 |
|        567 |        1 |
|        844 |        1 |
|       1160 |        1 |
|       1189 |        1 |
+------------+----------+



-- Now let's account fo rduplciate probe entries

-- exon
No change in counts

-- gene
-- sizes unchanged but quite a few duplicates removed.

+------------+----------+
| name_count | count(*) |
+------------+----------+
|          1 |     4565 |
|          2 |      207 |
|          3 |      157 |
|          4 |      505 |
|          5 |       84 |
|          6 |       80 |
|          7 |       63 |
|          8 |      182 |
|          9 |      189 |
|         10 |      173 |
|         11 |      176 |
|         12 |      230 |
|         13 |      205 |
|         14 |      235 |
|         15 |      270 |
|         16 |      261 |
|         17 |      260 |
|         18 |      333 |
|         19 |      371 |
|         20 |      499 |
|         21 |      632 |
|         22 |      740 |
|         23 |      946 |
|         24 |     2047 |
|         25 |     4713 |
|         26 |     1744 |
|         27 |     1373 |
|         28 |     1508 |
|         29 |      928 |
|         30 |     1269 |
|         31 |      815 |
|         32 |     1064 |
|         33 |      797 |
|         34 |      701 |
|         35 |      482 |
|         36 |      733 |
|         37 |      349 |
|         38 |      422 |
|         39 |      281 |
|         40 |      405 |
|         41 |      265 |
|         42 |      319 |
|         43 |      208 |
|         44 |      291 |
|         45 |      196 |
|         46 |      203 |
|         47 |      149 |
|         48 |      135 |
|         49 |       40 |
|         50 |       66 |
|         51 |       54 |
|         52 |       67 |
|         53 |       48 |
|         54 |       59 |
|         55 |       37 |
|         56 |       42 |
|         57 |       41 |
|         58 |       37 |
|         59 |       27 |
|         60 |       35 |
|         61 |       27 |
|         62 |       34 |
|         63 |       14 |
|         64 |       29 |
|         65 |       29 |
|         66 |       19 |
|         67 |       19 |
|         68 |       20 |
|         69 |       12 |
|         70 |       24 |
|         71 |       17 |
|         72 |       11 |
|         73 |       17 |
|         74 |       11 |
|         75 |       11 |
|         76 |        5 |
|         77 |       15 |
|         78 |        3 |
|         79 |       13 |
|         80 |        5 |
|         81 |        6 |
|         82 |        5 |
|         83 |        6 |
|         84 |        8 |
|         85 |        6 |
|         86 |       11 |
|         87 |        5 |
|         88 |        6 |
|         89 |        5 |
|         90 |        5 |
|         91 |        4 |
|         92 |        2 |
|         93 |        7 |
|         94 |        3 |
|         95 |        3 |
|         96 |        2 |
|         97 |        1 |
|         98 |        4 |
|         99 |        1 |
|        100 |        2 |
|        101 |        1 |
|        102 |        1 |
|        104 |        2 |
|        105 |        2 |
|        106 |        3 |
|        107 |        1 |
|        108 |        2 |
|        109 |        3 |
|        110 |        2 |
|        111 |        3 |
|        112 |        3 |
|        113 |        3 |
|        114 |        1 |
|        115 |        1 |
|        116 |        1 |
|        117 |        4 |
|        118 |        1 |
|        120 |        1 |
|        121 |        1 |
|        123 |        1 |
|        124 |        1 |
|        125 |        3 |
|        132 |        1 |
|        133 |        3 |
|        134 |        1 |
|        136 |        1 |
|        140 |        1 |
|        141 |        1 |
|        144 |        1 |
|        149 |        2 |
|        151 |        1 |
|        153 |        1 |
|        156 |        2 |
|        160 |        1 |
|        162 |        1 |
|        180 |        1 |
|        184 |        1 |
|        187 |        1 |
|        203 |        1 |
|        204 |        1 |
|        210 |        1 |
|        285 |        1 |
|        300 |        1 |
|        315 |        1 |
|        347 |        1 |
|        348 |        1 |
|        350 |        1 |
|        355 |        1 |
|        380 |        1 |
|        406 |        1 |
|        407 |        1 |
|        417 |        1 |
|        419 |        1 |
|        430 |        1 |
|        431 |        1 |
|        442 |        1 |
|        454 |        1 |
|        456 |        1 |
|        458 |        2 |
|        463 |        1 |
|        468 |        1 |
|        470 |        1 |
|        472 |        1 |
|        474 |        1 |
|        475 |        2 |
|        477 |        1 |
|        479 |        1 |
|        480 |        3 |
|        481 |        1 |
|        483 |        2 |
|        484 |        3 |
|        485 |        1 |
|        489 |        1 |

-- So let's take a look at some unmapped objects
-- Is this because we have reduced the extension?
-- ST Gene probes certainly seem designed against the 3' UTR, and quite a bit downstream.
-- Is this case for re-introducing an extension.


-- We should really build stats on this when we run the mapping, add upsteam/downstream calls?
-- This would mean we'd have to bring back all the oligo_features for a probeset, even if they don't overlap.
-- We want to build these stats anyway for ProbeFeature view.
-- Do we need to refactor the probemapping pipeline such that we do things on a probeset centric basis, rather than a transcript centric basis. Then we can parallelise on the farm!!!!
-- What caches do we need to account for? Are there any which we could not build this way?

-- We would have to be mindful about the potential for the same probeset on different arrays having slightly different probes.
-- Can't do this easily with core due to lack of ProbeSet API, currently eFG ProbeSet not directly linked to array, so hard to pull back based on array. Need link table, probeset_array_chip or probeset_array(in case we have probesets spread across different array_chips) ProbeSet would have to have functionality to figure out if it varied across arrays. This would be dependant on probes being collapsed into the same record properly.

-- Do we need different rules for st arrays? We should follow the biological model, not the array model, so should be same for all arrays.

-- Maybe this is an exon boundary issue.

-- let's take a closer look at some that have no transcript mappings at all.


mysql> select * from unmapped_object where unmapped_reason_id =750 limit 10;
+--------------------+------------------+-------------+----------------+------------+--------------------+-------------+--------------+------------+---------
------------+--------+
| unmapped_object_id | type             | analysis_id | external_db_id | identifier | unmapped_reason_id | query_score | target_score | ensembl_id | ensembl_
object_type | parent |
+--------------------+------------------+-------------+----------------+------------+--------------------+-------------+--------------+------------+---------
------------+--------+
|           24794856 | probe2transcript |        5253 |           NULL | 3937834    |                750 |        NULL |         NULL |       NULL | NULL    
            | NULL   |
|           24794857 | probe2transcript |        5253 |           NULL | 2360807    |                750 |        NULL |         NULL |       NULL | NULL    
            | NULL   |
|           24794858 | probe2transcript |        5253 |           NULL | 3685949    |                750 |        NULL |         NULL |       NULL | NULL    
            | NULL   |
|           24794859 | probe2transcript |        5253 |           NULL | 2607882    |                750 |        NULL |         NULL |       NULL | NULL    
            | NULL   |
|           24794860 | probe2transcript |        5253 |           NULL | 2405078    |                750 |        NULL |         NULL |       NULL | NULL    
            | NULL   |
|           24794861 | probe2transcript |        5253 |           NULL | 3815020    |                750 |        NULL |         NULL |       NULL | NULL    
            | NULL   |
|           24794862 | probe2transcript |        5253 |           NULL | 2344228    |                750 |        NULL |         NULL |       NULL | NULL    
            | NULL   |
|           24794863 | probe2transcript |        5253 |           NULL | 3927333    |                750 |        NULL |         NULL |       NULL | NULL    
            | NULL   |
|           24794864 | probe2transcript |        5253 |           NULL | 3490673    |                750 |        NULL |         NULL |       NULL | NULL    
            | NULL   |
|           24794865 | probe2transcript |        5253 |           NULL | 3379823    |                750 |        NULL |         NULL |       NULL | NULL    
            | NULL   |
+--------------------+------------------+-------------+----------------+------------+--------------------+-------------+--------------+------------+---------
------------+--------+
-- these are all exon probeset?


-- query are very slow on ensdb-archive, have tables been analysed/optimised?

mysql>  select of.* from oligo_feature of, oligo_probe op where op.probeset=3937834 and op.oligo_probe_id=of.oligo_probe_id;
+------------------+---------------+------------------+----------------+-------------------+------------+----------------+-------------+
| oligo_feature_id | seq_region_id | seq_region_start | seq_region_end | seq_region_strand | mismatches | oligo_probe_id | analysis_id |
+------------------+---------------+------------------+----------------+-------------------+------------+----------------+-------------+
|         35821107 |        226046 |         19662069 |       19662093 |                -1 |          0 |        7465409 |        5249 |
|         36779766 |        226046 |         19662067 |       19662091 |                -1 |          0 |       11373167 |        5249 |
|         39789796 |        226046 |         19662071 |       19662095 |                -1 |          0 |        9199101 |        5249 |
|         40293127 |        226046 |         19662066 |       19662090 |                -1 |          0 |       10878970 |        5249 |
+------------------+---------------+------------------+----------------+-------------------+------------+----------------+-------------+

-- THese are overlapping mappings by almost all of the 25mer probes

--Let's try another one.

mysql>  select of.* from oligo_feature of, oligo_probe op where op.probeset=3379823 and op.oligo_probe_id=of.oligo_probe_id;
+------------------+---------------+------------------+----------------+-------------------+------------+----------------+-------------+
| oligo_feature_id | seq_region_id | seq_region_start | seq_region_end | seq_region_strand | mismatches | oligo_probe_id | analysis_id |
+------------------+---------------+------------------+----------------+-------------------+------------+----------------+-------------+
|         36058894 |        226044 |         68625271 |       68625295 |                 1 |          0 |        9215078 |        5249 |
|         39052659 |        226044 |         68625261 |       68625285 |                 1 |          0 |       10721610 |        5249 |
|         41343148 |        226044 |         68625268 |       68625292 |                 1 |          0 |       11313286 |        5249 |
|         44551932 |        226044 |         68625266 |       68625290 |                 1 |          0 |       10158700 |        5249 |
+------------------+---------------+------------------+----------------+-------------------+------------+----------------+-------------+
4 rows in set (1 min 29.73 sec)

-- again overlapping mappings.

-- let's look at the sequences for these.

affy@bc-9-1-01 njst_homo_sapiens_core_51_36m>grep -A 1 'ProbeSetID=3379823' ~/data/affy/Human/HuEx-1_0-st-v2.probe.fa 
>probe:HuEx-1_0-st-v2:4418184;2183:1725; ProbeID=4418184; ProbeSetID=3379823; Assembly=build-34/hg16; Seqname=chr11; Start=68644045; Stop=68644069; Strand=-; Sense; category=main
TGGGTCCCACATCGGTGCTGGGTCC
>probe:HuEx-1_0-st-v2:4082681;2040:1594; ProbeID=4082681; ProbeSetID=3379823; Assembly=build-34/hg16; Seqname=chr11; Start=68644050; Stop=68644074; Strand=-; Sense; category=main
CCCACATCGGTGCTGGGTCCCAATG
>probe:HuEx-1_0-st-v2:4773675;1834:1864; ProbeID=4773675; ProbeSetID=3379823; Assembly=build-34/hg16; Seqname=chr11; Start=68644052; Stop=68644076; Strand=-; Sense; category=main
CACATCGGTGCTGGGTCCCAATGCC
>probe:HuEx-1_0-st-v2:2580085;2164:1007; ProbeID=2580085; ProbeSetID=3379823; Assembly=build-34/hg16; Seqname=chr11; Start=68644055; Stop=68644079; Strand=-; Sense; category=main
ATCGGTGCTGGGTCCCAATGCCCAC


-- These do overlap

-- Try a different probeset

affy@bc-9-1-01 njst_homo_sapiens_core_51_36m>grep -A 1 'ProbeSetID=3379822' ~/data/affy/Human/HuEx-1_0-st-v2.probe.fa 

>probe:HuEx-1_0-st-v2:3027620;1699:1182; ProbeID=3027620; ProbeSetID=3379822; Assembly=build-34/hg16; Seqname=chr11; Start=68643900; Stop=68643924; Strand=-; Sense; category=main
GAGGGGACACAGTCTTGGATATGAG
>probe:HuEx-1_0-st-v2:4031276;1835:1574; ProbeID=4031276; ProbeSetID=3379822; Assembly=build-34/hg16; Seqname=chr11; Start=68643902; Stop=68643926; Strand=-; Sense; category=main
GGGGACACAGTCTTGGATATGAGCC
>probe:HuEx-1_0-st-v2:2218776;1815:866; ProbeID=2218776; ProbeSetID=3379822; Assembly=build-34/hg16; Seqname=chr11; Start=68643903; Stop=68643927; Strand=-; Sense; category=main
GGGACACAGTCTTGGATATGAGCCC
>probe:HuEx-1_0-st-v2:4345854;1533:1697; ProbeID=4345854; ProbeSetID=3379822; Assembly=build-34/hg16; Seqname=chr11; Start=68643904; Stop=68643928; Strand=-; Sense; category=main
GGACACAGTCTTGGATATGAGCCCA

-- These are worse!! This can only be targetting one exon of an alt trans. Okay, but why is this not mapping? 


-- Let's take a look at the transcript neighbourhood.


select * from transcript where seq_region_start <= 68625295 and seq_region_end >=  68625261 and seq_region_id=226044 and seq_region_strand =1;

--nothing, so pad by 2kb
select * from transcript where seq_region_start <= 68627295 and seq_region_end >=  68623261 and seq_region_id=226044 and seq_region_strand =1;
--nothing, pad by 10kb
select * from transcript where seq_region_start <= 68635295 and seq_region_end >=  68614261 and seq_region_id=226044 and seq_region_strand =1;




-- No mappings at all for gene array

-- Need to implement counts for each array in QC logs!



-- clean up procedure

delete of from oligo_feature of, oligo_probe op, oligo_array oa where oa.name like '%_st_%' and oa.oligo_array_id=op.oligo_array_id and op.oligo_probe_id=of.oligo_probe_id;

delete op, oa from oligo_probe op, oligo_array oa where oa.name like '%_st_%' and oa.oligo_array_id=op.oligo_array_id;



-- remember to optimize/analyze all tables after update
