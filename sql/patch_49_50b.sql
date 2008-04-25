--insert into feature_type values(NULL, 'TSS 2.5KB', 'Pseudo', 'Feature is within +/-2.5KB of the Transcript Start Site');
i--nsert into feature_type values(NULL, 'TTS 2.5KB', 'Pseudo', 'Feature is within +/-2.5KB of the Transcript Termination Site');


--Human only data patch
--This may change if we do another build for v50
--update meta set meta_value='9,8,23,24,25,26,11,27,12,6,5,4,28,29,13,30,31,14,32,33,34,10,35,36,8,9,11,12,4,13,14,10,178400,178401' where meta_key='regulatory_string_feature_type_id';


--remove MT reg features
-- delete from regulatory_feature where seq_region_id =137;
--End human only data patch
