#!/bin/sh

#. ~/src/ensembl-efg/scripts/.efg

USER=$1
PASS=$2

if [ ! $USER ] || [ ! $PASS ]; then
 echo 'Usage: run_stable_id_mapping.sh $USER $PASS [ -other args ]'; 
 exit 1;
fi

shift
shift




time perl -w $EFG_SRC/scripts/regulatory_build/stable_id_mapper.pl\
	-ndbname homo_sapiens_funcgen_55_37 \
	-nhost ens-genomics1 \
	-nuser $USER \
   	-npass $PASS \
	-dnadb_name homo_sapiens_core_55_37\
	-clobber \
	-out_dir $HOME/data/efg/reg_build/homo_sapiens_funcgen_55_37/stable_id_mapping\
	-old_fset_name RegulatoryFeatures_v5 \
	-new_fset_name RegulatoryFeatures \
	-old_assembly NCBI36\
	-new_assembly GRCh37\
	$@
