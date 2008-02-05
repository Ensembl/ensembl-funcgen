

#$1 is password

$EFG_SRC/scripts/load_pf.pl \
    -file           /lustre/work1/ensembl/prf1/ecs4/work6/prf1/Sanger/H3K4me3-GM06990.hg18.update.hits \
    -cell_type      GM06990\
    -feature_set_id 1\
    -feature_type   H3K4me3 \
    -cdbname        homo_sapiens_core_42_36d \
    -dbname         homo_sapiens_funcgen_44_36f \
    -dbhost         ens-genomics1 \
    -pass           $1
