#!/bin/sh

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


USER=$1
shift



dbname=homo_sapiens_funcgen_77_38
#dbname=mus_musculus_funcgen_76_38
dbhost=ens-staging1
dnadb_host=ens-staging1

dnadb_user=$USER

nfs_root='/nfs/ensnfs-dev/staging';
ftp_root='/lustre/scratch109/ensembl/funcgen/output/ftp'


if [ ! $USER ]; then
	echo "Must provide a user argument"; exit; 
fi


job_cmd="perl -w $EFG_SRC/scripts/release/create_ftp_dir_links.pl\
  -nfs_root $nfs_root\
  -ftp_root $ftp_root\
	-dbhost $dbhost\
	-dbuser $USER \
	-dbname $dbname\
	-dnadb_host $dnadb_host\
	-dnadb_user $dnadb_user $@"

$job_cmd







