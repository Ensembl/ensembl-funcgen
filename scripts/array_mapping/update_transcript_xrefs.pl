#
#/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITH$OUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

update_transcript_xrefs.pl

=head1 SYNOPSIS

This script performs probe(set) to transcript mapping based on a few simple parameters. Overlap analysis 
of ProbeFeatures is performed and annotations are stored as xrefs for individual ProbeFeatures, Probes or 
ProbeSets as a whole. Any probe(set)s which fail the mapping procedure are by default stored in the 
UnmappedObject tables and a logfile is also written.

e.g. perl update_xref_transcripts.pl --species $SPECIES --transcript_dbname $DNADB_NAME --transcript_host $DNADB_HOST --transcript_port $DNADB_PORT --transcript_user $DNADB_USER --xref_host $DB_HOST --xref_dbname $DB_NAME --xref_user $DB_USER --xref_pass $DB_PASS 


=head1 DESCRIPTION


=head1 OPTIONS

Mandatory
-species    Latin name as used in DB name or meta table e.g. homo_sapiens

-reg_verbose    Turns on verbose output when loading the registry

-reg_host
-reg_port
-reg_user
-reg_pass

or

-reg_file

or

-transcript_host          Mandatory
-transcript_user          Mandatory
-transcript_port    
-transcript_pass          Mandatory
-transcript_dbname        Mandatory
-transcript_multi_species
-transcript_species_id

-xref_host          Mandatory
-xref_user          Mandatory
-xref_port
-xref_pass          Mandatory
-xref_dbname        Mandatory
-xref_multi_species
-xref_species_id

Other options:
-import_edb         Automatically imports the external_db record if not present
-tee                Tees output to STD$OUT
-filename           Sets name to be used in output and logfile, default is xref_dbname_probe2transcript.log|out
-help               Prints this POD documentation and exits


=head1 EXAMPLE


=head1 SEE ALSO

ensembl-funcgen/scripts/environments/arrays.env

=cut

#To do

# 1. Reimpliment validate arrays, see old script?
# 2. Add unannotated UTR clipping dependant on nearest neighbour
# 3. Extend UTRs to default length is they are less than defaults, so long as they don't overlap neighbour, 
#    then use annotated if present or clip to neighbour start/end if not, also accounting for default UTRs 
#    in the neighbour.
# 4. Separate UTR multipliers for 3' and 5'?
# 5. Implement incremental update from list of stable IDs. Consider unmapped probe changes etc. 
# 6. Parallelise by probeset chunks, can't do this by chromosome slices as we need to know genomewide 
#    counts for a given probeset. Calc UTRs then submit chunks jobs to farm
#    Chunk by retrieving all probesets and sorting an array of probeset names, then splice the array 
#    according to the number of chunks. We're still going to have retrieve all the transcripts and retrieve 
#    all probes for each, so we are really not gaining anything!! The only gain we can make is by chunking 
#    by slice, but then we need to know how many times something has mapped. Can we do some clean up afterwards? 
#    Let's add a clean up mode which simply deletes all probe sets which map too many times. We would need to 
#    ignore this threshold as we were mapping!!! So we don't delete and then mess up the counts for post run 
#    clean up.
# 7. There is no reason to have separate probe and xref DBs???
# 8. Validate array format against arrays specified? May want to just use an array format as a template???
# 9. Add mismatch filter for ProbeTranscriptAlign xrefs as match rules can differ between alignment and 
#    annotation
# 10.Handle ProbeAlign mismatch vs overlap mis match. Currently the overlap calculation is naive to the 
#    presence of alignment mis-matches.  Which means there is a possiblity of including probes with a total 
#    sequence mismatch of (align mismatch + overlap mismatch). This has always been the case.
# 11.Move ProbeAlign unmapped object storage to write_output, then this will not get written in test mode and 
#    we won't get duplication should the job fail halfway through. This is because hceck existing only check oxs, not uos.
# 12.Enable probesets to have different sizes on different arrays, see notes in cache_arrays_per_object
# 13.Collect warning into summary repoprt to list at very end.
# 14 Reduce max_transcripts as this is never being hit due to alignment threshold
# 15 Why can't we omit -arrays if we have -format?
# 16 Add UTR only overlap  in range registry.
# 17 Check for ProbeFeature xrefs and UOs in check_existing_and_exit?
# 18 PostAlign/PreXref processing
#    Remove duplicated ProbeFeatures(from ProbeTranscriptAlign) and redirect Xrefs
#    Being careful to make sure cigarlines are valid for both.
#    Remove ProbeTranscriptAlign ProbeFeaturess which have been called promiscuous 
#    by ProbeAlign, and update to promiscuous if sum of ProbeAlign and 
#    ProbeTranscriptAlign features render a Probe promiscuous

#Ensembl Genomes stuff
# TEST Registry usage required as species will come from same DB
# In which case we need to take a species param for each of the transcript, array and xref DBs
# Or can we force delete to be species specific? We would need to do this anyway to support updating of species asynchronously
# We probably need to think about this for the 1st stage too, 
# but will be easy as we just need to dump the correct top level sequence
# Validate species against registry alias and use this to generate species_core_Gene DB rather than ensembl_core_Gene
# patch other efg DBs and alter External parsers accordingly.
# Can't rely on Registry as species aliases may not be present or loaded

# Issues
# Cannot account for running non-linked arrays which may use the same probe/set name.  This may cause failure if the probeset sizes are different. Xrefs and counts should be unaffected as we base these on the probe_set_ids not the names. This is not really an issue as unlinked arrays should not be run together
# Cannot currently handle probesets with different sizes between arrays, defaults to lowest probeset size to be permissive. See todo 12.

use strict;

use Pod::Usage;
use Getopt::Long;
use File::Temp qw/tempfile/;
use IO::Handle;

use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::UnmappedObject;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Mapper::RangeRegistry;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (median mean get_date);

$| = 1; # auto flush stdout

#Helper params
$main::_log_file = undef;
$main::_tee      = 0;    
our $Helper = new Bio::EnsEMBL::Funcgen::Utils::Helper;
my $debug = 0;

main();

sub main {
	my $options = get_options();

	$Helper->log("Setting global constant variables", 0, 'append_date');
	my ($transcript_db, $xref_db) = @{get_databases($options)};
	my $schema_build = $xref_db->_get_schema_build($transcript_db);
	my $transc_edb_name = "$options->{species}_core_Transcript";
	my $transc_edb_id = get_external_db_id($xref_db, $options->{species}, $transc_edb_name, $schema_build, $options->{import_edb});
	my $transcripts = get_transcripts($transcript_db, $options->{test_slice}, $options->{test_transcript_sid});
	my $new_transcripts = get_unknown_transcript_list($transcripts, $xref_db, $transc_edb_id);
	if (! scalar @$new_transcripts) {
		$Helper->log("No transcripts need adding\n", 0, 1);
	} else {
		$Helper->log("Found ".scalar @$new_transcripts ." transcripts to add\n", 0, 1);
		load_new_transcript_xrefs ($new_transcripts, $xref_db, $transc_edb_id);
	}
}

# ----------------------------------------------------------------------

sub usage {

	print << "EOF";
Updates transcript xref ids.

perl $0 {options}

Options ([..] indicates optional):

READING TRANSCRIPTS:
--transcript_host            The database server to read transcripts from.
[--transcript_port]          The port to use for reading transcripts. Defaults to 3306.
--transcript_user            Database username for reading transcripts.
--transcript_pass            Password for transcript_user, if required.
--transcript_dbname          Database name to read transcripts from.
[--transcript_multi_species] Indicates that the transcript database is multi-species
[--transcript_species_id]    Species ID to use top access multi-species data

WRITING XREFS:
--xref_host            The database server to write xrefs to.
[--xref_port]          The port to use for writing xrefs.. Defaults to 3306.
--xref_user            Database username for xrefs. Must allow writing.
--xref_pass            Password for xref_user, if required.
--xref_dbname          Database name to write xrefs to.
[--xref_multi_species] Indicates that the transcript database is multi-species
[--xref_species_id]    Species ID to use top access multi-species data
EOF

	exit(0);

}

# ----------------------------------------------------------------------

sub get_options {
	my $options = {};

# Default options
	$options->{reg_verbose} = 0;
	$options->{transcript_port} = 3306; 
	$options->{xref_port} = 3306;

	my @tmp_args = @ARGV;

	GetOptions(
		'transcript_host=s'      => \$options->{transcript_host},
	'transcript_user=s'      => \$options->{transcript_user},
	'transcript_port=i'      => \$options->{transcript_port},
	'transcript_pass=s'      => \$options->{transcript_pass},
	'transcript_dbname=s'    => \$options->{transcript_dbname},
	'transcript_species_id=i' => \$options->{transcript_species_id},
	'transcript_multi_species' => \$options->{transcript_multi_species},
	'xref_host=s'            => \$options->{xref_host},
	'xref_user=s'            => \$options->{xref_user},
	'xref_port=i'            => \$options->{xref_port},
	'xref_pass=s'            => \$options->{xref_pass},
	'xref_dbname=s'          => \$options->{xref_dbname},
	'xref_species_id=i' => \$options->{xref_species_id},
	'xref_multi_species' => \$options->{xref_multi_species},
	'reg_file=s'             => \$options->{reg_file},
	'reg_host=s'             => \$options->{reg_host},
	'reg_user=s'             => \$options->{reg_user},
	'reg_pass=s'             => \$options->{reg_pass},
	'reg_port=i'             => \$options->{reg_port},
	'reg_verbose'            => \$options->{reg_verbose},
	'species=s'              => \$options->{species},
#Helper params
	'tee'                    => \$main::_tee,
	'filename'               => \$main::_log_file,
#add a reduced log to minimize memory usage?
	'help'                   => sub { pos2usage(-exitval => 0, -message => "Params are:\t@tmp_args"); }
	) or pod2usage(
		-exitval => 1,
		-message => "Params are:\t@tmp_args"
	);

#Set log type so we are no over writing to the same files for different 
#format, or custom formats
	$options->{log_type} = $options->{format} || $$;
	$options->{filename} ||= "$options->{xref_dbname}_$options->{log_type}_probe2transcript";
	$main::_log_file ||=  "./$options->{filename}.log";
	$options->{hostname} = `hostname`;
	chomp($options->{hostname});
	$Helper->log_header('Running on probe2transcript.pl on: '.$options->{hostname}, 0, 'append_date');
	$Helper->log("Params are:\t@tmp_args");

	if(! $options->{species}) {
		die("Must provide a -species");
	}

	if($options->{reg_host} && ! ($options->{reg_user} && $options->{reg_pass})) {
		die('Must provide at least a -reg_user -reg_pass (optional -reg_port -reg_verbose) if loading from db');
	}

	if (!$options->{reg_host}) {
		if (!$options->{transcript_user} || !$options->{transcript_dbname} || !$options->{transcript_host}) {
			die("You must specify a -transcript_user -transcript_dbname -transcript_host\n");
		} 
		if(!$options->{xref_user} || !$options->{xref_dbname} || !$options->{xref_host}) {
			die("You must specify a -xref_user -xref_dbname and -xref_host\n");
		}
	}

	return $options;
}

# ----------------------------------------------------------------------

sub get_databases {
	my ($options) = @_;
	my ($transcript_db, $xref_db);
	if($options->{reg_file} || $options->{reg_host}) {
		my $reg = 'Bio::EnsEMBL::Registry';
		if($options->{reg_file}) {
			$Helper->log("Loading registry from:\t".$options->{reg_file});
			$reg->load_all($options->{reg_file}, $options->{reg_verbose});
		} else{
			$reg->load_registry_from_db(
				-host    => $options->{reg_host},
				-port    => $options->{reg_port} || 3306,
				-user    => $options->{reg_user},
				-pass    => $options->{reg_pass},
				-verbose => $options->{reg_verbose},
			);
		}        

		$transcript_db = $reg->get_DBAdaptor($options->{species}, 'Core');
		$xref_db       = $reg->get_DBAdaptor($options->{species}, 'Funcgen');
	} else{
		$transcript_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
			-host    => $options->{transcript_host},
			-port    => $options->{transcript_port},
			-user    => $options->{transcript_user},
			-pass    => $options->{transcript_pass},
			-dbname  => $options->{transcript_dbname},
			-species => $options->{species},
			-multispecies_db => $options->{transcript_multi_species},
			-species_id => $options->{transcript_species_id}
		);

		$xref_db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
			-host => $options->{xref_host},
			-port   => $options->{xref_port},
			-user   => $options->{xref_user},
			-pass   => $options->{xref_pass},
			-dbname => $options->{xref_dbname},
			-species => $options->{species},
			-multispecies_db => $options->{xref_multi_species},
			-species_id => $options->{xref_species_id}
		);
	}

#Test the DBs here before starting
	$transcript_db->dbc->db_handle;
	$xref_db->dbc->db_handle;

#Allow automatic reconnection
	$transcript_db->dbc->disconnect_if_idle(1);
	$xref_db->dbc->disconnect_if_idle(1);

	return [$transcript_db, $xref_db];
}

# ----------------------------------------------------------------------
# Params:
# - xref_db: Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor object
# - species: string name
# - edb_name: MySQL database name 
# - schema build: string
# - import_edb: boolean
#
# Returns: ID number of external id in external_db table (xref system)

sub get_external_db_id {
	my($xref_db, $species, $edb_name, $schema_build, $import_edb) = @_;
#Check for external_db records for species DBs
	my ($transc_edb_display_name, $edb_display);
	$transc_edb_display_name = "EnsemblTranscript";
	$edb_display             = $transc_edb_display_name;

	my $sql = "SELECT external_db_id, db_release from external_db where db_name='$edb_name'";
	my @versions = @{$xref_db->dbc->db_handle->selectall_arrayref($sql)};
	my @tmp;

	foreach my $row(@versions) {
		my ($edb_id, $version) = @$row;
		push @tmp, $version;

		if($schema_build eq $version) { 
			return $edb_id;
		}
	}

	$sql = 'INSERT into external_db(db_name, db_release, status, dbprimary_acc_linkable, priority, db_display_name, type) values('.
	"'${edb_name}', '${schema_build}', 'KNOWNXREF', 1, 5, '$edb_display', 'MISC')";

	if(! $import_edb) {
		die("Could not find current external_db $edb_name $schema_build from available versions:\t @tmp\nMaybe you have mis-spelt the -trans-species or you may need to manually add the external_db to the table and master file:\n\n$sql\n\n");
	} 

	$Helper->log("Importing external_db using: $sql");
	$xref_db->dbc->db_handle->do($sql);
	return $xref_db->last_insert_id();
}


# ----------------------------------------------------------------------
# Params:
# - transcript_db: Bio::EnsEMBL::DBSQL::DBAdaptor object
# - test_slice: string describing a slice (see Bio::EnsEMBL::DBSQL::SliceAdaptor::fetch_by_name)
# - test_transcript_sid: Ensembl transcript stable identifier
#
# Returns: Array ref of Bio::EnsEMBL::Transcript objects 

sub get_transcripts {
	my ($transcript_db, $test_slice, $test_transcript_sid) = @_;
	my $transcript_adaptor = $transcript_db->get_TranscriptAdaptor();
	my $slice_adaptor = $transcript_db->get_SliceAdaptor();
	my $transcripts;
	if($test_slice) {
		$Helper->log("Running in test mode with slice:\t$test_slice\n".
			"WARNING:\tPromiscuous probesets will not be caught! Calculated UTRs will be wrong!");
#Need to add better text here when we have implemented parallel runs
		my $slice = $slice_adaptor->fetch_by_name($test_slice);
		if (! defined $slice) {
			die("Could not get slice from the DB:\t$slice");
		}
		$transcripts = $transcript_adaptor->fetch_all_by_Slice($slice);
	} elsif($test_transcript_sid) {
		$Helper->log("Running test mode with transcript:\t$test_transcript_sid\n".
			"WARNING:\tPromiscuous probeset will not be caught!\n".
			"WARNING:\t--calc_utrs will not work");

		$transcripts = $transcript_adaptor->fetch_by_stable_id($test_transcript_sid);
	} else{
		$transcripts = $transcript_adaptor->fetch_all();
	}

	if (scalar @$transcripts == 0) {
		die('Could not find any transcripts');
	}
	$Helper->log("Identified ".scalar(@$transcripts)." transcripts for probe mapping");

	my @final = sort {$a->stable_id cmp $b->stable_id} @$transcripts;
	return \@final;
}


# Get internal xref_id for the transcripts in the xref table
sub get_transcript_xref_ids {
	my ($xref_db, $transc_edb_id) = @_;
	my $hash = {};
	my $sql = "SELECT dbprimary_acc, xref_id FROM xref WHERE external_db_id = $transc_edb_id;";
	foreach my $row (@{$xref_db->dbc->db_handle->selectall_arrayref($sql)}) {
		$hash->{$row->[0]} = $row->[1];
	}
	return $hash;
}

sub get_unknown_transcript_list {
	my ($transcripts, $xref_db, $transc_edb_id) = @_;
	my $transcript_xref_id = get_transcript_xref_ids($xref_db, $transc_edb_id);
	my %seen = ();
	my @new_transcripts = grep {!$seen{$_->stable_id}++ && !$transcript_xref_id->{$_->stable_id}} @$transcripts;
	return \@new_transcripts;
}

# Load unknown transcripts into the xref table
sub load_new_transcript_xrefs {
	my ($new_transcripts, $xref_db, $transc_edb_id) = @_;
	my ($fh, $filename) = tempfile(DIR=>'/nfs/users/nfs_d/dz1/lustre2/');
	foreach my $transcript (@$new_transcripts) {
		print $fh join("\t", ('/N', $transc_edb_id, $transcript->stable_id, $transcript->display_id, $transcript->version, '\N', 'MISC','TRANSCRIPT'))."\n";
	}
	chmod 0644, $filename;
	$fh->autoflush;
	my $sql = "LOAD DATA LOCAL INFILE \"$filename\" INTO TABLE xref";
	$xref_db->dbc->db_handle->do($sql);
	close $fh;
}

