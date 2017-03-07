#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

e.g. perl update_transcript_xrefs.pl --out_dir $WORK_DIR --species $SPECIES --transcript_dbname $DNADB_NAME --transcript_host $DNADB_HOST --transcript_port $DNADB_PORT --transcript_user $DNADB_USER --xref_host $DB_HOST --xref_dbname $DB_NAME --xref_user $DB_USER --xref_pass $DB_PASS 


=head1 DESCRIPTION


=head1 OPTIONS

Note: [--parameter] Denotes optional parameter.

[--out_dir]  Default is current working directory.
--species    Latin name as used in DB name or meta table e.g. homo_sapiens

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

OR USING A REGISTRY:
[--reg_verbose]    Turns on verbose output when loading the registry

--reg_host
--reg_port
--reg_user
[--reg_pass]
or
[--reg_file]


Other options:
--tee                Tees output to STDOUT
--log_file           Default is $out_dir/xref_dbname_update_transcript_xrefs.log. Or if --out_dir is
                     not specified, the default Helper log directory.
--help               Prints this POD documentation and exits


=head1 EXAMPLE


=head1 SEE ALSO

ensembl-funcgen/scripts/environments/arrays.env
ensembl-funcgen/scripts/probe2transcript.pl

=cut

#To do

# 1. Make this do the UTR calculations once, so it's not done redundantly in probe2transcript.
#    Move this code to a separate module, so we can keep functionatlity in both scripts.
#    Will need to bsub this, and write/flow output somewhere.

# 2. Remove Helper and print to std filename directly. Currently log_file is never used as Helper
#    is initialised too early. Hence some output going may go to STDOUT and some to the default log dir

use strict;

use Pod::Usage;
use Getopt::Long;
use File::Temp qw/tempfile/;
use IO::Handle;

use Bio::EnsEMBL::DBSQL::Driver; # For over-riding connect_params method
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_backtick_cmd );

$| = 1; # auto flush stdout

# Helper params
my $Helper = new Bio::EnsEMBL::Funcgen::Utils::Helper;  
# This needs initisalising later, see To do 2 above
$main::_log_file = undef;
$main::_tee      = 0;    
my $debug = 0;

# print "$0 @ARGV\n"; 

main();

sub main {
	my $options = get_options();

	$Helper->log("Setting global constant variables", 0, 'append_date');
	my ($transcript_db, $xref_db) = @{get_databases($options)};
	my $schema_build = $xref_db->_get_schema_build($transcript_db);
	my $transc_edb_name = "$options->{species}_core_Transcript";
	$Helper->log("Getting or creating external DB id\n");
	my $transc_edb_id = get_external_db_id($xref_db, $options->{species}, $transc_edb_name, $schema_build);
	$Helper->log("Extracting transcripts from core database\n");
	my $transcripts = get_transcripts($transcript_db, $options->{test_slice}, $options->{test_transcript_sid});
	$Helper->log("Count unknown transcripts\n");
	my $new_transcripts = get_unknown_transcript_list($transcripts, $xref_db, $transc_edb_id);
	if (! scalar @$new_transcripts) {
		$Helper->log("No transcripts need adding\n", 0, 1);
	} else {
		$Helper->log("Found ".scalar @$new_transcripts ." transcripts to add\n", 0, 1);
		load_new_transcript_xrefs ($new_transcripts, $xref_db, $transc_edb_id, $options);
	}
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
    'transcript_host=s'        => \$options->{transcript_host},
    'transcript_user=s'        => \$options->{transcript_user},
    'transcript_port=i'        => \$options->{transcript_port},
    'transcript_pass=s'        => \$options->{transcript_pass},
    'transcript_dbname=s'      => \$options->{transcript_dbname},
    'transcript_species_id=i'  => \$options->{transcript_species_id},
    'transcript_multi_species' => \$options->{transcript_multi_species},
    'xref_host=s'              => \$options->{xref_host},
    'xref_user=s'              => \$options->{xref_user},
    'xref_port=i'              => \$options->{xref_port},
    'xref_pass=s'              => \$options->{xref_pass},
    'xref_dbname=s'            => \$options->{xref_dbname},
    'xref_species_id=i'        => \$options->{xref_species_id},
    'xref_multi_species'       => \$options->{xref_multi_species},
    'reg_file=s'               => \$options->{reg_file},
    'reg_host=s'               => \$options->{reg_host},
    'reg_user=s'               => \$options->{reg_user},
    'reg_pass=s'               => \$options->{reg_pass},
    'reg_port=i'               => \$options->{reg_port},
    'reg_verbose'              => \$options->{reg_verbose},
    'species=s'                => \$options->{species},
    'out_dir=s'                => \$options->{out_dir},
    # Helper params
    'tee'                    => \$main::_tee,
    'log_file'               => \$main::_log_file,
    'help'                   => sub { pos2usage(-exitval => 0, -message => "Params are:\t@tmp_args"); }
	) or pod2usage(
		-exitval => 1,
		-message => "Params are:\t@tmp_args"
	);

#use Data::Dumper; print Dumper($options); exit;
  # This is not working correctly, and is currently writing log file to default log dir, not workdir
  # no_log should be set in caller

	#$options->{filename} ||= "$options->{xref_dbname}_$options->{log_type}_probe2transcript";


  if(defined $options->{out_dir}){

    if(! -d $options->{out_dir}){
      die("Parameter --out_dir is not a valid directory:\t".$options->{out_dir});
    }

    #append a slash is not already present
    if($options->{out_dir} !~ /\/$/o){
      $options->{out_dir} .= '/';
    }
  }
  else{ $options->{out_dir} = ''; }  # Avoid under concat warning

  # This is currently useless, see To do 2 above
  if(! defined $main::_log_file){
    $main::_log_file = $options->{out_dir}.$options->{xref_dbname}.'_update_transcript_xrefs.log';
  }

	$options->{hostname} = run_backtick_cmd('hostname');
	$Helper->log_header("Running on $0 on: ".$options->{hostname}, 0, 'append_date');
	$Helper->log("Params are:\t@tmp_args");

	if(! $options->{species}) {
		die('Must provide a -species');
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
			-host   => $options->{xref_host},
			-port   => $options->{xref_port},
			-user   => $options->{xref_user},
			-pass   => $options->{xref_pass},
			-dbname => $options->{xref_dbname},
			-species => $options->{species},
			-multispecies_db => $options->{xref_multi_species},
			-species_id => $options->{xref_species_id}
		);
	}

  # Test the DBs here before starting
 
	$transcript_db->dbc->db_handle;

  # Inject/over-write Driver::connect_params method to handle mysql_local_infile config
  # This has the advantage of maintaining the use of the Ensembl DBAdaptors, so we don't
  # have to write any more exception handling or other wrappers
  # The cost here is the risk of this (very stable method) getting out of sync
  no strict 'refs';
  
  *{'Bio::EnsEMBL::DBSQL::Driver::connect_params'} = sub { 
    my $self = shift;
    my $conn = shift;

    my $dbname = $conn->dbname();
    my $dbparam = ($dbname) ? "database=${dbname};" : q{};

    my $dsn = sprintf( "DBI:%s:%shost=%s;port=%s",
     $conn->driver(), $dbparam,
     $conn->host(),   $conn->port() );

    my $attrs = { 'RaiseError' => 1 };

    if ((exists $conn->{local_infile}) && $conn->{local_infile}){
      #$dsn .= ';mysql_local_infile=1';
      $attrs->{mysql_local_infile} = 1;
    }

    if ( $conn->{'disconnect_when_inactive'} ) {
      $conn->{'count'}++;
      if ( $conn->{'count'} > 1000 ) {
        sleep 1;
        $conn->{'count'} = 0;
      }
    }

    return {
      dsn        => $dsn,
      username   => $conn->username(),
      password   => $conn->password(),
      attributes => $attrs,
    };
  };
     
  use strict;

  $xref_db->{local_infile} = 1;
	$xref_db->dbc->db_handle;

  # Allow automatic reconnection
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
#
# Returns: ID number of external id in external_db table (xref system)

sub get_external_db_id {
	my($xref_db, $species, $edb_name, $schema_build) = @_;
#Check for external_db records for species DBs
	my ($transc_edb_display_name, $edb_display);
	$transc_edb_display_name = "EnsemblTranscript";
	$edb_display             = $transc_edb_display_name;

	# Look for pre-existing
	my $sql = "
	SELECT 
		external_db_id
	FROM
		external_db
	WHERE	
		db_name='$edb_name'
		AND db_release = '$schema_build'";
	my @ids = @{$xref_db->dbc->db_handle->selectall_arrayref($sql)};
  	if (scalar @ids > 0) {
		return $ids[0][0];
  	}

	$sql = 'INSERT into external_db(db_name, db_release, status, dbprimary_acc_linkable, priority, db_display_name, type) values('.
	"'${edb_name}', '${schema_build}', 'KNOWNXREF', 1, 5, '$edb_display', 'MISC')";

	$Helper->log("Importing external_db using: $sql");
	$xref_db->dbc->do($sql);
	return $xref_db->db_handle->last_insert_id();
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
	my ($new_transcripts, $xref_db, $transc_edb_id, $options) = @_;
	my ($fh, $filename) = tempfile();
	foreach my $transcript (@$new_transcripts) {
		print $fh join("\t", ('\N', $transc_edb_id, $transcript->stable_id, $transcript->display_id, $transcript->version, '\N', 'MISC','TRANSCRIPT'))."\n";
	}

	chmod 0644, $filename; # Just in case default means mysql can't read it
	$fh->autoflush;
	my $cmd = "mysql --local-infile -u $options->{xref_user} -h $options->{xref_host} -D $options->{xref_dbname} -e 'LOAD DATA LOCAL INFILE \"$filename\" INTO TABLE xref'";
	if (defined $options->{xref_port}) {
	$cmd .= " -P $options->{xref_port}";
	}
	if (defined $options->{xref_pass}) {
	$cmd .= " -p$options->{xref_pass}";
	}
	run($cmd);
	close $fh;
	unlink $filename;
}

########################################################
## System calls 
## Wrapper function for system calls
## Params:
## - Command line
## Actions:
## - Runs command, prints out error in case of failure
########################################################

sub run {
  my ($cmd) = @_;
  $Helper->log("Running $cmd\n", 0, 'append_date');
  my $exit_code = system($cmd);
  if ($exit_code != 0) {
    die("Failure when running command\n$cmd\n")
  }
}

