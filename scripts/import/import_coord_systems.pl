#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

use warnings;
use strict;
use Carp;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;
use Pod::Usage;

my $opts = {};
GetOptions($opts, qw(seq_regions|sr help|? verbose man)) or pod2usage(1);
pod2usage(1) if $opts->{help};
pod2usage(-exitstatus => 0, -verbose => 2) if $opts->{man};

#Currently we ask for this config for the adaptor information but a far better
#way of solving the problem would be to use command line options ... that said
#ProbeAlign is an integral part of the EFG mapping pipeline so for the moment
#it should be okay
use Bio::EnsEMBL::Analysis::Config::ProbeAlign;
my $d_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%{$PROBE_CONFIG->{DEFAULT}->{DNADB}});
my $f_dba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(%{$PROBE_CONFIG->{DEFAULT}->{OUTDB}}, -dnadb => $d_dba);

if($opts->{verbose}) {
	use Data::Dumper;
	print Dumper ($f_dba), "\n";
	print Dumper ($d_dba), "\n";
}

my $sa = $d_dba->get_SliceAdaptor();
my $pfa = $f_dba->get_ProbeFeatureAdaptor();
my $fcsa = $f_dba->get_FGCoordSystemAdaptor();

#Have to fetch slices then get toplevel coord systems
my $slices = $sa->fetch_all('toplevel');

foreach my $cs (sort { $a->rank() <=> $b->rank() } unique_cs($slices)) {
	store_cs($cs);
}

foreach my $slice (@{$slices}) {
	store_slice($slice);
}

sub unique_cs {
	my ($slices) = @_;
	my %seen;
	my @coord_systems =
		grep {
			my $cs = $_;
			my $key = join('_|_', $cs->name(), $cs->version());
			if($seen{$key}) {
				0;
			}
			else {
				$seen{$key} = 1;
				1;
			}
		}
		
		map {
			my $slice = $_;
			$slice->coord_system();
		}
		@{$slices};
	
	return @coord_systems;
}

sub fg_cs {
	my ($cs) = @_;
	return $fcsa->fetch_by_name($cs->name(), $cs->version());
}

sub store_cs {
	my ($cs) = @_;
	my $fg_cs = fg_cs($cs);
	if($fg_cs) {
		print STDERR 'Coord-system ', $cs->name(), ' already exists', "\n" if $opts->{verbose};
	}
	else {
		print STDERR 'Coord-system ', $cs->name(), ' needs to be inserted', "\n" if $opts->{verbose};
		$fcsa->validate_and_store_coord_system($cs);
	}
}

sub store_slice {
	my($slice) = @_;
	
	if($opts->{verbose}) {
		print STDERR 'Processing slice '.$slice->name()."\n";
	}
	
	$pfa->build_seq_region_cache($slice);
	
	my $seq_region_id = $pfa->get_seq_region_id_by_Slice($slice, undef, 1);

  if(! $seq_region_id) {
  	my $fg_cs = fg_cs($slice->coord_system());
		#check whether we have an equivalent seq_region_id
		$seq_region_id = $pfa->get_seq_region_id_by_Slice($slice, $fg_cs);
		my $schema_build = $f_dba->_get_schema_build($slice->adaptor->db());
		my $sql;
		my @args = ($slice->seq_region_name(), $fg_cs->dbID(), $slice->get_seq_region_id(), $schema_build);

		#Add to comparable seq_region		
		if($seq_region_id) {
			$sql = 'insert into seq_region(seq_region_id, name, coord_system_id, core_seq_region_id, schema_build) values (?,?,?,?,?)';
			unshift(@args, $seq_region_id);
		}
		#No compararble seq_region
		else{
			$sql = 'insert into seq_region(name, coord_system_id, core_seq_region_id, schema_build) values (?,?,?,?)';			
		}

		my $sth = $pfa->prepare($sql);
	
		#Need to eval this
		eval{$sth->execute(@args);};
	
		if(!$@){
		  $seq_region_id =  $sth->{'mysql_insertid'};
		  print STDERR "Just inserted new seq region with ID ${seq_region_id}\n" if $opts->{verbose};;
		}
		else {
			croak("Could not insert new sequence region: $@");
		}
  }
  else {
  	print STDERR 'Sequence region already represented', "\n"  if $opts->{verbose};
  }
}

__END__
=pod

=head1 NAME

	import_coord_systems.pl
	
=head1 SYNOPSYS

	./import_coord_systems.pl [-seq_regions] [-v] [-h] [-m]

=head1 DESCRIPTION

This class is used to initalise the coordinate systems for a given source database. The
code reuses the configuration used for the ProbeAlign module which gives us access
to all required database information. If you are not using this module then well there's going
to be some problems.

This script will then allow us to avoid the race condition which can appear during EFG mapping
when we have coordinate systems which are not chromosomes but are toplevel. The effect of this
occuring in an unmanaged manner is that database contraints will cause the pipeline to
restart. This script allows us to register the coordinate systems in a single process & therefore
far more managed than the old system.

The script can also import sequence regions by using a copy of the sequence
region insertion logic. Whilst this is not ideal it is required for genomes
where regions may be intentionally missed by a mapping procedure and then
later assumed they have been imported into the functional genomics schema.

=head1 OPTIONS

=over 8

=item B<-help>
 
Prints this message

=item B<-man>

Prints the manual version

=item B<-verbose>

Prints a bit of extra information about the program you are running. This is feedback about the
inserted coordindate systems & connection settings using Data::Dumper. This is crude feedback at
best but the script is temporary.

=item B<-seq_regions>

If specified this will attempt to write a feature to the probe feature table
and then delete it. The result of this is to trigger the insertion of 
sequence regions which may not have been persisted due to no probes being
mapped to them. This causes problems during probe2transcript.pl where all
slices from a core database are expected to be found in the target database.

=back

=cut
