#!/usr/local/bin/perl

use warnings;
use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;
use Pod::Usage;

my $opts = {};
GetOptions($opts, qw(help|? verbose man)) or pod2usage(1);
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
my $fcsa = $f_dba->get_FGCoordSystemAdaptor();

#Have to fetch slices then get toplevel coord systems
my $slices = $sa->fetch_all('toplevel');
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

foreach my $cs (sort { $a->rank() <=> $b->rank() } @coord_systems) {
	my $fg_cs = $fcsa->fetch_by_name($cs->name(), $cs->version());
	if($fg_cs) {
		print STDERR 'Coord-system ', $cs->name(), ' already exists', "\n" if $opts->{verbose};
	}
	else {
		print STDERR 'Coord-system ', $cs->name(), ' needs to be inserted', "\n" if $opts->{verbose};
		$fcsa->validate_and_store_coord_system($cs);
	}
}

__END__
=pod

=head1 NAME

	import_coord_systems.pl
	
=head1 SYNOPSYS

	./import_coord_systems.pl [-v] [-h] [-m]

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
 
=back

=cut