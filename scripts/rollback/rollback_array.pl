=pod

=head1 NAME

ensembl-efg rollback_array.pl

=head1 SYNOPSIS

rollback_array.pl [options]

Options:

Mandatory
  -arrays|a        List of Array names (space separated)
  -chips|c         List of ArrayChip design names to restrict the delete to
  -mode|m          Rollback mode supported by Helper::rollback_ArrayChip, default if full rollback i.e. 'probe'
  -pass|p          DB password
  -dbname|n        DB namme
  -port            DB port
  -host|h          DB host
  -user|u          DB user name
  -species         Latin name e.g. homo_sapiens
  -force|f         Performs a full delete, removing all records upto the specified level.  This is to specifically force the rollback of probe2transcript entries, as these can be associated with more than one array i.e. you may be deleting annotations which also belong to another array.
  -help            Brief help message
  -man             Full documentation

=head1 DESCRIPTION

B<This program> removes all the imported probe, probe_feature and/or probe2transcript annotation data given set of Arrays or ArrayChips.

=cut



use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
$| =1;

my ($array_name, $pass, $force, @chips, $vendor, %achips);
my ($host, $dbname, $species, $mode, $port, $user);
my @tmp_args = @ARGV;

GetOptions (
			"arrays|a=s"      => \$array_name,
			"chips|c=s{,}"    => \@chips,
			"vendor|v=s"      => \$vendor,
			"mode|m=s"        => \$mode,
			"dbpass|p=s"      => \$pass,
			"dbport=s"        => \$port,
			"dbname|n=s"      => \$dbname,
			"dbhost|h=s"      => \$host,
			"species=s"       => \$species,
			"dbuser|u=s"      => \$user,
			"force|f"         => \$force,
			"help|?"          => sub { pos2usage(-exitval => 0,  
												 -verbose => 2, 
												 -message => "Params are:\t@tmp_args");
									 },
		   ) or pod2usage(
						 -exitval => 1,
						 -message => "Params are:\t@tmp_args"
						);

die('Must define a -user parameter')    if ! $user;
die('Must define a -port parameter')    if ! $port;
die('Must define a -dbname parameter')  if ! $dbname;
die('Must define a -dbhost parameter')  if ! $host;
die('Must define a -pass parameter')    if ! $pass;
die('Must define a -species parameter') if ! $species;
$force = 'force' if $force;

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													  -host => $host,
													  -dbname => $dbname,
													  -user => $user,
													  -pass => $pass,
						   						      -port => $port,
													  -species => $species,
													 );
$db->dbc->db_handle;#Test the DB connection
my $array       = $db->get_ArrayAdaptor->fetch_by_name_vendor($array_name, $vendor);
die ("Could not retrieve $vendor $array_name Array") if ! $array;
my $Helper = new Bio::EnsEMBL::Funcgen::Utils::Helper(
													  no_log => 1,#tees automatically with no_log
													 );
my $arraychip_a = $db->get_ArrayChipAdaptor();
my %acs;
map {$acs{$_->design_id} = $_} @{$array->get_ArrayChips};

#do chips belong to experiment?
if(@chips){
  
  foreach my $chip_name(@chips){
	
	if(! exists $acs{$chip_name}){
	  die("$chip_name is not a valid ArrayChip design_id for the $vendor $array_name Array");
	}

	$achips{$chip_name} = $acs{$chip_name};	
  }
}
else{
  %achips = %acs;
}


foreach my $ac(values(%achips)){
  $Helper->rollback_ArrayChip($ac, $mode, $force);
}

1;
