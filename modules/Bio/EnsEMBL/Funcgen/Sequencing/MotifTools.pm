=pod

=head1 NAME

Bio::EnsEMBL::Funcgen::Sequencing::MotifTools

=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Funcgen::Sequencing::MotifTools;

use warnings;
use strict;

#use feature qw(say);
#use DBI     qw(:sql_types);

#use Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor;
use Bio::EnsEMBL::Funcgen::InputSubset;
use Bio::EnsEMBL::Funcgen::Experiment;
use Bio::EnsEMBL::Utils::SqlHelper;

#use File::Basename                        qw( fileparse );
use Bio::EnsEMBL::Utils::Scalar            qw( assert_ref );
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( dump_data
                                               run_system_cmd
                                               run_backtick_cmd
                                               write_checksum
                                               validate_path
                                               check_file
                                               validate_package_path
                                               open_file
                                               which_path );

use base qw( Exporter );
use vars qw( @EXPORT );

@EXPORT = qw(
  rev_comp_matrix 
);






# assumes a jaspar matrix pfm file with 4 rows A C G T
sub rev_comp_matrix{
  my $file    = shift;
  my $out_dir = shift; 

  my (@mat, $cols, $line, $rc_file, $suffix);
   
  #Assume whatever suffix is correct
  if($file =~ /\.([a-zA-Z]+)$/){
    $suffix = $1;
    ($rc_file = $file) =~ s/${suffix}$/rc.${suffix}/;       
  }
  else{
    throw("Could not identify matrix file suffix:\t$file");      
  }
  
  if(defined $out_dir){
    #change this file name format 
    #$rc_file = $matrix_id.'rc.'.$version.'.pfm'; 
    $rc_file =~ s/.*\///o;
    $rc_file = $out_dir.'/'.$rc_file;
  }
  
  #print $rc_file."\n";
  
  # jaspar IDs can be MA for core, CN for CNE, PB, PH and PF for PHYLOFACTS
  # our own PWM IDs are FG
  #my($matrix_id) = $file =~ /.*([MPCF][AFNG][0-9]+).pfm/;
  #my($matrix_id, $version) = $file =~ /.*([MPCF][AFNGBHL][0-9]+).([0-9]*).pfm/;
  #print $file."\n".$matrix_id." $version\n";
  #throw("Problem parsing matrix id/version from file name:\t$file") if ! defined $matrix_id;
  
  #Make this a global ReadOnly my %swap
  #so we don't have to instantiate it for every call
  #ReadOnly is not a core module though
  my %swap = (0 => 3,
              1 => 2,
              2 => 1,
              3 => 0);
  my $rows = 0;
  my $in_file = open_file($file, '<');

  #If we pass an invalid file, this will parse the whoel thing before dying
  #Probably better than slowing down the whole run by inloop or preloop validation

  while( ($line = $in_file->getline) && defined $line){
    #Need to tighter validation than this?
    next if $line !~ /[0-9]/;
    chomp $line;
        
    $line =~ s/^\s+(.*)/$1/; # remove leading whitespace
    $line =~ s/[\[\]]//g; #remove brackets if any
    #my @field = split(/\s+/,$line); # split on white space
    #print join("\t",@field)."\n";
    #print join("~",@field)."\n";
    #my @rev= reverse(@field);
    my @rev = reverse(split(/\s+/,$line));
    $mat[$swap{$rows}] = \@rev;
    #print join("\t",reverse(@field))."\n";
    $rows++;
  }
  
  close($in_file);  
  throw("Found incorrect number of rows($rows) in matrix file:\t$file") if $rows != 4;
 
  my $out_file = open_file($rc_file, '>');

  foreach my $idx(0..3){
    print $out_file join(" ",@{$mat[$idx]})."\n";
  }

  close($out_file);

    #NJ don't need this now, we just need to add support fo DB access where ever this is used
    #was this ever being used again

    # get the relevant line from the matrix_list.txt file
    # by grepping for the id at the start of the line
    # add rc to the ID and append the line to matrix_list.txt
    #my $res = &backtick("grep '^$matrix_id.$version' $work_dir".
 #                       "matrix_list.txt");
 #   print $res;
    #chop($res);
    #my @field = split("\t",$res);
    ##$field[0] .= 'rc';
    #$field[0] = $matrix_id.'rc.'.$version;
    #$res = join("\t",@field);
    ##print $res."\n";
    #open(OUT, ">> $work_dir"."matrix_list.txt") or 
    #    die "failed to open $work_dir"."matrix_list.txt for appending";
    #print OUT $res."\n" or die "failed to write to  $work_dir".
    #                           "matrix_list.txt";
    #close(OUT);

  return;

}


1;