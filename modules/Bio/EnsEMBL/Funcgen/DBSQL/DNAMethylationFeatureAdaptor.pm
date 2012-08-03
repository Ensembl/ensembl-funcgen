#
# Ensembl module for Bio::EnsEMBL::Funcgen::DNAMethylationFeatureAdaptor
#

=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.


=head1 NAME

Bio::EnsEMBL::Funcgen::DNAMethylationFeatureAdaptor - An Adaptor for fetching DNAMethylationFeatures from
BigBed format files.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
my $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor -> new (-user   => 'user', -pass => 'pass', -dbname => 'dbname',-host   => 'host', -species => 'species',-driver => 'mysql', -DNADB_USER   => 'dnadb_user', -DNADB_HOST   => 'dnadb_host', -DNADB_NAME   => 'dnadb_name', -DNADB_PORT => dnadb_port );

my $rsa =$efgdba->get_adaptor("resultset");
my @a =@{$rsa->fetch_all_by_name('ES_5mC_Stadler2011_PMID22170606')};

my $dnaa =$efgdba->get_adaptor("DNAMethylationFeature"); 
$dnaa->load_resultset($a[0]);

my $slice_adaptor = $efgdba->get_adaptor("slice");
my $slice = $slice_adaptor->fetch_by_region('chromosome',1,3010493,3011550); 

my $dna_meth_features= $dnaa ->get_DNAMethylationFeatures (-SLICE => $slice);

foreach my $df (@{$dna_meth_features})
{
print $df->methylated_reads;
print $df->total_reads;
print $df->percent_methylation;
print $df->context;
print $df->display_label;
print $df->cell_type->name;
print $df->feature_type->name;
print $df->analysis->logic_name;
.
.
.
}


=head1 DESCRIPTION

The Bio::EnsEMBL::Funcgen::DNAMethylationFeatureAdaptor uses Lincoln Stein's Bio::DB::BigFile interface to
BigBed files and creates Bio::EnsEMBL::Funcgen::DNAMethylationFeature objects. For information about 
BigBed files please see http://genome.ucsc.edu/FAQ/FAQformat.html. The fourth field of the BigBed file is expected
to contain information about cytosine context and total reads. Fifth field in the file is score which is percentage methylation multiplied by 10 (ranges from 0-1000)
An example line of the bed file that represents a cytosine in CG context with 33 methylated reads out of a total of 37 reads would be as under:

chr1    3010492 3010493 CG/37   892     +


=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DNAMethylationFeature

Bio::DB::BigFile


=cut

package Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFeatureAdaptor;
use strict;
use warnings;

use Bio::DB::BigFile;
use Bio::DB::BigFile::Constants;
use Bio::EnsEMBL::Funcgen::DNAMethylationFeature;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
#use Bio::EnsEMBL::Utils::Scalar qw(assert_ref);
#use Scalar::Util qw(weaken);
use Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor);


=head2 load_resultset

 
  Arg [1]                : Ref - Bio::EnsEMBL::Funcgen::ResultSet object
  
  Example                : $dnamethfeatadaptor->load_resultset ($resultset);


  Description            : Loads resultset into the adaptor.
  Returntype             : None
  Exceptions             : None
  Caller                 : General
  Status                 : At Risk

=cut






sub load_resultset 
		{
		   my ($self, $rs) = @_;
  		   my $bb_file_path= $rs->dbfile_data_dir . "/" . $rs->name . ".bb";
  		   $bb_file_path=~s/\/\//\//;
  		   my $bigbed_adaptor = $self -> bigBed_File_Adaptor($bb_file_path);
  		   $self->{big_bed} = $bigbed_adaptor;
  		   $self->{bb_file_path}=$bb_file_path;
  		   $self->{resultset}=$rs;
  		   #return $self;
		}  
		
		
		
=head2 bigBed_File_Path

 
  Example                : $dnamethfeatadaptor->bigBed_File_Path ("path/to/bigbed")

  Description            : returns the file path of Bigbed file used by the adaptor.
  Returntype             : String
  Exceptions             : None
  Caller                 : General
  Status                 : At Risk

=cut
		

sub bigBed_File_Path

       {
       my $self = @_;
       return $self->{bb_file_path};
       }
       		



=head2 bigBed_File_Adaptor

  Args                  : (optional) string - path to a BigBed file
  Example               : $dnamethfeatadaptor->bigBed_File_Adaptor ()
  Description           : getter/setter for the BigBed file adaptor
  Returntype            : Bio::DB::bbiFile
  Exceptions            : None
  Caller                : General
  Status                : At Risk

=cut


sub bigBed_File_Adaptor

	{
		 my ($self, $path) = @_;
		 
		 return $self->{big_bed} if (!defined $path && exists  $self->{big_bed});
		 
		 unless ((-e $path) && (-r $path))
		 {
		  throw("-BIG_BED argument ' $path ' is not a valid path or the file is not readable");
		 }
		 
		 if ((-B $path) && ($path=~/\.bb$/))
		 { 
		   my $bed = Bio::DB::BigFile->bigBedFileOpen($path);
		   if($bed) {
    				 if(!ref($bed) || !$bed->isa('Bio::DB::bbiFile'))
    						{
      						 throw("-BIG_BED argument ' $path ' could not be coerced into a Bio::DB::bbiFile object");
    						}
 					}
 					
 			return $bed;
 			}
 					
 		 else
 		 {
 		  throw("-BIG_BED argument ' $path ' does not seem to be a bigBed file");
 		 }
 		 
 	}
 	
 	
=head2 get_default_DNAMethylationFeatures

  Arg [-SLICE]           : Bio::EnsEMBL::Slice - The slice on which this feature is
  Arg [-NAME]            : (optional) string - string to be prepended to display label
  Arg [-STRAND]          : int - The orientation of this feature. Valid values are 1, -1 and 0
  
  Example                : $dnamethfeatadaptor -> get_default_DNAMethylationFeatures (-SLICE => $slice , -NAME => "ES_Biseq" , -STRAND => 0);
  
  
  Description            : Returns an ARRAYREF of Bio::EnsEMBL::Funcgen::DNAMethylationFeature objects for this slice with a default coverage
                           of greater than or equal to 10 reads. Slice objects with chromosome coordinates should be only supplied 
  Returntype             : ARRAYREF of Bio::EnsEMBL::Funcgen::DNAMethylationFeature objects
  Exceptions             : None
  Caller                 : General
  Status                 : At Risk

=cut


sub get_default_DNAMethylationFeatures

		{

		my $self = shift;

		my $refdnamethylationfeatures = $self->get_DNAMethylationFeatures( -MIN_READ_DEPTH => 10,@_);

		return $refdnamethylationfeatures;

		}


 


=head2 get_DNAMethylationFeatures

  Arg [-SLICE]           : Bio::EnsEMBL::Slice - The slice on which this feature is
  Arg [-MIN_READ_DEPTH]  : (optional) int - The minimum read depth to filter the features. If not specified, all features are returned
  
  Example                : $dnamethfeatadaptor -> get_DNAMethylationFeatures (-SLICE => $slice , -MIN_READ_DEPTH => 10);
  
  
  Description            : fetches DNAMethylationFeatures for this slice that have a total coverage of greater or equal to MIN_READ_DEPTH
  Returntype             : ARRAYREF of Bio::EnsEMBL::Funcgen::DNAMethylationFeature objects. Slice objects with chromosome coordinates should be supplied only
  Exceptions             : None
  Caller                 : General
  Status                 : At Risk

=cut 
 

#slice with toplevel coordinates only

sub get_DNAMethylationFeatures

		{
		my $self = shift;
		
		if(!ref($self) ||!$self->isa('Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFeatureAdaptor'))
    						{
      						 throw("Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFeatureAdaptor object required");
    						}
		my ($slice, $min_read_depth) = rearrange(['SLICE','MIN_READ_DEPTH'],@_);
		
		
		#hard coded for now
		my $chr="chr". $slice->seq_region_name ();
		#$name.="_" if defined $name;
		my $resultset=$self->{resultset};
		my $slice_start=$slice->start();
		my $slice_end=$slice->end();
		
		#set minimum read depth to zero if not defined
		$min_read_depth=0 unless defined $min_read_depth;
		
		my $list_head = $self->bigBed_File_Adaptor-> bigBedIntervalQuery($chr ,$slice_start - 1, $slice_end + 1);
		my @dnamethylationfeatures;
		
		for (my $i=$list_head->head;$i;$i=$i->next) 
		{
		#split the white space delimited line into a temporary array
		#the first element of the array contains "methylated reds/total reads/context" and the second contains score (0-1000)
		my @temp =split /\s+/,$i->rest;
		my $total_reads=(split /\//, $temp[0])[-1];
		
			if ($total_reads >= $min_read_depth)
		
			{
			
			     #my $percent_methylation =$temp[1]/10;
		         #my $methylated_reads= sprintf "%d", (($percent_methylation/100)*$total_reads);
		         #my $score=$temp[1];
		
		
				#my $ratio = $methylated_reads/$total_reads;
		
				#if no context defined in bigBed then the default context is set as CG
				#$context="CG" unless defined $context;
		
				#my $display_label = defined $name ? "${name}_${context}_${methylated_reads}/$total_reads" : "${context}_${methylated_reads}/$total_reads";
				#set strand if defined in the bb file as +1 for '+' and -1 for '-' or 0 if not defined or the field itself if anything other than + or -
				#my $strand= $temp[2] eq "+" ? 1 : ($temp[2] eq "-" ? -1 : (! defined $temp[2] ? 0 : $temp[2]));

	
		my $dnamethylationfeature = Bio::EnsEMBL::Funcgen::DNAMethylationFeature -> new
		
		(
		-SLICE=>$slice,
		-SET=>$resultset,
		-FILE_DATA=>[$chr,$i->start,$i->end,@temp],
		-ADAPTOR=>$self,
		);
		
		
		
		
=i			
				my $dnamethylationfeature = Bio::EnsEMBL::Funcgen::DNAMethylationFeature -> new 
		
																				(
																				
																				-START => (($i->start) +1) - $slice_start +1,
																				-END =>  $i->end - $slice_start +1,
																				-METHLATED_READS => $methylated_reads,
																				-TOTAL_READS => $total_reads,
																				-PERCENT_METHYLATION => $percent_methylation,
																				#-DISPLAY_LABEL => $display_label,
																				-CONTEXT => $context,
																				-ADAPTOR => $self,
																				-SLICE => $slice,
																				-STRAND=> $strand,
																				-SET=>$resultset,
																				);
=cut																				
				push (@dnamethylationfeatures, $dnamethylationfeature);
			}
		
		}

		return \@dnamethylationfeatures;
		}
		
		
sub get_DNAMethylationFeature_params_from_file_data

       {
           my $self = shift;
       	   if(!ref($self) ||!$self->isa('Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFeatureAdaptor'))
    						{
      						 throw("Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFeatureAdaptor object required");
    						}
    	   my ($bed_data,$slice)=@_;
    	   my $slice_start=$slice->start();
		   my $slice_end=$slice->end();
		   my $start= ($bed_data->[1]+1)-($slice_start+1);
		   my $end= ($bed_data->[2])-($slice_start+1);
		   my ($context,$total_reads)=split /\//, $bed_data->[3];
		   $context="CG" unless defined $context;
		   my $percent_methylation =($bed_data->[4])/10;
		   my $methylated_reads= sprintf "%d", (($percent_methylation/100)*$total_reads);
		   #my $display_label = defined $name ? "${name}_${context}_${methylated_reads}/$total_reads" : "${context}_${methylated_reads}/$total_reads";
		   #set strand if defined in the bb file as +1 for '+' and -1 for '-' or 0 if not defined or the field itself if anything other than + or -
		   my $strand= $bed_data->[5] eq "+" ? 1 : ($bed_data->[5] eq "-" ? -1 : (! defined $bed_data->[5] ? 0 : $bed_data->[5]));
		   
		   return ($start,$end,$strand,$methylated_reads,$total_reads,$percent_methylation,$context);
       
       }

1; 		 
