#
# Ensembl module for Bio::EnsEMBL::Funcgen::BindingMatrix
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::BindingMatrix - A module to represent a BindingMatrix. 
In EFG this represents the binding affinities of a Transcription Factor to DNA.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::BindingMatrix;

my $matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new(
                                                        -name  => "MA0122.1",
                                                        -type => "Jaspar",
                                                        -description => "Nkx3-2 Jaspar Matrix",
                                                       );
$matrix->frequencies("A  [ 4  1 13 24  0  0  6  4  9 ]
C  [ 7  4  1  0  0  0  0  6  7 ]
G  [ 4  5  7  0 24  0 18 12  5 ]
T  [ 9 14  3  0  0 24  0  2  3 ]");

print $matrix->compare_to_optimal_site("TGGCCACCA")."\n";
                                                                
=head1 DESCRIPTION

This is a simple class to represent information about a BindingMatrix, containing the name 
(e.g. the Jaspar ID, or an internal name), its type ('Jaspar' or 'Inferred')

 and description. 

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor

=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.


=cut


use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::BindingMatrix;

use Bio::EnsEMBL::Utils::Argument qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::Storable;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Storable);

#SOME CONSTANTS
my %available_matrix_types = (
							  jaspar   => 'Jaspar',
							  inferred => 'Inferred',
							 );

=head2 new

  Arg [-name]: string - name of Matrix
  Arg [-type]: string - type of Matrix ('Jaspar' or 'Inferred')
  Arg [-frequencies]: (optional) string - frequencies representing the binding affinities of a Matrix
  Arg [-description]: (optional) string - descriptiom of Matrix
  Example    : my $matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new(
                                                               -name  => "MA0122.1",
                                                               -type => "Jaspar",
                                                               -description => "Jaspar Matrix",
                                                                );
  Description: Constructor method for BindingMatrix class
  Returntype : Bio::EnsEMBL::Funcgen::BindingMatrix
  Exceptions : Throws if name or/and type not defined
  Caller     : General
  Status     : Medium risk

=cut

sub new {
  my $caller = shift;

  my $obj_class = ref($caller) || $caller;
  my $self = $obj_class->SUPER::new(@_);
  
  my ( $name, $type, $freq, $desc, $ftype, $assoc_ftypes ) = rearrange
	( [
	   'NAME', 'TYPE', 'FREQUENCIES', 'DESCRIPTION', 'FEATURE_TYPE',
	  ], @_);
  
  
  if(! defined $name){
    throw("Must supply a name\n");
  }

  if(! defined $type){
    throw("Must supply a type\n");
  }
  
  if(! (ref($ftype) && $ftype->isa('Bio::EnsEMBL::Funcgen::FeatureType'))){
	throw("You must define a valid Bio::EnsEMBL::Funcgen::FeatureType");
	#leave is stored test to adaptor
  }

  $self->name($name);
  $self->type($type);
  $self->{feature_type} = $ftype;
  $self->frequencies($freq) if $freq;
  $self->description($desc) if $desc;

  return $self;
}



=head2 feature_type

  Example    : my $ft_name = $efeature->feature_type()->name();
  Description: Getter for the feature_type attribute for this feature.
  Returntype : Bio::EnsEMBL::Funcgen:FeatureType
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub feature_type{
  my $self = shift;
  
  return $self->{'feature_type'};
}


=head2 name

  Arg [1]    : (optional) string - name
  Example    : my $name = $matrix->name();
  Description: Getter and setter of name attribute
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub name {
    my $self = shift;
    $self->{'name'} = shift if @_;
    return $self->{'name'};
}

=head2 description

  Arg [1]    : (optional) string - description
  Example    : my $desc = $matrix->description();
  Description: Getter and setter of description attribute 
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Low Risk

=cut

sub description {
    my $self = shift;
    $self->{'description'} = shift if @_;
    return $self->{'description'};
}


=head2 type

  Arg [1]    : (optional) string - type
  Example    : $matrix->type('Jaspar');
  Description: Getter and setter of type attribute
  Returntype : string
  Exceptions : Throws if type not present in %available_matrix_types
  Caller     : General
  Status     : At Risk

=cut

sub type {
  my ($self, $type) = @_;
  
  if(defined $type){
	
  	if(exists $available_matrix_types{lc($type)}){ 
	  $self->{'type'} =  $available_matrix_types{lc($type)};
  	} else {
	  throw("$type is not a valid type: ".join(', ', values(%available_matrix_types)));	
  	}
  }

  return $self->{'type'};
}


=head2 frequencies

  Arg [1]    : (optional) string - frequencies
  Example    : $matrix->frequencies($frequencies_string);
  Description: Getter and setter of frequencies attribute

  The attribute is a string representing the matrix binding
  affinities in the Jaspar format. E.g. 
  	">
  	[ ]
  	"
  
  Returntype : string
  Exceptions : Throws if the string attribute is not a properly
  formed matrix in the Jaspar format
  Caller     : General
  Status     : At Risk

=cut

sub frequencies {
  my $self = shift;
 
  #TODO Make Base Backgroung more realistic than equiprobability
 
  my $frequencies = shift if @_; 
  if($frequencies){
  	$self->_weights($frequencies);
  	$self->{'frequencies'} = $frequencies;  	
  }
  return $self->{'frequencies'};
}

=head2 frequencies_revcomp

  Example    : $matrix->frequencies_revcomp();
  Description: Getter for the reverse complement frequencies attribute

  The attribute represents the reverse complement of frequencies
  
  Returntype : string
  Caller     : General
  Status     : At Risk

=cut

sub frequencies_revcomp {
  my $self = shift;
 
  return $self->{'frequencies_revcomp'};
}


=head2 relative_affinity

  Arg [1]    : string - Binding Site Sequence
  Arg [2]    : (optional) boolean - 1 if results are to be in linear scale (default is log scale)
  Example    : $matrix->relative_affinity($sequence);
  Description: Calculates the binding affinity of a given sequence
	relative to the optimal site for the matrix
	The site is taken as if it were in the proper orientation
        Considers a purely random background p(A)=p(C)=p(G)=p(T)
  Returntype : double
  Exceptions : Throws if the sequence length does not have the matrix length
  or if the sequence has unclear bases (N is not accepted)
  Caller     : General
  Status     : At Risk

=cut

sub relative_affinity {
  my ($self, $sequence, $linear) = (shift, shift, shift);
  $sequence =~ s/^\s+//;
  $sequence =~ s/\s+$//;
  
  throw "No sequence given" if !$sequence;
  $sequence = uc($sequence);
  if($sequence =~ /[^ACGT]/){
    throw "Sequence $sequence contains invalid characters: Only Aa Cc Gg Tt accepted";	
  }
  
  my $weight_matrix = $self->_weights;
  my $matrix_length = scalar(@{$weight_matrix->{'A'}});
  if(length($sequence) != $matrix_length){
    throw "Sequence $sequence does not have length $matrix_length";
  }
  
  my $log_odds = 0;
  my @bases = split(//,$sequence);
  for(my $i=0;$i<$matrix_length;$i++){
    $log_odds += $weight_matrix->{$bases[$i]}->[$i];	
  }
  
  #This log scale may be quite unrealistic... but usefull just for comparisons...
  if(!$linear){
    return ($log_odds - $self->_min_bind) / ($self->_max_bind - $self->_min_bind);
  } else {
    return (exp($log_odds) - exp($self->_min_bind)) / (exp($self->_max_bind) - exp($self->_min_bind));
  }
  #Linearizing it may be too strict, on the other hand...
  #return (exp($log_odds) - exp($self->_min_bind)) / (exp($self->_max_bind) - exp($self->_min_bind));
  #Comparing just to maximum may be useful...
  #return $log_odds / $self->_max_bind;
  #Otherwise return just the raw result
  #return $log_odds;
}

=head2 length

  Example    : $bm->length();
  Description: Returns the length of the the matrix (e.g. 19bp long)
  Returntype : int with the length of this binding matrix 
  Exceptions : none
  Caller     : General
  Status     : At Risk

=cut

sub length {
  my $self = shift;

  my $weight_matrix = $self->_weights;

  return scalar(@{$weight_matrix->{'A'}});
}

=head2 _weights

  Arg [1]    : (optional) string - frequencies
  Example    : _weights($frequencies);
  Description: Private Getter Setter for the weight matrix based on frequencies
  Returntype : HASHREF with the weights of this binding matrix 
  Exceptions : Throws if the frequencies attribute string does not correspond 
  to 4 rows of an equal number of integer numbers.
  Caller     : Self
  Status     : At Risk

=cut

sub _weights {
	my $self = shift;
	
 	#for the moment use equiprobability and constant pseudo-count
 	#TODO allow for it to be passed as parameters?
  	my $frequencies = shift if @_; 
  	if($frequencies){
  		$frequencies =~ s/^(>.*?\n)//;
		my $header = $1;

  		my ($a,$c,$g,$t) = split(/\n/,$frequencies);
  		my @As = split(/\s+/,_parse_matrix_line('[A\[\]]',$a));
  		my @Cs = split(/\s+/,_parse_matrix_line('[C\[\]]',$c));
  		my @Gs = split(/\s+/,_parse_matrix_line('[G\[\]]',$g));
  		my @Ts = split(/\s+/,_parse_matrix_line('[T\[\]]',$t));
		if((scalar(@As)!=scalar(@Cs)) || (scalar(@As)!=scalar(@Gs)) || (scalar(@As)!=scalar(@Ts)) ){
			throw "Frequencies provided are not a valid frequency matrix"		
		}
		
		#Create the reverse complement
		my @revT = reverse(@As);
		my @revA = reverse(@Ts);
		my @revC = reverse(@Gs);
		my @revG = reverse(@Cs);
		my $revcomp = $header;
		$revcomp.= "A [ ".join("\t",@revA)." ]\n";
		$revcomp.= "C [ ".join("\t",@revC)." ]\n";
		$revcomp.= "G [ ".join("\t",@revG)." ]\n";
		$revcomp.= "T [ ".join("\t",@revT)." ]\n";
		$self->{'frequencies_revcomp'} = $revcomp;

  		my @totals;
  		for(my $i=0;$i<scalar(@As);$i++){ 
  			$totals[$i]=$As[$i]+$Cs[$i]+$Gs[$i]+$Ts[$i];
  		}
  		
		my %weights;			
		#We can allow distinct background per nucleotide, instead of 0.25 for all... pass as parameter
		#But if the matrix was obtained using in-vivo data, it shouldn't matter the organism nucleotide bias..
		#We're using 0.1 as pseudo-count... the matrix cannot have very few elements... (e.g. <30 not good)
		my @was; for(my $i=0;$i<scalar(@As);$i++){ $was[$i] = log((($As[$i] + 0.1) / ($totals[$i]+0.4)) / 0.25); };
		$weights{'A'} = \@was;
		my @wcs; for(my $i=0;$i<scalar(@Cs);$i++){ $wcs[$i] = log((($Cs[$i] + 0.1) / ($totals[$i]+0.4)) / 0.25); };
		$weights{'C'} = \@wcs;
		my @wgs; for(my $i=0;$i<scalar(@Gs);$i++){ $wgs[$i] = log((($Gs[$i] + 0.1) / ($totals[$i]+0.4)) / 0.25); };
		$weights{'G'} = \@wgs;
		my @wts; for(my $i=0;$i<scalar(@Ts);$i++){ $wts[$i] = log((($Ts[$i] + 0.1) / ($totals[$i]+0.4)) / 0.25); };	     
		$weights{'T'} = \@wts;
	
		$self->{'weights'} = \%weights;

		my $max = 0; my $min = 0;
		for(my $i=0;$i<scalar(@As);$i++){
			my $col = [ $was[$i], $wcs[$i], $wgs[$i], $wts[$i] ];
			$min += _min($col);
			$max += _max($col);
		}

		#Log scale
		$self->_max_bind($max);
		$self->_min_bind($min);
	}
	
	return $self->{'weights'};	

}

sub _parse_matrix_line {
	 my ($pat,$line) = (shift,shift);
	 $line=~s/$pat//g; $line=~s/^\s+//; $line=~s/\s+$//;	
	 return $line;
}

sub _max { return _min_max(shift, 0); }

sub _min { return _min_max(shift, 1); }

sub _min_max {
	my ($list,$min) = (shift, shift);
	my $min_max = $list->[0];
	map { if($min ? $_ < $min_max : $_ > $min_max){ $min_max = $_; } } @$list;
	return $min_max;
}


=head2 _max_bind

  Arg [1]    : (optional) double - maximum binding affinity
  Example    : $matrix->_max_bind(10.2);
  Description: Private Getter and setter of max_bind attribute (not to be called directly)
  Returntype : float with the maximum binding affinity of the matrix 
  Exceptions : None
  Caller     : Self
  Status     : At Risk

=cut

sub _max_bind {
  my $self = shift;
  
  $self->{'max_bind'} = shift if @_;

  return $self->{'max_bind'};
}

=head2 _min_bind

  Arg [1]    : (optional) double - minimum binding affinity
  Example    : $matrix->_min_bind(-10.2);
  Description: Private Getter and setter of min_bind attribute (not to be called directly)
  Returntype : float with the minimum binding affinity of the matrix 
  Exceptions : None
  Caller     : Self
  Status     : At Risk

=cut

sub _min_bind {
  my $self = shift;
  
  $self->{'min_bind'} = shift if @_;

  return $self->{'min_bind'};
}

1;
