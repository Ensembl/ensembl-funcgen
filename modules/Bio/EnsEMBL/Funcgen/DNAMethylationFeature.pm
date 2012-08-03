#
# Ensembl module for Bio::EnsEMBL::Funcgen::DNAMethylationFeature
#

=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Funcgen::DNAMethylationFeature - A data class to represent DNA methylation.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::DNAMethylationFeature;

my $df = Bio::EnsEMBL::Funcgen::DNAMethylationFeature->new

  (
    -SLICE     => $slice,
    -SET       => $resultset,
    -FILE_DATA => \@bedline,
    -ADAPTOR   =>
      $dnamethfeatadaptor # Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFeatureAdaptor object
  );
  
print "Methylated reads" . $df->methylated_reads . "\n";
print "Total reads" .$df->total_reads . "\n";
print "Percent Methylation" . $df->percent_methylation . "\n";
print "Context" . $df->context . "\n";
print "Display label" . $df->display_label . "\n";
print "Cell Type" . $df->cell_type->name . "\n";
print "Feature Type" . $df->feature_type->name . "\n";
print "Analysis Method" . $df->analysis->logic_name . "\n";



=head1 DESCRIPTION


The Bio::EnsEMBL::Funcgen::DNAMethylationFeature class represents the methylation status of a single base e.g. cytosine methylation. 

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFeatureAdaptor

Bio::EnsEMBL::Funcgen::SetFeature

=cut

package Bio::EnsEMBL::Funcgen::DNAMethylationFeature;
use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );

#use parent recommended but available from perl 5.10.1
use base qw(Bio::EnsEMBL::Funcgen::SetFeature);


=head2 new

 
  Arg [-FILE_DATA]                     : Arrayref of parsed bed line. If this argument is passed, all the rest of the arguments
                                         in the alternative arguments section become unnecessary and the constructor will throw
                                         if any of these is defined along with FILE_DATA
  #alternative arguments of FILE_DATA                                       
  Arg [-START]                         : Int - start coordinate of feature on slice. Throws if FILE_DATA defined
  Arg [-END]                           : Int - end coordinate of feature on slice. Throws if FILE_DATA defined
  Arg [-METHLATED_READS]               : Int - Number of reads indicative of methylation for the cytosine. Throws if FILE_DATA defined
  Arg [-TOTAL_READS]                   : Int - Total read coverage for the cytosine. Throws if FILE_DATA defined
  Arg [-PERCENT_METHYLATION]           : Float - Percentage of methylation (methylated reads/total reads). Throws if FILE_DATA defined
  Arg [-CONTEXT]                       : String - sequence context for the cytosine (CG, CHG or CHH). Throws if FILE_DATA defined
  Arg [-STRAND]                        : Int - The orientation of this feature. Valid values are 1, -1 and 0. Throws if FILE_DATA defined
  
  #optional argument
  Arg [-DISPLAY_LABEL]                 : String (optional) - a display label for the feature.
  
  #other arguments for super class
  Arg [-SLICE]                         : Bio::EnsEMBL::Slice - The slice on which this feature is
  Arg [-SET]                           : Bio::EnsEMBL::Funcgen::ResultSet object
  Arg [-ADAPTOR]                       : Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFeatureAdaptor object
  
  Examples:

  my $dnamethylationfeature = Bio::EnsEMBL::Funcgen::DNAMethylationFeature->new

  (
    -SLICE               => $slice,
    -STRAND              => 0,
    -START               => 1,
    -END                 => 2,
    -METHLATED_READS     => 33,
    -TOTAL_READS         => 37,
    -PERCENT_METHYLATION => 89,
    -DISPLAY_LABEL       => $display_label,
    -CONTEXT             => $context,
    -ADAPTOR             => $dnamethfeatadaptor,
    -SET                 => $resultset,
  );

my $dnamethylationfeature = Bio::EnsEMBL::Funcgen::DNAMethylationFeature->new

  (
    -SLICE         => $slice,
    -SET           => $resultset,
    -FILE_DATA     => \@bedline,
    -ADAPTOR       => $dnamethfeatadaptor,
    -DISPLAY_LABEL => $display_label,
  );

  Description            : Constructor for DNAMethylationFeature objects.
  Returntype             : Bio::EnsEMBL::Funcgen::DNAMethylationFeature
  Exceptions             : throws if appropriate arguments are not passed
  Caller                 : General
  Status                 : At Risk

=cut


sub new {
    my $caller = shift;
    my $class = ref($caller) || $caller;
    my (
        $start,                           $end,
        $strand,                          $file_data,
        $methylated_reads,                $total_reads,
        $percent_methylation,             $context,
        $dna_methylation_feature_adaptor, $display_label,
        $slice,                           $resultset
      )
      = rearrange(
        [
            'START',               'END',
            'STRAND',              'FILE_DATA',
            'METHLATED_READS',     'TOTAL_READS',
            'PERCENT_METHYLATION', 'CONTEXT',
            'ADAPTOR',             'DISPLAY_LABEL',
            'SLICE',               'SET'
        ],
        @_
      );

    #if FILE_DATA is defined, alternative arguments must be undefined

    if (   defined $file_data
        && defined $slice
        && defined $resultset
        && defined $dna_methylation_feature_adaptor )
    {

        if (   defined $start
            || defined $end
            || defined $strand
            || defined $methylated_reads
            || defined $total_reads
            || defined $percent_methylation
            || defined $context )
        {
            throw(
               "One or more alternative arguments[START|END|STRAND|METHLATED_READS|TOTAL_READS|CONTEXT] are defined along with FILE_DATA"
            );
        }
        
        #set values from the FILE_DATA

        else {
            (
                $start, $end, $strand, $methylated_reads, $total_reads,
                $percent_methylation, $context
              )
              = $dna_methylation_feature_adaptor
              ->get_DNAMethylationFeature_params_from_file_data( $file_data,
                $slice );
        }
    }
    
    #if FILE_DATA is not defined, all the alternative arguments must be defined

    elsif ( !defined $file_data
        && defined $slice
        && defined $resultset
        && defined $dna_methylation_feature_adaptor )
    {
        my $check_params =
          !defined $start ? 'START'
          : (
            !defined $end ? 'END'
            : (
                !defined $strand ? 'STRAND'
                : (
                    !defined $methylated_reads ? 'METHLATED_READS'
                    : (
                        !defined $total_reads ? 'TOTAL_READS'
                        : (
                            !defined $percent_methylation
                            ? 'PERCENT_METHYLATION'
                            : ( !defined $context ? 'CONTEXT' : undef )
                        )
                    )
                )
            )
          );
        throw(
             "Either define FILE_DATA or $check_params parameter needs to be passed"
        ) if ( defined $check_params );
    }
    
    #also throw if eithe FILE_DATA or alternative arguments are defined but one or
    #more of the arguments for super class are missing

    else {
        my $check_params =
          !defined $slice ? 'SLICE'
          : (
            !defined $resultset ? 'SET'
            : ( !defined $dna_methylation_feature_adaptor ? 'ADAPTOR' : undef )
          );
        throw("$check_params is not defined") if ( defined $check_params );
    }

    my $self = $class->SUPER::new(
        -SLICE   => $slice,
        -START   => $start,
        -END     => $end,
        -STRAND  => $strand,
        -SET     => $resultset,
        -ADAPTOR => $dna_methylation_feature_adaptor,
    );

    $self->{methylated_reads}    = $methylated_reads;
    $self->{total_reads}         = $total_reads;
    $self->{percent_methylation} = $percent_methylation;
    $self->{context}             = $context;

    $self->{display_label} = $display_label if defined($display_label);

    return $self;
}



=head2 methylated_reads

  Example    : print "Methylated reads\t". $dnamethylationfeature->methylated_reads ."\n";
  Description: Gets number of methylated reads for this DNAMethylationFeature object
  Returntype : Int
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub methylated_reads {
    return $_[0]->{methylated_reads};
}


=head2 unmethylated_reads

  Example    : print "Unmethylated reads\t". $dnamethylationfeature->unmethylated_reads . "\n";
  Description: Getter/setter of unmethylated reads for this DNAMethylationFeature object
  Returntype : Int
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut


sub unmethylated_reads {
    my $self = shift;
    return $self->{'unmethylated_reads'}
      if ( ( exists $self->{'unmethylated_reads'} ) && !defined @_ );

    if (@_) {
        $self->{'unmethylated_reads'} = shift;
    }
    else {
        $self->{'unmethylated_reads'} =
          ( $self->{'total_reads'} - $self->{'methylated_reads'} );
    }
    return $self->{'unmethylated_reads'};
}


=head2 total_reads

  Example    : print "Total reads\t". $dnamethylationfeature->total_reads . "\n";
  Description: Gets total read coverage for this DNAMethylationFeature object
  Returntype : Int
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub total_reads {
    return $_[0]->{total_reads};
}


=head2 percent_methylation

  Example    : print "Percent methylation\t". $dnamethylationfeature->percent_methylation . "\n";
  Description: Gets the percent of methylated reads for this DNAMethylationFeature object
  Returntype : Float
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub percent_methylation {
     return $_[0]->{percent_methylation};
}

=head2 context

  Example    : print "Sequence context\t". $dnamethylationfeature->context . "\n"; #(context is one of CG, CHG or CHH)
  Description: returns the sequence context for this DNAMethylationFeature object
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub context {
     return $_[0]->{context};
}

=head2 display_label

  Example    : print "Display label\t". $dnamethylationfeature->display_label . "\n";
  Description: Getter for the display_label attribute for this DNAMethylationFeature.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub display_label {
    my $self = shift;

    #$self->{'display_label'} = shift if(@_);

    if ( !$self->{'display_label'} ) {

        $self->{'display_label'} = (
            $self->cell_type->display_label()
              ||

              $self->cell_type()->name()
        );

        $self->{'display_label'} .= " " . $self->feature_type->name();

    }
    return $self->{'display_label'};
}

1;
