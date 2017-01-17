#!/usr/bin/env perl

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);

use FindBin;
use Getopt::Long;
use Config::Tiny;
use DateTime;
use Data::Printer;
use JSON;
use HTTP::Request;
use LWP::UserAgent;
use Term::ReadKey;

use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::ApiVersion;

main();

sub main {

    # -----------------
    # initialize logger
    # -----------------
    my $logger = Bio::EnsEMBL::Utils::Logger->new();

    # $logger->init_log();

    # ----------------------------
    # read command line parameters
    # ----------------------------
    my ( $relco, $release, $password, $help, $tickets_tsv, $config );

    GetOptions(
        'relco=s'    => \$relco,
        'release=i'  => \$release,
        'password=s' => \$password,
        'p=s'        => \$password,
        'tickets=s'  => \$tickets_tsv,
        'config=s'   => \$config,
        'c=s'        => \$config,
        'help'       => \$help,
        'h'          => \$help,
    );

    # ------------
    # display help
    # ------------
    if ($help) {
        usage();
    }

    # ---------------------------------
    # deal with command line parameters
    # ---------------------------------
    ( $relco, $release, $password, $tickets_tsv, $config )
        = set_parameters( $relco, $release, $password, $tickets_tsv, $config,
        $logger );

    # ---------------------------
    # read config file parameters
    # ---------------------------
    my $parameters = Config::Tiny->read($config);

    # check_dates($parameters);

    # integrate command line parameters to parameters object
    $parameters->{relco}    = $relco;
    $parameters->{password} = $password;
    $parameters->{release}  = $release;

    # ------------------
    # parse tickets file
    # ------------------
    my $tickets = parse_tickets_file( $parameters, $tickets_tsv, $logger );

    # --------------------
    # validate JIRA fields
    # --------------------
    for my $ticket ( @{$tickets} ) {

        # validate_fields( $ticket, $parameters, $logger );

        # die;
    }

    # -------------------
    # create JIRA tickets
    # -------------------
    for my $ticket ( @{$tickets} ) {
        create_ticket( $ticket, $parameters, $logger );
    }

}

sub set_parameters {
    my ( $relco, $release, $password, $tickets_tsv, $config, $logger ) = @_;

    $relco = $ENV{'USER'} if !$relco;
    $relco = validate_relco_name( $relco, $logger );

    $release = Bio::EnsEMBL::ApiVersion->software_version() if !$release;

    $tickets_tsv = $FindBin::Bin . '/jira_recurrent_tickets.tsv'
        if !$tickets_tsv;

    if ( !-e $tickets_tsv ) {
        $logger->error(
            'Tickets file '
                . $tickets_tsv
                . ' not found! Please specify one using the -tickets option!',
            0, 0
        );
    }

    $config = $FindBin::Bin . '/jira.conf' if !$config;

    if ( !-e $config ) {
        $logger->error(
            'Config file '
                . $config
                . ' not found! Please specify one using the -config option!',
            0, 0
        );
    }

    printf( "\trelco: %s\n\trelease: %i\n\ttickets: %s\n\tconfig: %s\n",
        $relco, $release, $tickets_tsv, $config );
    print "Are the above parameters correct? (y,N) : ";
    my $response = readline();
    chomp $response;
    if ( $response ne 'y' ) {
        $logger->error(
            'Aborted by user. Please rerun with correct parameters.',
            0, 0 );
    }

    if ( !$password ) {
        print 'Please type your JIRA password:';

        ReadMode('noecho');    # make password invisible on terminal
        $password = ReadLine(0);
        chomp $password;
        ReadMode(0);           # restore terminal
        print "\n";
    }

    return ( $relco, $release, $password, $tickets_tsv, $config );
}

=head2 validate_relco_name

  Arg[1]      : String $relco - the coordinator of the current release
  Arg[2]      : Bio::EnsEMBL::Utils::Logger $logger - object used for logging
  Example     : my $valid_relco = validate_relco_name($relco, $logger)
  Description : Checks if the provided relco name is valid, returns valid JIRA
                username
  Return type : String
  Exceptions  : none
  Caller      : general

=cut

sub validate_relco_name {
    my ( $relco, $logger ) = @_;

    my %valid_relco_names = (
        'ilavidas' => 'ilavidas',
        'ilias'    => 'ilavidas',
        'Ilias'    => 'ilavidas',
        'juettema' => 'juettema',
        'thomas'   => 'juettema',
        'Thomas'   => 'juettema',
        'kostadim' => 'kostadim',
        'myrto'    => 'kostadim',
        'Myrto'    => 'kostadim',
        'mnuhn'    => 'mnuhn',
        'michael'  => 'mnuhn',
        'Michael'  => 'mnuhn',
    );

    if ( exists $valid_relco_names{$relco} ) {
        return $valid_relco_names{$relco};
    }
    else {
        my $valid_names = join( "\n", sort keys %valid_relco_names );
        $logger->error(
            "Relco name not valid! Here is a list of valid names:\n"
                . $valid_names,
            0, 0
        );
    }
}

# sub check_dates {
# }

sub parse_tickets_file {

    my ( $parameters, $tickets_tsv, $logger ) = @_;
    my @tickets;

    open my $tsv, '<', $tickets_tsv;

    my $header = readline $tsv;
    chomp $header;

    #    my @jira_fields = split /\t/, $header;

    while ( readline $tsv ) {
        my $line = $_;
        chomp $line;

        $line = replace_placeholders( $line, $parameters );

        my ($project,            $issue_type, $summary,
            $reporter,           $assignee,   $priority,
            $fix_version_string, $due_date,   $component_string,
            $description
        ) = split /\t/, $line;

        # if ( scalar @jira_fields != scalar @jira_values ) {
        #     $logger->error(
        #         "Mismatch between header and values. Check JIRA ticket:\n"
        #             . $line,
        #         0, 0
        #     );
        # }
        if ($assignee) {
            $assignee = validate_user_name( $assignee, $logger );
        }

        my @fix_versions;
        my @fvs = split /,/, $fix_version_string;

        for my $fv (@fvs) {
            push @fix_versions, { 'name' => $fv };
        }

        my @components;
        my @comps = split /,/, $component_string;

        for my $comp (@comps) {
            push @components, { 'name' => $comp };
        }

        my %ticket = (
            'project'     => { 'key'  => $project },
            'issuetype'   => { 'name' => $issue_type },
            'summary'     => $summary,
            'reporter'    => { 'name' => $reporter },
            'assignee'    => { 'name' => $assignee },
            'priority'    => { 'name' => $priority },
            'fixVersions' => \@fix_versions,
            'duedate'     => $due_date,
            'components'  => \@components,
            'description' => $description,
        );

        for my $key ( keys %ticket ) {

            #delete empty fields
            if ( !$ticket{$key} ) {
                delete $ticket{$key};
            }

            # # store multiple values to arrayrefs
            # if ( lc $key eq 'fixversion' || lc $key eq 'component' ) {
            #     if ( $ticket{$key} =~ /,/ ) {
            #         $ticket{$key} =~ s/\s//g;
            #         my @values = split /,/, $ticket{$key};
            #         $ticket{$key} = \@values;
            #     }
            # }
        }

        push @tickets, \%ticket;
    }

    return \@tickets;
}

sub replace_placeholders {
    my ( $line, $parameters ) = @_;

    my $today = DateTime->now();

    $line =~ s/<RelCo>/$parameters->{relco}/g;
    $line =~ s/<version>/$parameters->{release}/g;
    $line =~ s/<preHandover_date>/$parameters->{dates}->{preHandover}/g;
    $line =~ s/<handover_date>/$parameters->{dates}->{handover}/g;
    $line =~ s/<codeBranching_date>/$parameters->{dates}->{codeBranching}/g;
    $line =~ s/<release_date>/$parameters->{dates}->{release}/g;
    $line =~ s/<date_created>/$today/g;

    return $line;
}

sub validate_fields {
    my ( $ticket, $parameters, $logger ) = @_;

    my %fields_to_be_validated = (
        'project'   => 1,
        'issuetype' => 1,
        'reporter'  => 1,
        'priority'  => 1,

        # 'fixversion' => 1,
        # 'component'  => 1,
    );

    my $endpoint = 'rest/api/latest/search';

    for my $key ( keys %{$ticket} ) {
        my $value = $ticket->{$key};
        my ( $response, $content );

        if ( $fields_to_be_validated{ lc $key } ) {
            if ( ref($value) ne 'ARRAY' ) {
                $content = { "jql" => "$key=$value", "maxResults" => 1 };
                $response = post_request( $endpoint, $content, $parameters,
                    $logger );
                say $key . "\t" . $value;

            }
            else {
                for my $element ( @{$value} ) {
                    $content
                        = { "jql" => "$key=$element", "maxResults" => 1 };
                    $response
                        = post_request( $endpoint, $content, $parameters,
                        $logger );
                    say $key . "\t" . $element;
                }
            }
        }

    }

}

sub create_ticket {
    my ( $ticket, $parameters, $logger ) = @_;
    my $endpoint = 'rest/api/latest/issue';

    my $content = { 'fields' => $ticket };
    my $response = post_request( $endpoint, $content, $parameters, $logger );
}

=head2 post_request

  Arg[1]      : String $endpoint - the request's endpoint
  Arg[2]      : Hashref $content - the request's content
  Arg[3]      : Hashref $parameters - parameters used for authorization
  Arg[4]      : Bio::EnsEMBL::Utils::Logger $logger - object used for logging
  Example     : my $response = post_request( $endpoint, $content, $parameters, $logger )
  Description : Sends a POST request to the JIRA server
  Return type : HTTP::Response object
  Exceptions  : none
  Caller      : general

=cut

sub post_request {
    my ( $endpoint, $content, $parameters, $logger ) = @_;

    my $host = 'https://www.ebi.ac.uk/panda/jira/';
    my $url  = $host . $endpoint;

    my $json_content = encode_json($content);

    my $request = HTTP::Request->new( 'POST', $url );

    $request->authorization_basic( $parameters->{relco},
        $parameters->{password} );
    $request->header( 'Content-Type' => 'application/json' );
    $request->content($json_content);

    my $agent    = LWP::UserAgent->new();
    my $response = $agent->request($request);

    if ( $response->code() == 401 ) {
        $logger->error( 'Your JIRA password is not correct. Please try again',
            0, 0 );
    }

    if ( $response->code() == 403 ) {
        $logger->error(
            'Your do not have permission to submit JIRA tickets programmatically',
            0, 0
        );
    }

    if ( !$response->is_success() ) {
        # p $response;

        # my $error_message
        #     = decode_json( $response->content() )->{errorMessages}->[0];

        my $error_message = $response->code() . ' ' . $response->message();

        $logger->error( $error_message, 0, 0 );
    }

    return $response;
}

sub usage {
    exit;
}

#TODO check date timeline
