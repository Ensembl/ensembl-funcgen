#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $regfeat_adaptor = $registry->get_adaptor('Human', 'funcgen', 'regulatoryfeature');

# Regulatory Features: associated expermental evidence
# Print out the display_label, start/end values of all the evidence features for Human regulatory feature "ENSR00000623613".
# Compare with the start/end values of the regulatory feature itself.
my $rf = $regfeat_adaptor->fetch_by_stable_id('ENSR00000623613');
print $rf->stable_id.": \n";
print "\t".$rf->seq_region_name.":".$rf->bound_start."..".$rf->start."-".$rf->end."..".$rf->bound_end."\n";
print "\tEvidence Features: \n";
map { print_feature($_) } @{$rf->regulatory_attributes()};

sub print_feature {
	my $af = shift;
	print "\t\tDisplay Label: ".$af->display_label.";";
	print " Position: ".$af->seq_region_name.":".$af->start."-".$af->end.";";
	print "\n";
}


__END__

>perl regulatory_features_2.pl

ENSR00000623613:
        7:27178438..27178438-27185061..27185061
        Evidence Features:
                Display Label: DNase1 - HUVEC Enriched Site; Position: 7:27178438-27178886;
                Display Label: DNase1 - HUVEC Enriched Site; Position: 7:27178845-27179343;
                Display Label: DNase1 - HMEC Enriched Site; Position: 7:27178845-27179108;
                Display Label: DNase1 - H1ESC Enriched Site; Position: 7:27178986-27179243;
                Display Label: DNase1 - H1ESC Enriched Site; Position: 7:27178989-27179416;
                Display Label: DNase1 - H1ESC Enriched Site; Position: 7:27179299-27179639;
                Display Label: DNase1 - HMEC Enriched Site; Position: 7:27179426-27180165;
                Display Label: DNase1 - HUVEC Enriched Site; Position: 7:27179469-27181331;
                Display Label: DNase1 - H1ESC Enriched Site; Position: 7:27179494-27179718;
                Display Label: DNase1 - HUVEC Enriched Site; Position: 7:27180115-27180563;
                Display Label: DNase1 - HUVEC Enriched Site; Position: 7:27180681-27181089;
                Display Label: DNase1 - NHEK Enriched Site; Position: 7:27180774-27181071;
                Display Label: DNase1 - HMEC Enriched Site; Position: 7:27180886-27181213;
                Display Label: DNase1 - HUVEC Enriched Site; Position: 7:27181264-27181548;
                Display Label: DNase1 - HMEC Enriched Site; Position: 7:27181433-27181703;
                Display Label: DNase1 - HUVEC Enriched Site; Position: 7:27181453-27181969;
                Display Label: DNase1 - H1ESC Enriched Site; Position: 7:27181526-27181828;
                Display Label: DNase1 - HUVEC Enriched Site; Position: 7:27181552-27182185;
                Display Label: CTCF - HUVEC Enriched Site; Position: 7:27181590-27181895;
                Display Label: CTCF - HUVEC Enriched Site; Position: 7:27181598-27181960;
                Display Label: DNase1 - H1ESC Enriched Site; Position: 7:27181605-27181940;
                Display Label: DNase1 - K562 Enriched Site; Position: 7:27181699-27181971;
                Display Label: DNase1 - HeLa-S3 Enriched Site; Position: 7:27181881-27182143;
                Display Label: DNase1 - H1ESC Enriched Site; Position: 7:27181951-27182238;
                Display Label: DNase1 - HeLa-S3 Enriched Site; Position: 7:27182161-27183752;
                Display Label: DNase1 - HMEC Enriched Site; Position: 7:27182497-27182794;
                Display Label: Ini1 - HeLa-S3 Enriched Site; Position: 7:27182514-27183224;
                Display Label: DNase1 - IMR90 Enriched Site; Position: 7:27182547-27182888;
                Display Label: DNase1 - H1ESC Enriched Site; Position: 7:27182583-27183151;
                Display Label: DNase1 - HUVEC Enriched Site; Position: 7:27182589-27183983;
                Display Label: DNase1 - HMEC Enriched Site; Position: 7:27182813-27183810;
                Display Label: DNase1 - IMR90 Enriched Site; Position: 7:27182966-27183691;
                Display Label: DNase1 - IMR90 Enriched Site; Position: 7:27183012-27183757;
                Display Label: Max - HeLa-S3 Enriched Site; Position: 7:27183014-27183632;
                Display Label: DNase1 - HeLa-S3 Enriched Site; Position: 7:27183050-27183625;
                Display Label: DNase1 - H1ESC Enriched Site; Position: 7:27183053-27183732;
                Display Label: DNase1 - NHEK Enriched Site; Position: 7:27183061-27183626;
                Display Label: DNase1 - NH-A Enriched Site; Position: 7:27183095-27183640;
                Display Label: DNase1 - NHEK Enriched Site; Position: 7:27183113-27183569;
                Display Label: DNase1 - HUVEC Enriched Site; Position: 7:27183127-27184353;
                Display Label: DNase1 - H1ESC Enriched Site; Position: 7:27183170-27183753;
                Display Label: USF1 - H1ESC Enriched Site; Position: 7:27183191-27183596;
                Display Label: CTCF - NHEK Enriched Site; Position: 7:27183192-27183781;
                Display Label: Cmyc - HeLa-S3 Enriched Site; Position: 7:27183192-27183659;
                Display Label: Srf - HepG2 Enriched Site; Position: 7:27183197-27183525;
                Display Label: USF1 - HepG2 Enriched Site; Position: 7:27183204-27183558;
                Display Label: Srf - H1ESC Enriched Site; Position: 7:27183205-27183516;
                Display Label: CTCF - NHEK Enriched Site; Position: 7:27183238-27183747;
                Display Label: CTCF - H1ESC Enriched Site; Position: 7:27183247-27183697;
                Display Label: Max - HUVEC Enriched Site; Position: 7:27183250-27183672;     
                Display Label: CTCF - HepG2 Enriched Site; Position: 7:27183260-27183684;
                Display Label: CTCF - HeLa-S3 Enriched Site; Position: 7:27183261-27183766;
                Display Label: CTCF - HepG2 Enriched Site; Position: 7:27183263-27183681;
                Display Label: DNase1 - HepG2 Enriched Site; Position: 7:27183264-27183743;
                Display Label: Cfos - HeLa-S3 Enriched Site; Position: 7:27183265-27183596;
                Display Label: CTCF - HeLa-S3 Enriched Site; Position: 7:27183266-27183704;
                Display Label: CTCF - HMEC Enriched Site; Position: 7:27183274-27183793;
                Display Label: CTCF - HUVEC Enriched Site; Position: 7:27183275-27183887;
                Display Label: CTCF - HUVEC Enriched Site; Position: 7:27183276-27183796;
                Display Label: CTCF - HUVEC Enriched Site; Position: 7:27183285-27183696;
                Display Label: CTCF - HMEC Enriched Site; Position: 7:27183289-27183763;
                Display Label: CTCF - HepG2 Enriched Site; Position: 7:27183292-27183755;
                Display Label: CTCF - HepG2 Enriched Site; Position: 7:27183315-27183683;
                Display Label: CTCF - H1ESC Enriched Site; Position: 7:27183316-27183740;
                Display Label: DNase1 - HepG2 Enriched Site; Position: 7:27183321-27183671;
                Display Label: CTCF - HeLa-S3 Enriched Site; Position: 7:27183331-27183797;
                Display Label: Rad21 - HepG2 Enriched Site; Position: 7:27183333-27183686;
                Display Label: Rad21 - H1ESC Enriched Site; Position: 7:27183335-27183714;
                Display Label: CTCF - H1ESC Enriched Site; Position: 7:27183365-27183681;
                Display Label: CTCF - NHEK Enriched Site; Position: 7:27183370-27183789;
                Display Label: CTCF - HSMM Enriched Site; Position: 7:27183401-27183715;
                Display Label: DNase1 - HMEC Enriched Site; Position: 7:27183910-27185048;
                Display Label: DNase1 - IMR90 Enriched Site; Position: 7:27183949-27184276;
                Display Label: DNase1 - HUVEC Enriched Site; Position: 7:27184294-27184981;
                Display Label: DNase1 - H1ESC Enriched Site; Position: 7:27184403-27184592;
                Display Label: DNase1 - HUVEC Enriched Site; Position: 7:27184517-27185061;
                Display Label: Srf:MA0083.1; Position: 7:27183353-27183364;
                Display Label: Srf:PB0078.1; Position: 7:27183353-27183366;
                Display Label: Srf:MA0083.1; Position: 7:27183355-27183366;
                Display Label: Max:PB0043.1; Position: 7:27183374-27183389;
                Display Label: Max:PL0007.1; Position: 7:27183375-27183390;
                Display Label: Max:PL0007.1; Position: 7:27183375-27183390;
                Display Label: mxl-1::mdl-1:PL0014.1; Position: 7:27183375-27183390;
                Display Label: mxl-1::mdl-1:PL0014.1; Position: 7:27183375-27183390;
                Display Label: Max:PB0043.1; Position: 7:27183376-27183391;
                Display Label: Max:MA0058.1; Position: 7:27183377-27183386;
                Display Label: Cmyc:MA0059.1; Position: 7:27183377-27183387;
                Display Label: Cmyc:MA0059.1; Position: 7:27183378-27183388;
                Display Label: USF1:MA0093.1; Position: 7:27183379-27183385;
                Display Label: Max:MA0058.1; Position: 7:27183379-27183388;
                Display Label: USF1:MA0093.1; Position: 7:27183380-27183386;
                Display Label: CTCF:MA0139.1; Position: 7:27183414-27183432;
                Display Label: CTCF:MA0139.1; Position: 7:27183531-27183549;
