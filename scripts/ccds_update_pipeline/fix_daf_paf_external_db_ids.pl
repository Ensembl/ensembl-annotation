#!/usr/bin/env perl

# $Source: /cvsroot/ensembl/ensembl-personal/genebuilders/ccds/scripts/fix_daf_paf_external_db_ids.pl,v $
# $Revision: 1.2 $

# This is a patching script and should not normally be needed

use warnings;
use strict;

use Getopt::Long;
use Bio::EnsEMBL::DBSQL::BaseAlignFeatureAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

#this database will be updated - you should make a backup first:
my $ccds_host   = '';
my $ccds_port   = '3306';
my $ccds_user   = 'ensro';
my $ccds_pass   = undef;
my $ccds_dbname = '';

#this is thr database from which ccds_db was copied
my $ref_host   = '';
my $ref_port   = '3306';
my $ref_user   = 'ensro';
my $ref_dbname = '';

#This is a file with first column=hit_name
#and second column=dna_align_feature or protein_align_feature eg.
#Q9Z1N8.1 protein_align_feature
#AW555552.2     dna_align_feature
my $hit_file;

# ccds logic_name
my $logic_name = 'ccds';

&GetOptions( 'ccds_host=s'   => \$ccds_host,
             'ccds_port=s'   => \$ccds_port,
             'ccds_user=s'   => \$ccds_user,
             'ccds_pass=s'   => \$ccds_pass,
             'ccds_dbname=s' => \$ccds_dbname,
             'ref_host=s'    => \$ref_host,
             'ref_port=s'    => \$ref_port,
             'ref_user=s'    => \$ref_user,
             'ref_dbname=s'  => \$ref_dbname,
             'file=s'        => \$hit_file,
             'logic=s'       => \$logic_name, );

if ( !$hit_file ) {
  throw("Please enter a -file");
}
if ( !$logic_name ) {
  throw("Please enter a -logic");
}

# db connections
my $ref_db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $ref_host, -user => $ref_user, -port => $ref_port, -dbname => $ref_dbname );

my $ccds_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new( -host   => $ccds_host,
                                                   -user   => $ccds_user,
                                                   -pass   => $ccds_pass,
                                                   -port   => $ccds_port,
                                                   -dbname => $ccds_dbname );

# db adaptors
my $ref_dafAdaptor = $ref_db->get_DnaAlignFeatureAdaptor;
my $ref_pafAdaptor = $ref_db->get_ProteinAlignFeatureAdaptor;

my $ccds_dafAdaptor = $ccds_db->get_DnaAlignFeatureAdaptor;
my $ccds_pafAdaptor = $ccds_db->get_ProteinAlignFeatureAdaptor;

# read in the file
open( READ, $hit_file ) or die "Can't open $hit_file $!\n";
my %hits;
while (<READ>) {
  chomp;
  my ( $hit_name, $align_feature_type ) = split /\s+/, $_;
  push @{ $hits{$align_feature_type} }, $hit_name;
}
close READ;

# now do the update
foreach my $align_feature_type ( keys %hits ) {
  print STDERR "Doing $align_feature_type...\n";

  my $ref_adaptor;
  my $ccds_adaptor;
  if ( $align_feature_type eq 'dna_align_feature' ) {
    $ref_adaptor  = $ref_db->get_DnaAlignFeatureAdaptor;
    $ccds_adaptor = $ccds_db->get_DnaAlignFeatureAdaptor;
  }
  elsif ( $align_feature_type eq 'protein_align_feature' ) {
    $ref_adaptor  = $ref_db->get_ProteinAlignFeatureAdaptor;
    $ccds_adaptor = $ccds_db->get_ProteinAlignFeatureAdaptor;
  }
  else {
    throw("align_feature_type not recognised. Must be 'dna_align_feature' or 'protein_align_feature'.\n");
  }

  # get hits
  foreach my $hit_name ( @{ $hits{$align_feature_type} } ) {
    my @ref_features  = @{ $ref_adaptor->fetch_all_by_hit_name($hit_name) };
    my @ccds_features = @{ $ccds_adaptor->fetch_all_by_hit_name($hit_name) };
    #    print STDERR "Have ".scalar(@ref_features)." ref features and ".scalar(@ccds_features)." ccds features for ".$hit_name."\n";
    update_ccds_features( \@ref_features, \@ccds_features, $align_feature_type );
  }
} ## end foreach my $align_feature_type...

sub update_ccds_features {
  my ( $ref, $ccds, $align_feature_type ) = @_;
  my @ref_features  = @{$ref};
  my @ccds_features = @{$ccds};

  # loop through
CCDS: foreach my $ccds_hit (@ccds_features) {
    # only look at the ones we need to fix
    if ( $ccds_hit->external_db_id ) {
      next CCDS;
    }
    else {
      #      print STDERR "No external_db_id for ".$ccds_hit->hseqname."\n";
    }
    if ( $ccds_hit->analysis->logic_name ne $logic_name ) {
      next CCDS;
    }
    foreach my $ref_hit (@ref_features) {
      if ( $ccds_hit->start == $ref_hit->start &&
           $ccds_hit->end == $ref_hit->end       &&
           $ccds_hit->strand == $ref_hit->strand &&
           $ccds_hit->hstart == $ref_hit->hstart &&
           $ccds_hit->hend == $ref_hit->hend     &&
           $ccds_hit->analysis->logic_name eq $ref_hit->analysis->logic_name )
      {

        #        print STDERR "Updating ".$ccds_hit->hseqname.":\n<".$ccds_hit->slice->coord_system_name."> <".$ccds_hit->start.
        #              "> <".$ccds_hit->end."> <".$ccds_hit->strand."> <".$ccds_hit->hstart."> <".
        #              $ccds_hit->hend."> <".$ccds_hit->analysis->logic_name.">\n";
        print "mysql -u" . $ccds_db->dbc->username . " -h" . $ccds_db->dbc->host . " -P" . $ccds_db->dbc->port .
          " -p" . $ccds_db->dbc->password . " -D" . $ccds_db->dbc->dbname . " -e \"update $align_feature_type set external_db_id = " .
          $ref_hit->external_db_id . " where " . "hit_start = " . $ccds_hit->hstart . " and hit_end = " . $ccds_hit->hend .
          " and hit_name = '" . $ccds_hit->hseqname . "' and seq_region_strand = " . $ccds_hit->strand . " and seq_region_start = " .
          $ccds_hit->start . " and seq_region_end = " . $ccds_hit->end . " and analysis_id = " . $ccds_hit->analysis->dbID . "\"\n\n";

      }
      else {
        #        print STDERR "Not updating ".$ccds_hit->hseqname.":\nCCDS: <".$ccds_hit->slice->coord_system_name."> <".$ccds_hit->start.
        #              "> <".$ccds_hit->end."> <".$ccds_hit->strand."> <".$ccds_hit->hstart."> <".
        #              $ccds_hit->hend.">\nREF:<".$ref_hit->slice->coord_system_name."> <".$ref_hit->start.
        #              "> <".$ref_hit->end."> <".$ref_hit->strand."> <".$ref_hit->hstart."> <".
        #              $ref_hit->hend.">\n";
      }
    }    # foreach my $ref_hit (@ref_features) {
  }    # CCDS: foreach my $ccds_hit (@ccds_features) {
} ## end sub update_ccds_features
