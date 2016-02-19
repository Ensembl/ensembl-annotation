#!/usr/bin/env perl

# $Source: /cvsroot/ensembl/ensembl-personal/genebuilders/ccds/scripts/make_enst_to_ccds.pl,v $
# $Revision: 1.2 $

use warnings;
use strict;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $ccds_host   = '';
my $ccds_port   = '3306';
my $ccds_user   = 'ensro';
my $ccds_pass   = undef;
my $ccds_dbname = '';

my $host   = '';
my $port   = '';
my $user   = 'ensro';
my $pass   = undef;
my $dbname = '';

my $dna_host   = '';
my $dna_port   = '';
my $dna_user   = '';
my $dna_pass   = undef;
my $dna_dbname = '';

my $path = 'NCBI36';

&GetOptions( 'ccds_host=s'   => \$ccds_host,
             'ccds_port=s'   => \$ccds_port,
             'ccds_user=s'   => \$ccds_user,
             'ccds_pass=s'   => \$ccds_pass,
             'ccds_dbname=s' => \$ccds_dbname,
             'host=s'        => \$host,
             'port=s'        => \$port,
             'user=s'        => \$user,
             'pass=s'        => \$pass,
             'dbname=s'      => \$dbname,
             'dna_host=s'    => \$dna_host,
             'dna_port=s'    => \$dna_port,
             'dna_user=s'    => \$dna_user,
             'dna_pass=s'    => \$dna_pass,
             'dna_dbname=s'  => \$dna_dbname,
             'path=s'        => \$path, );

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $host, -user => $user, -pass => $pass, -port => $port, -dbname => $dbname );

my $ccds_db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $ccds_host,
                                                  -user   => $ccds_user,
                                                  -pass   => $ccds_pass,
                                                  -port   => $ccds_port,
                                                  -dbname => $ccds_dbname );

my $dnadb;
if ($dna_dbname) {
  $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dna_host,
                                               -user   => $dna_user,
                                               -pass   => $dna_pass,
                                               -port   => $dna_port,
                                               -dbname => $dna_dbname );

  $db->dnadb($dnadb);
  $ccds_db->dnadb($dnadb);
}

my $ccds_sa = $ccds_db->get_SliceAdaptor;
my $sa      = $db->get_SliceAdaptor;

foreach my $chr ( @{ $ccds_sa->fetch_all('chromosome') } ) {

  #print "on chr ".$chr->name."\n";

  foreach my $ccds_gene ( @{ $chr->get_all_Genes() } ) {

    $ccds_gene = $ccds_gene->transform( 'chromosome', $path );

    foreach my $ccds_trans ( @{ $ccds_gene->get_all_Transcripts() } ) {

      #find the ccds_id
      my $ccds_id;
      my @db_entries = @{ $ccds_trans->get_all_DBEntries('CCDS') };
      my %xref_hash;

      foreach my $entry (@db_entries) {
        $xref_hash{ $entry->display_id() } = 1;
      }

      if ( scalar keys %xref_hash != 1 ) {
        print "something odd going on " . scalar( keys %xref_hash ) . " xrefs\n";
        foreach my $entry ( keys %xref_hash ) {
          print "xref " . $entry . "\n";
        }
      }
      else {
        foreach my $entry ( keys %xref_hash ) {
          $ccds_id = $entry;
          #print "on ccds $ccds_id\n";
        }
      }

      my $chr_name = $ccds_trans->slice->seq_region_name;
      my $start    = $ccds_trans->start();
      my $end      = $ccds_trans->end();

      my $slice = $sa->fetch_by_region( 'chromosome', $chr_name, $start, $end, '1', $path );

      #print "slice name ".$slice->name."\n";

      my @ccds_exons = @{ $ccds_trans->get_all_translateable_Exons() };

      #print "have ".@ccds_exons." ccds exons\n";

      foreach my $gene ( @{ $slice->get_all_Genes() } ) {

        next if ( $gene->biotype ne "protein_coding" );

        #print "on gene ".$gene->display_id."\n";
        $gene = $gene->transform( 'chromosome', $path );

        foreach my $trans ( @{ $gene->get_all_Transcripts } ) {

          #print "on trans ".$trans->display_id."\n";

          my @exons = @{ $trans->get_all_translateable_Exons() };

          #print "have ".@exons." exons\n";

          my $match = 0;

          if ( scalar @exons == scalar @ccds_exons ) {

            for ( my $i = 0; $i < @exons; $i++ ) {

              if ( $ccds_exons[$i]->start == $exons[$i]->start && $ccds_exons[$i]->end == $exons[$i]->end ) {

                $match++;
              }
              #else{
              #  print "no match ".$ccds_exons[$i]->start." != ".$exons[$i]->start." or ".
              #	$ccds_exons[$i]->end." != ".$exons[$i]->end."\n";
              #}
            }

            if ( $match == scalar @exons ) {
              print $trans->stable_id . "\t" . $ccds_id . "\n";
            }
          }
        } ## end foreach my $trans ( @{ $gene...})
      } ## end foreach my $gene ( @{ $slice...})
    } ## end foreach my $ccds_trans ( @{...})
  } ## end foreach my $ccds_gene ( @{ ...})
} ## end foreach my $chr ( @{ $ccds_sa...})

