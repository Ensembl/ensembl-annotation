#!/usr/bin/env perl

# # #
# Making a status_file:
#
#   grep "no locations" /lustre/work1/ensembl/ba1/mouse_CCDS/load_mouse_ccds_transcripts_using_api_12mar09.out | perl -e 'while (<>) {my @fields = split(/\s+/, $_); $fields[3] =~ /(\d+)\.(\d+)/; print "mysql -NB -uensro -hgenebuild7 -Dba1_cdstrack_22apr09 -e \"select CcdsUids.ccds_uid,ccds_version ,ccds_status from CcdsUids,Groups,GroupVersions,CcdsStatusVals where CcdsUids.group_uid =GroupVersions.group_uid and Groups.group_uid = CcdsUids.group_uid and GroupVersions.ccds_status_val_uid =
#   CcdsStatusVals.ccds_status_val_uid and Groups.current_version = GroupVersions.version and CcdsUids.ccds_uid = $1 and ccds_version = $2\"\n"}' > /lustre/work1/ensembl/ba1/mouse_CCDS/no_locations.sql
#     chmod +x /lustre/work1/ensembl/ba1/mouse_CCDS/no_locations.sql
#       /lustre/work1/ensembl/ba1/mouse_CCDS/no_locations.sql > /lustre/work1/ensembl/ba1/mouse_CCDS/no_locations.sql.out
#
#
# # #

use warnings;
use strict;

use Getopt::Long;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my $host   = '';
my $port   = '3306';
my $user   = 'ensro';
my $pass   = undef;
my $dbname = 'core';

my $outhost   = '';
my $outport   = '3306';
my $outuser   = '';
my $outpass   = '';
my $outdbname = 'ccds';

my $dna_host   = '';
my $dna_port   = '3306';
my $dna_user   = 'ensro';
my $dna_pass   = undef;
my $dna_dbname = 'core';

my $analtype = "ccds_gene";
my $file     = undef;
my $status_file;

&GetOptions(
  'host=s'        => \$host,
  'port=s'        => \$port,
  'user=s'        => \$user,
  'pass=s'        => \$pass,
  'dbname=s'      => \$dbname,
  'outhost=s'     => \$outhost,
  'outport=s'     => \$outport,
  'outuser=s'     => \$outuser,
  'outpass=s'     => \$outpass,
  'outdbname=s'   => \$outdbname,
  'dna_host=s'    => \$dna_host,
  'dna_port=s'    => \$dna_port,
  'dna_user=s'    => \$dna_user,
  'dna_pass=s'    => \$dna_pass,
  'dna_dbname=s'  => \$dna_dbname,
  'analtype=s'    => \$analtype,
  'file:s'        => \$file,
  'status_file=s' => \$status_file,

);
if ( !$status_file ) {
  throw("No status file has been entered. Transcripts will have no description");
}
else {
  print "Using status file $file \n";
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $host, -user => $user, -pass => $pass, -port => $port, -dbname => $dbname );

my $outdb =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $outhost, -user => $outuser, -pass => $outpass, -port => $outport, -dbname => $outdbname );

my $dnadb;
$dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dna_host,
                                             -user   => $dna_user,
                                             -pass   => $dna_pass,
                                             -port   => $dna_port,
                                             -dbname => $dna_dbname );

$db->dnadb($dnadb);
$outdb->dnadb($dnadb);

my $ta       = $db->get_TranscriptAdaptor();
my $analysis = new Bio::EnsEMBL::Analysis( -logic_name => $analtype );
my $dea      = $db->get_DBEntryAdaptor;
my $outdea   = $outdb->get_DBEntryAdaptor;

# IN file
open( IN, $file ) or die "Can't open $file\n";

TRANSCRIPT: while (<IN>) {

  chomp;

  my @fields              = split( /\s+/, $_ );
  my $ens_trans_stable_id = $fields[0];
  my $ccdsid              = $fields[1];

  print "on ENSEMBL $ens_trans_stable_id, CCDS $ccdsid\n";
  # fetch status
  my $status = fetch_status( $ccdsid, $status_file );

  # fetch transcript from current core database
  my $ccds_trans = $ta->fetch_by_stable_id($ens_trans_stable_id);
  # we only want coding exons for CCDS
  my @ccds_exons = @{ $ccds_trans->get_all_translateable_Exons() };

  print "Got " . scalar(@ccds_exons) . " exons for " . $ccds_trans->stable_id . "\n";

  # make a new Transcript object for storing in our CCDS model db
  my $trans = new Bio::EnsEMBL::Transcript;
  $trans->strand( $ccds_trans->strand );
  $trans->analysis($analysis);
  $trans->version(1);
  if ($ccdsid) {
    $trans->stable_id($ccdsid);
  }

  # set the description ie. status
  $trans->description( $status . "_NoLocation" );
  my $slice = $outdb->get_SliceAdaptor->fetch_by_region( 'chromosome', $ccds_trans->slice->seq_region_name );

  # add all exons to the transcript
  foreach my $ccds_exon (@ccds_exons) {
    my $exon = new Bio::EnsEMBL::Exon();
    $exon->start( $ccds_exon->start );
    $exon->end( $ccds_exon->end );
    $exon->strand( $ccds_exon->strand );
    $exon->slice($slice);
    $trans->add_Exon($exon);
  }
  my @exons = @{ $trans->get_all_Exons };

  # set phases of exons
  my $phase     = 0;
  my $end_phase = undef;
  foreach my $e (@exons) {
    $e->phase($phase);
    $end_phase = ( $e->length + $phase ) % 3;
    $e->end_phase($end_phase);
    $phase = $end_phase;
  }

  # and now add a translation
  my $translation = new Bio::EnsEMBL::Translation;
  $trans->translation($translation);
  $translation->start_Exon( $exons[0] );
  $translation->end_Exon( $exons[-1] );
  $translation->start(1);
  $translation->end( $exons[-1]->length );
  $translation->version(1);
  $translation->stable_id($ccdsid);

  # finally add this all to a gene
  my $gene = new Bio::EnsEMBL::Gene;
  $gene->analysis($analysis);
  $gene->biotype($analtype);
  $gene->add_Transcript($trans);

  # store gene:
  $outdb->get_GeneAdaptor->store($gene);

  #store the CCDS_id as an xref:
  #my ($ccds_id, $ccds_version) = split/\./, $gv->get_ccds_id();
  my @db_entries = @{ $dea->fetch_all_by_Transcript( $ccds_trans, 'CCDS' ) };
  if ( scalar(@db_entries) ) {
    # sometime no db entries are returned
  DBENTRY: foreach my $dbe (@db_entries) {
      if ( $ccdsid ne $dbe->display_id ) {
        warning( "Not storing CCDS xref " . $dbe->display_id . " for $ens_trans_stable_id = $ccdsid" );
        next DBENTRY;
      }
      my ( $ccds_id, $ccds_version ) = split /\./, $dbe->display_id();
      my @trans = @{ $gene->get_all_Transcripts() };
      foreach my $t (@trans) {
        my $entry = new Bio::EnsEMBL::DBEntry( -adaptor    => $outdea,
                                               -primary_id => $ccds_id,
                                               -display_id => $ccds_id . "\." . $ccds_version,
                                               -version    => $ccds_version,
                                               -dbname     => 'CCDS',
                                               -release    => '1', );
        $entry->status("XREF");
        $outdea->store( $entry, $t->dbID, 'Transcript' );
      }
    }
  }

  my @refseqdna_dbes = @{ $dea->fetch_all_by_Transcript( $ccds_trans, 'RefSeq_mRNA' ) };
  if ( scalar(@refseqdna_dbes) ) {
    foreach my $dbe (@refseqdna_dbes) {
      my @trans = @{ $gene->get_all_Transcripts() };
      foreach my $t (@trans) {
        my $entry = new Bio::EnsEMBL::DBEntry( -adaptor    => $outdea,
                                               -primary_id => $dbe->display_id(),
                                               -display_id => $dbe->display_id(),
                                               -version    => 0,
                                               -dbname     => 'RefSeq_mRNA',
                                               -release    => '1', );
        $outdea->store( $entry, $t->dbID, 'Transcript' );
      }
    }
  }
} ## end TRANSCRIPT: while (<IN>)
close STATUS;
close IN;

sub fetch_status {
  my ( $ccdsid, $status_file ) = @_;

  open( STATUS, $status_file ) or die "Can't open $status_file\n";

  my $status;
  # columns = ccds, version, status
  while (<STATUS>) {
    my @fields = split( /\t/, $_ );
    if ( "CCDS" . $fields[0] . "." . $fields[1] eq $ccdsid ) {
      $status = $fields[2];
      $status =~ s/\s+$//;
      print "Status of $ccdsid is " . $fields[2] . "\n";
      last;
    }
  }

  close(STATUS);

  if ( !defined $status ) {
    throw("Status not defined for $ccdsid in status file");
  }

  return $status;
} ## end sub fetch_status

