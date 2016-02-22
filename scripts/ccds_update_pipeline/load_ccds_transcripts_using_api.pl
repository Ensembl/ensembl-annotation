#!/usr/bin/env perl

use warnings;
use strict;

use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::SeqEdit;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::ExternalData::CDSTrack::GroupVersion;
use Bio::EnsEMBL::ExternalData::CDSTrack::Accession;
use Bio::EnsEMBL::ExternalData::CDSTrack::CcdsStatus;
use Bio::EnsEMBL::ExternalData::CDSTrack::Ccds;
use Bio::EnsEMBL::ExternalData::CDSTrack::Location;
use Bio::EnsEMBL::ExternalData::CDSTrack::DBSQL::DBAdaptor;
use Bio::EnsEMBL::ExternalData::Mole::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SeqEdit;

my $cds_host   = 'genebuild4';
my $cds_port   = '3306';
my $cds_user   = 'ensro';
my $cds_pass   = undef;
my $cds_dbname = 'sd3_cdstrack_mar27';

my $host   = 'genebuild6';
my $port   = '3306';
my $user   = 'ensadmin';
my $pass   = 'ensembl';
my $dbname = 'sd3_human_36_ccds_apr07';

my $bad_host   = '';
my $bad_port   = '';
my $bad_user   = '';
my $bad_pass   = undef;
my $bad_dbname = '';

my $dna_host   = 'genebuild4';
my $dna_port   = '3306';
my $dna_user   = 'ensro';
my $dna_pass   = undef;
my $dna_dbname = 'sd3_homo_sapiens_core_44_36f';

my $mole_host   = 'cbi5d';
my $mole_port   = '3306';
my $mole_user   = 'genero';
my $mole_dbname = 'refseq_36';

my $release = undef;

my $tax_id            = 9606;
my $ncbi_build_number = 'GRCh37';
my $analtype          = "ccds_import";

my $path = 'GRCh37';

my %withdrawn_status =
  ( 'Candidate' => 1, 'Withdrawn' => 1, 'Reviewed, withdrawal pending' => 1, 'Withdrawn, inconsistent annotation' => 1, );

my %bad_status = ( 'Candidate' => 1, 'Withdrawn' => 1, 'Reviewed, withdrawal pending' => 1, );

&GetOptions( 'cds_host=s'          => \$cds_host,
             'cds_port=s'          => \$cds_port,
             'cds_user=s'          => \$cds_user,
             'cds_pass=s'          => \$cds_pass,
             'cds_dbname=s'        => \$cds_dbname,
             'host=s'              => \$host,
             'port=s'              => \$port,
             'user=s'              => \$user,
             'pass=s'              => \$pass,
             'dbname=s'            => \$dbname,
             'bad_host=s'          => \$bad_host,
             'bad_port=s'          => \$bad_port,
             'bad_user=s'          => \$bad_user,
             'bad_pass=s'          => \$bad_pass,
             'bad_dbname=s'        => \$bad_dbname,
             'dna_host=s'          => \$dna_host,
             'dna_port=s'          => \$dna_port,
             'dna_user=s'          => \$dna_user,
             'dna_pass=s'          => \$dna_pass,
             'dna_dbname=s'        => \$dna_dbname,
             'mole_host=s'         => \$mole_host,
             'mole_port=s'         => \$mole_port,
             'mole_user=s'         => \$mole_user,
             'mole_dbname=s'       => \$mole_dbname,
             'tax_id=s'            => \$tax_id,
             'ncbi_build_number=s' => \$ncbi_build_number,
             'analtype=s'          => \$analtype,
             'path=s'              => \$path,
             'release=s'           => \$release, );

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $host, -user => $user, -pass => $pass, -port => $port, -dbname => $dbname );

my $bad_db;
if ($bad_host) {
  $bad_db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $bad_host,
                                                -user   => $bad_user,
                                                -pass   => $bad_pass,
                                                -port   => $bad_port,
                                                -dbname => $bad_dbname );
}

my $mole_db = new Bio::EnsEMBL::ExternalData::Mole::DBSQL::DBAdaptor( -host   => $mole_host,
                                                                      -user   => $mole_user,
                                                                      -port   => $mole_port,
                                                                      -dbname => $mole_dbname );

my $cds_db =
  Bio::EnsEMBL::ExternalData::CDSTrack::DBSQL::DBAdaptor->new( -host   => $cds_host,
                                                               -user   => $cds_user,
                                                               -pass   => $cds_pass,
                                                               -port   => $cds_port,
                                                               -dbname => $cds_dbname );

my $dnadb;
if ($dna_dbname) {
  $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dna_host,
                                               -user   => $dna_user,
                                               -pass   => $dna_pass,
                                               -port   => $dna_port,
                                               -dbname => $dna_dbname );

  $db->dnadb($dnadb);
  if ($bad_host) {
    $bad_db->dnadb($dnadb);
  }
}
my $time = time();

#my $aa   = $cds_db->get_AccessionAdaptor;
#my $csva = $cds_db->get_CcdsStatusAdaptor;
#my $cua  = $cds_db->get_CcdsAdaptor;
my $gva = $cds_db->get_GroupVersionAdaptor;
#my $la   = $cds_db->get_LocationAdaptor;

my $analysis = new Bio::EnsEMBL::Analysis( -logic_name => $analtype );
my $dea = $db->get_DBEntryAdaptor;

#to get the current set
my @gv = @{ $gva->fetch_all_current( $tax_id, $ncbi_build_number ) };

#to retrieve by status - preliminary ones have version > current_version
#my @gv;
#my @status = ('Public','Reviewed, update pending','Under review, update','Under review, withdrawal');
#my @status = ('Preliminary');
#foreach my $st (@status){
#  print "on status $st\n";
#  push @gv,  @{$gva->fetch_all_by_status($st, $tax_id, $ncbi_build_number)};
#}

GROUPVERSION: foreach my $gv (@gv) {

  print "Doing group " . $gv->group_id . " version " . $gv->dbID . "\n";
  #check the structure does not have a 'withdrawn' status
  #my $status = $csva->fetch_by_status_id($gv->ccds_status_val_id)->ccds_status;
  my $status = $gv->get_status();

  if ( !exists $withdrawn_status{$status} ) {

    my $trans = new Bio::EnsEMBL::Transcript;
    $trans->strand( $gv->strand );
    $trans->analysis($analysis);
    $trans->description($status);
    $trans->version(1);
    $trans->stable_id( $gv->get_ccds_id );
    $trans->modified_date($time);
    $trans->created_date($time);
    my $slice = $db->get_SliceAdaptor->fetch_by_region( 'chromosome', $gv->chromosome, undef, undef, undef, $path );

    #my @loc = @{$la->fetch_all_by_GroupVersion($gv)};
    my @loc = @{ $gv->get_all_Locations() };

    if ( !@loc ) {
      print "no locations for group version with group_id " . $gv->group_id . " version " . $gv->dbID . "\n";
    }
    next GROUPVERSION unless (@loc);
    print "In load script, have " . scalar(@loc) . " locations\n";

    my $count;
    foreach my $loc (@loc) {
      $count++;
      my $exon = new Bio::EnsEMBL::Exon();
      $exon->start( $loc->exon_start );
      $exon->end( $loc->exon_end );
      $exon->strand( $gv->strand );
      $exon->slice($slice);
      $exon->stable_id( $gv->get_ccds_id . "" . $count );
      $trans->add_Exon($exon);
    }

    my @exons = @{ $trans->get_all_Exons };

    #set phases
    my $phase     = 0;
    my $end_phase = undef;

    foreach my $e (@exons) {
      $e->phase($phase);
      $end_phase = ( $e->length + $phase ) % 3;
      $e->end_phase($end_phase);
      $phase = $end_phase;
    }

    my $translation = new Bio::EnsEMBL::Translation;
    $trans->translation($translation);
    $translation->version(1);
    $translation->stable_id( $gv->get_ccds_id );
    $translation->modified_date($time);
    $translation->created_date($time);

    $translation->start_Exon( $exons[0] );
    $translation->end_Exon( $exons[-1] );
    $translation->start(1);
    $translation->end( $exons[-1]->length );

    # now check that the change did the right thing
    # fetch the protein sequence
    # fetch amino acid at position X
    # check it against old
    # check transcript's codon at this position

    my $gene = new Bio::EnsEMBL::Gene;
    $gene->analysis($analysis);
    $gene->biotype($analtype);
    $gene->add_Transcript($trans);
    $gene->version(1);
    $gene->stable_id( $gv->get_ccds_id );
    $gene->modified_date($time);
    $gene->created_date($time);

    #store:
    $db->get_GeneAdaptor->store($gene);

    # # #
    # now look at translation attributes
    # and public notes
    # # #
    #  fetch gene just stored out of db
    my $fetched_gene = $db->get_GeneAdaptor->fetch_by_stable_id( $gene->stable_id );
    foreach my $transcript ( @{ $fetched_gene->get_all_Transcripts } ) {
      # check for translation attribs
      print STDERR "    Checking for translation exceptions\n";
      my @translation_exce = @{ $gv->get_all_Interpretations("Translation exception") };
      print STDERR "    Checking for public notes\n";
      my $public_note = $gv->get_Public_Note();

      if (@translation_exce) {
        # print the transcript sequence
        print "TRANSCRIPT_SEQ_BEFORE: " . $transcript->seq->seq . "\n";
        print "TRANSLATION_SEQ_BEFORE:" . $transcript->translate->seq . "\n";
      }

      foreach my $tlne (@translation_exce) {
        my $comment = $tlne->comment();
        print "COMMENT: $comment\n";
        my ( $code, $name, $description, $start, $end, $old_seq, $old_codon, $alt_seq ) = parse_comment( $tlne->comment() );
        print "Changing at start $start end $end\n";
        # now we need to parse the comment
        # add them in
        my $seq_edit = Bio::EnsEMBL::SeqEdit->new( -CODE    => $code,           # '_selenocysteine',
                                                   -NAME    => $name,           # 'Selenocysteine',
                                                   -DESC    => $description,    # 'Selenocysteine',
                                                   -START   => $start,          # 10,
                                                   -END     => $end,            # 10,
                                                   -ALT_SEQ => $alt_seq         # 'U'
        );

        # check the codon being replaced
        if ( substr( $transcript->seq->seq, 3*( $start - 1 ), 3 ) ne $old_codon ) {
          throw( "Substring to be replaced (" . $old_codon . ") does not match transcript sequence (" .
                 substr( $transcript->seq->seq, 3*$start, 3 ) . ")" );
        }

        my $attribute = $seq_edit->get_Attribute();
        #$translation->add_Attributes($attribute);
        #my $attribute_adaptor = $db->get_AttributeAdaptor();
        $db->get_AttributeAdaptor()->store_on_Translation( $transcript->translation, [$attribute] );

        print "SUBSTRING: " . substr( $transcript->seq->seq, 3*( $start - 1 ), 3 ) . "\n";
        # now re-fetch
        my $refetched_transcript = $db->get_TranscriptAdaptor->fetch_by_stable_id( $transcript->stable_id );
        print "TRANSCRIPT_SEQ_AFTER: " . $refetched_transcript->seq->seq . "\n";
        print "TRANSLATION_SEQ_AFTER:" . $refetched_transcript->translate->seq . "\n";

        print STDERR "get ccds id after tln exception " . $gv->get_ccds_id . "\n";
      }    # end translation exception
           #start public notes
      if ( !$public_note ) {
        print "No public note for " . $transcript->stable_id . "\n";
      }
      elsif ( $public_note && $public_note->comment =~ /\w+/ ) {
        print "Found public note for " . $transcript->stable_id . "\n";
        my $tln_attrib = Bio::EnsEMBL::Attribute->new(
                                          -CODE        => 'CCDS_PublicNote',
                                          -NAME        => 'CCDS Public Note',
                                          -DESCRIPTION => 'Public Note for CCDS identifier, provided by http://www.ncbi.nlm.nih.gov/CCDS',
                                          -VALUE       => $public_note->comment() );
        $db->get_AttributeAdaptor()->store_on_Translation( $transcript->translation, [$tln_attrib] );
        #
      }    # close public notes
    }    # close transcript loop

    #store the CCDS_id as an xref:
    #my $ccds_id = $cua->fetch_by_GroupVersion($gv)->ccds_id;
    my ( $ccds_id, $ccds_version ) = split /\./, $gv->get_ccds_id();

    my @trans = @{ $gene->get_all_Transcripts() };
    if ( scalar(@trans) != 1 ) {
      throw("Gene has more than one transcript");
    }

    #also wanted to store all the accessions as xrefs...
    #my @acc = @{$aa->fetch_all_by_GroupVersion($gv)};
    my @acc = @{ $gv->get_all_Accessions() };

  ACCESSION: foreach my $acc (@acc) {
      # we are reading objects out of the Acessions table in cdstrack
      my $dbname = undef;
      if ( $acc->organization eq "EBI,WTSI" ) {
        if ( $tax_id == 9606 ) {
          $dbname = "Ens_Hs_transcript";
        }
        elsif ( $tax_id == 10090 ) {
          $dbname = "Ens_Mm_transcript";
        }
        else {
          throw("Tax id $tax_id not recognised");
        }
      }
      elsif ( $acc->organization eq "NCBI" ) {
        $dbname = "RefSeq_mRNA";
      }
      else {
        throw( "Don't recognize dbname " . $acc->organization . "\n" );
      }
      print "dbname $dbname\t";

      foreach my $trans (@trans) {
        print "trans_stable_id " . $trans->stable_id . " gv_ccds_id " . $gv->get_ccds_id . "\n";
        next if ( $trans->stable_id ne $gv->get_ccds_id );
        my $entry = new Bio::EnsEMBL::DBEntry(
                                   -adaptor    => $dea,
                                   -primary_id => $acc->transcript_stable_id,
                                   -display_id => $acc->transcript_version ? $acc->transcript_stable_id . "." . $acc->transcript_version :
                                     $acc->transcript_stable_id,
                                   -version => $acc->transcript_version,
                                   -dbname  => $dbname,
                                   -release => $release, );
        $entry->status("XREF");
        $dea->store( $entry, $trans->dbID, 'Transcript' );
        # and now the protein acc
        if ( $dbname eq 'RefSeq_mRNA' ) {
          $dbname = 'RefSeq_peptide';
        }
        elsif ( $dbname eq 'Ens_Hs_transcript' ) {
          $dbname = 'Ens_Hs_translation';
        }
        elsif ( $dbname eq 'Ens_Mm_transcript' ) {
          $dbname = 'Ens_Mm_translation';
        }
        my $transl_entry = new Bio::EnsEMBL::DBEntry(
                                -adaptor    => $dea,
                                -primary_id => $acc->translation_stable_id,
                                -display_id => $acc->translation_version ? $acc->translation_stable_id . "." . $acc->translation_version :
                                  $acc->translation_stable_id,
                                -version => $acc->translation_version,
                                -dbname  => $dbname,
                                -release => $release, );
        $transl_entry->status("XREF");
        $dea->store( $transl_entry, $trans->translation->dbID, 'Translation' );
      }    # foreach ccds transcript
    }    # ACCESSION

    # now do some checks
    foreach my $trans (@trans) {
      if ( !defined $trans->translation ) {
        throw( "All transcripts require a translation: " . $trans->stable_id );
      }
      elsif ( defined $trans->translation && $trans->translation->seq =~ /\*/ ) {
        my @stops;
        # transcript is defined but has stops
        warning( "Translation has a stop codon: transcript db ID " . $trans->dbID . ", see CCDS $ccds_id" );
        print STDERR $trans->translation->seq . "\n";
        if ( defined $trans->translation && $trans->translation->seq =~ /[A-Z]*\*$/ ) {
          warning( "Translation has a stop codon at the end of the translation for transcript db ID " .
                   $trans->dbID . ", see CCDS $ccds_id. OK to continue..." );
        }
        else {
          warning( "Translation has a stop codon in the middle of the translation for transcript db ID " .
                   $trans->dbID . ", see CCDS $ccds_id. NOT OK to continue..." );
        }
      }
      if ( defined $trans->translation && $trans->translation->seq !~ /^M/ ) {
        warning( "Translation starts with non-methionine amino acid (aa = " .
                 substr( $trans->translation->seq, 0, 1 ) . " and dna = " .
                 substr( $trans->seq->seq,         0, 3 ) . "): transcript db ID " . $trans->dbID . ", see CCDS $ccds_id" );
      }

      next if ( $trans->stable_id ne $gv->get_ccds_id );

      my $entry = new Bio::EnsEMBL::DBEntry( -adaptor    => $dea,
                                             -primary_id => $ccds_id,
                                             -display_id => $ccds_id . "\." . $ccds_version,
                                             -version    => $ccds_version,
                                             -dbname     => 'CCDS',
                                             -release    => undef, );

      $entry->status("XREF");
      $dea->store( $entry, $trans->dbID, 'Transcript' );
      $dea->store( $entry, $gene->dbID,  'Gene' );

      # Update gene and transcript entries so they have a display xref
      $gene->display_xref($entry);
      $db->get_GeneAdaptor->store($gene);
      $trans->display_xref($entry);
      $db->get_TranscriptAdaptor->store($trans);
    } ## end foreach my $trans (@trans)

  } ## end if ( !exists $withdrawn_status...)
} ## end GROUPVERSION: foreach my $gv (@gv)

sub get_protein_seq_from_mole {
  my ( $db, $acc_vers ) = @_;

  my $sql = "select sequence from sequence s, entry e where e.entry_id = s.entry_id and e.accession_version = '$acc_vers'";
  print "$sql\n";
  my $sth = $db->dbc->prepare($sql);
  $sth->execute();
  my $sequence = $sth->fetchrow;
  $sth->finish();

  if ( !defined $sequence ) {
    $sql = "select sequence from sequence s, entry e where e.entry_id = s.entry_id and e.name = '$acc_vers'";
    print "$sql\n";
    $sth = $db->dbc->prepare($sql);
    $sth->execute();
    $sequence = $sth->fetchrow;
    $sth->finish();
  }
  if ( !defined $sequence ) {
    warning("Cannot fetch sequence for $acc_vers");
  }
  return $sequence;
}

sub get_first_position_of_char {
  my ( $sequence, $char, $start_index ) = @_;

  my $pos = index( $sequence, $char, $start_index );

  if ( !defined $pos ) {
    warning("Char $char not found in sequence $sequence at or after position $start_index");
  }
  return $pos;
}

sub compare_sequences {
  my ( $seq1_with_stop, $seq2 ) = @_;
  print STDERR "Comparing...\n";
  my $selenos_agree;

  # seq1, get stop position
  my @stops;
  my $transl_subseq_index = 0;
  print STDERR "Finding stops\n";
  while ( substr( $seq1_with_stop, $transl_subseq_index ) =~ /\*/ ) {
    # somewhere downstream there is a stop codon
    my $tmp_index = get_first_position_of_char( $seq1_with_stop, '*', $transl_subseq_index );
    push @stops, $tmp_index;
    $transl_subseq_index = $tmp_index + 1;
  }

  # seq2, get seleno position
  my @selenos;
  $transl_subseq_index = 0;
  print STDERR "Finding selenos\n";
  while ( substr( $seq2, $transl_subseq_index ) =~ /U/ ) {
    my $tmp_index = get_first_position_of_char( $seq2, 'U', $transl_subseq_index );
    push @selenos, $tmp_index;
    print STDERR "tmp_index $tmp_index, transl_subseq_index $transl_subseq_index\n";
    $transl_subseq_index = $tmp_index + 1;
  }

  # now check that the arrays contain the same values
  for ( my $i = 0; $i < scalar(@stops); $i++ ) {
    print STDERR "    ccds seq1 found stop at $i\n";
  }
  for ( my $j = 0; $j < scalar(@selenos); $j++ ) {
    print STDERR "    refseq seq2 found stop at $j\n";
  }

  # now check that the arrays contain the same values
  if ( scalar(@stops) == scalar(@selenos) ) {
    for ( my $i = 0; $i < scalar(@stops); $i++ ) {
      if ( $stops[$i] == $selenos[$i] ) {
        $selenos_agree++;
      }
    }
  }

  if ( scalar(@stops) && $selenos_agree != scalar(@stops) ) {
    $selenos_agree = undef;
  }
  return $selenos_agree;
} ## end sub compare_sequences

sub parse_comment {
  my ($comment) = @_;
  # eg. replace the symbol 'L' (codon CTG) with 'M' at amino acid position 0

  my ( $attrib_code,   $attrib_name, $attrib_desc );
  my ( $old_seq,       $old_codon );
  my ( $seqedit_start, $seqedit_end, $new_seq );

  if ( $comment =~ /^replace the symbol '([A-Z*])' \(codon ([ATGC]{3})\) with '([A-Z])' at amino acid position (\d+)$/ ) {
    $old_seq       = $1;
    $old_codon     = $2;
    $new_seq       = $3;
    $seqedit_start = $4;
    $seqedit_end   = $4;
  }
  else {
    throw("Comment not recognised");
  }

  # the ncbi numbers start counting at 0 but we start at 1
  #$seqedit_start += 1;
  #$seqedit_end += 1;

  # now decide which type of exception it is
  if ( $old_seq eq '*' && $new_seq eq 'U' ) {
    # this is a seleno
    $attrib_code = '_selenocysteine';
    $attrib_name = 'Selenocysteine';
    $attrib_desc = 'The UGA codon at this position translates to a selenocysteine and not a stop codon';

  }
  elsif ( $old_seq eq '*' && $new_seq =~ /[A-Z]/ ) {
    # stop codon should be an amino acid
    $attrib_code = 'amino_acid_sub';            # UGA exception
    $attrib_name = 'Amino acid substitution';
    $attrib_desc =
'Some translations have been manually curated for amino acid substitiutions. For example a stop codon may be changed to an amino acid in order to prevent premature truncation, or one amino acid can be substituted for another.';

  }
  elsif ( $old_seq =~ /[A-Z]/ && $new_seq eq 'M' && $seqedit_start == 1 ) {
    # COMMENT: replace the symbol 'I' (codon ATT) with 'M' at amino acid position 1
    $attrib_code = 'amino_acid_sub';            # initial methionine
    $attrib_name = 'Amino acid substitution';
    $attrib_desc =
'Some translations have been manually curated for amino acid substitiutions. For example a stop codon may be changed to an amino acid in order to prevent premature truncation, or one amino acid can be substituted for another.';

  }
  else {
    $attrib_code = $old_seq . "_" . $old_codon . "_" . $new_seq;
    $attrib_name = $comment;
    $attrib_desc = $comment;
  }

  return ( $attrib_code, $attrib_name, $attrib_desc, $seqedit_start, $seqedit_end, $old_seq, $old_codon, $new_seq );
} ## end sub parse_comment
