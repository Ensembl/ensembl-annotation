#!/usr/bin/env perl

# $Source: /cvsroot/ensembl/ensembl-personal/genebuilders/ccds/scripts/find_unique_exons.pl,v $
# $Revision: 1.4 $

# A script designed to find cases where Ensembl
# has predicted a coding exons, but Havana
# has not predicted a coding exon.
#
# Unclustered ensembl genes, and ensembl genes
# that do not cluster with havana genes, are not
# considered.

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;

# ensembl genes before merge
my $host   = 'ens-staging1';
my $user   = 'ensro';
my $dbname = 'homo_sapiens_core_db';
my $port   = 3306;

# havana genes before merge
my $havanahost   = 'genebuild';
my $havanauser   = 'ensro';
my $havanadbname = 'some_vega_db';
my $havanaport   = 3306;

# output
# note that the file name must be of the format '*.html'
# to see the links, type 'file:///nfs/acari/ba1/WWW/*.html' in the address bar
my $file = 'stdout';

# biotypes to be considered
my @def_set_genetypes       = ('protein_coding');
my @genetypes               = ();
my $ensembl_source          = undef;
my @defhavana_set_genetypes = ('protein_coding');
my @havana_genetypes        = ();
my $havana_gene_source      = undef;                # = 'vega';
my @havana_new_genetypes    = ('havana');

# genome bits
my $coordsystem          = 'chromosome';
my $coord_system_version = 'GRCh37';
my @seq_region_names;
my $all;

# should not need to change varibales below here
my $text;
my $ens_setname    = 'ensembl';
my $havana_setname = 'havana';

$| = 1;

&GetOptions( 'host:s'                 => \$host,
             'user:s'                 => \$user,
             'dbname:s'               => \$dbname,
             'port:n'                 => \$port,
             'havanahost:s'           => \$havanahost,
             'havanauser:s'           => \$havanauser,
             'havanadbname:s'         => \$havanadbname,
             'havanaport:n'           => \$havanaport,
             'havana_genetypes:s'     => \@havana_genetypes,
             'genetypes:s'            => \@genetypes,
             'file:s'                 => \$file,
             'coord_system_version:s' => \$coord_system_version,
             'text'                   => \$text,
             'source:s'               => \$ensembl_source,
             'havana_source:s'        => \$havana_gene_source,
             'chromosomes:s'          => \@seq_region_names,
             'all'                    => \$all, );

# some checks etc...
if ( scalar(@havana_genetypes) ) {
  @havana_genetypes = split( /,/, join( ',', @havana_genetypes ) );
}
else {
  @havana_genetypes = @defhavana_set_genetypes;
}

if ( scalar(@genetypes) ) {
  @genetypes = split( /,/, join( ',', @genetypes ) );
}
else {
  @genetypes = @def_set_genetypes;
}

print STDERR "Finding unique exons in database $dbname (host $host)" . "compared to database $havanadbname (host $havanahost)\n";

my $fh;
if ( $file && $file ne "stdout" ) {
  print STDERR "Printing to file $file\n";
  open FH, ">$file" or die "couldn't open file " . $file . " $!";
  $fh = \*FH;
}
else {
  print STDERR "Printing to STDOUT\n";
  $fh = \*STDOUT;
}

# connect to dbs and get adaptors
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $host, -user => $user, -port => $port, -dbname => $dbname );
$db->assembly_type($coord_system_version);
my $sa = $db->get_SliceAdaptor();
my $ga = $db->get_GeneAdaptor();

# havana db
my $havanadb =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $havanahost, -user => $havanauser, -port => $havanaport, -dbname => $havanadbname );
$havanadb->assembly_type($coord_system_version);
my $havanasa = $havanadb->get_SliceAdaptor();
my $havanaga = $havanadb->get_GeneAdaptor();

# # #
# OK, now we begin to do stuff
# # #
my %unique_ens_trans;
print $fh "Count\tGenomic_location\tEnsembl_Transcript_stable_IDs\t" .
  "Havana_Transcript_stable_IDs\tEnsembl_Exon_stable_ID\tExon_length\tExon_number (from 5')\n";

# choose which chromosomes to do
my @slices;
if ($all) {
  if ( scalar(@seq_region_names) == 0 ) {
    # we're ok, do not filter
    foreach my $sl ( @{ $sa->fetch_all( 'toplevel', undef, undef, undef ) } ) {
      push @slices, $sl;
    }
  }
  else {
    throw("You have entered -all and -chr.");
  }
}
else {
  if ( scalar(@seq_region_names) == 0 ) {
    throw("Need either -all or chr name(s).");
  }
  else {
    foreach my $sr_name (@seq_region_names) {
      push @slices, $sa->fetch_all( 'toplevel', $sr_name );
    }
  }
}

foreach my $chr (@slices) {
  my $slicename     = $chr->name;
  my $havana_slice  = $havanasa->fetch_by_name($slicename);
  my $ensembl_slice = $sa->fetch_by_name($slicename);

  my $count;

  # # #
  # now we only fetch thr havana gene biotypes that we are interested in
  # # #
  print STDERR "Fetching havana genes for slice " . $havana_slice->name . "\n";
  my @havanagenes;
  foreach my $genetype (@havana_genetypes) {
    my @tmp = @{ $havana_slice->get_all_Genes( undef, undef, 1, $havana_gene_source, $genetype ) };
    print STDERR "Got " . scalar(@tmp) . " havana $genetype genes\n";
    foreach my $g (@tmp) {
      $g->biotype($havana_setname);
      push @havanagenes, $g;
    }
  }
  my @havana_biotype_genes;
  foreach my $h (@havanagenes) {
    $h->biotype('havana');
    push @havana_biotype_genes, $h;
  }
  # # #
  # now we only fetch the ensembl gene biotypes that we are interested in
  # # #
  print STDERR "Fetching ensembl genes for slice " . $ensembl_slice->name . "\n";
  my @genes;
  foreach my $genetype (@genetypes) {
    print STDERR "fetching $ensembl_source,$genetype\n";
    my @tmp = @{ $ensembl_slice->get_all_Genes( undef, undef, 1, $ensembl_source, $genetype ) };
    print STDERR "Got " . scalar(@tmp) . " ensembl $genetype genes\n";
    push @genes, @tmp;
  }

  # # #
  # Now look for unique ensembl exons
  # This will be done by clustering all the ensembl and havana genes
  # and then looking at each cluster individually. We want to find
  # ensembl coding exons that do not overlap with any havana coding exon.
  # We might also like to find cases where we have 2 genes and havana has
  # one gene (or vice versa).
  # # #
  my %enstypes_hash;
  my %types_hash;

  # define the types to be clustered
  $enstypes_hash{$ens_setname} = \@genetypes;
  $types_hash{$havana_setname} = \@havana_new_genetypes;
  $types_hash{$ens_setname}    = \@genetypes;
  print STDERR "    Number of havana  genes before = " . scalar(@havana_biotype_genes) . "\n";
  print STDERR "    Number of ens coding genes before    = " . scalar(@genes) . "\n";

  # Now look for useful stuff
  my @allthegenes = ( @genes, @havana_biotype_genes );
  my ( $clusters, $unclustered ) = cluster_Genes( \@allthegenes, \%types_hash, 'compare_Genes_all_Exons' );
  print STDERR "  Have " . scalar( @{$clusters} ) . " clustered and " . scalar( @{$unclustered} ) . " unclustered\n";

  # groupings
  my @single_havana;
  my @single_ensembl;

  #  foreach my $uncl (@$unclustered) {
  #   my @ensembl_genes = $uncl->get_Genes_by_Set($ens_setname);
  #   my @comp_genes = $uncl->get_Genes_by_Set($havana_setname);
  #    if (scalar(@comp_genes)==0 && scalar(@ensembl_genes)>0) {
  #      # keep this one
  #      $unique_ens_trans{'chr'} = $chr;
  #      my @transc =  @{$ensembl_genes[0]->get_all_Transcripts};
  #      $unique_ens_trans{$transc[0]->dbID} = $transc[0];
  #    }
  #  }

  # # #
  # loop thru the clusters and find interesting stuff
  # # #
  foreach my $c (@$clusters) {
    my @inc_sets = @{ $c->get_sets_included };

    #print FP "Doing cluster on chr $chr start ".$c->start." end ".$c->end."\n";
    if ( ( scalar(@inc_sets) == 1 ) && ( $inc_sets[0] eq $havana_setname ) ) {
      # clusters that contain only Havana genes are kept
      print STDERR "havana only\n";
      push @single_havana, @{ $c->get_Genes };

    }
    elsif ( ( scalar(@inc_sets) == 1 ) && ( $inc_sets[0] eq $ens_setname ) ) {
      # cluster that contains only ensembl genes
      print STDERR "ensembl only\n";
      push @single_ensembl, @{ $c->get_Genes };

    }
    else {
      my @ensembl_genes = @{ $c->get_Genes_by_Set($ens_setname) };
      my @comp_genes    = @{ $c->get_Genes_by_Set($havana_setname) };
      #print STDERR "Have ".scalar(@ensembl_genes)." ensembl and ".scalar(@comp_genes)." havana genes\n";

      # # #
      # now we want to cluster the exons
      # # #
      my @exonclusters = @{ $c->get_coding_exon_clustering_from_gene_cluster() };
    ECLUSTER: foreach my $eclust (@exonclusters) {
        print "ECLUSTER\n";
        my $eclust_got_ensembl;
        my $eclust_got_havana;
        my @exons_to_use;

        # look for a cluster that contains only e! exons
        #print "Have exon cluster start ".$eclust->start." end ".$eclust->end."\n";
      EXON: foreach my $ex ( @{ $eclust->get_all_Exons_in_ExonCluster } ) {
          my $ex_bio = $eclust->get_biotype_of_Exon($ex);
          my $exon_is_ensembl;
          my $exon_is_havana;
          # print "exon ".$ex->dbID." biotype $ex_bio\n";

          # # #
          # loop thru the 2 sets' biotypes to work out with
          # set this exon belongs to
          # # #
          foreach my $set_biotype (@genetypes) {
            # ensembl
            print "  checking $ex_bio vs $set_biotype\n";
            if ( $ex_bio eq $set_biotype ) {
              $eclust_got_ensembl = 1;
              $exon_is_ensembl    = 1;
            }
          }
          foreach my $set_biotype (@havana_new_genetypes) {
            # havana
            print "  checking $ex_bio vs $set_biotype\n";
            if ( $ex_bio eq $set_biotype ) {
              $eclust_got_havana = 1;
              $exon_is_havana    = 1;
              print "  can't use this exoncluster\n";
              next ECLUSTER;
            }
          }
          if ( !$eclust_got_ensembl && !$eclust_got_havana ) {
            throw("Biotype '$ex_bio' not recognised\n");
          }
          elsif ( $eclust_got_ensembl && $eclust_got_havana ) {
            throw("Biotype cannot be ensembl and havana");
          }
          elsif ( !$exon_is_havana && $exon_is_ensembl ) {
            print "  might use dbID " . $ex->dbID . "\n";
            push @exons_to_use, $ex;
          }
          elsif ( $exon_is_havana && !$exon_is_ensembl ) {
            print "  can't use this exoncluster, exon is havana";
          }
          else {
            throw("can't use this exoncluster");
          }
        }    # exon

        # # #
        # now that we know the biotype of the exon and thus which
        # set (ensembl or havana) it belongs to, we can collect
        # info
        # # #
        if ( !$eclust_got_havana && $eclust_got_ensembl ) {
          # this is a case of interest
          # and we want to print the results
          print "USING ";
          foreach my $e (@exons_to_use) {
            print $e->dbID . " ";
          }
          print "\n";

          # # #
          # We're in an exoncluster
          # so let's only print 1 line per exoncluster
          # # #
          # print $fh "Count\tGenomic_location\tEnsembl_Transcript_stable_IDs\t".
          #          "Havana_Transcript_stable_IDs\tEnsembl_Exon_stable_ID\tExon_length\tExon_number (from 5')\n";

          $count++;
          my %included_transcripts;
          foreach my $exon_to_use (@exons_to_use) {
            # get all transcripts with this exon
            my @included_transcripts = @{ $eclust->get_transcripts_having_this_Exon_in_ExonCluster($exon_to_use) };
            print "Have " . scalar(@included_transcripts) . " transcripts for exon " . $exon_to_use->dbID . "\n";
            foreach my $t (@included_transcripts) {
              $included_transcripts{$t} = $t;
            }
          }

          # get genomic location of cluster
          my @genomic_fields = split( ":", $ensembl_slice->name );
          my $genomic_location =
            $genomic_fields[0] . ":" . $genomic_fields[1] . ":" . $genomic_fields[2] . ":" . $c->start . ":" . $c->end . ":" . $c->strand;

          # ensembl transcript stable id string
          my @transcript_sis;
          foreach my $et ( keys %included_transcripts ) {
            if ( $included_transcripts{$et}->stable_id ) {
              push @transcript_sis, $included_transcripts{$et}->stable_id;
            }
            else {
              push @transcript_sis, $included_transcripts{$et}->dbID;
            }
          }
          my $transcript_str = join( ',', @transcript_sis );

          #havana stable id string
          my @havana_sis;
          my $strand;
          foreach my $hg (@comp_genes) {
            if ( !$strand ) {
              $strand = $hg->strand;
            }
            foreach my $ht ( @{ $hg->get_all_Transcripts } ) {
              if ( $ht->stable_id ) {
                push @havana_sis, $ht->stable_id;
              }
              else {
                push @havana_sis, $ht->dbID;
              }
            }
          }
          my $havana_si_string = join( ',', @havana_sis );

          # # #
          # and now print...
          # # #
          my @exon_ids;
          my @exon_lengths;
          foreach my $ex (@exons_to_use) {
            if ( $ex->stable_id ) {
              push @exon_ids, $ex->stable_id;
            }
            else {
              push @exon_ids, "dbID_" . $ex->dbID . ":start_" . $ex->start . ":end_" . $ex->end;
            }
            push @exon_lengths, $ex->length;
          }
          my $exon_str   = join( ',', @exon_ids );
          my $length_str = join( ',', @exon_lengths );

          print $fh $count . "\t" . $genomic_location . "\t$transcript_str\t" . "\t$havana_si_string\t$exon_str\t$length_str\t.\n";

          #          my @transcript_exons = @{$included_transcripts[0]->get_all_translateable_Exons};
          #          my $printed;
          #          for (my $x = 0; $x < scalar(@transcript_exons); $x++) {
          #            #if ($transcript_exons[$x]->stable_id eq $exon_to_use->stable_id) {
          #            if ($transcript_exons[$x]->start == $exon_to_use->start &&
          #                 $transcript_exons[$x]->end == $exon_to_use->end &&
          #                 $transcript_exons[$x]->strand == $exon_to_use->strand ) {
          #              # only print the coding length of the exon
          #              $printed = 1;
          #              print $fh ($count)."\t".$genomic_location."\t$transcript_str\t".
          #                         "\t$havana_si_string\t";
          #              if ($exon_to_use->stable_id) {
          #                print $fh $exon_to_use->stable_id."\t".$transcript_exons[$x]->length."\t".
          #                         ($x+1)." of ".scalar(@transcript_exons)."";
          #              } else {
          #                print $fh "exon_dbID".$exon_to_use->dbID.":".$exon_to_use->start.":".
          #                          $exon_to_use->end."\t".$transcript_exons[$x]->length."\t".
          #                          ($x+1)." of ".scalar(@transcript_exons)."";
          #              }
          #
          #              if ($x == 0) {
          #                print $fh " (first exon)\n";
          #              } elsif ($x == (scalar(@transcript_exons) -1)) {
          #                print $fh " (last exon)\n";
          #              } else {
          #                print $fh "\n";
          #              }
          #            }
          #          }
          #
          #          if (!$printed) {
          #            print STDERR "ERROR: Cannot use exon\n";
          #          }
        } ## end if ( !$eclust_got_havana...)
        else {
          print "NOT_USING ";
          foreach my $etu (@exons_to_use) {
            print $etu->stable_id . " ";
          }
        }    # unique
      }    # exon cluster
    }    # have ensembl and havana  in gene cluster
  }    # foreach cluster
}    # foreach chr

# # #
# we have gathered all the data
# now we just need to print it all out sensibly
# # #
# first, sort:
my @unsorted;
foreach my $transcript_dbid ( keys %unique_ens_trans ) {
  push @unsorted, $unique_ens_trans{$transcript_dbid};
}
my @sorted = sort { $a->start <=> $b->start ? $a->start <=> $b->start : $b->end <=> $a->end } @unsorted;

# now, print:
#print_html($file, \@sorted, \%unique_ens_trans, $text);

sub print_html {
  my ( $file, $sorted_transcripts, $transcripthash, $text ) = @_;
  my @sorted           = @{$sorted_transcripts};
  my %unique_ens_trans = %{$transcripthash};

  my $count;

  if ( $file ne "stdout" ) {
    open FH, ">$file";
  }
  else {
    open FH, ">-";
  }

  if ( !$text ) {
    print FH "<HTML>\n<TITLE>Ensembl-Havana Comparison</TITLE>\n<BODY>\n<P>\n";
    print FH "<H3><font color=\"red\">Comparison: Ensembl coding exons vs Havana coding exons</font></H3>\n<P>\n<P>\n";
    print FH "<TABLE border=\"1\">\n";
    print FH "<CAPTION><EM>List of ENST transcripts that have a coding exon not predicted by Havana:</EM></CAPTION>\n";
    print FH
"<TR><TH>Count</TH><TH>Genomic location</TH><TH>Ensembl</TH><TH>Havana</TH><TH>Exon</TH><TH>Exon length</TH><TH>Exon number (from 5')</TH></TR>\n";
  }
  else {
    print FH "Ensembl-Havana Comparison\n";
    print FH "Comparison: Ensembl coding exons vs Havana coding exons\n";
    print FH "List of ENST transcripts that have a coding exon not predicted by Havana:\n";
    print FH
"Count\tGenomic location\tEnsembl Transcript stable ID\tHavana Transcripts\tEnesembl Exon stable ID\tExon length\tExon number (from 5')\n";
  }

  # now loop thru and print:
  foreach my $transcript_obj (@sorted) {
    $count++;
    my %unique_havana_display_id = ();
    my %unique_havana_primary_id = ();

    my $transcript_dbid = $transcript_obj->dbID;
    #my $transcript_obj = $unique_ens_trans{$transcript_dbid};
    my $gene = $db->get_GeneAdaptor->fetch_by_transcript_id($transcript_dbid);

    # print the slice
    if ($text) {
#print FH "$count\t$coordsystem:$coord_system_version:".$unique_ens_trans{$transcript_dbid}{'chr'}.":".$transcript_obj->start.":".$transcript_obj->end.":".$transcript_obj->strand."\t";
    }
    else {
# print FH "<TR><TH>$count</TH><TH>$coordsystem:$coord_system_version:".$unique_ens_trans{$transcript_dbid}{'chr'}.":".$transcript_obj->start.
      ":" . $transcript_obj->end . ":" . $transcript_obj->strand . "</TH><TH>\t";
    }

    # print clickable link for transcript
    if ( $transcript_obj->stable_id ) {
      if ($text) {
        print FH $transcript_obj->stable_id . "\t";
      }
      else {
        print FH "<A HREF=\"http://www.ensembl.org/Homo_sapiens/transview?transcript=" .
          $transcript_obj->stable_id . "\">" . $transcript_obj->stable_id . "</A></TH><TH>\n";
      }
    }
    else {
      if ($text) {
        print FH "Unknown (may have been removed)\t";
      }
      else {
        print FH "Unknown (may have been removed)</TH><TH>\n";
      }
    }

    # print out the OTT and gene name
    foreach my $ht ( @{ $unique_ens_trans{$transcript_dbid}{'havana_transcripts'} } ) {
      my @ht_dbe = @{ $ht->get_all_DBEntries };
      foreach my $hd (@ht_dbe) {
        $unique_havana_display_id{ $hd->display_id } = 1;
        $unique_havana_primary_id{ $hd->primary_id } = 1;
      }
    }
    if ($text) {
      my @tmp_dids;
      foreach my $display_id ( keys %unique_havana_display_id ) {
        if ( $display_id =~ m/^OTT/ ) {
          push @tmp_dids, $display_id;
        }
      }
      my $string = join( ",", @tmp_dids );
      foreach my $primary_id ( keys %unique_havana_primary_id ) {
        if ( !exists $unique_havana_display_id{$primary_id} && $primary_id =~ m/^OTT/ ) {
          $string .= "," . $primary_id;
        }
      }
      print FH "$string\t";
    }
    else {
      foreach my $display_id ( keys %unique_havana_display_id ) {
        print FH "$display_id\t";
      }
      foreach my $primary_id ( keys %unique_havana_primary_id ) {
        if ( !exists $unique_havana_display_id{$primary_id} ) {
          print FH "$primary_id\t";
        }
      }
    }
    # clickable link to contigview
    if ($text) {
      print FH $unique_ens_trans{$transcript_dbid}{'exon'}->stable_id . "\t";
    }
    else {
      print FH "</TH><TH><A HREF=\"http://www.ensembl.org/Homo_sapiens/contigview?l=" . $unique_ens_trans{$transcript_dbid}{'chr'} .
        ":" . $unique_ens_trans{$transcript_dbid}{'exon'}->start . "-" . $unique_ens_trans{$transcript_dbid}{'exon'}->end .
        ";context=100\">" . $unique_ens_trans{$transcript_dbid}{'exon'}->stable_id . "</A></TH>\n";
    }
    # print exon length and position
    if ($text) {
      print FH $unique_ens_trans{$transcript_dbid}{'exon'}->length . "\t";
    }
    else {
      print FH "<TH>" . $unique_ens_trans{$transcript_dbid}{'exon'}->length . "</TH>";
    }
    my @transcript_exons = @{ $transcript_obj->get_all_translateable_Exons };
    for ( my $x = 0; $x < scalar(@transcript_exons); $x++ ) {
      if ( $transcript_exons[$x]->stable_id eq $unique_ens_trans{$transcript_dbid}{'exon'}->stable_id ) {
        if ($text) {
          print FH ( $x + 1 ) . " of " . scalar(@transcript_exons) . "";
        }
        else {
          print FH"<TH>" . ( $x + 1 ) . " of " . scalar(@transcript_exons) . "";
        }
        if ( $x == 0 ) {
          print FH " (first exon)";
        }
        elsif ( $x == ( scalar(@transcript_exons) - 1 ) ) {
          print FH " (last exon)";
        }
        if ( !$text ) {
          print FH "</TH>";
        }
      }
    }
    # close off the row
    if ($text) {
      print FH "\n";
    }
    else {
      print FH "</TR>\n\n";
    }
  } ## end foreach my $transcript_obj ...

  if ( !$text ) {
    print FH "</TABLE>\n";
    print FH "<P>\n";
    print FH "<P>\n";

    print FH "<ADDRESS>For more information please contact\n";
    print FH "<A HREF=\"MAILTO:ba1\@sanger.ac.uk\">Bronwen</A> at Ensembl.<ADDRESS>\n";
    print FH "</BODY>\n";
    print FH "</HTML>\n";
  }
} ## end sub print_html
