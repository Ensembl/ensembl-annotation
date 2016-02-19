
use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

# ensembl genes and dna
my $ensemblhost   = 'ens-staging1';
my $ensembluser   = 'ensro';
my $ensembldbname = 'homo_sapiens_core_75_37';
my $ensemblport   = 3306;
my $ensemblpass   = '';

# ccds genes
my $ccdshost   = 'ens-livemirror';
my $ccdsuser   = 'ensro';
my $ccdsdbname = 'ccds_human_75';
my $ccdsport   = 3306;
my $ccdspass   = '';

# dna genes
my $dnahost   = 'ens-staging1';
my $dnauser   = 'ensro';
my $dnadbname = 'homo_sapiens_core_75_37';
my $dnaport   = 3306;

# output
my $outfile = 'stdout';
my $verbose;

# assembly
my $coord_system_version = 'GRCh37';
my $coord_system_name    = 'chromosome';
my $ccds_logicname       = 'ccds_gene';

# region
my $chr_num;
#my $chr_num = '22';

#  ~~~~~~
#  No changes required below this point
#  ~~~~~~

# where do we get the ens structures from:
$| = 1;

&GetOptions( 'ensemblhost:s'          => \$ensemblhost,
             'ensembluser:s'          => \$ensembluser,
             'ensembldbname:s'        => \$ensembldbname,
             'ensemblport:n'          => \$ensemblport,
             'ccdshost:s'             => \$ccdshost,
             'ccdsuser:s'             => \$ccdsuser,
             'ccdsport:s'             => \$ccdsport,
             'ccdsdbname:s'           => \$ccdsdbname,
             'dnahost:s'              => \$dnahost,
             'dnauser:s'              => \$dnauser,
             'dnadbname:s'            => \$dnadbname,
             'dnaport:s'              => \$dnaport,
             'coord_system_version:s' => \$coord_system_version,
             'ccds_logicname:s'       => \$ccds_logicname,
             'outfile:s'              => \$outfile,
             'chr_num:s'              => \$chr_num, );

# print out what we are working with
print STDERR "ENSEMBL database: name $ensembldbname host $ensemblhost port $ensemblport\n";
print STDERR "CCDS database: name $ccdsdbname host $ccdshost port $ccdsport\n";
print STDERR " CCDS logic: $ccds_logicname\n";
print STDERR "\n\n\n";

# connect to dbs and get adaptors
my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $dnahost, -user => $dnauser, -port => $dnaport, -dbname => $dnadbname, );
my $dnasa = $dnadb->get_SliceAdaptor();

my $ensembldb = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $ensemblhost,
                                                    -user   => $ensembluser,
                                                    -port   => $ensemblport,
                                                    -pass   => $ensemblpass,
                                                    -dbname => $ensembldbname );
$ensembldb->dnadb($dnadb);
my $ensemblsa = $ensembldb->get_SliceAdaptor();
my $ensemblga = $ensembldb->get_GeneAdaptor();

my $ccdsdb = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $ccdshost,
                                                 -user   => $ccdsuser,
                                                 -port   => $ccdsport,
                                                 -pass   => $ccdspass,
                                                 -dbname => $ccdsdbname, );
$ccdsdb->dnadb($dnadb);
my $ccdssa = $ccdsdb->get_SliceAdaptor();
my $ccdsga = $ccdsdb->get_GeneAdaptor();

# open outfile
my $fh;
if ( $outfile && $outfile ne "stdout" ) {
  open FH, ">$outfile" or die "couldn't open file " . $outfile . " $!";
  $fh = \*FH;
}
else {
  $fh = \*STDOUT;
}

# # #
# OK, now we begin to do stuff
# # #

# fetch slices
my $ccds_slice = $ccdssa->fetch_by_region( $coord_system_name, $chr_num, undef, undef, undef, $coord_system_version );
my $ensembl_slice = $ensemblsa->fetch_by_name( $ccds_slice->name );

# # #
# now we only fetch the ccds gene biotypes that we are interested in
# # #
print STDERR "Fetching CCDS genes...\n";
my @ccds_genes = @{ do_gene_fetch( $ccds_slice, $ccds_logicname ) };
print STDERR "Got " . scalar(@ccds_genes) . " genes on " . $ccds_slice->name . "\n";

my $num_ccds_found = 0;
foreach my $ccds_gene (@ccds_genes) {
  # get ensembl genes
  my $ensembl_slice = $ensemblsa->fetch_by_name( get_genomic_location($ccds_gene) );
  my $ensembl_genes = $ensemblga->fetch_all_by_Slice($ensembl_slice);
  my $found         = confirm_ccds_in_ensembl( $ccds_gene, $ensembl_genes );
  if ( $found != 1 ) {
    print "ERROR cannot find ccds dbID " . $ccds_gene->dbID . " stable_id " . $ccds_gene->stable_id . "\n";
  }
  $num_ccds_found += $found;
}
print "Found $num_ccds_found of " . ( scalar(@ccds_genes) ) . " ccds\nDONE\n";

sub confirm_ccds_in_ensembl {
  my ( $ccds_gene, $ensembl_genes ) = @_;
  my $matched = 0;

  my $ccds_exons;
  foreach my $ccds_transc ( @{ $ccds_gene->get_all_Transcripts } ) {
    $ccds_exons = $ccds_transc->get_all_translateable_Exons;
  }
  foreach my $ensembl_gene ( @{$ensembl_genes} ) {
    next if ( $ensembl_gene->strand != $ccds_gene->strand );
    foreach my $ensembl_transc ( @{ $ensembl_gene->get_all_Transcripts } ) {
      my $ensembl_exons = $ensembl_transc->get_all_translateable_Exons;
      next if ( scalar(@$ccds_exons) != scalar(@$ensembl_exons) );
      # now compare
      my $num_exons_matched = 0;
      for ( my $i = 0; $i < scalar(@$ccds_exons); $i++ ) {
        if ( $ccds_exons->[$i]->seq_region_start == $ensembl_exons->[$i]->seq_region_start &&
             $ccds_exons->[$i]->seq_region_end == $ensembl_exons->[$i]->seq_region_end )
        {
          $num_exons_matched++;
          #print "Matched exon start ".$ccds_exons->[$i]->seq_region_start." end ".$ccds_exons->[$i]->seq_region_end." \n";
        }
      }
      if ( $num_exons_matched == scalar(@$ccds_exons) ) {
        $matched = 1;
      }
    }
  }
  return $matched;
} ## end sub confirm_ccds_in_ensembl

=head2 do_gene_fetch

  Example    : my $genes = do_gene_fetch($slice, undef, \@biotypes, undef);
  Description: Fetches genes on a slice having a specific biotype
  Returntype : Arrayref of Gene objects
  Exceptions : none

=cut

sub do_gene_fetch {
  my ( $slice, $logic ) = @_;
  my $genes = $slice->get_all_Genes($logic);
  return $genes;
}

=head2 get_genomic_location

  Example    : my $genomic_location = get_genomic_location($slice->name, $cluster);
  Description: Gets the genomic location of a cluster, as a slice name
  Returntype : String
  Exceptions : none

=cut

sub get_genomic_location {
  my ($gene) = @_;

  my @genomic_fields = split( ":", $gene->slice->name );
  my $genomic_location = $genomic_fields[0] . ":" . $genomic_fields[1] .
    ":" . $genomic_fields[2] . ":" . $gene->seq_region_start . ":" . $gene->seq_region_end . ":" . $genomic_fields[5];

  return $genomic_location;
}

