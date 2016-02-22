#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $int_file  = "/path/to/withdrawn_interpretations.txt";
my $ccds_file = "/path/to/ccds_to_enst.txt";
my $species   = undef;                                      # 'Homo_sapiens' # 'Mus_musculus'

&GetOptions( 'int_file:s' => \$int_file, 'ccds_file:s' => \$ccds_file, 'species:s' => \$species, );

if ( !$int_file || !$ccds_file || !$species ) {
  die("Need a -ccds_file, an -int_file and a -species");
}

if ( $species ne 'Homo_sapiens' && $species ne 'Mus_musculus' ) {
  die("Species should be 'Homo_sapiens' or 'Mus_musculus'");
}

# read in the ccds file (It's a list of ENST to CCDS ids)
open( IN, $ccds_file ) or die "Can't open $ccds_file $!\n";
my %c2e;
while (<IN>) {
  chomp;
  my ( $ccds, $enst ) = split /\s+/, $_;
  push @{ $c2e{$ccds} }, $enst;
  #print "Got ccds *$ccds* with enst *$enst*\n";
}
close IN;

# read in int_file (List of ccds_id, date_time, comment, val_description, char_val, interpretation_type, interpretation_subtype, name)
my %seen;
open( IN, $int_file ) or die "Can't open $int_file $!\n";
while (<IN>) {
  chomp;
  my @cols = split /\t/, $_;
  if ( exists $c2e{ $cols[0] } ) {
    if ( !exists $seen{ $cols[0] } ) {
      foreach my $c ( @{ $c2e{ $cols[0] } } ) {
        print "http://www.ensembl.org/" . $species . "/geneview?db=core;gene=" . $c . "\n";
        $seen{ $cols[0] } = 1;
        print "$_\n";
      }
    }
    else {
      # print STDERR "Already printed url\n";
    }
  }
  else {
    print STDERR "Not found in -ccds_file : '" . $cols[0] . "'\n";
  }
}
