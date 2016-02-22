#!/usr/bin/env perl

use warnings;
use strict;

use Getopt::Long;

my $file = "/ensembl-personal/genebuilders/ccds/ccds_withdrawal_list.txt";

&GetOptions( 'file=s' => \$file, );

$/ = "";

open( IN, $file ) or die "Can't open $file $! \n";

while (<IN>) {
  my @lines = split /\n/, $_;

  if ( $lines[0] =~ /^CCDS/ || $lines[0] =~ /^http/ ) {
    print "expecting comment as first line - seems to be confused $_\n";
    exit;
  }

  my $action  = undef;
  my $ccds    = undef;
  my $comment = "";

  if ( $lines[0] =~ /^Keep/ ) {
    $action = 'Keep';
  }
  elsif ( $lines[0] =~ /^Already/ ) {
    $action = 'Remove transcript';
  }
  elsif ( $lines[0] =~ /^Remove transcript/ ) {
    $action = 'Remove transcript';
  }
  elsif ( $lines[0] =~ /^Remove gene/ ) {
    $action = 'Remove gene';
  }
  elsif ( $lines[0] =~ /^Remove/ ) {
    $action = 'Remove transcript';
  }
  else {
    print "Line should start with 'Keep' or 'Remove transcript|gene' don't understand $_\n";
    exit;
  }

  my @tmp = split /[-.]/, $lines[0];
  if ( @tmp > 1 ) {
    shift @tmp;
    $comment = join "-", @tmp;
  }
  else {
    $comment = $lines[0];
    #print "no comment? using $comment\n";
  }

  $comment =~ s/\t/ /g;
  $comment =~ s/^\s+//;

  foreach my $line (@lines) {

    if ( $line =~ /^CCDS(\d+)\.\d+/ ) {
      my $this_ccds = $1;

      if ( defined $ccds ) {
        if ( $ccds ne $this_ccds ) {
          print "Have more than 1 CCDS listed here $ccds and $this_ccds \n";
        }
      }
      else {
        $ccds = $this_ccds;
      }
    }
  }

  #  print "CCDS: $ccds\n";
  #  print "Action: $action\n";
  #  print "comment:$comment\n";
  #  print "***\n";

  print "$ccds\t$action\t$comment\n";

} ## end while (<IN>)
close IN;
