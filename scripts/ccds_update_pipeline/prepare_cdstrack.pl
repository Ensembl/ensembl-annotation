#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;

my $dir = undef;

&GetOptions( 'dir:s' => \$dir, );

unless ( defined $dir ) {
  die "need to define -dir \n";
}

#3 files to edit:

my $sql  = "sql/";
my $data = "data/";

my $create_tables   = "createTables.sql";
my $interpretations = "Interpretations.txt";
my $create_keys     = "createKeys.sql";

#fix create_tables

open( IN, $dir . "/" . $sql . $create_tables ) or die "can't open $dir/" . $sql . $create_tables . "\n";
open( OUT, ">" . $dir . "/" . $sql . "new_" . $create_tables ) or die "can't open $dir/" . $sql . "new_" . $create_tables . "\n";

$/ = "\n\n";
my $line = "";

while (<IN>) {

  next if ( $_ =~ /ALTER TABLE/ );

  my $entry = $_;
  $entry =~ s/\[dbo\]\.//g;
  $entry =~ s/CREATE\s+TABLE\s+dbo\./CREATE TABLE /g;
  $entry =~ s/\bGO\b/;/gi;
  $entry =~ s/\[gap_count\] \[int\] NOT NULL ,/gap_count int NOT NULL/;
  $entry =~ s/\[//g;
  $entry =~ s/\]//g;
  $entry =~ s/IDENTITY\s*\(1,\s*1\)/UNIQUE AUTO_INCREMENT/g;
  $entry =~ s/General_//g;
  $entry =~ s/NOT FOR REPLICATION //g;
  $entry =~ s/ ON PRIMARY TEXTIMAGE_ON PRIMARY/ COLLATE=latin1_swedish_ci ENGINE=MyISAM/g;
  $entry =~ s/ ON PRIMARY/ COLLATE=latin1_swedish_ci ENGINE=MyISAM/g;

  $line .= "$entry\n";

}
$line =~ s/,\s*\n*\s*\)\s*ENGINE/\n) ENGINE/g;
print OUT $line;
$line = undef;

close IN;

#Add entry for extra table detailing Withdrawals:
print OUT "\n\nCREATE TABLE EnsemblWithdrawals (\n" .
  "\tccds_uid int NOT NULL ,\n" . "\taction ENUM( 'Keep', 'Remove transcript', 'Remove gene' ), \n" .
  "\tcomment text COLLATE Latin1_BIN NULL\n" . ") COLLATE=latin1_swedish_ci ENGINE=MyISAM\n" . ";\n";

close OUT;

#fix interpretations

open( IN, $dir . "/" . $data . $interpretations ) or die "can't open $dir/" . $data . $interpretations . "\n";
open( OUT, ">" . $dir . "/" . $data . "new_" . $interpretations ) or die "can't open $dir/" . $data . "new_" . $interpretations . "\n";

$/ = "\n";

while (<IN>) {

  my $entry = $_;
  $entry =~ s/\t\t/\t\\N\t/g;
  $entry =~ s/\t\t/\t\\N\t/g;
  $entry =~ s/\t$/\t\\N/g;
  while ( $entry =~ /\r/ ) {
    $entry =~ s/\r//g;
  }
  print OUT "$entry";

}

close IN;
close OUT;

#fix create_keys

open( IN, $dir . "/" . $sql . $create_keys ) or die "can't open $dir/" . $sql . $create_keys . "\n";
open( OUT, ">" . $dir . "/" . $sql . "new_" . $create_keys ) or die "can't open $dir/" . $sql . "new_" . $create_keys . "\n";

$/ = "\n\n";

while (<IN>) {

  my $entry = $_;

  next if ( $entry =~ /---------------------- Primary keys ------------------------/ );
  next if ( $entry =~ /---------------------- Foreign keys ------------------------/ );

  $entry =~ s/\[dbo\]\.//g;
  $entry =~ s/\bGO\b/;/gi;
  $entry =~ s/\[//g;
  $entry =~ s/\]//g;
  $entry =~ s/ WITH  FILLFACTOR = 90  ON PRIMARY//g;

  #hack to add chromosome as a primary key...
  #  if ($entry =~/PK_locations_groupVersions/){
  #    my @tmp = split/\n/, $entry;
  #    $tmp[4] = join "", $tmp[4], ",\n\t\tchromosome";
  #    $entry = join "\n", @tmp, "\n";
  #  }

  #separate out the alter_table clauses

  if ( $entry =~ /FOREIGN KEY/ ) {
    my @rows = split /CONSTRAINT/, $entry;
    my $first_line = shift @rows;
    foreach my $constraint (@rows) {
      $constraint = join " CONSTRAINT ", $first_line, $constraint;
      $constraint =~ s/\),/\);\n/;
      print OUT $constraint . "\n";
    }
  }
  else {
    print OUT "$entry";
  }

} ## end while (<IN>)

close IN;

#Add entry for Steve's extra table:
print OUT "\n\nALTER TABLE EnsemblWithdrawals ADD\n" .
  "\tCONSTRAINT PK_ensemblWithdrawals PRIMARY KEY  CLUSTERED\n" . "\t(\n" . "\t\tccds_uid\n" . "\t)\n" . ";\n";

print OUT "\n\nALTER TABLE EnsemblWithdrawals ADD\n" . "\tCONSTRAINT  FK_EnsemblWithdrawals_CcdsUids FOREIGN KEY\n" .
  "\t(\n" . "\t\tccds_uid\n" . "\t) REFERENCES CcdsUids (\n" . "\t\tccds_uid\n" . ");\n";

close OUT;

