#!/usr/bin/perl -w

# Usage:
# ./make_metopb_iasi.pl < rtcoef_metop_2_iasi.dat
#   where rtcoef_metop_2_iasi.dat is any Metop-A IASI coefficient file.
#
# The output is a Metop-B IASI coefficient file rtcoef_metop_1_iasi.dat
# in the current directory.

use strict;

open OUT, ">rtcoef_metop_1_iasi.dat";

foreach (<>) {
  if (m/metop-2/)   { s/metop-2/metop-1/; } 
  elsif (m/sat_id/) { s/2/1/; }
  print OUT $_;
}

close OUT;
