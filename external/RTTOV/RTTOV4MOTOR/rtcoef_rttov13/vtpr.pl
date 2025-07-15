#!/usr/bin/perl -w
#
# This script updates a RTTOV coefficient file with a 
# new satellite number and a new satellite height.
#
# This is usefull for VTPR for which we have a single set 
# of spectral response functions without knowing the correspondig satellites.
#
#rttov_lbl_import_satellite.exe --platform noaa --altitude 1450 --serial 2
#rttov_lbl_import_satellite.exe --platform noaa --altitude 1516 --serial 3
#rttov_lbl_import_satellite.exe --platform noaa --altitude 1460 --serial 4
# ATTENTION NOAA-5 1530km et TIROS-N 870km aussi appele noaa5
#rttov_lbl_import_satellite.exe --platform noaa --altitude 1530 --serial 5
#
#
# P. Brunel April 2012
#
use strict;
use FileHandle;

( @ARGV == 4 ) or die "Usage: $0 input_file new_sat new_alt output_file \n";

my ( $I1, $newsat, $newalt, $O ) = @ARGV;

my $fh1 = 'FileHandle'->new( "<$I1" );
my $fh2 = 'FileHandle'->new( ">$O" );

# On recopie le fichier 1 du debut a l'identificateur 
while( my $line = <$fh1> ) {
  $fh2->print( $line );
  last if ( $line =~ m/^IDENTIFICATION/ );
}

$fh2->print ( << "EOF" );
 ! 
  1  $newsat  7          ! platform sat_id instrument
 noaa-$newsat   vtpr1
EOF

# On saute les lignes du fichier 1
while( my $line = <$fh1> ) {
  last if ( $line =~ m/noaa.*vtpr1/ );
}

# On recopie le fichier 1 jusqu'a la ligne precendete de l'altitude
while( my $line = <$fh1> ) {
  $fh2->print( $line );
  last if ( $line =~ m/! Planck constants/ );
}

$fh2->print ( << "EOF" );
   $newalt                    ! nominal satellite height (km)
EOF

# On saute les lignes du fichier 1
while( my $line = <$fh1> ) {
  last if ( $line =~ m/! nominal satellite height / );
}

# On recopie le fichier 1 jusque la fin
while( my $line = <$fh1> ) {
  $fh2->print( $line );
}

