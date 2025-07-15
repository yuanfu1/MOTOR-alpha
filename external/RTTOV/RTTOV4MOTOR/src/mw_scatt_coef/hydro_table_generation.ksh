#!/usr/bin/ksh
#
#  To be run from the hydro table source directory (src/mw_scatt_coef/) interactively:
#
#    ./hydro_table_generation.ksh <channels.dat filename>
#
#  When compiled with OpenMP and run a modern multi-core linux workstation, the full 
#  set of hydrotables (everything in channels.dat_all) can be generated in 5 - 10 minutes.
#

# Set number of threads if using RTTOV compiled with OpenMP
export OMP_NUM_THREADS=4

# Executable name
EXEC=../../bin/rttov_scatt_make_coef.exe

# Run the executable
$EXEC $1
