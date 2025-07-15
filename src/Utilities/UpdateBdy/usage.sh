#!/bin/bash

#SBATCH -J DA
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH -p MOTOR
#SBATCH -o log
#SBATCH -e log

###
### Prepare 3 data files + namelist.input before your run: 
### Note start time and end time parameters in the namelist file
### the easiest way to do this is by setting both $INPUT_DIR and $OUTPUT_DIR as `pwd`
### $INPUT_BDY/namelist.input (provide grid info)
### $INPUT_BDY/grapesinput  (from SI)
### $INPUT_BDY/grapesbdy (from SI)
### $INPUT_BDY/grapesinput-MOTORDA
###
### The updated bdy file will be generatead as $OUTPUT_BDY/grapesbdy_update
###

CURRDIR=`pwd`
cd /public/home/wuyl/run/motor/MOTOR02/MOTOR
source setEnv.sh
cd $CURRDIR

export INPUT_BDY=`pwd`
echo $INPUT_BDY
export OUTPUT_BDY=`pwd`
echo $OUTPUT_BDY
export NEW_INPUT="grapesinput-MOTORDA"

/public/software/openmpi-4.0.2/bin/mpirun -n 1 -mca btl_tcp_if_include ib0 ./Test_updatebdy.exe
#cp $STATIC_DIR/../src/Utilities/UpdateBdy/grads/*ctl .

