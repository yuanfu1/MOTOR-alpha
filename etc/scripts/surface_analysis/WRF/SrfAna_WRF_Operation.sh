#!/bin/bash
#
#--------------------------------------------------------------
# Script to run MOTOR-DA as an operation system.
#
# Author: Jiongming Pang (pang.j.m@hotmail.com)
#         @GBA-MWF/SIMI, 2021-11-11
#
#---------------------------------------------------------------

#SBATCH -J Opr-IDA
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --ntasks-per-node 4
#SBATCH -p YJZX
#SBATCH -o log.%j
#SBATCH -e log.%j

  export NNum=1
  export node=4
  export nNum=`expr ${node} \* ${NNum}`

  #source /public/home/qzl/spack/share/spack/setup-env.sh
  #spack load gcc@9.3.0
  #spack load openmpi@4.1.1
  #spack load cmake@3.20.5
  #spack load netcdf-fortran@4.5.3


# 1. TIME AND DIRCTORY SETTING (ONLY THIS PART CAN BE MODIFIED BY USERS):
# ==================================================================================
#
  # intitial time setting
  export INTERVAL_SRFANA=6    # in minutes
  export INTERVAL_OBSFILE=5   # in minutes
  export INTERVAL_BAKFILE=60  # in minutes
  export TIME_WINDOW=90       # in minutes
  export TIME_ZONE=8          # in hours
  export TIME_DELAY=0         # in minutes, operation time delay

  # intitial directory setting, can be modified by users
  export TOPDIR=/public/home/srf/data/surface_ana/motor/MOTOR/
  export OBSDIRsource=/public/publicdata/rawData/AWS
  export BAKDIRsource=/public/home/haps/HAPS_ENS/recentWRF3km
  cd ${TOPDIR}
  source ./pathEnv.sh
  cd ${TOPDIR}/run


# 2. GET THE TIME WINDOWS OF ANALYSIS, OBSFILES, and BAKFILES
# ==================================================================================
#
  # intitial time setting
  export BAKDIR=${TOPDIR}/grid/input/
  export OBSDIR=${TOPDIR}/data/
  export OUTDIR=${TOPDIR}/grid/output/
  
  # get the time window for analysis
  export tmp=`date +%s`
  export tmp=`expr ${tmp} - ${TIME_ZONE} \* 60 \* 60 - ${TIME_DELAY} \* 60`
  export FINL_UNIXTIME_SRFANA=$(printf "%.0f" `expr $tmp / $INTERVAL_SRFANA / 60`)
  export FINL_UNIXTIME_SRFANA=`expr $FINL_UNIXTIME_SRFANA \* $INTERVAL_SRFANA \* 60`
  export INIT_UNIXTIME_SRFANA=`expr $FINL_UNIXTIME_SRFANA - $TIME_WINDOW \* 60`

  # get the time window for observation files
  export tmp=$(printf "%.0f" `expr $FINL_UNIXTIME_SRFANA / $INTERVAL_OBSFILE / 60`)
  export FINL_UNIXTIME_OBSFILE=`expr \( $tmp + 1 \) \* $INTERVAL_OBSFILE \* 60`
  export tmp=$(printf "%.0f" `expr $INIT_UNIXTIME_SRFANA / $INTERVAL_OBSFILE / 60`)
  export INIT_UNIXTIME_OBSFILE=`expr \( $tmp \) \* $INTERVAL_OBSFILE \* 60`

  # get the time window for background files from model output
  export tmp=$(printf "%.0f" `expr $FINL_UNIXTIME_SRFANA / $INTERVAL_BAKFILE / 60`)
  export FINL_UNIXTIME_BAKFILE=`expr \( $tmp + 1 \) \* $INTERVAL_BAKFILE \* 60`
  export tmp=$(printf "%.0f" `expr $INIT_UNIXTIME_SRFANA / $INTERVAL_BAKFILE / 60`)
  export INIT_UNIXTIME_BAKFILE=`expr \( $tmp \) \* $INTERVAL_BAKFILE \* 60`


# 3. GET THE NAMES OF ALL OBS FILES AND BAK FILES, AND PREPARE INPUT FILES:
# ==================================================================================
#
  export OBSDIR=${TOPDIR}/data/
  export BAKDIR=${TOPDIR}/grid/input/
  export NMLDIR=${TOPDIR}/static/Application

  # get the names of all obs files
  export NUM_OBSFILE=`expr \( $FINL_UNIXTIME_OBSFILE - $INIT_UNIXTIME_OBSFILE \) / $INTERVAL_OBSFILE / 60 + 1 `
  export timetmp=${INIT_UNIXTIME_OBSFILE}
  for ((i=1;i<=${NUM_OBSFILE};i++)) ; do
    NAME_OBSFILE[$i]=`date "+%Y%m%d_%H%M" -d @${timetmp}`
    ln -sf ${OBSDIRsource}/${NAME_OBSFILE[$i]} ${OBSDIR}/
    if [ $i -lt ${NUM_OBSFILE} ] ; then
      NAME_OBSFILE[$i]=${NAME_OBSFILE[$i]}","
    fi
    timetmp=`expr ${timetmp} + ${INTERVAL_OBSFILE} \* 60 `
  done
  #echo ${NAME_OBSFILE[*]}

  # get the names of all bak files
  export NUM_BAKFILE=`expr \( $FINL_UNIXTIME_BAKFILE - $INIT_UNIXTIME_BAKFILE \) / $INTERVAL_BAKFILE / 60 + 1 `
  export timetmp=${INIT_UNIXTIME_BAKFILE}
  for ((i=1;i<=$NUM_BAKFILE;i++)) ; do
    NAME_BAKFILE[$i]=wrfout_d01_`date "+%Y-%m-%d_%H:%M:%S" -d @${timetmp}`
    ln -sf ${BAKDIRsource}/${NAME_BAKFILE[$i]} ${BAKDIR}/
    if [ $i -lt ${NUM_BAKFILE} ] ; then
      NAME_BAKFILE[$i]=${NAME_BAKFILE[$i]}","
    fi
    timetmp=`expr ${timetmp} + ${INTERVAL_BAKFILE} \* 60 `
  done
  #echo ${NAME_BAKFILE[*]}

  # creat namelist for surface analysis
  cp ${NMLDIR}/App_3DVAR_MultiGrid.nl.tmp ${NMLDIR}/App_3DVAR_MultiGrid.nl
  sed -i.bak "s/INIT_UNIXTIME_SRFANA/${INIT_UNIXTIME_SRFANA}/g" ${NMLDIR}/App_3DVAR_MultiGrid.nl
  sed -i.bak "s/FINL_UNIXTIME_SRFANA/${FINL_UNIXTIME_SRFANA}/g" ${NMLDIR}/App_3DVAR_MultiGrid.nl
  sed -i.bak "s/NAME_BAKFILE/${NAME_BAKFILE[*]}/g" ${NMLDIR}/App_3DVAR_MultiGrid.nl
  sed -i.bak "s/NAME_OBSFILE/${NAME_OBSFILE[*]}/g" ${NMLDIR}/App_3DVAR_MultiGrid.nl


# 4. RUN EXE
# ==================================================================================
#
  export TIME_ANA=`date "+%Y%m%d_%H%M" -d @${FINL_UNIXTIME_SRFANA}`
  export OUTDIR=${TOPDIR}/grid/output/
  export SAVEDIR=${OUTDIR}/${TIME_ANA}

  # run
  rm log.*
  mpirun -n 4 ./App_3DVAR_MultiGrid.exe >& ./log_${TIME_ANA}
  #mpirun -n 4 ./App_3DVAR_MultiGrid.exe
  #rm log.* hosts
  #srun hostname | sort -u > hosts
  #mpirun -machinefile hosts -n 4 ./App_3DVAR_MultiGrid.exe >& ./log_${TIME_ANA}

  # move files
  mkdir -p ${SAVEDIR}
  mv ${OUTDIR}/bak/bak.nc ${SAVEDIR}/bak_${TIME_ANA}.nc
  mv ${OUTDIR}/ana/ana.nc ${SAVEDIR}/ana_${TIME_ANA}.nc
  mv ./log_${TIME_ANA} ${SAVEDIR}/log_${TIME_ANA}


