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
#SBATCH -n 16
#SBATCH --ntasks-per-node 16
#SBATCH -p YJZX
#SBATCH -o log.%j
#SBATCH -e log.%j

  export NNum=1
  export node=16
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
  export TIME_DELAY=5         # in minutes, operation time delay

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
  export INDATADIR=${TOPDIR}/input/
  export OUTDATADIR=${TOPDIR}/output/
  
  # get the time window for analysis
  export tmp=`date +%s`
  export tmp=`expr ${tmp} - ${TIME_ZONE} \* 60 \* 60 - ${TIME_DELAY} \* 60`
  export FINL_UNIXTIME_SRFANA=$(printf "%.0f" `expr $tmp / $INTERVAL_SRFANA / 60`)
  export FINL_UNIXTIME_SRFANA=`expr $FINL_UNIXTIME_SRFANA \* $INTERVAL_SRFANA \* 60`
  export INIT_UNIXTIME_SRFANA=`expr $FINL_UNIXTIME_SRFANA - $TIME_WINDOW \* 60`

  export sYear=`date "+%_Y" -d @${INIT_UNIXTIME_SRFANA}`
  export sMon=`date "+%_m" -d @${INIT_UNIXTIME_SRFANA}`
  export sDay=`date "+%_d" -d @${INIT_UNIXTIME_SRFANA}`
  export sHor=`date "+%_H" -d @${INIT_UNIXTIME_SRFANA}`
  export sMin=`date "+%_M" -d @${INIT_UNIXTIME_SRFANA}`
  export sSec=`date "+%_S" -d @${INIT_UNIXTIME_SRFANA}`

  export eYear=`date "+%_Y" -d @${FINL_UNIXTIME_SRFANA}`
  export eMon=`date "+%_m" -d @${FINL_UNIXTIME_SRFANA}`
  export eDay=`date "+%_d" -d @${FINL_UNIXTIME_SRFANA}`
  export eHor=`date "+%_H" -d @${FINL_UNIXTIME_SRFANA}`
  export eMin=`date "+%_M" -d @${FINL_UNIXTIME_SRFANA}`
  export eSec=`date "+%_S" -d @${FINL_UNIXTIME_SRFANA}`

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
  export NMLDIR=${TOPDIR}/static/Application

  # get the names of all obs files
  export NUM_OBSFILE=`expr \( $FINL_UNIXTIME_OBSFILE - $INIT_UNIXTIME_OBSFILE \) / $INTERVAL_OBSFILE / 60 + 1 `
  export timetmp=${INIT_UNIXTIME_OBSFILE}
  for ((i=1;i<=${NUM_OBSFILE};i++)) ; do
    NAME_OBSFILE[$i]=`date "+%Y%m%d_%H%M" -d @${timetmp}`
    ln -sf ${OBSDIRsource}/${NAME_OBSFILE[$i]} ${INDATADIR}/
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
    ln -sf ${BAKDIRsource}/${NAME_BAKFILE[$i]} ${INDATADIR}/
    if [ $i -lt ${NUM_BAKFILE} ] ; then
      NAME_BAKFILE[$i]=${NAME_BAKFILE[$i]}","
    fi
    timetmp=`expr ${timetmp} + ${INTERVAL_BAKFILE} \* 60 `
  done
  #echo ${NAME_BAKFILE[*]}

  # creat namelist for surface analysis
  #cp ${NMLDIR}/App_3DVAR_MultiGrid.nl.tmp ${NMLDIR}/App_3DVAR_MultiGrid.nl
  cp ${TOPDIR}/run/App_3DVAR_MultiGrid_WRF2.nl.tmp ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl
  #sed -i.bak "s/INIT_UNIXTIME_SRFANA/${INIT_UNIXTIME_SRFANA}/g" ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl
  #sed -i.bak "s/FINL_UNIXTIME_SRFANA/${FINL_UNIXTIME_SRFANA}/g" ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl
  sed -i.bak "s/sYear/${sYear}/g" ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl
  sed -i.bak "s/sMon/${sMon}/g" ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl
  sed -i.bak "s/sDay/${sDay}/g" ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl
  sed -i.bak "s/sHor/${sHor}/g" ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl
  sed -i.bak "s/sMin/${sMin}/g" ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl
  sed -i.bak "s/sSec/${sSec}/g" ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl
  sed -i.bak "s/eYear/${eYear}/g" ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl
  sed -i.bak "s/eMon/${eMon}/g" ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl
  sed -i.bak "s/eDay/${eDay}/g" ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl
  sed -i.bak "s/eHor/${eHor}/g" ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl
  sed -i.bak "s/eMin/${eMin}/g" ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl
  sed -i.bak "s/eSec/${eSec}/g" ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl
  sed -i.bak "s/NAME_BAKFILE/${NAME_BAKFILE[*]}/g" ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl
  sed -i.bak "s/NAME_OBSFILE/${NAME_OBSFILE[*]}/g" ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl


# 4. RUN EXE
# ==================================================================================
#
  export TIME_ANA=`date "+%Y%m%d_%H%M" -d @${FINL_UNIXTIME_SRFANA}`
  export MON_ANA=`date "+%Y%m" -d @${FINL_UNIXTIME_SRFANA}`
  export DAY_ANA=`date "+%d" -d @${FINL_UNIXTIME_SRFANA}`
  export HM_ANA=`date "+%H%M" -d @${FINL_UNIXTIME_SRFANA}`
  export OUTDIR=${TOPDIR}/grid/output/
  export SAVEDIR=${OUTDIR}/${MON_ANA}/${DAY_ANA}

  # run
  rm log.*
  mpirun -n ${nNum} ./App_3DVAR_MultiGrid_WRF.exe >& ./log_${TIME_ANA}
  #srun hostname | sort -u > hosts
  #mpirun -machinefile hosts -n 4 ./App_3DVAR_MultiGrid.exe >& ./log_${TIME_ANA}

  # move files
  #mkdir -p ${SAVEDIR}
  #mv ${OUTDIR}/bak/bak.nc ${SAVEDIR}/bak_${TIME_ANA}.nc
  #mv ${OUTDIR}/ana/ana.nc ${SAVEDIR}/ana_${TIME_ANA}.nc
  mv ./log_${TIME_ANA} ${OUTDATADIR}/log/log_${TIME_ANA}
  cp ${NMLDIR}/App_3DVAR_MultiGrid_WRF.nl ${OUTDATADIR}/log/nl_${TIME_ANA}

  # check status
  if grep "Time cost" ${OUTDATADIR}/log/log_${TIME_ANA} ; then
    echo "`date` Done ${TIME_ANA} surface analysis" >> ${TOPDIR}/run/Operation2.log
  else
    echo "!!!!!! ERROR `date` surface analysis in ${TIME_ANA} failed ..." >> ${TOPDIR}/run/Operation2.log
  fi







