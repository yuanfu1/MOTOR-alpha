#!/bin/bash
#
#--------------------------------------------------------------
# Script to run MOTOR-DA with MultiDA as an operation system.
#
# Author: Jiongming Pang (pang.j.m@hotmail.com)
#         @GBA-MWF/SIMI, 2022-02-10
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
  # initial time setting
  export INTERVAL_ANA=6          # in minutes, real-time-analysis' interval 
  export INTERVAL_OBSSurface=5   # in minutes, Obs_Surface files' interval
  export INTERVAL_OBSSound=180   # in minutes, Obs_Sound files' interval
  export INTERVAL_OBSVwpw=360    # in minutes, Obs_Vwpw files' interval
  export INTERVAL_BAKFILE=60     # in minutes, background outputs's interval
  export TIME_WINDOW=96          # in minutes
  export TIME_ZONE=8             # in hours
  export TIME_DELAY=0            # in minutes, operation time delay

  # initial directory setting, can be modified by users
  export TOPDIR=/public/home/srf/data/surface_ana/motor/MOTOR/
  export BAKDIR=/public/home/haps/HAPS_ENS/recentWRF3km
  export OBSDIR=/public/publicdata/rawData

  export INDATADIR=${TOPDIR}/input/
  export OUTDATADIR=${TOPDIR}/output/

  # set parts of options in namelist file
  export VLEVEL=60   # vLevel in analysis
  export ANAVAR='temp, uwnd, vwnd'   # variables for analysis
  export TSLOTS3=$(( ${TIME_WINDOW} / ${INTERVAL_ANA} + 1 ))
  export TSLOTS2=$(( ( ${TSLOTS3} - 1 ) / 2 + 1 ))
  export TSLOTS1=$(( ( ${TSLOTS2} - 1 ) / 2 + 1 ))

  cd ${TOPDIR}
  source ./pathEnv.sh
  cd ${TOPDIR}/run

# 2. GET THE TIME WINDOWS OF ANALYSIS, BAKFILES, and OBSFILES
# ==================================================================================
#
  # initial time setting
  tmpA=`date +%s`
  tmpA=$(( ${tmpA} - ${TIME_DELAY} * 60 ))

  # get the time window for analysis
  FINL_UNIXTIME_ANA=$(( $tmpA / $INTERVAL_ANA / 60 ))
  FINL_UNIXTIME_ANA=$(( $FINL_UNIXTIME_ANA * $INTERVAL_ANA * 60 - ${TIME_ZONE} * 60 * 60 ))
  INIT_UNIXTIME_ANA=$(( $FINL_UNIXTIME_ANA - $TIME_WINDOW * 60 ))

  sYear=`date "+%_Y" -d @${INIT_UNIXTIME_ANA}`
  sMon=`date "+%_m" -d @${INIT_UNIXTIME_ANA}`
  sDay=`date "+%_d" -d @${INIT_UNIXTIME_ANA}`
  sHor=`date "+%_H" -d @${INIT_UNIXTIME_ANA}`
  sMin=`date "+%_M" -d @${INIT_UNIXTIME_ANA}`
  sSec=`date "+%_S" -d @${INIT_UNIXTIME_ANA}`

  eYear=`date "+%_Y" -d @${FINL_UNIXTIME_ANA}`
  eMon=`date "+%_m" -d @${FINL_UNIXTIME_ANA}`
  eDay=`date "+%_d" -d @${FINL_UNIXTIME_ANA}`
  eHor=`date "+%_H" -d @${FINL_UNIXTIME_ANA}`
  eMin=`date "+%_M" -d @${FINL_UNIXTIME_ANA}`
  eSec=`date "+%_S" -d @${FINL_UNIXTIME_ANA}`

  # get the time window for background files from model output
  tmp=$(bc <<< "scale=4; ( ( $FINL_UNIXTIME_ANA + ${TIME_ZONE} * 60 * 60 ) / $INTERVAL_BAKFILE / 60 )")
  tmp=`echo ${tmp}|awk '{print int($1)==$1?int($1):int(int($1*10/10+1))}'`
  FINL_UNIXTIME_BAKFILE=$(( $tmp * $INTERVAL_BAKFILE * 60  - ${TIME_ZONE} * 60 * 60 ))
  tmp=$(( ( $INIT_UNIXTIME_ANA + ${TIME_ZONE} * 60 * 60 ) / $INTERVAL_BAKFILE / 60 ))
  INIT_UNIXTIME_BAKFILE=$(( $tmp * $INTERVAL_BAKFILE * 60 - ${TIME_ZONE} * 60 * 60 ))

  # get the time window for obs_surface files
  tmp=$(bc <<< "scale=4; ( ( $FINL_UNIXTIME_ANA + ${TIME_ZONE} * 60 * 60 ) / $INTERVAL_OBSSurface / 60 )")
  tmp=`echo ${tmp}|awk '{print int($1)==$1?int($1):int(int($1*10/10+1))}'`
  FINL_UNIXTIME_OBSSurface=$(( $tmp * $INTERVAL_OBSSurface * 60 - ${TIME_ZONE} * 60 * 60 ))
  tmp=$(( ( $INIT_UNIXTIME_ANA + ${TIME_ZONE} * 60 * 60 ) / $INTERVAL_OBSSurface / 60 ))
  INIT_UNIXTIME_OBSSurface=$(( $tmp * $INTERVAL_OBSSurface * 60 - ${TIME_ZONE} * 60 * 60 ))

  # get the time window for obs_sound files
  tmp=$(bc <<< "scale=4; ( ( $FINL_UNIXTIME_ANA + ${TIME_ZONE} * 60 * 60 ) / $INTERVAL_OBSSound / 60 )")
  tmp=`echo ${tmp}|awk '{print int($1)==$1?int($1):int(int($1*10/10+1))}'`
  FINL_UNIXTIME_OBSSound=$(( $tmp * $INTERVAL_OBSSound * 60 - ${TIME_ZONE} * 60 * 60 ))
  tmp=$(( ( $INIT_UNIXTIME_ANA + ${TIME_ZONE} * 60 * 60 ) / $INTERVAL_OBSSound / 60 ))
  INIT_UNIXTIME_OBSSound=$(( $tmp * $INTERVAL_OBSSound * 60 - ${TIME_ZONE} * 60 * 60 ))

  # get the time window for obs_vwpw files
  tmp=$(bc <<< "scale=4; ( ( $FINL_UNIXTIME_ANA + ${TIME_ZONE} * 60 * 60 ) / $INTERVAL_OBSVwpw / 60 )")
  tmp=`echo ${tmp}|awk '{print int($1)==$1?int($1):int(int($1*10/10+1))}'`
  FINL_UNIXTIME_OBSVwpw=$(( $tmp * $INTERVAL_OBSVwpw * 60 - ${TIME_ZONE} * 60 * 60 ))
  tmp=$(( ( $INIT_UNIXTIME_ANA + ${TIME_ZONE} * 60 * 60 ) / $INTERVAL_OBSVwpw / 60 ))
  INIT_UNIXTIME_OBSVwpw=$(( $tmp * $INTERVAL_OBSVwpw * 60 - ${TIME_ZONE} * 60 * 60 ))


# 3. GET THE NAMES OF ALL BAK FILES AND OBS FILES, AND PREPARE INPUT FILES:
# ==================================================================================
#
  export NMLDIR=${TOPDIR}/static/Application

  # get the names of all bak files
  NUM_BAKFILE=$(( ( $FINL_UNIXTIME_BAKFILE - $INIT_UNIXTIME_BAKFILE ) / $INTERVAL_BAKFILE / 60 + 1 ))
  timetmp=${INIT_UNIXTIME_BAKFILE}
  j=0
  for ((i=1;i<=$NUM_BAKFILE;i++)) ; do
    nametmp=wrfout_d01_`date "+%Y-%m-%d_%H:%M:%S" -d @${timetmp}`
    if [ -e ${BAKDIR}/${nametmp} ] ; then
      j=$(( ${j} + 1 ))
      ln -sf ${BAKDIR}/${nametmp} ${INDATADIR}/
      NAME_BAKFILE[$j]=${nametmp}","
    fi
    timetmp=$(( ${timetmp} + ${INTERVAL_BAKFILE} * 60 ))
  done
  NAME_BAKFILE[$j]=`echo ${NAME_BAKFILE[$j]%,*}`
  #echo ${NAME_BAKFILE[*]}

  # get the names of all obs_surface files
  NUM_OBSSurface=$(( ( $FINL_UNIXTIME_OBSSurface - $INIT_UNIXTIME_OBSSurface ) / $INTERVAL_OBSSurface / 60 + 1 ))
  timetmp=${INIT_UNIXTIME_OBSSurface}
  j=0
  for ((i=1;i<=${NUM_OBSSurface};i++)) ; do
    nametmp=`date "+%Y%m%d_%H%M" -d @${timetmp}`
    if [ -e ${OBSDIR}/AWS/${nametmp} ] ; then
      j=$(( ${j} + 1 ))
      ln -sf ${OBSDIR}/AWS/${nametmp} ${INDATADIR}/
      NAME_OBSFILE_Surface[$j]=${nametmp}","
    fi
    timetmp=$(( ${timetmp} + ${INTERVAL_OBSSurface} * 60 ))
  done
  NAME_OBSFILE_Surface[$j]=`echo ${NAME_OBSFILE_Surface[$j]%,*}`
  #echo ${NAME_OBSFILE_Surface[*]}
  if [ $j = 0 ] ; then
    FlagSurface=0
  else
    FlagSurface=1
  fi

  # get the names of all obs_sound files
  NUM_OBSSound=$(( ( $FINL_UNIXTIME_OBSSound - $INIT_UNIXTIME_OBSSound ) / $INTERVAL_OBSSound / 60 + 1 ))
  timetmp=${INIT_UNIXTIME_OBSSound}
  j=0
  for ((i=1;i<=${NUM_OBSSound};i++)) ; do
    yyyymm=`date "+%Y%m" -d @${timetmp}`
    nametmp=rec_RTEMP_`date "+%Y%m%d%H%M%S" -d @${timetmp}`_g.dat
    if [ -e ${OBSDIR}/rapid/RTEMP/${yyyymm}/${nametmp} ] ; then
      j=$(( ${j} + 1 ))
      ln -sf ${OBSDIR}/rapid/RTEMP/${yyyymm}/${nametmp} ${INDATADIR}/
      NAME_OBSFILE_Sound[$j]=${nametmp}","
    fi
    timetmp=$(( ${timetmp} + ${INTERVAL_OBSSound} * 60 ))
  done
  NAME_OBSFILE_Sound[$j]=`echo ${NAME_OBSFILE_Sound[$j]%,*}`
  #echo ${NAME_OBSFILE_Sound[*]}
  if [ $j = 0 ] ; then
    FlagSound=0
  else
    FlagSound=1
  fi

  # get the names of all obs_vwpw files
  export NUM_OBSVwpw=$(( ( $FINL_UNIXTIME_OBSVwpw - $INIT_UNIXTIME_OBSVwpw ) / $INTERVAL_OBSVwpw / 60 + 1 ))
  export timetmp=${INIT_UNIXTIME_OBSVwpw}
  j=0
  for ((i=1;i<=${NUM_OBSVwpw};i++)) ; do
    yyyymm=`date "+%Y%m" -d @${timetmp}`
    nametmp=vwpw`date "+%y%m%d%H" -d @${timetmp}`.dat
    if [ -e ${OBSDIR}/vwpw/${yyyymm}/${nametmp} ] ; then
      j=$(( ${j} + 1 ))
      ln -sf ${OBSDIR}/vwpw/${yyyymm}/${nametmp} ${INDATADIR}/
      NAME_OBSFILE_Vwpw[$j]=${nametmp}","
    fi
    timetmp=$(( ${timetmp} + ${INTERVAL_OBSVwpw} * 60 ))
  done
  NAME_OBSFILE_Vwpw[$j]=`echo ${NAME_OBSFILE_Vwpw[$j]%,*}`
  #echo ${NAME_OBSFILE_Vwpw[*]}
  if [ $j = 0 ] ; then
    FlagVwpw=0
  else
    FlagVwpw=1
  fi

  # creat namelist for surface analysis
  cp ${TOPDIR}/run/App_3DVAR_MultiGrid_MultiDA_WRF.nl.tmp ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/sYear/${sYear}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/sMon/${sMon}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/sDay/${sDay}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/sHor/${sHor}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/sMin/${sMin}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/sSec/${sSec}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/eYear/${eYear}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/eMon/${eMon}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/eDay/${eDay}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/eHor/${eHor}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/eMin/${eMin}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/eSec/${eSec}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/TSLOTS1/${TSLOTS2}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/TSLOTS2/${TSLOTS3}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/TSLOTS3/${TSLOTS3}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/VLEVEL/${VLEVEL}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/ANAVAR/${ANAVAR}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/NAME_BAKFILE/${NAME_BAKFILE[*]}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/FlagSurface/${FlagSurface}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/FlagSound/${FlagSound}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/FlagVwpw/${FlagVwpw}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/NAME_OBSFILE_Surface/${NAME_OBSFILE_Surface[*]}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/NAME_OBSFILE_Sound/${NAME_OBSFILE_Sound[*]}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl
  sed -i.bak "s/NAME_OBSFILE_Vwpw/${NAME_OBSFILE_Vwpw[*]}/g" ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl


# 4. RUN EXE
# ==================================================================================
#
  TIME_ANA=`date "+%Y%m%d_%H%M" -d @${FINL_UNIXTIME_ANA}`
  MON_ANA=`date "+%Y%m" -d @${FINL_UNIXTIME_ANA}`
  DAY_ANA=`date "+%d" -d @${FINL_UNIXTIME_ANA}`
  HM_ANA=`date "+%H%M" -d @${FINL_UNIXTIME_ANA}`

  export OUTDATADIR=${TOPDIR}/output/
  export infoDIR=${OUTDATADIR}/log/${MON_ANA}/${DAY_ANA}

  # move the files created two days ago
  SAVE_UNIXTIME=$(( ${FINL_UNIXTIME_ANA} - 24 * 60 * 60 ))
  TIME_SAVE=`date "+%Y%m%d" -d @${SAVE_UNIXTIME}`
  MON_SAVE=`date "+%Y%m" -d @${SAVE_UNIXTIME}`
  DAY_SAVE=`date "+%d" -d @${SAVE_UNIXTIME}`
  export anaDIR=${OUTDATADIR}/ana/${MON_SAVE}/${DAY_SAVE}
  export bakDIR=${OUTDATADIR}/bak/${MON_SAVE}/${DAY_SAVE}

  # run
  rm log.*
  #mpirun -n ${nNum} ./App_3DVAR_MultiGrid_WRF.exe >& ./log_${TIME_ANA}
  mpirun -n ${nNum} ./App_3DVAR_MultiGrid_WRF_MultiObs.exe >& ./log_${TIME_ANA}

  # move info-files
  mkdir -p ${infoDIR}
  mv ./log_${TIME_ANA} ${infoDIR}/log_${TIME_ANA}
  cp ${NMLDIR}/App_3DVAR_MultiDA_WRF.nl ${infoDIR}/nl_${TIME_ANA}

  # restore files
  mkdir -p ${anaDIR}
  mkdir -p ${bakDIR}
  mv ${OUTDATADIR}/ana_${TIME_SAVE}*.nc ${anaDIR}/
  mv ${OUTDATADIR}/bak_${TIME_SAVE}*.nc ${bakDIR}/
  

  # check status
  if grep "Time cost" ${infoDIR}/log_${TIME_ANA} ; then
    echo "`date` Done ${TIME_ANA} surface analysis" >> ${TOPDIR}/run/Operation.log
    if [ $FlagSurface = 1 ] ; then
      echo "use Obs_Surface" >> ${TOPDIR}/run/Operation.log
    fi
    if [ $FlagSound = 1 ] ; then
      echo "use Obs_Sound" >> ${TOPDIR}/run/Operation.log
    fi
    if [ $FlagVwpw = 1 ] ; then
      echo "use Obs_Vwpw" >> ${TOPDIR}/run/Operation.log
    fi
  else
    echo "!!!!!! ERROR `date` surface analysis in ${TIME_ANA} failed ..." >> ${TOPDIR}/run/Operation.log
  fi







