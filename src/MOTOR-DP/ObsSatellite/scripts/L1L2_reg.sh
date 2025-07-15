
#rm ./*HDF
#rm ./*NC
ROOTDIR=`pwd`
for start_date in 2022052623 2022052700
#2021050408 2021050409 2021050410 2021050411 2021050412 2021050413 2021050414 2021050415 2021050416
do

  ss1=`${BIN_DIR}/da_advance_time.exe ${start_date} -1 -f ccyymmddhh`
  ss2=`${BIN_DIR}/da_advance_time.exe ${start_date} +0 -f ccyymmddhh`
  start_date_curr=`${BIN_DIR}/da_advance_time.exe ${ss1} +15min -f ccyymmddhhnnss`
  
  while [ ${ss1} -le $ss2 ]; do
    
    end_date_curr=`${BIN_DIR}/da_advance_time.exe ${start_date_curr} +4min17s -f ccyymmddhhnnss`
    echo 'start_date=' $start_date_curr
    echo 'end_date=' $end_date_curr
    echo 'ss1=' $ss1
    echo 'ss2=' $ss2
    yyyy=`echo $start_date_curr|cut -c 1-4`
    yyyymmdd=`echo $start_date_curr|cut -c 1-8`
    hh=`echo $start_date_curr|cut -c 9-10`
    yyyymmddhh=`echo $start_date_curr|cut -c 1-10`
    if [ $hh -gt 15 ]; then
      date_path=`${BIN_DIR}/da_advance_time.exe ${start_date_curr} +1d -f ccyymmddhhnnss`
      yyyymmdd_path=`echo $date_path|cut -c 1-8`
    else
      yyyymmdd_path=$yyyymmdd
    fi
    echo "WORK ON " $hh $yyyymmdd_path $start_date_curr $end_date_curr
    
    OBSDIR=/public/publicdata/fy4a/
    cd $OBSDIR/SATE_FY4A_L1/REG/4KM/$yyyymmdd_path/
    #cp Z_SATE*FDI-_MULT_NOM_${start_date_curr}_${end_date_curr}_4000M_V0001.HDF  ${ROOTDIR}/
    cp Z_SATE*FDI-_MULT_NOM_${ss1}*_4000M_V0001.HDF  ${ROOTDIR}/
    cd $OBSDIR/SATE_FY4A_L1/REG/4KMGEO/$yyyymmdd_path/
    #cp Z_SATE*GEO-_MULT_NOM_${start_date_curr}_${end_date_curr}_4000M_V0001.HDF ${ROOTDIR}/
    cp Z_SATE*GEO-_MULT_NOM_${ss1}*_4000M_V0001.HDF ${ROOTDIR}/
    cd $OBSDIR/SATE_FY4A_L2/REG/CLM/$yyyymmdd_path/
    #cp Z_SATE*CLM-_MULT_NOM_${start_date_curr}_${end_date_curr}_4000M_V0001.NC  ${ROOTDIR}/
    cp Z_SATE*CLM-_MULT_NOM_${ss1}*_4000M_V0001.NC  ${ROOTDIR}/
    
    cd $ROOTDIR
    
    echo "Renaming..."
    # rename
    fnames=`ls Z_SATE*FDI-_MULT_NOM_${yyyymmddhh}*_4000M_V0001.HDF`
    for fname in $fnames
    do
      tstring=${fname##*NOM}
      mv $fname L1_REG_FDI${tstring}
    done
    #
    # rename
    fnames=`ls Z_SATE*GEO-_MULT_NOM_${yyyymmddhh}*_4000M_V0001.HDF`
    for fname in $fnames
    do
      tstring=${fname##*NOM}
      mv $fname L1_REG_GEO${tstring}
    done
    #
    # rename
    fnames=`ls Z_SATE*CLM-_MULT_NOM_${yyyymmddhh}*_4000M_V0001.NC`
    for fname in $fnames
    do
      tstring=${fname##*NOM}
      mv $fname L2_REG_CLM${tstring}
    done
    #
    
    start_date_curr=`${BIN_DIR}/da_advance_time.exe ${start_date_curr} +4min18s -f ccyymmddhhnnss`
    ss1=`${BIN_DIR}/da_advance_time.exe ${start_date_curr} +4min18s -f ccyymmddhhnnss`
  done
#start_date=`${BIN_DIR}/da_advance_time.exe ${start_date} +24`
done

