
#rm ./*HDF
#rm ./*NC
ROOTDIR=`pwd`
for start_date in 2022052623 2022052700
#2021050408 2021050409 2021050410 2021050411 2021050412 2021050413 2021050414 2021050415 2021050416
do

ss1=`${BIN_DIR}/da_advance_time.exe ${start_date} -3 -f ccyymmddhhnn`
ss2=`${BIN_DIR}/da_advance_time.exe ${start_date} +7 -f ccyymmddhhnn`

while [ ${ss1} -le $ss2 ]; do
start_date_curr=`${BIN_DIR}/da_advance_time.exe ${ss1} +0min -f ccyymmddhhnnss`
end_date_curr=`${BIN_DIR}/da_advance_time.exe ${start_date_curr} +14min59s -f ccyymmddhhnnss`
echo 'start_date=' $start_date_curr
echo 'end_date=' $end_date_curr
yyyy=`echo $start_date_curr|cut -c 1-4`
yyyymmdd=`echo $start_date_curr|cut -c 1-8`
hh=`echo $start_date_curr|cut -c 9-10`
if [ $hh -gt 15 ]; then
  date_path=`${BIN_DIR}/da_advance_time.exe ${start_date_curr} +1d -f ccyymmddhhnnss`
  yyyymmdd_path=`echo $date_path|cut -c 1-8`
else
  yyyymmdd_path=$yyyymmdd
fi
echo "WORK ON " $hh $yyyymmdd_path $start_date_curr $end_date_curr

OBSDIR=/public/publicdata/fy4a/
cd $OBSDIR/SATE_FY4A_L1/GLB/4KM/$yyyymmdd_path/
cp Z_SATE*FDI-_MULT_NOM_${start_date_curr}_${end_date_curr}_4000M_V0001.HDF  ${ROOTDIR}/
cd $OBSDIR/SATE_FY4A_L1/GLB/4KMGEO/$yyyymmdd_path/
cp Z_SATE*GEO-_MULT_NOM_${start_date_curr}_${end_date_curr}_4000M_V0001.HDF ${ROOTDIR}/
cd $OBSDIR/SATE_FY4A_L2/GLB/CLM/$yyyymmdd_path/
cp Z_SATE*CLM-_MULT_NOM_${start_date_curr}_${end_date_curr}_4000M_V0001.NC  ${ROOTDIR}/

cd $ROOTDIR

# rename
fnames=`ls Z_SATE*FDI-_MULT_NOM_${start_date_curr}_${end_date_curr}_4000M_V0001.HDF`
for fname in $fnames
do
  tstring=${fname##*NOM}
  mv $fname L1_GLB_FDI${tstring}
done
#
# rename
fnames=`ls Z_SATE*GEO-_MULT_NOM_${start_date_curr}_${end_date_curr}_4000M_V0001.HDF`
for fname in $fnames
do
  tstring=${fname##*NOM}
  mv $fname L1_GLB_GEO${tstring}
done
#
# rename
fnames=`ls Z_SATE*CLM-_MULT_NOM_${start_date_curr}_${end_date_curr}_4000M_V0001.NC`
for fname in $fnames
do
  tstring=${fname##*NOM}
  mv $fname L2_GLB_CLM${tstring}
done
#

ss1=`${BIN_DIR}/da_advance_time.exe ${ss1} +15min -f ccyymmddhhnn`
done
#start_date=`${BIN_DIR}/da_advance_time.exe ${start_date} +24`
done

