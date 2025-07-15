for start_date in 2022052623

do

yyyymmddhh=`echo $start_date|cut -c 1-10`
end_date=`${BIN_DIR}/da_advance_time.exe ${start_date} +1 -f ccyymmddhh`

echo "WORK ON " $hh $yyyymmdd_path $start_date_curr $end_date_curr

echo "Cat filenames to txt ..."

fnames_fdi=`ls L1_REG_FDI_${yyyymmddhh}*_4000M_V0001.HDF L1_REG_FDI_${end_date}*_4000M_V0001.HDF L1_GLB_FDI_${yyyymmddhh}*_4000M_V0001.HDF L1_GLB_FDI_${end_date}*_4000M_V0001.HDF`
ls $fnames_fdi &>FDI.txt
#len_fdi=${#fnames[@]}
size_fdi=`ls $fnames_fdi |wc -l`

fnames_geo=`ls L1_REG_GEO_${yyyymmddhh}*_4000M_V0001.HDF L1_REG_GEO_${end_date}*_4000M_V0001.HDF L1_GLB_GEO_${yyyymmddhh}*_4000M_V0001.HDF L1_GLB_GEO_${end_date}*_4000M_V0001.HDF`
ls $fnames_geo &>GEO.txt
size_geo=`ls $fnames_geo |wc -l`

fnames_clm=`ls L2_REG_CLM_${yyyymmddhh}*_4000M_V0001.NC L2_REG_CLM_${end_date}*_4000M_V0001.NC L2_GLB_CLM_${yyyymmddhh}*_4000M_V0001.NC L2_GLB_CLM_${end_date}*_4000M_V0001.NC`
ls $fnames_clm &>CLM.txt
size_clm=`ls $fnames_clm |wc -l`

echo "LENS = " $size_fdi $size_geo $size_clm

# count=0
# while [ $count -le $size_fdi ]; do 
#   echo $count
#   tstring0=${fnames_clm[$count]}
#   echo $tstring0
#   tstring1=${tstring0##*CLM}
#   tstring=${tstring1%NC*}
#   echo $tstring 
#   count=`expr $count + 1`
# done

done

