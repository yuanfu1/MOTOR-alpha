! Description:
!> @file
!!   Subroutines for University of Wisconsin IR emissivity atlas
!
!> @brief
!!   Subroutines for University of Wisconsin IR emissivity atlas
!!
!! @details
!!   It is intended that this atlas be used via the RTTOV interface
!!   rather than by calling these subroutines directly.
!!
!!   Borbas, E. E. and B. C. Ruston, 2010.
!!   The RTTOV UWiremis IR land surface emissivity module.
!!   NWP SAF report: http://nwpsaf.eu/vs_reports/nwpsaf-mo-vs-042.pdf
!!
!!   Borbas, E, 2014.
!!   The RTTOV UWiremis module Investigation into the angular dependence of IR surface emissivity.
!!   NWP SAF report: http://nwpsaf.eu/vs_reports/nwpsaf-mo-vs-050.pdf
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 7 December 2016, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, DWD and MeteoFrance.
!
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
MODULE mod_uwiremis_atlas

  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0      02/06/2010  Based on UW IR atlas code (E. Borbas, B. Ruston, J. Hocking)
  !  1.1      02/06/2014  satellite zenith angle correction added (E. Borbas)

#include "throw.h"

#ifdef _RTTOV_HDF
#undef _RTTOV_NETCDF
#endif

  USE parkind1, ONLY : &
    jpim, &   ! 32-bit int
    jpis, &   ! 16-bit
    jpit, &   ! 8-bit
    jprb, &   ! 64-bit real
    jprm, &   ! 32-bit
    jplm      ! logical

#ifdef _RTTOV_HDF
  USE rttov_hdf_mod, ONLY : &
    open_hdf,             &
    close_hdf,            &
    is_hdf_open,          &
    is_hdf_64bit_reals
#endif

  IMPLICIT NONE

#include "rttov_errorreport.interface"

#ifdef _RTTOV_NETCDF
  INCLUDE 'netcdf.inc'
#endif

  PRIVATE
  PUBLIC :: uwiremis_atlas_data,        &
            rttov_uwiremis_init,        &
            rttov_uwiremis,             &
            rttov_uwiremis_close_atlas, &
            bfemis_gridres, numpcs      ! Useful for external programs

  ! Atlas constants

  INTEGER(KIND=jpim), PARAMETER :: numpcs = 6
  INTEGER(KIND=jpim), PARAMETER :: numwave = 416

  INTEGER(KIND=jpim), PARAMETER :: db_ver_year = 2007

  INTEGER(KIND=jpim), PARAMETER :: seaice_flag = 70           ! flag value returned for sea-ice
  REAL(KIND=jprb),    PARAMETER :: default_std = 0.05_jprb    ! default standard deviation

  INTEGER(KIND=jpim), PARAMETER :: bfemis_gridres = 100       ! 0.1 deg
  INTEGER(KIND=jpim), PARAMETER :: bfemis_ygrid1 = 89950      ! 89.95 deg
  INTEGER(KIND=jpim), PARAMETER :: bfemis_xgrid1 = -179950    ! -179.95 deg

  INTEGER(KIND=jpim), PARAMETER :: cov_emis_gridres = 500     ! 0.5 deg
  INTEGER(KIND=jpim), PARAMETER :: cov_emis_ygrid1 = 89750    ! 89.75 deg
  INTEGER(KIND=jpim), PARAMETER :: cov_emis_xgrid1 = -179750  ! -179.75 deg

  INTEGER(KIND=jpim), PARAMETER :: igbp_gridres = 50          ! 0.05 deg
  INTEGER(KIND=jpim), PARAMETER :: igbp_ygrid1 = 89950        ! 89.95 deg
  INTEGER(KIND=jpim), PARAMETER :: igbp_xgrid1 = -179950      ! -179.95 deg

  REAL(KIND=jprb),    PARAMETER :: angcorrminzen = 5._jprb
  REAL(KIND=jprb),    PARAMETER :: angcorrterm = 85._jprb     ! Solar zenith angle of terminator

  REAL(KIND=jprb),    PARAMETER :: hkod = -999._jprb          ! Missing data


  !> Data type for UWIRemis atlas data
  TYPE uwiremis_atlas_data
    PRIVATE

    ! Flags
    LOGICAL(KIND=jplm) :: single_inst
    LOGICAL(KIND=jplm) :: std_init
    LOGICAL(KIND=jplm) :: do_ang_corr

    ! Dimensions
    INTEGER(KIND=jpim) :: nb_lats
    INTEGER(KIND=jpim) :: nb_lons
    INTEGER(KIND=jpim) :: nb_pack

    INTEGER(KIND=jpim) :: cv_lats
    INTEGER(KIND=jpim) :: cv_lons
    INTEGER(KIND=jpim) :: cv_pack

    INTEGER(KIND=jpim) :: igbp_lats
    INTEGER(KIND=jpim) :: igbp_lons
    INTEGER(KIND=jpim) :: nb_igbp

    ! Atlas data loaded by initialisation routine
    INTEGER(KIND=jpis), POINTER :: bfemis_flag(:,:)  ! dims are (nb_lats,nb_lons)
    INTEGER(KIND=jpim), POINTER :: bfemis_lut(:,:)   ! dims are (nb_lats,nb_lons)
    REAL(KIND=jprm),    POINTER :: pca_coef(:,:)     ! dims are (nb_pack,numpcs)

    INTEGER(KIND=jpim), POINTER :: cov_emis_lut(:,:) ! dims are (cv_lats,cv_lons)
    INTEGER(KIND=jpis), POINTER :: cov_emis(:,:)     ! dims are (cv_pack,numwave)

    INTEGER(KIND=jpit), POINTER :: igbp(:,:)         ! dims are (igbp_lats,igbp_lons)

    REAL(KIND=jprm),    POINTER :: p1d(:,:)          ! dims are (numwave,nb_igbp)
    REAL(KIND=jprm),    POINTER :: p2d(:,:)          ! dims are (numwave,nb_igbp)
    REAL(KIND=jprm),    POINTER :: p3d(:,:)          ! dims are (numwave,nb_igbp)
    REAL(KIND=jprm),    POINTER :: p1n(:,:)          ! dims are (numwave,nb_igbp)
    REAL(KIND=jprm),    POINTER :: p2n(:,:)          ! dims are (numwave,nb_igbp)
    REAL(KIND=jprm),    POINTER :: p3n(:,:)          ! dims are (numwave,nb_igbp)

    REAL(KIND=jprm),    POINTER :: pcu(:,:)          ! pcu need to be read in from flat file

    ! Data to allow more efficient memory usage (convert jprm to jprb only at point of use)
    REAL(KIND=jprm),    POINTER :: pca_sfac(:)       ! PCA coef scale factors
    REAL(KIND=jprm),    POINTER :: pca_offs(:)       ! PCA coef offsets
    REAL(KIND=jprm) :: cov_sfac                      ! Stdev scale factor

    ! Arrays to hold hsr data interpolated onto channel wavenumbers
    INTEGER(KIND=jpim) :: platform_id
    INTEGER(KIND=jpim) :: sat_id
    INTEGER(KIND=jpim) :: inst_id
    INTEGER(KIND=jpim) :: ncoefchans                 ! Number of channels in coef file
    REAL(KIND=jprb), POINTER :: cov_emis_int(:,:)    ! dims are (cv_pack,nchannels)
    REAL(KIND=jprb), POINTER :: pcu_int(:,:,:)       ! dims are (numpcs,nchannels,1/2)
    REAL(KIND=jprb), POINTER :: pcm_int(:,:)
    REAL(KIND=jprb), POINTER :: sice_em_int(:), snow_em_int(:)
    REAL(KIND=jprm), POINTER :: p1d_int(:,:,:)
    REAL(KIND=jprm), POINTER :: p2d_int(:,:,:)
    REAL(KIND=jprm), POINTER :: p3d_int(:,:,:)
    REAL(KIND=jprm), POINTER :: p1n_int(:,:,:)
    REAL(KIND=jprm), POINTER :: p2n_int(:,:,:)
    REAL(KIND=jprm), POINTER :: p3n_int(:,:,:)
  END TYPE uwiremis_atlas_data


  ! Atlas data contained in this file - shared, does not change

  REAL(KIND=jprb) :: sice_em(numwave), snow_em(numwave)
  REAL(KIND=jprb) :: sice_stdv, snow_stdv
  REAL(KIND=jprb) :: pcm(numwave)
  REAL(KIND=jprb) :: hsr_wavenum(numwave)

  DATA hsr_wavenum / &
    699.3_jprb,  704.3_jprb,  709.3_jprb,  714.3_jprb,  719.3_jprb,  724.3_jprb,  729.3_jprb,  734.3_jprb,  &
    739.3_jprb,  744.3_jprb,  749.3_jprb,  754.3_jprb,  759.3_jprb,  764.3_jprb,  769.3_jprb,  774.3_jprb,  &
    779.3_jprb,  784.3_jprb,  789.3_jprb,  794.3_jprb,  799.3_jprb,  804.3_jprb,  809.3_jprb,  814.3_jprb,  &
    819.3_jprb,  824.3_jprb,  829.3_jprb,  834.3_jprb,  839.3_jprb,  844.3_jprb,  849.3_jprb,  854.3_jprb,  &
    859.3_jprb,  864.3_jprb,  869.3_jprb,  874.3_jprb,  879.3_jprb,  884.3_jprb,  889.3_jprb,  894.3_jprb,  &
    899.3_jprb,  904.3_jprb,  909.3_jprb,  914.3_jprb,  919.3_jprb,  924.3_jprb,  929.3_jprb,  934.3_jprb,  &
    939.3_jprb,  944.3_jprb,  949.3_jprb,  954.3_jprb,  959.3_jprb,  964.3_jprb,  969.3_jprb,  974.3_jprb,  &
    979.3_jprb,  984.3_jprb,  989.3_jprb,  994.3_jprb,  999.3_jprb, 1004.3_jprb, 1009.3_jprb, 1014.3_jprb,  &
    1019.3_jprb, 1024.3_jprb, 1029.3_jprb, 1034.3_jprb, 1039.3_jprb, 1044.3_jprb, 1049.3_jprb, 1054.3_jprb,  &
    1059.3_jprb, 1064.3_jprb, 1069.3_jprb, 1074.3_jprb, 1079.3_jprb, 1084.3_jprb, 1089.3_jprb, 1094.3_jprb,  &
    1099.3_jprb, 1104.3_jprb, 1109.3_jprb, 1114.3_jprb, 1119.3_jprb, 1124.3_jprb, 1129.3_jprb, 1134.3_jprb,  &
    1139.3_jprb, 1144.3_jprb, 1149.3_jprb, 1154.3_jprb, 1159.3_jprb, 1164.3_jprb, 1169.3_jprb, 1174.3_jprb,  &
    1179.3_jprb, 1184.3_jprb, 1189.3_jprb, 1194.3_jprb, 1199.3_jprb, 1204.3_jprb, 1209.3_jprb, 1214.3_jprb,  &
    1219.3_jprb, 1224.3_jprb, 1229.3_jprb, 1234.3_jprb, 1239.3_jprb, 1244.3_jprb, 1249.3_jprb, 1254.3_jprb,  &
    1259.3_jprb, 1264.3_jprb, 1269.3_jprb, 1274.3_jprb, 1279.3_jprb, 1284.3_jprb, 1289.3_jprb, 1294.3_jprb,  &
    1299.3_jprb, 1304.3_jprb, 1309.3_jprb, 1314.3_jprb, 1319.3_jprb, 1324.3_jprb, 1329.3_jprb, 1334.3_jprb,  &
    1339.3_jprb, 1344.3_jprb, 1349.3_jprb, 1354.3_jprb, 1359.3_jprb, 1364.3_jprb, 1369.3_jprb, 1374.3_jprb,  &
    1379.3_jprb, 1384.3_jprb, 1389.3_jprb, 1394.3_jprb, 1399.3_jprb, 1404.3_jprb, 1409.3_jprb, 1414.3_jprb,  &
    1419.3_jprb, 1424.3_jprb, 1429.3_jprb, 1434.3_jprb, 1439.3_jprb, 1444.3_jprb, 1449.3_jprb, 1454.3_jprb,  &
    1459.3_jprb, 1464.3_jprb, 1469.3_jprb, 1474.3_jprb, 1479.3_jprb, 1484.3_jprb, 1489.3_jprb, 1494.3_jprb,  &
    1499.3_jprb, 1504.3_jprb, 1509.3_jprb, 1514.3_jprb, 1519.3_jprb, 1524.3_jprb, 1529.3_jprb, 1534.3_jprb,  &
    1539.3_jprb, 1544.3_jprb, 1549.3_jprb, 1554.3_jprb, 1559.3_jprb, 1564.3_jprb, 1569.3_jprb, 1574.3_jprb,  &
    1579.3_jprb, 1584.3_jprb, 1589.3_jprb, 1594.3_jprb, 1599.3_jprb, 1604.3_jprb, 1609.3_jprb, 1614.3_jprb,  &
    1619.3_jprb, 1624.3_jprb, 1629.3_jprb, 1634.3_jprb, 1639.3_jprb, 1644.3_jprb, 1649.3_jprb, 1654.3_jprb,  &
    1659.3_jprb, 1664.3_jprb, 1669.3_jprb, 1674.3_jprb, 1679.3_jprb, 1684.3_jprb, 1689.3_jprb, 1694.3_jprb,  &
    1699.3_jprb, 1704.3_jprb, 1709.3_jprb, 1714.3_jprb, 1719.3_jprb, 1724.3_jprb, 1729.3_jprb, 1734.3_jprb,  &
    1739.3_jprb, 1744.3_jprb, 1749.3_jprb, 1754.3_jprb, 1759.3_jprb, 1764.3_jprb, 1769.3_jprb, 1774.3_jprb,  &
    1779.3_jprb, 1784.3_jprb, 1789.3_jprb, 1794.3_jprb, 1799.3_jprb, 1804.3_jprb, 1809.3_jprb, 1814.3_jprb,  &
    1819.3_jprb, 1824.3_jprb, 1829.3_jprb, 1834.3_jprb, 1839.3_jprb, 1844.3_jprb, 1849.3_jprb, 1854.3_jprb,  &
    1859.3_jprb, 1864.3_jprb, 1869.3_jprb, 1874.3_jprb, 1879.3_jprb, 1884.3_jprb, 1889.3_jprb, 1894.3_jprb,  &
    1899.3_jprb, 1904.3_jprb, 1909.3_jprb, 1914.3_jprb, 1919.3_jprb, 1924.3_jprb, 1929.3_jprb, 1934.3_jprb,  &
    1939.3_jprb, 1944.3_jprb, 1949.3_jprb, 1954.3_jprb, 1959.3_jprb, 1964.3_jprb, 1969.3_jprb, 1974.3_jprb,  &
    1979.3_jprb, 1984.3_jprb, 1989.3_jprb, 1994.3_jprb, 1999.3_jprb, 2004.3_jprb, 2009.3_jprb, 2014.3_jprb,  &
    2019.3_jprb, 2024.3_jprb, 2029.3_jprb, 2034.3_jprb, 2039.3_jprb, 2044.3_jprb, 2049.3_jprb, 2054.3_jprb,  &
    2059.3_jprb, 2064.3_jprb, 2069.3_jprb, 2074.3_jprb, 2079.3_jprb, 2084.3_jprb, 2089.3_jprb, 2094.3_jprb,  &
    2099.3_jprb, 2104.3_jprb, 2109.3_jprb, 2114.3_jprb, 2119.3_jprb, 2124.3_jprb, 2129.3_jprb, 2134.3_jprb,  &
    2139.3_jprb, 2144.3_jprb, 2149.3_jprb, 2154.3_jprb, 2159.3_jprb, 2164.3_jprb, 2169.3_jprb, 2174.3_jprb,  &
    2179.3_jprb, 2184.3_jprb, 2189.3_jprb, 2194.3_jprb, 2199.3_jprb, 2204.3_jprb, 2209.3_jprb, 2214.3_jprb,  &
    2219.3_jprb, 2224.3_jprb, 2229.3_jprb, 2234.3_jprb, 2239.3_jprb, 2244.3_jprb, 2249.3_jprb, 2254.3_jprb,  &
    2259.3_jprb, 2264.3_jprb, 2269.3_jprb, 2274.3_jprb, 2279.3_jprb, 2284.3_jprb, 2289.3_jprb, 2294.3_jprb,  &
    2299.3_jprb, 2304.3_jprb, 2309.3_jprb, 2314.3_jprb, 2319.3_jprb, 2324.3_jprb, 2329.3_jprb, 2334.3_jprb,  &
    2339.3_jprb, 2344.3_jprb, 2349.3_jprb, 2354.3_jprb, 2359.3_jprb, 2364.3_jprb, 2369.3_jprb, 2374.3_jprb,  &
    2379.3_jprb, 2384.3_jprb, 2389.3_jprb, 2394.3_jprb, 2399.3_jprb, 2404.3_jprb, 2409.3_jprb, 2414.3_jprb,  &
    2419.3_jprb, 2424.3_jprb, 2429.3_jprb, 2434.3_jprb, 2439.3_jprb, 2444.3_jprb, 2449.3_jprb, 2454.3_jprb,  &
    2459.3_jprb, 2464.3_jprb, 2469.3_jprb, 2474.3_jprb, 2479.3_jprb, 2484.3_jprb, 2489.3_jprb, 2494.3_jprb,  &
    2499.3_jprb, 2504.3_jprb, 2509.3_jprb, 2514.3_jprb, 2519.3_jprb, 2524.3_jprb, 2529.3_jprb, 2534.3_jprb,  &
    2539.3_jprb, 2544.3_jprb, 2549.3_jprb, 2554.3_jprb, 2559.3_jprb, 2564.3_jprb, 2569.3_jprb, 2574.3_jprb,  &
    2579.3_jprb, 2584.3_jprb, 2589.3_jprb, 2594.3_jprb, 2599.3_jprb, 2604.3_jprb, 2609.3_jprb, 2614.3_jprb,  &
    2619.3_jprb, 2624.3_jprb, 2629.3_jprb, 2634.3_jprb, 2639.3_jprb, 2644.3_jprb, 2649.3_jprb, 2654.3_jprb,  &
    2659.3_jprb, 2664.3_jprb, 2669.3_jprb, 2674.3_jprb, 2679.3_jprb, 2684.3_jprb, 2689.3_jprb, 2694.3_jprb,  &
    2699.3_jprb, 2704.3_jprb, 2709.3_jprb, 2714.3_jprb, 2719.3_jprb, 2724.3_jprb, 2729.3_jprb, 2734.3_jprb,  &
    2739.3_jprb, 2744.3_jprb, 2749.3_jprb, 2754.3_jprb, 2759.3_jprb, 2764.3_jprb, 2769.3_jprb, 2774.3_jprb /

  DATA pcm / &
        0.9782182_jprb,   0.9770744_jprb,   0.9763290_jprb,   0.9763215_jprb,   0.9760258_jprb,  &
        0.9763704_jprb,   0.9767076_jprb,   0.9763077_jprb,   0.9758835_jprb,   0.9753462_jprb,  &
        0.9748067_jprb,   0.9734465_jprb,   0.9721510_jprb,   0.9717180_jprb,   0.9714773_jprb,  &
        0.9706340_jprb,   0.9710826_jprb,   0.9722888_jprb,   0.9731166_jprb,   0.9732918_jprb,  &
        0.9736975_jprb,   0.9751787_jprb,   0.9770049_jprb,   0.9773170_jprb,   0.9765164_jprb,  &
        0.9759824_jprb,   0.9750199_jprb,   0.9746831_jprb,   0.9738413_jprb,   0.9731615_jprb,  &
        0.9720387_jprb,   0.9716908_jprb,   0.9708628_jprb,   0.9705366_jprb,   0.9697853_jprb,  &
        0.9694459_jprb,   0.9688896_jprb,   0.9688236_jprb,   0.9689180_jprb,   0.9692774_jprb,  &
        0.9693237_jprb,   0.9692513_jprb,   0.9689918_jprb,   0.9686664_jprb,   0.9684489_jprb,  &
        0.9681804_jprb,   0.9672847_jprb,   0.9667084_jprb,   0.9661347_jprb,   0.9655386_jprb,  &
        0.9650131_jprb,   0.9641176_jprb,   0.9628995_jprb,   0.9620982_jprb,   0.9605948_jprb,  &
        0.9590283_jprb,   0.9572537_jprb,   0.9552648_jprb,   0.9529146_jprb,   0.9505763_jprb,  &
        0.9486620_jprb,   0.9468448_jprb,   0.9446425_jprb,   0.9428397_jprb,   0.9415421_jprb,  &
        0.9398234_jprb,   0.9378662_jprb,   0.9358756_jprb,   0.9338515_jprb,   0.9317511_jprb,  &
        0.9296144_jprb,   0.9274116_jprb,   0.9248639_jprb,   0.9219664_jprb,   0.9197029_jprb,  &
        0.9187206_jprb,   0.9195539_jprb,   0.9211251_jprb,   0.9227578_jprb,   0.9242273_jprb,  &
        0.9256495_jprb,   0.9265392_jprb,   0.9276078_jprb,   0.9279289_jprb,   0.9282181_jprb,  &
        0.9284544_jprb,   0.9289097_jprb,   0.9299400_jprb,   0.9314128_jprb,   0.9329405_jprb,  &
        0.9349486_jprb,   0.9377099_jprb,   0.9380918_jprb,   0.9354525_jprb,   0.9330018_jprb,  &
        0.9316696_jprb,   0.9308965_jprb,   0.9296793_jprb,   0.9282659_jprb,   0.9273711_jprb,  &
        0.9268156_jprb,   0.9265846_jprb,   0.9264724_jprb,   0.9278417_jprb,   0.9298262_jprb,  &
        0.9342009_jprb,   0.9397170_jprb,   0.9451398_jprb,   0.9501663_jprb,   0.9547508_jprb,  &
        0.9586911_jprb,   0.9618842_jprb,   0.9649577_jprb,   0.9675525_jprb,   0.9696881_jprb,  &
        0.9708689_jprb,   0.9717879_jprb,   0.9722518_jprb,   0.9724457_jprb,   0.9728941_jprb,  &
        0.9731293_jprb,   0.9731925_jprb,   0.9730867_jprb,   0.9733831_jprb,   0.9735166_jprb,  &
        0.9740434_jprb,   0.9742066_jprb,   0.9746855_jprb,   0.9748268_jprb,   0.9749292_jprb,  &
        0.9751188_jprb,   0.9752902_jprb,   0.9751062_jprb,   0.9751985_jprb,   0.9752622_jprb,  &
        0.9750626_jprb,   0.9755121_jprb,   0.9755228_jprb,   0.9760818_jprb,   0.9759580_jprb,  &
        0.9758280_jprb,   0.9755163_jprb,   0.9754220_jprb,   0.9750829_jprb,   0.9743836_jprb,  &
        0.9745844_jprb,   0.9742978_jprb,   0.9740397_jprb,   0.9744191_jprb,   0.9745796_jprb,  &
        0.9749123_jprb,   0.9750853_jprb,   0.9746974_jprb,   0.9747824_jprb,   0.9746920_jprb,  &
        0.9735873_jprb,   0.9733123_jprb,   0.9725510_jprb,   0.9718717_jprb,   0.9713586_jprb,  &
        0.9706160_jprb,   0.9701124_jprb,   0.9698699_jprb,   0.9698430_jprb,   0.9694992_jprb,  &
        0.9691019_jprb,   0.9690002_jprb,   0.9678345_jprb,   0.9668854_jprb,   0.9659764_jprb,  &
        0.9666998_jprb,   0.9669611_jprb,   0.9665817_jprb,   0.9679645_jprb,   0.9695909_jprb,  &
        0.9711555_jprb,   0.9724632_jprb,   0.9737635_jprb,   0.9746142_jprb,   0.9748497_jprb,  &
        0.9752109_jprb,   0.9752749_jprb,   0.9754022_jprb,   0.9753313_jprb,   0.9746057_jprb,  &
        0.9745884_jprb,   0.9747860_jprb,   0.9752877_jprb,   0.9753085_jprb,   0.9759305_jprb,  &
        0.9752344_jprb,   0.9748027_jprb,   0.9757417_jprb,   0.9751943_jprb,   0.9748128_jprb,  &
        0.9743713_jprb,   0.9741939_jprb,   0.9725359_jprb,   0.9723988_jprb,   0.9716700_jprb,  &
        0.9708291_jprb,   0.9705051_jprb,   0.9699901_jprb,   0.9689955_jprb,   0.9683419_jprb,  &
        0.9684200_jprb,   0.9672046_jprb,   0.9660766_jprb,   0.9658424_jprb,   0.9648336_jprb,  &
        0.9640325_jprb,   0.9642861_jprb,   0.9636880_jprb,   0.9638920_jprb,   0.9638573_jprb,  &
        0.9641714_jprb,   0.9648057_jprb,   0.9648220_jprb,   0.9639065_jprb,   0.9635883_jprb,  &
        0.9626419_jprb,   0.9616417_jprb,   0.9600965_jprb,   0.9587714_jprb,   0.9576451_jprb,  &
        0.9557189_jprb,   0.9545730_jprb,   0.9550443_jprb,   0.9551759_jprb,   0.9560625_jprb,  &
        0.9576327_jprb,   0.9587138_jprb,   0.9594474_jprb,   0.9598546_jprb,   0.9601094_jprb,  &
        0.9601356_jprb,   0.9597549_jprb,   0.9590299_jprb,   0.9581512_jprb,   0.9572046_jprb,  &
        0.9557602_jprb,   0.9538486_jprb,   0.9521495_jprb,   0.9503905_jprb,   0.9491790_jprb,  &
        0.9485527_jprb,   0.9479896_jprb,   0.9475234_jprb,   0.9468080_jprb,   0.9469628_jprb,  &
        0.9469683_jprb,   0.9465806_jprb,   0.9468755_jprb,   0.9466828_jprb,   0.9471480_jprb,  &
        0.9470276_jprb,   0.9470209_jprb,   0.9468378_jprb,   0.9464890_jprb,   0.9462101_jprb,  &
        0.9459322_jprb,   0.9449111_jprb,   0.9435923_jprb,   0.9416961_jprb,   0.9401403_jprb,  &
        0.9387150_jprb,   0.9374595_jprb,   0.9347988_jprb,   0.9319339_jprb,   0.9295776_jprb,  &
        0.9268476_jprb,   0.9243815_jprb,   0.9224647_jprb,   0.9208075_jprb,   0.9195780_jprb,  &
        0.9183103_jprb,   0.9171674_jprb,   0.9164810_jprb,   0.9160877_jprb,   0.9151877_jprb,  &
        0.9148492_jprb,   0.9142842_jprb,   0.9142084_jprb,   0.9138089_jprb,   0.9137760_jprb,  &
        0.9137531_jprb,   0.9141592_jprb,   0.9136598_jprb,   0.9125727_jprb,   0.9108481_jprb,  &
        0.9093652_jprb,   0.9080561_jprb,   0.9062355_jprb,   0.9046820_jprb,   0.9028210_jprb,  &
        0.9018152_jprb,   0.9008504_jprb,   0.9000632_jprb,   0.8995758_jprb,   0.8989593_jprb,  &
        0.8987811_jprb,   0.8992507_jprb,   0.8999549_jprb,   0.9013391_jprb,   0.9020863_jprb,  &
        0.9025120_jprb,   0.9023982_jprb,   0.9015658_jprb,   0.9008633_jprb,   0.8996401_jprb,  &
        0.8981582_jprb,   0.8969440_jprb,   0.8946483_jprb,   0.8925536_jprb,   0.8906261_jprb,  &
        0.8889833_jprb,   0.8870751_jprb,   0.8845615_jprb,   0.8825631_jprb,   0.8811586_jprb,  &
        0.8796447_jprb,   0.8779839_jprb,   0.8765292_jprb,   0.8754975_jprb,   0.8739760_jprb,  &
        0.8725729_jprb,   0.8714029_jprb,   0.8706908_jprb,   0.8710466_jprb,   0.8699325_jprb,  &
        0.8697992_jprb,   0.8718969_jprb,   0.8713725_jprb,   0.8701416_jprb,   0.8695096_jprb,  &
        0.8698574_jprb,   0.8700698_jprb,   0.8694080_jprb,   0.8693934_jprb,   0.8693246_jprb,  &
        0.8698239_jprb,   0.8696592_jprb,   0.8681608_jprb,   0.8656288_jprb,   0.8654716_jprb,  &
        0.8640761_jprb,   0.8639477_jprb,   0.8635154_jprb,   0.8630069_jprb,   0.8623275_jprb,  &
        0.8623751_jprb,   0.8627441_jprb,   0.8630516_jprb,   0.8638958_jprb,   0.8644919_jprb,  &
        0.8655882_jprb,   0.8666160_jprb,   0.8676174_jprb,   0.8692035_jprb,   0.8695340_jprb,  &
        0.8703975_jprb,   0.8714244_jprb,   0.8715467_jprb,   0.8713564_jprb,   0.8712272_jprb,  &
        0.8714187_jprb,   0.8701625_jprb,   0.8697796_jprb,   0.8688766_jprb,   0.8682391_jprb,  &
        0.8680181_jprb,   0.8676605_jprb,   0.8672657_jprb,   0.8679592_jprb,   0.8675538_jprb,  &
        0.8686572_jprb,   0.8682060_jprb,   0.8688578_jprb,   0.8693632_jprb,   0.8689557_jprb,  &
        0.8681611_jprb,   0.8684876_jprb,   0.8680010_jprb,   0.8675498_jprb,   0.8675414_jprb,  &
        0.8677824_jprb,   0.8665875_jprb,   0.8668503_jprb,   0.8665696_jprb,   0.8671130_jprb,  &
        0.8669835_jprb,   0.8671956_jprb,   0.8683699_jprb,   0.8685648_jprb,   0.8682314_jprb,  &
        0.8683055_jprb,   0.8694246_jprb,   0.8689486_jprb,   0.8693868_jprb,   0.8694460_jprb,  &
        0.8701811_jprb,   0.8704424_jprb,   0.8709887_jprb,   0.8712862_jprb,   0.8721344_jprb,  &
        0.8724745_jprb,   0.8727338_jprb,   0.8740577_jprb,   0.8748575_jprb,   0.8747587_jprb,  &
        0.8762293_jprb,   0.8772818_jprb,   0.8779803_jprb,   0.8791369_jprb,   0.8807610_jprb,  &
        0.8813813_jprb/

  DATA sice_stdv / 0.015_jprb /
  DATA sice_em / &
      0.9370_jprb, 0.9370_jprb, 0.9370_jprb, 0.9370_jprb, 0.9367_jprb, 0.9367_jprb, 0.9366_jprb, 0.9365_jprb, &
      0.9365_jprb, 0.9365_jprb, 0.9365_jprb, 0.9366_jprb, 0.9367_jprb, 0.9370_jprb, 0.9374_jprb, 0.9381_jprb, &
      0.9386_jprb, 0.9393_jprb, 0.9401_jprb, 0.9408_jprb, 0.9415_jprb, 0.9427_jprb, 0.9440_jprb, 0.9452_jprb, &
      0.9464_jprb, 0.9481_jprb, 0.9496_jprb, 0.9511_jprb, 0.9525_jprb, 0.9544_jprb, 0.9563_jprb, 0.9582_jprb, &
      0.9602_jprb, 0.9620_jprb, 0.9640_jprb, 0.9658_jprb, 0.9678_jprb, 0.9702_jprb, 0.9725_jprb, 0.9748_jprb, &
      0.9770_jprb, 0.9792_jprb, 0.9814_jprb, 0.9836_jprb, 0.9856_jprb, 0.9872_jprb, 0.9885_jprb, 0.9897_jprb, &
      0.9905_jprb, 0.9911_jprb, 0.9913_jprb, 0.9913_jprb, 0.9912_jprb, 0.9910_jprb, 0.9907_jprb, 0.9904_jprb, &
      0.9901_jprb, 0.9897_jprb, 0.9893_jprb, 0.9889_jprb, 0.9885_jprb, 0.9880_jprb, 0.9876_jprb, 0.9871_jprb, &
      0.9867_jprb, 0.9864_jprb, 0.9861_jprb, 0.9858_jprb, 0.9854_jprb, 0.9852_jprb, 0.9849_jprb, 0.9846_jprb, &
      0.9844_jprb, 0.9842_jprb, 0.9840_jprb, 0.9838_jprb, 0.9836_jprb, 0.9834_jprb, 0.9832_jprb, 0.9831_jprb, &
      0.9829_jprb, 0.9828_jprb, 0.9826_jprb, 0.9824_jprb, 0.9822_jprb, 0.9821_jprb, 0.9820_jprb, 0.9819_jprb, &
      0.9817_jprb, 0.9816_jprb, 0.9814_jprb, 0.9813_jprb, 0.9811_jprb, 0.9810_jprb, 0.9808_jprb, 0.9807_jprb, &
      0.9805_jprb, 0.9804_jprb, 0.9803_jprb, 0.9801_jprb, 0.9799_jprb, 0.9797_jprb, 0.9796_jprb, 0.9794_jprb, &
      0.9792_jprb, 0.9791_jprb, 0.9789_jprb, 0.9787_jprb, 0.9786_jprb, 0.9785_jprb, 0.9784_jprb, 0.9784_jprb, &
      0.9783_jprb, 0.9782_jprb, 0.9782_jprb, 0.9782_jprb, 0.9781_jprb, 0.9781_jprb, 0.9781_jprb, 0.9781_jprb, &
      0.9781_jprb, 0.9781_jprb, 0.9781_jprb, 0.9780_jprb, 0.9780_jprb, 0.9780_jprb, 0.9780_jprb, 0.9780_jprb, &
      0.9780_jprb, 0.9780_jprb, 0.9779_jprb, 0.9779_jprb, 0.9778_jprb, 0.9778_jprb, 0.9777_jprb, 0.9777_jprb, &
      0.9777_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, &
      0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9775_jprb, 0.9775_jprb, 0.9775_jprb, &
      0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, 0.9777_jprb, 0.9777_jprb, 0.9777_jprb, 0.9777_jprb, &
      0.9777_jprb, 0.9777_jprb, 0.9777_jprb, 0.9777_jprb, 0.9777_jprb, 0.9776_jprb, 0.9776_jprb, 0.9776_jprb, &
      0.9775_jprb, 0.9775_jprb, 0.9774_jprb, 0.9773_jprb, 0.9773_jprb, 0.9773_jprb, 0.9773_jprb, 0.9773_jprb, &
      0.9774_jprb, 0.9774_jprb, 0.9775_jprb, 0.9776_jprb, 0.9777_jprb, 0.9778_jprb, 0.9779_jprb, 0.9780_jprb, &
      0.9781_jprb, 0.9782_jprb, 0.9783_jprb, 0.9785_jprb, 0.9786_jprb, 0.9788_jprb, 0.9790_jprb, 0.9792_jprb, &
      0.9793_jprb, 0.9795_jprb, 0.9797_jprb, 0.9799_jprb, 0.9801_jprb, 0.9802_jprb, 0.9803_jprb, 0.9805_jprb, &
      0.9806_jprb, 0.9807_jprb, 0.9808_jprb, 0.9809_jprb, 0.9810_jprb, 0.9811_jprb, 0.9811_jprb, 0.9811_jprb, &
      0.9810_jprb, 0.9810_jprb, 0.9810_jprb, 0.9809_jprb, 0.9808_jprb, 0.9808_jprb, 0.9807_jprb, 0.9807_jprb, &
      0.9806_jprb, 0.9805_jprb, 0.9805_jprb, 0.9804_jprb, 0.9803_jprb, 0.9802_jprb, 0.9802_jprb, 0.9801_jprb, &
      0.9800_jprb, 0.9799_jprb, 0.9798_jprb, 0.9797_jprb, 0.9797_jprb, 0.9795_jprb, 0.9795_jprb, 0.9794_jprb, &
      0.9793_jprb, 0.9792_jprb, 0.9791_jprb, 0.9791_jprb, 0.9789_jprb, 0.9789_jprb, 0.9788_jprb, 0.9787_jprb, &
      0.9786_jprb, 0.9785_jprb, 0.9785_jprb, 0.9783_jprb, 0.9783_jprb, 0.9782_jprb, 0.9781_jprb, 0.9781_jprb, &
      0.9780_jprb, 0.9779_jprb, 0.9779_jprb, 0.9778_jprb, 0.9777_jprb, 0.9777_jprb, 0.9776_jprb, 0.9775_jprb, &
      0.9774_jprb, 0.9774_jprb, 0.9773_jprb, 0.9772_jprb, 0.9771_jprb, 0.9771_jprb, 0.9770_jprb, 0.9769_jprb, &
      0.9769_jprb, 0.9768_jprb, 0.9767_jprb, 0.9766_jprb, 0.9765_jprb, 0.9765_jprb, 0.9764_jprb, 0.9764_jprb, &
      0.9763_jprb, 0.9762_jprb, 0.9762_jprb, 0.9761_jprb, 0.9761_jprb, 0.9760_jprb, 0.9759_jprb, 0.9758_jprb, &
      0.9757_jprb, 0.9757_jprb, 0.9756_jprb, 0.9756_jprb, 0.9755_jprb, 0.9755_jprb, 0.9754_jprb, 0.9754_jprb, &
      0.9754_jprb, 0.9754_jprb, 0.9754_jprb, 0.9753_jprb, 0.9753_jprb, 0.9753_jprb, 0.9753_jprb, 0.9752_jprb, &
      0.9752_jprb, 0.9752_jprb, 0.9753_jprb, 0.9754_jprb, 0.9755_jprb, 0.9756_jprb, 0.9757_jprb, 0.9757_jprb, &
      0.9758_jprb, 0.9758_jprb, 0.9759_jprb, 0.9759_jprb, 0.9760_jprb, 0.9760_jprb, 0.9760_jprb, 0.9760_jprb, &
      0.9761_jprb, 0.9761_jprb, 0.9761_jprb, 0.9761_jprb, 0.9761_jprb, 0.9761_jprb, 0.9761_jprb, 0.9761_jprb, &
      0.9761_jprb, 0.9760_jprb, 0.9760_jprb, 0.9759_jprb, 0.9759_jprb, 0.9758_jprb, 0.9758_jprb, 0.9757_jprb, &
      0.9757_jprb, 0.9757_jprb, 0.9756_jprb, 0.9756_jprb, 0.9755_jprb, 0.9754_jprb, 0.9753_jprb, 0.9753_jprb, &
      0.9752_jprb, 0.9751_jprb, 0.9751_jprb, 0.9750_jprb, 0.9750_jprb, 0.9749_jprb, 0.9749_jprb, 0.9748_jprb, &
      0.9747_jprb, 0.9746_jprb, 0.9746_jprb, 0.9746_jprb, 0.9745_jprb, 0.9744_jprb, 0.9743_jprb, 0.9742_jprb, &
      0.9742_jprb, 0.9741_jprb, 0.9740_jprb, 0.9739_jprb, 0.9739_jprb, 0.9739_jprb, 0.9738_jprb, 0.9737_jprb, &
      0.9736_jprb, 0.9735_jprb, 0.9735_jprb, 0.9734_jprb, 0.9733_jprb, 0.9732_jprb, 0.9731_jprb, 0.9731_jprb, &
      0.9730_jprb, 0.9729_jprb, 0.9728_jprb, 0.9727_jprb, 0.9726_jprb, 0.9725_jprb, 0.9724_jprb, 0.9723_jprb, &
      0.9723_jprb, 0.9722_jprb, 0.9721_jprb, 0.9720_jprb, 0.9719_jprb, 0.9718_jprb, 0.9717_jprb, 0.9716_jprb, &
      0.9715_jprb, 0.9714_jprb, 0.9713_jprb, 0.9712_jprb, 0.9711_jprb, 0.9709_jprb, 0.9708_jprb, 0.9706_jprb, &
      0.9705_jprb, 0.9704_jprb, 0.9703_jprb, 0.9702_jprb, 0.9700_jprb, 0.9699_jprb, 0.9698_jprb, 0.9696_jprb, &
      0.9695_jprb, 0.9693_jprb, 0.9691_jprb, 0.9690_jprb, 0.9688_jprb, 0.9686_jprb, 0.9685_jprb, 0.9683_jprb, &
      0.9682_jprb, 0.9681_jprb, 0.9679_jprb, 0.9677_jprb, 0.9676_jprb, 0.9674_jprb, 0.9671_jprb, 0.9669_jprb/

  DATA snow_stdv / 0.015_jprb /
  DATA snow_em / &
      0.9716_jprb, 0.9716_jprb, 0.9716_jprb, 0.9716_jprb, 0.9713_jprb, 0.9710_jprb, 0.9708_jprb, 0.9706_jprb, &
      0.9705_jprb, 0.9705_jprb, 0.9705_jprb, 0.9703_jprb, 0.9701_jprb, 0.9700_jprb, 0.9699_jprb, 0.9700_jprb, &
      0.9702_jprb, 0.9703_jprb, 0.9705_jprb, 0.9707_jprb, 0.9710_jprb, 0.9714_jprb, 0.9717_jprb, 0.9722_jprb, &
      0.9728_jprb, 0.9734_jprb, 0.9740_jprb, 0.9746_jprb, 0.9753_jprb, 0.9759_jprb, 0.9765_jprb, 0.9771_jprb, &
      0.9778_jprb, 0.9784_jprb, 0.9792_jprb, 0.9798_jprb, 0.9806_jprb, 0.9814_jprb, 0.9824_jprb, 0.9833_jprb, &
      0.9842_jprb, 0.9852_jprb, 0.9863_jprb, 0.9873_jprb, 0.9882_jprb, 0.9891_jprb, 0.9901_jprb, 0.9908_jprb, &
      0.9914_jprb, 0.9920_jprb, 0.9925_jprb, 0.9926_jprb, 0.9928_jprb, 0.9927_jprb, 0.9926_jprb, 0.9926_jprb, &
      0.9923_jprb, 0.9920_jprb, 0.9918_jprb, 0.9916_jprb, 0.9915_jprb, 0.9913_jprb, 0.9911_jprb, 0.9907_jprb, &
      0.9905_jprb, 0.9903_jprb, 0.9902_jprb, 0.9900_jprb, 0.9897_jprb, 0.9896_jprb, 0.9894_jprb, 0.9892_jprb, &
      0.9890_jprb, 0.9889_jprb, 0.9886_jprb, 0.9884_jprb, 0.9883_jprb, 0.9884_jprb, 0.9885_jprb, 0.9885_jprb, &
      0.9884_jprb, 0.9883_jprb, 0.9881_jprb, 0.9880_jprb, 0.9880_jprb, 0.9880_jprb, 0.9880_jprb, 0.9879_jprb, &
      0.9879_jprb, 0.9879_jprb, 0.9879_jprb, 0.9879_jprb, 0.9879_jprb, 0.9879_jprb, 0.9878_jprb, 0.9877_jprb, &
      0.9876_jprb, 0.9876_jprb, 0.9877_jprb, 0.9876_jprb, 0.9875_jprb, 0.9874_jprb, 0.9873_jprb, 0.9873_jprb, &
      0.9873_jprb, 0.9874_jprb, 0.9875_jprb, 0.9875_jprb, 0.9874_jprb, 0.9874_jprb, 0.9874_jprb, 0.9874_jprb, &
      0.9874_jprb, 0.9874_jprb, 0.9874_jprb, 0.9874_jprb, 0.9874_jprb, 0.9874_jprb, 0.9874_jprb, 0.9873_jprb, &
      0.9873_jprb, 0.9873_jprb, 0.9873_jprb, 0.9873_jprb, 0.9873_jprb, 0.9872_jprb, 0.9871_jprb, 0.9872_jprb, &
      0.9871_jprb, 0.9870_jprb, 0.9870_jprb, 0.9870_jprb, 0.9870_jprb, 0.9869_jprb, 0.9868_jprb, 0.9868_jprb, &
      0.9867_jprb, 0.9866_jprb, 0.9866_jprb, 0.9865_jprb, 0.9865_jprb, 0.9865_jprb, 0.9866_jprb, 0.9866_jprb, &
      0.9865_jprb, 0.9865_jprb, 0.9865_jprb, 0.9866_jprb, 0.9866_jprb, 0.9866_jprb, 0.9866_jprb, 0.9867_jprb, &
      0.9868_jprb, 0.9868_jprb, 0.9868_jprb, 0.9867_jprb, 0.9867_jprb, 0.9867_jprb, 0.9866_jprb, 0.9867_jprb, &
      0.9867_jprb, 0.9867_jprb, 0.9867_jprb, 0.9867_jprb, 0.9867_jprb, 0.9868_jprb, 0.9868_jprb, 0.9868_jprb, &
      0.9869_jprb, 0.9869_jprb, 0.9870_jprb, 0.9872_jprb, 0.9873_jprb, 0.9874_jprb, 0.9874_jprb, 0.9874_jprb, &
      0.9875_jprb, 0.9875_jprb, 0.9875_jprb, 0.9875_jprb, 0.9876_jprb, 0.9876_jprb, 0.9876_jprb, 0.9876_jprb, &
      0.9877_jprb, 0.9877_jprb, 0.9877_jprb, 0.9877_jprb, 0.9878_jprb, 0.9879_jprb, 0.9879_jprb, 0.9879_jprb, &
      0.9878_jprb, 0.9878_jprb, 0.9878_jprb, 0.9879_jprb, 0.9879_jprb, 0.9879_jprb, 0.9879_jprb, 0.9878_jprb, &
      0.9877_jprb, 0.9876_jprb, 0.9876_jprb, 0.9877_jprb, 0.9877_jprb, 0.9876_jprb, 0.9876_jprb, 0.9876_jprb, &
      0.9876_jprb, 0.9876_jprb, 0.9876_jprb, 0.9875_jprb, 0.9874_jprb, 0.9873_jprb, 0.9873_jprb, 0.9872_jprb, &
      0.9870_jprb, 0.9869_jprb, 0.9869_jprb, 0.9869_jprb, 0.9869_jprb, 0.9869_jprb, 0.9868_jprb, 0.9867_jprb, &
      0.9867_jprb, 0.9866_jprb, 0.9866_jprb, 0.9865_jprb, 0.9865_jprb, 0.9865_jprb, 0.9864_jprb, 0.9863_jprb, &
      0.9862_jprb, 0.9862_jprb, 0.9862_jprb, 0.9862_jprb, 0.9862_jprb, 0.9861_jprb, 0.9860_jprb, 0.9860_jprb, &
      0.9860_jprb, 0.9859_jprb, 0.9859_jprb, 0.9859_jprb, 0.9858_jprb, 0.9858_jprb, 0.9857_jprb, 0.9857_jprb, &
      0.9857_jprb, 0.9856_jprb, 0.9855_jprb, 0.9855_jprb, 0.9854_jprb, 0.9854_jprb, 0.9853_jprb, 0.9853_jprb, &
      0.9852_jprb, 0.9852_jprb, 0.9852_jprb, 0.9852_jprb, 0.9852_jprb, 0.9852_jprb, 0.9851_jprb, 0.9850_jprb, &
      0.9850_jprb, 0.9850_jprb, 0.9850_jprb, 0.9851_jprb, 0.9850_jprb, 0.9850_jprb, 0.9849_jprb, 0.9849_jprb, &
      0.9850_jprb, 0.9851_jprb, 0.9851_jprb, 0.9850_jprb, 0.9850_jprb, 0.9849_jprb, 0.9849_jprb, 0.9849_jprb, &
      0.9849_jprb, 0.9848_jprb, 0.9848_jprb, 0.9848_jprb, 0.9849_jprb, 0.9849_jprb, 0.9849_jprb, 0.9849_jprb, &
      0.9849_jprb, 0.9849_jprb, 0.9849_jprb, 0.9848_jprb, 0.9848_jprb, 0.9849_jprb, 0.9849_jprb, 0.9850_jprb, &
      0.9850_jprb, 0.9850_jprb, 0.9851_jprb, 0.9851_jprb, 0.9852_jprb, 0.9852_jprb, 0.9853_jprb, 0.9853_jprb, &
      0.9854_jprb, 0.9854_jprb, 0.9854_jprb, 0.9854_jprb, 0.9854_jprb, 0.9854_jprb, 0.9854_jprb, 0.9854_jprb, &
      0.9855_jprb, 0.9856_jprb, 0.9856_jprb, 0.9856_jprb, 0.9856_jprb, 0.9856_jprb, 0.9856_jprb, 0.9856_jprb, &
      0.9856_jprb, 0.9855_jprb, 0.9855_jprb, 0.9855_jprb, 0.9854_jprb, 0.9853_jprb, 0.9853_jprb, 0.9853_jprb, &
      0.9853_jprb, 0.9853_jprb, 0.9853_jprb, 0.9853_jprb, 0.9853_jprb, 0.9853_jprb, 0.9853_jprb, 0.9853_jprb, &
      0.9852_jprb, 0.9851_jprb, 0.9851_jprb, 0.9850_jprb, 0.9849_jprb, 0.9849_jprb, 0.9849_jprb, 0.9848_jprb, &
      0.9848_jprb, 0.9848_jprb, 0.9848_jprb, 0.9848_jprb, 0.9848_jprb, 0.9848_jprb, 0.9847_jprb, 0.9846_jprb, &
      0.9846_jprb, 0.9846_jprb, 0.9847_jprb, 0.9846_jprb, 0.9845_jprb, 0.9844_jprb, 0.9844_jprb, 0.9843_jprb, &
      0.9842_jprb, 0.9842_jprb, 0.9842_jprb, 0.9842_jprb, 0.9841_jprb, 0.9841_jprb, 0.9840_jprb, 0.9839_jprb, &
      0.9838_jprb, 0.9838_jprb, 0.9837_jprb, 0.9837_jprb, 0.9837_jprb, 0.9836_jprb, 0.9836_jprb, 0.9835_jprb, &
      0.9835_jprb, 0.9834_jprb, 0.9833_jprb, 0.9832_jprb, 0.9832_jprb, 0.9832_jprb, 0.9831_jprb, 0.9831_jprb, &
      0.9830_jprb, 0.9829_jprb, 0.9828_jprb, 0.9828_jprb, 0.9827_jprb, 0.9827_jprb, 0.9826_jprb, 0.9826_jprb, &
      0.9825_jprb, 0.9824_jprb, 0.9824_jprb, 0.9823_jprb, 0.9821_jprb, 0.9821_jprb, 0.9820_jprb, 0.9821_jprb, &
      0.9820_jprb, 0.9820_jprb, 0.9818_jprb, 0.9817_jprb, 0.9817_jprb, 0.9816_jprb, 0.9815_jprb, 0.9815_jprb, &
      0.9815_jprb, 0.9814_jprb, 0.9813_jprb, 0.9813_jprb, 0.9812_jprb, 0.9811_jprb, 0.9811_jprb, 0.9810_jprb/

CONTAINS

!------------------------------------------
! Routines for initialising database
!------------------------------------------

  !> Initialise a UWIRemis atlas data structure. By default the atlas data are
  !! general and can be used for any sensor with IR channels.
  !! If the last four optional arguments are specified then the atlas data are
  !! initialised for use with the specified instrument only: this makes calls
  !! to obtain emissivities from the atlas much faster.
  !! @param[in]       path            path to atlas data files
  !! @param[in]       imonth          month of data to read (1-12)
  !! @param[in]       verbose         flag to turn verbose output on/off
  !! @param[in]       std_init        flag to load covariance data
  !! @param[in]       do_ang_corr     flag to include zenith angle correction
  !! @param[in,out]   atlas           UWIRemis atlas data structure to initialise
  !! @param[out]      err             status on exit
  !! @param[in]       instr_wavenum   channel wavenumber list, optional
  !! @param[in]       platform_id     platform ID from associated rtcoef structure, optional
  !! @param[in]       sat_id          satellite ID from associated rtcoef structure, optional
  !! @param[in]       inst_id         instrument ID from associated rtcoef structure, optional
  SUBROUTINE rttov_uwiremis_init(    &
                      path,          &
                      imonth,        &
                      verbose,       &
                      std_init,      &
                      do_ang_corr,   &
                      atlas,         &
                      err,           &
                      instr_wavenum, &
                      platform_id,   &
                      sat_id,        &
                      inst_id)

    ! Description:
    ! initialize the rttov_uwiremis algorithm by (1) reading in the UW BF IR Global
    ! Emissivty data and (2) the eigenvectors of the laboratory spectra, and (2) make some
    ! precalculations for the PCA regression
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.9    03/31/2009   origianl code E. Borbas UW-Madison/CIMSS
    !  1.0    03/31/2009  New F90 code with structures (E Borbas B Ruston)

    USE rttov_const, ONLY : &
      errorstatus_success, &
      errorstatus_fatal

    CHARACTER(LEN=*),          INTENT(IN)           :: path
    INTEGER(KIND=jpim),        INTENT(IN)           :: imonth
    LOGICAL(KIND=jplm),        INTENT(IN)           :: verbose
    LOGICAL(KIND=jplm),        INTENT(IN)           :: std_init
    LOGICAL(KIND=jplm),        INTENT(IN)           :: do_ang_corr
    TYPE(uwiremis_atlas_data), INTENT(INOUT)        :: atlas
    INTEGER(KIND=jpim),        INTENT(OUT)          :: err
    REAL(KIND=jprb),           INTENT(IN), OPTIONAL :: instr_wavenum(:)
    INTEGER(KIND=jpim),        INTENT(IN), OPTIONAL :: platform_id
    INTEGER(KIND=jpim),        INTENT(IN), OPTIONAL :: sat_id
    INTEGER(KIND=jpim),        INTENT(IN), OPTIONAL :: inst_id

    CHARACTER(LEN=300) :: fn
    CHARACTER(LEN=4)   :: cyear
    CHARACTER(LEN=2)   :: cmonth
    LOGICAL(KIND=jplm) :: file_exists
    CHARACTER(LEN=128) :: errmsg

#ifdef _RTTOV_HDF
    LOGICAL(KIND=jplm) :: hdf_was_open, hdf_was_64bit_reals
#endif

    TRY

    atlas%single_inst = PRESENT(instr_wavenum)
    atlas%std_init    = std_init
    atlas%do_ang_corr = do_ang_corr

    IF (atlas%single_inst) THEN
      IF (.NOT. PRESENT(platform_id) .OR. &
          .NOT. PRESENT(sat_id) .OR. &
          .NOT. PRESENT(inst_id)) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0, 'platform_id, sat_id and inst_id required if instr_wavenum supplied')
      ENDIF
      atlas%platform_id = platform_id
      atlas%sat_id      = sat_id
      atlas%inst_id     = inst_id
    ENDIF

    CALL rttov_uwiremis_nullify_pointers(atlas)

    WRITE(cyear, '(i4)') db_ver_year
    WRITE(cmonth, '(i2.2)') imonth

#ifdef _RTTOV_HDF
    hdf_was_open = is_hdf_open
    hdf_was_64bit_reals = is_hdf_64bit_reals
    CALL open_hdf(.TRUE._jplm, err)
    THROW(err.NE.0)
    WRITE(errmsg,'(a)') ', RTTOV was compiled with HDF5 so HDF5 atlas files are required'
#else
#ifdef _RTTOV_NETCDF
    WRITE(errmsg,'(a)') ', RTTOV was compiled with netCDF so netCDF atlas files are required'
#else
    err = errorstatus_fatal
    THROWM(err.NE.0,'RTTOV must be compiled with HDF5 to use the UWIRemis atlas')
#endif
#endif
    !----------------------------------------------------------------------------
    ! reading the 0.1 degree resolution PCA Coefficients of UW BF IR Land Surface Global Emissivity
    !----------------------------------------------------------------------------
    fn = TRIM(path)//'UWirbfemis_COEF_V2.1_0.1deg_'//cyear//cmonth//'_mask'
#ifdef _RTTOV_HDF
    fn = TRIM(fn)//'.H5'
#else
    fn = TRIM(fn)//'.nc'
#endif
    INQUIRE(FILE=fn, EXIST=file_exists)
    IF (file_exists) THEN
      CALL rttov_uwiremis_read_coefs(err, TRIM(fn), atlas)
      THROW(err.NE.0)
    ELSE
      err = errorstatus_fatal
      THROWM(err .NE. errorstatus_success, 'UWiremis PCA coefs file not found: '//TRIM(fn)//TRIM(errmsg))
    ENDIF
    IF (verbose) INFO('Using UWiremis coefs: '//TRIM(fn))

    !----------------------------------------------------------------------------
    ! reading the 0.5 degree resolution UW IR Land Surface Global Emissivty STD DEV
    !----------------------------------------------------------------------------
    IF (atlas%std_init) THEN
      fn = TRIM(path)//'UWiremis_hsremis_covmat_V1.0_deg0.5_month'//cmonth//'_mask'
#ifdef _RTTOV_HDF
      fn = TRIM(fn)//'.H5'
#else
      fn = TRIM(fn)//'.nc'
#endif
      INQUIRE(FILE=fn, EXIST=file_exists)
      IF (file_exists) THEN
        CALL rttov_uwiremis_read_cov(err, TRIM(fn), atlas)
        THROW(err.NE.0)
      ELSE
        err = errorstatus_fatal
        THROWM(err .NE. errorstatus_success, 'UWiremis covariances file not found: '//TRIM(fn)//TRIM(errmsg))
      ENDIF
      IF (verbose) INFO('Using UWiremis covariances: '//TRIM(fn))
    ENDIF

    !-------------------------------------------------------------------
    ! reading the eigienvectors of the 128 selected laboratory spectra
    !-------------------------------------------------------------------
    ALLOCATE(atlas%pcu(numpcs,numwave))
    fn = TRIM(path)//'UWiremis_labeigvects'
#ifdef _RTTOV_HDF
    fn = TRIM(fn)//'.H5'
#else
    fn = TRIM(fn)//'.nc'
#endif
    INQUIRE(FILE=fn, EXIST=file_exists)
    IF (file_exists) THEN
      CALL rttov_uwiremis_read_labevecs(err, TRIM(fn), atlas)
      THROW(err.NE.0)
    ELSE
      err = errorstatus_fatal
      THROWM(err .NE. errorstatus_success, 'UWiremis eigenvector file not found: '//TRIM(fn)//TRIM(errmsg))
    ENDIF
    IF (verbose) INFO('Using UWiremis eigenvector file: '//TRIM(fn))


    IF (atlas%do_ang_corr) THEN
      !-------------------------------------------------------------------
      ! reading the IGBP ecosystem map file
      !-------------------------------------------------------------------
      fn = TRIM(path)//'uwiremis_igbpmap'
#ifdef _RTTOV_HDF
      fn = TRIM(fn)//'.H5'
#else
      fn = TRIM(fn)//'.nc'
#endif
      INQUIRE(FILE=fn, EXIST=file_exists)
      IF (file_exists) THEN
        CALL rttov_uwiremis_read_igbp(err, TRIM(fn), atlas)
        THROW(err.NE.0)
      ELSE
        err = errorstatus_fatal
        THROWM(err .NE. errorstatus_success, 'UWiremis IGBP atlas file not found: '//TRIM(fn)//TRIM(errmsg))
      ENDIF
      IF (verbose) INFO('Using UWiremis IGBP atlas file: '//TRIM(fn))

      !-------------------------------------------------------------------
      ! reading the angular correcction coefficient file
      !-------------------------------------------------------------------
      IF (imonth == 12 .OR. imonth == 1 .OR. imonth == 2) THEN
        fn = TRIM(path)//'uwiremis_angcorr_201301_V1.5'
      ELSE IF (imonth >= 3 .AND. imonth <= 5) THEN
        fn = TRIM(path)//'uwiremis_angcorr_201204_V1.5'
      ELSE IF (imonth >= 6 .AND. imonth <= 8) THEN
        fn = TRIM(path)//'uwiremis_angcorr_201207_V1.5'
      ELSE
        fn = TRIM(path)//'uwiremis_angcorr_201210_V1.5'
      ENDIF

#ifdef _RTTOV_HDF
      fn = TRIM(fn)//'.H5'
#else
      fn = TRIM(fn)//'.nc'
#endif
      INQUIRE(FILE=fn, EXIST=file_exists)
      IF (file_exists) THEN
        CALL rttov_uwiremis_read_angfunc(err, TRIM(fn), atlas)
        THROW(err.NE.0)
      ELSE
        err = errorstatus_fatal
        THROWM(err .NE. errorstatus_success, 'UWiremis angle correction file not found: '//TRIM(fn)//TRIM(errmsg))
      ENDIF
      IF (verbose) INFO('Using UWiremis angle correction file: '//TRIM(fn))
    ENDIF

#ifdef _RTTOV_HDF
    ! If HDF5 lib was open before rttov_read_coefs was called, make sure the real
    ! kind is the same as it was previously. Otherwise close the library.
    IF (hdf_was_open) THEN
      CALL open_hdf(hdf_was_64bit_reals, err)
      THROW(err.NE.0)
    ELSE
      CALL close_hdf(err)
      THROW(err.NE.0)
    ENDIF
#endif

    IF (atlas%single_inst) THEN
      ! If a channel wavenumber list is supplied for a particular instrument
      ! we can retrieve emissivities much faster.

      ! Interpolate hsr data onto the channel wavenumbers
      atlas%ncoefchans = SIZE(instr_wavenum)
      IF (atlas%do_ang_corr) THEN
        ALLOCATE(atlas%pcu_int(numpcs,atlas%ncoefchans,2), atlas%pcm_int(atlas%ncoefchans,2))
        ALLOCATE(atlas%p1d_int(atlas%ncoefchans,atlas%nb_igbp,2), &
                 atlas%p2d_int(atlas%ncoefchans,atlas%nb_igbp,2), &
                 atlas%p3d_int(atlas%ncoefchans,atlas%nb_igbp,2), &
                 atlas%p1n_int(atlas%ncoefchans,atlas%nb_igbp,2), &
                 atlas%p2n_int(atlas%ncoefchans,atlas%nb_igbp,2), &
                 atlas%p3n_int(atlas%ncoefchans,atlas%nb_igbp,2))
      ELSE
        ALLOCATE(atlas%pcu_int(numpcs,atlas%ncoefchans,1), atlas%pcm_int(atlas%ncoefchans,1))
      ENDIF
      ALLOCATE(atlas%sice_em_int(atlas%ncoefchans), atlas%snow_em_int(atlas%ncoefchans))
      IF (atlas%std_init) ALLOCATE(atlas%cov_emis_int(atlas%cv_pack,atlas%ncoefchans))

      CALL rttov_uwiremis_hsr_interp(instr_wavenum, atlas)

      ! HSR data is no longer required
      DEALLOCATE(atlas%pcu)
      IF (atlas%std_init) DEALLOCATE(atlas%cov_emis)
      IF (atlas%do_ang_corr) DEALLOCATE(atlas%p1d, atlas%p2d, atlas%p3d, atlas%p1n, atlas%p2n, atlas%p3n)
    ENDIF

    CATCH
  END SUBROUTINE rttov_uwiremis_init


#ifdef _RTTOV_HDF
  !> Read UWIRemis atlas covariance data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   UWIRemis atlas data structure
  SUBROUTINE rttov_uwiremis_read_cov(err, fn, atlas)

    ! Description:
    ! read the 0.1 degree resolution UW BF IR Global Emissivty data
    ! from the HDF5 file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.1      03/06/2014 original code J. Vidot

    INTEGER(KIND=jpim),        INTENT(OUT)   :: err
    CHARACTER(LEN=*),          INTENT(IN)    :: fn
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpim)          :: indexlut
    INTEGER(KIND=jpim)          :: i,j
    INTEGER(KIND=jpis), POINTER :: emis_cov(:,:)
    INTEGER(KIND=jpit), POINTER :: pack_cov(:,:)

#include "rttov_hdf_load.interface"

    TRY

    CALL RTTOV_HDF_LOAD(err, fn, "/EMIS", SNAME="MASK", PIT2=pack_cov)
    THROWM(ERR.NE.0,"Cannot load Emissivity Covariance flag from "//TRIM(FN))

    atlas%cv_lats = SIZE(pack_cov,1)
    atlas%cv_lons = SIZE(pack_cov,2)
    ALLOCATE(atlas%cov_emis_lut(atlas%cv_lats,atlas%cv_lons))

    ! Generate the look-up table into the covariance data
    atlas%cov_emis_lut(:,:) = -1
    indexlut = 1
    DO i = 1, atlas%cv_lons
      DO j = 1, atlas%cv_lats
        IF (pack_cov(j,i) > 0) THEN
          atlas%cov_emis_lut(j,i) = indexlut
          indexlut = indexlut + 1
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(pack_cov)
    NULLIFY(pack_cov)

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="DIAGCOV", PIS2=emis_cov)
    THROWM(ERR.NE.0,"Cannot get covariance values from "//TRIM(fn))

    IF (SIZE(emis_cov,1) .NE. numwave) THEN
      DEALLOCATE(emis_cov)
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in covariance data: numwave")
    ENDIF
    atlas%cv_pack = SIZE(emis_cov,2)
    ALLOCATE(atlas%cov_emis(atlas%cv_pack,numwave))

    atlas%cov_emis = TRANSPOSE(emis_cov)
    DEALLOCATE(emis_cov)
    NULLIFY(emis_cov)

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="DIAGCOV_SFAC", RM0=atlas%cov_sfac)
    THROWM(ERR.NE.0,"Cannot get covariance scale factor from "//TRIM(fn))

    CATCH
  END SUBROUTINE rttov_uwiremis_read_cov
#else
#ifdef _RTTOV_NETCDF
  !> Read UWIRemis atlas covariance data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   UWIRemis atlas data structure
  SUBROUTINE  rttov_uwiremis_read_cov(err, fn, atlas)

    ! Description:
    ! read the 0.1 degree resolution UW BF IR Global Emissivty data
    ! from the netCDF file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0    03/31/2009   origianl code B. Ruston

    INTEGER(KIND=jpim),        INTENT(OUT)   :: err
    CHARACTER(LEN=*),          INTENT(IN)    :: fn
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpim) :: nvars     ! number of variables
    INTEGER(KIND=jpim) :: ndims     ! number of dimensions
    INTEGER(KIND=jpim) :: errstat   ! error code
    INTEGER(KIND=jpim) :: recdim    ! record dimension
    INTEGER(KIND=jpim) :: nc_dim(4)

    INTEGER(KIND=jpim) :: i, j
    INTEGER(KIND=jpim) :: ncid, ngatts, nrecs, varid

    INTEGER(KIND=jpis), ALLOCATABLE :: emis_cov(:,:)
    INTEGER(KIND=jpit), ALLOCATABLE :: pack_cov(:,:)

    CHARACTER(LEN=1024) :: strbuf ! string buffer for var
    INTEGER(KIND=jpim)  :: indexlut

    TRY

    ! Open netCDF file.
    errstat = nf_open(TRIM(fn), nf_nowrite, ncid)

    ! Get info on the record dimension for this file.
    errstat = nf_inq(ncid, ndims, nvars, ngatts, recdim)

    DO recdim = 1, ndims
      errstat = nf_inq_dim(ncid, recdim, strbuf, nrecs)
      nc_dim(recdim) = nrecs
    ENDDO

    ! Retrieve database of diagonal of covariance of MODIS emissivity values
    atlas%cv_lats = nc_dim(2)
    atlas%cv_lons = nc_dim(3)
    ALLOCATE(pack_cov(atlas%cv_lats,atlas%cv_lons), atlas%cov_emis_lut(atlas%cv_lats,atlas%cv_lons))
    errstat = nf_inq_varid (ncid, 'mask', varid)
    errstat = nf_get_var_int1(ncid, varid, pack_cov)

    ! Generate the look-up table into the covariance data
    atlas%cov_emis_lut(:,:) = -1
    indexlut = 1
    DO i = 1, nc_dim(3)
      DO j = 1, nc_dim(2)
        IF (pack_cov(j,i) > 0) THEN
          atlas%cov_emis_lut(j,i) = indexlut
          indexlut = indexlut + 1
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(pack_cov)

    ! Retrieve database of diagonal of covariance of MODIS emissivity values
    IF (nc_dim(1) .NE. numwave) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in covariance data: numwave")
    ENDIF
    atlas%cv_pack = nc_dim(4)
    ALLOCATE(emis_cov(numwave,atlas%cv_pack), atlas%cov_emis(atlas%cv_pack,numwave))
    errstat = nf_inq_varid (ncid, 'emis_diagCov', varid)
    errstat = nf_get_var_int2(ncid, varid, emis_cov)
    errstat = nf_get_att_real(ncid, varid, 'scale_factor', atlas%cov_sfac)

    atlas%cov_emis = TRANSPOSE(emis_cov)
    DEALLOCATE(emis_cov)

    errstat = nf_close(ncid)

    CATCH
  END SUBROUTINE rttov_uwiremis_read_cov
#else
  !> Read UWIRemis atlas covariance data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   UWIRemis atlas data structure
  SUBROUTINE rttov_uwiremis_read_cov(err, fn, atlas)
    ! Stub
    INTEGER(KIND=jpim),        INTENT(OUT)   :: err
    CHARACTER(LEN=*),          INTENT(IN)    :: fn
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas
    err = 0
  END SUBROUTINE rttov_uwiremis_read_cov
#endif
#endif

#ifdef _RTTOV_HDF
  !> Read UWIRemis atlas data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   UWIRemis atlas data structure
  SUBROUTINE rttov_uwiremis_read_coefs(err, fn, atlas)

    ! Description:
    ! read the 0.1 degree resolution UW BF IR Global Emissivty data
    ! from the HDF5 file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.1      03/06/2014 original code J. Vidot

    INTEGER(KIND=jpim),        INTENT(OUT)   :: err
    CHARACTER(LEN=*),          INTENT(IN)    :: fn
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpim)          :: indexlut
    INTEGER(KIND=jpim)          :: i, j

#include "rttov_hdf_load.interface"
    TRY

    CALL RTTOV_HDF_LOAD(ERR, fn, "/EMIS", SNAME="FLAG", PIS2=atlas%bfemis_flag)
    THROWM(ERR.NE.0,"Cannot load emissivity flag from "//TRIM(fn))

    atlas%nb_lats = SIZE(atlas%bfemis_flag,1)
    atlas%nb_lons = SIZE(atlas%bfemis_flag,2)
    ALLOCATE(atlas%bfemis_lut(atlas%nb_lats,atlas%nb_lons))

    ! Generate the look-up table into the emissivity data
    atlas%bfemis_lut(:,:) = -1_jpim
    indexlut = 1_jpim
    DO i = 1, atlas%nb_lons
      DO j = 1, atlas%nb_lats
        IF (atlas%bfemis_flag(j,i) > 0) THEN
          atlas%bfemis_lut(j,i) = indexlut
          indexlut = indexlut + 1_jpim
        ENDIF
      ENDDO
    ENDDO

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="COEF", PRM2=atlas%pca_coef)
    THROWM(ERR.NE.0,"Cannot get PCA coeffcients from "//TRIM(fn))

    atlas%nb_pack = SIZE(atlas%pca_coef,1)
    IF (SIZE(atlas%pca_coef,2) .NE. numpcs) THEN
      DEALLOCATE(atlas%pca_coef)
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in emissivity data: numpcs")
    ENDIF

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="COEF_SFAC", PRM1=atlas%pca_sfac)
    THROWM(ERR.NE.0,"Cannot get PCA coeffcients scale factor from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="COEF_OFFS", PRM1=atlas%pca_offs)
    THROWM(ERR.NE.0,"Cannot get PCA coeffcients offset from "//TRIM(fn))

    CATCH
  END SUBROUTINE rttov_uwiremis_read_coefs
#else
#ifdef _RTTOV_NETCDF
  !> Read UWIRemis atlas data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   UWIRemis atlas data structure
  SUBROUTINE  rttov_uwiremis_read_coefs(err, fn, atlas)

    ! Description:
    ! read the 0.1 degree resolution UW BF IR Global Emissivty data
    ! from the netCDF file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.1    11/30/2012   modified for coef files by E. Borbas
    !  1.0    03/31/2009   original code B. Ruston

    INTEGER(KIND=jpim),        INTENT(OUT)   :: err
    CHARACTER(LEN=*),          INTENT(IN)    :: fn
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpim)  :: nvars     ! number of variables
    INTEGER(KIND=jpim)  :: ndims     ! number of dimensions
    INTEGER(KIND=jpim)  :: errstat   ! error code
    INTEGER(KIND=jpim)  :: recdim    ! record dimension
    INTEGER(KIND=jpim)  :: nc_dim(4) ! hng_pnt, lats, lons, pack_len

    INTEGER(KIND=jpim)  :: i, j, k
    INTEGER(KIND=jpim)  :: ncid, ngatts, nrecs, varid

    CHARACTER(LEN=1024) :: strbuf ! string buffer for var
    CHARACTER(LEN=6)    :: cfld
    INTEGER(KIND=jpim)  :: indexlut

    TRY

    ! Open netCDF file.
    errstat = nf_open(TRIM(fn), nf_nowrite, ncid)

    ! Get info on the record dimension for this file.
    errstat = nf_inq(ncid, ndims, nvars, ngatts, recdim)

    DO recdim = 1, ndims
      errstat = nf_inq_dim(ncid, recdim, strbuf, nrecs)
      nc_dim(recdim) = nrecs
    ENDDO

    atlas%nb_lats = nc_dim(2)
    atlas%nb_lons = nc_dim(3)
    ALLOCATE(atlas%bfemis_flag(atlas%nb_lats,atlas%nb_lons), atlas%bfemis_lut(atlas%nb_lats,atlas%nb_lons))

    ! Retrieve emissivity database flag value
    errstat = nf_inq_varid(ncid, 'emis_flag', varid)
    errstat = nf_get_var_int2(ncid, varid, atlas%bfemis_flag)

    ! Generate the look-up table into the emissivity data
    atlas%bfemis_lut(:,:) = -1_jpim
    indexlut = 1_jpim
    DO i = 1, atlas%nb_lons
      DO j = 1, atlas%nb_lats
        IF (atlas%bfemis_flag(j,i) > 0) THEN
          atlas%bfemis_lut(j,i) = indexlut
          indexlut = indexlut + 1_jpim
        ENDIF
      ENDDO
    ENDDO

    ! Retrieve database of 6 coefs
    atlas%nb_pack = nc_dim(4)
    IF (nc_dim(1) .NE. numpcs) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in emissivity data: numpcs")
    ENDIF
    ALLOCATE(atlas%pca_coef(atlas%nb_pack,numpcs),atlas%pca_sfac(numpcs),atlas%pca_offs(numpcs))
    DO k = 1, numpcs
      IF (k < 10) THEN
        WRITE(cfld, '("coef",i1," ")') k
      ELSE
        WRITE(cfld, '("coef",i2.2)') k
      ENDIF
      errstat = nf_inq_varid (ncid, cfld, varid)
      errstat = nf_get_var_real(ncid, varid, atlas%pca_coef(:,k))
      errstat = nf_get_att_real(ncid, varid, 'scale_factor', atlas%pca_sfac(k))
      errstat = nf_get_att_real(ncid, varid, 'add_offset', atlas%pca_offs(k))
    ENDDO

    errstat = nf_close(ncid)

    CATCH
  END SUBROUTINE rttov_uwiremis_read_coefs
#else
  !> Read UWIRemis atlas data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   UWIRemis atlas data structure
  SUBROUTINE rttov_uwiremis_read_coefs(err, fn, atlas)
    ! Stub
    INTEGER(KIND=jpim),        INTENT(OUT)   :: err
    CHARACTER(LEN=*),          INTENT(IN)    :: fn
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas
    err = 0
  END SUBROUTINE rttov_uwiremis_read_coefs
#endif
#endif

#ifdef _RTTOV_HDF
  !> Read UWIRemis atlas eigenvector data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   UWIRemis atlas data structure
  SUBROUTINE rttov_uwiremis_read_labevecs(err, fn, atlas)

    ! Description:
    ! read the eigenvectors of selected laboratory measurements
    ! (created by E borbas) from the HDF5 file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.1      03/06/2014 original code J. Vidot

    INTEGER(KIND=jpim),        INTENT(OUT)   :: err
    CHARACTER(LEN=*),          INTENT(IN)    :: fn
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas

    REAL(KIND=jprm),  POINTER    :: pcu_x(:,:)

#include "rttov_hdf_load.interface"

    TRY

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="PC_SCORES", PRM2=pcu_x)
    THROWM(ERR.NE.0,"Cannot load PC scores from "//TRIM(fn))

    IF (SIZE(pcu_x,1) < numpcs) THEN
      DEALLOCATE(pcu_x)
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in eigenvector data: numpcs")
    ENDIF
    IF (SIZE(pcu_x,2) .NE. numwave) THEN
      DEALLOCATE(pcu_x)
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in eigenvector data: numwave")
    ENDIF
    atlas%pcu = pcu_x(1:numpcs,1:numwave)
    DEALLOCATE(pcu_x)

    CATCH
  END SUBROUTINE rttov_uwiremis_read_labevecs
#else
#ifdef _RTTOV_NETCDF
  !> Read UWIRemis atlas eigenvector data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   UWIRemis atlas data structure
  SUBROUTINE  rttov_uwiremis_read_labevecs(err, fn, atlas)

    ! Description:
    ! read the eigenvectors of selected laboratory measurements
    ! (created by E borbas) from the netCDF file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0    03/31/2009   origianl code B. Ruston

    INTEGER(KIND=jpim),        INTENT(OUT)   :: err
    CHARACTER(LEN=*),          INTENT(IN)    :: fn
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpim) :: errstat
    INTEGER(KIND=jpim) :: ncid, varid
    INTEGER(KIND=jpim) :: ndims, nvars, ngatts, recdim, nrecs
    INTEGER(KIND=jpim) :: nc_dim(2)

    CHARACTER (len=1024) :: strbuf ! string buffer for var

    TRY

    ! Open netCDF file.
    errstat = nf_open(TRIM(fn), nf_nowrite, ncid)

    ! Get info on the record dimension for this file.
    errstat = nf_inq(ncid, ndims, nvars, ngatts, recdim)

    DO recdim = 1, ndims
      errstat = nf_inq_dim(ncid, recdim, strbuf, nrecs)
      nc_dim(recdim) = nrecs
    ENDDO

    IF (nc_dim(1) < numpcs) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in eigenvector data: numpcs")
    ENDIF
    IF (nc_dim(2) .NE. numwave) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in eigenvector data: numwave")
    ENDIF

    ! Read the laboratory eigenvalues into array pcu
    errstat = nf_inq_varid(ncid, 'PC_scores', varid)

    ! Specify integer kinds for compatibility
    errstat = nf_get_vara_real(ncid, varid, INT((/1,1/), KIND=jpim), INT((/numpcs,numwave/), KIND=jpim), atlas%pcu)

    errstat = nf_close(ncid)

    CATCH
  END SUBROUTINE rttov_uwiremis_read_labevecs
#else
  !> Read UWIRemis atlas eigenvector data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   UWIRemis atlas data structure
  SUBROUTINE  rttov_uwiremis_read_labevecs(err, fn, atlas)
    ! Stub
    INTEGER(KIND=jpim),        INTENT(OUT)   :: err
    CHARACTER(LEN=*),          INTENT(IN)    :: fn
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas
    err = 0
  END SUBROUTINE rttov_uwiremis_read_labevecs
#endif
#endif

#ifdef _RTTOV_HDF
  !> Read UWIRemis atlas angular correction data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   UWIRemis atlas data structure
  SUBROUTINE  rttov_uwiremis_read_angfunc(err, fn, atlas)

    ! Description:
    ! read the three day-night coefficients for angular correcton
    ! (created by E borbas) from the HDF5 file into memory.

    INTEGER(KIND=jpim),        INTENT(OUT)   :: err
    CHARACTER(LEN=*),          INTENT(IN)    :: fn
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas

#include "rttov_hdf_load.interface"

    TRY

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p1d", PRM2=atlas%p1d)
    THROWM(ERR.NE.0,"Cannot load angular correction day 1 coefs data from "//TRIM(fn))

    IF (SIZE(atlas%p1d,1) .NE. numwave) THEN
      DEALLOCATE(atlas%p1d)
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in angular correction data: numwave")
    ENDIF

    atlas%nb_igbp = SIZE(atlas%p1d,2)

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p2d", PRM2=atlas%p2d)
    THROWM(ERR.NE.0,"Cannot load angular correction day 2 coefs data from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p3d", PRM2=atlas%p3d)
    THROWM(ERR.NE.0,"Cannot load angular correction day 3 coefs data from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p1n", PRM2=atlas%p1n)
    THROWM(ERR.NE.0,"Cannot load angular correction night 1 coefs data from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p2n", PRM2=atlas%p2n)
    THROWM(ERR.NE.0,"Cannot load angular correction night 2 coefs data from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p3n", PRM2=atlas%p3n)
    THROWM(ERR.NE.0,"Cannot load angular correction night 3 coefs data from "//TRIM(fn))

    CATCH
  END SUBROUTINE rttov_uwiremis_read_angfunc
#else
#ifdef _RTTOV_NETCDF
  !> Read UWIRemis atlas angular correction data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   UWIRemis atlas data structure
  SUBROUTINE  rttov_uwiremis_read_angfunc(err, fn, atlas)

    ! Description:
    ! read the three day-night coefficients for angular correcton
    ! (created by E borbas) from the netCDF file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0    02/04/2014   Eva Borbas

    INTEGER(KIND=jpim),        INTENT(OUT)   :: err
    CHARACTER(LEN=*),          INTENT(IN)    :: fn
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpim) :: nvars     ! number of variables
    INTEGER(KIND=jpim) :: ndims     ! number of dimensions
    INTEGER(KIND=jpim) :: errstat   ! error code
    INTEGER(KIND=jpim) :: recdim    ! record dimension
    INTEGER(KIND=jpim) :: nc_dim(2)

    INTEGER(KIND=jpim) :: ncid, ngatts, nrecs, varid
    CHARACTER (len=1024) :: strbuf ! string buffer for var

    TRY

   ! Open netCDF file
    errstat = nf_open(TRIM(fn), nf_nowrite, ncid)

    ! Get info on the record dimension for this file.
    errstat = nf_inq(ncid, ndims, nvars, ngatts, recdim)

    DO recdim = 1, ndims
      errstat = nf_inq_dim(ncid, recdim, strbuf, nrecs)
      nc_dim(recdim) = nrecs
    ENDDO

    IF (nc_dim(1) .NE. numwave) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in angular correction data: numwave")
    ENDIF
    atlas%nb_igbp = nc_dim(2)
    ALLOCATE(atlas%p1d(numwave,atlas%nb_igbp))
    ALLOCATE(atlas%p2d(numwave,atlas%nb_igbp))
    ALLOCATE(atlas%p3d(numwave,atlas%nb_igbp))
    ALLOCATE(atlas%p1n(numwave,atlas%nb_igbp))
    ALLOCATE(atlas%p2n(numwave,atlas%nb_igbp))
    ALLOCATE(atlas%p3n(numwave,atlas%nb_igbp))

    ! Retrieve p(3) for day and night
    errstat = nf_inq_varid (ncid, 'p1d', varid)
    errstat = nf_get_var_real(ncid, varid, atlas%p1d)
    errstat = nf_inq_varid (ncid, 'p2d', varid)
    errstat = nf_get_var_real(ncid, varid, atlas%p2d)
    errstat = nf_inq_varid (ncid, 'p3d', varid)
    errstat = nf_get_var_real(ncid, varid, atlas%p3d)
    errstat = nf_inq_varid (ncid, 'p1n', varid)
    errstat = nf_get_var_real(ncid, varid, atlas%p1n)
    errstat = nf_inq_varid (ncid, 'p2n', varid)
    errstat = nf_get_var_real(ncid, varid, atlas%p2n)
    errstat = nf_inq_varid (ncid, 'p3n', varid)
    errstat = nf_get_var_real(ncid, varid, atlas%p3n)

    errstat = nf_close(ncid)

    CATCH
  END SUBROUTINE rttov_uwiremis_read_angfunc
#else
  !> Read UWIRemis atlas angular correction data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   UWIRemis atlas data structure
  SUBROUTINE  rttov_uwiremis_read_angfunc(err, fn, atlas)
    ! Stub
    INTEGER(KIND=jpim),        INTENT(OUT)   :: err
    CHARACTER(LEN=*),          INTENT(IN)    :: fn
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas
    err = 0
  END SUBROUTINE rttov_uwiremis_read_angfunc
#endif
#endif

#ifdef _RTTOV_HDF
  !> Read UWIRemis atlas IGBP data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   UWIRemis atlas data structure
  SUBROUTINE  rttov_uwiremis_read_igbp(err, fn, atlas)

    ! Description:
    ! read in the IGBP ecosystem map interpolated for the 0.05x0.05 degree grid and
    ! modified for IR emissivity application
    ! (created by E borbas) from the HDF5 file into memory.

    INTEGER(KIND=jpim),        INTENT(OUT)   :: err
    CHARACTER(LEN=*),          INTENT(IN)    :: fn
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas

#include "rttov_hdf_load.interface"

    TRY

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="IGBP", PIT2=atlas%igbp)
    THROWM(ERR.NE.0,"Cannot load IGBP data from "//TRIM(fn))

    atlas%igbp_lats = SIZE(atlas%igbp,1)
    atlas%igbp_lons = SIZE(atlas%igbp,2)

    CATCH
  END SUBROUTINE rttov_uwiremis_read_igbp
#else
#ifdef _RTTOV_NETCDF
  !> Read UWIRemis atlas IGBP data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   UWIRemis atlas data structure
  SUBROUTINE  rttov_uwiremis_read_igbp(err, fn, atlas)

    ! Description:
    ! read in the IGBP ecosystem map interpolated for the 0.05x0.05 degree grid and
    ! modified for IR emissivity application
    ! (created by E borbas) from the netCDF file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0    04/04/2014   Eva Borbas

    INTEGER(KIND=jpim),        INTENT(OUT)   :: err
    CHARACTER(LEN=*),          INTENT(IN)    :: fn
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpim) :: nvars     ! number of variables
    INTEGER(KIND=jpim) :: ndims     ! number of dimensions
    INTEGER(KIND=jpim) :: errstat   ! error code
    INTEGER(KIND=jpim) :: recdim    ! record dimension
    INTEGER(KIND=jpim) :: nc_dim(2)

    INTEGER(KIND=jpim) :: ncid, varid, ngatts, nrecs
    CHARACTER(LEN=1024) :: strbuf ! string buffer for var

    TRY

    ! Open netCDF file
    errstat = nf_open(TRIM(fn), nf_nowrite, ncid)

    ! Get info on the record dimension for this file.
    errstat = nf_inq(ncid, ndims, nvars, ngatts, recdim)

    DO recdim = 1, ndims
      errstat = nf_inq_dim(ncid, recdim, strbuf, nrecs)
      nc_dim(recdim) = nrecs
    ENDDO

    atlas%igbp_lats = nc_dim(1)
    atlas%igbp_lons = nc_dim(2)
    ALLOCATE(atlas%igbp(atlas%igbp_lats,atlas%igbp_lons))

    ! Retrieve igbp
    errstat = nf_inq_varid (ncid, 'IGBP', varid)
    errstat = nf_get_var_int1(ncid, varid, atlas%igbp)

    errstat = nf_close(ncid)

    CATCH
  END SUBROUTINE rttov_uwiremis_read_igbp
#else
  !> Read UWIRemis atlas IGBP data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   UWIRemis atlas data structure
  SUBROUTINE  rttov_uwiremis_read_igbp(err, fn, atlas)
    ! Stub
    INTEGER(KIND=jpim),        INTENT(OUT)   :: err
    CHARACTER(LEN=*),          INTENT(IN)    :: fn
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas
    err = 0
  END SUBROUTINE rttov_uwiremis_read_igbp
#endif
#endif

!------------------------------------------
! Routines for returning emissivity values
!------------------------------------------

  !> Main subroutine to return emissivities and, optionally, emissivity error
  !! estimates for the given location, surface type, satellite geometry, and
  !! wavenumbers. Note that the platform_id, sat_id, inst_id and ncoefchans
  !! arguments are only used for safety checks when the atlas data were
  !! initialised for a specific instrument.
  !! @param[out]      err                 status on exit
  !! @param[in]       verbose             flag to turn verbose output on/off
  !! @param[in]       nchs                number of channels for which to return emissivities
  !! @param[in]       lat                 latitude for which to return emissivities
  !! @param[in]       lon                 longitude for which to return emissivities
  !! @param[in]       satzen              viewing zenith angle
  !! @param[in]       solzen              solar zenith angle
  !! @param[in]       surfacetype         surface type as defined in RTTOV profile
  !! @param[in]       snowfrac            snow fraction
  !! @param[in]       instr_wavenum       channel wavenumber list for required emissivities
  !! @param[in]       channels            channel list for required emissivities (for single-inst init)
  !! @param[in]       platform_id         platform ID from associated rtcoef structure
  !! @param[in]       sat_id              satellite ID from associated rtcoef structure
  !! @param[in]       inst_id             instrument ID from associated rtcoef structure
  !! @param[in]       ncoefchans          number of channels in associated rtcoef structure
  !! @param[in]       atlas               UWIRemis atlas data structure
  !! @param[out]      instr_emis          calculated emissivities values
  !! @param[out]      instr_emis_cov      estimated errors in emissivities
  !! @param[out]      instr_emis_flag     quality flags for returned emissivities
  !! @param[out]      instr_coeff         PC eigenvalues for returned emissivities, optional
  !! @param[out]      instr_pcu           PC eigenvectors for returned emissivities, optional
  !! @param[out]      instr_pcm           PC constants for returned emissivities, optional
  SUBROUTINE rttov_uwiremis( &
          err,               &
          verbose,           &
          nchs,              &
          lat,               &
          lon,               &
          satzen,            &
          solzen,            &
          surfacetype,       &
          snowfrac,          &
          instr_wavenum,     &
          channels,          &
          platform_id,       &
          sat_id,            &
          inst_id,           &
          ncoefchans,        &
          atlas,             &
          instr_emis,        &
          instr_emis_cov,    &
          instr_emis_flag,   &
          instr_coeff,       &
          instr_pcu,         &
          instr_pcm)

    ! Description:
    ! To compute IR emissivty for a given location and frequency
    ! from the 0.1 degree resolution UW BF IR Global Emissivty data
    ! (at 10 hinge points) (http://cimss.ssec.wisc.edu/iremis/)
    ! and labratory measurements using principal component analyses
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.9       03/31/2009  Original code E Borbas UW-Madison/CIMSS
    !  1.0       03/31/2009  New F90 code with structures (E Borbas B Ruston)

    USE rttov_const, ONLY : surftype_land, surftype_seaice

    INTEGER(KIND=jpim),        INTENT(OUT) :: err
    LOGICAL(KIND=jplm),        INTENT(IN)  :: verbose
    INTEGER(KIND=jpim),        INTENT(IN)  :: nchs
    INTEGER(KIND=jpim),        INTENT(IN)  :: surfacetype
    REAL(KIND=jprb),           INTENT(IN)  :: lat, lon
    REAL(KIND=jprb),           INTENT(IN)  :: satzen, solzen
    REAL(KIND=jprb),           INTENT(IN)  :: snowfrac
    REAL(KIND=jprb),           INTENT(IN)  :: instr_wavenum(nchs)
    INTEGER(KIND=jpim),        INTENT(IN)  :: channels(nchs)
    INTEGER(KIND=jpim),        INTENT(IN)  :: platform_id
    INTEGER(KIND=jpim),        INTENT(IN)  :: sat_id
    INTEGER(KIND=jpim),        INTENT(IN)  :: inst_id
    INTEGER(KIND=jpim),        INTENT(IN)  :: ncoefchans
    TYPE(uwiremis_atlas_data), INTENT(IN)  :: atlas
    REAL(KIND=jprb),           INTENT(OUT) :: instr_emis(nchs)
    REAL(KIND=jprb),           INTENT(OUT) :: instr_emis_cov(nchs)
    INTEGER(KIND=jpim),        INTENT(OUT) :: instr_emis_flag
    REAL(KIND=jprb), OPTIONAL, INTENT(OUT) :: instr_coeff(numpcs)
    REAL(KIND=jprb), OPTIONAL, INTENT(OUT) :: instr_pcu(numpcs,nchs)
    REAL(KIND=jprb), OPTIONAL, INTENT(OUT) :: instr_pcm(nchs)

    REAL(KIND=jprb)    :: coeff(numpcs)

    REAL(KIND=jprb)    :: angcorr(numwave)
    REAL(KIND=jprb)    :: angcorrchn(nchs,2)
    REAL(KIND=jprb)    :: instr_emis_angcorr(nchs,2)
    REAL(KIND=jprb)    :: instr_pcu_angcorr(numpcs,nchs,2)
    REAL(KIND=jprb)    :: instr_pcm_angcorr(nchs,2)

    REAL(KIND=jprb)    :: hsremis(numwave)
    REAL(KIND=jprb)    :: emis_cov(numwave)

    REAL(KIND=jprb)    :: long
    INTEGER(KIND=jpim) :: gridy, gridx, rnd_x, rnd_y, i, ilat, ilon, j

    INTEGER(KIND=jpim) :: igbp_type
    INTEGER(KIND=jpim) :: grid05y, grid05x

    REAL(KIND=jprb)    :: hsr_ir_emis(nchs)
    REAL(KIND=jprb)    :: hsr_ir_emis_cov(nchs)
    REAL(KIND=jprb)    :: cov_buff(numwave)

    LOGICAL(KIND=jplm) :: get_pc
    CHARACTER(LEN=128) :: msg

    TRY

    instr_emis(:)     = hkod
    instr_emis_cov(:) = hkod

    get_pc = PRESENT(instr_coeff) .AND. PRESENT(instr_pcu) .AND. PRESENT(instr_pcm)

    IF (get_pc) THEN
      instr_coeff = 0._jprb
      instr_pcu   = 0._jprb
      instr_pcm   = 0._jprb
    ENDIF

    IF (atlas%single_inst) THEN
      IF (atlas%platform_id /= platform_id .OR. &
          atlas%inst_id /= inst_id .OR. &
          atlas%ncoefchans /= ncoefchans) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0,'UWIRemis atlas called for different coefs to that with which it was initialised')
      ENDIF
      IF (atlas%sat_id /= sat_id) THEN
        IF (verbose) &
          WARN('WARNING: UWIRemis atlas called for instrument with different sat_id to that with which it was initialised')
      ENDIF
    ENDIF

    IF (surfacetype == surftype_land) THEN

    !------------------------------------------------------------
    ! find the closest grid point from the uwiremis database
    !------------------------------------------------------------

      long = MODULO(lon, 360.0_jprb)
      IF (long >= 180.0) THEN
        long = long - 360.0_jprb
      ENDIF

      ilat = NINT(lat * 1000._jprb, KIND=jpim)
      ilon = NINT(long * 1000._jprb, KIND=jpim)

      gridy = NINT(ABS(bfemis_ygrid1 - ilat) * 1._jprb / bfemis_gridres, KIND=jpim) + 1_jpim
      gridx = NINT(ABS(bfemis_xgrid1 - ilon) * 1._jprb / bfemis_gridres, KIND=jpim) + 1_jpim
      gridy = MAX(MIN(gridy, atlas%nb_lats), 1)
      gridx = MAX(MIN(gridx, atlas%nb_lons), 1)

      IF (atlas%do_ang_corr) THEN
        grid05y = NINT(ABS(igbp_ygrid1 - ilat) * 1._jprb / igbp_gridres, KIND=jpim) + 1_jpim
        grid05x = NINT(ABS(igbp_xgrid1 - ilon) * 1._jprb / igbp_gridres, KIND=jpim) + 1_jpim
        grid05y = MAX(MIN(grid05y, atlas%igbp_lats), 1)
        grid05x = MAX(MIN(grid05x, atlas%igbp_lons), 1)
        igbp_type = atlas%igbp(grid05y,grid05x)
      ENDIF

      IF (atlas%std_init) THEN
        rnd_y = NINT(ABS(cov_emis_ygrid1 - ilat) * 1._jprb / cov_emis_gridres, KIND=jpim) + 1_jpim
        rnd_x = NINT(ABS(cov_emis_xgrid1 - ilon) * 1._jprb / cov_emis_gridres, KIND=jpim) + 1_jpim
        rnd_y = MAX(MIN(rnd_y, atlas%cv_lats), 1)
        rnd_x = MAX(MIN(rnd_x, atlas%cv_lons), 1)
      ENDIF

      instr_emis_flag = atlas%bfemis_flag(gridy,gridx)

    !------------------------------
    ! check if it is a land pixel
    !------------------------------

      IF (instr_emis_flag > 0) THEN

        ! Find the emissivity or coefs and covariances

        IF (atlas%bfemis_lut(gridy,gridx) > 0) THEN
          coeff(:) = REAL(atlas%pca_coef(atlas%bfemis_lut(gridy,gridx),:), KIND=jprb) * atlas%pca_sfac(:) + atlas%pca_offs(:)
          IF (get_pc) instr_coeff(:) = coeff(:)
        ELSE
          RETURN
        ENDIF

        IF (atlas%std_init) THEN
          IF (atlas%single_inst) THEN
            IF (atlas%cov_emis_lut(rnd_y,rnd_x) > 0) THEN
              ! Note atlas%cov_emis_int contains the standard deviations (i.e. sqrt is already taken)
              instr_emis_cov(:) = atlas%cov_emis_int(atlas%cov_emis_lut(rnd_y,rnd_x),channels(:))
            ELSE
              instr_emis_cov(:) = default_std
            ENDIF
          ELSE
            IF (atlas%cov_emis_lut(rnd_y,rnd_x) > 0) THEN
              emis_cov(:) = SQRT(REAL(atlas%cov_emis(atlas%cov_emis_lut(rnd_y,rnd_x),:), KIND=jprb) * atlas%cov_sfac)
            ELSE
              emis_cov(:) = default_std
            ENDIF
          ENDIF
        ENDIF

        IF (atlas%single_inst) THEN
          !--------------------------------------------------------------------------------------------------------
          ! compute the emissivity from the PCs
          !-------------------------------------------------------------------------------------------------------

          IF (atlas%do_ang_corr) THEN

            ! The angular correction is not linear so we must apply it to the HSR
            ! emissivities before they are linearly interpolated to the channel wavenumbers

            ! Reconstruct the two nearest HSR emissivities
            DO j = 1, 2
              CALL rttov_uwiremis_recon_emis( &
                  atlas,                      &
                  coeff,                      &
                  channels,                   &
                  instr_emis_angcorr(:,j), j)
            ENDDO

            IF (get_pc) THEN
              IF (coeff(1) /= -999._jprb) THEN
                instr_pcu_angcorr = atlas%pcu_int(:,channels(:),:)
                instr_pcm_angcorr = atlas%pcm_int(channels(:),:)
              ELSE
                instr_pcu_angcorr = 0._jprb
                instr_pcm_angcorr = 0._jprb
              ENDIF
            ENDIF

            IF (satzen > angcorrminzen) THEN

              ! Calculate the corresponding angular correction factors
              DO j = 1, 2
                CALL rttov_uwiremis_angcorr( &
                      atlas%p1d_int(:,:,j), atlas%p2d_int(:,:,j), atlas%p3d_int(:,:,j), &
                      atlas%p1n_int(:,:,j), atlas%p2n_int(:,:,j), atlas%p3n_int(:,:,j), &
                      solzen,                                         &
                      satzen,                                         &
                      igbp_type,                                      &
                      angcorrchn(:,j))
              ENDDO

              ! Apply angular correction
              instr_emis_angcorr = instr_emis_angcorr * angcorrchn
              IF (get_pc) THEN
                DO i = 1, numpcs
                  instr_pcu_angcorr(i,:,:) = instr_pcu_angcorr(i,:,:) * angcorrchn
                ENDDO
                instr_pcm_angcorr = instr_pcm_angcorr * angcorrchn
              ENDIF

            ENDIF

            ! Interpolate to channel wavenumbers (note the reconstructed emissivities
            ! and PCs already include the relevant interpolation weights)
            instr_emis = instr_emis_angcorr(:,1) + instr_emis_angcorr(:,2)
            IF (get_pc) THEN
              instr_pcu = instr_pcu_angcorr(:,:,1) + instr_pcu_angcorr(:,:,2)
              instr_pcm = instr_pcm_angcorr(:,1) + instr_pcm_angcorr(:,2)
            ENDIF

          ELSE

            ! With no angular correction the emissivity PCs have been interpolated
            ! to the channel central wavenumbers so we can just reconstruct the
            ! emissivities directly

            CALL rttov_uwiremis_recon_emis( &
                atlas,                      &
                coeff,                      &! in
                channels,                   &! in
                instr_emis)                  ! out

            IF (get_pc) THEN
              IF (coeff(1) /= -999._jprb) THEN
                instr_pcu = atlas%pcu_int(:,channels(:),1)
                instr_pcm = atlas%pcm_int(channels(:),1)
              ELSE
                instr_pcu = 0._jprb
                instr_pcm = 0._jprb
              ENDIF
            ENDIF

          ENDIF

          !---------------------------------------------------
          ! Linearly blend avg snow emis with snow cover frac
          !---------------------------------------------------
          IF (snowfrac > 0.0_jprb) THEN

            IF (snowfrac > 1.0_jprb) THEN
              instr_emis(:) = atlas%snow_em_int(channels(:))
              IF (atlas%std_init) instr_emis_cov(:) = snow_stdv
            ELSE
              instr_emis(:) = snowfrac * atlas%snow_em_int(channels(:)) + (1.0_jprb - snowfrac) * instr_emis(:)
              IF (atlas%std_init) instr_emis_cov(:) = &
                    (/ (snowfrac * snow_stdv, i = 1, nchs) /) + (1.0_jprb - snowfrac) * instr_emis_cov(:)
            ENDIF

            IF (get_pc) THEN
              ! Currently the PCs are not returned for snowy cases - overwrite values from above
              instr_coeff = 0._jprb
              instr_pcu = 0._jprb
              instr_pcm = 0._jprb
            ENDIF

          ENDIF  ! snow chk

        ELSE

          !--------------------------------------------------------------------------------------------------------
          ! compute the hsr emissivity spectra at 416 wavenumber points from the 10 BF emissivity hinge points
          !-------------------------------------------------------------------------------------------------------

          CALL rttov_uwiremis_recon_hsremis( &
              atlas,                         &! in
              coeff,                         &! in
              hsremis)                        ! out

          !--------------------------------------------------------------------------------------------------------
          ! apply angular correction to the hsr emissivity spectra at 416 wavenumber points
          !-------------------------------------------------------------------------------------------------------

          IF (atlas%do_ang_corr .AND. satzen > angcorrminzen) THEN
            CALL rttov_uwiremis_angcorr( &
                  atlas%p1d, atlas%p2d, atlas%p3d,     &
                  atlas%p1n, atlas%p2n, atlas%p3n,     &
                  solzen,            &! in
                  satzen,            &! in
                  igbp_type,         &! in
                  angcorr)            ! out
            hsremis = hsremis * angcorr
          ELSEIF (get_pc) THEN
            angcorr = 1._jprb
          ENDIF

          !--------------------------------------------------------------------------------
          ! create instrument specific emis/stdv by finidng the closest wavenumber value
          !--------------------------------------------------------------------------------

          IF (get_pc) THEN
            CALL rttov_uwiremis_select_wavenum( &
                atlas,                          &! in
                hsremis,                        &! in
                emis_cov,                       &! in
                nchs,                           &! in
                instr_wavenum(1:nchs),          &! in
                instr_emis,                     &! out
                instr_emis_cov,                 &! out
                instr_pcu,                      &! out
                instr_pcm,                      &! out
                angcorr)                         ! in
          ELSE
            CALL rttov_uwiremis_select_wavenum( &
                atlas,                          &! in
                hsremis,                        &! in
                emis_cov,                       &! in
                nchs,                           &! in
                instr_wavenum(1:nchs),          &! in
                instr_emis,                     &! out
                instr_emis_cov)                  ! out
          ENDIF

          !---------------------------------------------------
          ! Linearly blend avg snow emis with snow cover frac
          !---------------------------------------------------
          IF (snowfrac > 0.0_jprb) THEN

            ! snow_stdv is a fixed stdv
            cov_buff(:) = snow_stdv
            CALL rttov_uwiremis_select_wavenum( &
                      atlas,                    &! in
                      snow_em,                  &! in
                      cov_buff,                 &! in
                      nchs,                     &! in
                      instr_wavenum,            &! in
                      hsr_ir_emis,              &! out
                      hsr_ir_emis_cov)           ! out

            IF (snowfrac > 1.0_jprb) THEN
              instr_emis(:) = hsr_ir_emis(:)
              ! Note: a stdv was passed into hsr_ir_emis_cov -- no sqrt
              IF (atlas%std_init) instr_emis_cov(:) = hsr_ir_emis_cov(:)
            ELSE
              instr_emis(:) = snowfrac * hsr_ir_emis(:) + (1._jprb - snowfrac) * instr_emis(:)
              ! Note: a stdv was passed into hsr_ir_emis_cov -- no sqrt
              IF (atlas%std_init) instr_emis_cov(:) = &
                    snowfrac * hsr_ir_emis_cov(:) + (1._jprb - snowfrac) * instr_emis_cov(:)
            ENDIF

            IF (get_pc) THEN
              ! Currently the PCs are not returned for snowy cases - overwrite values from above
              instr_coeff = 0._jprb
              instr_pcu = 0._jprb
              instr_pcm = 0._jprb
            ENDIF

          ENDIF  ! snow chk

        ENDIF ! single_inst

      ENDIF  ! emis flag > 0

    ELSE IF (surfacetype == surftype_seaice) THEN

      !---------------------------------------
      ! Return emissivity and stdv for seaice
      !---------------------------------------

      IF (atlas%single_inst) THEN

        instr_emis(:) = atlas%sice_em_int(channels(:))
        instr_emis_cov(:) = sice_stdv
        instr_emis_flag = seaice_flag

      ELSE

        ! sice_stdv is a fixed stdv
        cov_buff(:) = sice_stdv
        CALL rttov_uwiremis_select_wavenum( &
                    atlas,                  &! in
                    sice_em,                &! in
                    cov_buff,               &! in
                    nchs,                   &! in
                    instr_wavenum,          &! out
                    instr_emis,             &! out
                    instr_emis_cov)          ! out

        instr_emis_flag = seaice_flag

      ENDIF

    ELSE
      IF (verbose) THEN
        WRITE (msg, '(a)') 'Warning: IR emissivity atlas should only be called for land and seaice surface types'
        CALL rttov_errorreport(errorstatus_success, msg)
      ENDIF
    ENDIF

    ! Cap the final emissivities here for consistency between single-
    ! and multi-instrument initialisation.
    DO i = 1, nchs
      instr_emis(i) = MIN(instr_emis(i), 1._jprb)
    END DO

    CATCH
  END SUBROUTINE rttov_uwiremis

  !> Interpolate PCs onto instrument channel wavenumbers. This is called
  !! during initialisation when single-instrument init is selected.
  !! @param[in]      instr_wavenum   instrument channel wavenumbers
  !! @param[in,out]  atlas           UWIRemis atlas data structure
  SUBROUTINE rttov_uwiremis_hsr_interp(instr_wavenum, atlas)

    ! Description:
    ! Initialisation for a single instrument.
    ! Interpolate PC data onto a specific set of wavenumbers:
    ! this can be precomputed for a given instrument during
    ! atlas initialisation to enable very rapid emissivity
    ! calculations.

    REAL(KIND=jprb),           INTENT(IN) :: instr_wavenum(:)
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpim) :: j, k, nchs

    REAL(KIND=jprb) :: dist(numwave)
    REAL(KIND=jprb) :: mindist
    INTEGER(KIND=jpim) :: ind_mindist

    REAL(KIND=jprb) :: dwvnum1, dwvnum2
    REAL(KIND=jprb) :: pcu1(numpcs), pcu2(numpcs), pcm1, pcm2
    REAL(KIND=jprb) :: sice_em1, sice_em2, snow_em1, snow_em2
    REAL(KIND=jprb) :: cov_emis1(atlas%cv_pack), cov_emis2(atlas%cv_pack)

    !---------------------------------------------------------------
    ! finding the closest frequency from the hsr emissivity spectr
    !--------------------------------------------------------------

    nchs = SIZE(instr_wavenum)

    DO j = 1, nchs

      ! The angular correction is not linear so in this case we need to
      ! store data for the nearest two HSR wavenumbers so that the
      ! corresponding angcorr values can be applied to the reconstructed
      ! HSR emissivities and then the linearinterpolation to channel
      ! wavenumber can be done.

      IF (instr_wavenum(j) <= hsr_wavenum(1)) THEN

        atlas%pcu_int(:,j,1)   = REAL(atlas%pcu(:,1), KIND=jprb)
        atlas%pcm_int(j,1)     = pcm(1)
        IF (atlas%do_ang_corr) THEN
          atlas%pcu_int(:,j,2) = 0._jprb
          atlas%pcm_int(j,2)   = 0._jprb
          atlas%p1d_int(j,:,1) = atlas%p1d(1,:)
          atlas%p2d_int(j,:,1) = atlas%p2d(1,:)
          atlas%p3d_int(j,:,1) = atlas%p3d(1,:)
          atlas%p1n_int(j,:,1) = atlas%p1n(1,:)
          atlas%p2n_int(j,:,1) = atlas%p2n(1,:)
          atlas%p3n_int(j,:,1) = atlas%p3n(1,:)
          atlas%p1d_int(j,:,2) = 0._jprb
          atlas%p2d_int(j,:,2) = 0._jprb
          atlas%p3d_int(j,:,2) = 0._jprb
          atlas%p1n_int(j,:,2) = 0._jprb
          atlas%p2n_int(j,:,2) = 0._jprb
          atlas%p3n_int(j,:,2) = 0._jprb
        ENDIF
        atlas%sice_em_int(j) = sice_em(1)
        atlas%snow_em_int(j) = snow_em(1)
        IF (atlas%std_init) atlas%cov_emis_int(:,j) = SQRT(atlas%cov_emis(:,1) * atlas%cov_sfac)

      ELSE IF (instr_wavenum(j) >= hsr_wavenum(numwave)) THEN

        atlas%pcu_int(:,j,1)   = REAL(atlas%pcu(:,numwave), KIND=jprb)
        atlas%pcm_int(j,1)     = pcm(numwave)
        IF (atlas%do_ang_corr) THEN
          atlas%pcu_int(:,j,2) = 0._jprb
          atlas%pcm_int(j,2)   = 0._jprb
          atlas%p1d_int(j,:,1) = atlas%p1d(numwave,:)
          atlas%p2d_int(j,:,1) = atlas%p2d(numwave,:)
          atlas%p3d_int(j,:,1) = atlas%p3d(numwave,:)
          atlas%p1n_int(j,:,1) = atlas%p1n(numwave,:)
          atlas%p2n_int(j,:,1) = atlas%p2n(numwave,:)
          atlas%p3n_int(j,:,1) = atlas%p3n(numwave,:)
          atlas%p1d_int(j,:,2) = 0._jprb
          atlas%p2d_int(j,:,2) = 0._jprb
          atlas%p3d_int(j,:,2) = 0._jprb
          atlas%p1n_int(j,:,2) = 0._jprb
          atlas%p2n_int(j,:,2) = 0._jprb
          atlas%p3n_int(j,:,2) = 0._jprb
        ENDIF
        atlas%sice_em_int(j) = sice_em(numwave)
        atlas%snow_em_int(j) = snow_em(numwave)
        IF (atlas%std_init) atlas%cov_emis_int(:,j) = SQRT(atlas%cov_emis(:,numwave) * atlas%cov_sfac)

      ELSE ! within wavenumber compute range

        mindist = 100._jprb
        ind_mindist = 100000_jpim

        DO k = 1, numwave
          ! calculate distances between the instr freq and hsr emissivities
          dist(k) = ABS(instr_wavenum(j) - hsr_wavenum(k))

          ! finding the closest frequency from the hsr emissivity
          IF (dist(k) <= mindist) THEN
            mindist = dist(k)
            ind_mindist = k
          ENDIF
        ENDDO

  !--------------------------------
  ! Interpolate values
  !--------------------------------

  ! Bilinear mean of the two closest spectral points
        k = 1
        IF (atlas%std_init) THEN
          cov_emis1(:) = 0._jprb
          cov_emis2(:) = 0._jprb
        ENDIF

        IF (instr_wavenum(j) <= hsr_wavenum(ind_mindist)) k = -1

        dwvnum1 = dist(ind_mindist) / (dist(ind_mindist) + dist(ind_mindist + k))
        dwvnum2 = 1._jprb - dwvnum1

        pcu1(:)  = dwvnum1 * REAL(atlas%pcu(:,ind_mindist + k), KIND=jprb)
        pcu2(:)  = dwvnum2 * REAL(atlas%pcu(:,ind_mindist), KIND=jprb)
        pcm1     = dwvnum1 * pcm(ind_mindist + k)
        pcm2     = dwvnum2 * pcm(ind_mindist)
        sice_em1 = dwvnum1 * sice_em(ind_mindist + k)
        sice_em2 = dwvnum2 * sice_em(ind_mindist)
        snow_em1 = dwvnum1 * snow_em(ind_mindist + k)
        snow_em2 = dwvnum2 * snow_em(ind_mindist)

        IF (atlas%do_ang_corr) THEN
          atlas%pcu_int(:,j,1) = pcu1(:)
          atlas%pcm_int(j,1)   = pcm1
          atlas%pcu_int(:,j,2) = pcu2(:)
          atlas%pcm_int(j,2)   = pcm2
          atlas%p1d_int(j,:,1) = atlas%p1d(ind_mindist + k,:)
          atlas%p2d_int(j,:,1) = atlas%p2d(ind_mindist + k,:)
          atlas%p3d_int(j,:,1) = atlas%p3d(ind_mindist + k,:)
          atlas%p1n_int(j,:,1) = atlas%p1n(ind_mindist + k,:)
          atlas%p2n_int(j,:,1) = atlas%p2n(ind_mindist + k,:)
          atlas%p3n_int(j,:,1) = atlas%p3n(ind_mindist + k,:)
          atlas%p1d_int(j,:,2) = atlas%p1d(ind_mindist,:)
          atlas%p2d_int(j,:,2) = atlas%p2d(ind_mindist,:)
          atlas%p3d_int(j,:,2) = atlas%p3d(ind_mindist,:)
          atlas%p1n_int(j,:,2) = atlas%p1n(ind_mindist,:)
          atlas%p2n_int(j,:,2) = atlas%p2n(ind_mindist,:)
          atlas%p3n_int(j,:,2) = atlas%p3n(ind_mindist,:)
        ELSE
          atlas%pcu_int(:,j,1) = pcu1(:) + pcu2(:)
          atlas%pcm_int(j,1)   = pcm1 + pcm2
        ENDIF
        atlas%sice_em_int(j) = sice_em1 + sice_em2
        atlas%snow_em_int(j) = snow_em1 + snow_em2

        IF (atlas%std_init) THEN
          cov_emis1(:) = dwvnum1 * SQRT(REAL(atlas%cov_emis(:,ind_mindist + k), KIND=jprb) * atlas%cov_sfac)
          cov_emis2(:) = dwvnum2 * SQRT(REAL(atlas%cov_emis(:,ind_mindist), KIND=jprb) * atlas%cov_sfac)
          atlas%cov_emis_int(:,j) = cov_emis1(:) + cov_emis2(:)
        ENDIF

      ENDIF

    ENDDO
  END SUBROUTINE rttov_uwiremis_hsr_interp

  !> Reconstruct emissivities at instrument wavenumbers from interpolated PCs
  !! (used when atlas is initialised for a specific instrument)
  !! @param[in]   atlas           UWIRemis atlas data structure
  !! @param[in]   coef            PC coefficients
  !! @param[in]   channels        channels for which emissivities are required
  !! @param[out]  emis            reconstructed emissivities
  !! @param[in]   ind             selects which of the two interpolated sets of PCs to use
  !!                                (used with angular correction), optional
  SUBROUTINE rttov_uwiremis_recon_emis(atlas, coef, channels, emis, ind)

    ! Description:
    ! Used with single-instrument initialisation.
    ! Reconstruct the emissivities at the instrument wavenumbers
    ! from interpolated Principal Components.

    TYPE(uwiremis_atlas_data), INTENT(IN)           :: atlas
    REAL(KIND=jprb),           INTENT(IN)           :: coef(numpcs)
    INTEGER(KIND=jpim),        INTENT(IN)           :: channels(:)
    REAL(KIND=jprb),           INTENT(OUT)          :: emis(SIZE(channels))
    INTEGER(KIND=jpim),        INTENT(IN), OPTIONAL :: ind

    INTEGER(KIND=jpim) :: j, k, nchn

    !-----------------------------------
    ! apply regcoef to get the emissivities
    !-----------------------------------

    IF (coef(1) /= -999._jprb) THEN
      j = 1
      IF (PRESENT(ind)) j = ind
      nchn = SIZE(emis)
      DO k = 1, nchn
        emis(k) = SUM(coef(:) * atlas%pcu_int(:,channels(k),j)) + atlas%pcm_int(channels(k),j)
      ENDDO
    ELSE
      emis = hkod
    ENDIF

  END SUBROUTINE rttov_uwiremis_recon_emis

  !> Reconstruct high-spectral-resolution emissivities from PCs
  !! @param[in]   atlas           UWIRemis atlas data structure
  !! @param[in]   coef            PC coefficients
  !! @param[out]  hsremis         reconstructed emissivities
  SUBROUTINE rttov_uwiremis_recon_hsremis(atlas, coef, hsremis)

    ! Description:
    ! Used with multiple-instrument initialisation.
    ! To creates high spectra resolution emissivties at 416 wavenumbers
    ! from the PCA Coefficitents of the UW BF IR Global Emissivty data
    ! and labratory measurements using principal component analyses
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.9       03/31/2009  Original code E Borbas UW-Madison/CIMSS
    !  1.0       03/31/2009  New F90 code with structures (E Borbas B Ruston)
    !  1.1       11/30/2012  Removed the coef calcualtion part (E Borbas )

    TYPE(uwiremis_atlas_data), INTENT(IN)  :: atlas
    REAL(KIND=jprb),           INTENT(IN)  :: coef(numpcs)
    REAL(KIND=jprb),           INTENT(OUT) :: hsremis(numwave)

    INTEGER(KIND=jpim) :: k

    !-----------------------------------
    ! apply regcoef to get the hsr dataset
    !-----------------------------------

    IF (coef(1) /= -999._jprb) THEN
      DO k = 1, numwave
        hsremis(k) = SUM(coef(:) * REAL(atlas%pcu(:,k), KIND=jprb)) + pcm(k)
      ENDDO
    ELSE
      hsremis = hkod
    ENDIF

  END SUBROUTINE rttov_uwiremis_recon_hsremis

  !> Linear interpolation of high-spectral-resolution emissivity data onto
  !! instrument channel wavenumbers
  !! @param[in]      atlas           UWIRemis atlas data structure
  !! @param[in]      hsremis         high-spectral-resolution emissivity data
  !! @param[in]      emis_cov        high-spectral-resolution emissivity covariance data
  !! @param[in]      nchs            number of instrument channels
  !! @param[in]      instr_wavenum   channel central wavenumbers
  !! @param[out]     instr_emis      interpolated emissivity values
  !! @param[out]     instr_emis_cov  interpolated emissivity error values
  !! @param[in,out]  instr_pcu       interpolated PC eigenvectors, optional
  !! @param[in,out]  instr_pcm       interpolated PC constants, optional
  !! @param[in]      angcorr         angular correction factor, optional, required for pcu/pcm output

  SUBROUTINE rttov_uwiremis_select_wavenum ( &
          atlas,                             &
          hsremis,                           &
          emis_cov,                          &
          nchs,                              &
          instr_wavenum,                     &
          instr_emis,                        &
          instr_emis_cov,                    &
          instr_pcu,                         &
          instr_pcm,                         &
          angcorr)

    ! Description:
    ! Used with multiple-instrument initialisation.
    ! Subroutine to find the closest wavenumber from the UW HSR emissivity spectra
    ! for the instrument frequency and assign the instrument emissivity by choosing the
    ! closest spectral point value or bilinear interpolating  between the two
    ! closest spectral point values
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.9       03/31/2009  Original code E Borbas UW-Madison/CIMSS
    !  1.0       03/31/2009  New F90 code with structures (E Borbas B Ruston)

    TYPE(uwiremis_atlas_data), INTENT(IN)    :: atlas
    REAL(KIND=jprb),           INTENT(IN)    :: hsremis(numwave)
    REAL(KIND=jprb),           INTENT(IN)    :: emis_cov(numwave)
    INTEGER(KIND=jpim),        INTENT(IN)    :: nchs
    REAL(KIND=jprb),           INTENT(IN)    :: instr_wavenum(nchs)
    REAL(KIND=jprb),           INTENT(OUT)   :: instr_emis(nchs)
    REAL(KIND=jprb),           INTENT(OUT)   :: instr_emis_cov(nchs)
    REAL(KIND=jprb), OPTIONAL, INTENT(INOUT) :: instr_pcu(numpcs,nchs)
    REAL(KIND=jprb), OPTIONAL, INTENT(INOUT) :: instr_pcm(nchs)
    REAL(KIND=jprb), OPTIONAL, INTENT(IN)    :: angcorr(numwave)

    INTEGER(KIND=jpim) :: j, k

    REAL(KIND=jprb) :: dist(numwave)
    REAL(KIND=jprb) :: mindist
    INTEGER(KIND=jpim) :: ind_mindist

    REAL(KIND=jprb) :: dwvnum1, dwvnum2, dwvsum
    REAL(KIND=jprb) :: hsremis1, hsremis2, emis_cov1, emis_cov2
    REAL(KIND=jprb) :: pcu1(numpcs), pcu2(numpcs), pcm1, pcm2
    LOGICAL(KIND=jplm) :: lcpu_emis, lcpu_cov, lcpu_pc

    ! initialize instrument emissivity

    !---------------------------------------------------------------
    ! finding the closest frequency from the hsr emissivity spectra
    !---------------------------------------------------------------
    lcpu_emis = .NOT. ALL(hsremis == hsremis(1))
    IF (atlas%std_init) THEN
      lcpu_cov = .NOT. ALL(emis_cov == emis_cov(1))
    ELSE
      lcpu_cov = .FALSE.
    ENDIF
    lcpu_pc = PRESENT(instr_pcu) .AND. PRESENT(instr_pcm) .AND. PRESENT(angcorr)

    IF (lcpu_emis .OR. lcpu_cov) THEN
      instr_emis(:)       = hkod
      instr_emis_cov(:)   = hkod
      DO j = 1, nchs

        IF (instr_wavenum(j) <= hsr_wavenum(1)) THEN

          instr_emis(j) = hsremis(1)
          IF (atlas%std_init) instr_emis_cov(j) = emis_cov(1)
          IF (lcpu_pc) THEN
            instr_pcu(:,j) = atlas%pcu(:,1) * angcorr(1)
            instr_pcm(j) = pcm(1) * angcorr(1)
          ENDIF

        ELSEIF (instr_wavenum(j) >= hsr_wavenum(numwave)) THEN

          instr_emis(j) = hsremis(numwave)
          IF (atlas%std_init) instr_emis_cov(j) = emis_cov(numwave)
          IF (lcpu_pc) THEN
            instr_pcu(:,j) = atlas%pcu(:,numwave) * angcorr(numwave)
            instr_pcm(j) = pcm(numwave) * angcorr(numwave)
          ENDIF

        ELSE ! within wavenumber compute range

          mindist = 100._jprb
          ind_mindist = 100000_jpim

          DO k = 1, numwave

            ! calculate distances between the instr freq end hsr emissivities
            dist(k) = ABS(instr_wavenum(j) - hsr_wavenum(k))

            ! finding the closest frequency from the hsr emissivity
            IF (dist(k) <=  mindist) THEN
              mindist = dist(k)
              ind_mindist = k
            ENDIF
          ENDDO

    !--------------------------------
    ! assign instrument emissivity
    !--------------------------------
    !  closest spectral point
    !                       instr_emis(j)=hsremis(ind_mindist)
    !                       instr_emis_cov(j)=atlas%cov_emis(ind_mindist)

    ! or bilinear mean of the two closest spectral points

          k = 1
          IF (instr_wavenum(j) <= hsr_wavenum(ind_mindist)) k = -1

          dwvnum1 = dist(ind_mindist)
          dwvnum2 = dist(ind_mindist + k)
          dwvsum = dwvnum1 + dwvnum2

          IF (lcpu_emis) THEN
            hsremis1 = dwvnum1 * hsremis(ind_mindist + k)
            hsremis2 = dwvnum2 * hsremis(ind_mindist)
            instr_emis(j) = (hsremis1 + hsremis2) / dwvsum
          ELSE
            instr_emis(j) = hsremis(1)
          ENDIF

          IF (atlas%std_init) THEN
            IF (lcpu_cov) THEN
              emis_cov1 = dwvnum1 * emis_cov(ind_mindist + k)
              emis_cov2 = dwvnum2 * emis_cov(ind_mindist)
              instr_emis_cov(j) = (emis_cov1 + emis_cov2) / dwvsum
            ELSE
              instr_emis_cov(j) = emis_cov(1)
            ENDIF
          ENDIF

          IF (lcpu_pc) THEN
            pcu1 = dwvnum1 * atlas%pcu(:,ind_mindist + k) * angcorr(ind_mindist + k)
            pcu2 = dwvnum2 * atlas%pcu(:,ind_mindist) * angcorr(ind_mindist)
            instr_pcu(:,j) = (pcu1 + pcu2) / dwvsum
            pcm1 = dwvnum1 * pcm(ind_mindist + k) * angcorr(ind_mindist + k)
            pcm2 = dwvnum2 * pcm(ind_mindist) * angcorr(ind_mindist)
            instr_pcm(j) = (pcm1 + pcm2) / dwvsum
          ENDIF

        ENDIF    !==  (instr_wavenum(j) <= hsr_wavenum(1))

      ENDDO

    ELSE  ! all logical computes (lcpu_emis, lcpu_cov) are false
      instr_emis(:) = hsremis(1)
      IF (atlas%std_init) instr_emis_cov(:) = emis_cov(1)
    ENDIF

  END SUBROUTINE rttov_uwiremis_select_wavenum

  !> Calculate zenith angle correction
  !! @param[in]       p1d         Daytime angular correction coefficients
  !! @param[in]       p2d         Daytime angular correction coefficients
  !! @param[in]       p3d         Daytime angular correction coefficients
  !! @param[in]       p1n         Nighttime angular correction coefficients
  !! @param[in]       p2n         Nighttime angular correction coefficients
  !! @param[in]       p3n         Nighttime angular correction coefficients
  !! @param[in]       solzen      Solar zenith angle
  !! @param[in]       satzen      Satellite zenith angle
  !! @param[in]       igbp_type   IGBP surface type classification
  !! @param[out]      angcorr     calculated angular correction factor
  SUBROUTINE rttov_uwiremis_angcorr(  &
        p1d, p2d, p3d, p1n, p2n, p3n, &
        solzen,                       &
        satzen,                       &
        igbp_type,                    &
        angcorr)

    ! Description:
    ! Subroutine to calculate the satellite zenith angle correction
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.9       01/22/2014  Original code E Borbas UW-Madison/CIMSS

    REAL(KIND=jprm),    INTENT(IN)  :: p1d(:,:), p2d(:,:), p3d(:,:)
    REAL(KIND=jprm),    INTENT(IN)  :: p1n(:,:), p2n(:,:), p3n(:,:)
    REAL(KIND=jprb),    INTENT(IN)  :: satzen, solzen
    INTEGER(KIND=jpim), INTENT(IN)  :: igbp_type
    REAL(KIND=jprb),    INTENT(OUT) :: angcorr(:)

    INTEGER(KIND=jpim) :: indclst
    REAL(KIND=jprb)    :: pv(SIZE(angcorr)), dpv(SIZE(angcorr)), p1(SIZE(angcorr))
!     CHARACTER(LEN=128) :: msg

    IF (igbp_type > 0) THEN

      indclst = NINT(ABS(satzen))

      ! Note from E Borbas:
      ! We have the nonlinear fitting function (with three coefficients) and then we
      ! calculate the differences of the values between the nadir point (p3) and at
      ! the giving viewing angle. In the future we may develop the method further
      ! and p1 may not be equal to p3, but also added some bias.
      ! Therefore current code (where p3 is not actually used) is kept for clarity.

      IF (solzen <= angcorrterm) THEN
        pv(:) = p1d(:,igbp_type) * indclst**2 + p2d(:,igbp_type) * indclst + p3d(:,igbp_type)
        p1(:) = p3d(:,igbp_type)
      ELSE
        pv(:) = p1n(:,igbp_type) * indclst**2 + p2n(:,igbp_type) * indclst + p3n(:,igbp_type)
        p1(:) = p3n(:,igbp_type)
      ENDIF

      dpv(:) = p1(:) - pv(:)
      angcorr(:) = (1._jprb - dpv(:))

    ELSE
!       IF (verbose) &
!         WRITE (msg,'(a)') &
!           'Warning: IGBP type = ocean water or coast line: no IR emissivity angular correction is applied'
!         CALL rttov_errorreport(errorstatus_success, msg)
      angcorr(:) = 1._jprb
    ENDIF

  END SUBROUTINE rttov_uwiremis_angcorr

!------------------------------------
! Routine to deallocate atlas arrays
!------------------------------------

  !> Deallocate data in UWIRemis atlas data structure
  !! @param[in,out]   atlas   UWIRemis atlas data structure to deallocate
  SUBROUTINE rttov_uwiremis_close_atlas(atlas)
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas

    IF ( ASSOCIATED(atlas%bfemis_flag)  ) DEALLOCATE(atlas%bfemis_flag)
    IF ( ASSOCIATED(atlas%bfemis_lut)   ) DEALLOCATE(atlas%bfemis_lut)
    IF ( ASSOCIATED(atlas%pca_coef)     ) DEALLOCATE(atlas%pca_coef)
    IF ( ASSOCIATED(atlas%pca_offs)     ) DEALLOCATE(atlas%pca_offs)
    IF ( ASSOCIATED(atlas%pca_sfac)     ) DEALLOCATE(atlas%pca_sfac)
    IF ( ASSOCIATED(atlas%pcu)          ) DEALLOCATE(atlas%pcu)
    IF ( ASSOCIATED(atlas%cov_emis_lut) ) DEALLOCATE(atlas%cov_emis_lut)
    IF ( ASSOCIATED(atlas%cov_emis)     ) DEALLOCATE(atlas%cov_emis)
    IF ( ASSOCIATED(atlas%pcu_int)      ) DEALLOCATE(atlas%pcu_int)
    IF ( ASSOCIATED(atlas%pcm_int)      ) DEALLOCATE(atlas%pcm_int)
    IF ( ASSOCIATED(atlas%sice_em_int)  ) DEALLOCATE(atlas%sice_em_int)
    IF ( ASSOCIATED(atlas%snow_em_int)  ) DEALLOCATE(atlas%snow_em_int)
    IF ( ASSOCIATED(atlas%cov_emis_int) ) DEALLOCATE(atlas%cov_emis_int)

    IF ( ASSOCIATED(atlas%igbp)      ) DEALLOCATE(atlas%igbp)
    IF ( ASSOCIATED(atlas%p1d)       ) DEALLOCATE(atlas%p1d)
    IF ( ASSOCIATED(atlas%p2d)       ) DEALLOCATE(atlas%p2d)
    IF ( ASSOCIATED(atlas%p3d)       ) DEALLOCATE(atlas%p3d)
    IF ( ASSOCIATED(atlas%p1n)       ) DEALLOCATE(atlas%p1n)
    IF ( ASSOCIATED(atlas%p2n)       ) DEALLOCATE(atlas%p2n)
    IF ( ASSOCIATED(atlas%p3n)       ) DEALLOCATE(atlas%p3n)
    IF ( ASSOCIATED(atlas%p1d_int)   ) DEALLOCATE(atlas%p1d_int)
    IF ( ASSOCIATED(atlas%p2d_int)   ) DEALLOCATE(atlas%p2d_int)
    IF ( ASSOCIATED(atlas%p3d_int)   ) DEALLOCATE(atlas%p3d_int)
    IF ( ASSOCIATED(atlas%p1n_int)   ) DEALLOCATE(atlas%p1n_int)
    IF ( ASSOCIATED(atlas%p2n_int)   ) DEALLOCATE(atlas%p2n_int)
    IF ( ASSOCIATED(atlas%p3n_int)   ) DEALLOCATE(atlas%p3n_int)

    CALL rttov_uwiremis_nullify_pointers(atlas)

  END SUBROUTINE rttov_uwiremis_close_atlas

  !> Nullify pointers in UWIRemis atlas data structure
  !! @param[in,out]   atlas   UWIRemis atlas data structure to nullify
  SUBROUTINE rttov_uwiremis_nullify_pointers(atlas)
    TYPE(uwiremis_atlas_data), INTENT(INOUT) :: atlas

    NULLIFY(atlas%bfemis_flag)
    NULLIFY(atlas%bfemis_lut)
    NULLIFY(atlas%pca_coef)
    NULLIFY(atlas%pca_offs)
    NULLIFY(atlas%pca_sfac)
    NULLIFY(atlas%pcu)
    NULLIFY(atlas%cov_emis_lut)
    NULLIFY(atlas%cov_emis)
    NULLIFY(atlas%pcu_int)
    NULLIFY(atlas%pcm_int)
    NULLIFY(atlas%sice_em_int)
    NULLIFY(atlas%snow_em_int)
    NULLIFY(atlas%cov_emis_int)

    NULLIFY(atlas%igbp)
    NULLIFY(atlas%p1d)
    NULLIFY(atlas%p2d)
    NULLIFY(atlas%p3d)
    NULLIFY(atlas%p1n)
    NULLIFY(atlas%p2n)
    NULLIFY(atlas%p3n)
    NULLIFY(atlas%p1d_int)
    NULLIFY(atlas%p2d_int)
    NULLIFY(atlas%p3d_int)
    NULLIFY(atlas%p1n_int)
    NULLIFY(atlas%p2n_int)
    NULLIFY(atlas%p3n_int)

  END SUBROUTINE rttov_uwiremis_nullify_pointers

END MODULE mod_uwiremis_atlas
