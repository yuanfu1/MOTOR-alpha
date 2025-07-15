! Description:
!> @file
!!   Subroutines for visible/near-IR BRDF atlas
!
!> @brief
!!   Subroutines for visible/near-IR BRDF atlas
!!
!! @details
!!   It is intended that this atlas be used via the RTTOV interface
!!   rather than by calling these subroutines directly.
!!
!!   Vidot, J. and E. Borbas, 2013
!!   Land surface VIS/NIR BRDF atlas for RTTOV-11: Model and Validation against
!!   SEVIRI Land SAF Albedo product. Q.J.R.M.S. 140, 2186-2196
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
MODULE mod_brdf_atlas

  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0      03/01/2012  Based on UW IR atlas code
  !  1.1      12/05/2015  -Add hdf5 atlas files read capability (J Vidot)
  !                       -Remove snow reflectance spectra and
  !                        replace by PC calculation for snow surfaces
  !                       -Replace all NetCDF KINDs by RTTOV KINDs

#include "throw.h"

#ifdef _RTTOV_HDF
#undef _RTTOV_NETCDF
#endif

  USE parkind1, ONLY : &
    jpim, &   ! 32-bit int
    jpis, &   ! 16-bit
    jprb, &   ! 64-bit real
    jprm, &   ! 32-bit
    jplm      ! logical

  USE rttov_solar_refl_mod, ONLY : &
    coastal_waters_ref, &
    ocean_waters_ref

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
  PUBLIC :: brdf_atlas_data,              &
            rttov_visnirbrdf_init,        &
            rttov_visnirbrdf,             &
            rttov_visnirbrdf_close_atlas, &
            brdfmodis_gridres               ! This can be useful for external programs

  ! Atlas constants

  INTEGER(KIND=jpim), PARAMETER :: numpcs = 6
  INTEGER(KIND=jpim), PARAMETER :: hngpnts = 7
  INTEGER(KIND=jpim), PARAMETER :: numwave = 2101

  INTEGER(KIND=jpim), PARAMETER :: db_ver_year = 2007

  REAL(KIND=jprb),    PARAMETER :: default_std = 0.05_jprb ! default standard deviation

  INTEGER(KIND=jpim), PARAMETER :: brdfmodis_gridres = 100_jprb     ! 0.1 deg
  INTEGER(KIND=jpim), PARAMETER :: brdfmodis_ygrid1 = 89950_jprb    ! 89.95 deg
  INTEGER(KIND=jpim), PARAMETER :: brdfmodis_xgrid1 = -179950_jprb  ! -179.95 deg

  REAL(KIND=jprb),    PARAMETER :: hkod = -999._jprb ! Missing data

  !> Data type for BRDF atlas data
  TYPE brdf_atlas_data
    PRIVATE

    LOGICAL(KIND=jplm) :: single_inst

    INTEGER(KIND=jpim) :: nb_lats
    INTEGER(KIND=jpim) :: nb_lons
    INTEGER(KIND=jpim) :: nb_pack

    ! Atlas data loaded by initialisation routine

    INTEGER(KIND=jpis), POINTER :: brdfmodis_flag(:,:)  ! dims are (nb_lats,nb_lons)
    INTEGER(KIND=jpim), POINTER :: brdfmodis_lut(:,:)   ! dims are (nb_lats,nb_lons)
    INTEGER(KIND=jpis), POINTER :: brdfmodis(:,:,:)     ! dims are (nb_kernel,nb_pack,hngpnts)

    REAL(KIND=jprb),    POINTER :: D(:,:), D_snow(:,:)
    REAL(KIND=jprm),    POINTER :: pcu(:,:), pcu_snow(:,:) ! pcu need to be read in from flat file

    ! Data to allow more efficient memory usage (convert integers to reals only at point of use)
    ! size is hngpnts

    REAL(KIND=jprm), POINTER :: sfac_fiso(:) ! iso parameter scale factors
    REAL(KIND=jprm), POINTER :: offs_fiso(:) ! iso parameter offsets
    REAL(KIND=jprm), POINTER :: sfac_fvol(:) ! vol parameter scale factors
    REAL(KIND=jprm), POINTER :: offs_fvol(:) ! vol parameter offsets
    REAL(KIND=jprm), POINTER :: sfac_fgeo(:) ! geo parameter scale factors
    REAL(KIND=jprm), POINTER :: offs_fgeo(:) ! geo parameter offsets

    ! Arrays to hold hsr data interpolated onto instrument wavenumbers

    INTEGER(KIND=jpim) :: platform_id
    INTEGER(KIND=jpim) :: sat_id
    INTEGER(KIND=jpim) :: inst_id
    INTEGER(KIND=jpim) :: ncoefchans
    REAL(KIND=jprb), POINTER :: coastal_waters_ref_int(:), ocean_waters_ref_int(:)
    REAL(KIND=jprb), POINTER :: pcu_int(:,:), pcu_int_snow(:,:)
    REAL(KIND=jprb), POINTER :: pcm_int(:), pcm_int_snow(:)

  END TYPE brdf_atlas_data


  ! Atlas data contained in this file - shared, does not change

  INTEGER(KIND=jpim) :: i
  REAL(KIND=jprb) :: pcm(numwave),pcm_snow(numwave)
  REAL(KIND=jprb) :: hsr_wavenum(numwave)
  REAL(KIND=jprb) :: pcu_modres(hngpnts,hngpnts),pcu_modres_snow(hngpnts,hngpnts)
  REAL(KIND=jprb) :: pcm_modres(hngpnts),pcm_modres_snow(hngpnts)

  DATA pcm_modres / &
  39.8051262_jprb,  49.6526299_jprb,  52.5679855_jprb,  49.9688873_jprb,  39.3762932_jprb,  &
  36.0293274_jprb,  28.9127731_jprb /

  DATA pcm_modres_snow / &
    4.4675913_jprb,   8.9177074_jprb,  46.8297043_jprb,  86.3503723_jprb,  92.3513489_jprb,  &
   92.8217773_jprb,  92.4939270_jprb/

!   number of PC is the first dimension and spectra is the second for pcu_modres (eigenvectors)
  DATA pcu_modres / &
  -0.0153509_jprb,  -0.0174285_jprb,  -0.0550101_jprb,  -0.0456830_jprb,  -0.0492308_jprb,  &
  -0.0334872_jprb,   0.0433589_jprb,  -0.0207704_jprb,   0.0020112_jprb,  -0.0397478_jprb,  &
  -0.0045270_jprb,   0.0223504_jprb,   0.0344748_jprb,   0.0374810_jprb,  -0.0216471_jprb,  &
   0.0288548_jprb,  -0.0150199_jprb,   0.0039704_jprb,   0.0347403_jprb,  -0.0112829_jprb,  &
   0.0245785_jprb,  -0.0213760_jprb,   0.0332810_jprb,   0.0133695_jprb,  -0.0124754_jprb,  &
  -0.0246968_jprb,   0.0084204_jprb,  -0.0061599_jprb,  -0.0256441_jprb,  -0.0138075_jprb,  &
  -0.0049740_jprb,   0.0402420_jprb,  -0.0040104_jprb,  -0.0184281_jprb,  -0.0085442_jprb,  &
  -0.0232838_jprb,  -0.0140538_jprb,   0.0088426_jprb,   0.0090562_jprb,  -0.0227179_jprb,  &
   0.0317618_jprb,   0.0141060_jprb,  -0.0217532_jprb,  -0.0192157_jprb,   0.0172922_jprb,  &
  -0.0152554_jprb,   0.0092314_jprb,   0.0020664_jprb,  -0.0013151_jprb /

  DATA pcu_modres_snow / &
   -0.0072172_jprb,  -0.0072680_jprb,  -0.0307175_jprb,  -0.0153461_jprb,   0.0000093_jprb,  &
   -0.0321871_jprb,  -0.0121987_jprb,  -0.0135923_jprb,  -0.0134654_jprb,  -0.0483233_jprb,  &
   -0.0191789_jprb,   0.0000021_jprb,  -0.0306654_jprb,  -0.0061526_jprb,  -0.0399308_jprb,  &
   -0.0359425_jprb,  -0.0134985_jprb,   0.0239141_jprb,  -0.0000455_jprb,   0.0520646_jprb,  &
    0.0099703_jprb,  -0.0202572_jprb,  -0.0062896_jprb,   0.0219874_jprb,  -0.0219007_jprb,  &
    0.0001583_jprb,  -0.0178375_jprb,   0.0249178_jprb,  -0.0165172_jprb,   0.0146778_jprb,  &
    0.0041047_jprb,  -0.0320968_jprb,  -0.0002225_jprb,   0.0182142_jprb,  -0.0091505_jprb,  &
   -0.0173791_jprb,   0.0217582_jprb,  -0.0016791_jprb,  -0.0093627_jprb,  -0.0008085_jprb,  &
    0.0106018_jprb,  -0.0318982_jprb,  -0.0189624_jprb,   0.0269639_jprb,  -0.0053196_jprb,  &
    0.0186084_jprb,  -0.0014913_jprb,  -0.0077672_jprb,  -0.0041959_jprb/

  DATA (pcm(i), i=1, 1000) / &
  24.8277779_jprb,  25.0408955_jprb,  25.3756046_jprb,  25.7037067_jprb,  26.0310898_jprb,  &
  26.3604908_jprb,  26.7347832_jprb,  27.1387291_jprb,  27.5377216_jprb,  27.9030590_jprb,  &
  28.2610722_jprb,  28.6590042_jprb,  29.0834980_jprb,  29.4976254_jprb,  29.8795204_jprb,  &
  30.1842937_jprb,  30.4282646_jprb,  30.6161823_jprb,  30.7734871_jprb,  30.9610100_jprb,  &
  31.2550640_jprb,  31.5511723_jprb,  31.6682377_jprb,  31.6455383_jprb,  31.3637829_jprb,  &
  31.0997677_jprb,  30.9445992_jprb,  30.9195175_jprb,  31.0114574_jprb,  31.1538925_jprb,  &
  31.3271103_jprb,  31.5918617_jprb,  31.9364128_jprb,  32.5347023_jprb,  33.1924286_jprb,  &
  33.7807274_jprb,  34.3048897_jprb,  34.8109322_jprb,  35.2768326_jprb,  35.6295013_jprb,  &
  35.8979759_jprb,  36.0899544_jprb,  36.2127190_jprb,  36.3239937_jprb,  36.4798775_jprb,  &
  36.7381706_jprb,  37.0174065_jprb,  37.2550430_jprb,  37.3141747_jprb,  37.1758308_jprb,  &
  36.7321091_jprb,  36.0366364_jprb,  35.4615135_jprb,  35.1427155_jprb,  35.2301254_jprb,  &
  35.6141052_jprb,  35.9188957_jprb,  36.1277504_jprb,  36.3305817_jprb,  36.5077209_jprb,  &
  36.6811256_jprb,  36.8497658_jprb,  37.0669022_jprb,  37.4394073_jprb,  37.8410454_jprb,  &
  38.2695389_jprb,  38.6767044_jprb,  38.9945450_jprb,  39.2578239_jprb,  39.4313698_jprb,  &
  39.5493889_jprb,  39.6527290_jprb,  39.7408524_jprb,  39.8403931_jprb,  39.9821548_jprb,  &
  40.1356277_jprb,  40.2465286_jprb,  40.3448677_jprb,  40.4148369_jprb,  40.4683533_jprb,  &
  40.5087280_jprb,  40.5428581_jprb,  40.5748482_jprb,  40.5822258_jprb,  40.5920639_jprb,  &
  40.5929642_jprb,  40.5849724_jprb,  40.5702782_jprb,  40.5456543_jprb,  40.5149002_jprb,  &
  40.4415970_jprb,  40.3617287_jprb,  40.2473755_jprb,  40.1151505_jprb,  39.9589920_jprb,  &
  39.7548714_jprb,  39.5531769_jprb,  39.3287430_jprb,  39.0908279_jprb,  38.8686562_jprb,  &
  38.6769829_jprb,  38.4995613_jprb,  38.3691139_jprb,  38.2617950_jprb,  38.1423302_jprb,  &
  38.0185242_jprb,  37.8886223_jprb,  37.7509346_jprb,  37.6172295_jprb,  37.4692764_jprb,  &
  37.3169022_jprb,  37.1636925_jprb,  37.0170670_jprb,  36.8825645_jprb,  36.7567787_jprb,  &
  36.6402168_jprb,  36.5216675_jprb,  36.3995590_jprb,  36.2845993_jprb,  36.1637955_jprb,  &
  36.0548248_jprb,  35.9530907_jprb,  35.8764610_jprb,  36.0405540_jprb,  36.3078537_jprb,  &
  36.6830826_jprb,  37.5555115_jprb,  38.4996796_jprb,  39.5082893_jprb,  40.6039925_jprb,  &
  41.6787682_jprb,  42.6532631_jprb,  43.5256271_jprb,  44.3674774_jprb,  45.0891762_jprb,  &
  45.7204056_jprb,  46.3044815_jprb,  46.6977043_jprb,  47.0319710_jprb,  47.3322449_jprb,  &
  47.5031433_jprb,  47.6438637_jprb,  47.7786484_jprb,  47.8841782_jprb,  47.9842186_jprb,  &
  48.0831642_jprb,  48.1158714_jprb,  48.1428604_jprb,  48.1731644_jprb,  48.2052116_jprb,  &
  48.2262611_jprb,  48.2345009_jprb,  48.2307587_jprb,  48.2203102_jprb,  48.1980476_jprb,  &
  48.1617775_jprb,  48.1188393_jprb,  48.0772324_jprb,  48.0221825_jprb,  47.9725037_jprb,  &
  47.8647423_jprb,  47.7546425_jprb,  47.6628494_jprb,  47.6102982_jprb,  47.5757103_jprb,  &
  47.5605965_jprb,  47.5695610_jprb,  47.5893745_jprb,  47.6571159_jprb,  47.7293015_jprb,  &
  47.8073463_jprb,  47.8942604_jprb,  47.9880371_jprb,  48.0909271_jprb,  48.1990356_jprb,  &
  48.3001366_jprb,  48.3927689_jprb,  48.4773979_jprb,  48.5763664_jprb,  48.6816711_jprb,  &
  48.7880669_jprb,  48.8990593_jprb,  49.0132561_jprb,  49.1233025_jprb,  49.2419930_jprb,  &
  49.3548508_jprb,  49.4367981_jprb,  49.5028229_jprb,  49.5608177_jprb,  49.6137123_jprb,  &
  49.6512108_jprb,  49.6875305_jprb,  49.7134285_jprb,  49.7348328_jprb,  49.7571220_jprb,  &
  49.7906685_jprb,  49.8226013_jprb,  49.8518486_jprb,  49.8824654_jprb,  49.9062805_jprb,  &
  49.9288368_jprb,  49.9525146_jprb,  49.9775429_jprb,  49.9752960_jprb,  49.9725151_jprb,  &
  49.9668655_jprb,  49.9510574_jprb,  49.9350204_jprb,  49.9062233_jprb,  49.8749199_jprb,  &
  49.8329430_jprb,  49.7877769_jprb,  49.7458458_jprb,  49.7070007_jprb,  49.6658325_jprb,  &
  49.6287918_jprb,  49.5949745_jprb,  49.5527306_jprb,  49.5147476_jprb,  49.4577179_jprb,  &
  49.4103966_jprb,  49.3577080_jprb,  49.3042564_jprb,  49.2495346_jprb,  49.1900711_jprb,  &
  49.1349220_jprb,  49.0783005_jprb,  49.0104485_jprb,  48.9424515_jprb,  48.8712959_jprb,  &
  48.7976608_jprb,  48.7201958_jprb,  48.6492691_jprb,  48.5846863_jprb,  48.5148926_jprb,  &
  48.4375496_jprb,  48.3640823_jprb,  48.2843971_jprb,  48.2191467_jprb,  48.1516533_jprb,  &
  48.0809593_jprb,  48.0056038_jprb,  47.9305763_jprb,  47.8318100_jprb,  47.7247353_jprb,  &
  47.6254768_jprb,  47.5232315_jprb,  47.4170990_jprb,  47.3067856_jprb,  47.1959763_jprb,  &
  47.1025696_jprb,  47.0158081_jprb,  46.9364243_jprb,  46.8515549_jprb,  46.7669792_jprb,  &
  46.6809807_jprb,  46.5889931_jprb,  46.4778671_jprb,  46.3589363_jprb,  46.2331390_jprb,  &
  46.1061630_jprb,  45.9807281_jprb,  45.8536797_jprb,  45.7196045_jprb,  45.5729408_jprb,  &
  45.4256516_jprb,  45.2711411_jprb,  45.1085663_jprb,  44.9622421_jprb,  44.8323746_jprb,  &
  44.6884499_jprb,  44.5197906_jprb,  44.3317146_jprb,  44.0937538_jprb,  43.8393288_jprb,  &
  43.5828819_jprb,  43.3454704_jprb,  43.2157364_jprb,  43.1692314_jprb,  43.1494408_jprb,  &
  43.1466866_jprb,  43.1388931_jprb,  43.0929909_jprb,  43.0425720_jprb,  42.9774513_jprb,  &
  42.9205513_jprb,  42.8812256_jprb,  42.8628807_jprb,  42.8499832_jprb,  42.8666153_jprb,  &
  42.8766289_jprb,  42.8651009_jprb,  42.8294296_jprb,  42.7809982_jprb,  42.6996040_jprb,  &
  42.6092682_jprb,  42.4630814_jprb,  42.2992706_jprb,  42.1220589_jprb,  41.9337540_jprb,  &
  41.7471809_jprb,  41.5761147_jprb,  41.4013367_jprb,  41.2029266_jprb,  40.9854431_jprb,  &
  40.7861519_jprb,  40.6723328_jprb,  40.6026917_jprb,  40.7962875_jprb,  41.2006149_jprb,  &
  41.6306000_jprb,  42.0552292_jprb,  42.4793282_jprb,  42.8819809_jprb,  43.2628098_jprb,  &
  43.6393890_jprb,  43.9990768_jprb,  44.3819275_jprb,  44.8599586_jprb,  45.4965630_jprb,  &
  46.1477737_jprb,  46.7114563_jprb,  47.2458572_jprb,  47.7549171_jprb,  48.1832123_jprb,  &
  48.5828934_jprb,  48.9394531_jprb,  49.2532997_jprb,  49.5449371_jprb,  49.7617188_jprb,  &
  49.9556465_jprb,  50.1164284_jprb,  50.2043686_jprb,  50.2867813_jprb,  50.3782043_jprb,  &
  50.4792099_jprb,  50.5746956_jprb,  50.6758156_jprb,  50.7880707_jprb,  50.8856621_jprb,  &
  50.9989624_jprb,  51.1025085_jprb,  51.2065086_jprb,  51.3081245_jprb,  51.4042892_jprb,  &
  51.4893990_jprb,  51.5658188_jprb,  51.6339455_jprb,  51.7216644_jprb,  51.8101959_jprb,  &
  51.8976097_jprb,  51.9952774_jprb,  52.0927658_jprb,  52.1877747_jprb,  52.2839622_jprb,  &
  52.3732262_jprb,  52.4584694_jprb,  52.5331650_jprb,  52.6083336_jprb,  52.6867561_jprb,  &
  52.7561340_jprb,  52.8231392_jprb,  52.8769035_jprb,  52.9340477_jprb,  52.9769707_jprb,  &
  53.0189934_jprb,  53.0571861_jprb,  53.1026077_jprb,  53.1264534_jprb,  53.1512184_jprb,  &
  53.1757889_jprb,  53.1957359_jprb,  53.2113342_jprb,  53.2238350_jprb,  53.2309532_jprb,  &
  53.2288017_jprb,  53.2214279_jprb,  53.2135010_jprb,  53.2051048_jprb,  53.1969490_jprb,  &
  53.1745415_jprb,  53.1487122_jprb,  53.1208115_jprb,  53.0874443_jprb,  53.0537109_jprb,  &
  53.0197678_jprb,  52.9902916_jprb,  52.9605255_jprb,  52.9354286_jprb,  52.9177933_jprb,  &
  52.9032173_jprb,  52.8898621_jprb,  52.8802643_jprb,  52.8731766_jprb,  52.8647499_jprb,  &
  52.8530045_jprb,  52.8380737_jprb,  52.8203926_jprb,  52.8020592_jprb,  52.7679634_jprb,  &
  52.7333336_jprb,  52.7019920_jprb,  52.6730728_jprb,  52.6412010_jprb,  52.6088791_jprb,  &
  52.5843277_jprb,  52.5531616_jprb,  52.5271034_jprb,  52.4952164_jprb,  52.4646378_jprb,  &
  52.4378052_jprb,  52.4118347_jprb,  52.3874474_jprb,  52.3652267_jprb,  52.3381195_jprb,  &
  52.3142014_jprb,  52.2914925_jprb,  52.2736092_jprb,  52.2451439_jprb,  52.2222404_jprb,  &
  52.1886978_jprb,  52.1623955_jprb,  52.1310196_jprb,  52.1028595_jprb,  52.0730438_jprb,  &
  52.0445671_jprb,  52.0113258_jprb,  51.9817390_jprb,  51.9564095_jprb,  51.9228745_jprb,  &
  51.8925934_jprb,  51.8659782_jprb,  51.8516502_jprb,  51.8358459_jprb,  51.8187065_jprb,  &
  51.8114624_jprb,  51.8038597_jprb,  51.7975960_jprb,  51.7957802_jprb,  51.7862206_jprb,  &
  51.7749062_jprb,  51.7648163_jprb,  51.7499771_jprb,  51.7300301_jprb,  51.7053070_jprb,  &
  51.6875572_jprb,  51.6764526_jprb,  51.6724510_jprb,  51.6580887_jprb,  51.6512489_jprb,  &
  51.6510277_jprb,  51.6447945_jprb,  51.6359177_jprb,  51.6288338_jprb,  51.6105728_jprb,  &
  51.5943718_jprb,  51.5859909_jprb,  51.5827827_jprb,  51.5809479_jprb,  51.5743904_jprb,  &
  51.5832405_jprb,  51.5928879_jprb,  51.6025124_jprb,  51.6131096_jprb,  51.6350632_jprb,  &
  51.6606941_jprb,  51.6791229_jprb,  51.7093048_jprb,  51.7432480_jprb,  51.7822037_jprb,  &
  51.8244705_jprb,  51.8735008_jprb,  51.9337273_jprb,  52.0021095_jprb,  52.0662193_jprb,  &
  52.1219902_jprb,  52.1815948_jprb,  52.2372551_jprb,  52.2917328_jprb,  52.3450737_jprb,  &
  52.3958054_jprb,  52.4485970_jprb,  52.4963264_jprb,  52.5285416_jprb,  52.5456848_jprb,  &
  52.5682220_jprb,  52.5838470_jprb,  52.5915985_jprb,  52.6025658_jprb,  52.6067963_jprb,  &
  52.6079178_jprb,  52.6029167_jprb,  52.5942574_jprb,  52.5815468_jprb,  52.5732155_jprb,  &
  52.5633888_jprb,  52.5515327_jprb,  52.5406914_jprb,  52.5288696_jprb,  52.5188370_jprb,  &
  52.5111122_jprb,  52.5029030_jprb,  52.4930000_jprb,  52.4852066_jprb,  52.4742317_jprb,  &
  52.4616852_jprb,  52.4498520_jprb,  52.4321480_jprb,  52.4159622_jprb,  52.3998108_jprb,  &
  52.3804893_jprb,  52.3675957_jprb,  52.3594437_jprb,  52.3474426_jprb,  52.3394012_jprb,  &
  52.3224907_jprb,  52.3066368_jprb,  52.2871361_jprb,  52.2686691_jprb,  52.2519913_jprb,  &
  52.2347488_jprb,  52.2113647_jprb,  52.1953735_jprb,  52.1752243_jprb,  52.1597404_jprb,  &
  52.1403046_jprb,  52.1238899_jprb,  52.1051445_jprb,  52.0915756_jprb,  52.0763893_jprb,  &
  52.0601921_jprb,  52.0443954_jprb,  52.0285988_jprb,  52.0057831_jprb,  51.9835510_jprb,  &
  51.9618835_jprb,  51.9426765_jprb,  51.9239616_jprb,  51.9005165_jprb,  51.8794785_jprb,  &
  51.8587990_jprb,  51.8369713_jprb,  51.8125839_jprb,  51.7946358_jprb,  51.7740822_jprb,  &
  51.7539864_jprb,  51.7282219_jprb,  51.7064743_jprb,  51.6843719_jprb,  51.6619110_jprb,  &
  51.6393089_jprb,  51.6160889_jprb,  51.5870743_jprb,  51.5575371_jprb,  51.5361748_jprb,  &
  51.5095062_jprb,  51.4877815_jprb,  51.4607735_jprb,  51.4419670_jprb,  51.4202919_jprb,  &
  51.4010620_jprb,  51.3809471_jprb,  51.3631134_jprb,  51.3402901_jprb,  51.3124466_jprb,  &
  51.2802010_jprb,  51.2503929_jprb,  51.2237434_jprb,  51.1952133_jprb,  51.1688309_jprb,  &
  51.1414261_jprb,  51.1203613_jprb,  51.0967331_jprb,  51.0638733_jprb,  51.0370522_jprb,  &
  51.0052795_jprb,  50.9706154_jprb,  50.9383621_jprb,  50.9037666_jprb,  50.8628426_jprb,  &
  50.8253746_jprb,  50.7880478_jprb,  50.7480202_jprb,  50.7141953_jprb,  50.6797256_jprb,  &
  50.6510429_jprb,  50.6217155_jprb,  50.5882530_jprb,  50.5581360_jprb,  50.5372543_jprb,  &
  50.5106659_jprb,  50.4889908_jprb,  50.4598312_jprb,  50.4428482_jprb,  50.4221764_jprb,  &
  50.3889351_jprb,  50.3576393_jprb,  50.3225365_jprb,  50.2929764_jprb,  50.2579727_jprb,  &
  50.2222214_jprb,  50.1976433_jprb,  50.1746941_jprb,  50.1478615_jprb,  50.1312180_jprb,  &
  50.1038666_jprb,  50.0832787_jprb,  50.0704613_jprb,  50.0520058_jprb,  50.0398445_jprb,  &
  50.0267334_jprb,  50.0163155_jprb,  50.0069122_jprb,  49.9936562_jprb,  49.9829483_jprb,  &
  49.9692192_jprb,  49.9593430_jprb,  49.9469185_jprb,  49.9347458_jprb,  49.9212189_jprb,  &
  49.9120140_jprb,  49.9000587_jprb,  49.8828354_jprb,  49.8688316_jprb,  49.8582878_jprb,  &
  49.8416824_jprb,  49.8188858_jprb,  49.8003883_jprb,  49.7844124_jprb,  49.7631149_jprb,  &
  49.7411499_jprb,  49.7240791_jprb,  49.7277603_jprb,  49.7255249_jprb,  49.7225990_jprb,  &
  49.7235184_jprb,  49.7244034_jprb,  49.7272987_jprb,  49.7290535_jprb,  49.7340431_jprb,  &
  49.7380943_jprb,  49.7374840_jprb,  49.7377243_jprb,  49.7411957_jprb,  49.7446899_jprb,  &
  49.7512856_jprb,  49.7585487_jprb,  49.7601433_jprb,  49.7730370_jprb,  49.7852173_jprb,  &
  49.7938004_jprb,  49.8022423_jprb,  49.8105507_jprb,  49.8190346_jprb,  49.8340187_jprb,  &
  49.8599854_jprb,  49.8791466_jprb,  49.8962135_jprb,  49.9132423_jprb,  49.9302406_jprb,  &
  49.9448166_jprb,  49.9586143_jprb,  49.9748116_jprb,  49.9914742_jprb,  50.0030594_jprb,  &
  50.0160561_jprb,  50.0209503_jprb,  50.0245743_jprb,  50.0328445_jprb,  50.0443077_jprb,  &
  50.0511284_jprb,  50.0550766_jprb,  50.0624504_jprb,  50.0567436_jprb,  50.0464134_jprb,  &
  50.0464554_jprb,  50.0418129_jprb,  50.0353546_jprb,  50.0275078_jprb,  50.0182991_jprb,  &
  50.0115814_jprb,  49.9986649_jprb,  49.9969940_jprb,  49.9932785_jprb,  49.9854317_jprb,  &
  49.9782829_jprb,  49.9752884_jprb,  49.9710770_jprb,  49.9630280_jprb,  49.9549904_jprb,  &
  49.9537354_jprb,  49.9515114_jprb,  49.9452438_jprb,  49.9397163_jprb,  49.9324188_jprb,  &
  49.9341507_jprb,  49.9386711_jprb,  49.9422760_jprb,  49.9458580_jprb,  49.9523811_jprb,  &
  49.9523811_jprb,  49.9547119_jprb,  49.9586601_jprb,  49.9626236_jprb,  49.9655609_jprb,  &
  49.9674263_jprb,  49.9661255_jprb,  49.9658737_jprb,  49.9669991_jprb,  49.9691620_jprb,  &
  49.9715080_jprb,  49.9721489_jprb,  49.9716530_jprb,  49.9718475_jprb,  49.9727516_jprb,  &
  49.9733696_jprb,  49.9708176_jprb,  49.9688606_jprb,  49.9712830_jprb,  49.9716225_jprb,  &
  49.9664421_jprb,  49.9607506_jprb,  49.9544907_jprb,  49.9515419_jprb,  49.9500313_jprb,  &
  49.9477959_jprb,  49.9434090_jprb,  49.9423103_jprb,  49.9394264_jprb,  49.9348564_jprb,  &
  49.9301224_jprb,  49.9295044_jprb,  49.9240608_jprb,  49.9237061_jprb,  49.9297562_jprb,  &
  49.9353065_jprb,  49.9357147_jprb,  49.9364967_jprb,  49.9389839_jprb,  49.9436493_jprb,  &
  49.9436493_jprb,  49.9500465_jprb,  49.9561234_jprb,  49.9587288_jprb,  49.9587288_jprb,  &
  49.9617653_jprb,  49.9663429_jprb,  49.9707603_jprb,  49.9713631_jprb,  49.9714279_jprb,  &
  49.9734993_jprb,  49.9766426_jprb,  49.9784393_jprb,  49.9791260_jprb,  49.9800758_jprb,  &
  49.9815826_jprb,  49.9839630_jprb,  49.9846802_jprb,  49.9863243_jprb,  49.9880943_jprb,  &
  49.9880943_jprb,  49.9896965_jprb,  49.9894905_jprb,  49.9895668_jprb,  49.9966049_jprb,  &
  49.9956169_jprb,  49.9964828_jprb,  50.0000000_jprb,  50.0000000_jprb,  49.9968452_jprb,  &
  49.9953918_jprb,  49.9976425_jprb,  49.9982185_jprb,  49.9976425_jprb,  49.9966583_jprb,  &
  49.9960823_jprb,  49.9860840_jprb,  49.9757996_jprb,  49.9660645_jprb,  49.9563522_jprb,  &
  49.9436417_jprb,  49.9303856_jprb,  49.9178734_jprb,  49.9114265_jprb,  49.9006920_jprb,  &
  49.8869209_jprb,  49.8750648_jprb,  49.8820724_jprb,  49.8952179_jprb,  49.9087067_jprb,  &
  49.9226913_jprb,  49.9384460_jprb,  49.9551392_jprb,  49.9662971_jprb,  49.9789619_jprb,  &
  49.9923134_jprb,  50.0056419_jprb,  50.0198593_jprb,  50.0290070_jprb,  50.0301132_jprb,  &
  50.0306892_jprb,  50.0315247_jprb,  50.0325394_jprb,  50.0325394_jprb,  50.0366821_jprb,  &
  50.0424576_jprb,  50.0446434_jprb,  50.0444450_jprb,  50.0435028_jprb,  50.0429611_jprb,  &
  50.0410995_jprb,  50.0381851_jprb,  50.0338593_jprb,  50.0290680_jprb,  50.0268478_jprb,  &
  50.0284653_jprb,  50.0255661_jprb,  50.0255165_jprb,  50.0292664_jprb,  50.0288467_jprb,  &
  50.0309105_jprb,  50.0357094_jprb,  50.0378380_jprb,  50.0396957_jprb,  50.0412788_jprb,  &
  50.0423393_jprb,  50.0442009_jprb,  50.0468445_jprb,  50.0500031_jprb,  50.0515862_jprb,  &
  50.0515862_jprb,  50.0526237_jprb,  50.0519371_jprb,  50.0493240_jprb,  50.0447273_jprb,  &
  50.0386925_jprb,  50.0308914_jprb,  50.0273399_jprb,  50.0226707_jprb,  50.0154266_jprb,  &
  50.0147209_jprb,  50.0150795_jprb,  50.0150795_jprb,  50.0154724_jprb,  50.0163231_jprb,  &
  50.0183716_jprb,  50.0180206_jprb,  50.0166054_jprb,  50.0160980_jprb,  50.0155869_jprb,  &
  50.0150795_jprb,  50.0094986_jprb,  50.0067902_jprb,  50.0078011_jprb,  50.0048714_jprb,  &
  50.0004501_jprb,  49.9939041_jprb,  49.9873734_jprb,  49.9813309_jprb,  49.9783249_jprb,  &
  49.9720001_jprb,  49.9630051_jprb,  49.9660683_jprb,  49.9688110_jprb,  49.9703026_jprb,  &
  49.9748764_jprb,  49.9796295_jprb,  49.9811134_jprb,  49.9752426_jprb,  49.9638939_jprb,  &
  49.9525757_jprb,  49.9417839_jprb,  49.9319649_jprb,  49.9225235_jprb,  49.9138374_jprb,  &
  49.9099312_jprb,  49.9088898_jprb,  49.9108391_jprb,  49.9056740_jprb,  49.8997993_jprb,  &
  49.8954391_jprb,  49.8925095_jprb,  49.8905754_jprb,  49.8859024_jprb,  49.8807907_jprb,  &
  49.8750191_jprb,  49.8676910_jprb,  49.8595352_jprb,  49.8509216_jprb,  49.8418274_jprb,  &
  49.8318062_jprb,  49.8179283_jprb,  49.8022957_jprb,  49.7942238_jprb,  49.7859993_jprb,  &
  49.7774734_jprb,  49.7667694_jprb,  49.7549629_jprb,  49.7385483_jprb,  49.7234764_jprb,  &
  49.7117271_jprb,  49.7011414_jprb,  49.6912956_jprb,  49.6767616_jprb,  49.6633415_jprb,  &
  49.6554146_jprb,  49.6472588_jprb,  49.6388893_jprb,  49.6279030_jprb,  49.6166306_jprb,  &
  49.6083031_jprb,  49.5962639_jprb,  49.5782700_jprb,  49.5603027_jprb,  49.5423546_jprb,  &
  49.5231552_jprb,  49.5037918_jprb,  49.4841080_jprb,  49.4644470_jprb,  49.4448128_jprb,  &
  49.4255524_jprb,  49.4064255_jprb,  49.3850861_jprb,  49.3621750_jprb,  49.3345413_jprb,  &
  49.3073959_jprb,  49.2807388_jprb,  49.2501335_jprb,  49.2181396_jprb,  49.1827126_jprb,  &
  49.1461792_jprb,  49.1058388_jprb,  49.0671272_jprb,  49.0304718_jprb,  48.9971504_jprb,  &
  48.9654961_jprb,  48.9351006_jprb,  48.9032974_jprb,  48.8576202_jprb,  48.8112907_jprb,  &
  48.7635269_jprb,  48.7192497_jprb,  48.6782303_jprb,  48.6360245_jprb,  48.5933533_jprb,  &
  48.5451736_jprb,  48.4958916_jprb,  48.4416351_jprb,  48.3870430_jprb,  48.3316231_jprb,  &
  48.2760696_jprb,  48.2203636_jprb,  48.1568756_jprb,  48.0891685_jprb,  48.0233879_jprb,  &
  47.9580269_jprb,  47.8983803_jprb,  47.8384514_jprb,  47.7764587_jprb,  47.7113800_jprb,  &
  47.6392517_jprb,  47.5664330_jprb,  47.4927902_jprb,  47.4156647_jprb,  47.3362427_jprb,  &
  47.2579308_jprb,  47.1800041_jprb,  47.1000214_jprb,  47.0197868_jprb,  46.9422836_jprb,  &
  46.8630714_jprb,  46.7731094_jprb,  46.6837540_jprb,  46.5961456_jprb,  46.5114517_jprb,  &
  46.4316444_jprb,  46.3515358_jprb,  46.2711067_jprb,  46.1839142_jprb,  46.0920143_jprb,  &
  46.0002861_jprb,  45.9086380_jprb,  45.8165321_jprb,  45.7243080_jprb,  45.6323013_jprb,  &
  45.5403290_jprb,  45.4414749_jprb,  45.3423462_jprb,  45.2438774_jprb,  45.1460304_jprb,  &
  45.0531769_jprb,  44.9611778_jprb,  44.8730888_jprb,  44.7848625_jprb,  44.6962013_jprb,  &
  44.6074944_jprb,  44.5186729_jprb,  44.4296494_jprb,  44.3402634_jprb,  44.2506561_jprb /

  DATA (pcm(i), i=1001, 2101) / &
  44.1607094_jprb,  44.0731812_jprb,  43.9887390_jprb,  43.9055557_jprb,  43.8237686_jprb,  &
  43.7400627_jprb,  43.6544876_jprb,  43.5656395_jprb,  43.4738922_jprb,  43.3855782_jprb,  &
  43.3000793_jprb,  43.2144966_jprb,  43.1288414_jprb,  43.0444603_jprb,  42.9610329_jprb,  &
  42.8788757_jprb,  42.7976685_jprb,  42.7186127_jprb,  42.6411896_jprb,  42.5643272_jprb,  &
  42.4879112_jprb,  42.4143944_jprb,  42.3432999_jprb,  42.2731361_jprb,  42.2038078_jprb,  &
  42.1262093_jprb,  42.0402870_jprb,  41.9597473_jprb,  41.8853607_jprb,  41.8125687_jprb,  &
  41.7418900_jprb,  41.6727829_jprb,  41.6061783_jprb,  41.5405464_jprb,  41.4768372_jprb,  &
  41.4141502_jprb,  41.3540878_jprb,  41.2948189_jprb,  41.2383842_jprb,  41.1822052_jprb,  &
  41.1274643_jprb,  41.0734787_jprb,  41.0284386_jprb,  40.9833984_jprb,  40.9400063_jprb,  &
  40.8966179_jprb,  40.8550797_jprb,  40.8137131_jprb,  40.7774200_jprb,  40.7422295_jprb,  &
  40.7087517_jprb,  40.6759262_jprb,  40.6424599_jprb,  40.6086006_jprb,  40.5786285_jprb,  &
  40.5523720_jprb,  40.5258446_jprb,  40.4989166_jprb,  40.4718018_jprb,  40.4442024_jprb,  &
  40.4175949_jprb,  40.3967705_jprb,  40.3760033_jprb,  40.3600616_jprb,  40.3441162_jprb,  &
  40.3246727_jprb,  40.3047104_jprb,  40.2845230_jprb,  40.2642479_jprb,  40.2461090_jprb,  &
  40.2295685_jprb,  40.2119942_jprb,  40.1929398_jprb,  40.1752472_jprb,  40.1617279_jprb,  &
  40.1482124_jprb,  40.1347351_jprb,  40.1212540_jprb,  40.1068344_jprb,  40.0923080_jprb,  &
  40.0796089_jprb,  40.0676575_jprb,  40.0519829_jprb,  40.0328445_jprb,  40.0155716_jprb,  &
  40.0022507_jprb,  39.9889793_jprb,  39.9760551_jprb,  39.9631348_jprb,  39.9518661_jprb,  &
  39.9407730_jprb,  39.9260254_jprb,  39.9096107_jprb,  39.8932228_jprb,  39.8768578_jprb,  &
  39.8610001_jprb,  39.8468132_jprb,  39.8326263_jprb,  39.8258629_jprb,  39.8191414_jprb,  &
  39.8177261_jprb,  39.8180771_jprb,  39.8131714_jprb,  39.8029785_jprb,  39.7918663_jprb,  &
  39.7778473_jprb,  39.7638283_jprb,  39.7553101_jprb,  39.7469215_jprb,  39.7385559_jprb,  &
  39.7301941_jprb,  39.7209587_jprb,  39.7105408_jprb,  39.6998329_jprb,  39.6873703_jprb,  &
  39.6749039_jprb,  39.6656952_jprb,  39.6570625_jprb,  39.6512833_jprb,  39.6478424_jprb,  &
  39.6415672_jprb,  39.6264725_jprb,  39.6113739_jprb,  39.5991821_jprb,  39.5872116_jprb,  &
  39.5752602_jprb,  39.5633278_jprb,  39.5532265_jprb,  39.5477867_jprb,  39.5423470_jprb,  &
  39.5340309_jprb,  39.5255585_jprb,  39.5175095_jprb,  39.5097389_jprb,  39.5021439_jprb,  &
  39.4950714_jprb,  39.4879951_jprb,  39.4809418_jprb,  39.4738884_jprb,  39.4654083_jprb,  &
  39.4556961_jprb,  39.4457703_jprb,  39.4347534_jprb,  39.4237404_jprb,  39.4153824_jprb,  &
  39.4077301_jprb,  39.4006233_jprb,  39.3943214_jprb,  39.3879433_jprb,  39.3757095_jprb,  &
  39.3634720_jprb,  39.3536835_jprb,  39.3454399_jprb,  39.3374062_jprb,  39.3301773_jprb,  &
  39.3229485_jprb,  39.3178062_jprb,  39.3132172_jprb,  39.3081551_jprb,  39.3022766_jprb,  &
  39.2963982_jprb,  39.2932854_jprb,  39.2903557_jprb,  39.2858047_jprb,  39.2796402_jprb,  &
  39.2733307_jprb,  39.2645950_jprb,  39.2558594_jprb,  39.2483025_jprb,  39.2415314_jprb,  &
  39.2345772_jprb,  39.2265396_jprb,  39.2185020_jprb,  39.2100487_jprb,  39.2013931_jprb,  &
  39.1928711_jprb,  39.1848831_jprb,  39.1768951_jprb,  39.1696091_jprb,  39.1626015_jprb,  &
  39.1544037_jprb,  39.1420212_jprb,  39.1296387_jprb,  39.1166000_jprb,  39.1033058_jprb,  &
  39.0901566_jprb,  39.0775375_jprb,  39.0649185_jprb,  39.0551796_jprb,  39.0466881_jprb,  &
  39.0380936_jprb,  39.0289993_jprb,  39.0199089_jprb,  39.0114403_jprb,  39.0033150_jprb,  &
  38.9953842_jprb,  38.9891548_jprb,  38.9829254_jprb,  38.9763603_jprb,  38.9695282_jprb,  &
  38.9626999_jprb,  38.9484596_jprb,  38.9342232_jprb,  38.9233665_jprb,  38.9168854_jprb,  &
  38.9104080_jprb,  38.9007607_jprb,  38.8906097_jprb,  38.8807297_jprb,  38.8715324_jprb,  &
  38.8623390_jprb,  38.8516769_jprb,  38.8403702_jprb,  38.8293648_jprb,  38.8211403_jprb,  &
  38.8129158_jprb,  38.8053017_jprb,  38.7983170_jprb,  38.7913361_jprb,  38.7822495_jprb,  &
  38.7728691_jprb,  38.7637253_jprb,  38.7552795_jprb,  38.7468300_jprb,  38.7378387_jprb,  &
  38.7285194_jprb,  38.7192001_jprb,  38.7160988_jprb,  38.7131042_jprb,  38.7093086_jprb,  &
  38.7039299_jprb,  38.6985550_jprb,  38.6919746_jprb,  38.6848297_jprb,  38.6776543_jprb,  &
  38.6681557_jprb,  38.6586609_jprb,  38.6500778_jprb,  38.6432724_jprb,  38.6364708_jprb,  &
  38.6283264_jprb,  38.6194801_jprb,  38.6106377_jprb,  38.6043701_jprb,  38.5981979_jprb,  &
  38.5919533_jprb,  38.5855103_jprb,  38.5790672_jprb,  38.5736046_jprb,  38.5689354_jprb,  &
  38.5642624_jprb,  38.5462418_jprb,  38.5255814_jprb,  38.5055542_jprb,  38.4907646_jprb,  &
  38.4759712_jprb,  38.4659157_jprb,  38.4638939_jprb,  38.4618683_jprb,  38.4560661_jprb,  &
  38.4480019_jprb,  38.4399338_jprb,  38.4268913_jprb,  38.4131088_jprb,  38.3994522_jprb,  &
  38.3868599_jprb,  38.3742676_jprb,  38.3613091_jprb,  38.3476181_jprb,  38.3339272_jprb,  &
  38.3180656_jprb,  38.3004417_jprb,  38.2828140_jprb,  38.2699661_jprb,  38.2586327_jprb,  &
  38.2472992_jprb,  38.2403336_jprb,  38.2335548_jprb,  38.2267380_jprb,  38.2197037_jprb,  &
  38.2126656_jprb,  38.2033386_jprb,  38.1893158_jprb,  38.1752892_jprb,  38.1601753_jprb,  &
  38.1439590_jprb,  38.1277466_jprb,  38.1148758_jprb,  38.1037292_jprb,  38.0925865_jprb,  &
  38.0762939_jprb,  38.0587997_jprb,  38.0413055_jprb,  38.0249214_jprb,  38.0085945_jprb,  &
  37.9921532_jprb,  37.9745026_jprb,  37.9568520_jprb,  37.9401245_jprb,  37.9269295_jprb,  &
  37.9137383_jprb,  37.9009972_jprb,  37.8892174_jprb,  37.8774376_jprb,  37.8631363_jprb,  &
  37.8453941_jprb,  37.8276482_jprb,  37.8100739_jprb,  37.7926598_jprb,  37.7752495_jprb,  &
  37.7593307_jprb,  37.7444191_jprb,  37.7295074_jprb,  37.7139015_jprb,  37.6979599_jprb,  &
  37.6820183_jprb,  37.6686821_jprb,  37.6562958_jprb,  37.6439095_jprb,  37.6292305_jprb,  &
  37.6139336_jprb,  37.5986366_jprb,  37.5867119_jprb,  37.5754814_jprb,  37.5642471_jprb,  &
  37.5502930_jprb,  37.5359001_jprb,  37.5215073_jprb,  37.5048103_jprb,  37.4878120_jprb,  &
  37.4708138_jprb,  37.4529190_jprb,  37.4349174_jprb,  37.4169197_jprb,  37.4013443_jprb,  &
  37.3860474_jprb,  37.3707466_jprb,  37.3536263_jprb,  37.3362770_jprb,  37.3189278_jprb,  &
  37.3034554_jprb,  37.2882614_jprb,  37.2730713_jprb,  37.2579231_jprb,  37.2427826_jprb,  &
  37.2276421_jprb,  37.2088013_jprb,  37.1890335_jprb,  37.1692619_jprb,  37.1514854_jprb,  &
  37.1343765_jprb,  37.1172676_jprb,  37.1028595_jprb,  37.0896873_jprb,  37.0765114_jprb,  &
  37.0611534_jprb,  37.0444221_jprb,  37.0276909_jprb,  37.0130348_jprb,  37.0002098_jprb,  &
  36.9873848_jprb,  36.9724579_jprb,  36.9548225_jprb,  36.9371872_jprb,  36.9202461_jprb,  &
  36.9047089_jprb,  36.8891716_jprb,  36.8728256_jprb,  36.8535385_jprb,  36.8342476_jprb,  &
  36.8150330_jprb,  36.7965698_jprb,  36.7781105_jprb,  36.7596474_jprb,  36.7407684_jprb,  &
  36.7218666_jprb,  36.7029648_jprb,  36.6841202_jprb,  36.6652870_jprb,  36.6464500_jprb,  &
  36.6266975_jprb,  36.6064301_jprb,  36.5861626_jprb,  36.5639420_jprb,  36.5395088_jprb,  &
  36.5150719_jprb,  36.4915161_jprb,  36.4701500_jprb,  36.4487839_jprb,  36.4279404_jprb,  &
  36.4123459_jprb,  36.3967514_jprb,  36.3811531_jprb,  36.3616753_jprb,  36.3416977_jprb,  &
  36.3217163_jprb,  36.2965317_jprb,  36.2687569_jprb,  36.2409821_jprb,  36.2143250_jprb,  &
  36.1891022_jprb,  36.1638756_jprb,  36.1390076_jprb,  36.1155815_jprb,  36.0921555_jprb,  &
  36.0687294_jprb,  36.0453835_jprb,  36.0220413_jprb,  35.9986992_jprb,  35.9755783_jprb,  &
  35.9525642_jprb,  35.9295502_jprb,  35.9059792_jprb,  35.8815994_jprb,  35.8572197_jprb,  &
  35.8324814_jprb,  35.8050613_jprb,  35.7776413_jprb,  35.7502251_jprb,  35.7244530_jprb,  &
  35.6990509_jprb,  35.6736488_jprb,  35.6487770_jprb,  35.6244240_jprb,  35.6000710_jprb,  &
  35.5758171_jprb,  35.5520287_jprb,  35.5282364_jprb,  35.5044479_jprb,  35.4815331_jprb,  &
  35.4587784_jprb,  35.4360237_jprb,  35.4116631_jprb,  35.3856850_jprb,  35.3597031_jprb,  &
  35.3340645_jprb,  35.3105316_jprb,  35.2870026_jprb,  35.2634735_jprb,  35.2387543_jprb,  &
  35.2136688_jprb,  35.1885834_jprb,  35.1640663_jprb,  35.1404762_jprb,  35.1168823_jprb,  &
  35.0932922_jprb,  35.0679588_jprb,  35.0425949_jprb,  35.0172272_jprb,  34.9915161_jprb,  &
  34.9655495_jprb,  34.9395828_jprb,  34.9134903_jprb,  34.8866959_jprb,  34.8598976_jprb,  &
  34.8331032_jprb,  34.8068695_jprb,  34.7808609_jprb,  34.7548561_jprb,  34.7286949_jprb,  &
  34.7020950_jprb,  34.6754990_jprb,  34.6489029_jprb,  34.6229324_jprb,  34.5971184_jprb,  &
  34.5713043_jprb,  34.5450020_jprb,  34.5176888_jprb,  34.4903793_jprb,  34.4630699_jprb,  &
  34.4367981_jprb,  34.4107246_jprb,  34.3846474_jprb,  34.3579903_jprb,  34.3302078_jprb,  &
  34.3024216_jprb,  34.2746353_jprb,  34.2480545_jprb,  34.2217255_jprb,  34.1953964_jprb,  &
  34.1690941_jprb,  34.1428642_jprb,  34.1166306_jprb,  34.0904007_jprb,  34.0639038_jprb,  &
  34.0373268_jprb,  34.0107460_jprb,  33.9838638_jprb,  33.9558258_jprb,  33.9277878_jprb,  &
  33.8997536_jprb,  33.8666687_jprb,  33.8307571_jprb,  33.7948456_jprb,  33.7593575_jprb,  &
  33.7304268_jprb,  33.7014999_jprb,  33.6725693_jprb,  33.6459503_jprb,  33.6219673_jprb,  &
  33.5979881_jprb,  33.5740089_jprb,  33.5482216_jprb,  33.5221443_jprb,  33.4960632_jprb,  &
  33.4684982_jprb,  33.4361839_jprb,  33.4038696_jprb,  33.3715553_jprb,  33.3411522_jprb,  &
  33.3120079_jprb,  33.2828636_jprb,  33.2537193_jprb,  33.2248878_jprb,  33.1960678_jprb,  &
  33.1672516_jprb,  33.1380882_jprb,  33.1080856_jprb,  33.0780830_jprb,  33.0480843_jprb,  &
  33.0174904_jprb,  32.9865303_jprb,  32.9555664_jprb,  32.9246025_jprb,  32.8961258_jprb,  &
  32.8678513_jprb,  32.8395767_jprb,  32.8113747_jprb,  32.7834244_jprb,  32.7554741_jprb,  &
  32.7275200_jprb,  32.6999435_jprb,  32.6727371_jprb,  32.6455307_jprb,  32.6183243_jprb,  &
  32.5931396_jprb,  32.5685806_jprb,  32.5440216_jprb,  32.5194664_jprb,  32.4952126_jprb,  &
  32.4709587_jprb,  32.4467049_jprb,  32.4211998_jprb,  32.3919983_jprb,  32.3627968_jprb,  &
  32.3335953_jprb,  32.3058357_jprb,  32.2796783_jprb,  32.2535172_jprb,  32.2273598_jprb,  &
  32.2004166_jprb,  32.1731071_jprb,  32.1458015_jprb,  32.1184921_jprb,  32.0921860_jprb,  &
  32.0660248_jprb,  32.0398598_jprb,  32.0139084_jprb,  31.9917793_jprb,  31.9696484_jprb,  &
  31.9475174_jprb,  31.9254494_jprb,  31.9036102_jprb,  31.8817711_jprb,  31.8599319_jprb,  &
  31.8373051_jprb,  31.8132954_jprb,  31.7892838_jprb,  31.7652740_jprb,  31.7417164_jprb,  &
  31.7186184_jprb,  31.6955204_jprb,  31.6724224_jprb,  31.6495037_jprb,  31.6266994_jprb,  &
  31.6038971_jprb,  31.5810928_jprb,  31.5579185_jprb,  31.5345955_jprb,  31.5112743_jprb,  &
  31.4879513_jprb,  31.4651890_jprb,  31.4425659_jprb,  31.4199409_jprb,  31.3973179_jprb,  &
  31.3761768_jprb,  31.3552513_jprb,  31.3343277_jprb,  31.3134022_jprb,  31.2921829_jprb,  &
  31.2709408_jprb,  31.2496967_jprb,  31.2284546_jprb,  31.2088394_jprb,  31.1892796_jprb,  &
  31.1697178_jprb,  31.1501560_jprb,  31.1285019_jprb,  31.1068325_jprb,  31.0851631_jprb,  &
  31.0634918_jprb,  31.0434933_jprb,  31.0234928_jprb,  31.0034924_jprb,  30.9834938_jprb,  &
  30.9657249_jprb,  30.9479752_jprb,  30.9302254_jprb,  30.9124756_jprb,  30.8942242_jprb,  &
  30.8759556_jprb,  30.8576870_jprb,  30.8394184_jprb,  30.8213997_jprb,  30.8034000_jprb,  &
  30.7854004_jprb,  30.7674026_jprb,  30.7479382_jprb,  30.7282581_jprb,  30.7085762_jprb,  &
  30.6888962_jprb,  30.6694336_jprb,  30.6500263_jprb,  30.6306190_jprb,  30.6112118_jprb,  &
  30.5930939_jprb,  30.5755043_jprb,  30.5579147_jprb,  30.5403271_jprb,  30.5217361_jprb,  &
  30.5024853_jprb,  30.4832344_jprb,  30.4639854_jprb,  30.4447708_jprb,  30.4256001_jprb,  &
  30.4064274_jprb,  30.3872547_jprb,  30.3688240_jprb,  30.3518295_jprb,  30.3348351_jprb,  &
  30.3178406_jprb,  30.3011398_jprb,  30.2857361_jprb,  30.2703323_jprb,  30.2549286_jprb,  &
  30.2395267_jprb,  30.2243767_jprb,  30.2092247_jprb,  30.1940746_jprb,  30.1789227_jprb,  &
  30.1636677_jprb,  30.1483898_jprb,  30.1331120_jprb,  30.1178341_jprb,  30.1022491_jprb,  &
  30.0864716_jprb,  30.0706921_jprb,  30.0549126_jprb,  30.0396652_jprb,  30.0252609_jprb,  &
  30.0108566_jprb,  29.9964523_jprb,  29.9821110_jprb,  29.9681377_jprb,  29.9541645_jprb,  &
  29.9401932_jprb,  29.9262199_jprb,  29.9103298_jprb,  29.8941879_jprb,  29.8780460_jprb,  &
  29.8619061_jprb,  29.8450203_jprb,  29.8276520_jprb,  29.8102837_jprb,  29.7929153_jprb,  &
  29.7757988_jprb,  29.7592392_jprb,  29.7426777_jprb,  29.7261181_jprb,  29.7095585_jprb,  &
  29.6996498_jprb,  29.6897545_jprb,  29.6798592_jprb,  29.6699657_jprb,  29.6582737_jprb,  &
  29.6456814_jprb,  29.6330910_jprb,  29.6205006_jprb,  29.6072922_jprb,  29.5927544_jprb,  &
  29.5782166_jprb,  29.5636806_jprb,  29.5491428_jprb,  29.5346622_jprb,  29.5201855_jprb,  &
  29.5057106_jprb,  29.4912338_jprb,  29.4785175_jprb,  29.4671650_jprb,  29.4558125_jprb,  &
  29.4444599_jprb,  29.4330292_jprb,  29.4211864_jprb,  29.4093437_jprb,  29.3974991_jprb,  &
  29.3856564_jprb,  29.3738499_jprb,  29.3620567_jprb,  29.3502617_jprb,  29.3384686_jprb,  &
  29.3269520_jprb,  29.3160992_jprb,  29.3052444_jprb,  29.2943897_jprb,  29.2835350_jprb,  &
  29.2743473_jprb,  29.2654877_jprb,  29.2566280_jprb,  29.2477703_jprb,  29.2384205_jprb,  &
  29.2281876_jprb,  29.2179546_jprb,  29.2077217_jprb,  29.1974869_jprb,  29.1856327_jprb,  &
  29.1735115_jprb,  29.1613884_jprb,  29.1492653_jprb,  29.1374588_jprb,  29.1262627_jprb,  &
  29.1150646_jprb,  29.1038685_jprb,  29.0926723_jprb,  29.0812340_jprb,  29.0697365_jprb,  &
  29.0582390_jprb,  29.0467434_jprb,  29.0351295_jprb,  29.0231609_jprb,  29.0111923_jprb,  &
  28.9992237_jprb,  28.9872551_jprb,  28.9762497_jprb,  28.9657135_jprb,  28.9551773_jprb,  &
  28.9446430_jprb,  28.9341106_jprb,  28.9236202_jprb,  28.9131298_jprb,  28.9026394_jprb,  &
  28.8921471_jprb,  28.8813610_jprb,  28.8702316_jprb,  28.8591003_jprb,  28.8479710_jprb,  &
  28.8368397_jprb,  28.8234921_jprb,  28.8096809_jprb,  28.7958698_jprb,  28.7820587_jprb,  &
  28.7683182_jprb,  28.7549057_jprb,  28.7414951_jprb,  28.7280827_jprb,  28.7146721_jprb,  &
  28.7011204_jprb,  28.6874294_jprb,  28.6737366_jprb,  28.6600456_jprb,  28.6463528_jprb,  &
  28.6348705_jprb,  28.6239319_jprb,  28.6129913_jprb,  28.6020527_jprb,  28.5911160_jprb,  &
  28.5802250_jprb,  28.5693321_jprb,  28.5584393_jprb,  28.5475483_jprb,  28.5352097_jprb,  &
  28.5201931_jprb,  28.5051765_jprb,  28.4901619_jprb,  28.4751453_jprb,  28.4604645_jprb,  &
  28.4460125_jprb,  28.4315605_jprb,  28.4171085_jprb,  28.4026566_jprb,  28.3889313_jprb,  &
  28.3753681_jprb,  28.3618069_jprb,  28.3482456_jprb,  28.3346806_jprb,  28.3210144_jprb,  &
  28.3073463_jprb,  28.2936783_jprb,  28.2800121_jprb,  28.2669716_jprb,  28.2564793_jprb,  &
  28.2459869_jprb,  28.2354946_jprb,  28.2250023_jprb,  28.2137108_jprb,  28.2009792_jprb,  &
  28.1882477_jprb,  28.1755161_jprb,  28.1627846_jprb,  28.1492786_jprb,  28.1349773_jprb,  &
  28.1206779_jprb,  28.1063766_jprb,  28.0920753_jprb,  28.0768261_jprb,  28.0609722_jprb,  &
  28.0451164_jprb,  28.0292606_jprb,  28.0134068_jprb,  27.9966908_jprb,  27.9796181_jprb,  &
  27.9625454_jprb,  27.9454708_jprb,  27.9283981_jprb,  27.9122562_jprb,  27.8963795_jprb,  &
  27.8805046_jprb,  27.8646297_jprb,  27.8487549_jprb,  27.8316154_jprb,  27.8142128_jprb,  &
  27.7968121_jprb,  27.7794113_jprb,  27.7620106_jprb,  27.7442665_jprb,  27.7264652_jprb,  &
  27.7086658_jprb,  27.6908646_jprb,  27.6730633_jprb,  27.6553326_jprb,  27.6376114_jprb,  &
  27.6198902_jprb,  27.6021690_jprb,  27.5844479_jprb,  27.5658588_jprb,  27.5471153_jprb,  &
  27.5283699_jprb,  27.5096264_jprb,  27.4908829_jprb,  27.4715710_jprb,  27.4521275_jprb,  &
  27.4326820_jprb,  27.4132385_jprb,  27.3937931_jprb,  27.3746490_jprb,  27.3556042_jprb,  &
  27.3365593_jprb,  27.3175144_jprb,  27.2984695_jprb,  27.2805119_jprb,  27.2631073_jprb,  &
  27.2457027_jprb,  27.2282982_jprb,  27.2108917_jprb,  27.1937027_jprb,  27.1766872_jprb,  &
  27.1596699_jprb,  27.1426525_jprb,  27.1256351_jprb,  27.1083889_jprb,  27.0908337_jprb,  &
  27.0732765_jprb,  27.0557213_jprb,  27.0381641_jprb,  27.0203800_jprb,  27.0019836_jprb,  &
  26.9835873_jprb,  26.9651909_jprb,  26.9467945_jprb,  26.9284649_jprb,  26.9107628_jprb,  &
  26.8930607_jprb,  26.8753586_jprb,  26.8576565_jprb,  26.8399544_jprb,  26.8213730_jprb,  &
  26.8026886_jprb,  26.7840023_jprb,  26.7653179_jprb,  26.7466335_jprb,  26.7292213_jprb,  &
  26.7124367_jprb,  26.6956520_jprb,  26.6788673_jprb,  26.6620827_jprb,  26.6450748_jprb,  &
  26.6277657_jprb,  26.6104546_jprb,  26.5931454_jprb,  26.5758343_jprb,  26.5585136_jprb,  &
  26.5411339_jprb,  26.5237522_jprb,  26.5063725_jprb,  26.4889908_jprb,  26.4716110_jprb,  &
  26.4569683_jprb,  26.4428005_jprb,  26.4286308_jprb,  26.4144630_jprb,  26.4002934_jprb,  &
  26.3871841_jprb,  26.3750114_jprb,  26.3628368_jprb,  26.3506641_jprb,  26.3384914_jprb,  &
  26.3263550_jprb,  26.3143845_jprb,  26.3024158_jprb,  26.2904472_jprb,  26.2784767_jprb,  &
  26.2665081_jprb,  26.2512283_jprb,  26.2351971_jprb,  26.2191639_jprb,  26.2031307_jprb,  &
  26.1870995_jprb,  26.1714039_jprb,  26.1561775_jprb,  26.1409512_jprb,  26.1257267_jprb,  &
  26.1105003_jprb,  26.0952740_jprb,  26.0812817_jprb,  26.0672913_jprb,  26.0533028_jprb,  &
  26.0393124_jprb,  26.0253239_jprb,  26.0126495_jprb,  26.0010452_jprb,  25.9894428_jprb,  &
  25.9778385_jprb,  25.9662361_jprb,  25.9545670_jprb,  25.9421520_jprb,  25.9297352_jprb,  &
  25.9173203_jprb,  25.9049053_jprb,  25.8924885_jprb,  25.8806934_jprb,  25.8693409_jprb,  &
  25.8579903_jprb,  25.8466377_jprb,  25.8352871_jprb,  25.8239040_jprb,  25.8120346_jprb,  &
  25.8001652_jprb,  25.7882938_jprb,  25.7764244_jprb,  25.7645531_jprb,  25.7520580_jprb,  &
  25.7389641_jprb,  25.7258682_jprb,  25.7127724_jprb,  25.6996765_jprb,  25.6865807_jprb,  &
  25.6744690_jprb,  25.6624279_jprb,  25.6503849_jprb,  25.6383438_jprb,  25.6263008_jprb,  &
  25.6137161_jprb,  25.6000385_jprb,  25.5863628_jprb,  25.5726852_jprb,  25.5590076_jprb,  &
  25.5453300_jprb,  25.5319958_jprb,  25.5188046_jprb,  25.5056133_jprb,  25.4924221_jprb,  &
  25.4792290_jprb,  25.4660549_jprb,  25.4532051_jprb,  25.4403553_jprb,  25.4275074_jprb,  &
  25.4146576_jprb,  25.4018078_jprb,  25.3888283_jprb,  25.3756237_jprb,  25.3624191_jprb,  &
  25.3492126_jprb,  25.3360081_jprb,  25.3228035_jprb,  25.3096390_jprb,  25.2964973_jprb,  &
  25.2833557_jprb,  25.2702122_jprb,  25.2570705_jprb,  25.2439289_jprb,  25.2310963_jprb,  &
  25.2182941_jprb,  25.2054901_jprb,  25.1926861_jprb,  25.1798820_jprb,  25.1671085_jprb,  &
  25.1545029_jprb,  25.1418972_jprb,  25.1292915_jprb,  25.1166859_jprb,  25.1040802_jprb,  &
  25.0913010_jprb,  25.0782108_jprb,  25.0651207_jprb,  25.0520287_jprb,  25.0389385_jprb,  &
  25.0258484_jprb,  25.0122795_jprb,  24.9983025_jprb,  24.9843254_jprb,  24.9703484_jprb,  &
  24.9563713_jprb,  24.9423923_jprb,  24.9293041_jprb,  24.9166107_jprb,  24.9039154_jprb,  &
  24.8912220_jprb,  24.8785286_jprb,  24.8658333_jprb,  24.8535175_jprb,  24.8412895_jprb,  &
  24.8290596_jprb,  24.8168297_jprb,  24.8045998_jprb,  24.7923717_jprb,  24.7742329_jprb,  &
  24.7555084_jprb,  24.7367859_jprb,  24.7180614_jprb,  24.6993389_jprb,  24.6806145_jprb,  &
  24.6643181_jprb,  24.6480808_jprb,  24.6318455_jprb,  24.6156082_jprb,  24.5993710_jprb,  &
  24.5831604_jprb,  24.5689888_jprb,  24.5548191_jprb,  24.5406475_jprb,  24.5264759_jprb,  &
  24.5123062_jprb,  24.4981022_jprb,  24.4822884_jprb,  24.4664726_jprb,  24.4506588_jprb,  &
  24.4348431_jprb,  24.4190292_jprb,  24.4032135_jprb,  24.3889141_jprb,  24.3746204_jprb,  &
  24.3603249_jprb,  24.3460312_jprb,  24.3317356_jprb,  24.3174419_jprb,  24.3038292_jprb,  &
  24.2902584_jprb,  24.2766857_jprb,  24.2631149_jprb,  24.2495422_jprb,  24.2359695_jprb,  &
  24.2227917_jprb,  24.2096767_jprb,  24.1965618_jprb,  24.1834450_jprb,  24.1703300_jprb,  &
  24.1572151_jprb,  24.1433697_jprb,  24.1292858_jprb,  24.1152000_jprb,  24.1011162_jprb,  &
  24.0870323_jprb,  24.0729465_jprb,  24.0578842_jprb,  24.0421982_jprb,  24.0265121_jprb,  &
  24.0108261_jprb,  23.9951401_jprb,  23.9794540_jprb,  23.9644184_jprb,  23.9502163_jprb,  &
  23.9360161_jprb,  23.9218140_jprb,  23.9076138_jprb,  23.8934116_jprb,  23.8783321_jprb,  &
  23.8603821_jprb,  23.8424339_jprb,  23.8244839_jprb,  23.8065338_jprb,  23.7885838_jprb,  &
  23.7706356_jprb /

  DATA (pcm_snow(i), i=1, 1000) / &
    4.0553737_jprb,   4.1246114_jprb,   4.1928368_jprb,   4.2636189_jprb,   4.3337398_jprb,  &
    4.4063163_jprb,   4.4800682_jprb,   4.5561152_jprb,   4.6343346_jprb,   4.7189765_jprb,  &
    4.7996798_jprb,   4.8793116_jprb,   4.9580832_jprb,   5.0318923_jprb,   5.1114373_jprb,  &
    5.1895919_jprb,   5.2723780_jprb,   5.3648539_jprb,   5.4638696_jprb,   5.5692205_jprb,  &
    5.6860938_jprb,   5.8161249_jprb,   5.9647784_jprb,   6.1386642_jprb,   6.3357239_jprb,  &
    6.5664897_jprb,   6.8242984_jprb,   7.1092768_jprb,   7.4305677_jprb,   7.7802258_jprb,  &
    8.1575098_jprb,   8.5571251_jprb,   8.9828501_jprb,   9.4280252_jprb,   9.8945770_jprb,  &
   10.3629398_jprb,  10.8178520_jprb,  11.2583809_jprb,  11.6762304_jprb,  12.0569487_jprb,  &
   12.3957186_jprb,  12.6877918_jprb,  12.9178181_jprb,  13.0827341_jprb,  13.1847124_jprb,  &
   13.2117481_jprb,  13.1602411_jprb,  13.0392017_jprb,  12.8595781_jprb,  12.6281862_jprb,  &
   12.3657923_jprb,  12.0876589_jprb,  11.7963610_jprb,  11.4928522_jprb,  11.1861620_jprb,  &
   10.8751917_jprb,  10.5504141_jprb,  10.2162085_jprb,   9.8714495_jprb,   9.5135803_jprb,  &
    9.1376371_jprb,   8.7611217_jprb,   8.3838892_jprb,   8.0066757_jprb,   7.6149139_jprb,  &
    7.2238712_jprb,   6.8367028_jprb,   6.4565711_jprb,   6.0764236_jprb,   5.7040277_jprb,  &
    5.3444228_jprb,   4.9989600_jprb,   4.6669912_jprb,   4.3473830_jprb,   4.0495453_jprb,  &
    3.7781017_jprb,   3.5276525_jprb,   3.2919629_jprb,   3.0769436_jprb,   2.8832572_jprb,  &
    2.7075057_jprb,   2.5467250_jprb,   2.4020565_jprb,   2.2785034_jprb,   2.1708882_jprb,  &
    2.0732737_jprb,   1.9852085_jprb,   1.9068877_jprb,   1.8385501_jprb,   1.7804509_jprb,  &
    1.7317245_jprb,   1.6880711_jprb,   1.6488868_jprb,   1.6147453_jprb,   1.5856470_jprb,  &
    1.5590069_jprb,   1.5332441_jprb,   1.5093766_jprb,   1.4903805_jprb,   1.4769915_jprb,  &
    1.4673688_jprb,   1.4614617_jprb,   1.4587358_jprb,   1.4627540_jprb,   1.4743274_jprb,  &
    1.4915078_jprb,   1.5155025_jprb,   1.5466309_jprb,   1.5892223_jprb,   1.6462225_jprb,  &
    1.7172153_jprb,   1.8015628_jprb,   1.9023141_jprb,   2.0237238_jprb,   2.1780450_jprb,  &
    2.3626380_jprb,   2.5735414_jprb,   2.8147416_jprb,   3.1017640_jprb,   3.4467602_jprb,  &
    3.8346391_jprb,   4.2658944_jprb,   4.7574492_jprb,   5.3149104_jprb,   5.9249001_jprb,  &
    6.5867810_jprb,   7.3210101_jprb,   8.0712986_jprb,   8.8337259_jprb,   9.6074133_jprb,  &
   10.3865108_jprb,  11.1650133_jprb,  11.9357920_jprb,  12.6365261_jprb,  13.3111420_jprb,  &
   13.9507036_jprb,  14.5209055_jprb,  15.0041761_jprb,  15.4113321_jprb,  15.7224703_jprb,  &
   15.9541159_jprb,  16.1331406_jprb,  16.2587776_jprb,  16.3334980_jprb,  16.3743229_jprb,  &
   16.3829975_jprb,  16.3609028_jprb,  16.3169346_jprb,  16.2537804_jprb,  16.1662960_jprb,  &
   16.0674839_jprb,  15.9593048_jprb,  15.8483410_jprb,  15.7382011_jprb,  15.6288462_jprb,  &
   15.5307112_jprb,  15.4354954_jprb,  15.3437147_jprb,  15.2536201_jprb,  15.1666365_jprb,  &
   15.0853271_jprb,  15.0052948_jprb,  14.9247894_jprb,  14.8409700_jprb,  14.7579441_jprb,  &
   14.6753168_jprb,  14.5896807_jprb,  14.5002909_jprb,  14.4058676_jprb,  14.3099689_jprb,  &
   14.2128115_jprb,  14.1119270_jprb,  14.0089111_jprb,  13.9043102_jprb,  13.7987461_jprb,  &
   13.6917839_jprb,  13.5815182_jprb,  13.4652328_jprb,  13.3454618_jprb,  13.2192478_jprb,  &
   13.0904217_jprb,  12.9557676_jprb,  12.8172388_jprb,  12.6761875_jprb,  12.5346012_jprb,  &
   12.3934193_jprb,  12.2540188_jprb,  12.1158447_jprb,  11.9789352_jprb,  11.8428793_jprb,  &
   11.7071428_jprb,  11.5731201_jprb,  11.4398155_jprb,  11.3082771_jprb,  11.1756840_jprb,  &
   11.0421667_jprb,  10.9088335_jprb,  10.7755632_jprb,  10.6446409_jprb,  10.5136204_jprb,  &
   10.3820734_jprb,  10.2512646_jprb,  10.1215572_jprb,   9.9978018_jprb,   9.8782005_jprb,  &
    9.7669287_jprb,   9.6581450_jprb,   9.5563879_jprb,   9.4560547_jprb,   9.3646936_jprb,  &
    9.2752285_jprb,   9.1915960_jprb,   9.1095400_jprb,   9.0300608_jprb,   8.9497738_jprb,  &
    8.8686943_jprb,   8.7854099_jprb,   8.7007570_jprb,   8.6093712_jprb,   8.5153236_jprb,  &
    8.4138937_jprb,   8.3106766_jprb,   8.2011461_jprb,   8.0907774_jprb,   7.9711218_jprb,  &
    7.8509483_jprb,   7.7281280_jprb,   7.6049399_jprb,   7.4730115_jprb,   7.3407402_jprb,  &
    7.2039762_jprb,   7.0665069_jprb,   6.9215403_jprb,   6.7764301_jprb,   6.6299238_jprb,  &
    6.4834390_jprb,   6.3372664_jprb,   6.1918778_jprb,   6.0593410_jprb,   5.9270425_jprb,  &
    5.8052850_jprb,   5.6835237_jprb,   5.5714068_jprb,   5.4595499_jprb,   5.3591290_jprb,  &
    5.2598019_jprb,   5.1625376_jprb,   5.0656695_jprb,   4.9755502_jprb,   4.8876719_jprb,  &
    4.8039341_jprb,   4.7224364_jprb,   4.6461449_jprb,   4.5743880_jprb,   4.5054483_jprb,  &
    4.4406414_jprb,   4.3787889_jprb,   4.3252454_jprb,   4.2719593_jprb,   4.2208495_jprb,  &
    4.1697426_jprb,   4.1248293_jprb,   4.0803547_jprb,   4.0390439_jprb,   3.9988046_jprb,  &
    3.9672081_jprb,   3.9427361_jprb,   3.9270284_jprb,   3.9288402_jprb,   3.9338572_jprb,  &
    3.9666328_jprb,   3.9994082_jprb,   4.0758743_jprb,   4.1595087_jprb,   4.2763886_jprb,  &
    4.4154940_jprb,   4.5744205_jprb,   4.7746387_jprb,   4.9768190_jprb,   5.2348313_jprb,  &
    5.4928441_jprb,   5.8082228_jprb,   6.1450505_jprb,   6.5089293_jprb,   6.9123740_jprb,  &
    7.3205743_jprb,   7.7934623_jprb,   8.2663517_jprb,   8.7814236_jprb,   9.3133192_jprb,  &
    9.8769264_jprb,  10.5007563_jprb,  11.1245985_jprb,  11.8259106_jprb,  12.5305901_jprb,  &
   13.3271646_jprb,  14.1972742_jprb,  15.0918217_jprb,  16.1419392_jprb,  17.1920719_jprb,  &
   18.3762035_jprb,  19.6179638_jprb,  20.8995876_jprb,  22.3036671_jprb,  23.7077599_jprb,  &
   25.2310696_jprb,  26.7891235_jprb,  28.3656197_jprb,  29.9883385_jprb,  31.6110497_jprb,  &
   33.0576324_jprb,  34.4545746_jprb,  35.7899628_jprb,  36.9466515_jprb,  38.1033592_jprb,  &
   39.0609474_jprb,  39.9389038_jprb,  40.7734146_jprb,  41.3697205_jprb,  41.9660301_jprb,  &
   42.3929214_jprb,  42.6948853_jprb,  42.9968529_jprb,  43.2506599_jprb,  43.5034218_jprb,  &
   43.7357178_jprb,  43.9320259_jprb,  44.1282349_jprb,  44.2968140_jprb,  44.4549942_jprb,  &
   44.6106644_jprb,  44.7339745_jprb,  44.8572998_jprb,  44.9692268_jprb,  45.0640945_jprb,  &
   45.1589546_jprb,  45.2526360_jprb,  45.3458481_jprb,  45.4390488_jprb,  45.5309601_jprb,  &
   45.6228790_jprb,  45.7109528_jprb,  45.7881355_jprb,  45.8653107_jprb,  45.9374924_jprb,  &
   46.0048065_jprb,  46.0721474_jprb,  46.1243095_jprb,  46.1709137_jprb,  46.2175026_jprb,  &
   46.2389030_jprb,  46.2585907_jprb,  46.2762337_jprb,  46.2787933_jprb,  46.2813301_jprb,  &
   46.2840157_jprb,  46.2870560_jprb,  46.2900772_jprb,  46.2862282_jprb,  46.2730827_jprb,  &
   46.2599297_jprb,  46.2464981_jprb,  46.2329025_jprb,  46.2192917_jprb,  46.2048759_jprb,  &
   46.1900978_jprb,  46.1753464_jprb,  46.1570129_jprb,  46.1372871_jprb,  46.1175499_jprb,  &
   46.0926819_jprb,  46.0662804_jprb,  46.0398712_jprb,  46.0122604_jprb,  45.9842682_jprb,  &
   45.9563255_jprb,  45.9261971_jprb,  45.8955345_jprb,  45.8648911_jprb,  45.8247223_jprb,  &
   45.7820396_jprb,  45.7393150_jprb,  45.7026825_jprb,  45.6680565_jprb,  45.6334190_jprb,  &
   45.6128998_jprb,  45.5988083_jprb,  45.5847473_jprb,  45.5800705_jprb,  45.5816345_jprb,  &
   45.5831871_jprb,  45.5929413_jprb,  45.6113014_jprb,  45.6296692_jprb,  45.6625137_jprb,  &
   45.7224007_jprb,  45.7822914_jprb,  45.8454170_jprb,  45.9225960_jprb,  45.9997559_jprb,  &
   46.0769157_jprb,  46.1710892_jprb,  46.2652397_jprb,  46.3593407_jprb,  46.4769173_jprb,  &
   46.6008644_jprb,  46.7247314_jprb,  46.8686371_jprb,  47.0290909_jprb,  47.1895103_jprb,  &
   47.3553658_jprb,  47.5350533_jprb,  47.7147560_jprb,  47.8944702_jprb,  48.1136856_jprb,  &
   48.3332062_jprb,  48.5527573_jprb,  48.8041611_jprb,  49.0710106_jprb,  49.3378868_jprb,  &
   49.6293259_jprb,  49.9710693_jprb,  50.3128357_jprb,  50.6546364_jprb,  51.0576210_jprb,  &
   51.4635010_jprb,  51.8693962_jprb,  52.2835007_jprb,  52.7042732_jprb,  53.1250267_jprb,  &
   53.5450249_jprb,  53.9595299_jprb,  54.3739891_jprb,  54.7884636_jprb,  55.1968155_jprb,  &
   55.6021461_jprb,  56.0074348_jprb,  56.4066925_jprb,  56.7784119_jprb,  57.1500473_jprb,  &
   57.5217438_jprb,  57.8758430_jprb,  58.2210732_jprb,  58.5662804_jprb,  58.9133949_jprb,  &
   59.2743225_jprb,  59.6353149_jprb,  59.9962349_jprb,  60.3680267_jprb,  60.7488861_jprb,  &
   61.1297340_jprb,  61.5106163_jprb,  61.9026833_jprb,  62.2955627_jprb,  62.6883430_jprb,  &
   63.0848885_jprb,  63.4901581_jprb,  63.8953972_jprb,  64.3006287_jprb,  64.7092133_jprb,  &
   65.1198273_jprb,  65.5303726_jprb,  65.9409790_jprb,  66.2812576_jprb,  66.6163712_jprb,  &
   66.9514542_jprb,  67.2714767_jprb,  67.5308609_jprb,  67.7902451_jprb,  68.0496292_jprb,  &
   68.2836075_jprb,  68.4844742_jprb,  68.6854172_jprb,  68.8863449_jprb,  69.0521469_jprb,  &
   69.1980209_jprb,  69.3439102_jprb,  69.4898376_jprb,  69.5849686_jprb,  69.6678772_jprb,  &
   69.7508545_jprb,  69.8337936_jprb,  69.9457245_jprb,  70.0594254_jprb,  70.1730804_jprb,  &
   70.2869720_jprb,  70.4062042_jprb,  70.5254135_jprb,  70.6446381_jprb,  70.7624969_jprb,  &
   70.8699646_jprb,  70.9774857_jprb,  71.0849991_jprb,  71.1882095_jprb,  71.2658463_jprb,  &
   71.3434601_jprb,  71.4211121_jprb,  71.4972916_jprb,  71.5645370_jprb,  71.6318588_jprb,  &
   71.6991272_jprb,  71.7625046_jprb,  71.7904663_jprb,  71.8184128_jprb,  71.8463745_jprb,  &
   71.8739548_jprb,  71.8849564_jprb,  71.8959656_jprb,  71.9069443_jprb,  71.9179993_jprb,  &
   71.9147415_jprb,  71.9100113_jprb,  71.9052963_jprb,  71.9006348_jprb,  71.8937454_jprb,  &
   71.8862000_jprb,  71.8787155_jprb,  71.8711472_jprb,  71.8516769_jprb,  71.8232269_jprb,  &
   71.7946777_jprb,  71.7661743_jprb,  71.7344437_jprb,  71.6963654_jprb,  71.6582108_jprb,  &
   71.6200790_jprb,  71.5808792_jprb,  71.5286636_jprb,  71.4764557_jprb,  71.4241791_jprb,  &
   71.3719559_jprb,  71.3068314_jprb,  71.2377625_jprb,  71.1687393_jprb,  71.0996628_jprb,  &
   71.0190201_jprb,  70.9224243_jprb,  70.8258591_jprb,  70.7291946_jprb,  70.6319809_jprb,  &
   70.5155182_jprb,  70.3990860_jprb,  70.2825775_jprb,  70.1660919_jprb,  70.0471268_jprb,  &
   69.9265747_jprb,  69.8059311_jprb,  69.6853561_jprb,  69.5649261_jprb,  69.4454193_jprb,  &
   69.3259659_jprb,  69.2064819_jprb,  69.0870361_jprb,  68.9788666_jprb,  68.8777313_jprb,  &
   68.7765121_jprb,  68.6753311_jprb,  68.5757523_jprb,  68.5003510_jprb,  68.4248657_jprb,  &
   68.3494568_jprb,  68.2740402_jprb,  68.2172852_jprb,  68.1823196_jprb,  68.1474152_jprb,  &
   68.1124191_jprb,  68.0774918_jprb,  68.0608749_jprb,  68.0483780_jprb,  68.0358887_jprb,  &
   68.0233536_jprb,  68.0169373_jprb,  68.0525742_jprb,  68.0882111_jprb,  68.1238251_jprb,  &
   68.1594849_jprb,  68.2064667_jprb,  68.2710724_jprb,  68.3357010_jprb,  68.4003754_jprb,  &
   68.4649582_jprb,  68.5593338_jprb,  68.6727905_jprb,  68.7863159_jprb,  68.8998108_jprb,  &
   69.0133286_jprb,  69.1550064_jprb,  69.3047562_jprb,  69.4545670_jprb,  69.6043396_jprb,  &
   69.7541351_jprb,  69.9359207_jprb,  70.1212616_jprb,  70.3065567_jprb,  70.4918518_jprb,  &
   70.6772003_jprb,  70.8711700_jprb,  71.0653915_jprb,  71.2596054_jprb,  71.4538651_jprb,  &
   71.6480942_jprb,  71.8503494_jprb,  72.0526047_jprb,  72.2548294_jprb,  72.4571075_jprb,  &
   72.6594009_jprb,  72.8439865_jprb,  73.0280838_jprb,  73.2122192_jprb,  73.3963394_jprb,  &
   73.5804596_jprb,  73.7547989_jprb,  73.9280090_jprb,  74.1012497_jprb,  74.2745361_jprb,  &
   74.4477844_jprb,  74.6139755_jprb,  74.7782288_jprb,  74.9423676_jprb,  75.1065292_jprb,  &
   75.2706680_jprb,  75.4324722_jprb,  75.5926361_jprb,  75.7528229_jprb,  75.9130096_jprb,  &
   76.0731659_jprb,  76.2323685_jprb,  76.3897858_jprb,  76.5472488_jprb,  76.7047043_jprb,  &
   76.8621674_jprb,  77.0198212_jprb,  77.1802063_jprb,  77.3406372_jprb,  77.5010910_jprb,  &
   77.6615143_jprb,  77.8219376_jprb,  77.9791870_jprb,  78.1352539_jprb,  78.2913284_jprb,  &
   78.4473801_jprb,  78.6034622_jprb,  78.7522812_jprb,  78.8864288_jprb,  79.0206223_jprb,  &
   79.1548157_jprb,  79.2889557_jprb,  79.4231262_jprb,  79.5382690_jprb,  79.6508408_jprb,  &
   79.7633972_jprb,  79.8759689_jprb,  79.9884491_jprb,  80.0916138_jprb,  80.1785507_jprb,  &
   80.2654190_jprb,  80.3524017_jprb,  80.4393158_jprb,  80.5262604_jprb,  80.5989761_jprb,  &
   80.6681671_jprb,  80.7373657_jprb,  80.8065720_jprb,  80.8757401_jprb,  80.9411316_jprb,  &
   80.9877930_jprb,  81.0344467_jprb,  81.0811081_jprb,  81.1277237_jprb,  81.1743622_jprb,  &
   81.2213593_jprb,  81.2685699_jprb,  81.3158112_jprb,  81.3630295_jprb,  81.4102631_jprb,  &
   81.4575500_jprb,  81.4991837_jprb,  81.5387726_jprb,  81.5783615_jprb,  81.6179581_jprb,  &
   81.6575470_jprb,  81.6971359_jprb,  81.7357712_jprb,  81.7741928_jprb,  81.8127136_jprb,  &
   81.8512344_jprb,  81.8897400_jprb,  81.9278870_jprb,  81.9594727_jprb,  81.9909897_jprb,  &
   82.0225449_jprb,  82.0540314_jprb,  82.0855789_jprb,  82.1171112_jprb,  82.1484833_jprb,  &
   82.1798477_jprb,  82.2112961_jprb,  82.2427063_jprb,  82.2740479_jprb,  82.3052139_jprb,  &
   82.3340073_jprb,  82.3628769_jprb,  82.3917389_jprb,  82.4205856_jprb,  82.4493484_jprb,  &
   82.4782257_jprb,  82.5075684_jprb,  82.5369644_jprb,  82.5663605_jprb,  82.5956879_jprb,  &
   82.6250458_jprb,  82.6544266_jprb,  82.6924973_jprb,  82.7315140_jprb,  82.7706299_jprb,  &
   82.8097229_jprb,  82.8487930_jprb,  82.8878860_jprb,  82.9339218_jprb,  82.9832077_jprb,  &
   83.0324478_jprb,  83.0817032_jprb,  83.1309967_jprb,  83.1802368_jprb,  83.2350311_jprb,  &
   83.2979584_jprb,  83.3609085_jprb,  83.4238205_jprb,  83.4867249_jprb,  83.5495987_jprb,  &
   83.6132660_jprb,  83.6882324_jprb,  83.7631989_jprb,  83.8380737_jprb,  83.9130478_jprb,  &
   83.9879608_jprb,  84.0629425_jprb,  84.1447754_jprb,  84.2305450_jprb,  84.3162079_jprb,  &
   84.4020157_jprb,  84.4877930_jprb,  84.5734558_jprb,  84.6593552_jprb,  84.7461472_jprb,  &
   84.8328857_jprb,  84.9196548_jprb,  85.0064163_jprb,  85.0931854_jprb,  85.1799316_jprb,  &
   85.2636566_jprb,  85.3450470_jprb,  85.4263840_jprb,  85.5078049_jprb,  85.5892029_jprb,  &
   85.6706314_jprb,  85.7520447_jprb,  85.8229599_jprb,  85.8929901_jprb,  85.9631271_jprb,  &
   86.0332336_jprb,  86.1033630_jprb,  86.1734543_jprb,  86.2419739_jprb,  86.3045654_jprb,  &
   86.3671875_jprb,  86.4298172_jprb,  86.4923782_jprb,  86.5550232_jprb,  86.6175232_jprb,  &
   86.6770248_jprb,  86.7317276_jprb,  86.7864914_jprb,  86.8412628_jprb,  86.8959885_jprb,  &
   86.9507523_jprb,  87.0056305_jprb,  87.0547256_jprb,  87.0987473_jprb,  87.1428757_jprb,  &
   87.1870270_jprb,  87.2311172_jprb,  87.2752380_jprb,  87.3193436_jprb,  87.3561859_jprb,  &
   87.3873978_jprb,  87.4185333_jprb,  87.4497910_jprb,  87.4809647_jprb,  87.5122147_jprb,  &
   87.5434265_jprb,  87.5711517_jprb,  87.5956497_jprb,  87.6201706_jprb,  87.6447830_jprb,  &
   87.6693420_jprb,  87.6938782_jprb,  87.7184677_jprb,  87.7383118_jprb,  87.7511520_jprb,  &
   87.7639084_jprb,  87.7767105_jprb,  87.7895126_jprb,  87.8022385_jprb,  87.8150635_jprb,  &
   87.8263550_jprb,  87.8309479_jprb,  87.8354568_jprb,  87.8400803_jprb,  87.8446426_jprb,  &
   87.8491898_jprb,  87.8537521_jprb,  87.8583527_jprb,  87.8612671_jprb,  87.8639984_jprb,  &
   87.8666916_jprb,  87.8694916_jprb,  87.8721695_jprb,  87.8749313_jprb,  87.8776245_jprb,  &
   87.8815384_jprb,  87.8867569_jprb,  87.8919983_jprb,  87.8971176_jprb,  87.9023590_jprb,  &
   87.9075470_jprb,  87.9127655_jprb,  87.9179764_jprb,  87.9217453_jprb,  87.9254150_jprb,  &
   87.9291229_jprb,  87.9327927_jprb,  87.9365234_jprb,  87.9402390_jprb,  87.9440155_jprb,  &
   87.9503479_jprb,  87.9601746_jprb,  87.9700089_jprb,  87.9797745_jprb,  87.9896088_jprb,  &
   87.9994354_jprb,  88.0092697_jprb,  88.0189972_jprb,  88.0332565_jprb,  88.0487823_jprb,  &
   88.0643463_jprb,  88.0800018_jprb,  88.0955505_jprb,  88.1110916_jprb,  88.1266861_jprb,  &
   88.1422272_jprb,  88.1619720_jprb,  88.1816483_jprb,  88.2013245_jprb,  88.2210770_jprb,  &
   88.2407379_jprb,  88.2604752_jprb,  88.2801743_jprb,  88.3005447_jprb,  88.3250198_jprb,  &
   88.3495941_jprb,  88.3739853_jprb,  88.3985214_jprb,  88.4230042_jprb,  88.4473724_jprb,  &
   88.4719315_jprb,  88.4972305_jprb,  88.5262299_jprb,  88.5551453_jprb,  88.5841141_jprb,  &
   88.6131058_jprb,  88.6421127_jprb,  88.6710358_jprb,  88.7000351_jprb,  88.7289810_jprb,  &
   88.7578278_jprb,  88.7866898_jprb,  88.8154755_jprb,  88.8443146_jprb,  88.8731384_jprb,  &
   88.9019775_jprb,  88.9308243_jprb,  88.9595947_jprb,  88.9865341_jprb,  89.0133209_jprb,  &
   89.0401001_jprb,  89.0668182_jprb,  89.0936813_jprb,  89.1203308_jprb,  89.1471405_jprb,  &
   89.1739426_jprb,  89.2010345_jprb,  89.2283630_jprb,  89.2557449_jprb,  89.2830200_jprb,  &
   89.3104019_jprb,  89.3376999_jprb,  89.3650284_jprb,  89.3923721_jprb,  89.4203033_jprb,  &
   89.4499741_jprb,  89.4795303_jprb,  89.5091324_jprb,  89.5386505_jprb,  89.5683060_jprb,  &
   89.5979080_jprb,  89.6274719_jprb,  89.6570206_jprb,  89.6859818_jprb,  89.7148743_jprb,  &
   89.7435989_jprb,  89.7724228_jprb,  89.8012772_jprb,  89.8301010_jprb,  89.8588409_jprb,  &
   89.8877182_jprb,  89.9160919_jprb,  89.9434509_jprb,  89.9706726_jprb,  89.9979324_jprb,  &
   90.0251389_jprb,  90.0523834_jprb,  90.0796890_jprb,  90.1068878_jprb,  90.1341019_jprb,  &
   90.1612167_jprb,  90.1880951_jprb,  90.2150574_jprb,  90.2418518_jprb,  90.2687225_jprb,  &
   90.2955780_jprb,  90.3224335_jprb,  90.3493347_jprb,  90.3761673_jprb,  90.3990555_jprb,  &
   90.4200516_jprb,  90.4409943_jprb,  90.4618683_jprb,  90.4828949_jprb,  90.5038223_jprb,  &
   90.5246506_jprb,  90.5456314_jprb,  90.5664978_jprb,  90.5854874_jprb,  90.6035309_jprb,  &
   90.6216278_jprb,  90.6396713_jprb,  90.6577377_jprb,  90.6758194_jprb,  90.6939316_jprb,  &
   90.7119446_jprb,  90.7300262_jprb,  90.7448502_jprb,  90.7578583_jprb,  90.7708282_jprb,  &
   90.7838745_jprb,  90.7968674_jprb,  90.8099747_jprb,  90.8229218_jprb,  90.8360062_jprb,  &
   90.8489914_jprb,  90.8605347_jprb,  90.8703842_jprb,  90.8801651_jprb,  90.8900223_jprb,  &
   90.8997955_jprb,  90.9096680_jprb,  90.9195709_jprb,  90.9292831_jprb,  90.9391479_jprb,  &
   90.9489365_jprb,  90.9585342_jprb,  90.9682770_jprb,  90.9779434_jprb,  90.9875565_jprb,  &
   90.9971924_jprb,  91.0068741_jprb,  91.0164948_jprb,  91.0260849_jprb,  91.0356827_jprb,  &
   91.0468674_jprb,  91.0588837_jprb,  91.0707474_jprb,  91.0827026_jprb,  91.0945206_jprb,  &
   91.1064453_jprb,  91.1183853_jprb,  91.1303482_jprb,  91.1422348_jprb,  91.1538696_jprb,  &
   91.1633759_jprb,  91.1728897_jprb,  91.1824265_jprb,  91.1919403_jprb,  91.2014313_jprb,  &
   91.2109909_jprb,  91.2204819_jprb,  91.2299881_jprb,  91.2394791_jprb,  91.2499695_jprb/

  DATA (pcm_snow(i), i=1001, 2101) / &
   91.2617645_jprb,  91.2736130_jprb,  91.2854080_jprb,  91.2972412_jprb,  91.3090820_jprb,  &
   91.3209076_jprb,  91.3327637_jprb,  91.3445511_jprb,  91.3564377_jprb,  91.3688507_jprb,  &
   91.3819656_jprb,  91.3950958_jprb,  91.4081726_jprb,  91.4212570_jprb,  91.4343033_jprb,  &
   91.4473953_jprb,  91.4604416_jprb,  91.4734955_jprb,  91.4866638_jprb,  91.4991531_jprb,  &
   91.5113602_jprb,  91.5235901_jprb,  91.5357132_jprb,  91.5478668_jprb,  91.5600586_jprb,  &
   91.5722122_jprb,  91.5843658_jprb,  91.5965195_jprb,  91.6087418_jprb,  91.6203842_jprb,  &
   91.6311264_jprb,  91.6420593_jprb,  91.6529388_jprb,  91.6638336_jprb,  91.6747208_jprb,  &
   91.6855316_jprb,  91.6964645_jprb,  91.7073593_jprb,  91.7181854_jprb,  91.7289886_jprb,  &
   91.7391357_jprb,  91.7492752_jprb,  91.7593689_jprb,  91.7694168_jprb,  91.7795868_jprb,  &
   91.7896729_jprb,  91.7997742_jprb,  91.8098145_jprb,  91.8199234_jprb,  91.8300018_jprb,  &
   91.8388672_jprb,  91.8470764_jprb,  91.8553543_jprb,  91.8636398_jprb,  91.8718872_jprb,  &
   91.8802109_jprb,  91.8884811_jprb,  91.8966675_jprb,  91.9050674_jprb,  91.9132919_jprb,  &
   91.9213257_jprb,  91.9279251_jprb,  91.9344482_jprb,  91.9410782_jprb,  91.9476776_jprb,  &
   91.9543228_jprb,  91.9608917_jprb,  91.9675446_jprb,  91.9741135_jprb,  91.9807358_jprb,  &
   91.9872818_jprb,  91.9934845_jprb,  91.9991302_jprb,  92.0047531_jprb,  92.0103531_jprb,  &
   92.0159225_jprb,  92.0216522_jprb,  92.0271378_jprb,  92.0327454_jprb,  92.0383835_jprb,  &
   92.0440292_jprb,  92.0495987_jprb,  92.0548172_jprb,  92.0597534_jprb,  92.0645676_jprb,  &
   92.0695038_jprb,  92.0744095_jprb,  92.0792236_jprb,  92.0842361_jprb,  92.0890656_jprb,  &
   92.0939178_jprb,  92.0987701_jprb,  92.1036224_jprb,  92.1082535_jprb,  92.1125641_jprb,  &
   92.1168137_jprb,  92.1210632_jprb,  92.1254425_jprb,  92.1297073_jprb,  92.1340332_jprb,  &
   92.1382599_jprb,  92.1425018_jprb,  92.1468964_jprb,  92.1511917_jprb,  92.1555557_jprb,  &
   92.1601181_jprb,  92.1647797_jprb,  92.1693420_jprb,  92.1739655_jprb,  92.1785965_jprb,  &
   92.1831818_jprb,  92.1877823_jprb,  92.1924210_jprb,  92.1970901_jprb,  92.2016602_jprb,  &
   92.2062988_jprb,  92.2111435_jprb,  92.2162476_jprb,  92.2212143_jprb,  92.2261887_jprb,  &
   92.2311478_jprb,  92.2361984_jprb,  92.2411270_jprb,  92.2461243_jprb,  92.2511292_jprb,  &
   92.2560883_jprb,  92.2610703_jprb,  92.2659225_jprb,  92.2704163_jprb,  92.2749329_jprb,  &
   92.2792664_jprb,  92.2837830_jprb,  92.2882614_jprb,  92.2926559_jprb,  92.2971191_jprb,  &
   92.3016281_jprb,  92.3061066_jprb,  92.3105774_jprb,  92.3150253_jprb,  92.3191605_jprb,  &
   92.3231735_jprb,  92.3272247_jprb,  92.3311996_jprb,  92.3353271_jprb,  92.3392868_jprb,  &
   92.3433228_jprb,  92.3472137_jprb,  92.3512650_jprb,  92.3553162_jprb,  92.3593063_jprb,  &
   92.3633499_jprb,  92.3673172_jprb,  92.3713760_jprb,  92.3754501_jprb,  92.3794556_jprb,  &
   92.3835297_jprb,  92.3875732_jprb,  92.3916855_jprb,  92.3956604_jprb,  92.3997269_jprb,  &
   92.4037857_jprb,  92.4078217_jprb,  92.4118805_jprb,  92.4151077_jprb,  92.4173508_jprb,  &
   92.4196320_jprb,  92.4218826_jprb,  92.4242096_jprb,  92.4264526_jprb,  92.4287872_jprb,  &
   92.4310226_jprb,  92.4332962_jprb,  92.4356079_jprb,  92.4379272_jprb,  92.4401932_jprb,  &
   92.4426575_jprb,  92.4460220_jprb,  92.4493713_jprb,  92.4527512_jprb,  92.4561005_jprb,  &
   92.4594421_jprb,  92.4627914_jprb,  92.4661255_jprb,  92.4694748_jprb,  92.4728317_jprb,  &
   92.4761200_jprb,  92.4795761_jprb,  92.4828339_jprb,  92.4866333_jprb,  92.4905472_jprb,  &
   92.4944153_jprb,  92.4984589_jprb,  92.5023422_jprb,  92.5063248_jprb,  92.5101852_jprb,  &
   92.5142136_jprb,  92.5180893_jprb,  92.5221176_jprb,  92.5259476_jprb,  92.5298843_jprb,  &
   92.5338211_jprb,  92.5378952_jprb,  92.5419693_jprb,  92.5460587_jprb,  92.5500717_jprb,  &
   92.5542374_jprb,  92.5583115_jprb,  92.5623703_jprb,  92.5664215_jprb,  92.5705414_jprb,  &
   92.5745773_jprb,  92.5786819_jprb,  92.5827103_jprb,  92.5867538_jprb,  92.5902863_jprb,  &
   92.5938416_jprb,  92.5972977_jprb,  92.6009445_jprb,  92.6043701_jprb,  92.6079559_jprb,  &
   92.6115189_jprb,  92.6150131_jprb,  92.6185455_jprb,  92.6220627_jprb,  92.6256180_jprb,  &
   92.6291885_jprb,  92.6327667_jprb,  92.6370697_jprb,  92.6415482_jprb,  92.6461105_jprb,  &
   92.6505737_jprb,  92.6549911_jprb,  92.6594696_jprb,  92.6639786_jprb,  92.6683807_jprb,  &
   92.6729279_jprb,  92.6774063_jprb,  92.6818924_jprb,  92.6863174_jprb,  92.6908417_jprb,  &
   92.6940384_jprb,  92.6965179_jprb,  92.6990891_jprb,  92.7017670_jprb,  92.7041779_jprb,  &
   92.7068024_jprb,  92.7093430_jprb,  92.7120132_jprb,  92.7144699_jprb,  92.7170486_jprb,  &
   92.7195892_jprb,  92.7220840_jprb,  92.7247009_jprb,  92.7271194_jprb,  92.7286530_jprb,  &
   92.7301483_jprb,  92.7317123_jprb,  92.7331924_jprb,  92.7347336_jprb,  92.7362823_jprb,  &
   92.7377930_jprb,  92.7392807_jprb,  92.7408218_jprb,  92.7423019_jprb,  92.7438126_jprb,  &
   92.7453232_jprb,  92.7468872_jprb,  92.7482605_jprb,  92.7494278_jprb,  92.7506027_jprb,  &
   92.7517319_jprb,  92.7528992_jprb,  92.7539825_jprb,  92.7552032_jprb,  92.7563629_jprb,  &
   92.7575226_jprb,  92.7586670_jprb,  92.7598114_jprb,  92.7609406_jprb,  92.7620544_jprb,  &
   92.7632904_jprb,  92.7643280_jprb,  92.7651291_jprb,  92.7659836_jprb,  92.7667923_jprb,  &
   92.7676926_jprb,  92.7684479_jprb,  92.7692871_jprb,  92.7701874_jprb,  92.7709122_jprb,  &
   92.7717514_jprb,  92.7725601_jprb,  92.7734451_jprb,  92.7743073_jprb,  92.7750931_jprb,  &
   92.7758942_jprb,  92.7763901_jprb,  92.7768555_jprb,  92.7773819_jprb,  92.7778778_jprb,  &
   92.7783890_jprb,  92.7787552_jprb,  92.7793198_jprb,  92.7798309_jprb,  92.7802963_jprb,  &
   92.7808380_jprb,  92.7812576_jprb,  92.7817688_jprb,  92.7822189_jprb,  92.7828293_jprb,  &
   92.7833786_jprb,  92.7840576_jprb,  92.7847214_jprb,  92.7854767_jprb,  92.7861328_jprb,  &
   92.7868881_jprb,  92.7875290_jprb,  92.7882614_jprb,  92.7889481_jprb,  92.7896805_jprb,  &
   92.7903214_jprb,  92.7910233_jprb,  92.7917328_jprb,  92.7924194_jprb,  92.7930756_jprb,  &
   92.7938843_jprb,  92.7946930_jprb,  92.7954330_jprb,  92.7962723_jprb,  92.7970200_jprb,  &
   92.7978058_jprb,  92.7985382_jprb,  92.7993317_jprb,  92.8001480_jprb,  92.8008652_jprb,  &
   92.8017426_jprb,  92.8024750_jprb,  92.8032532_jprb,  92.8040314_jprb,  92.8047943_jprb,  &
   92.8051987_jprb,  92.8055344_jprb,  92.8057709_jprb,  92.8061218_jprb,  92.8064117_jprb,  &
   92.8067398_jprb,  92.8069763_jprb,  92.8072968_jprb,  92.8076630_jprb,  92.8079224_jprb,  &
   92.8081436_jprb,  92.8085556_jprb,  92.8088379_jprb,  92.8091431_jprb,  92.8094940_jprb,  &
   92.8097610_jprb,  92.8100357_jprb,  92.8103561_jprb,  92.8107681_jprb,  92.8110580_jprb,  &
   92.8114090_jprb,  92.8118591_jprb,  92.8121338_jprb,  92.8124237_jprb,  92.8128510_jprb,  &
   92.8131866_jprb,  92.8134689_jprb,  92.8137970_jprb,  92.8141632_jprb,  92.8144608_jprb,  &
   92.8148727_jprb,  92.8151932_jprb,  92.8154831_jprb,  92.8157272_jprb,  92.8160477_jprb,  &
   92.8163910_jprb,  92.8166733_jprb,  92.8169861_jprb,  92.8173370_jprb,  92.8176422_jprb,  &
   92.8179703_jprb,  92.8182602_jprb,  92.8185425_jprb,  92.8188782_jprb,  92.8191757_jprb,  &
   92.8194275_jprb,  92.8197556_jprb,  92.8200455_jprb,  92.8202972_jprb,  92.8205948_jprb,  &
   92.8208923_jprb,  92.8211517_jprb,  92.8214417_jprb,  92.8217087_jprb,  92.8220291_jprb,  &
   92.8222809_jprb,  92.8225250_jprb,  92.8229141_jprb,  92.8231583_jprb,  92.8233490_jprb,  &
   92.8236618_jprb,  92.8239517_jprb,  92.8241653_jprb,  92.8240662_jprb,  92.8238602_jprb,  &
   92.8238068_jprb,  92.8236160_jprb,  92.8234863_jprb,  92.8233719_jprb,  92.8232193_jprb,  &
   92.8231201_jprb,  92.8230209_jprb,  92.8228531_jprb,  92.8226929_jprb,  92.8226089_jprb,  &
   92.8224869_jprb,  92.8223114_jprb,  92.8222198_jprb,  92.8220901_jprb,  92.8222198_jprb,  &
   92.8223190_jprb,  92.8224106_jprb,  92.8226471_jprb,  92.8227615_jprb,  92.8229446_jprb,  &
   92.8230591_jprb,  92.8231812_jprb,  92.8232803_jprb,  92.8234787_jprb,  92.8236237_jprb,  &
   92.8238144_jprb,  92.8239059_jprb,  92.8240433_jprb,  92.8241730_jprb,  92.8243103_jprb,  &
   92.8244858_jprb,  92.8244095_jprb,  92.8242722_jprb,  92.8242645_jprb,  92.8242645_jprb,  &
   92.8241348_jprb,  92.8240433_jprb,  92.8239899_jprb,  92.8239365_jprb,  92.8238907_jprb,  &
   92.8237686_jprb,  92.8237228_jprb,  92.8237839_jprb,  92.8236847_jprb,  92.8236160_jprb,  &
   92.8235016_jprb,  92.8234482_jprb,  92.8233795_jprb,  92.8231659_jprb,  92.8228989_jprb,  &
   92.8226776_jprb,  92.8225021_jprb,  92.8222809_jprb,  92.8220749_jprb,  92.8218155_jprb,  &
   92.8215256_jprb,  92.8213120_jprb,  92.8211517_jprb,  92.8207932_jprb,  92.8206711_jprb,  &
   92.8203964_jprb,  92.8202362_jprb,  92.8199539_jprb,  92.8196945_jprb,  92.8195724_jprb,  &
   92.8189545_jprb,  92.8184204_jprb,  92.8177643_jprb,  92.8172760_jprb,  92.8166885_jprb,  &
   92.8160553_jprb,  92.8154068_jprb,  92.8147430_jprb,  92.8142166_jprb,  92.8136139_jprb,  &
   92.8129883_jprb,  92.8124161_jprb,  92.8118286_jprb,  92.8111954_jprb,  92.8106232_jprb,  &
   92.8100128_jprb,  92.8093719_jprb,  92.8088608_jprb,  92.8083191_jprb,  92.8077774_jprb,  &
   92.8072662_jprb,  92.8067322_jprb,  92.8062057_jprb,  92.8056793_jprb,  92.8050919_jprb,  &
   92.8046341_jprb,  92.8040314_jprb,  92.8036499_jprb,  92.8030930_jprb,  92.8025436_jprb,  &
   92.8019791_jprb,  92.8015060_jprb,  92.8009033_jprb,  92.8004532_jprb,  92.7999268_jprb,  &
   92.7992630_jprb,  92.7985535_jprb,  92.7977142_jprb,  92.7970200_jprb,  92.7961426_jprb,  &
   92.7953568_jprb,  92.7946091_jprb,  92.7938690_jprb,  92.7930603_jprb,  92.7922821_jprb,  &
   92.7915115_jprb,  92.7907562_jprb,  92.7899170_jprb,  92.7891312_jprb,  92.7883377_jprb,  &
   92.7875977_jprb,  92.7868118_jprb,  92.7860107_jprb,  92.7851715_jprb,  92.7843552_jprb,  &
   92.7834930_jprb,  92.7826843_jprb,  92.7818222_jprb,  92.7810135_jprb,  92.7800827_jprb,  &
   92.7792282_jprb,  92.7784195_jprb,  92.7775497_jprb,  92.7766342_jprb,  92.7759323_jprb,  &
   92.7750168_jprb,  92.7741318_jprb,  92.7732468_jprb,  92.7724457_jprb,  92.7716370_jprb,  &
   92.7707748_jprb,  92.7699356_jprb,  92.7690125_jprb,  92.7680435_jprb,  92.7670441_jprb,  &
   92.7661133_jprb,  92.7651138_jprb,  92.7640686_jprb,  92.7631760_jprb,  92.7622070_jprb,  &
   92.7612076_jprb,  92.7601471_jprb,  92.7592545_jprb,  92.7582550_jprb,  92.7572327_jprb,  &
   92.7562714_jprb,  92.7552490_jprb,  92.7543488_jprb,  92.7533493_jprb,  92.7523651_jprb,  &
   92.7513962_jprb,  92.7503738_jprb,  92.7493820_jprb,  92.7484436_jprb,  92.7474365_jprb,  &
   92.7464447_jprb,  92.7454834_jprb,  92.7444992_jprb,  92.7435150_jprb,  92.7425003_jprb,  &
   92.7415543_jprb,  92.7405319_jprb,  92.7396011_jprb,  92.7385483_jprb,  92.7376251_jprb,  &
   92.7365570_jprb,  92.7356796_jprb,  92.7346497_jprb,  92.7336349_jprb,  92.7326736_jprb,  &
   92.7316284_jprb,  92.7306213_jprb,  92.7295074_jprb,  92.7284851_jprb,  92.7274323_jprb,  &
   92.7263641_jprb,  92.7253265_jprb,  92.7242279_jprb,  92.7232361_jprb,  92.7221680_jprb,  &
   92.7211380_jprb,  92.7199936_jprb,  92.7189636_jprb,  92.7179337_jprb,  92.7168121_jprb,  &
   92.7158432_jprb,  92.7147598_jprb,  92.7137070_jprb,  92.7125931_jprb,  92.7115479_jprb,  &
   92.7105103_jprb,  92.7093201_jprb,  92.7081070_jprb,  92.7069397_jprb,  92.7056198_jprb,  &
   92.7044220_jprb,  92.7033234_jprb,  92.7021179_jprb,  92.7009583_jprb,  92.6996536_jprb,  &
   92.6984558_jprb,  92.6972733_jprb,  92.6960678_jprb,  92.6949539_jprb,  92.6936951_jprb,  &
   92.6924896_jprb,  92.6912842_jprb,  92.6901016_jprb,  92.6888580_jprb,  92.6876984_jprb,  &
   92.6864548_jprb,  92.6851730_jprb,  92.6839600_jprb,  92.6826401_jprb,  92.6812973_jprb,  &
   92.6800461_jprb,  92.6787109_jprb,  92.6774063_jprb,  92.6761169_jprb,  92.6748199_jprb,  &
   92.6735687_jprb,  92.6722260_jprb,  92.6709442_jprb,  92.6695938_jprb,  92.6683578_jprb,  &
   92.6670761_jprb,  92.6657486_jprb,  92.6644287_jprb,  92.6631775_jprb,  92.6618881_jprb,  &
   92.6605530_jprb,  92.6592255_jprb,  92.6579361_jprb,  92.6566544_jprb,  92.6554108_jprb,  &
   92.6539612_jprb,  92.6527100_jprb,  92.6514053_jprb,  92.6501007_jprb,  92.6488037_jprb,  &
   92.6475067_jprb,  92.6461945_jprb,  92.6448212_jprb,  92.6435394_jprb,  92.6422653_jprb,  &
   92.6409760_jprb,  92.6396255_jprb,  92.6383209_jprb,  92.6370392_jprb,  92.6357117_jprb,  &
   92.6343307_jprb,  92.6329803_jprb,  92.6317139_jprb,  92.6302338_jprb,  92.6287384_jprb,  &
   92.6272812_jprb,  92.6257172_jprb,  92.6242371_jprb,  92.6227570_jprb,  92.6212845_jprb,  &
   92.6198349_jprb,  92.6183395_jprb,  92.6168213_jprb,  92.6154251_jprb,  92.6138000_jprb,  &
   92.6124268_jprb,  92.6108627_jprb,  92.6094055_jprb,  92.6079254_jprb,  92.6064758_jprb,  &
   92.6049805_jprb,  92.6034241_jprb,  92.6019363_jprb,  92.6005554_jprb,  92.5988617_jprb,  &
   92.5972672_jprb,  92.5956573_jprb,  92.5939255_jprb,  92.5922852_jprb,  92.5907135_jprb,  &
   92.5890121_jprb,  92.5874786_jprb,  92.5857315_jprb,  92.5841370_jprb,  92.5825043_jprb,  &
   92.5808487_jprb,  92.5792999_jprb,  92.5776062_jprb,  92.5759583_jprb,  92.5743332_jprb,  &
   92.5726852_jprb,  92.5709839_jprb,  92.5693436_jprb,  92.5677414_jprb,  92.5661392_jprb,  &
   92.5644455_jprb,  92.5628967_jprb,  92.5612946_jprb,  92.5596466_jprb,  92.5582275_jprb,  &
   92.5566711_jprb,  92.5550613_jprb,  92.5535126_jprb,  92.5519333_jprb,  92.5503693_jprb,  &
   92.5488739_jprb,  92.5472717_jprb,  92.5457458_jprb,  92.5441437_jprb,  92.5426178_jprb,  &
   92.5410385_jprb,  92.5394592_jprb,  92.5379333_jprb,  92.5363235_jprb,  92.5348206_jprb,  &
   92.5332642_jprb,  92.5316467_jprb,  92.5300980_jprb,  92.5284958_jprb,  92.5267868_jprb,  &
   92.5251007_jprb,  92.5233459_jprb,  92.5216141_jprb,  92.5199203_jprb,  92.5182037_jprb,  &
   92.5164719_jprb,  92.5147858_jprb,  92.5130539_jprb,  92.5113525_jprb,  92.5096130_jprb,  &
   92.5078812_jprb,  92.5061874_jprb,  92.5044708_jprb,  92.5027390_jprb,  92.5010452_jprb,  &
   92.4993134_jprb,  92.4975815_jprb,  92.4959183_jprb,  92.4942017_jprb,  92.4923859_jprb,  &
   92.4907684_jprb,  92.4888840_jprb,  92.4870071_jprb,  92.4851761_jprb,  92.4833221_jprb,  &
   92.4813614_jprb,  92.4794998_jprb,  92.4776077_jprb,  92.4757080_jprb,  92.4737701_jprb,  &
   92.4719849_jprb,  92.4700012_jprb,  92.4681702_jprb,  92.4663010_jprb,  92.4643631_jprb,  &
   92.4624786_jprb,  92.4606628_jprb,  92.4586868_jprb,  92.4568405_jprb,  92.4549713_jprb,  &
   92.4531097_jprb,  92.4511108_jprb,  92.4492188_jprb,  92.4473190_jprb,  92.4454803_jprb,  &
   92.4435272_jprb,  92.4414749_jprb,  92.4394608_jprb,  92.4375687_jprb,  92.4355621_jprb,  &
   92.4335251_jprb,  92.4314575_jprb,  92.4295349_jprb,  92.4275208_jprb,  92.4255829_jprb,  &
   92.4235229_jprb,  92.4215622_jprb,  92.4195099_jprb,  92.4175949_jprb,  92.4155807_jprb,  &
   92.4135895_jprb,  92.4115982_jprb,  92.4095840_jprb,  92.4076080_jprb,  92.4055634_jprb,  &
   92.4035950_jprb,  92.4015656_jprb,  92.3995895_jprb,  92.3976288_jprb,  92.3957291_jprb,  &
   92.3937683_jprb,  92.3919067_jprb,  92.3899384_jprb,  92.3880539_jprb,  92.3861465_jprb,  &
   92.3841858_jprb,  92.3822403_jprb,  92.3803406_jprb,  92.3784256_jprb,  92.3764877_jprb,  &
   92.3745728_jprb,  92.3726578_jprb,  92.3707657_jprb,  92.3687820_jprb,  92.3668900_jprb,  &
   92.3649597_jprb,  92.3630524_jprb,  92.3611221_jprb,  92.3591843_jprb,  92.3572464_jprb,  &
   92.3553467_jprb,  92.3534470_jprb,  92.3516006_jprb,  92.3495407_jprb,  92.3475266_jprb,  &
   92.3455505_jprb,  92.3436584_jprb,  92.3416443_jprb,  92.3396454_jprb,  92.3376541_jprb,  &
   92.3357315_jprb,  92.3337479_jprb,  92.3316956_jprb,  92.3298111_jprb,  92.3277512_jprb,  &
   92.3259201_jprb,  92.3238297_jprb,  92.3219376_jprb,  92.3199234_jprb,  92.3179779_jprb,  &
   92.3159943_jprb,  92.3140717_jprb,  92.3119736_jprb,  92.3100281_jprb,  92.3081055_jprb,  &
   92.3060455_jprb,  92.3041077_jprb,  92.3021469_jprb,  92.3002090_jprb,  92.2982788_jprb,  &
   92.2963562_jprb,  92.2944641_jprb,  92.2925720_jprb,  92.2905426_jprb,  92.2886505_jprb,  &
   92.2867203_jprb,  92.2847824_jprb,  92.2828903_jprb,  92.2809677_jprb,  92.2790375_jprb,  &
   92.2771149_jprb,  92.2752762_jprb,  92.2733307_jprb,  92.2713852_jprb,  92.2694778_jprb,  &
   92.2675095_jprb,  92.2655869_jprb,  92.2636414_jprb,  92.2617340_jprb,  92.2597885_jprb,  &
   92.2579575_jprb,  92.2559738_jprb,  92.2540359_jprb,  92.2521667_jprb,  92.2503128_jprb,  &
   92.2482986_jprb,  92.2464676_jprb,  92.2445297_jprb,  92.2426300_jprb,  92.2407532_jprb,  &
   92.2388840_jprb,  92.2370300_jprb,  92.2350540_jprb,  92.2331619_jprb,  92.2312469_jprb,  &
   92.2293091_jprb,  92.2274551_jprb,  92.2256088_jprb,  92.2236786_jprb,  92.2217789_jprb,  &
   92.2199249_jprb,  92.2179794_jprb,  92.2161026_jprb,  92.2141647_jprb,  92.2122269_jprb,  &
   92.2103577_jprb,  92.2085190_jprb,  92.2065582_jprb,  92.2046204_jprb,  92.2028122_jprb,  &
   92.2009964_jprb,  92.1991119_jprb,  92.1973267_jprb,  92.1955948_jprb,  92.1938095_jprb,  &
   92.1919708_jprb,  92.1902008_jprb,  92.1884308_jprb,  92.1865311_jprb,  92.1847610_jprb,  &
   92.1829529_jprb,  92.1811829_jprb,  92.1793976_jprb,  92.1775284_jprb,  92.1757050_jprb,  &
   92.1738739_jprb,  92.1721725_jprb,  92.1703415_jprb,  92.1684952_jprb,  92.1667557_jprb,  &
   92.1649094_jprb,  92.1630936_jprb,  92.1612930_jprb,  92.1595154_jprb,  92.1576309_jprb,  &
   92.1559296_jprb,  92.1542740_jprb,  92.1529999_jprb,  92.1516876_jprb,  92.1503143_jprb,  &
   92.1490707_jprb,  92.1477203_jprb,  92.1463776_jprb,  92.1450882_jprb,  92.1437988_jprb,  &
   92.1424332_jprb,  92.1410370_jprb,  92.1398315_jprb,  92.1384354_jprb,  92.1371460_jprb,  &
   92.1357956_jprb,  92.1345062_jprb,  92.1331711_jprb,  92.1318893_jprb,  92.1305618_jprb,  &
   92.1291962_jprb,  92.1278610_jprb,  92.1265640_jprb,  92.1253662_jprb,  92.1239777_jprb,  &
   92.1226120_jprb,  92.1213455_jprb,  92.1200485_jprb,  92.1186295_jprb,  92.1167526_jprb,  &
   92.1148758_jprb,  92.1130295_jprb,  92.1110992_jprb,  92.1092682_jprb,  92.1074524_jprb,  &
   92.1055603_jprb,  92.1037369_jprb,  92.1017303_jprb,  92.0999222_jprb,  92.0979919_jprb,  &
   92.0961685_jprb,  92.0943680_jprb,  92.0924301_jprb,  92.0905762_jprb,  92.0887222_jprb,  &
   92.0868301_jprb,  92.0849457_jprb,  92.0830917_jprb,  92.0812302_jprb,  92.0793839_jprb,  &
   92.0775452_jprb,  92.0755615_jprb,  92.0737457_jprb,  92.0718918_jprb,  92.0699997_jprb,  &
   92.0681534_jprb,  92.0663605_jprb,  92.0644760_jprb,  92.0626450_jprb,  92.0608444_jprb,  &
   92.0589905_jprb,  92.0571976_jprb,  92.0553589_jprb,  92.0534668_jprb,  92.0516815_jprb,  &
   92.0498810_jprb,  92.0480804_jprb,  92.0462036_jprb,  92.0444107_jprb,  92.0425415_jprb,  &
   92.0407410_jprb,  92.0389023_jprb,  92.0371246_jprb,  92.0353012_jprb,  92.0334549_jprb,  &
   92.0316925_jprb,  92.0297623_jprb,  92.0279846_jprb,  92.0261383_jprb,  92.0243225_jprb,  &
   92.0225220_jprb,  92.0207443_jprb,  92.0188522_jprb,  92.0170670_jprb,  92.0152588_jprb,  &
   92.0134048_jprb,  92.0115891_jprb,  92.0098343_jprb,  92.0079880_jprb,  92.0062180_jprb,  &
   92.0043411_jprb,  92.0025635_jprb,  92.0007782_jprb,  91.9990082_jprb,  91.9971542_jprb,  &
   91.9954147_jprb,  91.9935532_jprb,  91.9917831_jprb,  91.9899445_jprb,  91.9881744_jprb,  &
   91.9863815_jprb,  91.9845886_jprb,  91.9827499_jprb,  91.9809875_jprb,  91.9791718_jprb,  &
   91.9773026_jprb,  91.9755478_jprb,  91.9737320_jprb,  91.9719315_jprb,  91.9701462_jprb,  &
   91.9683456_jprb,  91.9665070_jprb,  91.9647598_jprb,  91.9628830_jprb,  91.9611359_jprb,  &
   91.9592056_jprb,  91.9571228_jprb,  91.9552460_jprb,  91.9532852_jprb,  91.9513397_jprb,  &
   91.9493942_jprb,  91.9474411_jprb,  91.9454575_jprb,  91.9434738_jprb,  91.9415817_jprb,  &
   91.9396591_jprb,  91.9376602_jprb,  91.9356613_jprb,  91.9337463_jprb,  91.9317780_jprb,  &
   91.9298553_jprb,  91.9279404_jprb,  91.9259109_jprb,  91.9239960_jprb,  91.9220276_jprb,  &
   91.9200516_jprb,  91.9181290_jprb,  91.9160690_jprb,  91.9142227_jprb,  91.9122086_jprb,  &
   91.9102325_jprb,  91.9083099_jprb,  91.9064026_jprb,  91.9044418_jprb,  91.9024277_jprb,  &
   91.9006958_jprb,  91.8988953_jprb,  91.8971863_jprb,  91.8954239_jprb,  91.8936615_jprb,  &
   91.8919144_jprb,  91.8901672_jprb,  91.8884201_jprb,  91.8866730_jprb,  91.8848953_jprb,  &
   91.8831253_jprb,  91.8813553_jprb,  91.8796234_jprb,  91.8778915_jprb,  91.8761215_jprb,  &
   91.8743668_jprb,  91.8726883_jprb,  91.8709183_jprb,  91.8691864_jprb,  91.8673935_jprb,  &
   91.8656235_jprb,  91.8638992_jprb,  91.8621445_jprb,  91.8603897_jprb,  91.8586349_jprb,  &
   91.8568649_jprb,  91.8550720_jprb,  91.8534164_jprb,  91.8516464_jprb,  91.8498154_jprb,  &
   91.8481064_jprb/

CONTAINS

!------------------------------------------
! Routines for initialising database
!------------------------------------------

  !> Initialise a BRDF atlas data structure. By default the atlas data are
  !! general and can be used for any sensor with visible/near-IR channels.
  !! If the last four optional arguments are specified then the atlas data are
  !! initialised for use with the specified instrument only: this makes calls
  !! to obtain BRDFs from the atlas much faster.
  !! @param[in]       path            path to atlas data files
  !! @param[in]       imonth          month of data to read (1-12)
  !! @param[in]       verbose         flag to turn verbose output on/off
  !! @param[in,out]   atlas           BRDF atlas data structure to initialise
  !! @param[out]      err             status on exit
  !! @param[in]       instr_wavenum   channel wavenumber list, optional
  !! @param[in]       platform_id     platform ID from associated rtcoef structure, optional
  !! @param[in]       sat_id          satellite ID from associated rtcoef structure, optional
  !! @param[in]       inst_id         instrument ID from associated rtcoef structure, optional
  SUBROUTINE rttov_visnirbrdf_init( &
                    path,          &
                    imonth,        &
                    verbose,       &
                    atlas,         &
                    err,           &
                    instr_wavenum, &
                    platform_id,   &
                    sat_id,        &
                    inst_id)

    ! Description:
    ! initialize the rttov vis/nir brdf algorithm by (1) reading in the MODIS VIS/NIR Global
    ! BRDF kernel parameters data and (2) the eigenvectors of the laboratory spectra, and (2) make some
    ! precalculations for the PCA regression
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0     03/01/2012 original code J. Vidot

    USE rttov_lapack_mod, ONLY : dgesv

    CHARACTER(LEN=*),      INTENT(IN)           :: path
    INTEGER(KIND=jpim),    INTENT(IN)           :: imonth
    LOGICAL(KIND=jplm),    INTENT(IN)           :: verbose
    TYPE(brdf_atlas_data), INTENT(INOUT)        :: atlas
    INTEGER(KIND=jpim),    INTENT(OUT)          :: err
    REAL(KIND=jprb),       INTENT(IN), OPTIONAL :: instr_wavenum(:)
    INTEGER(KIND=jpim),    INTENT(IN), OPTIONAL :: platform_id
    INTEGER(KIND=jpim),    INTENT(IN), OPTIONAL :: sat_id
    INTEGER(KIND=jpim),    INTENT(IN), OPTIONAL :: inst_id

    INTEGER(KIND=jpim) :: i
    INTEGER            :: ipiv(numpcs), info                  ! LAPACK arguments, no KIND
    INTEGER, PARAMETER :: npc = numpcs

    REAL(KIND=jprb)  :: A(numpcs,hngpnts), AT(hngpnts,numpcs)
    DOUBLE PRECISION :: B(numpcs,numpcs), C(numpcs,numpcs)    ! LAPACK arguments, no KIND

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

    CALL rttov_visnirbrdf_nullify(atlas)

    WRITE(cyear,'(i4)') db_ver_year
    WRITE(cmonth,'(i2.2)') imonth
    WRITE(errmsg,'(a)') ''

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
    THROWM(err.NE.0,'RTTOV must be compiled with HDF5 to use the BRDF atlas')
#endif
#endif

    !-------------------------------------------------------------------
    !  reading the eigienvectors of the 126 selected laboratory spectra
    !-------------------------------------------------------------------
    ALLOCATE(atlas%pcu(numpcs,numwave), atlas%pcu_snow(numpcs,numwave))
    fn=TRIM(path)//'CMSsolarbrdf_labeigvects_V1.1'
#ifdef _RTTOV_HDF
    fn = TRIM(fn)//'.H5'
#else
    fn = TRIM(fn)//'.nc'
#endif
    INQUIRE(FILE=fn, EXIST=file_exists)
    IF (file_exists) THEN
      CALL rttov_brdf_read_labeigvects(err, TRIM(fn), atlas)
      THROW(err.NE.0)
    ELSE
      err = errorstatus_fatal
      THROWM(err .NE. errorstatus_success, 'BRDF eigenvector file not found: '//TRIM(fn)//TRIM(errmsg))
    ENDIF
    IF (verbose) INFO('Using BRDF eigenvector file: '//TRIM(fn))

    !----------------------------------------------------------------------------
    ! reading the 0.1 degree resolution MODIS BRDF kernel parameters
    !----------------------------------------------------------------------------
    fn=TRIM(path)//'CMSsolarbrdf_V1.0_0.1deg_'//cyear//cmonth//'_mask'
#ifdef _RTTOV_HDF
    fn = TRIM(fn)//'.H5'
#else
    fn = TRIM(fn)//'.nc'
#endif
    INQUIRE(FILE=fn, EXIST=file_exists)
    IF (file_exists) THEN
      CALL rttov_brdf_read_brdf(err, TRIM(fn), atlas)
      THROW(err.NE.0)
    ELSE
      err = errorstatus_fatal
      THROWM(err .NE. errorstatus_success, 'BRDF data file not found: '//TRIM(fn)//TRIM(errmsg))
    ENDIF
    IF (verbose) INFO('Using BRDF coefs: '//TRIM(fn))

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

    !--------------------------------------------------
    ! create hsr_wavenum array
    !--------------------------------------------------
    !  numwave = 2101, hsr_wavenum(4000:25000:10)

    DO i = 1, numwave
      hsr_wavenum(i) = 4000._jprb + (i - 1) * 10._jprb
    ENDDO

    !---------------------------------------------------
    ! compute D matrix for HSR brdf construction
    !---------------------------------------------------

    ALLOCATE(atlas%D(hngpnts,numpcs), atlas%D_snow(hngpnts,numpcs))

    ! define matrix A
    A = pcu_modres(1:numpcs,1:hngpnts)

    ! calculate transpose matrix A : AT [6X4]
    AT = TRANSPOSE(A)

    !  B=A*A' [4X6]X[6X4]   first: row second: column
    B(:,:) = 0.0_jprb
    B = MATMUL(A,AT)

    ! computing inverse of B
    C(:,:) = 0.0_jprb
    DO i = 1, numpcs
      C(i,i) = 1.0_jprb
    ENDDO
    CALL DGESV(npc, npc, B, npc, ipiv, C, npc, info)

    ! compute D=A'*inv(A*A')=A'*inv(B)= AT * C [6X4] * [4X4]
    atlas%D = MATMUL(AT,C)


    ! define matrix A
    A = pcu_modres_snow(1:numpcs,1:hngpnts)

    ! calculate transpose matrix A : AT [6X4]
    AT = TRANSPOSE(A)

    !  B=A*A' [4X6]X[6X4]   first: row second: column
    B(:,:) = 0.0_jprb
    B = MATMUL(A,AT)

    ! computing inverse of B
    C(:,:) = 0.0_jprb
    DO i = 1, numpcs
      C(i,i) = 1.0_jprb
    ENDDO
    CALL DGESV(npc, npc, B, npc, IPIV, C, npc, info)

    ! compute D=A'*inv(A*A')=A'*inv(B)= AT * C [6X4] * [4X4]
    atlas%D_snow = MATMUL(AT,C)

    IF (atlas%single_inst) THEN
      ! If a channel wavenumber list is supplied for a particular instrument
      ! we can retrieve BRDFs much faster.

      ! Interpolate hsr data onto the channel wavenumbers
      atlas%ncoefchans = SIZE(instr_wavenum)
      ALLOCATE(atlas%pcu_int(numpcs,atlas%ncoefchans),         &
               atlas%pcm_int(atlas%ncoefchans),                &
               atlas%pcu_int_snow(numpcs,atlas%ncoefchans),    &
               atlas%pcm_int_snow(atlas%ncoefchans),           &
               atlas%coastal_waters_ref_int(atlas%ncoefchans), &
               atlas%ocean_waters_ref_int(atlas%ncoefchans))
      CALL rttov_brdf_hsr_interp(instr_wavenum(:), atlas)
    ENDIF

    CATCH
  END SUBROUTINE rttov_visnirbrdf_init

#ifdef _RTTOV_HDF
  !> Read BRDF atlas data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   BRDF atlas data structure
  SUBROUTINE rttov_brdf_read_brdf(err, fn, atlas)

    ! Description:
    ! read the 0.1 degree resolution MODIS BRDF kernel parameters data
    ! from the HDF5 file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.1      03/06/2014 original code J. Vidot
    !

    INTEGER(KIND=jpim),    INTENT(OUT)   :: err
    CHARACTER(LEN=*),      INTENT(IN)    :: fn    ! filename including full path
    TYPE(brdf_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpis) , POINTER :: brdf_x(:,:)
    INTEGER(KIND=jpim)           :: indexlut
    INTEGER(KIND=jpim)           :: i,j

#include "rttov_hdf_load.interface"

    TRY

    CALL RTTOV_HDF_LOAD( ERR, fn, "/BRDF", SNAME="FLAG", PIS2=atlas%brdfmodis_flag)
    THROWM(ERR.NE.0,"Cannot load BRDF flag from "//TRIM(FN))

    atlas%nb_lats = SIZE(atlas%brdfmodis_flag,1)
    atlas%nb_lons = SIZE(atlas%brdfmodis_flag,2)

    ALLOCATE(atlas%brdfmodis_lut(atlas%nb_lats,atlas%nb_lons))

    ! Generate the look-up table into the emissivity data
    atlas%brdfmodis_lut(:,:) = -1_jpim
    indexlut = 1_jpim
    DO i = 1, atlas%nb_lons
      DO j = 1, atlas%nb_lats
        IF (atlas%brdfmodis_flag(j,i) > 0) THEN
          atlas%brdfmodis_lut(j,i) = indexlut
          indexlut = indexlut + 1_jpim
        ENDIF
      ENDDO
    ENDDO

    CALL RTTOV_HDF_LOAD( ERR, fn, "/BRDF", SNAME="FISO", PIS2=brdf_x)
    THROWM(ERR.NE.0,"Cannot BRDF fiso kernel from "//TRIM(FN))

    atlas%nb_pack = SIZE(brdf_x,1)
    IF (SIZE(brdf_x,2) .NE. hngpnts) THEN
      DEALLOCATE(brdf_x)
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in BRDF data: hngpnts")
    ENDIF
    ALLOCATE(atlas%brdfmodis(3,atlas%nb_pack,hngpnts))

    atlas%brdfmodis(1,:,:) = brdf_x
    DEALLOCATE(brdf_x)

    CALL RTTOV_HDF_LOAD( ERR, fn, "/BRDF", SNAME="FVOL", PIS2=brdf_x)
    THROWM(ERR.NE.0,"Cannot BRDF fvol kernel from "//TRIM(FN))
    atlas%brdfmodis(2,:,:) = brdf_x
    DEALLOCATE(brdf_x)

    CALL RTTOV_HDF_LOAD( ERR, fn, "/BRDF", SNAME="FGEO", PIS2=brdf_x)
    THROWM(ERR.NE.0,"Cannot BRDF fgeo kernel from "//TRIM(FN))
    atlas%brdfmodis(3,:,:) = brdf_x
    DEALLOCATE(brdf_x)
    NULLIFY(brdf_x)

!     CALL RTTOV_HDF_LOAD( ERR, fn, "/BRDF", SNAME="KERNEL", PIS3=atlas%brdfmodis)
!     THROWM(ERR.NE.0,"Cannot BRDF kernel from "//TRIM(FN))
!
!     atlas%nb_pack = SIZE(atlas%brdfmodis,2)
!     IF (SIZE(atlas%brdfmodis,3) .NE. hngpnts) THEN
!       DEALLOCATE(atlas%brdfmodis)
!       err = errorstatus_fatal
!       THROWM(err.NE.0,"Inconsistent dimensions in BRDF data: hngpnts")
!     ENDIF

    CALL RTTOV_HDF_LOAD( ERR, fn, "/BRDF", SNAME="FISO_SFAC", PRM1=atlas%sfac_fiso)
    THROWM(ERR.NE.0,"Cannot BRDF fiso kernel offset from "//TRIM(FN))

    CALL RTTOV_HDF_LOAD( ERR, fn, "/BRDF", SNAME="FISO_OFFS", PRM1=atlas%offs_fiso)
    THROWM(ERR.NE.0,"Cannot BRDF fiso kernel scale factor from "//TRIM(FN))

    CALL RTTOV_HDF_LOAD( ERR, fn, "/BRDF", SNAME="FVOL_SFAC", PRM1=atlas%sfac_fvol)
    THROWM(ERR.NE.0,"Cannot BRDF fvol kernel offset from "//TRIM(FN))

    CALL RTTOV_HDF_LOAD( ERR, fn, "/BRDF", SNAME="FVOL_OFFS", PRM1=atlas%offs_fvol)
    THROWM(ERR.NE.0,"Cannot BRDF fvol kernel scale factor from "//TRIM(FN))

    CALL RTTOV_HDF_LOAD( ERR, fn, "/BRDF", SNAME="FGEO_SFAC", PRM1=atlas%sfac_fgeo)
    THROWM(ERR.NE.0,"Cannot BRDF fgeo kernel offset from "//TRIM(FN))

    CALL RTTOV_HDF_LOAD( ERR, fn, "/BRDF", SNAME="FGEO_OFFS", PRM1=atlas%offs_fgeo)
    THROWM(ERR.NE.0,"Cannot BRDF fgeo kernel scale factor from "//TRIM(FN))

    CATCH
  END SUBROUTINE rttov_brdf_read_brdf
#else
#ifdef _RTTOV_NETCDF
  !> Read BRDF atlas data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   BRDF atlas data structure
  SUBROUTINE rttov_brdf_read_brdf(err, fn, atlas)

    ! Description:
    ! read the 0.1 degree resolution MODIS BRDF kernel parameters data
    ! from the netCDF file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0      03/01/2012 original code J. Vidot
    !

    INTEGER(KIND=jpim),    INTENT(OUT)   :: err
    CHARACTER(LEN=*),      INTENT(IN)    :: fn
    TYPE(brdf_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpim) :: nvars     ! number of variables
    INTEGER(KIND=jpim) :: ndims     ! number of dimensions
    INTEGER(KIND=jpim) :: errstat   ! error code
    INTEGER(KIND=jpim) :: recdim    ! record dimension
    INTEGER(KIND=jpim) :: nc_dim(4) ! hng_pnt, lats, lons, pack_len

    INTEGER(KIND=jpim) :: i,j,k
    INTEGER(KIND=jpim) :: ncid, ngatts, nrecs, varid

    INTEGER(KIND=jpis), ALLOCATABLE :: brdf_ch(:)

    CHARACTER (LEN=1024) :: strbuf ! string buffer for var
    CHARACTER (LEN=6)    :: cfld
    INTEGER(KIND=jpim)   :: indexlut

    TRY

    ! Open netCDF file.
    errstat = nf_open(TRIM(fn), nf_nowrite,ncid)

    ! Get info on the record dimension for this file.
    errstat = nf_inq(ncid, ndims, nvars, ngatts, recdim)

    DO recdim=1,ndims
      errstat = nf_inq_dim(ncid, recdim, strbuf, nrecs)
      nc_dim(recdim) = nrecs
    ENDDO

    IF (nc_dim(1) .NE. hngpnts) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in BRDF data: hngpnts")
    ENDIF
    atlas%nb_lats = nc_dim(2)
    atlas%nb_lons = nc_dim(3)
    atlas%nb_pack = nc_dim(4)

    ALLOCATE(atlas%brdfmodis_flag(atlas%nb_lats,atlas%nb_lons), &
             atlas%brdfmodis_lut(atlas%nb_lats,atlas%nb_lons),  &
             atlas%brdfmodis(3,atlas%nb_pack,hngpnts),    &
             brdf_ch(atlas%nb_pack))

    ALLOCATE(atlas%sfac_fiso(hngpnts), &
             atlas%offs_fiso(hngpnts), &
             atlas%sfac_fvol(hngpnts), &
             atlas%offs_fvol(hngpnts), &
             atlas%sfac_fgeo(hngpnts), &
             atlas%offs_fgeo(hngpnts))

    ! Retrieve brdf database flag value
    errstat = nf_inq_varid(ncid, 'brdf_flag', varid)
    errstat = nf_get_var_int2(ncid, varid, atlas%brdfmodis_flag)

    ! Generate the look-up table into the emissivity data
    atlas%brdfmodis_lut(:,:) = -1_jpim
    indexlut = 1_jpim
    DO i = 1, atlas%nb_lons
      DO j = 1, atlas%nb_lats
        IF (atlas%brdfmodis_flag(j,i) > 0) THEN
          atlas%brdfmodis_lut(j,i) = indexlut
          indexlut = indexlut + 1_jpim
        ENDIF
      ENDDO
    ENDDO

    ! Retrieve database of 7 MODIS hinge-point FISO values
    Do k=1,nc_dim(1)
      WRITE(cfld,'("fiso",i1," ")') k
      errstat = nf_inq_varid (ncid,cfld,varid)
      errstat = nf_get_var_int2(ncid,varid,brdf_ch)
      errstat = nf_get_att_real(ncid,varid,'scale_factor',atlas%sfac_fiso(k))
      errstat = nf_get_att_real(ncid,varid,'add_offset',atlas%offs_fiso(k))
      atlas%brdfmodis(1,:,k) = brdf_ch
    ENDDO

    ! Retrieve database of 7 MODIS hinge-point FVOL values
    Do k=1,nc_dim(1)
      WRITE(cfld,'("fvol",i1," ")') k
      errstat = nf_inq_varid (ncid,cfld,varid)
      errstat = nf_get_var_int2(ncid,varid,brdf_ch)
      errstat = nf_get_att_real(ncid,varid,'scale_factor',atlas%sfac_fvol(k))
      errstat = nf_get_att_real(ncid,varid,'add_offset',atlas%offs_fvol(k))
      atlas%brdfmodis(2,:,k) = brdf_ch
    ENDDO

    ! Retrieve database of 7 MODIS hinge-point FGEO values
    Do k=1,nc_dim(1)
      WRITE(cfld,'("fgeo",i1," ")') k
      errstat = nf_inq_varid (ncid,cfld,varid)
      errstat = nf_get_var_int2(ncid,varid,brdf_ch)
      errstat = nf_get_att_real(ncid,varid,'scale_factor',atlas%sfac_fgeo(k))
      errstat = nf_get_att_real(ncid,varid,'add_offset',atlas%offs_fgeo(k))
      atlas%brdfmodis(3,:,k) = brdf_ch
    ENDDO

    DEALLOCATE(brdf_ch)

    errstat = nf_close(ncid)

    CATCH
  END SUBROUTINE rttov_brdf_read_brdf
#else
  !> Read BRDF atlas data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   BRDF atlas data structure
  SUBROUTINE rttov_brdf_read_brdf(err, fn, atlas)
    ! Stub
    INTEGER(KIND=jpim),    INTENT(OUT)   :: err
    CHARACTER(LEN=*),      INTENT(IN)    :: fn
    TYPE(brdf_atlas_data), INTENT(INOUT) :: atlas
    err = 0
  END SUBROUTINE rttov_brdf_read_brdf
#endif
#endif

#ifdef _RTTOV_HDF
  !> Read BRDF atlas eigenvector file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   BRDF atlas data structure
  SUBROUTINE rttov_brdf_read_labeigvects(err, fn, atlas)

    ! Description:
    ! read the eigenvectors of the 126 selected laboratory spectra
    !(created by E borbas) from the HDF5 file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.1      03/06/2014 original code J. Vidot
    !

    INTEGER(KIND=jpim),    INTENT(OUT)   :: err
    CHARACTER(LEN=*),      INTENT(IN)    :: fn
    TYPE(brdf_atlas_data), INTENT(INOUT) :: atlas

    REAL(KIND=jprm), POINTER :: pcu_x(:,:)

#include "rttov_hdf_load.interface"

    TRY

    CALL RTTOV_HDF_LOAD(err, fn, "/BRDF", SNAME="PC_SCORES_LAND", PRM2=pcu_x)
    THROWM(ERR.NE.0,"Cannot load PC scores from "//TRIM(FN))

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

    CALL RTTOV_HDF_LOAD(err, fn, "/BRDF", SNAME="PC_SCORES_SNOW", PRM2=pcu_x)
    THROWM(ERR.NE.0,"Cannot load PC scores from "//TRIM(FN))

    atlas%pcu_snow = pcu_x(1:numpcs,1:numwave)
    DEALLOCATE(pcu_x)
    NULLIFY(pcu_x)

    CATCH
  END SUBROUTINE rttov_brdf_read_labeigvects
#else
#ifdef _RTTOV_NETCDF
  !> Read BRDF atlas eigenvector file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   BRDF atlas data structure
  SUBROUTINE rttov_brdf_read_labeigvects(err, fn, atlas)

    ! Description:
    ! read the eigenvectors of the 126 selected laboratory spectra
    !(created by E borbas) from the netCDF file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0     03/01/2012 original code J. Vidot
    !

    INTEGER(KIND=jpim),    INTENT(OUT)   :: err
    CHARACTER(LEN=*),      INTENT(IN)    :: fn
    TYPE(brdf_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpim) :: errstat   ! error code
    INTEGER(KIND=jpim) :: ncid, varid
    INTEGER(KIND=jpim) :: ndims, nvars, ngatts, recdim, nrecs
    INTEGER(KIND=jpim) :: nc_dim(2)

    CHARACTER(LEN=1024) :: strbuf ! string buffer for var

    TRY

    ! Open netCDF file.
    errstat = nf_open(TRIM(fn), nf_nowrite,ncid)

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

    ! This file has only one dimension: data is numwave x numwave
    IF (nc_dim(1) .NE. numwave) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in eigenvector data: numwave")
    ENDIF

    ! Read the laboratory eigenvalues into array pcu
    errstat = nf_inq_varid(ncid, 'PC_scores_land', varid)

    ! Specify integer KINDs for compatibility
    errstat = nf_get_vara_real(ncid, varid, INT((/1,1/), KIND=jpim), INT((/numpcs,numwave/), KIND=jpim), atlas%pcu)

    ! Read the snow eigenvalues into array pcu
    errstat = nf_inq_varid(ncid, 'PC_scores_snow', varid)

    ! Specify integer KINDs for compatibility
    errstat = nf_get_vara_real(ncid, varid, INT((/1,1/), KIND=jpim), INT((/numpcs,numwave/), KIND=jpim), atlas%pcu_snow)

    errstat = nf_close(ncid)

    CATCH
  END SUBROUTINE rttov_brdf_read_labeigvects
#else
  !> Read BRDF atlas eigenvector file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   BRDF atlas data structure
  SUBROUTINE rttov_brdf_read_labeigvects(err, fn, atlas)
    ! Stub
    INTEGER(KIND=jpim),    INTENT(OUT)   :: err
    CHARACTER(LEN=*),      INTENT(IN)    :: fn
    TYPE(brdf_atlas_data), INTENT(INOUT) :: atlas
    err = 0
  END SUBROUTINE rttov_brdf_read_labeigvects
#endif
#endif

!------------------------------------------
! Routines for returning emissivity values
!------------------------------------------

  !> Main subroutine to return BRDFs and, optionally, directional-hemipsherical albedos
  !! for the given location, surface type, satellite and solar geometry, and
  !! wavenumbers. Note that the platform_id, sat_id, inst_id and ncoefchans
  !! arguments are only used for safety checks when the atlas data were
  !! initialised for a specific instrument.
  !! @param[out]      err                 status on exit
  !! @param[in]       verbose             flag to turn verbose output on/off
  !! @param[in]       nchs                number of channels for which to return BRDFs
  !! @param[in]       lat                 latitude for which to return BRDFs
  !! @param[in]       lon                 longitude for which to return BRDFs
  !! @param[in]       surfacetype         surface type as defined in RTTOV profile
  !! @param[in]       watertype           water type as defined in RTTOV profile
  !! @param[in]       vzangle             viewing zenith angle
  !! @param[in]       vaangle             viewing azimuth angle
  !! @param[in]       szangle             solar zenith angle
  !! @param[in]       saangle             solar azimuth angle
  !! @param[in]       instr_wavenum       channel wavenumber list for required BRDFs
  !! @param[in]       channels            channel list for required BRDFs (for single-inst init)
  !! @param[in]       platform_id         platform ID from associated rtcoef structure
  !! @param[in]       sat_id              satellite ID from associated rtcoef structure
  !! @param[in]       inst_id             instrument ID from associated rtcoef structure
  !! @param[in]       ncoefchans          number of channels in associated rtcoef structure
  !! @param[in]       atlas               BRDF atlas data structure
  !! @param[out]      instr_brdf          calculated BRDF values
  !! @param[out]      instr_brdf_flag     quality flags for returned BRDFs
  !! @param[out]      instr_bs_albedo     calculated directional-hemipsherical albedo values, optional
  SUBROUTINE rttov_visnirbrdf( &
          err,               &
          verbose,           &
          nchs,              &
          lat,               &
          lon,               &
          surfacetype,       &
          watertype,         &
          vzangle,           &
          vaangle,           &
          szangle,           &
          saangle,           &
          instr_wavenum,     &
          channels,          &
          platform_id,       &
          sat_id,            &
          inst_id,           &
          ncoefchans,        &
          atlas,             &
          instr_brdf,        &
          instr_brdf_flag,   &
          instr_bs_albedo)

    ! Description:
    ! To compute BRDF for a given location, frequency and geometry
    ! from the 0.1 degree resolution MODIS BRDF kernel parameters data
    ! (at 7 hinge points)
    ! and labratory measurements using principal component analyses
    !
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0     03/01/2012 original code J. Vidot
    !

    USE rttov_const, ONLY : pi, surftype_land, surftype_sea, watertype_fresh_water !need to be verify

    INTEGER(KIND=jpim),    INTENT(OUT)           :: err
    LOGICAL(KIND=jplm),    INTENT(IN)            :: verbose
    INTEGER(KIND=jpim),    INTENT(IN)            :: nchs
    REAL(KIND=jprb),       INTENT(IN)            :: lat, lon
    INTEGER(KIND=jpim),    INTENT(IN)            :: surfacetype, watertype
    REAL(KIND=jprb),       INTENT(IN)            :: vzangle, vaangle, szangle, saangle
    REAL(KIND=jprb),       INTENT(IN)            :: instr_wavenum(nchs)
    INTEGER(KIND=jpim),    INTENT(IN)            :: channels(nchs)
    INTEGER(KIND=jpim),    INTENT(IN)            :: platform_id
    INTEGER(KIND=jpim),    INTENT(IN)            :: sat_id
    INTEGER(KIND=jpim),    INTENT(IN)            :: inst_id
    INTEGER(KIND=jpim),    INTENT(IN)            :: ncoefchans
    TYPE(brdf_atlas_data), INTENT(IN)            :: atlas
    REAL(KIND=jprb),       INTENT(OUT)           :: instr_brdf(nchs)
    INTEGER(KIND=jpim),    INTENT(OUT)           :: instr_brdf_flag
    REAL(KIND=jprb),       INTENT(OUT), OPTIONAL :: instr_bs_albedo(nchs)

    REAL(KIND=jprb) :: modisfiso(hngpnts), modisfvol(hngpnts), modisfgeo(hngpnts)
    REAL(KIND=jprb) :: modisbrdf(hngpnts), modisalbedo(hngpnts)
    REAL(KIND=jprb) :: hsrbrdf(numwave), albedo(numwave)

    REAL(KIND=jprb) :: long
    INTEGER(KIND=jpim) :: gridy, gridx, ilat, ilon

    TRY

    instr_brdf(:) = hkod
    instr_brdf_flag = hkod

    IF (atlas%single_inst) THEN
      IF (atlas%platform_id /= platform_id .OR. &
          atlas%inst_id /= inst_id .OR. &
          atlas%ncoefchans /= ncoefchans) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0,'BRDF atlas called for different coefs to that with which it was initialised')
        RETURN
      ENDIF
      IF (atlas%sat_id /= sat_id) THEN
        IF (verbose) &
          WARN('WARNING: BRDF atlas called for instrument with different sat_id to that with which it was initialised')
      ENDIF
    ENDIF

    IF (PRESENT(instr_bs_albedo)) THEN
      instr_bs_albedo(:) = hkod
    ENDIF

    IF (surfacetype == surftype_land) THEN

    !------------------------------------------------------------
    ! find the closest grid point from the MODIS BRDF kernel parameters database
    !------------------------------------------------------------

      long = MODULO(lon, 360.0_jprb)
      IF (long >= 180.0) THEN
        long = long - 360.0_jprb
      ENDIF

      ilat = NINT(lat * 1000._jprb, KIND=jpim)
      ilon = NINT(long * 1000._jprb, KIND=jpim)

      gridy = NINT(ABS(brdfmodis_ygrid1-ilat) * 1._jprb / brdfmodis_gridres, KIND=jpim) + 1_jpim
      gridx = NINT(ABS(brdfmodis_xgrid1-ilon) * 1._jprb / brdfmodis_gridres, KIND=jpim) + 1_jpim
      gridy = MAX(MIN(gridy, atlas%nb_lats), 1)
      gridx = MAX(MIN(gridx, atlas%nb_lons), 1)

      instr_brdf_flag = atlas%brdfmodis_flag(gridy,gridx)

      IF (instr_brdf_flag == 6) RETURN  ! No BRDF data for this location

      !------------------------------
      ! check if it is a land pixel
      !------------------------------

      IF (instr_brdf_flag > 0) THEN

        ! Find the brdf

        IF (atlas%brdfmodis_lut(gridy,gridx) > 0) THEN
          modisfiso(:) = REAL(atlas%brdfmodis(1,atlas%brdfmodis_lut(gridy,gridx),:), KIND=jprb) * &
                         atlas%sfac_fiso(:) + atlas%offs_fiso(:)
          modisfvol(:) = REAL(atlas%brdfmodis(2,atlas%brdfmodis_lut(gridy,gridx),:), KIND=jprb) * &
                         atlas%sfac_fvol(:) + atlas%offs_fvol(:)
          modisfgeo(:) = REAL(atlas%brdfmodis(3,atlas%brdfmodis_lut(gridy,gridx),:), KIND=jprb) * &
                         atlas%sfac_fgeo(:) + atlas%offs_fgeo(:)
        ELSE
          modisfiso(:) = 0.0_jprb
          modisfvol(:) = 0.0_jprb
          modisfgeo(:) = 0.0_jprb
        ENDIF

        CALL rttov_recon_modisbrdf(&
            vzangle,                            &! in
            vaangle,                            &! in
            szangle,                            &! in
            saangle,                            &! in
            modisfiso,                          &! in
            modisfvol,                          &! in
            modisfgeo,                          &! in
            modisbrdf,                          &! out
            modisalbedo,                        &! out
            PRESENT(instr_bs_albedo))            ! in

        IF (atlas%single_inst) THEN

          !-----------------------------------------------
          ! compute the brdf at the instrument wavenumbers
          !-----------------------------------------------

          IF (instr_brdf_flag == 4_jpim) THEN
            CALL rttov_brdf_recon_brdf_snow( &
                atlas,                  &! in
                modisbrdf,              &! in
                channels,               &! in
                instr_brdf)              ! out
          ELSE
            CALL rttov_brdf_recon_brdf( &
                atlas,                  &! in
                modisbrdf,              &! in
                channels,               &! in
                instr_brdf)              ! out
          ENDIF

          WHERE (instr_wavenum < hsr_wavenum(1) .OR. &
                 instr_wavenum > hsr_wavenum(numwave))
            instr_brdf = hkod
          ENDWHERE

          IF (PRESENT(instr_bs_albedo)) THEN

            IF (instr_brdf_flag == 4_jpim) THEN
              instr_bs_albedo = hkod
            ELSE
              CALL rttov_brdf_recon_brdf( &
                  atlas,                  &! in
                  modisalbedo,            &! in
                  channels,               &! in
                  instr_bs_albedo)         ! out

              WHERE (instr_wavenum < hsr_wavenum(1) .OR. &
                     instr_wavenum > hsr_wavenum(numwave))
                instr_bs_albedo = hkod
              ENDWHERE
            ENDIF

          ENDIF

        ELSE

          !------------------------------------------------------------------------------------------
          ! compute the hsr brdf spectra at 2101 wavenumber points from the 7 MODIS BRDF hinge points
          !------------------------------------------------------------------------------------------

          IF (instr_brdf_flag == 4_jpim) THEN
            CALL rttov_brdf_recon_hsrbrdf_snow( &
                atlas,                     &! in
                modisbrdf,                 &! in
                hsrbrdf)                    ! out
          ELSE
            CALL rttov_brdf_recon_hsrbrdf( &
                atlas,                     &! in
                modisbrdf,                 &! in
                hsrbrdf)                    ! out
          ENDIF

          IF (PRESENT(instr_bs_albedo)) THEN
            CALL rttov_brdf_recon_hsrbrdf( &
                atlas,                     &! in
                modisalbedo,               &! in
                albedo)                     ! out
          ENDIF

          !--------------------------------------------------------------------------------
          ! create instrument specific brdf by finding the closest wavenumber value
          !--------------------------------------------------------------------------------

          CALL rttov_brdf_select_wavenum( &
              hsrbrdf,                    &! in
              nchs,                       &! in
              instr_wavenum(1:nchs),      &! in
              instr_brdf)                  ! out

          !--------------------------------------------------------------------------------
          ! create instrument specific directional-hemispherical albedo by finding the closest wavenumber value
          !--------------------------------------------------------------------------------

          If (PRESENT(instr_bs_albedo)) THEN

            If (instr_brdf_flag == 4_jpim) THEN
              instr_bs_albedo = hkod
            ELSE
              CALL rttov_brdf_select_wavenum( &
                  albedo,                     &! in
                  nchs,                       &! in
                  instr_wavenum(1:nchs),      &! in
                  instr_bs_albedo)             ! out
            ENDIF

          ENDIF

        ENDIF ! single_inst

      ENDIF  ! brdf flag = 0

    ENDIF ! surftype_land

    !--------------------------------------------------------------------------------
    ! if water surface, then used default spectra for open ocean or coastal waters
    !--------------------------------------------------------------------------------

    IF (surfacetype == surftype_sea) THEN

      IF (atlas%single_inst) THEN

        IF (watertype == watertype_fresh_water) THEN
          instr_brdf = atlas%coastal_waters_ref_int(channels(:))
        ELSE
          instr_brdf = atlas%ocean_waters_ref_int(channels(:))
        ENDIF

        WHERE (instr_wavenum < hsr_wavenum(1) .OR. &
               instr_wavenum > hsr_wavenum(numwave))
          instr_brdf = hkod
        ENDWHERE

      ELSE

        IF (watertype == watertype_fresh_water) THEN
          hsrbrdf = coastal_waters_ref
        ELSE
          hsrbrdf = ocean_waters_ref
        ENDIF

        CALL rttov_brdf_select_wavenum( &
              hsrbrdf,                  &! in
              nchs,                     &! in
              instr_wavenum(1:nchs),    &! in
              instr_brdf)                ! out

      ENDIF

    ENDIF

    ! Cap the final BRDFs here for consistency between single-
    ! and multi-instrument initialisation.
    WHERE (instr_brdf > hkod)
      ! The atlas gives BRFs, but we want the output to be a BRDF so normalise by pi
      instr_brdf = instr_brdf / pi

! Physical BRDFs may exceed 1. (e.g. Nicodemus et al 1977) so do not apply any upper cap
      instr_brdf = MAX(instr_brdf, 0._jprb)
!       instr_brdf = MIN(instr_brdf, 1._jprb)
    ENDWHERE

    IF (PRESENT(instr_bs_albedo)) THEN
      WHERE (instr_bs_albedo > hkod)
        instr_bs_albedo = MAX(instr_bs_albedo, 0._jprb)
        instr_bs_albedo = MIN(instr_bs_albedo, 1._jprb)
      ENDWHERE
    ENDIF

    CATCH
  END SUBROUTINE rttov_visnirbrdf


  !> Compute BRDF and directional-hemipsherical albedo from the semi-empirical model of
  !! Ross-Li from the three kernel parameters following the paper of Lucht et 
  !! al (2000) Eq 37 for BRDF and Eq. 46 for albedo.
  !! @param[in]      vzangle         viewing zenith angle
  !! @param[in]      vaangle         viewing azimuth angle
  !! @param[in]      szangle         solar zenith angle
  !! @param[in]      saangle         solar azimuth angle
  !! @param[in]      modisfiso       BRDF model coefficients
  !! @param[in]      modisfvol       BRDF model coefficients
  !! @param[in]      modisfgeo       BRDF model coefficients
  !! @param[out]     modisbrdf       calculated BRDF
  !! @param[out]     modisalbedo     calculated directional-hemipsherical albedo
  !! @param[in]      calc_albedo     if TRUE also calculate albedo
  SUBROUTINE rttov_recon_modisbrdf( &
            vzangle,                            &
            vaangle,                            &
            szangle,                            &
            saangle,                            &
            modisfiso,                          &
            modisfvol,                          &
            modisfgeo,                          &
            modisbrdf,                          &
            modisalbedo,                        &
            calc_albedo)

    ! Description:
    !   Compute the brdf and the directional-hemispherical albedo from the semi
    !   empirical model of Ross-Li from the three kernel parameters
    !   following the paper of Lucht et al.(2000) Eq 37 for BRDF
    !   and Eq. 46 for albedo.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0     03/01/2012 original code J. Vidot

    USE rttov_const, ONLY : deg2rad, pi

    REAL(KIND=jprb),    INTENT(IN)  :: vzangle, vaangle, szangle, saangle
    REAL(KIND=jprb),    INTENT(IN)  :: modisfiso(hngpnts)
    REAL(KIND=jprb),    INTENT(IN)  :: modisfvol(hngpnts)
    REAL(KIND=jprb),    INTENT(IN)  :: modisfgeo(hngpnts)
    REAL(KIND=jprb),    INTENT(OUT) :: modisbrdf(hngpnts)
    REAL(KIND=jprb),    INTENT(OUT) :: modisalbedo(hngpnts)
    LOGICAL(KIND=jplm), INTENT(IN)  :: calc_albedo

    INTEGER(KIND=jpim) :: j
    REAL(KIND=jprb)    :: cossz, cosvz, tansz, tanvz, cosva_sa
    REAL(KIND=jprb)    :: cos_ksi,sin_ksi,ksi
    REAL(KIND=jprb)    :: kernel_vol,kernel_geo
    REAL(KIND=jprb)    :: delta_kgeo,sec_szangle,sec_vzangle,sec,cos_t,sin_t,t_kgeo

    REAL(KIND=jprb)    :: bs_alb_param(6)

    !--------------------------------------------------------------------------------
    ! see Lucht et al., 2000, IEEE vol 38, no 2, pp 977-998 for equations and
    ! for coefficients to compute directional-hemispherical (or black-sky) albedo in Table 1
    !--------------------------------------------------------------------------------

    DATA bs_alb_param / &
    -0.007574_jprb, -0.070987_jprb, 0.307588_jprb, -1.284909_jprb, -0.166314_jprb, 0.041840_jprb /

    cossz = COS(szangle*deg2rad)
    cosvz = COS(vzangle*deg2rad)
    tansz = TAN(szangle*deg2rad)
    tanvz = TAN(vzangle*deg2rad)
    cosva_sa = COS((vaangle-saangle)*deg2rad)

    cos_ksi = cossz*cosvz + SIN(szangle*deg2rad)*SIN(vzangle*deg2rad)*cosva_sa
    sin_ksi = SQRT(1._jprb-cos_ksi**2_jpim)
    ksi = ACOS(cos_ksi)
    sec_szangle = 1._jprb/cossz
    sec_vzangle = 1._jprb/cosvz
    sec = sec_szangle + sec_vzangle
    kernel_vol = (((pi/2._jprb - ksi)*cos_ksi + sin_ksi) / &
        (cossz + cosvz)) - PI/4._jprb  !! Eq 38

    delta_kgeo = SQRT(tansz**2_jpim + tanvz**2_jpim - 2._jprb * tansz*tanvz*cosva_sa)

    cos_t = (2._jprb*SQRT(delta_kgeo**2_jpim + &
        (tansz*tanvz)**2_jpim * (1._jprb - cosva_sa**2_jpim))/sec)

    IF (cos_t < -1._jprb) cos_t = -1._jprb
    IF (cos_t > 1._jprb) cos_t = 1._jprb

    sin_t = SQRT(1._jprb - cos_t**2_jpim)
    t_kgeo = ACOS(cos_t)

    kernel_geo =  ((t_kgeo - sin_t*cos_t)*sec/pi - sec) + (1._jprb + cos_ksi)*sec_vzangle*sec_szangle/2._jprb !!eq 39

    DO j = 1, hngpnts
      modisbrdf(j) = modisfiso(j) + modisfvol(j) * kernel_vol + modisfgeo(j) * kernel_geo !! Eq 37
    ENDDO

    IF (calc_albedo) THEN
      DO j = 1, hngpnts
        modisalbedo(j) = modisfiso(j) + &
          modisfvol(j) * ( bs_alb_param(1) + bs_alb_param(2) * ((szangle*deg2rad)**2_jpim) + &
          bs_alb_param(3) * ((szangle*deg2rad)**3_jpim)) + &
          modisfgeo(j) * ( bs_alb_param(4) + bs_alb_param(5) * ((szangle*deg2rad)**2_jpim) + &
          bs_alb_param(6) * ((szangle*deg2rad)**3_jpim))
      ENDDO
    ENDIF

  END SUBROUTINE rttov_recon_modisbrdf


  !> Interpolate PCs onto instrument channel wavenumbers. This is called
  !! during initialisation when single-instrument init is selected.
  !! @param[in]      instr_wavenum   instrument channel wavenumbers
  !! @param[in,out]  atlas           BRDF atlas data structure
  SUBROUTINE rttov_brdf_hsr_interp(instr_wavenum, atlas)

    REAL(KIND=jprb),       INTENT(IN)    :: instr_wavenum(:)
    TYPE(brdf_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpim) :: j, k, nchs

    REAL(KIND=jprb) :: dist(numwave)
    REAL(KIND=jprb) :: mindist
    INTEGER(KIND=jpim) :: ind_mindist

    REAL(KIND=jprb) :: dwvnum1,dwvnum2,dwvsum
    REAL(KIND=jprb) :: pcu1(numpcs),pcu2(numpcs),pcm1,pcm2
    REAL(KIND=jprb) :: pcu1_snow(numpcs),pcu2_snow(numpcs),pcm1_snow,pcm2_snow
    REAL(KIND=jprb) :: coastal_waters_ref1, coastal_waters_ref2
    REAL(KIND=jprb) :: ocean_waters_ref1, ocean_waters_ref2

    !---------------------------------------------------------------
    ! finding the closest frequency from the hsr emissivity spectr
    !--------------------------------------------------------------

    nchs = SIZE(instr_wavenum)

    DO j = 1, nchs

      IF (instr_wavenum(j) <= hsr_wavenum(1)) THEN

        atlas%pcu_int(:,j)              = REAL(atlas%pcu(:,1), KIND=jprb)
        atlas%pcm_int(j)                = pcm(1)
        atlas%pcu_int_snow(:,j)         = REAL(atlas%pcu_snow(:,1), KIND=jprb)
        atlas%pcm_int_snow(j)           = pcm_snow(1)
        atlas%coastal_waters_ref_int(j) = coastal_waters_ref(1)
        atlas%ocean_waters_ref_int(j)   = ocean_waters_ref(1)

      ELSE IF (instr_wavenum(j) >= hsr_wavenum(numwave)) THEN

        atlas%pcu_int(:,j)              = REAL(atlas%pcu(:,numwave), KIND=jprb)
        atlas%pcm_int(j)                = pcm(numwave)
        atlas%pcu_int_snow(:,j)         = REAL(atlas%pcu_snow(:,numwave), KIND=jprb)
        atlas%pcm_int_snow(j)           = pcm_snow(numwave)
        atlas%coastal_waters_ref_int(j) = coastal_waters_ref(numwave)
        atlas%ocean_waters_ref_int(j)   = ocean_waters_ref(numwave)

      ELSE ! within wavenumber compute range

        mindist = 100._jprb
        ind_mindist = 100000_jpim

        DO k = 1, numwave

          ! calculate distances between the instr freq and hsr emissivities
          dist(k) = ABS(instr_wavenum(j) - hsr_wavenum(k))

          ! finding the closest frequency from the hsr emissivity
          IF (dist(k) <=  mindist) THEN

            mindist = dist(k)
            ind_mindist = k

          ENDIF
        ENDDO

  !--------------------------------
  ! Interpolate values
  !--------------------------------

  ! Bilinear mean of the two closest spectral points

        k = 1
        IF (instr_wavenum(j) <= hsr_wavenum(ind_mindist)) k = -1

        dwvnum1 = dist( ind_mindist )
        dwvnum2 = dist( ind_mindist + k )
        dwvsum = dwvnum1 + dwvnum2

        pcu1(:)   = dwvnum1 * REAL(atlas%pcu(:,ind_mindist+k), KIND=jprb)
        pcu2(:)   = dwvnum2 * REAL(atlas%pcu(:,ind_mindist), KIND=jprb)
        pcm1      = dwvnum1 * pcm(ind_mindist+k)
        pcm2      = dwvnum2 * pcm(ind_mindist)
        pcu1_snow(:)   = dwvnum1 * REAL(atlas%pcu_snow(:,ind_mindist+k), KIND=jprb)
        pcu2_snow(:)   = dwvnum2 * REAL(atlas%pcu_snow(:,ind_mindist), KIND=jprb)
        pcm1_snow      = dwvnum1 * pcm_snow(ind_mindist+k)
        pcm2_snow      = dwvnum2 * pcm_snow(ind_mindist)
        coastal_waters_ref1 = dwvnum1 * coastal_waters_ref(ind_mindist+k)
        coastal_waters_ref2 = dwvnum2 * coastal_waters_ref(ind_mindist)
        ocean_waters_ref1   = dwvnum1 * ocean_waters_ref(ind_mindist+k)
        ocean_waters_ref2   = dwvnum2 * ocean_waters_ref(ind_mindist)

        atlas%pcu_int(:,j)    = (pcu1(:) + pcu2(:)) / dwvsum
        atlas%pcm_int(j)      = (pcm1 + pcm2) / dwvsum
        atlas%pcu_int_snow(:,j)    = (pcu1_snow(:) + pcu2_snow(:)) / dwvsum
        atlas%pcm_int_snow(j)      = (pcm1_snow + pcm2_snow) / dwvsum
        atlas%coastal_waters_ref_int(j) = (coastal_waters_ref1 + coastal_waters_ref2) / dwvsum
        atlas%ocean_waters_ref_int(j)   = (ocean_waters_ref1 + ocean_waters_ref2) / dwvsum

      ENDIF

    ENDDO

  END SUBROUTINE rttov_brdf_hsr_interp

  !> Reconstruct BRDFs at instrument wavenumbers from interpolated PCs (used
  !! when atlas is initialised for a specific instrument)
  !! @param[in]   atlas           BRDF atlas data structure
  !! @param[in]   modisbrdf       MODIS kernel parameters
  !! @param[in]   channels        channels for which BRDFs are required
  !! @param[out]  brdf            reconstructed BRDFs
  SUBROUTINE rttov_brdf_recon_brdf( &
          atlas,                    &
          modisbrdf,                &
          channels,                 &
          brdf)

    ! Description:
    ! Reconstruct BRDFs at instrument wavenumbers from interpolated PCs
    !
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0     03/01/2012 original code J. Vidot
    !

    TYPE(brdf_atlas_data), INTENT(IN)  :: atlas
    REAL(KIND=jprb),       INTENT(IN)  :: modisbrdf(hngpnts)
    INTEGER(KIND=jpim),    INTENT(IN)  :: channels(:)
    REAL(KIND=jprb),       INTENT(OUT) :: brdf(SIZE(channels))

    INTEGER(KIND=jpim) :: k, nchn

    REAL(KIND=jprb) :: col(hngpnts)
    REAL(KIND=jprb) :: coef(numpcs)

    !-------------------------
    ! calculate the regcoef
    !-------------------------

    IF ( ALL(modisbrdf(:) > 0._jprb) ) THEN

      col(:) = modisbrdf(:) * 100._jprb - pcm_modres(:)

      Do k = 1, numpcs
        coef(k) = SUM(col(:) * atlas%D(:,k))
      ENDDO

      !-----------------------------------
      ! apply regcoef to get the hsr dataset
      !-----------------------------------

      nchn = SIZE(brdf)
      DO k = 1, nchn

        brdf(k) = SUM(coef(:) * REAL(atlas%pcu_int(:,channels(k)), KIND=jprb)) + atlas%pcm_int(channels(k))
        brdf(k) = brdf(k) / 100._jprb

        ! This is done in the calling subroutine
!         brdf(k) = MAX(brdf(k), 0._jprb)
!         brdf(k) = MIN(brdf(k), 1._jprb)

      ENDDO

    ELSE

      brdf = hkod

    ENDIF

  END SUBROUTINE rttov_brdf_recon_brdf


  !> Reconstruct snow BRDFs at instrument wavenumbers from interpolated PCs
  !! (used when atlas is initialised for a specific instrument)
  !! @param[in]   atlas           BRDF atlas data structure
  !! @param[in]   modisbrdf       MODIS kernel parameters
  !! @param[in]   channels        channels for which BRDFs are required
  !! @param[out]  brdf            reconstructed snow BRDFs
  SUBROUTINE rttov_brdf_recon_brdf_snow( &
          atlas,                    &! in
          modisbrdf,                &! in
          channels,                 &! in
          brdf)                      ! out

    ! Description:
    ! Reconstruct BRDFs at instrument wavenumbers from interpolated PCs
    !
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.1     03/06/2014 original code J. Vidot

    TYPE(brdf_atlas_data), INTENT(IN)  :: atlas
    REAL(KIND=jprb),       INTENT(IN)  :: modisbrdf(hngpnts)
    INTEGER(KIND=jpim),    INTENT(IN)  :: channels(:)
    REAL(KIND=jprb),       INTENT(OUT) :: brdf(SIZE(channels))

    INTEGER(KIND=jpim) :: k, nchn

    REAL(KIND=jprb) :: col(hngpnts)
    REAL(KIND=jprb) :: coef(numpcs)

    !-------------------------
    ! calculate the regcoef
    !-------------------------

    IF (ALL(modisbrdf(:) > 0._jprb)) THEN

      col(:) = modisbrdf(:) * 100._jprb - pcm_modres_snow(:)

      DO k = 1, numpcs
        coef(k) = SUM(col(:) * atlas%D_snow(:,k))
      ENDDO

      !-----------------------------------
      ! apply regcoef to get the hsr dataset
      !-----------------------------------

      nchn = SIZE(brdf)
      DO k = 1, nchn

        brdf(k) = SUM(coef(:) * REAL(atlas%pcu_int_snow(:,channels(k)), KIND=jprb)) + atlas%pcm_int_snow(channels(k))
        brdf(k) = brdf(k) / 100._jprb

      ENDDO

    ELSE

      brdf = hkod

    ENDIF

  END SUBROUTINE rttov_brdf_recon_brdf_snow


  !> Reconstruct high-spectral-resolution BRDFs from PCs
  !! @param[in]   atlas           BRDF atlas data structure
  !! @param[in]   modisbrdf       MODIS kernel parameters
  !! @param[out]  hsrbrdf         reconstructed BRDFs
  SUBROUTINE rttov_brdf_recon_hsrbrdf( &
          atlas,                       &
          modisbrdf,                   &
          hsrbrdf)

    ! Description:
    ! To creates high spectra resolution brdf at 2101 wavenumbers
    ! from the MODIS BRDF kernel parameters data (at 7 hinge points)
    ! and labratory measurements using principal component analyses
    !
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0     03/01/2012 original code J. Vidot

    TYPE(brdf_atlas_data), INTENT(IN)  :: atlas
    REAL(KIND=jprb),       INTENT(IN)  :: modisbrdf(hngpnts)
    REAL(KIND=jprb),       INTENT(OUT) :: hsrbrdf(numwave)

    INTEGER(KIND=jpim) :: k

    REAL(KIND=jprb) :: col(hngpnts)
    REAL(KIND=jprb) :: coef(numpcs)

    !-------------------------
    ! calculate the regcoef
    !-------------------------

    IF (ALL(modisbrdf(:) > 0._jprb)) THEN

      col(:) = modisbrdf(:) * 100._jprb - pcm_modres(:)

      DO k = 1, numpcs
        coef(k) = SUM(col(:) * atlas%D(:,k))
      ENDDO

      !-----------------------------------
      ! apply regcoef to get the hsr dataset
      !-----------------------------------

      DO k = 1, numwave

        hsrbrdf(k) = SUM(coef(:) * REAL(atlas%pcu(:,k), KIND=jprb)) + pcm(k)
        hsrbrdf(k) = hsrbrdf(k) / 100._jprb

      ENDDO

    ELSE

      hsrbrdf = hkod

    ENDIF

  END SUBROUTINE rttov_brdf_recon_hsrbrdf


  !> Reconstruct high-spectral-resolution BRDFs for snow from PCs
  !! @param[in]   atlas           BRDF atlas data structure
  !! @param[in]   modisbrdf       MODIS kernel parameters
  !! @param[out]  hsrbrdf         reconstructed snow BRDFs
  SUBROUTINE rttov_brdf_recon_hsrbrdf_snow( &
          atlas,                       &
          modisbrdf,                   &
          hsrbrdf)

    ! Description:
    ! To creates high spectra resolution brdf at 2101 wavenumbers
    ! from the MODIS BRDF kernel parameters data (at 7 hinge points)
    ! and labratory measurements using principal component analyses
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0     03/06/2014 original code J. Vidot

    TYPE(brdf_atlas_data), INTENT(IN)  :: atlas
    REAL(KIND=jprb),       INTENT(IN)  :: modisbrdf(hngpnts)
    REAL(KIND=jprb),       INTENT(OUT) :: hsrbrdf(numwave)

    INTEGER(KIND=jpim) :: k

    REAL(KIND=jprb) :: col(hngpnts)
    REAL(KIND=jprb) :: coef(numpcs)

    !-------------------------
    ! calculate the regcoef
    !-------------------------

    IF (ALL(modisbrdf(:) > 0._jprb)) THEN

      col(:) = modisbrdf(:) * 100._jprb - pcm_modres_snow(:)

      DO k = 1, numpcs
        coef(k) = SUM(col(:) * atlas%D_snow(:,k))
      ENDDO

      !-----------------------------------
      ! apply regcoef to get the hsr dataset
      !-----------------------------------

      DO k = 1, numwave

        hsrbrdf(k) = SUM(coef(:) * REAL(atlas%pcu_snow(:,k), KIND=jprb)) + pcm_snow(k)
        hsrbrdf(k) = hsrbrdf(k) / 100._jprb

      ENDDO

    ELSE

      hsrbrdf = hkod

    ENDIF

  END SUBROUTINE rttov_brdf_recon_hsrbrdf_snow


  !> Linear interpolation of high-spectral-resolution BRDF data onto
  !! instrument channel wavenumbers
  !! @param[in]   hsrbrdf         high-spectral-resolution BRDF data
  !! @param[in]   nchs            number of instrument channels
  !! @param[in]   instr_wavenum   channel central wavenumbers
  !! @param[out]  instr_brdf      interpolated BRDF values
  SUBROUTINE rttov_brdf_select_wavenum( &
          hsrbrdf,                      &
          nchs,                         &
          instr_wavenum,                &
          instr_brdf)

    ! Description:
    ! subroutine to find the closest wavenumber from the HSR brdf spectra
    ! for the instrument frequency and assign the instrument brdf y by choosing the
    ! closest spectral point value or bilinear interpolating  between the two
    ! closest spectral point values
    !
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0     03/01/2012 original code J. Vidot
    !

    INTEGER(KIND=jpim), INTENT(IN)  :: nchs
    REAL(KIND=jprb),    INTENT(IN)  :: hsrbrdf(numwave)
    REAL(KIND=jprb),    INTENT(IN)  :: instr_wavenum(nchs)
    REAL(KIND=jprb),    INTENT(OUT) :: instr_brdf(nchs)

    INTEGER(KIND=jpim) :: j, k

    REAL(KIND=jprb) :: dist(numwave)
    REAL(KIND=jprb) :: mindist
    INTEGER(KIND=jpim) :: ind_mindist

    REAL(KIND=jprb) :: dwvnum1,dwvnum2,dwvsum
    REAL(KIND=jprb) :: hsrbrdf1,hsrbrdf2
    LOGICAL(KIND=jplm) :: lcpu_brdf

    ! initialize instrument brdf

    !---------------------------------------------------------------
    ! finding the closest frequency from the hsr brdf spectr
    !--------------------------------------------------------------

    lcpu_brdf = .NOT. ALL(hsrbrdf == hsrbrdf(1))

    IF (lcpu_brdf) THEN
      instr_brdf(:) = hkod
      DO j = 1, nchs

        IF (instr_wavenum(j) >= hsr_wavenum(1) .AND. &
            instr_wavenum(j) <= hsr_wavenum(numwave)) THEN

          ! within wavenumber compute range

          mindist = 100._jprb
          ind_mindist = 100000_jpim

          DO k = 1, numwave

            ! calucalte distances between the instr freq end hsr brdf
            dist(k) = ABS(instr_wavenum(j) - hsr_wavenum(k))

            ! finding the closest frequency from the hsr brdf
            IF (dist(k) <=  mindist) THEN

              mindist = dist(k)
              ind_mindist = k

            ENDIF
          ENDDO

    !--------------------------------
    ! assign instrument brdf
    !--------------------------------
    !  closest spectral point
    !                       instr_brdf(j)=hsrbrdf(ind_mindist)

    ! or bilinear mean of the two closest spectral points

          k = 1
          IF (instr_wavenum(j) <= hsr_wavenum(ind_mindist) ) k=-1

          dwvnum1 = dist(ind_mindist)
          dwvnum2 = dist(ind_mindist + k)
          dwvsum = dwvnum1 + dwvnum2

          IF (lcpu_brdf) THEN
            hsrbrdf1 = dwvnum1 * hsrbrdf(ind_mindist + k)
            hsrbrdf2 = dwvnum2 * hsrbrdf(ind_mindist)
            instr_brdf(j) = (hsrbrdf1 + hsrbrdf2) / dwvsum
          ELSE
            instr_brdf(j) = hsrbrdf(1)
          ENDIF

        ENDIF    !==  (instr_wavenum(j) <= hsr_wavenum(1))

      ENDDO

    ELSE  ! all logical computes (lcpu_brdf) are false
      instr_brdf(:)=hsrbrdf(1)
    ENDIF

  END SUBROUTINE rttov_brdf_select_wavenum

!------------------------------------
! Routine to deallocate atlas arrays
!------------------------------------

  !> Deallocate data in BRDF atlas data structure
  !! @param[in,out]   atlas   BRDF atlas data structure to deallocate
  SUBROUTINE rttov_visnirbrdf_close_atlas(atlas)
    TYPE(brdf_atlas_data), INTENT(INOUT) :: atlas

    IF (ASSOCIATED(atlas%brdfmodis_flag))         DEALLOCATE(atlas%brdfmodis_flag)
    IF (ASSOCIATED(atlas%brdfmodis_lut))          DEALLOCATE(atlas%brdfmodis_lut)
    IF (ASSOCIATED(atlas%brdfmodis))              DEALLOCATE(atlas%brdfmodis)
    IF (ASSOCIATED(atlas%sfac_fiso))              DEALLOCATE(atlas%sfac_fiso)
    IF (ASSOCIATED(atlas%offs_fiso))              DEALLOCATE(atlas%offs_fiso)
    IF (ASSOCIATED(atlas%sfac_fvol))              DEALLOCATE(atlas%sfac_fvol)
    IF (ASSOCIATED(atlas%offs_fvol))              DEALLOCATE(atlas%offs_fvol)
    IF (ASSOCIATED(atlas%sfac_fgeo))              DEALLOCATE(atlas%sfac_fgeo)
    IF (ASSOCIATED(atlas%offs_fgeo))              DEALLOCATE(atlas%offs_fgeo)
    IF (ASSOCIATED(atlas%D))                      DEALLOCATE(atlas%D)
    IF (ASSOCIATED(atlas%D_snow))                 DEALLOCATE(atlas%D_snow)
    IF (ASSOCIATED(atlas%pcu))                    DEALLOCATE(atlas%pcu)
    IF (ASSOCIATED(atlas%pcu_snow))               DEALLOCATE(atlas%pcu_snow)
    IF (ASSOCIATED(atlas%pcu_int))                DEALLOCATE(atlas%pcu_int)
    IF (ASSOCIATED(atlas%pcm_int))                DEALLOCATE(atlas%pcm_int)
    IF (ASSOCIATED(atlas%pcu_int_snow))           DEALLOCATE(atlas%pcu_int_snow)
    IF (ASSOCIATED(atlas%pcm_int_snow))           DEALLOCATE(atlas%pcm_int_snow)
    IF (ASSOCIATED(atlas%coastal_waters_ref_int)) DEALLOCATE(atlas%coastal_waters_ref_int)
    IF (ASSOCIATED(atlas%ocean_waters_ref_int))   DEALLOCATE(atlas%ocean_waters_ref_int)

    CALL rttov_visnirbrdf_nullify(atlas)

  END SUBROUTINE rttov_visnirbrdf_close_atlas

  !> Nullify pointers in BRDF atlas data structure
  !! @param[in,out]   atlas   BRDF atlas data structure to nullify
  SUBROUTINE rttov_visnirbrdf_nullify(atlas)
    TYPE(brdf_atlas_data), INTENT(INOUT) :: atlas

    NULLIFY(atlas%brdfmodis_flag)
    NULLIFY(atlas%brdfmodis_lut)
    NULLIFY(atlas%brdfmodis)
    NULLIFY(atlas%sfac_fiso)
    NULLIFY(atlas%offs_fiso)
    NULLIFY(atlas%sfac_fvol)
    NULLIFY(atlas%offs_fvol)
    NULLIFY(atlas%sfac_fgeo)
    NULLIFY(atlas%offs_fgeo)
    NULLIFY(atlas%D)
    NULLIFY(atlas%D_snow)
    NULLIFY(atlas%pcu)
    NULLIFY(atlas%pcu_snow)
    NULLIFY(atlas%pcu_int)
    NULLIFY(atlas%pcm_int)
    NULLIFY(atlas%pcu_int_snow)
    NULLIFY(atlas%pcm_int_snow)
    NULLIFY(atlas%coastal_waters_ref_int)
    NULLIFY(atlas%ocean_waters_ref_int)

  END SUBROUTINE rttov_visnirbrdf_nullify

END MODULE mod_brdf_atlas
