! Description:
!> @file
!!   Subroutines for MEASURES CAMEL IR emissivity atlas
!
!> @brief
!!   Subroutines for MEASURES CAMEL IR emissivity atlas
!!
!! @details
!!   It is intended that this atlas be used via the RTTOV interface
!!   rather than by calling these subroutines directly.
!!
!!   Borbas, Eva, Glynn Hulley, Michelle Feltz, Robert Knuteson, and
!!   Simon Hook. "The Combined ASTER MODIS Emissivity over Land (CAMEL)
!!   Part 1: Methodology and High Spectral Resolution Application."
!!   Remote Sensing 10, no. 4 (2018): 643.

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
MODULE mod_camel_atlas

  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  0.6      06/06/2016  Based on MEASURES CAMEL IR atlas code (E. Borbas, G. Hulley, B. Ruston, J. Hocking)

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
  PUBLIC :: camel_atlas_data,        &
            rttov_camel_init,        &
            rttov_camel,             &
            rttov_camel_close_atlas, &
            numpcs  ! Useful for external programs

  ! Atlas constants

  INTEGER(KIND=jpim), PARAMETER :: numpcs = 9
  INTEGER(KIND=jpim), PARAMETER :: numwave = 417

  INTEGER(KIND=jpim), PARAMETER :: db_ver_year = 2007
  CHARACTER(LEN=8),   PARAMETER :: prd_version = 'v01r02'     ! Release date: 2016-11-22
  CHARACTER(LEN=4),   PARAMETER :: prd_version_cov = 'V002'   ! Release date: 2018-11-04 (CAMEL covariance data)

  INTEGER(KIND=jpim), PARAMETER :: seaice_flag = 70           ! flag value returned for sea-ice
  REAL(KIND=jprb),    PARAMETER :: default_std = 0.05_jprb    ! default standard deviation

  INTEGER(KIND=jpim), PARAMETER :: cov_emis_gridres = 250     ! 0.25 deg
  INTEGER(KIND=jpim), PARAMETER :: cov_emis_ygrid1 = 89875    ! 89.875 deg
  INTEGER(KIND=jpim), PARAMETER :: cov_emis_xgrid1 = -179875  ! -179.875 deg

  INTEGER(KIND=jpim), PARAMETER :: igbp_gridres = 50          ! 0.05 deg
  INTEGER(KIND=jpim), PARAMETER :: igbp_ygrid1 = 89950        ! 89.95 deg
  INTEGER(KIND=jpim), PARAMETER :: igbp_xgrid1 = -179950      ! -179.95 deg

  REAL(KIND=jprb),    PARAMETER :: angcorrminzen = 5._jprb
  REAL(KIND=jprb),    PARAMETER :: angcorrterm = 85._jprb     ! Solar zenith angle of terminator

  REAL(KIND=jprb),    PARAMETER :: hkod = -999._jprb          ! Missing data


  !> Data type for CAMEL atlas data
  TYPE camel_atlas_data
    PRIVATE

    ! Flags
    LOGICAL(KIND=jplm) :: single_inst
    LOGICAL(KIND=jplm) :: std_init
    LOGICAL(KIND=jplm) :: do_ang_corr

    ! Dimensions
    INTEGER(KIND=jpim) :: nb_lats        ! Dimensions of atlas data
    INTEGER(KIND=jpim) :: nb_lons
    INTEGER(KIND=jpim) :: nmask

    INTEGER(KIND=jpim) :: cv_lats        ! Dimensions of covariance data
    INTEGER(KIND=jpim) :: cv_lons
    INTEGER(KIND=jpim) :: cv_pack

    INTEGER(KIND=jpim) :: igbp_lats      ! Dimensions of IGBP data
    INTEGER(KIND=jpim) :: igbp_lons
    INTEGER(KIND=jpim) :: nb_igbp

    INTEGER(KIND=jpim) :: camel_gridres  ! Atlas grid resolution x1000
    INTEGER(KIND=jpim) :: camel_ygrid1   ! Atlas grid reference latitude x1000
    INTEGER(KIND=jpim) :: camel_xgrid1   ! Atlas grid reference longitude x1000

!     INTEGER(KIND=jpim) :: cov_emis_gridres  ! Cov grid resolution x1000
!     INTEGER(KIND=jpim) :: cov_emis_ygrid1   ! Cov grid reference latitude x1000
!     INTEGER(KIND=jpim) :: cov_emis_xgrid1   ! Cov grid reference longitude x1000

    ! Atlas data loaded by initialisation routine
    INTEGER(KIND=jpis), POINTER :: camel_qflag(:,:)  ! dims are (nb_lons,nb_lats)
    INTEGER(KIND=jpim), POINTER :: coef_lut(:,:)     ! dims are (nb_lons,nb_lats)
    REAL(KIND=jprm),    POINTER :: pca_coef(:,:)     ! dims are (numpcs,nmask)

    INTEGER(KIND=jpis), POINTER :: pc_labvs(:)
    INTEGER(KIND=jpis), POINTER :: pc_npcs(:)

    INTEGER(KIND=jpim), POINTER :: cov_emis_lut(:,:) ! dims are (cv_lats,cv_lons)
    INTEGER(KIND=jpis), POINTER :: cov_emis(:,:)     ! dims are (cv_pack,numwave)

    INTEGER(KIND=jpit), POINTER :: igbp(:,:)         ! dims are (igbp_lats,igbp_lons)

    REAL(KIND=jprm),    POINTER :: p1d(:,:)          ! dims are (numwave,nb_igbp)
    REAL(KIND=jprm),    POINTER :: p2d(:,:)          ! dims are (numwave,nb_igbp)
    REAL(KIND=jprm),    POINTER :: p3d(:,:)          ! dims are (numwave,nb_igbp)
    REAL(KIND=jprm),    POINTER :: p1n(:,:)          ! dims are (numwave,nb_igbp)
    REAL(KIND=jprm),    POINTER :: p2n(:,:)          ! dims are (numwave,nb_igbp)
    REAL(KIND=jprm),    POINTER :: p3n(:,:)          ! dims are (numwave,nb_igbp)

    REAL(KIND=jprb),    POINTER :: pcm_hsr_v8(:), pcm_hsr_v10(:), pcm_hsr_v12(:)
    REAL(KIND=jprb),    POINTER :: pcu_hsr_v8(:,:), pcu_hsr_v10(:,:), pcu_hsr_v12(:,:)

    ! Data to allow more efficient memory usage (convert jprm to jprb only at point of use)
    REAL(KIND=jprm) :: cov_sfac           ! Stdev scale factor

    ! Arrays to hold hsr data interpolated onto channel wavenumbers
    INTEGER(KIND=jpim) :: platform_id
    INTEGER(KIND=jpim) :: sat_id
    INTEGER(KIND=jpim) :: inst_id
    INTEGER(KIND=jpim) :: ncoefchans                    ! Number of channels in coef file
    REAL(KIND=jprb), POINTER :: cov_emis_int(:,:)       ! dims are (cv_pack,ncoefchans)
    REAL(KIND=jprb), POINTER :: pcu_int_v12(:,:,:)      ! dims are (numpcs,ncoefchans,1/2)
    REAL(KIND=jprb), POINTER :: pcm_int_v12(:,:)
    REAL(KIND=jprb), POINTER :: pcu_int_v10(:,:,:)
    REAL(KIND=jprb), POINTER :: pcm_int_v10(:,:)
    REAL(KIND=jprb), POINTER :: pcu_int_v8(:,:,:)
    REAL(KIND=jprb), POINTER :: pcm_int_v8(:,:)
    REAL(KIND=jprb), POINTER :: sice_em_int(:), snow_em_int(:)
    REAL(KIND=jprm), POINTER :: p1d_int(:,:,:)
    REAL(KIND=jprm), POINTER :: p2d_int(:,:,:)
    REAL(KIND=jprm), POINTER :: p3d_int(:,:,:)
    REAL(KIND=jprm), POINTER :: p1n_int(:,:,:)
    REAL(KIND=jprm), POINTER :: p2n_int(:,:,:)
    REAL(KIND=jprm), POINTER :: p3n_int(:,:,:)
  END TYPE camel_atlas_data


  ! Atlas data contained in this file - shared, does not change

  REAL(KIND=jprb) :: sice_em(numwave), snow_em(numwave)
  REAL(KIND=jprb) :: sice_stdv, snow_stdv

  REAL(KIND=jprb) :: hsr_wavenum(numwave)


  DATA sice_stdv / 0.015_jprb /
  DATA sice_em / & !0.9370_jprb, &
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
      0.9682_jprb, 0.9681_jprb, 0.9679_jprb, 0.9677_jprb, 0.9676_jprb, 0.9674_jprb, 0.9671_jprb, 0.9669_jprb, 0.9669_jprb/

  DATA snow_stdv / 0.015_jprb /
  DATA snow_em / & !0.9716_jprb, &
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
      0.9815_jprb, 0.9814_jprb, 0.9813_jprb, 0.9813_jprb, 0.9812_jprb, 0.9811_jprb, 0.9811_jprb, 0.9810_jprb, 0.9810_jprb/

CONTAINS

!------------------------------------------
! Routines for initialising database
!------------------------------------------

  !> Initialise a CAMEL atlas data structure. By default the atlas data are
  !! general and can be used for any sensor with IR channels.
  !! If the last four optional arguments are specified then the atlas data are
  !! initialised for use with the specified instrument only: this makes calls
  !! to obtain emissivities from the atlas much faster.
  !! @param[in]       path            path to atlas data files
  !! @param[in]       imonth          month of data to read (1-12)
  !! @param[in]       verbose         flag to turn verbose output on/off
  !! @param[in]       std_init        flag to load covariance data
  !! @param[in]       do_ang_corr     flag to include zenith angle correction
  !! @param[in,out]   atlas           CAMEL atlas data structure to initialise
  !! @param[out]      err             status on exit
  !! @param[in]       instr_wavenum   channel wavenumber list, optional
  !! @param[in]       platform_id     platform ID from associated rtcoef structure, optional
  !! @param[in]       sat_id          satellite ID from associated rtcoef structure, optional
  !! @param[in]       inst_id         instrument ID from associated rtcoef structure, optional
  SUBROUTINE rttov_camel_init(       &
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
    ! initialize the rttov_uwiremis algorithm by (1) reading in the MEASURES CAMEL Global
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

    CHARACTER(LEN=*),       INTENT(IN)           :: path
    INTEGER(KIND=jpim),     INTENT(IN)           :: imonth
    LOGICAL(KIND=jplm),     INTENT(IN)           :: verbose
    LOGICAL(KIND=jplm),     INTENT(IN)           :: std_init
    LOGICAL(KIND=jplm),     INTENT(IN)           :: do_ang_corr
    TYPE(camel_atlas_data), INTENT(INOUT)        :: atlas
    INTEGER(KIND=jpim),     INTENT(OUT)          :: err
    REAL(KIND=jprb),        INTENT(IN), OPTIONAL :: instr_wavenum(:)
    INTEGER(KIND=jpim),     INTENT(IN), OPTIONAL :: platform_id
    INTEGER(KIND=jpim),     INTENT(IN), OPTIONAL :: sat_id
    INTEGER(KIND=jpim),     INTENT(IN), OPTIONAL :: inst_id

    INTEGER(KIND=jpim) :: i, npc_int
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

    CALL rttov_camel_nullify_pointers(atlas)

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
    IF (std_init) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'CAMEL climatology standard deviation data only available in HDF5 format')
    ENDIF
    WRITE(errmsg,'(a)') ', RTTOV was compiled with netCDF so netCDF atlas files are required'
#else
    err = errorstatus_fatal
    THROWM(err.NE.0,'RTTOV must be compiled with HDF5 to use the CAMEL atlas')
#endif
#endif

    !----------------------------------------------------------------------------
    ! reading the PCA Coefficients of CAMEL IR Land Surface Global Emissivity
    !----------------------------------------------------------------------------
    fn = TRIM(path)//'CAMEL_coef_'//cyear//cmonth//'_'//TRIM(prd_version)
#ifdef _RTTOV_HDF
    fn = TRIM(fn)//'.H5'
#else
    fn = TRIM(fn)//'.nc'
#endif
    INQUIRE(FILE=fn, EXIST=file_exists)
    IF (file_exists) THEN
      CALL rttov_camel_read_coefs(err, TRIM(fn), atlas)
      THROW(err.NE.0)
    ELSE
      err = errorstatus_fatal
      THROWM(err .NE. errorstatus_success, 'MEASURES CAMEL PCA coefs file not found: '//TRIM(fn)//TRIM(errmsg))
    ENDIF
    IF (verbose) INFO('Using CAMEL coefs: '//TRIM(fn))

    !----------------------------------------------------------------------------
    ! reading the 0.25 degree CAMEL climatology IR Land Surface Global Emissivty STD DEV
    !----------------------------------------------------------------------------
    IF (atlas%std_init) THEN
      fn = TRIM(path)//'CAMEL_hsremis_covmat_deg0.25_2000-2016_month'//cmonth//'_'//TRIM(prd_version_cov)//'_mask.H5'
      INQUIRE(FILE=fn, EXIST=file_exists)
      IF (file_exists) THEN
        CALL rttov_camel_read_cov(err, TRIM(fn), atlas)
        THROW(err.NE.0)
      ELSE
        err = errorstatus_fatal
        THROWM(err .NE. errorstatus_success, 'CAMEL covariances file not found: '//TRIM(fn))
      ENDIF
      IF (verbose) INFO('Using CAMEL covariances: '//TRIM(fn))
    ENDIF

    !-------------------------------------------------------------------
    ! reading the eigenvectors
    !-------------------------------------------------------------------
    ALLOCATE(atlas%pcm_hsr_v8(numwave),  atlas%pcu_hsr_v8(numwave,numwave),  &
             atlas%pcm_hsr_v10(numwave), atlas%pcu_hsr_v10(numwave,numwave), &
             atlas%pcm_hsr_v12(numwave), atlas%pcu_hsr_v12(numwave,numwave))
    fn = TRIM(path)//'pchsr_v8.2'
#ifdef _RTTOV_HDF
    fn = TRIM(fn)//'.H5'
#else
    fn = TRIM(fn)//'.nc'
#endif
    INQUIRE(FILE=fn, EXIST=file_exists)
    IF (file_exists) THEN
      CALL rttov_camel_read_labeigvects(err, TRIM(fn), '8  ', atlas)
      THROW(err.NE.0)
    ELSE
      err = errorstatus_fatal
      THROWM(err .NE. errorstatus_success, 'Laboratory data file not found: '//TRIM(fn)//TRIM(errmsg))
    ENDIF
    IF (verbose) INFO('Using CAMEL laboratory data file: '//TRIM(fn))

    !----------------------------------------------------------------
    fn = TRIM(path)//'pchsr_v10.2'
#ifdef _RTTOV_HDF
    fn = TRIM(fn)//'.H5'
#else
    fn = TRIM(fn)//'.nc'
#endif
    INQUIRE(FILE=fn, EXIST=file_exists)
    IF (file_exists) THEN
      CALL rttov_camel_read_labeigvects(err, TRIM(fn), '10 ', atlas)
      THROW(err.NE.0)
    ELSE
      err = errorstatus_fatal
      THROWM(err .NE. errorstatus_success, 'Laboratory data file not found: '//TRIM(fn)//TRIM(errmsg))
    ENDIF
    IF (verbose) INFO('Using CAMEL laboratory data file: '//TRIM(fn))

    !----------------------------------------------------------------
    fn = TRIM(path)//'pchsr_v12.2'
#ifdef _RTTOV_HDF
    fn = TRIM(fn)//'.H5'
#else
    fn = TRIM(fn)//'.nc'
#endif
    INQUIRE(FILE=fn, EXIST=file_exists)
    IF (file_exists) THEN
      CALL rttov_camel_read_labeigvects(err, TRIM(fn), '12 ', atlas)
      THROW(err.NE.0)
    ELSE
      err = errorstatus_fatal
      THROWM(err .NE. errorstatus_success, 'Laboratory data file not found: '//TRIM(fn)//TRIM(errmsg))
    ENDIF
    IF (verbose) INFO('Using CAMEL laboratory data file: '//TRIM(fn))


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
        CALL rttov_camel_read_angfunc(err, TRIM(fn), atlas)
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

    !--------------------------------------------------
    ! create hsr_wavenum array
    !--------------------------------------------------
    !  numwave = 2778:-5:694 (698)

    DO i = 1, numwave
      hsr_wavenum(i) = 698 + (i - 1) * 5
    ENDDO

    IF (atlas%single_inst) THEN
      ! If a channel wavenumber list is supplied for a particular instrument
      ! we can retrieve emissivities much faster.

      ! Interpolate hsr data onto the channel wavenumbers
      atlas%ncoefchans = SIZE(instr_wavenum)
      npc_int = 1
      IF (atlas%do_ang_corr) THEN
        npc_int = 2
        ALLOCATE(atlas%p1d_int(atlas%ncoefchans,atlas%nb_igbp,2), &
                 atlas%p2d_int(atlas%ncoefchans,atlas%nb_igbp,2), &
                 atlas%p3d_int(atlas%ncoefchans,atlas%nb_igbp,2), &
                 atlas%p1n_int(atlas%ncoefchans,atlas%nb_igbp,2), &
                 atlas%p2n_int(atlas%ncoefchans,atlas%nb_igbp,2), &
                 atlas%p3n_int(atlas%ncoefchans,atlas%nb_igbp,2))
      ENDIF
      ALLOCATE(atlas%pcu_int_v12(numpcs,atlas%ncoefchans,npc_int), &
               atlas%pcm_int_v12(atlas%ncoefchans,npc_int),        &
               atlas%pcu_int_v10(numpcs,atlas%ncoefchans,npc_int), &
               atlas%pcm_int_v10(atlas%ncoefchans,npc_int),        &
               atlas%pcu_int_v8(numpcs,atlas%ncoefchans,npc_int),  &
               atlas%pcm_int_v8(atlas%ncoefchans,npc_int))
      ALLOCATE(atlas%sice_em_int(atlas%ncoefchans), atlas%snow_em_int(atlas%ncoefchans))
      IF (atlas%std_init) ALLOCATE(atlas%cov_emis_int(atlas%cv_pack,atlas%ncoefchans))
      CALL rttov_camel_hsr_interp('12 ', instr_wavenum(:), atlas)
      CALL rttov_camel_hsr_interp('10 ', instr_wavenum(:), atlas)
      CALL rttov_camel_hsr_interp('8  ', instr_wavenum(:), atlas)

      ! HSR data is no longer required
      DEALLOCATE(atlas%pcm_hsr_v8, atlas%pcm_hsr_v10, atlas%pcm_hsr_v12, &
                 atlas%pcu_hsr_v8, atlas%pcu_hsr_v10, atlas%pcu_hsr_v12)
      IF (atlas%std_init)    DEALLOCATE(atlas%cov_emis)
      IF (atlas%do_ang_corr) DEALLOCATE(atlas%p1d, atlas%p2d, atlas%p3d, atlas%p1n, atlas%p2n, atlas%p3n)
    ENDIF

    CATCH
  END SUBROUTINE rttov_camel_init


#ifdef _RTTOV_HDF
  !> Read CAMEL atlas covariance data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   CAMEL atlas data structure
  SUBROUTINE rttov_camel_read_cov(err, fn, atlas)

    ! Description:
    ! read the 0.25 degree resolution CAMEL climatology covariance data
    ! from the HDF5 file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.1      03/06/2014 original code J. Vidot
    !  2.0      11/16/2018 Modified for CAMEL Climatology

    INTEGER(KIND=jpim),     INTENT(OUT)   :: err
    CHARACTER(LEN=*),       INTENT(IN)    :: fn
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpim)          :: indexlut
    INTEGER(KIND=jpim)          :: i, j
    INTEGER(KIND=jpis), POINTER :: pack_cov(:,:)

#include "rttov_hdf_load.interface"

    TRY

    CALL RTTOV_HDF_LOAD(err, fn, "/EMIS", SNAME="NUMOBS", PIS2=pack_cov)
    THROWM(ERR.NE.0,"Cannot load Emissivity Covariance numobs from "//TRIM(FN))
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

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="DIAGCOV", PIS2=atlas%cov_emis)
    THROWM(ERR.NE.0,"Cannot get covariance values from "//TRIM(fn))

    IF (SIZE(atlas%cov_emis,2) .NE. numwave) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in covariance data: numwave")
    ENDIF
    atlas%cv_pack = SIZE(atlas%cov_emis,1)

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="DIAGCOV_SFAC", RM0=atlas%cov_sfac)
    THROWM(ERR.NE.0,"Cannot get covariance scale factor from "//TRIM(fn))

    ! CALL rttov_hdf_load(err, fn, "/EMIS", SNAME='COV_EMIS_YGRID1', I0=atlas%cov_emis_ygrid1)
    ! THROW(ERR.NE.0)
    ! CALL rttov_hdf_load(err, fn, "/EMIS", SNAME='COV_EMIS_XGRID1', I0=atlas%cov_emis_xgrid1)
    ! THROW(ERR.NE.0)
    ! CALL rttov_hdf_load(err, fn, "/EMIS", SNAME='COV_EMIS_GRIDRES', I0=atlas%cov_emis_gridres)
    ! THROW(ERR.NE.0)

    CATCH
  END SUBROUTINE rttov_camel_read_cov
#else
  !> Read CAMEL atlas covariance data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   CAMEL atlas data structure
  SUBROUTINE  rttov_camel_read_cov(err, fn, atlas)
    ! Stub
    INTEGER(KIND=jpim),     INTENT(OUT)   :: err
    CHARACTER(LEN=*),       INTENT(IN)    :: fn
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas
    err = 0
  END SUBROUTINE rttov_camel_read_cov
#endif

#ifdef _RTTOV_HDF
  !> Read CAMEL atlas data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   CAMEL atlas data structure
  SUBROUTINE rttov_camel_read_coefs(err, fn, atlas)

    ! Description:
    ! read the 0.05 degree resolution MEASURES CAMEL Global Emissivty data
    ! from the HDF5 file into memory.

    INTEGER(KIND=jpim),     INTENT(OUT)   :: err
    CHARACTER(LEN=*),       INTENT(IN)    :: fn
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpim)          :: indexlut
    INTEGER(KIND=jpim)          :: i, j

#include "rttov_hdf_load.interface"
    TRY

    CALL RTTOV_HDF_LOAD(ERR, fn, "/EMIS", SNAME="CAMEL_QFLAG", PIS2=atlas%camel_qflag)
    THROWM(ERR.NE.0,"Cannot load CAMEL flag from "//TRIM(fn))

    atlas%nb_lons = SIZE(atlas%camel_qflag,1)
    atlas%nb_lats = SIZE(atlas%camel_qflag,2)
    ALLOCATE(atlas%coef_lut(atlas%nb_lons,atlas%nb_lats))

    ! Generate the look-up table into the emissivity data
    atlas%coef_lut(:,:) = -1_jpim
    indexlut = 1_jpim
    DO i = 1, atlas%nb_lats
      DO j = 1, atlas%nb_lons
        IF (atlas%camel_qflag(j,i) > 0) THEN
          atlas%coef_lut(j,i) = indexlut
          indexlut = indexlut + 1_jpim
        ENDIF
      ENDDO
    ENDDO

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="COEF", PRM2=atlas%pca_coef)
    THROWM(ERR.NE.0,"Cannot get CAMEL PCA coeffcients from "//TRIM(fn))

    atlas%nmask = SIZE(atlas%pca_coef,2)
    IF (SIZE(atlas%pca_coef,1) .NE. numpcs) THEN
      DEALLOCATE(atlas%pca_coef)
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in emissivity data: numpcs")
    ENDIF

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="PC_LABVS", PIS1=atlas%pc_labvs)
    THROWM(ERR.NE.0,"Cannot get PC_LABVS from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="PC_NPCS", PIS1=atlas%pc_npcs)
    THROWM(ERR.NE.0,"Cannot get PC_NPCS from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="CAMEL_YGRID1", I0=atlas%camel_ygrid1)
    THROWM(ERR.NE.0,"Cannot get CAMEL_YGRID1 from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="CAMEL_XGRID1", I0=atlas%camel_xgrid1)
    THROWM(ERR.NE.0,"Cannot get CAMEL_XGRID1 from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="CAMEL_GRIDRES", I0=atlas%camel_gridres)
    THROWM(ERR.NE.0,"Cannot get CAMEL_GRIDRES from "//TRIM(fn))

    CATCH
  END SUBROUTINE rttov_camel_read_coefs
#else
#ifdef _RTTOV_NETCDF
  !> Read CAMEL atlas data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   CAMEL atlas data structure
  SUBROUTINE  rttov_camel_read_coefs(err, fn, atlas)

    ! Description:
    ! read the 0.05 degree resolution MEASURES CAMEL Global Emissivty data
    ! from the netCDF file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.6    06/07/2016   by E. Borbas

    INTEGER(KIND=jpim),     INTENT(OUT)   :: err
    CHARACTER(LEN=*),       INTENT(IN)    :: fn
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpim) :: nvars     ! number of variables
    INTEGER(KIND=jpim) :: ndims     ! number of dimensions
    INTEGER(KIND=jpim) :: errstat   ! error code
    INTEGER(KIND=jpim) :: recdim    ! record dimension
    INTEGER(KIND=jpim) :: nc_dim(4) ! lats, lons, hng_pnt, mask
    INTEGER(KIND=jpim) :: i, j
    INTEGER(KIND=jpim) :: ncid, ngatts, nrecs, varid

    CHARACTER(LEN=1024) :: strbuf ! string buffer for var
    INTEGER(KIND=jpim)  :: indexlut

    CHARACTER (LEN=256) :: atlas_res_str
    REAL(KIND=jprm)     :: lat_max, lon_min, atlas_res

    TRY

    ! Open netCDF file.
    errstat = nf_open(TRIM(fn), nf_nowrite, ncid)

    ! Get info on the record dimension for this file.
    errstat = nf_inq(ncid, ndims, nvars, ngatts, recdim)

    DO recdim = 1, ndims
      errstat = nf_inq_dim(ncid, recdim, strbuf, nrecs)
      nc_dim(recdim) = nrecs
    ENDDO

    atlas%nb_lats = nc_dim(1)
    atlas%nb_lons = nc_dim(2)
    atlas%nmask   = nc_dim(4)

    ALLOCATE(atlas%camel_qflag(atlas%nb_lons,atlas%nb_lats), atlas%coef_lut(atlas%nb_lons,atlas%nb_lats), &
             atlas%pc_npcs(atlas%nmask), atlas%pc_labvs(atlas%nmask), atlas%pca_coef(numpcs,atlas%nmask))

    ! Retrieve emissivity database flag value
    errstat = nf_inq_varid(ncid, 'camel_qflag', varid)
    errstat = nf_get_var_int2(ncid, varid, atlas%camel_qflag)

    ! Generate the look-up table into the emissivity data
    atlas%coef_lut(:,:) = -1_jpim
    indexlut = 1_jpim

    DO i = 1, atlas%nb_lats
      DO j = 1, atlas%nb_lons
        IF (atlas%camel_qflag(j,i) > 0) THEN
          atlas%coef_lut(j,i) = indexlut
          indexlut = indexlut + 1_jpim
        ENDIF
      ENDDO
    ENDDO

    ! Retrieve verson number of PCs of laboratory data
    errstat = nf_inq_varid(ncid, 'pc_labvs', varid)
    errstat = nf_get_var_int2(ncid, varid, atlas%pc_labvs)

    ! Retrieve number of PCs for each pixel
    errstat = nf_inq_varid(ncid, 'pc_npcs', varid)
    errstat = nf_get_var_int2(ncid, varid, atlas%pc_npcs)

    ! Retrieve ASTER Snow Fraction
!     errstat = nf_inq_varid (ncid, 'snow_fraction', varid)
!     errstat = nf_get_var_int2(ncid, varid, sf)

    ! Retrieve coef database
    errstat = nf_inq_varid(ncid, 'pc_coefs', varid)
    errstat = nf_get_var_real(ncid, varid, atlas%pca_coef)

    errstat = nf_get_att_real(ncid, NF_GLOBAL, 'geospatial_lat_max', lat_max)
    errstat = nf_get_att_real(ncid, NF_GLOBAL, 'geospatial_lon_min', lon_min)
    errstat = nf_get_att_text(ncid, NF_GLOBAL, 'spatial_resolution', atlas_res_str)
    READ(atlas_res_str,*) atlas_res

    atlas%camel_ygrid1 = NINT(lat_max * 1000._jprb)
    atlas%camel_xgrid1 = NINT(lon_min * 1000._jprb)
    atlas%camel_gridres = NINT(atlas_res * 1000._jprb)

    errstat = nf_close(ncid)

    CATCH
  END SUBROUTINE rttov_camel_read_coefs
#else
  !> Read CAMEL atlas data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   CAMEL atlas data structure
  SUBROUTINE  rttov_camel_read_coefs(err, fn, atlas)
    ! Stub
    INTEGER(KIND=jpim),     INTENT(OUT)   :: err
    CHARACTER(LEN=*),       INTENT(IN)    :: fn
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas
    err = 0
  END SUBROUTINE rttov_camel_read_coefs
#endif
#endif

#ifdef _RTTOV_HDF
  !> Read CAMEL atlas eigenvector data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in]      clabvs  eigenvector file identifier
  !! @param[in,out]  atlas   CAMEL atlas data structure
  SUBROUTINE rttov_camel_read_labeigvects(err, fn, clabvs, atlas)

    ! Description:
    ! read the eigenvectors of selected laboratory measurements
    ! (created by E borbas) from the HDF5 file into memory.

    INTEGER(KIND=jpim),     INTENT(OUT)   :: err
    CHARACTER(LEN=*),       INTENT(IN)    :: fn
    CHARACTER(LEN=*),       INTENT(IN)    :: clabvs
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas

    REAL(KIND=jprb), POINTER :: pcm_x(:)
    REAL(KIND=jprb), POINTER :: pcu_x(:,:)

#include "rttov_hdf_load.interface"

    TRY

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="EIGENVALUES", PR1=pcm_x)
    THROWM(ERR.NE.0,"Cannot load eigenvalues from "//TRIM(fn))

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="EIGENVECTORS", PR2=pcu_x)
    THROWM(ERR.NE.0,"Cannot load eigenvectors from "//TRIM(fn))

    IF (SIZE(pcu_x,1) .NE. numwave .OR. SIZE(pcu_x,2) .NE. numwave .OR. SIZE(pcm_x) .NE. numwave) THEN
      DEALLOCATE(pcm_x, pcu_x)
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in eigenvalue/eigenvector data: numwave")
    ENDIF
    IF (clabvs == '12 ') THEN
      atlas%pcu_hsr_v12 = pcu_x
      atlas%pcm_hsr_v12 = pcm_x
    ELSEIF (clabvs == '10 ') THEN
      atlas%pcu_hsr_v10 = pcu_x
      atlas%pcm_hsr_v10 = pcm_x
    ELSEIF (clabvs == '8  ') THEN
      atlas%pcu_hsr_v8 = pcu_x
      atlas%pcm_hsr_v8 = pcm_x
    ENDIF
    DEALLOCATE(pcm_x, pcu_x)

    CATCH
  END SUBROUTINE rttov_camel_read_labeigvects
#else
#ifdef _RTTOV_NETCDF
  !> Read CAMEL atlas eigenvector data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in]      clabvs  eigenvector file identifier
  !! @param[in,out]  atlas   CAMEL atlas data structure
  SUBROUTINE  rttov_camel_read_labeigvects(err, fn, clabvs, atlas)

    ! Description:
    ! read the eigenvectors of selected laboratory measurements
    ! (created by E borbas) from the netCDF file into memory.
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  1.0    03/31/2009   origianl code B. Ruston

    INTEGER(KIND=jpim),     INTENT(OUT)   :: err
    CHARACTER(LEN=*),       INTENT(IN)    :: fn
    CHARACTER(LEN=*),       INTENT(IN)    :: clabvs
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas

    CHARACTER(LEN=1024) :: strbuf ! string buffer for var

    INTEGER(KIND=jpim) :: nvars     ! number of variables
    INTEGER(KIND=jpim) :: ndims     ! number of dimensions
    INTEGER(KIND=jpim) :: errstat   ! error code
    INTEGER(KIND=jpim) :: recdim    ! record dimension
    INTEGER(KIND=jpim) :: nc_dim(2) ! hng_pnt, lats, lons, pack_len
    INTEGER(KIND=jpim) :: ncid, ngatts, nrecs, varid

    REAL(KIND=jprb) :: pcm(numwave), pcu(numwave,numwave)

    TRY

    ! Open netCDF file.
    errstat = nf_open(TRIM(fn), nf_nowrite, ncid)

    ! Get info on the record dimension for this file.
    errstat = nf_inq(ncid, ndims, nvars, ngatts, recdim)

    DO recdim = 1, ndims
      errstat = nf_inq_dim(ncid, recdim, strbuf, nrecs)
      nc_dim(recdim) = nrecs
    ENDDO

    IF (nc_dim(1) .NE. numwave) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in eigenvector data: numwave")
    ENDIF

    ! Retrieve eigenvalues
    errstat = nf_inq_varid(ncid, 'Eigenvalues', varid)
    errstat = nf_get_var_double(ncid, varid, pcm)

    ! Retrieve eigenvectors
    errstat = nf_inq_varid(ncid, 'Eigenvectors', varid)
    errstat = nf_get_var_double(ncid, varid, pcu)

    errstat = nf_close(ncid)

    IF (clabvs == '12 ') THEN
      atlas%pcu_hsr_v12 = pcu
      atlas%pcm_hsr_v12 = pcm
    ELSEIF (clabvs == '10 ') THEN
      atlas%pcu_hsr_v10 = pcu
      atlas%pcm_hsr_v10 = pcm
    ELSEIF (clabvs == '8  ') THEN
      atlas%pcu_hsr_v8 = pcu
      atlas%pcm_hsr_v8 = pcm
    ENDIF

    CATCH
  END SUBROUTINE rttov_camel_read_labeigvects
#else
  !> Read CAMEL atlas eigenvector data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in]      clabvs  eigenvector file identifier
  !! @param[in,out]  atlas   CAMEL atlas data structure
  SUBROUTINE  rttov_camel_read_labeigvects(err, fn, clabvs, atlas)
    ! Stub
    INTEGER(KIND=jpim),     INTENT(OUT)   :: err
    CHARACTER(LEN=*),       INTENT(IN)    :: fn
    CHARACTER(LEN=*),       INTENT(IN)    :: clabvs
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas
    err = 0
  END SUBROUTINE rttov_camel_read_labeigvects
#endif
#endif

#ifdef _RTTOV_HDF
  !> Read CAMEL atlas angular correction data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   CAMEL atlas data structure
  SUBROUTINE  rttov_camel_read_angfunc(err, fn, atlas)

    ! Description:
    ! read the three day-night coefficients for angular correcton
    ! (created by E borbas) from the HDF5 file into memory.

    INTEGER(KIND=jpim),     INTENT(OUT)   :: err
    CHARACTER(LEN=*),       INTENT(IN)    :: fn
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas

    REAL(KIND=jprm), POINTER :: pin(:,:)

#include "rttov_hdf_load.interface"

    TRY

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p1d", PRM2=pin)
    THROWM(ERR.NE.0,"Cannot load angular correction day 1 coefs data from "//TRIM(fn))

    IF (SIZE(pin,1) .NE. numwave-1) THEN
      DEALLOCATE(atlas%p1d)
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in angular correction data: numwave")
    ENDIF

    atlas%nb_igbp = SIZE(pin,2)

    ALLOCATE(atlas%p1d(numwave,atlas%nb_igbp))
    atlas%p1d(1:numwave-1,:) = pin
    atlas%p1d(numwave,:) = atlas%p1d(numwave-1,:)
    DEALLOCATE(pin)

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p2d", PRM2=pin)
    THROWM(ERR.NE.0,"Cannot load angular correction day 2 coefs data from "//TRIM(fn))
    ALLOCATE(atlas%p2d(numwave,atlas%nb_igbp))
    atlas%p2d(1:numwave-1,:) = pin
    atlas%p2d(numwave,:) = atlas%p2d(numwave-1,:)
    DEALLOCATE(pin)

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p3d", PRM2=pin)
    THROWM(ERR.NE.0,"Cannot load angular correction day 3 coefs data from "//TRIM(fn))
    ALLOCATE(atlas%p3d(numwave,atlas%nb_igbp))
    atlas%p3d(1:numwave-1,:) = pin
    atlas%p3d(numwave,:) = atlas%p3d(numwave-1,:)
    DEALLOCATE(pin)

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p1n", PRM2=pin)
    THROWM(ERR.NE.0,"Cannot load angular correction night 1 coefs data from "//TRIM(fn))
    ALLOCATE(atlas%p1n(numwave,atlas%nb_igbp))
    atlas%p1n(1:numwave-1,:) = pin
    atlas%p1n(numwave,:) = atlas%p1n(numwave-1,:)
    DEALLOCATE(pin)

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p2n", PRM2=pin)
    THROWM(ERR.NE.0,"Cannot load angular correction night 2 coefs data from "//TRIM(fn))
    ALLOCATE(atlas%p2n(numwave,atlas%nb_igbp))
    atlas%p2n(1:numwave-1,:) = pin
    atlas%p2n(numwave,:) = atlas%p2n(numwave-1,:)
    DEALLOCATE(pin)

    CALL rttov_hdf_load(err, fn, "/EMIS", SNAME="p3n", PRM2=pin)
    THROWM(ERR.NE.0,"Cannot load angular correction night 3 coefs data from "//TRIM(fn))
    ALLOCATE(atlas%p3n(numwave,atlas%nb_igbp))
    atlas%p3n(1:numwave-1,:) = pin
    atlas%p3n(numwave,:) = atlas%p3n(numwave-1,:)
    DEALLOCATE(pin)

    CATCH
  END SUBROUTINE rttov_camel_read_angfunc
#else
#ifdef _RTTOV_NETCDF
  !> Read CAMEL atlas angular correction data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   CAMEL atlas data structure
  SUBROUTINE  rttov_camel_read_angfunc(err, fn, atlas)

    ! Description:
    ! read the three day-night coefficients for angular correcton
    ! (created by E borbas) from the netCDF file into memory.

    INTEGER(KIND=jpim),     INTENT(OUT)   :: err
    CHARACTER(LEN=*),       INTENT(IN)    :: fn
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpim) :: nvars     ! number of variables
    INTEGER(KIND=jpim) :: ndims     ! number of dimensions
    INTEGER(KIND=jpim) :: errstat   ! error code
    INTEGER(KIND=jpim) :: recdim    ! record dimension
    INTEGER(KIND=jpim) :: nc_dim(2)

    INTEGER(KIND=jpim) :: ncid, ngatts, nrecs, varid
    CHARACTER (len=1024) :: strbuf ! string buffer for var
    REAL(KIND=jprm), POINTER :: pin(:,:)

    TRY

   ! Open netCDF file
    errstat = nf_open(TRIM(fn), nf_nowrite, ncid)

    ! Get info on the record dimension for this file.
    errstat = nf_inq(ncid, ndims, nvars, ngatts, recdim)

    DO recdim = 1, ndims
      errstat = nf_inq_dim(ncid, recdim, strbuf, nrecs)
      nc_dim(recdim) = nrecs
    ENDDO

    IF (nc_dim(1) .NE. numwave-1) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Inconsistent dimensions in angular correction data: numwave")
    ENDIF
    atlas%nb_igbp = nc_dim(2)
    ALLOCATE(pin(numwave-1,atlas%nb_igbp))
    ALLOCATE(atlas%p1d(numwave,atlas%nb_igbp))
    ALLOCATE(atlas%p2d(numwave,atlas%nb_igbp))
    ALLOCATE(atlas%p3d(numwave,atlas%nb_igbp))
    ALLOCATE(atlas%p1n(numwave,atlas%nb_igbp))
    ALLOCATE(atlas%p2n(numwave,atlas%nb_igbp))
    ALLOCATE(atlas%p3n(numwave,atlas%nb_igbp))

    ! Retrieve p(3) for day and night
    errstat = nf_inq_varid (ncid, 'p1d', varid)
    errstat = nf_get_var_real(ncid, varid, pin)
    atlas%p1d(1:numwave-1,:) = pin
    atlas%p1d(numwave,:) = atlas%p1d(numwave-1,:)
    errstat = nf_inq_varid (ncid, 'p2d', varid)
    errstat = nf_get_var_real(ncid, varid, pin)
    atlas%p2d(1:numwave-1,:) = pin
    atlas%p2d(numwave,:) = atlas%p2d(numwave-1,:)
    errstat = nf_inq_varid (ncid, 'p3d', varid)
    errstat = nf_get_var_real(ncid, varid, pin)
    atlas%p3d(1:numwave-1,:) = pin
    atlas%p3d(numwave,:) = atlas%p3d(numwave-1,:)
    errstat = nf_inq_varid (ncid, 'p1n', varid)
    errstat = nf_get_var_real(ncid, varid, pin)
    atlas%p1n(1:numwave-1,:) = pin
    atlas%p1n(numwave,:) = atlas%p1n(numwave-1,:)
    errstat = nf_inq_varid (ncid, 'p2n', varid)
    errstat = nf_get_var_real(ncid, varid, pin)
    atlas%p2n(1:numwave-1,:) = pin
    atlas%p2n(numwave,:) = atlas%p2n(numwave-1,:)
    errstat = nf_inq_varid (ncid, 'p3n', varid)
    errstat = nf_get_var_real(ncid, varid, pin)
    atlas%p3n(1:numwave-1,:) = pin
    atlas%p3n(numwave,:) = atlas%p3n(numwave-1,:)

    errstat = nf_close(ncid)

    CATCH
  END SUBROUTINE rttov_camel_read_angfunc
#else
  !> Read CAMEL atlas angular correction data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   CAMEL atlas data structure
  SUBROUTINE  rttov_camel_read_angfunc(err, fn, atlas)
    ! Stub
    INTEGER(KIND=jpim),     INTENT(OUT)   :: err
    CHARACTER(LEN=*),       INTENT(IN)    :: fn
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas
    err = 0
  END SUBROUTINE rttov_camel_read_angfunc
#endif
#endif

#ifdef _RTTOV_HDF
  !> Read CAMEL atlas IGBP data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   CAMEL atlas data structure
  SUBROUTINE  rttov_uwiremis_read_igbp(err, fn, atlas)

    ! Description:
    ! read in the IGBP ecosystem map interpolated for the 0.05x0.05 degree grid and
    ! modified for IR emissivity application
    ! (created by E borbas) from the HDF5 file into memory.

    INTEGER(KIND=jpim),     INTENT(OUT)   :: err
    CHARACTER(LEN=*),       INTENT(IN)    :: fn
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas

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
  !> Read CAMEL atlas IGBP data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   CAMEL atlas data structure
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

    INTEGER(KIND=jpim),     INTENT(OUT)   :: err
    CHARACTER(LEN=*),       INTENT(IN)    :: fn
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas

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
  !> Read CAMEL atlas IGBP data file
  !! @param[out]     err     status on exit
  !! @param[in]      fn      full path of file to read
  !! @param[in,out]  atlas   CAMEL atlas data structure
  SUBROUTINE  rttov_uwiremis_read_igbp(err, fn, atlas)
    ! Stub
    INTEGER(KIND=jpim),     INTENT(OUT)   :: err
    CHARACTER(LEN=*),       INTENT(IN)    :: fn
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas
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
  !! @param[in]       atlas               CAMEL atlas data structure
  !! @param[out]      instr_emis          calculated emissivities values
  !! @param[out]      instr_emis_cov      estimated errors in emissivities
  !! @param[out]      instr_emis_flag     quality flags for returned emissivities
  !! @param[out]      instr_coeff         PC eigenvalues for returned emissivities, optional
  !! @param[out]      instr_pcu           PC eigenvectors for returned emissivities, optional
  !! @param[out]      instr_pcm           PC constants for returned emissivities, optional
  SUBROUTINE rttov_camel( &
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
    ! from the 0.05 degree resolution MEASURES CAMEL Global Emissivity data
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
    TYPE(camel_atlas_data),    INTENT(IN)  :: atlas
    REAL(KIND=jprb),           INTENT(OUT) :: instr_emis(nchs)
    REAL(KIND=jprb),           INTENT(OUT) :: instr_emis_cov(nchs)
    INTEGER(KIND=jpim),        INTENT(OUT) :: instr_emis_flag
    REAL(KIND=jprb), OPTIONAL, INTENT(OUT) :: instr_coeff(numpcs)
    REAL(KIND=jprb), OPTIONAL, INTENT(OUT) :: instr_pcu(numpcs,nchs)
    REAL(KIND=jprb), OPTIONAL, INTENT(OUT) :: instr_pcm(nchs)

    REAL(KIND=jprb)    :: coeff(numpcs)
    INTEGER(KIND=jpim) :: npcs, labvs

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
      IF (atlas%do_ang_corr) THEN
        instr_pcu_angcorr = 0._jprb
        instr_pcm_angcorr = 0._jprb
      ENDIF
    ENDIF

    IF (atlas%single_inst) THEN
      IF (atlas%platform_id /= platform_id .OR. &
          atlas%inst_id /= inst_id .OR. &
          atlas%ncoefchans /= ncoefchans) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0,'CAMEL atlas called for different coefs to that with which it was initialised')
      ENDIF
      IF (atlas%sat_id /= sat_id) THEN
        IF (verbose) &
          WARN('WARNING: CAMEL atlas called for instrument with different sat_id to that with which it was initialised')
      ENDIF
    ENDIF

    IF (surfacetype == surftype_land) THEN

    !------------------------------------------------------------
    ! find the closest grid point from the camel database
    !------------------------------------------------------------

      long = MODULO(lon, 360.0_jprb)
      IF (long >= 180.0) THEN
        long = long - 360.0_jprb
      ENDIF

      ilat = NINT(lat * 1000._jprb, KIND=jpim)
      ilon = NINT(long * 1000._jprb, KIND=jpim)

      gridy = NINT(ABS(atlas%camel_ygrid1 - ilat) * 1._jprb / atlas%camel_gridres, KIND=jpim) + 1_jpim
      gridx = NINT(ABS(atlas%camel_xgrid1 - ilon) * 1._jprb / atlas%camel_gridres, KIND=jpim) + 1_jpim
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

      instr_emis_flag = atlas%camel_qflag(gridx,gridy)

    !------------------------------
    ! check if it is a land pixel
    !------------------------------

      IF (instr_emis_flag > 0) THEN

        ! Find the emissivity or coefs and covariances

        IF (atlas%coef_lut(gridx,gridy) > 0) THEN
          coeff(:) = REAL(atlas%pca_coef(:,atlas%coef_lut(gridx,gridy)), KIND=jprb)
          IF (get_pc) instr_coeff(:) = coeff(:)
          labvs = atlas%pc_labvs(atlas%coef_lut(gridx,gridy))
          npcs = atlas%pc_npcs(atlas%coef_lut(gridx,gridy))
        ELSE
          RETURN
        ENDIF

        IF (atlas%std_init) THEN
          IF (atlas%single_inst) THEN
            IF (atlas%cov_emis_lut(rnd_y,rnd_x) > 0) THEN
              ! Note cov_emis_int contains the standard deviations (i.e. sqrt is already taken)
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
              CALL rttov_camel_recon_emis( &
                  npcs,                       &
                  labvs,                      &
                  atlas,                      &
                  coeff,                      &
                  channels,                   &
                  instr_emis_angcorr(:,j),    &
                  get_pc,                     &
                  instr_pcu_angcorr(:,:,j),   &
                  instr_pcm_angcorr(:,j), j)
            ENDDO

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
                DO i = 1, npcs
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

            CALL rttov_camel_recon_emis( &
                npcs,                    &! in
                labvs,                   &! in
                atlas,                   &! in
                coeff,                   &! in
                channels,                &! in
                instr_emis,              &! out
                get_pc,                  &! in
                instr_pcu,               &! out
                instr_pcm)                ! out

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
          ! compute the hsr emissivity spectra at 417 wavenumber points from the 10 BF emissivity hinge points
          !-------------------------------------------------------------------------------------------------------

          CALL rttov_camel_recon_hsremis( &
              npcs,                          &! in
              labvs,                         &! in
              atlas,                         &! in
              coeff,                         &! in
              hsremis)                        ! out

          !--------------------------------------------------------------------------------------------------------
          ! apply angular correction to the hsr emissivity spectra at 417 wavenumber points
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
          ! create instrument specific emis/stdv by finding the closest wavenumber value
          !--------------------------------------------------------------------------------

          IF (get_pc) THEN
            CALL rttov_uwiremis_select_wavenum( &
                npcs,                           &! in
                labvs,                          &! in
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
                npcs,                           &! in
                labvs,                          &! in
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
                      npcs,                     &! in
                      labvs,                    &! in
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
                    npcs,                   &! in
                    labvs,                  &! in
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
  END SUBROUTINE rttov_camel


  !> Interpolate PCs onto instrument channel wavenumbers. This is called
  !! during initialisation when single-instrument init is selected.
  !! @param[in]      clabvs          eigenvector file identifier
  !! @param[in]      instr_wavenum   instrument channel wavenumbers
  !! @param[in,out]  atlas           CAMEL atlas data structure
  SUBROUTINE rttov_camel_hsr_interp(clabvs, instr_wavenum, atlas)

    ! Description:
    ! Initialisation for a single instrument.
    ! Interpolate PC data onto a specific set of wavenumbers:
    ! this can be precomputed for a given instrument during
    ! atlas initialisation to enable very rapid emissivity
    ! calculations.

    CHARACTER(LEN=3),       INTENT(IN)    :: clabvs
    REAL(KIND=jprb),        INTENT(IN)    :: instr_wavenum(:)
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas

    INTEGER(KIND=jpim) :: j, k, nchs

    REAL(KIND=jprb) :: dist(numwave)
    REAL(KIND=jprb) :: mindist
    INTEGER(KIND=jpim) :: ind_mindist

    REAL(KIND=jprb) :: dwvnum1, dwvnum2
    REAL(KIND=jprb) :: pcu1(numpcs), pcu2(numpcs), pcm1, pcm2
    REAL(KIND=jprb) :: sice_em1, sice_em2, snow_em1, snow_em2
    REAL(KIND=jprb) :: cov_emis1(atlas%cv_pack), cov_emis2(atlas%cv_pack)

    REAL(KIND=jprb) :: pcm(numwave), pcu(numpcs,numwave)
    REAL(KIND=jprb) :: pcm_int(atlas%ncoefchans,SIZE(atlas%pcu_int_v8,3))
    REAL(KIND=jprb) :: pcu_int(numpcs,atlas%ncoefchans,SIZE(atlas%pcu_int_v8,3))

    IF (clabvs == '12 ') THEN
      pcu = atlas%pcu_hsr_v12(1:numpcs,:)
      pcm = atlas%pcm_hsr_v12
    ELSEIF (clabvs == '10 ') THEN
      pcu = atlas%pcu_hsr_v10(1:numpcs,:)
      pcm = atlas%pcm_hsr_v10
    ELSEIF (clabvs == '8  ') THEN
      pcu = atlas%pcu_hsr_v8(1:numpcs,:)
      pcm = atlas%pcm_hsr_v8
    ENDIF

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

        pcu_int(:,j,1)   = REAL(pcu(:,1), KIND=jprb)
        pcm_int(j,1)     = pcm(1)
        IF (atlas%do_ang_corr) THEN
          pcu_int(:,j,2) = 0._jprb
          pcm_int(j,2)   = 0._jprb
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

        pcu_int(:,j,1)   = REAL(pcu(:,numwave), KIND=jprb)
        pcm_int(j,1)     = pcm(numwave)
        IF (atlas%do_ang_corr) THEN
          pcu_int(:,j,2) = 0._jprb
          pcm_int(j,2)   = 0._jprb
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

        pcu1(:)  = dwvnum1 * REAL(pcu(:,ind_mindist + k), KIND=jprb)
        pcu2(:)  = dwvnum2 * REAL(pcu(:,ind_mindist), KIND=jprb)
        pcm1     = dwvnum1 * pcm(ind_mindist + k)
        pcm2     = dwvnum2 * pcm(ind_mindist)
        sice_em1 = dwvnum1 * sice_em(ind_mindist + k)
        sice_em2 = dwvnum2 * sice_em(ind_mindist)
        snow_em1 = dwvnum1 * snow_em(ind_mindist + k)
        snow_em2 = dwvnum2 * snow_em(ind_mindist)

        IF (atlas%do_ang_corr) THEN
          pcu_int(:,j,1) = pcu1(:)
          pcm_int(j,1)   = pcm1
          pcu_int(:,j,2) = pcu2(:)
          pcm_int(j,2)   = pcm2
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
          pcu_int(:,j,1) = pcu1(:) + pcu2(:)
          pcm_int(j,1)   = pcm1 + pcm2
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

    IF (clabvs == '12 ') THEN
      atlas%pcu_int_v12 = pcu_int
      atlas%pcm_int_v12 = pcm_int
    ELSEIF (clabvs == '10 ') THEN
      atlas%pcu_int_v10 = pcu_int
      atlas%pcm_int_v10 = pcm_int
    ELSEIF (clabvs == '8  ') THEN
      atlas%pcu_int_v8 = pcu_int
      atlas%pcm_int_v8 = pcm_int
    ENDIF

  END SUBROUTINE rttov_camel_hsr_interp


  !> Reconstruct emissivities at instrument wavenumbers from interpolated PCs
  !! (used when atlas is initialised for a specific instrument)
  !! @param[in]      npcs            number of PC scores to use in reconstruction
  !! @param[in]      labvs           eigenvector set identifer
  !! @param[in]      atlas           CAMEL atlas data structure
  !! @param[in]      coef            PC coefficients
  !! @param[in]      channels        channels for which emissivities are required
  !! @param[out]     emis            reconstructed emissivities
  !! @param[in]      get_pc          flag to indicate if pcu/pcm should be output
  !! @param[in,out]  pcu             PC eigenvectors, optional
  !! @param[in,out]  pcm             PC constants, optional
  !! @param[in]      ind             selects which of the two interpolated sets of PCs to use
  !!                                   (used with angular correction), optional
  SUBROUTINE rttov_camel_recon_emis(npcs, labvs, atlas, coef, channels, emis, get_pc, pcu, pcm, ind)

    ! Description:
    ! Used with single-instrument initialisation.
    ! Reconstruct the emissivities at the instrument wavenumbers
    ! from interpolated Principal Components.

    INTEGER(KIND=jpim),     INTENT(IN)              :: npcs, labvs
    TYPE(camel_atlas_data), INTENT(IN)              :: atlas
    REAL(KIND=jprb),        INTENT(IN)              :: coef(numpcs)
    INTEGER(KIND=jpim),     INTENT(IN)              :: channels(:)
    REAL(KIND=jprb),        INTENT(OUT)             :: emis(SIZE(channels))
    LOGICAL(KIND=jplm),     INTENT(IN)              :: get_pc
    REAL(KIND=jprb),        INTENT(INOUT), OPTIONAL :: pcu(numpcs,SIZE(channels))
    REAL(KIND=jprb),        INTENT(INOUT), OPTIONAL :: pcm(SIZE(channels))
    INTEGER(KIND=jpim),     INTENT(IN),    OPTIONAL :: ind

    INTEGER(KIND=jpim) :: j, k, nchn
    REAL(KIND=jprb), POINTER :: pcu_int(:,:,:), pcm_int(:,:)

    IF (labvs == 12) THEN
      pcu_int => atlas%pcu_int_v12
      pcm_int => atlas%pcm_int_v12
    ELSEIF (labvs == 10) THEN
      pcu_int => atlas%pcu_int_v10
      pcm_int => atlas%pcm_int_v10
    ELSEIF (labvs == 8) THEN
      pcu_int => atlas%pcu_int_v8
      pcm_int => atlas%pcm_int_v8
    ENDIF

    !-----------------------------------
    ! apply regcoef to get the emissivities
    !-----------------------------------

    IF (coef(1) /= -999._jprb) THEN
      j = 1
      IF (PRESENT(ind)) j = ind
      nchn = SIZE(emis)
      DO k = 1, nchn
        emis(k) = SUM(coef(1:npcs) * pcu_int(1:npcs,channels(k),j)) + pcm_int(channels(k),j)
      ENDDO
      IF (get_pc) THEN
        pcu(1:npcs,:) = pcu_int(1:npcs,channels(:),j)
        pcm = pcm_int(channels(:),j)
      ENDIF
    ELSE
      emis = hkod
    ENDIF

  END SUBROUTINE rttov_camel_recon_emis

  !> Reconstruct high-spectral-resolution emissivities from PCs
  !! @param[in]   npcs            number of PC scores to use in reconstruction
  !! @param[in]   labvs           eigenvector set identifer
  !! @param[in]   atlas           CAMEL atlas data structure
  !! @param[in]   coef            PC coefficients
  !! @param[out]  hsremis         reconstructed emissivities
  SUBROUTINE rttov_camel_recon_hsremis(npcs, labvs, atlas, coef, hsremis)

    ! Description:
    ! Used with multiple-instrument initialisation.
    ! To creates high spectra resolution emissivties at 417 wavenumbers
    ! from the PCA Coefficitents of the MEASURES CAMEL Global Emissivty data
    ! and labratory measurements using principal component analyses
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.9       03/31/2009  Original code E Borbas UW-Madison/CIMSS
    !  1.0       03/31/2009  New F90 code with structures (E Borbas B Ruston)
    !  1.1       11/30/2012  Removed the coef calcualtion part (E Borbas )
    !  0.6       06/06/2016  Modified for CAMEL DB  (E Borbas)

    INTEGER(KIND=jpim),     INTENT(IN)  :: npcs, labvs
    TYPE(camel_atlas_data), INTENT(IN)  :: atlas
    REAL(KIND=jprb),        INTENT(IN)  :: coef(numpcs)
    REAL(KIND=jprb),        INTENT(OUT) :: hsremis(numwave)
    REAL(KIND=jprb), POINTER :: pcu(:,:), pcm(:)

    INTEGER(KIND=jpim) :: k

    IF (labvs == 12) THEN
      pcu => atlas%pcu_hsr_v12
      pcm => atlas%pcm_hsr_v12
    ELSEIF (labvs == 10) THEN
      pcu => atlas%pcu_hsr_v10
      pcm => atlas%pcm_hsr_v10
    ELSEIF (labvs == 8) THEN
      pcu => atlas%pcu_hsr_v8
      pcm => atlas%pcm_hsr_v8
    ENDIF

    !-----------------------------------
    ! apply regcoef to get the hsr dataset
    !-----------------------------------

    IF (coef(1) /= -999._jprb) THEN
      DO k = 1, numwave
        hsremis(k) = SUM(coef(1:npcs) * REAL(pcu(1:npcs,k), KIND=jprb)) + pcm(k)
      ENDDO
    ELSE
      hsremis = hkod
    ENDIF

  END SUBROUTINE rttov_camel_recon_hsremis


  !> Linear interpolation of high-spectral-resolution emissivity data onto
  !! instrument channel wavenumbers
  !! @param[in]      npcs            number of PC scores to use in reconstruction
  !! @param[in]      labvs           eigenvector set identifer
  !! @param[in]      atlas           CAMEL atlas data structure
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
          npcs,                              &
          labvs,                             &
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
    ! Subroutine to find the closest wavenumber from the MEASURES CAMEL HSR emissivity spectra
    ! for the instrument frequency and assign the instrument emissivity by choosing the
    ! closest spectral point value or bilinear interpolating  between the two
    ! closest spectral point values
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    !  0.9       03/31/2009  Original code E Borbas UW-Madison/CIMSS
    !  1.0       03/31/2009  New F90 code with structures (E Borbas B Ruston)

    INTEGER(KIND=jpim),        INTENT(IN)    :: npcs, labvs
    TYPE(camel_atlas_data),    INTENT(IN)    :: atlas
    REAL(KIND=jprb),           INTENT(IN)    :: hsremis(numwave)
    REAL(KIND=jprb),           INTENT(IN)    :: emis_cov(numwave)
    INTEGER(KIND=jpim),        INTENT(IN)    :: nchs
    REAL(KIND=jprb),           INTENT(IN)    :: instr_wavenum(nchs)
    REAL(KIND=jprb),           INTENT(OUT)   :: instr_emis(nchs)
    REAL(KIND=jprb),           INTENT(OUT)   :: instr_emis_cov(nchs)
    REAL(KIND=jprb), OPTIONAL, INTENT(INOUT) :: instr_pcu(numpcs,nchs)
    REAL(KIND=jprb), OPTIONAL, INTENT(INOUT) :: instr_pcm(nchs)
    REAL(KIND=jprb), OPTIONAL, INTENT(IN)    :: angcorr(numwave)

    REAL(KIND=jprb), POINTER :: pcu(:,:), pcm(:)
    INTEGER(KIND=jpim) :: j, k

    REAL(KIND=jprb) :: dist(numwave)
    REAL(KIND=jprb) :: mindist
    INTEGER(KIND=jpim) :: ind_mindist

    REAL(KIND=jprb) :: dwvnum1, dwvnum2, dwvsum
    REAL(KIND=jprb) :: hsremis1, hsremis2, emis_cov1, emis_cov2
    REAL(KIND=jprb) :: pcu1(npcs), pcu2(npcs), pcm1, pcm2
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

    IF (lcpu_pc) THEN
      IF (labvs == 12) THEN
        pcu => atlas%pcu_hsr_v12
        pcm => atlas%pcm_hsr_v12
      ELSEIF (labvs == 10) THEN
        pcu => atlas%pcu_hsr_v10
        pcm => atlas%pcm_hsr_v10
      ELSEIF (labvs == 8) THEN
        pcu => atlas%pcu_hsr_v8
        pcm => atlas%pcm_hsr_v8
      ENDIF
    ENDIF

    IF (lcpu_emis .OR. lcpu_cov) THEN
      instr_emis(:)       = hkod
      instr_emis_cov(:)   = hkod
      DO j = 1, nchs

        IF (instr_wavenum(j) <= hsr_wavenum(1)) THEN

          instr_emis(j) = hsremis(1)
          IF (atlas%std_init) instr_emis_cov(j) = emis_cov(1)
          IF (lcpu_pc) THEN
            instr_pcu(1:npcs,j) = pcu(1:npcs,1) * angcorr(1)
            instr_pcm(j) = pcm(1) * angcorr(1)
          ENDIF

        ELSEIF (instr_wavenum(j) >= hsr_wavenum(numwave)) THEN

          instr_emis(j) = hsremis(numwave)
          IF (atlas%std_init) instr_emis_cov(j) = emis_cov(numwave)
          IF (lcpu_pc) THEN
            instr_pcu(1:npcs,j) = pcu(1:npcs,numwave) * angcorr(numwave)
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
            pcu1 = dwvnum1 * pcu(1:npcs,ind_mindist + k) * angcorr(ind_mindist + k)
            pcu2 = dwvnum2 * pcu(1:npcs,ind_mindist) * angcorr(ind_mindist)
            instr_pcu(1:npcs,j) = (pcu1 + pcu2) / dwvsum
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
  !> Deallocate data in CAMEL atlas data structure
  !! @param[in,out]   atlas   CAMEL atlas data structure to deallocate
  SUBROUTINE rttov_camel_close_atlas(atlas)
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas

    IF ( ASSOCIATED(atlas%camel_qflag) ) DEALLOCATE(atlas%camel_qflag)
    IF ( ASSOCIATED(atlas%pc_labvs)    ) DEALLOCATE(atlas%pc_labvs)
    IF ( ASSOCIATED(atlas%pc_npcs)     ) DEALLOCATE(atlas%pc_npcs)
    IF ( ASSOCIATED(atlas%coef_lut)    ) DEALLOCATE(atlas%coef_lut)
    IF ( ASSOCIATED(atlas%pca_coef)    ) DEALLOCATE(atlas%pca_coef)
    IF ( ASSOCIATED(atlas%pcm_hsr_v8)  ) DEALLOCATE(atlas%pcm_hsr_v8)
    IF ( ASSOCIATED(atlas%pcm_hsr_v10) ) DEALLOCATE(atlas%pcm_hsr_v10)
    IF ( ASSOCIATED(atlas%pcm_hsr_v12) ) DEALLOCATE(atlas%pcm_hsr_v12)
    IF ( ASSOCIATED(atlas%pcu_hsr_v8)  ) DEALLOCATE(atlas%pcu_hsr_v8)
    IF ( ASSOCIATED(atlas%pcu_hsr_v10) ) DEALLOCATE(atlas%pcu_hsr_v10)
    IF ( ASSOCIATED(atlas%pcu_hsr_v12) ) DEALLOCATE(atlas%pcu_hsr_v12)
    IF ( ASSOCIATED(atlas%cov_emis_lut)) DEALLOCATE(atlas%cov_emis_lut)
    IF ( ASSOCIATED(atlas%cov_emis)    ) DEALLOCATE(atlas%cov_emis)
    IF ( ASSOCIATED(atlas%pcu_int_v12) ) DEALLOCATE(atlas%pcu_int_v12)
    IF ( ASSOCIATED(atlas%pcm_int_v12) ) DEALLOCATE(atlas%pcm_int_v12)
    IF ( ASSOCIATED(atlas%pcu_int_v10) ) DEALLOCATE(atlas%pcu_int_v10)
    IF ( ASSOCIATED(atlas%pcm_int_v10) ) DEALLOCATE(atlas%pcm_int_v10)
    IF ( ASSOCIATED(atlas%pcu_int_v8)  ) DEALLOCATE(atlas%pcu_int_v8)
    IF ( ASSOCIATED(atlas%pcm_int_v8)  ) DEALLOCATE(atlas%pcm_int_v8)
    IF ( ASSOCIATED(atlas%sice_em_int) ) DEALLOCATE(atlas%sice_em_int)
    IF ( ASSOCIATED(atlas%snow_em_int) ) DEALLOCATE(atlas%snow_em_int)
    IF ( ASSOCIATED(atlas%cov_emis_int)) DEALLOCATE(atlas%cov_emis_int)

    IF ( ASSOCIATED(atlas%igbp)        ) DEALLOCATE(atlas%igbp)
    IF ( ASSOCIATED(atlas%p1d)         ) DEALLOCATE(atlas%p1d)
    IF ( ASSOCIATED(atlas%p2d)         ) DEALLOCATE(atlas%p2d)
    IF ( ASSOCIATED(atlas%p3d)         ) DEALLOCATE(atlas%p3d)
    IF ( ASSOCIATED(atlas%p1n)         ) DEALLOCATE(atlas%p1n)
    IF ( ASSOCIATED(atlas%p2n)         ) DEALLOCATE(atlas%p2n)
    IF ( ASSOCIATED(atlas%p3n)         ) DEALLOCATE(atlas%p3n)
    IF ( ASSOCIATED(atlas%p1d_int)     ) DEALLOCATE(atlas%p1d_int)
    IF ( ASSOCIATED(atlas%p2d_int)     ) DEALLOCATE(atlas%p2d_int)
    IF ( ASSOCIATED(atlas%p3d_int)     ) DEALLOCATE(atlas%p3d_int)
    IF ( ASSOCIATED(atlas%p1n_int)     ) DEALLOCATE(atlas%p1n_int)
    IF ( ASSOCIATED(atlas%p2n_int)     ) DEALLOCATE(atlas%p2n_int)
    IF ( ASSOCIATED(atlas%p3n_int)     ) DEALLOCATE(atlas%p3n_int)

    CALL rttov_camel_nullify_pointers(atlas)

  END SUBROUTINE rttov_camel_close_atlas

  !> Nullify pointers in CAMEL atlas data structure
  !! @param[in,out]   atlas   CAMEL atlas data structure to nullify
  SUBROUTINE rttov_camel_nullify_pointers(atlas)
    TYPE(camel_atlas_data), INTENT(INOUT) :: atlas

    NULLIFY(atlas%camel_qflag)
    NULLIFY(atlas%pc_labvs)
    NULLIFY(atlas%pc_npcs)
    NULLIFY(atlas%coef_lut)
    NULLIFY(atlas%pca_coef)
    NULLIFY(atlas%pcm_hsr_v8)
    NULLIFY(atlas%pcm_hsr_v10)
    NULLIFY(atlas%pcm_hsr_v12)
    NULLIFY(atlas%pcu_hsr_v8)
    NULLIFY(atlas%pcu_hsr_v10)
    NULLIFY(atlas%pcu_hsr_v12)
    NULLIFY(atlas%cov_emis_lut)
    NULLIFY(atlas%cov_emis)
    NULLIFY(atlas%pcu_int_v12)
    NULLIFY(atlas%pcm_int_v12)
    NULLIFY(atlas%pcu_int_v10)
    NULLIFY(atlas%pcm_int_v10)
    NULLIFY(atlas%pcu_int_v8)
    NULLIFY(atlas%pcm_int_v8)
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

  END SUBROUTINE rttov_camel_nullify_pointers

END MODULE mod_camel_atlas
