! Description:
!> @file
!!   Read HTFRTC coefficient files into RTTOV coefficients structure.
!
!> @brief
!!   Read HTFRTC coefficient files into RTTOV coefficients structure.
!!
!! @details
!!   This subroutine reads HTFRTC coefficient data. HTFRTC requires two files:
!!   one static file which is used for all HTFRTC simulations, and one sensor-
!!   specific file. Both files are read by this subroutine.
!!
!!
!! @param[out]    err                status on exit
!! @param[in,out] coefs              RTTOV coefs structure
!! @param[in]     fname_coef         file name of static HTFRTC data
!! @param[in]     fname_sensor       file name of sensor-specific HTFRTC data
!! @param[in]     channels_rec       list of channels for which PC reconstructed radiances are required, optional
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
!    Copyright 2018, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_read_coefs_htfrtc(err, coefs, fname_coef, fname_sensor, channels_rec)

!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY: jpim
  USE rttov_types, ONLY : rttov_coefs
!INTF_OFF
USE parkind1, ONLY: jprb, jplm
#ifdef _RTTOV_NETCDF
  USE netcdf
#endif
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim)     , INTENT(OUT)          :: err
  TYPE(rttov_coefs) , INTENT(INOUT)        :: coefs
  CHARACTER(LEN=*)  , INTENT(IN)           :: fname_coef
  CHARACTER(LEN=*)  , INTENT(IN)           :: fname_sensor
  INTEGER(jpim)     , INTENT(IN), OPTIONAL :: channels_rec(:)
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_nullify_coefs.interface"

#ifdef _RTTOV_NETCDF
  INTEGER(jpim) :: dimid
  INTEGER(jpim) :: varid
#endif
  INTEGER(jpim) :: lun
  INTEGER(jpim) :: i_ch, j_ch
  CHARACTER(LEN=30)       :: aline
  INTEGER(jpim)           :: tmp_n_ch
  REAL(jprb), ALLOCATABLE :: tmp_sensor_freq(:)
  REAL(jprb), ALLOCATABLE :: tmp_ch_mean(:)
  REAL(jprb), ALLOCATABLE :: tmp_pc(:,:)
  INTEGER(jpim), ALLOCATABLE :: tmp_addch(:,:)
  INTEGER :: ios
  LOGICAL(jplm):: lopen

  TRY

  CALL rttov_nullify_coefs(coefs)

  ! Read coefficient file
  DO lun = 9, 99
    INQUIRE(lun, opened = lopen)
    IF (.NOT. lopen) EXIT
  ENDDO

  OPEN(lun, file=TRIM(fname_coef), iostat=err, status='old')
  THROW(err.NE.0)
  READ(lun, '(a30)', iostat=ios) aline(1:30)
  IF ((ios /= 0) .OR. aline(1:30) /= '*dim version                 :') THEN
    ! Not ASCII, assume NETCDF
    CLOSE(lun)

#ifndef _RTTOV_NETCDF
    err = errorstatus_fatal
    THROWM(err.NE.0,"HTFRTC requires RTTOV to be compiled against NetCDF")
#else
    err=NF90_OPEN(TRIM(fname_coef), NF90_NOWRITE, lun)
    THROW(err.NE.0)

    err=NF90_INQ_DIMID(lun, "linear_gas", dimid)
    THROW(err.NE.0)
    err=NF90_INQUIRE_DIMENSION(lun, dimid, LEN=coefs%coef_htfrtc%n_gas_l)
    THROW(err.NE.0)

    err=NF90_INQ_DIMID(lun, "pressure", dimid)
    THROW(err.NE.0)
    err=NF90_INQUIRE_DIMENSION(lun, dimid, LEN=coefs%coef_htfrtc%n_p)
    THROW(err.NE.0)

    err=NF90_INQ_DIMID(lun, "linear_coefficients", dimid)
    THROW(err.NE.0)
    err=NF90_INQUIRE_DIMENSION(lun, dimid, LEN=coefs%coef_htfrtc%n_val_l)
    THROW(err.NE.0)

    err=NF90_INQ_DIMID(lun, "planck", dimid)
    THROW(err.NE.0)
    err=NF90_INQUIRE_DIMENSION(lun, dimid, LEN=coefs%coef_htfrtc%n_b)
    THROW(err.NE.0)

    err=NF90_INQ_DIMID(lun, "wavenumber", dimid)
    THROW(err.NE.0)
    err=NF90_INQUIRE_DIMENSION(lun, dimid, LEN=coefs%coef_htfrtc%n_f)
    THROW(err.NE.0)

    err=NF90_INQ_DIMID(lun, "lintau", dimid)
    THROW(err.NE.0)
    err=NF90_INQUIRE_DIMENSION(lun, dimid, LEN=coefs%coef_htfrtc%n_lt)
    THROW(err.NE.0)

    err=NF90_INQ_DIMID(lun, "seasurfemparams", dimid)
    THROW(err.NE.0)
    err=NF90_INQUIRE_DIMENSION(lun, dimid, LEN=coefs%coef_htfrtc%n_ssemp)
    THROW(err.NE.0)

    err=NF90_INQ_DIMID(lun, "iremisparams", dimid)
    THROW(err.NE.0)
    err=NF90_INQUIRE_DIMENSION(lun, dimid, LEN=coefs%coef_htfrtc%n_iremis)
    THROW(err.NE.0)

    err=NF90_INQ_DIMID(lun, "pcs", dimid)
    THROW(err.NE.0)
    err=NF90_INQUIRE_DIMENSION(lun, dimid, LEN=coefs%coef_htfrtc%n_pc)
    THROW(err.NE.0)

    err=NF90_INQ_DIMID(lun, "pcs_oc", dimid)
    THROW(err.NE.0)
    err=NF90_INQUIRE_DIMENSION(lun, dimid, LEN=coefs%coef_htfrtc%n_pc_oc)
    THROW(err.NE.0)

    err=NF90_INQ_DIMID(lun, "mftlb", dimid)
    THROW(err.NE.0)
    err=NF90_INQUIRE_DIMENSION(lun, dimid, LEN=coefs%coef_htfrtc%n_mftlb)
    THROW(err.NE.0)

    err=NF90_INQ_DIMID(lun, "cont", dimid)
    THROW(err.NE.0)
    err=NF90_INQUIRE_DIMENSION(lun, dimid, LEN=coefs%coef_htfrtc%n_cont)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%p(coefs%coef_htfrtc%n_p))
    err=NF90_INQ_VARID(lun, "pressure", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%p)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%coef_l(coefs%coef_htfrtc%n_val_l, &
    coefs%coef_htfrtc%n_p,coefs%coef_htfrtc%n_gas_l, &
    coefs%coef_htfrtc%n_f))
    err=NF90_INQ_VARID(lun, "linear_mass_ext_coeff", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%coef_l)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%mixed_ref_frac(coefs%coef_htfrtc%n_p, &
    coefs%coef_htfrtc%n_gas_l))
    err=NF90_INQ_VARID(lun, "mixed_ref_fraction", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%mixed_ref_frac)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%val_b(coefs%coef_htfrtc%n_b))
    err=NF90_INQ_VARID(lun, "planck_temp", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%val_b)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%coef_b(coefs%coef_htfrtc%n_b, &
    coefs%coef_htfrtc%n_f))
    err=NF90_INQ_VARID(lun, "planck", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%coef_b)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%coef_ct(coefs%coef_htfrtc%n_cont, &
    coefs%coef_htfrtc%n_f))
    err=NF90_INQ_VARID(lun, "coef_ct", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%coef_ct)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%coef_ctt(coefs%coef_htfrtc%n_cont, &
    coefs%coef_htfrtc%n_b,coefs%coef_htfrtc%n_f))
    err=NF90_INQ_VARID(lun, "coef_ctt", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%coef_ctt)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%val_lt(coefs%coef_htfrtc%n_lt))
    err=NF90_INQ_VARID(lun, "val_lintau", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%val_lt)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%coef_lt(coefs%coef_htfrtc%n_lt))
    err=NF90_INQ_VARID(lun, "lintau", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%coef_lt)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%coef_ssemp(coefs%coef_htfrtc%n_ssemp, &
    coefs%coef_htfrtc%n_f))
    err=NF90_INQ_VARID(lun, "ssemp", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%coef_ssemp)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%coef_iremis(coefs%coef_htfrtc%n_iremis, &
    coefs%coef_htfrtc%n_f))
    err=NF90_INQ_VARID(lun, "iremis", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%coef_iremis)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%val_mean(coefs%coef_htfrtc%n_f))
    err=NF90_INQ_VARID(lun, "mean", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%val_mean)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%val_norm(coefs%coef_htfrtc%n_f))
    err=NF90_INQ_VARID(lun, "norm", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%val_norm)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%freq(coefs%coef_htfrtc%n_f))
    err=NF90_INQ_VARID(lun, "wavenumber", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%freq)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%mftlb(coefs%coef_htfrtc%n_mftlb))
    err=NF90_INQ_VARID(lun, "mftlb", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%mftlb)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%coef_pdt(coefs%coef_htfrtc%n_f, &
    coefs%coef_htfrtc%n_pc))
    err=NF90_INQ_VARID(lun, "pdt", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%coef_pdt)
    THROW(err.NE.0)

    ALLOCATE(coefs%coef_htfrtc%addf(coefs%coef_htfrtc%n_gas_l, &
    coefs%coef_htfrtc%n_f))
    err=NF90_INQ_VARID(lun, "addf", varid)
    THROW(err.NE.0)
    err=NF90_GET_VAR(lun, varid, coefs%coef_htfrtc%addf)
    THROW(err.NE.0)

    err=NF90_CLOSE(lun)
    THROW(err.NE.0)
#endif
  ELSE
    DO
      READ(lun, '(a30)', iostat=ios) aline(1:30)
      IF (ios < 0) EXIT
      SELECT CASE (aline(1:30))

        CASE ('*dim linear_gas              :')
          READ(lun,*,iostat=err) coefs%coef_htfrtc%n_gas_l
          THROW(err.NE.0)

        CASE ('*dim pressure                :')
          READ(lun,*,iostat=err) coefs%coef_htfrtc%n_p
          THROW(err.NE.0)

        CASE ('*dim linear_coefficients     :')
          READ(lun,*,iostat=err) coefs%coef_htfrtc%n_val_l
          THROW(err.NE.0)

        CASE ('*dim planck                  :')
          READ(lun,*,iostat=err) coefs%coef_htfrtc%n_b
          THROW(err.NE.0)

        CASE ('*dim wavenumber              :')
          READ(lun,*,iostat=err) coefs%coef_htfrtc%n_f
          THROW(err.NE.0)

        CASE ('*dim lintau                  :')
          READ(lun,*,iostat=err) coefs%coef_htfrtc%n_lt
          THROW(err.NE.0)

        CASE ('*dim seasurfemparams         :')
          READ(lun,*,iostat=err) coefs%coef_htfrtc%n_ssemp
          THROW(err.NE.0)

        CASE ('*dim iremisparams            :')
          READ(lun,*,iostat=err) coefs%coef_htfrtc%n_iremis
          THROW(err.NE.0)

        CASE ('*dim pcs                     :')
          READ(lun,*,iostat=err) coefs%coef_htfrtc%n_pc
          THROW(err.NE.0)

        CASE ('*dim pcs_oc                  :')
          READ(lun,*,iostat=err) coefs%coef_htfrtc%n_pc_oc
          THROW(err.NE.0)

        CASE ('*dim mftlb                   :')
          READ(lun,*,iostat=err) coefs%coef_htfrtc%n_mftlb
          THROW(err.NE.0)

        CASE ('*dim cont                    :')
          READ(lun,*,iostat=err) coefs%coef_htfrtc%n_cont
          THROW(err.NE.0)

        CASE ('*var pressure                :')
          ALLOCATE(coefs%coef_htfrtc%p(coefs%coef_htfrtc%n_p))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%p
          THROW(err.NE.0)

        CASE ('*var mixed_ref_fraction      :')
          ALLOCATE(coefs%coef_htfrtc%mixed_ref_frac(coefs%coef_htfrtc%n_p, &
          coefs%coef_htfrtc%n_gas_l))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%mixed_ref_frac
          THROW(err.NE.0)

        CASE ('*var planck_temp             :')
          ALLOCATE(coefs%coef_htfrtc%val_b(coefs%coef_htfrtc%n_b))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%val_b
          THROW(err.NE.0)

        CASE ('*var linear_mass_ext_coeff   :')
          ALLOCATE(coefs%coef_htfrtc%coef_l(coefs%coef_htfrtc%n_val_l, &
          coefs%coef_htfrtc%n_p,coefs%coef_htfrtc%n_gas_l, &
          coefs%coef_htfrtc%n_f))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%coef_l
          THROW(err.NE.0)

        CASE ('*var planck                  :')
          ALLOCATE(coefs%coef_htfrtc%coef_b(coefs%coef_htfrtc%n_b, &
          coefs%coef_htfrtc%n_f))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%coef_b
          THROW(err.NE.0)

        CASE ('*var coef_ct                 :')
          ALLOCATE(coefs%coef_htfrtc%coef_ct(2,coefs%coef_htfrtc%n_f))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%coef_ct
          THROW(err.NE.0)

        CASE ('*var coef_ctt                :')
          ALLOCATE(coefs%coef_htfrtc%coef_ctt(2,coefs%coef_htfrtc%n_b, &
          coefs%coef_htfrtc%n_f))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%coef_ctt
          THROW(err.NE.0)

        CASE ('*var lintau                  :')
          ALLOCATE(coefs%coef_htfrtc%coef_lt(coefs%coef_htfrtc%n_lt))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%coef_lt
          THROW(err.NE.0)

        CASE ('*var val_lintau              :')
          ALLOCATE(coefs%coef_htfrtc%val_lt(coefs%coef_htfrtc%n_lt))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%val_lt
          THROW(err.NE.0)

        CASE ('*var ssemp                   :')
          ALLOCATE(coefs%coef_htfrtc%coef_ssemp(coefs%coef_htfrtc%n_ssemp, &
          coefs%coef_htfrtc%n_f))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%coef_ssemp
          THROW(err.NE.0)

        CASE ('*var iremis                  :')
          ALLOCATE(coefs%coef_htfrtc%coef_iremis(coefs%coef_htfrtc%n_iremis, &
          coefs%coef_htfrtc%n_f))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%coef_iremis
          THROW(err.NE.0)

        CASE ('*var mean                    :')
          ALLOCATE(coefs%coef_htfrtc%val_mean(coefs%coef_htfrtc%n_f))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%val_mean
          THROW(err.NE.0)

        CASE ('*var norm                    :')
          ALLOCATE(coefs%coef_htfrtc%val_norm(coefs%coef_htfrtc%n_f))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%val_norm
          THROW(err.NE.0)

        CASE ('*var wavenumber              :')
          ALLOCATE(coefs%coef_htfrtc%freq(coefs%coef_htfrtc%n_f))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%freq
          THROW(err.NE.0)

        CASE ('*var mftlb                   :')
          ALLOCATE(coefs%coef_htfrtc%mftlb(coefs%coef_htfrtc%n_mftlb))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%mftlb
          THROW(err.NE.0)

        CASE ('*var pdt                     :')
          ALLOCATE(coefs%coef_htfrtc%coef_pdt(coefs%coef_htfrtc%n_f, &
          coefs%coef_htfrtc%n_pc))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%coef_pdt
          THROW(err.NE.0)

        CASE ('*var addf                    :')
          ALLOCATE(coefs%coef_htfrtc%addf(coefs%coef_htfrtc%n_gas_l, &
          coefs%coef_htfrtc%n_f))
          READ(lun,*,iostat=err) coefs%coef_htfrtc%addf
          THROW(err.NE.0)

      END SELECT
    END DO

    CLOSE(lun)

  ENDIF


  ! Read sensor file
  DO lun = 9, 99
    INQUIRE(lun, opened = lopen)
    IF (.NOT. lopen) EXIT
  ENDDO

  OPEN(lun, file=TRIM(fname_sensor), iostat=err, status='old')
  THROW(err.NE.0)
  READ(lun, '(a30)', iostat=ios) aline(1:30)
  IF ((ios /= 0) .OR. aline(1:30) /= '*dim version                 :') THEN
    ! Not ASCII, assume NETCDF
    CLOSE(lun)

#ifndef _RTTOV_NETCDF
    err = errorstatus_fatal
    THROWM(err.NE.0,"HTFRTC requires RTTOV to be compiled against NetCDF")
#else
    err=NF90_OPEN(TRIM(fname_sensor), NF90_NOWRITE, lun)
    THROW(err.NE.0)

    IF (PRESENT(channels_rec)) THEN
      coefs%coef_htfrtc%n_ch=size(channels_rec)
      ALLOCATE(coefs%coef_htfrtc%sensor_freq(coefs%coef_htfrtc%n_ch))
      err=NF90_INQ_VARID(lun, "wavenumber", varid)
      THROW(err.NE.0)
      DO i_ch=1,coefs%coef_htfrtc%n_ch
        j_ch=channels_rec(i_ch)
        err=NF90_GET_VAR(lun, varid, &
        coefs%coef_htfrtc%sensor_freq(i_ch), start=(/j_ch/))
        THROW(err.NE.0)
      ENDDO
      ALLOCATE(coefs%coef_htfrtc%ch_mean(coefs%coef_htfrtc%n_ch))
      err=NF90_INQ_VARID(lun, "mean", varid)
      THROW(err.NE.0)
      DO i_ch=1,coefs%coef_htfrtc%n_ch
        j_ch=channels_rec(i_ch)
        err=NF90_GET_VAR(lun, varid, &
        coefs%coef_htfrtc%ch_mean(i_ch), start=(/j_ch/))
        THROW(err.NE.0)
      ENDDO
      ALLOCATE(coefs%coef_htfrtc%pc(coefs%coef_htfrtc%n_pc, &
      coefs%coef_htfrtc%n_ch))
      err=NF90_INQ_VARID(lun, "eof", varid)
      THROW(err.NE.0)
      DO i_ch=1,coefs%coef_htfrtc%n_ch
        j_ch=channels_rec(i_ch)
        err=NF90_GET_VAR(lun, varid, &
        coefs%coef_htfrtc%pc(1:coefs%coef_htfrtc%n_pc,i_ch), start=(/1,j_ch/))
        THROW(err.NE.0)
      ENDDO
      ALLOCATE(coefs%coef_htfrtc%addch(coefs%coef_htfrtc%n_gas_l, &
      coefs%coef_htfrtc%n_ch))
      err=NF90_INQ_VARID(lun, "addch", varid)
      THROW(err.NE.0)
      DO i_ch=1,coefs%coef_htfrtc%n_ch
        j_ch=channels_rec(i_ch)
        err=NF90_GET_VAR(lun, varid, &
        coefs%coef_htfrtc%addch(1:coefs%coef_htfrtc%n_gas_l,i_ch), start=(/1,j_ch/))
        THROW(err.NE.0)
      ENDDO
    ELSE
      err=NF90_INQ_DIMID(lun, "wavenumber", dimid)
      THROW(err.NE.0)
      err=NF90_INQUIRE_DIMENSION(lun, dimid, LEN=coefs%coef_htfrtc%n_ch)
      THROW(err.NE.0)

      ALLOCATE(coefs%coef_htfrtc%sensor_freq(coefs%coef_htfrtc%n_ch))
      err=NF90_INQ_VARID(lun, "wavenumber", varid)
      THROW(err.NE.0)
      DO i_ch=1,coefs%coef_htfrtc%n_ch
        err=NF90_GET_VAR(lun, varid, &
        coefs%coef_htfrtc%sensor_freq(i_ch), start=(/i_ch/))
        THROW(err.NE.0)
      ENDDO
      ALLOCATE(coefs%coef_htfrtc%ch_mean(coefs%coef_htfrtc%n_ch))
      err=NF90_INQ_VARID(lun, "mean", varid)
      THROW(err.NE.0)
      DO i_ch=1,coefs%coef_htfrtc%n_ch
        err=NF90_GET_VAR(lun, varid, &
        coefs%coef_htfrtc%ch_mean(i_ch), start=(/i_ch/))
        THROW(err.NE.0)
      ENDDO
      ALLOCATE(coefs%coef_htfrtc%pc(coefs%coef_htfrtc%n_pc, &
      coefs%coef_htfrtc%n_ch))
      err=NF90_INQ_VARID(lun, "eof", varid)
      THROW(err.NE.0)
      DO i_ch=1,coefs%coef_htfrtc%n_ch
        err=NF90_GET_VAR(lun, varid, &
        coefs%coef_htfrtc%pc(1:coefs%coef_htfrtc%n_pc,i_ch), start=(/1,i_ch/))
        THROW(err.NE.0)
      ENDDO
      ALLOCATE(coefs%coef_htfrtc%addch(coefs%coef_htfrtc%n_gas_l, &
      coefs%coef_htfrtc%n_ch))
      err=NF90_INQ_VARID(lun, "addch", varid)
      THROW(err.NE.0)
      DO i_ch=1,coefs%coef_htfrtc%n_ch
        err=NF90_GET_VAR(lun, varid, &
        coefs%coef_htfrtc%addch(1:coefs%coef_htfrtc%n_gas_l,i_ch), start=(/1,i_ch/))
        THROW(err.NE.0)
      ENDDO
    ENDIF

    err=NF90_CLOSE(lun)
    THROW(err.NE.0)
#endif
  ELSE
    DO
      READ(lun, '(a30)', iostat=ios) aline(1:30)
      IF (ios < 0) EXIT
      SELECT CASE (aline(1:30))

        CASE ('*dim wavenumber              :')
          READ(lun,*) tmp_n_ch

        CASE ('*var wavenumber              :')
          ALLOCATE(tmp_sensor_freq(tmp_n_ch))
          READ(lun,*) tmp_sensor_freq

        CASE ('*var mean                    :')
          ALLOCATE(tmp_ch_mean(tmp_n_ch))
          READ(lun,*) tmp_ch_mean

        CASE ('*var eof                     :')
          ALLOCATE(tmp_pc(coefs%coef_htfrtc%n_pc,tmp_n_ch))
          READ(lun,*) tmp_pc

        CASE ('*var addch                   :')
          ALLOCATE(tmp_addch(coefs%coef_htfrtc%n_gas_l,tmp_n_ch))
          READ(lun,*) tmp_addch

      END SELECT
    END DO

    CLOSE(lun)

    IF (PRESENT(channels_rec)) THEN

      coefs%coef_htfrtc%n_ch=size(channels_rec)
      ALLOCATE(coefs%coef_htfrtc%sensor_freq(coefs%coef_htfrtc%n_ch))
      DO i_ch=1,coefs%coef_htfrtc%n_ch
        j_ch=channels_rec(i_ch)
        coefs%coef_htfrtc%sensor_freq(i_ch)=tmp_sensor_freq(j_ch)
      ENDDO
      ALLOCATE(coefs%coef_htfrtc%ch_mean(coefs%coef_htfrtc%n_ch))
      DO i_ch=1,coefs%coef_htfrtc%n_ch
        j_ch=channels_rec(i_ch)
        coefs%coef_htfrtc%ch_mean(i_ch)=tmp_ch_mean(j_ch)
      ENDDO
      ALLOCATE(coefs%coef_htfrtc%pc(coefs%coef_htfrtc%n_pc,coefs%coef_htfrtc%n_ch))
      DO i_ch=1,coefs%coef_htfrtc%n_ch
        j_ch=channels_rec(i_ch)
        coefs%coef_htfrtc%pc(1:coefs%coef_htfrtc%n_pc,i_ch)= &
        tmp_pc(1:coefs%coef_htfrtc%n_pc,j_ch)
      ENDDO
      ALLOCATE(coefs%coef_htfrtc%addch(coefs%coef_htfrtc%n_gas_l,coefs%coef_htfrtc%n_ch))
      DO i_ch=1,coefs%coef_htfrtc%n_ch
        j_ch=channels_rec(i_ch)
        coefs%coef_htfrtc%addch(1:coefs%coef_htfrtc%n_gas_l,i_ch)= &
        tmp_addch(1:coefs%coef_htfrtc%n_gas_l,j_ch)
      ENDDO

    ELSE

      coefs%coef_htfrtc%n_ch=tmp_n_ch
      ALLOCATE(coefs%coef_htfrtc%sensor_freq(coefs%coef_htfrtc%n_ch))
      DO i_ch=1,coefs%coef_htfrtc%n_ch
       coefs%coef_htfrtc%sensor_freq(i_ch)=tmp_sensor_freq(i_ch)
      ENDDO
      ALLOCATE(coefs%coef_htfrtc%ch_mean(coefs%coef_htfrtc%n_ch))
      DO i_ch=1,coefs%coef_htfrtc%n_ch
        coefs%coef_htfrtc%ch_mean(i_ch)=tmp_ch_mean(i_ch)
      ENDDO
      ALLOCATE(coefs%coef_htfrtc%pc(coefs%coef_htfrtc%n_pc,coefs%coef_htfrtc%n_ch))
      DO i_ch=1,coefs%coef_htfrtc%n_ch
        coefs%coef_htfrtc%pc(1:coefs%coef_htfrtc%n_pc,i_ch)= &
        tmp_pc(1:coefs%coef_htfrtc%n_pc,i_ch)
      ENDDO
      ALLOCATE(coefs%coef_htfrtc%addch(coefs%coef_htfrtc%n_gas_l,coefs%coef_htfrtc%n_ch))
      DO i_ch=1,coefs%coef_htfrtc%n_ch
        coefs%coef_htfrtc%addch(1:coefs%coef_htfrtc%n_gas_l,i_ch)= &
        tmp_addch(1:coefs%coef_htfrtc%n_gas_l,i_ch)
      ENDDO

    ENDIF

    IF (ALLOCATED(tmp_sensor_freq)) DEALLOCATE(tmp_sensor_freq)
    IF (ALLOCATED(tmp_ch_mean)) DEALLOCATE(tmp_ch_mean)
    IF (ALLOCATED(tmp_pc)) DEALLOCATE(tmp_pc)
    IF (ALLOCATED(tmp_addch)) DEALLOCATE(tmp_addch)

  ENDIF

  CATCH
END SUBROUTINE rttov_read_coefs_htfrtc
