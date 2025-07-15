! Description:
!> @file
!!   Read a v10/v11 HDF5 coefficient file into the v12 rttov_coef structure.
!!
!> @brief
!!   Read a v10/v11 HDF5 coefficient file into the v12 rttov_coef structure.
!!
!> @details
!!   There is no optional channel selection: do this before or after conversion.
!!
!! @param[out]    err             status on exit
!! @param[in]     file_coef       full path of input file
!! @param[in,out] coef            RTTOV optical depth coefficient structure
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
SUBROUTINE rttov11_read_hdf5_coef(err, file_coef, coef)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE parkind1, ONLY : jprb, jplm
#ifdef _RTTOV_HDF
  USE rttov_types, ONLY : rttov_fast_coef_hdf_io
  USE rttov_const, ONLY :  &
        gas_mass,               &
        mair,                   &
        gas_unit_specconc,      &
        gas_id_mixed,           &
        gas_id_watervapour,     &
        gas_id_ozone,           &
        gas_id_wvcont,          &
        gas_id_co2,             &
        gas_id_n2o,             &
        gas_id_co,              &
        gas_id_ch4,             &
        ngases_max,             &
        tscale_v11,             &
        gscale_v11
  USE rttov_fast_coef_utils_mod, ONLY : set_pointers
  USE hdf5
  USE rttov_hdf_mod
  USE rttov_hdf_rttov_fast_coef_io
  USE rttov_hdf_rttov_nlte_coef_io
#endif
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT)   :: err
  CHARACTER(LEN=*),   INTENT(IN)    :: file_coef
  TYPE(rttov_coef),   INTENT(INOUT) :: coef
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_nullify_coef.interface"

#ifdef _RTTOV_HDF
  INTEGER(HID_T)       :: file_id, g_id, g_id_sub
  INTEGER(KIND=jpim)   :: n, i
  CHARACTER(LEN=lensh) :: sname
  LOGICAL(KIND=jplm)   :: lext
  INTEGER(KIND=jpim), POINTER  :: tmp_gas_pos(:)
  INTEGER(KIND=jpim), POINTER  :: gaz_units(:)
  TYPE(rttov_fast_coef_hdf_io) :: FAST_COEF_temp
#endif

!- End of header --------------------------------------------------------
  TRY

#ifndef _RTTOV_HDF
  err = errorstatus_fatal
  THROWM(err.NE.0, "")
#else

  CALL H5FOPEN_F( FILE_COEF, H5F_ACC_RDONLY_F, FILE_ID, ERR )
  THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(FILE_COEF))

  CALL H5GOPEN_F( FILE_ID, "/COEF", G_ID, ERR )
  THROWM(ERR.NE.0,"CANNOT OPEN GROUP /COEF IN FILE "//TRIM(FILE_COEF))


  CALL rttov_nullify_coef(coef)

  coef%iremis_version = 0

  sname='ID_PLATFORM'
  CALL read_array_hdf(g_id,sname,err,i0=coef%ID_PLATFORM)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='ID_SAT'
  CALL read_array_hdf(g_id,sname,err,i0=coef%ID_SAT)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='ID_INST'
  CALL read_array_hdf(g_id,sname,err,i0=coef%ID_INST)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='ID_SENSOR'
  CALL read_array_hdf(g_id,sname,err,i0=coef%ID_SENSOR)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='ID_COMP_LVL'
  CALL read_array_hdf(g_id,sname,err,i0=coef%ID_COMP_LVL)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  IF (coef%id_comp_lvl < 10 .OR. coef%id_comp_lvl > 11) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, "Version of coefficient file is incompatible with this subroutine (must be RTTOV v10/v11)")
  ENDIF
  coef%id_comp_lvl = 12

  sname='ID_COMP_PC'
  CALL read_array_hdf(g_id,sname,err,i0=coef%ID_COMP_PC)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='ID_CREATION_DATE'
  CALL read_array_hdf(g_id,sname,err,i1=coef%ID_CREATION_DATE)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='ID_CREATION'
  CALL read_array_hdf(g_id,sname,err,c0=coef%ID_CREATION)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='ID_COMMON_NAME'
  CALL read_array_hdf(g_id,sname,err,c0=coef%ID_COMMON_NAME)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='LINE_BY_LINE'
  CALL read_array_hdf(g_id,sname,err,c1=coef%LINE_BY_LINE)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='FMV_MODEL_DEF'
  CALL read_array_hdf(g_id,sname,err,c0=coef%FMV_MODEL_DEF)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='FMV_MODEL_VER'
  CALL read_array_hdf(g_id,sname,err,i0=coef%FMV_MODEL_VER)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='FMV_CHN'
  CALL read_array_hdf(g_id,sname,err,i0=coef%FMV_CHN)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='FMV_GAS'
  CALL read_array_hdf(g_id,sname,err,i0=coef%FMV_GAS)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='FMV_GAS_ID'
  CALL read_array_hdf(g_id,sname,err,pi1=coef%FMV_GAS_ID)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='FMV_GAS_POS'
  CALL read_array_hdf(g_id,sname,err,pi1=tmp_gas_pos)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ALLOCATE(coef%FMV_GAS_POS(ngases_max))   ! ngases_max = 8 in RTTOV v11 but is 9 in RTTOV v12
  coef%FMV_GAS_POS = 0
  coef%FMV_GAS_POS(1:8) = tmp_gas_pos
  DEALLOCATE(tmp_gas_pos)

  sname='FMV_VAR'
  CALL read_array_hdf(g_id,sname,err,pi1=coef%FMV_VAR)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='FMV_COE'
  CALL read_array_hdf(g_id,sname,err,pi1=coef%FMV_COE)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='FMV_LVL'
  CALL read_array_hdf(g_id,sname,err,pi1=coef%FMV_LVL)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='INCZEEMAN'
  CALL read_array_hdf(g_id,sname,err,l0=coef%INCZEEMAN)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='FC_PLANCK_C1'
  CALL read_array_hdf(g_id,sname,err,r0=coef%FC_PLANCK_C1)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='FC_PLANCK_C2'
  CALL read_array_hdf(g_id,sname,err,r0=coef%FC_PLANCK_C2)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='FC_SAT_HEIGHT'
  CALL read_array_hdf(g_id,sname,err,r0=coef%FC_SAT_HEIGHT)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='GAZ_UNITS'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
    CALL read_array_hdf(g_id,sname,err,pi1=GAZ_UNITS)
    THROWM(err.ne.0,"CANNOT READ "//TRIM(sname))
  ELSE
    NULLIFY(gaz_units)
  ENDIF
  IF (.NOT. ASSOCIATED(gaz_units)) THEN
    ALLOCATE(gaz_units(coef%fmv_gas))
    gaz_units = 1  ! v11 value for kg/kg - the default if the section is missing
  ENDIF

  sname='README_SRF'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,c1=coef%README_SRF)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='FF_ORI_CHN'
  CALL read_array_hdf(g_id,sname,err,pi1=coef%FF_ORI_CHN)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='FF_VAL_CHN'
  CALL read_array_hdf(g_id,sname,err,pi1=coef%FF_VAL_CHN)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='FF_CWN'
  CALL read_array_hdf(g_id,sname,err,pr1=coef%FF_CWN)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='FF_BCO'
  CALL read_array_hdf(g_id,sname,err,pr1=coef%FF_BCO)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='FF_BCS'
  CALL read_array_hdf(g_id,sname,err,pr1=coef%FF_BCS)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='FF_GAM'
  CALL read_array_hdf(g_id,sname,err,pr1=coef%FF_GAM)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='TT_VAL_CHN'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,pi1=coef%TT_VAL_CHN)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='TT_A0'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,pr1=coef%TT_A0)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='TT_A1'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,pr1=coef%TT_A1)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='PW_VAL_CHN'
  CALL read_array_hdf(g_id,sname,err,pi1=coef%PW_VAL_CHN)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='SS_VAL_CHN'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,pi1=coef%SS_VAL_CHN)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='SS_SOLAR_SPECTRUM'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,pr1=coef%SS_SOLAR_SPECTRUM)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='WOC_WAOPC_OW'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf_cmplx(g_id,sname,err,pk1=coef%WOC_WAOPC_OW)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='WOC_WAOPC_FW'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf_cmplx(g_id,sname,err,pk1=coef%WOC_WAOPC_FW)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='WS_NOMEGA'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,i0=coef%WS_NOMEGA)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='WS_NPOINT'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,pr1=coef%WS_NPOINT)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='WS_K_OMEGA'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,pr1=coef%WS_K_OMEGA)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='FASTEM_POLAR'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,pi1=coef%FASTEM_POLAR)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='SSIREM_VER'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,i0=coef%SSIREM_VER)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='SSIREM_A0'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,pr1=coef%SSIREM_A0)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='SSIREM_A1'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,pr1=coef%SSIREM_A1)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='SSIREM_A2'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,pr1=coef%SSIREM_A2)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='SSIREM_XZN1'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,pr1=coef%SSIREM_XZN1)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='SSIREM_XZN2'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,pr1=coef%SSIREM_XZN2)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='REF_PRFL_P'
  CALL read_array_hdf(g_id,sname,err,pr1=coef%REF_PRFL_P)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='REF_PRFL_T'
  CALL read_array_hdf(g_id,sname,err,pr2=coef%REF_PRFL_T)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='REF_PRFL_MR'
  CALL read_array_hdf(g_id,sname,err,pr2=coef%REF_PRFL_MR)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='LIM_PRFL_P'
  CALL read_array_hdf(g_id,sname,err,pr1=coef%LIM_PRFL_P)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='LIM_PRFL_TMAX'
  CALL read_array_hdf(g_id,sname,err,pr1=coef%ENV_PRFL_TMAX)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='LIM_PRFL_TMIN'
  CALL read_array_hdf(g_id,sname,err,pr1=coef%ENV_PRFL_TMIN)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='LIM_PRFL_GMAX'
  CALL read_array_hdf(g_id,sname,err,pr2=coef%ENV_PRFL_GMAX)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))

  sname='LIM_PRFL_GMIN'
  CALL read_array_hdf(g_id,sname,err,pr2=coef%ENV_PRFL_GMIN)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))


  ! Convert limits read from file to profile envelope and take care of gas units (inc. background profile)
  ! Conversion of gas mixing ratio units (kg/kg wet -> ppmv dry) if needed
  DO n = 1, coef%fmv_gas
    IF (n /= gas_id_mixed) THEN
      IF (gaz_units(n) == gas_unit_specconc) THEN
        coef%ref_prfl_mr(:,n)   = (1.E06_jprb * Mair / gas_mass(coef%fmv_gas_id(n))) * &
                                  coef%ref_prfl_mr(:,n) / (1._jprb - coef%ref_prfl_mr(:,n))
        coef%env_prfl_gmin(:,n) = (1.E06_jprb * Mair / gas_mass(coef%fmv_gas_id(n))) * &
                                  coef%env_prfl_gmin(:,n) / (1._jprb - coef%env_prfl_gmin(:,n))
        coef%env_prfl_gmax(:,n) = (1.E06_jprb * Mair / gas_mass(coef%fmv_gas_id(n))) * &
                                  coef%env_prfl_gmax(:,n) / (1._jprb - coef%env_prfl_gmax(:,n))
      ENDIF
      coef%env_prfl_gmax(:,n) = coef%env_prfl_gmax(:,n) / (1._jprb + gscale_v11)
      coef%env_prfl_gmin(:,n) = coef%env_prfl_gmin(:,n) / (1._jprb - gscale_v11)
    ENDIF
  ENDDO

  coef%env_prfl_tmax = coef%env_prfl_tmax / (1._jprb + tscale_v11)
  coef%env_prfl_tmin = coef%env_prfl_tmin / (1._jprb - tscale_v11)

  ALLOCATE (coef%bkg_prfl_mr(coef%fmv_lvl(gas_id_mixed), coef%fmv_gas), STAT = ERR)
  THROWM(err.NE.0, "allocation of bkg_prfl_mr")
  coef%bkg_prfl_mr = coef%ref_prfl_mr

  IF (ASSOCIATED(gaz_units)) DEALLOCATE(gaz_units)



  sname='PMC_LENGTHCELL'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,r0=coef%PMC_LENGTHCELL)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='PMC_PNOMINAL'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,pr1=coef%PMC_PNOMINAL)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='PMC_TEMPCELL'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,r0=coef%PMC_TEMPCELL)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='PMC_BETAPLUS1'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,r0=coef%PMC_BETAPLUS1)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='PMC_NLAY'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,i0=coef%PMC_NLAY)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='PMC_NVAR'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,i0=coef%PMC_NVAR)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  sname='PMC_COEF'
  CALL h5lexists_f( g_id, sname, lext, err )
  IF (lext) THEN
  CALL read_array_hdf(g_id,sname,err,pr3=coef%PMC_COEF)
  THROWM(err.NE.0,"CANNOT READ "//TRIM(sname))
  ENDIF

  coef%fmv_ori_nchn = coef%fmv_chn
  coef%PMC_SHIFT = ASSOCIATED(coef%PMC_COEF)

  DO n = 1, coef%fmv_gas
    SELECT CASE (coef%fmv_gas_id(n))
    CASE (gas_id_mixed)
      coef%nmixed  = coef%fmv_var(n)
      coef%ncmixed = coef%fmv_coe(n)
      coef%nlevels = coef%fmv_lvl(n)
    CASE (gas_id_watervapour)
      coef%nwater  = coef%fmv_var(n)
      coef%ncwater = coef%fmv_coe(n)
    CASE (gas_id_ozone)
      coef%nozone  = coef%fmv_var(n)
      coef%ncozone = coef%fmv_coe(n)
    CASE (gas_id_wvcont)
      coef%nwvcont  = coef%fmv_var(n)
      coef%ncwvcont = coef%fmv_coe(n)
    CASE (gas_id_co2)
      coef%nco2  = coef%fmv_var(n)
      coef%ncco2 = coef%fmv_coe(n)
    CASE (gas_id_n2o)
      coef%nn2o  = coef%fmv_var(n)
      coef%ncn2o = coef%fmv_coe(n)
    CASE (gas_id_co)
      coef%nco  = coef%fmv_var(n)
      coef%ncco = coef%fmv_coe(n)
    CASE (gas_id_ch4)
      coef%nch4  = coef%fmv_var(n)
      coef%ncch4 = coef%fmv_coe(n)
    END SELECT
  END DO
  coef%NLAYERS = coef%NLEVELS - 1


  NULLIFY(FAST_COEF_temp%mixedgas)
  NULLIFY(FAST_COEF_temp%watervapour)
  NULLIFY(FAST_COEF_temp%ozone)
  NULLIFY(FAST_COEF_temp%wvcont)
  NULLIFY(FAST_COEF_temp%co2)
  NULLIFY(FAST_COEF_temp%n2o)
  NULLIFY(FAST_COEF_temp%co)
  NULLIFY(FAST_COEF_temp%ch4)

  CALL H5GOPEN_F( G_ID, "THERMAL", G_ID_SUB, ERR )
  THROWM(ERR.NE.0,"CANNOT OPEN GROUP THERMAL")

  CALL RTTOV_HDF_RTTOV_FAST_COEF_RH( FAST_COEF_temp, G_ID_SUB, ERR )
  THROWM(ERR.NE.0,"CANNOT READ COEF%THERMAL")

  ALLOCATE (coef%thermal(coef%fmv_chn), STAT = ERR)
  THROWM(ERR.NE.0, "allocation of thermal fast coefs")
  DO i = 1, coef%fmv_chn
    ALLOCATE (coef%thermal(i)%gasarray(coef%fmv_gas), STAT = ERR)
    THROWM(ERR.NE.0, "allocation of gasarray")
    NULLIFY (coef%thermal(i)%mixedgas)
    NULLIFY (coef%thermal(i)%watervapour)
    NULLIFY (coef%thermal(i)%ozone)
    NULLIFY (coef%thermal(i)%wvcont)
    NULLIFY (coef%thermal(i)%co2)
    NULLIFY (coef%thermal(i)%n2o)
    NULLIFY (coef%thermal(i)%co)
    NULLIFY (coef%thermal(i)%ch4)
  ENDDO

  IF (ASSOCIATED(FAST_COEF_temp%mixedgas)) THEN
    CALL move_hdf_fast_coef(coef, coef%thermal, FAST_COEF_temp%mixedgas, gas_id_mixed)
    DEALLOCATE(FAST_COEF_temp%mixedgas)    ; NULLIFY(FAST_COEF_temp%mixedgas)
  ENDIF
  IF (ASSOCIATED(FAST_COEF_temp%watervapour)) THEN
    CALL move_hdf_fast_coef(coef, coef%thermal, FAST_COEF_temp%watervapour, gas_id_watervapour)
    DEALLOCATE(FAST_COEF_temp%watervapour) ; NULLIFY(FAST_COEF_temp%watervapour)
  ENDIF
  IF (ASSOCIATED(FAST_COEF_temp%ozone)) THEN
    CALL move_hdf_fast_coef(coef, coef%thermal, FAST_COEF_temp%ozone, gas_id_ozone)
    DEALLOCATE(FAST_COEF_temp%ozone)       ; NULLIFY(FAST_COEF_temp%ozone)
  ENDIF
  IF (ASSOCIATED(FAST_COEF_temp%wvcont)) THEN
    CALL move_hdf_fast_coef(coef, coef%thermal, FAST_COEF_temp%wvcont, gas_id_wvcont)
    DEALLOCATE(FAST_COEF_temp%wvcont)      ; NULLIFY(FAST_COEF_temp%wvcont)
  ENDIF
  IF (ASSOCIATED(FAST_COEF_temp%co2)) THEN
    CALL move_hdf_fast_coef(coef, coef%thermal,FAST_COEF_temp%co2, gas_id_co2)
    DEALLOCATE(FAST_COEF_temp%co2)         ; NULLIFY(FAST_COEF_temp%co2)
  ENDIF
  IF (ASSOCIATED(FAST_COEF_temp%n2o)) THEN
    CALL move_hdf_fast_coef(coef, coef%thermal,FAST_COEF_temp%n2o, gas_id_n2o)
    DEALLOCATE(FAST_COEF_temp%n2o)         ; NULLIFY(FAST_COEF_temp%n2o)
  ENDIF
  IF (ASSOCIATED(FAST_COEF_temp%co)) THEN
    CALL move_hdf_fast_coef(coef, coef%thermal,FAST_COEF_temp%co, gas_id_co)
    DEALLOCATE(FAST_COEF_temp%co)          ; NULLIFY(FAST_COEF_temp%co)
  ENDIF
  IF (ASSOCIATED(FAST_COEF_temp%ch4)) THEN
    CALL move_hdf_fast_coef(coef, coef%thermal,FAST_COEF_temp%ch4, gas_id_ch4)
    DEALLOCATE(FAST_COEF_temp%ch4)         ; NULLIFY(FAST_COEF_temp%ch4)
  ENDIF

  CALL H5GCLOSE_F( G_ID_SUB, ERR )
  THROWM(ERR.NE.0,"CANNOT CLOSE GROUP THERMAL")


  coef%SOLARCOEF = .FALSE.
  NULLIFY( coef%SOLAR)

  CALL H5LEXISTS_F( G_ID, "SOLAR", LEXT, ERR )
  THROWM(ERR.NE.0,"CANNOT TEST SOLAR ")
  IF( LEXT ) THEN

    CALL H5GOPEN_F( G_ID, "SOLAR", G_ID_SUB, ERR )
    THROWM(ERR.NE.0,"CANNOT OPEN GROUP SOLAR")

    CALL RTTOV_HDF_RTTOV_FAST_COEF_RH( FAST_COEF_temp, G_ID_SUB, ERR )
    THROWM(ERR.NE.0,"CANNOT READ COEF%SOLAR")

    ALLOCATE (coef%solar(coef%fmv_chn), STAT = ERR)
    THROWM( ERR .NE. 0, "ALLOCATION OF SOLAR FAST COEFS")
    DO i = 1, coef%fmv_chn
      ALLOCATE (coef%solar(i)%gasarray(coef%fmv_gas), STAT = ERR)
      THROWM(ERR.NE.0, "allocation of gasarray")
      NULLIFY (coef%solar(i)%mixedgas)
      NULLIFY (coef%solar(i)%watervapour)
      NULLIFY (coef%solar(i)%ozone)
      NULLIFY (coef%solar(i)%wvcont)
      NULLIFY (coef%solar(i)%co2)
      NULLIFY (coef%solar(i)%n2o)
      NULLIFY (coef%solar(i)%co)
      NULLIFY (coef%solar(i)%ch4)
    ENDDO

    IF (ASSOCIATED(FAST_COEF_temp%mixedgas)) THEN
      CALL move_hdf_fast_coef(coef, coef%solar, FAST_COEF_temp%mixedgas, gas_id_mixed)
      DEALLOCATE(FAST_COEF_temp%mixedgas)    ; NULLIFY(FAST_COEF_temp%mixedgas)
    ENDIF
    IF (ASSOCIATED(FAST_COEF_temp%watervapour)) THEN
      CALL move_hdf_fast_coef(coef, coef%solar, FAST_COEF_temp%watervapour, gas_id_watervapour)
      DEALLOCATE(FAST_COEF_temp%watervapour) ; NULLIFY(FAST_COEF_temp%watervapour)
    ENDIF
    IF (ASSOCIATED(FAST_COEF_temp%ozone)) THEN
      CALL move_hdf_fast_coef(coef, coef%solar, FAST_COEF_temp%ozone, gas_id_ozone)
      DEALLOCATE(FAST_COEF_temp%ozone)       ; NULLIFY(FAST_COEF_temp%ozone)
    ENDIF
    IF (ASSOCIATED(FAST_COEF_temp%wvcont)) THEN
      CALL move_hdf_fast_coef(coef, coef%solar, FAST_COEF_temp%wvcont, gas_id_wvcont)
      DEALLOCATE(FAST_COEF_temp%wvcont)      ; NULLIFY(FAST_COEF_temp%wvcont)
    ENDIF
    IF (ASSOCIATED(FAST_COEF_temp%co2)) THEN
      CALL move_hdf_fast_coef(coef, coef%solar, FAST_COEF_temp%co2, gas_id_co2)
      DEALLOCATE(FAST_COEF_temp%co2)         ; NULLIFY(FAST_COEF_temp%co2)
    ENDIF
    IF (ASSOCIATED(FAST_COEF_temp%n2o)) THEN
      CALL move_hdf_fast_coef(coef, coef%solar, FAST_COEF_temp%n2o, gas_id_n2o)
      DEALLOCATE(FAST_COEF_temp%n2o)         ; NULLIFY(FAST_COEF_temp%n2o)
    ENDIF
    IF (ASSOCIATED(FAST_COEF_temp%co)) THEN
      CALL move_hdf_fast_coef(coef, coef%solar, FAST_COEF_temp%co, gas_id_co)
      DEALLOCATE(FAST_COEF_temp%co)          ; NULLIFY(FAST_COEF_temp%co)
    ENDIF
    IF (ASSOCIATED(FAST_COEF_temp%ch4)) THEN
      CALL move_hdf_fast_coef(coef, coef%solar, FAST_COEF_temp%ch4, gas_id_ch4)
      DEALLOCATE(FAST_COEF_temp%ch4)         ; NULLIFY(FAST_COEF_temp%ch4)
    ENDIF

    CALL H5GCLOSE_F( G_ID_SUB, ERR )
    THROWM(ERR.NE.0,"CANNOT CLOSE GROUP SOLAR")

    coef%SOLARCOEF = .TRUE.

  ENDIF


  coef%NLTECOEF = .FALSE.
  NULLIFY( coef%NLTE_COEF)

!   CALL H5LEXISTS_F( G_ID, "NLTE_COEF", LEXT, ERR )
!   THROWM(ERR.NE.0,"CANNOT TEST NLTE_COEF ")
!   IF( LEXT ) THEN
! 
!     CALL H5GOPEN_F( G_ID, "NLTE_COEF", G_ID_SUB, ERR )
!     THROWM(ERR.NE.0,"CANNOT OPEN GROUP NLTE_COEF")
! 
!     ALLOCATE (coef%NLTE_COEF, STAT = ERR)
!     THROWM( ERR .NE. 0, "ALLOCATION OF NLTE_COEF FAST COEFS")
! 
!     CALL RTTOV_HDF_RTTOV_NLTE_COEF_RH( coef%NLTE_COEF, G_ID_SUB, ERR )
!     THROWM(ERR.NE.0,"CANNOT READ COEF%NLTE_COEF")
! 
!     CALL H5GCLOSE_F( G_ID_SUB, ERR )
!     THROWM(ERR.NE.0,"CANNOT CLOSE GROUP NLTE_COEF")
! 
!     coef%NLTECOEF = .TRUE.
! 
!   ENDIF

  CALL H5GCLOSE_F( G_ID, ERR )
  THROWM(ERR.NE.0,"CANNOT CLOSE GROUP /COEF IN FILE "//TRIM(FILE_COEF))

  CALL H5FCLOSE_F( FILE_ID, ERR )
  THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(FILE_COEF))

#endif
  CATCH

CONTAINS

  SUBROUTINE move_hdf_fast_coef(coef, fast_coef, gas_in, gas_id)

    USE rttov_types, ONLY : rttov_fast_coef
    TYPE(rttov_coef), INTENT(INOUT) :: coef
    TYPE(rttov_fast_coef), INTENT(INOUT) :: fast_coef(:)
    REAL(jprb), INTENT(IN) :: gas_in(:,:,:)
    INTEGER(jpim), INTENT(IN) :: gas_id

    INTEGER(jpim) :: gas_pos, ncoef
    INTEGER(jpim) :: i

    gas_pos = coef%fmv_gas_pos(gas_id)
    ncoef = coef%fmv_coe(gas_pos)

    DO i = 1, coef%fmv_chn
      IF (ANY(gas_in(:, i, :) /= 0._jprb)) THEN
        ALLOCATE (fast_coef(i)%gasarray(gas_pos)%coef(ncoef, coef%nlayers))
        fast_coef(i)%gasarray(gas_pos)%coef = TRANSPOSE(gas_in(:, i , :))
        CALL set_pointers(fast_coef(i), gas_pos, gas_id)
      ELSE
        NULLIFY(fast_coef(i)%gasarray(gas_pos)%coef)
      ENDIF
    ENDDO

  END SUBROUTINE move_hdf_fast_coef

END SUBROUTINE rttov11_read_hdf5_coef
