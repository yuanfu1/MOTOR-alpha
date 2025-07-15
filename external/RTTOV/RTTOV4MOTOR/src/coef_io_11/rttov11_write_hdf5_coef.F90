! Description:
!> @file
!!   Write a v10/v11 HDF5 coefficient file from the v12 rttov_coef structure.
!!
!> @brief
!!   Write a v10/v11 HDF5 coefficient file from the v12 rttov_coef structure.
!!
!! @details
!!   Where necessary this subroutine converts the v12 coefficient structure to
!!   v10/v11 format. Suitable defaults are chosen for values which no longer
!!   exist in the v12 coefficient structure.
!!
!!   NB This subroutine should not be called for variable SO2 coefficients.
!!
!! @param[out]    err             status on exit
!! @param[in]     file_coef       full path of output file
!! @param[in]     coef            RTTOV optical depth coefficient structure
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
SUBROUTINE rttov11_write_hdf5_coef(err, file_coef, coef)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE parkind1, ONLY : jplm, jprb
#ifdef _RTTOV_HDF
  USE rttov_types, ONLY : rttov_fast_coef_hdf_io
  USE rttov_const, ONLY : &
         speedl,             &
         sensor_id_mw,       &
         sensor_id_po,       &
         gas_id_mixed,       &
         gas_id_watervapour, &
         gas_id_ozone,       &
         gas_id_wvcont,      &
         gas_id_co2,         &
         gas_id_n2o,         &
         gas_id_co,          &
         gas_id_ch4,         &
         tscale_v11, gscale_v11
  USE hdf5
  USE h5lt
  USE rttov_hdf_mod
  USE rttov_hdf_rttov_fast_coef_io
  USE rttov_hdf_rttov_nlte_coef_io
#endif
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim), INTENT(OUT) :: err
  CHARACTER(LEN=*),   INTENT(IN)  :: file_coef
  TYPE(rttov_coef)  , INTENT(IN)  :: coef
!INTF_END
#include "rttov_errorreport.interface"

#ifdef _RTTOV_HDF
  INTEGER(HID_T)       :: file_id, g_id, g_id_sub
  CHARACTER(LEN=lensh) :: sname
  INTEGER(KIND=jpim)   :: i
  INTEGER(KIND=jpim)   :: chanlist(coef%fmv_chn), gaz_units(coef%fmv_gas)
  REAL(KIND=jprb), ALLOCATABLE :: tlim(:), glim(:,:)
  TYPE(rttov_fast_coef_hdf_io) :: FAST_COEF_temp
#endif
!- End of header --------------------------------------------------------
  TRY

#ifndef _RTTOV_HDF
  err = errorstatus_fatal
  THROWM(err.NE.0, "")
#else

  IF (coef%fmv_gas > 8) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0,'Too many gases: cannot write variable SO2 files in v10/v11 format')
  ENDIF

  CALL RTTOV_HDF_OPEN_TRUNC( ERR, file_coef, FILE_ID)
  THROWM(ERR.NE.0,"CANNOT OPEN FILE "//TRIM(file_coef))

  CALL MKPAR( FILE_ID, "/COEF", G_ID, ERR )
  THROWM(ERR.NE.0,"CANNOT CREATE GROUP /COEF IN FILE "//TRIM(file_coef))

  chanlist = (/ (i, i = 1, coef%fmv_chn) /)

  sname='ID_PLATFORM'
  CALL write_array_hdf(g_id,sname,&
      'Platform ID, see documentation',&
      err,i0=coef%ID_PLATFORM  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='ID_SAT'
  CALL write_array_hdf(g_id,sname,&
      'Satellite ID, see documentation',&
      err,i0=coef%ID_SAT  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='ID_INST'
  CALL write_array_hdf(g_id,sname,&
      'Instrument ID, see documentation',&
      err,i0=coef%ID_INST  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='ID_SENSOR'
  CALL write_array_hdf(g_id,sname,&
      'Sensor ID,  1 = Infrared, 2 = Microwave, 3 = High resolution, 4 = Polarimetric',&
      err,i0=coef%ID_SENSOR  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='ID_COMP_LVL'
  CALL write_array_hdf(g_id,sname,&
      'RTTOV coefficient file version number',&
      err,i0=10_jpim  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='ID_COMP_PC'
  CALL write_array_hdf(g_id,sname,&
      'Principal component coefficient file version number',&
      err,i0=coef%ID_COMP_PC  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='ID_CREATION_DATE'
  CALL write_array_hdf(g_id,sname,&
      'Creation Date YYYY MM DD',&
      err,i1=coef%ID_CREATION_DATE  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='ID_CREATION'
  CALL write_array_hdf(g_id,sname,&
      'Creation comment',&
      err,c0=coef%ID_CREATION  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='ID_COMMON_NAME'
  CALL write_array_hdf(g_id,sname,&
      'Usual name of the satellite',&
      err,c0=coef%ID_COMMON_NAME  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='LINE_BY_LINE'
  CALL write_array_hdf(g_id,sname,&
      'Line-by-Line information',&
      err,c1=coef%LINE_BY_LINE  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='FMV_MODEL_DEF'
  CALL write_array_hdf(g_id,sname,&
      'Fast model predictors definition (RTTOV7 RTTOV8 RTTOV9)',&
      err,c0=coef%FMV_MODEL_DEF  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='FMV_MODEL_VER'
  CALL write_array_hdf(g_id,sname,&
      'Fast model version compatibility level',&
      err,i0=coef%FMV_MODEL_VER  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))


  sname='FMV_CHN'
  CALL write_array_hdf(g_id,sname,&
      'Number of channels read into coef structure',&
      err,i0=coef%FMV_CHN  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='FMV_GAS'
  CALL write_array_hdf(g_id,sname,&
      'Number of gases in file',&
      err,i0=coef%FMV_GAS  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='FMV_GAS_ID'
  IF (ASSOCIATED(coef%FMV_GAS_ID)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Gas ID numebr in gas_id list',&
        err,i1=coef%FMV_GAS_ID  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='FMV_GAS_POS'
  IF (ASSOCIATED(coef%FMV_GAS_POS)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Respective position of each gas of gas_id list',&
        err,i1=coef%FMV_GAS_POS(1:8)  )  ! ngases_max = 8 in RTTOV v11 but is 9 in RTTOV v12
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='FMV_VAR'
  IF (ASSOCIATED(coef%FMV_VAR)) THEN
    CALL write_array_hdf(g_id,sname,&
        'number of variables/predictors by gas',&
        err,i1=coef%FMV_VAR  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='FMV_COE'
  IF (ASSOCIATED(coef%FMV_COE)) THEN
    CALL write_array_hdf(g_id,sname,&
        'number of coefficients by gas',&
        err,i1=coef%FMV_COE  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='FMV_LVL'
  IF (ASSOCIATED(coef%FMV_LVL)) THEN
    CALL write_array_hdf(g_id,sname,&
        'number of levels(pres/absorber) by gas',&
        err,i1=coef%FMV_LVL  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='INCZEEMAN'
  CALL write_array_hdf(g_id,sname,&
      'Flag to include Zeeman effect for this sensor',&
      err,l0=coef%INCZEEMAN  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))


  sname='FC_SPEEDL'
  CALL write_array_hdf(g_id,sname,&
      'Speed of light',&
      err,r0=speedl , units = 'cm/s' )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='FC_PLANCK_C1'
  CALL write_array_hdf(g_id,sname,&
      'First radiation constant',&
      err,r0=coef%FC_PLANCK_C1 , units = 'mW/(m2*sr*cm-4)' )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='FC_PLANCK_C2'
  CALL write_array_hdf(g_id,sname,&
      'Second radiation constant',&
      err,r0=coef%FC_PLANCK_C2 , units = 'cm*K' )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='FC_SAT_HEIGHT'
  CALL write_array_hdf(g_id,sname,&
      'Satellite nominal altitude',&
      err,r0=coef%FC_SAT_HEIGHT , units = 'km' )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='GAZ_UNITS'
  gaz_units = 2  ! v11 value for ppmv dry
  CALL write_array_hdf(g_id,sname,&
      'Unit of gaz concentration for each gaz. Default value is specific conc. (kg/kg) value inside RTTOV calculations (ppmv)',&
      err,i1=GAZ_UNITS  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='README_SRF'
  CALL write_array_hdf(g_id,sname,&
      'About filter functions used',&
      err,c1=coef%README_SRF  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='FF_ORI_CHN'
  IF (ASSOCIATED(coef%FF_ORI_CHN)) THEN
    CALL write_array_hdf(g_id,sname,&
        'original chan number',&
        err,i1=coef%FF_ORI_CHN  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='FF_VAL_CHN'
  IF (ASSOCIATED(coef%FF_VAL_CHN)) THEN
    CALL write_array_hdf(g_id,sname,&
        'validity of the channel (1=OK)',&
        err,i1=coef%FF_VAL_CHN  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='FF_CWN'
  IF (ASSOCIATED(coef%FF_CWN)) THEN
    CALL write_array_hdf(g_id,sname,&
        'cental wave number (cm-1)',&
        err,r1=coef%FF_CWN , units = 'cm^-1' )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='FF_BCO'
  IF (ASSOCIATED(coef%FF_BCO)) THEN
    CALL write_array_hdf(g_id,sname,&
        'band correction offset (K)',&
        err,r1=coef%FF_BCO , units = 'K' )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='FF_BCS'
  IF (ASSOCIATED(coef%FF_BCS)) THEN
    CALL write_array_hdf(g_id,sname,&
        'band correction slope (K/K)',&
        err,r1=coef%FF_BCS , units = 'K/K' )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='FF_GAM'
  IF (ASSOCIATED(coef%FF_GAM)) THEN
    CALL write_array_hdf(g_id,sname,&
        'gamma factor transm. correction',&
        err,r1=coef%FF_GAM  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='TT_CHN'
  IF (ASSOCIATED(coef%TT_VAL_CHN)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Transmittance threshold chan number',&
        err,i1=chanlist  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='TT_VAL_CHN'
  IF (ASSOCIATED(coef%TT_VAL_CHN)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Transmittance threshold channel validity',&
        err,i1=coef%TT_VAL_CHN  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='TT_CWN'
  IF (ASSOCIATED(coef%TT_VAL_CHN)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Trasmittance threshold cental wave number (cm-1)',&
        err,r1=coef%FF_CWN , units = 'cm^-1' )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='TT_A0'
  IF (ASSOCIATED(coef%TT_A0)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Transmittance threshold A0',&
        err,r1=coef%TT_A0  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='TT_A1'
  IF (ASSOCIATED(coef%TT_A1)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Transmittance threshold A1',&
        err,r1=coef%TT_A1  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='PW_CHN'
  IF (ASSOCIATED(coef%PW_VAL_CHN)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Planck_weighted chan number',&
        err,i1=chanlist  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='PW_VAL_CHN'
  IF (ASSOCIATED(coef%PW_VAL_CHN)) THEN
    CALL write_array_hdf(g_id,sname,&
        '0 => non-PW thermal coefs, 1 => PW thermal coefs',&
        err,i1=coef%PW_VAL_CHN  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='SS_CHN'
  IF (ASSOCIATED(coef%SS_VAL_CHN)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Solar spectrum chan number',&
        err,i1=chanlist  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='SS_VAL_CHN'
  IF (ASSOCIATED(coef%SS_VAL_CHN)) THEN
    CALL write_array_hdf(g_id,sname,&
        '0 => thermal, 1 => thermal+solar , 2=> solar',&
        err,i1=coef%SS_VAL_CHN  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='SS_CWN'
  IF (ASSOCIATED(coef%SS_VAL_CHN)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Solar spectrum cental wave number (cm-1)',&
        err,r1=coef%FF_CWN , units = 'cm^-1' )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='SS_SOLAR_SPECTRUM'
  IF (ASSOCIATED(coef%SS_SOLAR_SPECTRUM)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Solar spectrum value',&
        err,r1=coef%SS_SOLAR_SPECTRUM  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='WOC_CHN'
  IF (ASSOCIATED(coef%WOC_WAOPC_OW)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Water optical constants chan number',&
        err,i1=chanlist  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='WOC_CWN'
  IF (ASSOCIATED(coef%WOC_WAOPC_OW)) THEN
    CALL write_array_hdf(g_id,sname,&
        'cental wave number (cm-1)',&
        err,r1=coef%FF_CWN , units = 'cm^-1' )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='WOC_WAOPC_OW'
  IF (ASSOCIATED(coef%WOC_WAOPC_OW)) THEN
    CALL write_array_hdf_cmplx(g_id,sname,&
        'Water optical constants ocean water',&
        err,k1=coef%WOC_WAOPC_OW  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='WOC_WAOPC_FW'
  IF (ASSOCIATED(coef%WOC_WAOPC_FW)) THEN
    CALL write_array_hdf_cmplx(g_id,sname,&
        'Water optical constants fresh water',&
        err,k1=coef%WOC_WAOPC_FW  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='WS_NOMEGA'
  CALL write_array_hdf(g_id,sname,&
      'Size of WAVE_SPECTRUM arrays',&
      err,i0=coef%WS_NOMEGA  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='WS_NPOINT'
  IF (ASSOCIATED(coef%WS_NPOINT)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Pulsation for waves at the surface of water',&
        err,r1=coef%WS_NPOINT , units = 'rad s-1' )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='WS_K_OMEGA'
  IF (ASSOCIATED(coef%WS_K_OMEGA)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Wavenumber (for water waves)',&
        err,r1=coef%WS_K_OMEGA , units = 'm-1' )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  IF (coef%ID_SENSOR == sensor_id_mw .OR. coef%ID_SENSOR == sensor_id_po) THEN
    sname='FASTEM_VER'
    CALL write_array_hdf(g_id,sname,&
        'Fastem version number',&
        err,i0=5_jpim  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

    sname='FASTEM_POLAR'
    IF (ASSOCIATED(coef%FASTEM_POLAR)) THEN
      CALL write_array_hdf(g_id,sname,&
          'polarisation of each channel',&
          err,i1=coef%FASTEM_POLAR , &
          units = '0=0.5 V+H  1=90-incident angle  2=incident angle  3=V  4=H  5=3rd Stokes comp  6=4th Stokes comp')
      THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
    ENDIF
  ENDIF

  sname='SSIREM_VER'
  CALL write_array_hdf(g_id,sname,&
      'ISEM version number',&
      err,i0=coef%SSIREM_VER  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='SSIREM_CHN'
  IF (ASSOCIATED(coef%SSIREM_A0)) THEN
    CALL write_array_hdf(g_id,sname,&
        'ISEM chan number',&
        err,i1=chanlist  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='SSIREM_A0'
  IF (ASSOCIATED(coef%SSIREM_A0)) THEN
    CALL write_array_hdf(g_id,sname,&
        'ISEM coef a0',&
        err,r1=coef%SSIREM_A0  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='SSIREM_A1'
  IF (ASSOCIATED(coef%SSIREM_A1)) THEN
    CALL write_array_hdf(g_id,sname,&
        'ISEM coef a1',&
        err,r1=coef%SSIREM_A1  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='SSIREM_A2'
  IF (ASSOCIATED(coef%SSIREM_A2)) THEN
    CALL write_array_hdf(g_id,sname,&
        'ISEM coef a2',&
        err,r1=coef%SSIREM_A2  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='SSIREM_XZN1'
  IF (ASSOCIATED(coef%SSIREM_XZN1)) THEN
    CALL write_array_hdf(g_id,sname,&
        'ISEM coef xzn1',&
        err,r1=coef%SSIREM_XZN1  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='SSIREM_XZN2'
  IF (ASSOCIATED(coef%SSIREM_XZN2)) THEN
    CALL write_array_hdf(g_id,sname,&
        'ISEM coef xzn2',&
        err,r1=coef%SSIREM_XZN2  )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF


  sname='REF_PRFL_P'
  IF (ASSOCIATED(coef%REF_PRFL_P)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Reference profile pressure',&
        err,r1=coef%REF_PRFL_P , units = 'hPa' )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='REF_PRFL_T'
  IF (ASSOCIATED(coef%REF_PRFL_T)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Reference profile temperature',&
        err,r2=coef%REF_PRFL_T , units = 'K' ,compress=.TRUE._jplm)
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='REF_PRFL_MR'
  IF (ASSOCIATED(coef%REF_PRFL_MR)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Reference profile gas concentration',&
        err,r2=coef%REF_PRFL_MR , units = 'ppmv' ,compress=.TRUE._jplm)
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='LIM_PRFL_P'
  IF (ASSOCIATED(coef%LIM_PRFL_P)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Envelope pressure',&
        err,r1=coef%LIM_PRFL_P , units = 'hPa' )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  ALLOCATE(tlim(coef%fmv_lvl(gas_id_mixed)), glim(coef%fmv_lvl(gas_id_mixed),coef%fmv_gas))

  tlim = coef%env_prfl_tmax * (1._jprb + tscale_v11)
  sname='LIM_PRFL_TMAX'
  CALL write_array_hdf(g_id,sname,&
      'Maximum temperature envelope',&
      err,r1=tlim , units = 'K' )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  tlim = coef%env_prfl_tmin * (1._jprb - tscale_v11)
  sname='LIM_PRFL_TMIN'
  CALL write_array_hdf(g_id,sname,&
      'Minimum temperature envelope',&
      err,r1=tlim , units = 'K' )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  glim(:,1)  = coef%env_prfl_gmax(:,1)
  glim(:,2:) = coef%env_prfl_gmax(:,2:) * (1._jprb + gscale_v11)
  sname='LIM_PRFL_GMAX'
  CALL write_array_hdf(g_id,sname,&
      'Maximum gas concentration envelope (levels, gases)',&
      err,r2=glim , units = 'ppmv' ,compress=.TRUE._jplm)
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  glim(:,1)  = coef%env_prfl_gmax(:,1)
  glim(:,2:) = coef%env_prfl_gmin(:,2:) * (1._jprb - gscale_v11)
  sname='LIM_PRFL_GMIN'
  CALL write_array_hdf(g_id,sname,&
      'Minimum gas concentration envelope (levels, gases)',&
      err,r2=glim , units = 'ppmv' ,compress=.TRUE._jplm)
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  DEALLOCATE(tlim, glim)

  sname='PMC_LENGTHCELL'
  CALL write_array_hdf(g_id,sname,&
      'Cell length',&
      err,r0=coef%PMC_LENGTHCELL , units = 'cm' )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='PMC_PNOMINAL'
  IF (ASSOCIATED(coef%PMC_PNOMINAL)) THEN
    CALL write_array_hdf(g_id,sname,&
        'Nominal cell pressure',&
        err,r1=coef%PMC_PNOMINAL , units = 'hPa' )
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF

  sname='PMC_TEMPCELL'
  CALL write_array_hdf(g_id,sname,&
      'Cell temperature',&
      err,r0=coef%PMC_TEMPCELL , units = 'K' )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='PMC_BETAPLUS1'
  CALL write_array_hdf(g_id,sname,&
      'CO2 band-average: self-HW/air-HW',&
      err,r0=coef%PMC_BETAPLUS1 , units = 'N/A' )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='PMC_NLAY'
  CALL write_array_hdf(g_id,sname,&
      'Number of layers used',&
      err,i0=coef%PMC_NLAY  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='PMC_NVAR'
  CALL write_array_hdf(g_id,sname,&
      'Number of variables used',&
      err,i0=coef%PMC_NVAR  )
  THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))

  sname='PMC_COEF'
  IF (ASSOCIATED(coef%PMC_COEF)) THEN
    CALL write_array_hdf(g_id,sname,&
        'PMC cell corrections',&
        err,r3=coef%PMC_COEF , units = 'N/A' ,compress=.TRUE._jplm)
    THROWM(err.NE.0,"CANNOT WRITE "//TRIM(sname))
  ENDIF



  CALL H5LTSET_ATTRIBUTE_STRING_F(G_ID, '.', "Description",   &
  "This is a RTTOV coefficient structure" // &
  CHAR(0), ERR )
  THROWM(ERR.NE.0,"CANNOT WRITE ATTRIBUTE")

  CALL MKPAR( G_ID, "THERMAL", G_ID_SUB, ERR )
  THROWM(ERR.NE.0,"CANNOT CREATE GROUP THERMAL")

  IF (coef%ncmixed > 0_jpim ) &
    CALL move_fast_coef_hdf(coef, coef%thermal(:), FAST_COEF_temp%mixedgas, gas_id_mixed)
  IF (coef%ncwater > 0_jpim ) &
    CALL move_fast_coef_hdf(coef, coef%thermal(:), FAST_COEF_temp%watervapour, gas_id_watervapour)
  IF (coef%ncozone > 0_jpim ) &
    CALL move_fast_coef_hdf(coef, coef%thermal(:), FAST_COEF_temp%ozone, gas_id_ozone)
  IF (coef%ncwvcont > 0_jpim ) &
    CALL move_fast_coef_hdf(coef, coef%thermal(:), FAST_COEF_temp%wvcont, gas_id_wvcont)
  IF (coef%ncco2 > 0_jpim ) &
    CALL move_fast_coef_hdf(coef, coef%thermal(:), FAST_COEF_temp%co2, gas_id_co2)
  IF (coef%ncn2o > 0_jpim ) &
    CALL move_fast_coef_hdf(coef, coef%thermal(:), FAST_COEF_temp%n2o, gas_id_n2o)
  IF (coef%ncco > 0_jpim ) &
    CALL move_fast_coef_hdf(coef, coef%thermal(:), FAST_COEF_temp%co, gas_id_co)
  IF (coef%ncch4 > 0_jpim ) &
    CALL move_fast_coef_hdf(coef, coef%thermal(:), FAST_COEF_temp%ch4, gas_id_ch4)

  CALL RTTOV_HDF_RTTOV_FAST_COEF_WH( FAST_COEF_temp, G_ID_SUB, ERR, COMPRESS=.TRUE._jplm )
  THROWM(ERR.NE.0,"CANNOT WRITE THERMAL COEFFICIENTS")

  IF (coef%ncmixed > 0_jpim ) &
    DEALLOCATE(FAST_COEF_temp%mixedgas)    ; NULLIFY(FAST_COEF_temp%mixedgas)
  IF (coef%ncwater > 0_jpim ) &
    DEALLOCATE(FAST_COEF_temp%watervapour) ; NULLIFY(FAST_COEF_temp%watervapour)
  IF (coef%ncozone > 0_jpim ) &
    DEALLOCATE(FAST_COEF_temp%ozone)       ; NULLIFY(FAST_COEF_temp%ozone)
  IF (coef%ncwvcont > 0_jpim ) &
    DEALLOCATE(FAST_COEF_temp%wvcont)      ; NULLIFY(FAST_COEF_temp%wvcont)
  IF (coef%ncco2 > 0_jpim ) &
    DEALLOCATE(FAST_COEF_temp%co2)         ; NULLIFY(FAST_COEF_temp%co2)
  IF (coef%ncn2o > 0_jpim ) &
    DEALLOCATE(FAST_COEF_temp%n2o)         ; NULLIFY(FAST_COEF_temp%n2o)
  IF (coef%ncco > 0_jpim ) &
    DEALLOCATE(FAST_COEF_temp%co)          ; NULLIFY(FAST_COEF_temp%co)
  IF (coef%ncch4 > 0_jpim ) &
    DEALLOCATE(FAST_COEF_temp%ch4)         ; NULLIFY(FAST_COEF_temp%ch4)

  CALL H5GCLOSE_F( G_ID_SUB, ERR )
  THROWM(ERR.NE.0,"CANNOT CLOSE GROUP THERMAL")

  IF( ASSOCIATED( coef%SOLAR ) .AND. coef%SOLARCOEF ) THEN
    CALL MKPAR( G_ID, "SOLAR", G_ID_SUB, ERR )
    THROWM(ERR.NE.0,"CANNOT CREATE GROUP SOLAR")

    IF (coef%ncmixed > 0_jpim ) &
      CALL move_fast_coef_hdf(coef, coef%solar(:), FAST_COEF_temp%mixedgas, gas_id_mixed)
    IF (coef%ncwater > 0_jpim ) &
      CALL move_fast_coef_hdf(coef, coef%solar(:), FAST_COEF_temp%watervapour, gas_id_watervapour)
    IF (coef%ncozone > 0_jpim ) &
      CALL move_fast_coef_hdf(coef, coef%solar(:), FAST_COEF_temp%ozone, gas_id_ozone)
    IF (coef%ncwvcont > 0_jpim ) &
      CALL move_fast_coef_hdf(coef, coef%solar(:), FAST_COEF_temp%wvcont, gas_id_wvcont)
    IF (coef%ncco2 > 0_jpim ) &
      CALL move_fast_coef_hdf(coef, coef%solar(:), FAST_COEF_temp%co2, gas_id_co2)
    IF (coef%ncn2o > 0_jpim ) &
      CALL move_fast_coef_hdf(coef, coef%solar(:), FAST_COEF_temp%n2o, gas_id_n2o)
    IF (coef%ncco > 0_jpim ) &
      CALL move_fast_coef_hdf(coef, coef%solar(:), FAST_COEF_temp%co, gas_id_co)
    IF (coef%ncch4 > 0_jpim ) &
      CALL move_fast_coef_hdf(coef, coef%solar(:), FAST_COEF_temp%ch4, gas_id_ch4)

    CALL RTTOV_HDF_RTTOV_FAST_COEF_WH( FAST_COEF_temp, G_ID_SUB, ERR, COMPRESS=.TRUE._jplm )
    THROWM(ERR.NE.0,"CANNOT WRITE SOLAR COEFFICIENTS")

    IF (coef%ncmixed > 0_jpim ) &
      DEALLOCATE(FAST_COEF_temp%mixedgas)    ; NULLIFY(FAST_COEF_temp%mixedgas)
    IF (coef%ncwater > 0_jpim ) &
      DEALLOCATE(FAST_COEF_temp%watervapour) ; NULLIFY(FAST_COEF_temp%watervapour)
    IF (coef%ncozone > 0_jpim ) &
      DEALLOCATE(FAST_COEF_temp%ozone)       ; NULLIFY(FAST_COEF_temp%ozone)
    IF (coef%ncwvcont > 0_jpim ) &
      DEALLOCATE(FAST_COEF_temp%wvcont)      ; NULLIFY(FAST_COEF_temp%wvcont)
    IF (coef%ncco2 > 0_jpim ) &
      DEALLOCATE(FAST_COEF_temp%co2)         ; NULLIFY(FAST_COEF_temp%co2)
    IF (coef%ncn2o > 0_jpim ) &
      DEALLOCATE(FAST_COEF_temp%n2o)         ; NULLIFY(FAST_COEF_temp%n2o)
    IF (coef%ncco > 0_jpim ) &
      DEALLOCATE(FAST_COEF_temp%co)          ; NULLIFY(FAST_COEF_temp%co)
    IF (coef%ncch4 > 0_jpim ) &
      DEALLOCATE(FAST_COEF_temp%ch4)         ; NULLIFY(FAST_COEF_temp%ch4)

    CALL H5GCLOSE_F( G_ID_SUB, ERR )
    THROWM(ERR.NE.0,"CANNOT CLOSE GROUP SOLAR")
  ENDIF

!   IF( coef%NLTECOEF ) THEN
!     CALL MKPAR( G_ID, "NLTE_COEF", G_ID_SUB, ERR )
!     THROWM(ERR.NE.0,"CANNOT CREATE GROUP NLTE_COEF")
! 
!     CALL RTTOV_HDF_RTTOV_NLTE_COEF_WH( coef%NLTE_COEF, G_ID_SUB, ERR, COMPRESS=.TRUE._jplm )
!     THROWM(ERR.NE.0,"CANNOT WRITE COEF%NLTE_COEF")
! 
!     CALL H5GCLOSE_F( G_ID_SUB, ERR )
!     THROWM(ERR.NE.0,"CANNOT CLOSE GROUP NLTE_COEF")
!   ENDIF


  CALL H5GCLOSE_F( G_ID, ERR )
  THROWM(ERR.NE.0,"CANNOT CLOSE GROUP /COEF IN FILE "//TRIM(file_coef))

  CALL H5FCLOSE_F( FILE_ID, ERR )
  THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(file_coef))

#endif
  CATCH

CONTAINS

  SUBROUTINE move_fast_coef_hdf(coef, fast_coef, gas_in, gas_id)

    USE rttov_types, ONLY : rttov_fast_coef
    TYPE(rttov_coef), INTENT(IN) :: coef
    TYPE(rttov_fast_coef), INTENT(IN) :: fast_coef(:)
    REAL(jprb), INTENT(INOUT), POINTER :: gas_in(:,:,:)
    INTEGER(jpim), INTENT(IN) :: gas_id

    INTEGER(jpim) :: i, gas_pos, ncoef

    gas_pos = coef%fmv_gas_pos(gas_id)
    ncoef = coef%fmv_coe(gas_pos)

    ALLOCATE(gas_in(coef%nlayers,coef%fmv_chn,ncoef))

    DO i = 1, coef%fmv_chn
      IF (ASSOCIATED(fast_coef(i)%gasarray(gas_pos)%coef)) THEN
        gas_in(:,i,:) = TRANSPOSE(fast_coef(i)%gasarray(gas_pos)%coef)
      ELSE
        gas_in(:,i,:) = 0._jprb
      ENDIF
    ENDDO

  END SUBROUTINE move_fast_coef_hdf

END SUBROUTINE rttov11_write_hdf5_coef
