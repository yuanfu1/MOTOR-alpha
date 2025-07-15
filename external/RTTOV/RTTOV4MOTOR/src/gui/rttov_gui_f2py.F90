! Description:
!> @file
!!   Subroutines for the GUI: this module defines the interface to run RTTOV
!!   from Python for the GUI.
!
!> @brief
!!   Subroutines for the GUI: this module defines the interface to run RTTOV
!!   from Python for the GUI.
!!
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
!    Copyright 2017, EUMETSAT, All Rights Reserved.
!
SUBROUTINE RTTOV_GUI_F2PY
! This subroutine is only here 
! to ensure the Makefile.PL will create dependencies
!INTF_END
USE PARKIND1
USE RTTOV_TYPES
USE RTTOV_CONST
USE RTTOV_HDF_MOD
USE RTTOV_GUI_HANDLE
USE RTTOV_GUI_CONTEXT

END SUBROUTINE RTTOV_GUI_F2PY


SUBROUTINE RTTOV_GUI_LOAD( CHANNELS, FILE_COEF, FILE_SCAER, FILE_SCCLD, FILE_MFASIS_CLD  ,FILE_PCCOEF, NCHANNELS, ERR )
!
#include "throw.h"

USE PARKIND1
USE RTTOV_GUI_HANDLE
USE RTTOV_HDF_MOD
USE RTTOV_CONST, ONLY: sensor_id_hi     

  IMPLICIT NONE

  INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
  INTEGER(KIND=JPIM), INTENT(IN)  :: CHANNELS(:)
  CHARACTER(LEN=*),   INTENT(IN)  :: FILE_COEF
  CHARACTER(LEN=*),   INTENT(IN)  :: FILE_SCAER
  CHARACTER(LEN=*),   INTENT(IN)  :: FILE_SCCLD
  CHARACTER(LEN=*),   INTENT(IN)  :: FILE_MFASIS_CLD
  CHARACTER(LEN=*),   INTENT(IN)  :: FILE_PCCOEF

  INTEGER(KIND=JPIM), INTENT(OUT) :: NCHANNELS
!INTF_END

#include "rttov_print_info.interface"

#include "throw.h"
  LOGICAL(KIND=jplm) :: verbose

TRY
  
  ERR = 0
  NCHANNELS = 0
  
  CALL OPEN_HDF( .FALSE., ERR ) ! 32 bits
  THROW(ERR.NE.0)

  IF( ASSOCIATED( THE_RTH ) ) THEN
    INFO("Drop coefficients: ")
    INFO(THE_RTH%FILE_COEF)
    IF( TRIM( THE_RTH%FILE_SCAER ) .NE. "" ) &
      INFO(THE_RTH%FILE_SCAER)
    IF( TRIM( THE_RTH%FILE_SCCLD ) .NE. "" ) &
      INFO(THE_RTH%FILE_SCCLD)
    IF( TRIM( THE_RTH%FILE_PCCOEF ) .NE. "" ) &
      INFO(THE_RTH%FILE_PCCOEF)

    CALL RTTOV_GUI_HANDLEDROP( THE_RTH, ERR )
    THROW(ERR.NE.0)
  ENDIF
  
  CALL RTTOV_GUI_HANDLELOAD(THE_RTH, CHANNELS, FILE_COEF, FILE_SCAER, FILE_SCCLD, FILE_MFASIS_CLD  ,FILE_PCCOEF, ERR )
  THROW(ERR.NE.0)
  
  INFO("Load coefficients: ")
  INFO(THE_RTH%FILE_COEF)
  IF( TRIM( THE_RTH%FILE_SCAER ) .NE. "" ) &
      INFO(THE_RTH%FILE_SCAER)
  IF( TRIM( THE_RTH%FILE_SCCLD ) .NE. "" ) &
      INFO(THE_RTH%FILE_SCCLD)
  IF( TRIM( THE_RTH%FILE_MFASIS_CLD ) .NE. "" ) &
      INFO(THE_RTH%FILE_MFASIS_CLD)
  IF( TRIM( THE_RTH%FILE_PCCOEF ) .NE. "" ) &
      INFO(THE_RTH%FILE_PCCOEF)

  NCHANNELS = THE_RTH%COEFS % COEF % FMV_CHN
  
  verbose= ( THE_RTH%COEFS % COEF %id_sensor .NE. sensor_id_hi )
  CALL RTTOV_PRINT_INFO (THE_RTH%COEFS, verbose=verbose)

  CALL CLOSE_HDF( ERR )
  THROW(ERR.NE.0)

CATCH
  
  CALL CLOSE_HDF( ERR )
  ERR = 1

  THE_RTH => NULL()  

END SUBROUTINE
SUBROUTINE RTTOV_GUI_DROP( ERR )
!
#include "throw.h"

USE PARKIND1
!INTF_OFF
USE RTTOV_GUI_HANDLE
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
!INTF_END
!
TRY

  ERR = 0

  IF( .NOT. ASSOCIATED( THE_RTH ) ) THEN
    INFO("No coefficients loaded")
    RETURN
  ENDIF
  !
  ! Note that handledrop drop coefficients and atlases
  !
  INFO("Drop coefficients: ")
  INFO(THE_RTH%FILE_COEF)
  IF( TRIM( THE_RTH%FILE_SCAER ) .NE. "" ) &
    INFO(THE_RTH%FILE_SCAER)
  IF( TRIM( THE_RTH%FILE_SCCLD ) .NE. "" ) &
    INFO(THE_RTH%FILE_SCCLD)
  IF( TRIM( THE_RTH%FILE_MFASIS_CLD ) .NE. "" ) &
      INFO(THE_RTH%FILE_MFASIS_CLD)
  IF( TRIM( THE_RTH%FILE_PCCOEF ) .NE. "" ) &
    INFO(THE_RTH%FILE_PCCOEF)
 
  CALL RTTOV_GUI_HANDLEDROP( THE_RTH, ERR )
  THROWM(ERR.NE.0,"An error occured while dropping coefficients")

  THE_RTH => NULL()
  
CATCH
  
THE_RTH => NULL()

END SUBROUTINE

SUBROUTINE RTTOV_GUI_GET_EMISBRDF( PATH_ATLAS, FILENAME_PROF, FILENAME_SURFACE, &
   & EMIS_ATLAS_ID, BRDF_ATLAS_ID, ERR )
!
#include "throw.h"

USE PARKIND1, ONLY : JPIM

!INTF_OFF
  USE RTTOV_TYPES
  USE RTTOV_CONST         
  USE RTTOV_HDF_MOD
  USE MOD_RTTOV_EMIS_ATLAS
  USE RTTOV_GUI_HANDLE
  USE RTTOV_GUI_CONTEXT
!INTF_ON
  
  IMPLICIT NONE
  CHARACTER(LEN=*),   INTENT(IN)  :: PATH_ATLAS, FILENAME_PROF, FILENAME_SURFACE
  INTEGER(KIND=JPIM), INTENT(IN)  :: EMIS_ATLAS_ID, BRDF_ATLAS_ID
  INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
!INTF_END

#include "rttov_hdf_load.interface"
#include "rttov_hdf_save.interface"

#include "rttov_setup_emis_atlas.interface"
#include "rttov_get_emis.interface"
#include "rttov_setup_brdf_atlas.interface"
#include "rttov_get_brdf.interface"

#include "rttov_deallocate_brdf_atlas.interface"
#include "rttov_deallocate_emis_atlas.interface"


!
  TYPE(RTTOV_CHANPROF), ALLOCATABLE :: CHANPROF(:)      ! PROFILE/CHANNEL INDICES (NCHANPROF)
  LOGICAL(KIND=JPLM) ,     POINTER  :: CALCEMIS(:)    => NULL()  ! SWITCHES FOR EMIS CALCS
  TYPE(RTTOV_EMISSIVITY),  POINTER  :: EMISSIVITY(:)  => NULL()  ! SURFACE EMIS
  LOGICAL(KIND=JPLM) ,     POINTER  :: CALCREFL(:)    => NULL()  ! SWITCHES FOR REFL CALCS
  TYPE(RTTOV_REFLECTANCE), POINTER  :: REFLECTANCE(:) => NULL()  ! SURFACE REFL
  TYPE(RTTOV_OPTIONS)               :: OPTS             ! RTTOV OPTIONS

  TYPE(RTTOVGUICONTEXT_TYPE)        :: CONTEXT
!
  INTEGER(KIND=JPIM) :: NCHANNELS
  INTEGER(KIND=JPIM) :: I

  INTEGER(KIND=JPIM) :: MONTH, ATLAS_TYPE
  LOGICAL(KIND=JPLM) :: DOTHERMAL
  LOGICAL(KIND=JPLM) :: DOSOLAR
  LOGICAL(KIND=JPLM) :: DOPURESOLAR

TRY

  IF( .NOT. ASSOCIATED( THE_RTH ) ) THEN
    ERR = 1
    THROWM(ERR.NE.0,"No coefficients file loaded")
  ENDIF
  
  CALL OPEN_HDF( .FALSE., ERR ) ! 32 bits
  THROW(ERR.NE.0)


! Read profile and options from HDF
  CALL RTTOVCONTEXTINITFROMPROFILE( CONTEXT, FILENAME_PROF, THE_RTH, ERR )
  THROW(ERR.NE.0)

  NCHANNELS = CONTEXT%RTH%COEFS%COEF%FMV_CHN
  MONTH     = CONTEXT%PROFILES(1)%DATE(2)
  
  ALLOCATE( CHANPROF  ( NCHANNELS ) ,STAT= ERR)
  THROW(ERR.NE.0)
! This interface is only valid for all channels, one profile
  CHANPROF(:)%CHAN = (/ (I, I=1, NCHANNELS) /)
  CHANPROF(:)%PROF = 1_JPIM

  DOTHERMAL = .TRUE.
  DOSOLAR   = .FALSE.
  DOPURESOLAR   = .FALSE.
    ! ss_val_chn = 0 => thermal-only
    ! ss_val_chn = 1 => thermal + solar
    ! ss_val_chn = 2 => solar-only
  IF (ASSOCIATED(CONTEXT%RTH%COEFS%COEF%SS_VAL_CHN)) THEN
    DOTHERMAL = ANY(CONTEXT%RTH%COEFS%COEF%SS_VAL_CHN < 2_JPIM)
    DOSOLAR   = ANY(CONTEXT%RTH%COEFS%COEF%SS_VAL_CHN > 0_JPIM)
    DOPURESOLAR   = ANY(CONTEXT%RTH%COEFS%COEF%SS_VAL_CHN > 1_JPIM)
  ENDIF
  
  IF( (CONTEXT%RTH%EMS_ATLAS_MONTH .NE. MONTH .OR.  &
    &  CONTEXT%RTH%EMIS_ATLAS_ID .NE. EMIS_ATLAS_ID ) .AND. &
    & DOTHERMAL ) THEN
    
    INFO("LOAD EMIS ATLAS")
    IF( CONTEXT%RTH%EMS_ATLAS_MONTH .GE. 1) THEN
      INFO("DROP PREVIOUS EMIS ATLAS")
      CALL RTTOV_DEALLOCATE_EMIS_ATLAS(CONTEXT%RTH%EMIS_ATLAS)
      CONTEXT%RTH%EMS_ATLAS_MONTH = -1_JPIM
      CONTEXT%RTH%EMIS_ATLAS_ID     = -1_JPIM
    ENDIF

    IF (CONTEXT%RTH%COEFS%COEF%ID_SENSOR == SENSOR_ID_MW .OR. &
        CONTEXT%RTH%COEFS%COEF%ID_SENSOR == SENSOR_ID_PO) THEN
      ATLAS_TYPE = ATLAS_TYPE_MW ! MW atlas
    ELSE
      ATLAS_TYPE = ATLAS_TYPE_IR ! IR atlas
    ENDIF
    CALL RTTOV_SETUP_EMIS_ATLAS(   &
             &  ERR,               &! OUT
             &  OPTS,              &! IN
             &  MONTH,             &! IN
             &  ATLAS_TYPE,        &! IN
             &  CONTEXT%RTH%EMIS_ATLAS, &
             &  ATLAS_ID = EMIS_ATLAS_ID, &
             &  PATH = TRIM(PATH_ATLAS)//"/emis_data", &
             &  COEFS = CONTEXT%RTH%COEFS, &
             &  IR_ATLAS_ANG_CORR = .TRUE.)
    THROW(ERR.NE.0)

    CONTEXT%RTH%EMS_ATLAS_MONTH = MONTH
    CONTEXT%RTH%EMIS_ATLAS_ID   = EMIS_ATLAS_ID
    INFO("LOAD EMIS ATLAS END")
  ENDIF

  IF( (CONTEXT%RTH%BRDF_ATLAS_MONTH .NE. MONTH .OR. &
    &  CONTEXT%RTH%BRDF_ATLAS_ID .NE. BRDF_ATLAS_ID ).AND. &
    & DOPURESOLAR ) THEN

    INFO("LOAD REFL ATLAS")
    IF( CONTEXT%RTH%BRDF_ATLAS_MONTH .GE. 1) THEN
      INFO("DROP PREVIOUS REFL ATLAS")
      CALL RTTOV_DEALLOCATE_BRDF_ATLAS(CONTEXT%RTH%BRDF_ATLAS)
      THROW(ERR.NE.0)
      CONTEXT%RTH%BRDF_ATLAS_MONTH = -1_jpim
      CONTEXT%RTH%BRDF_ATLAS_ID    = -1_JPIM
    ENDIF

    CALL RTTOV_SETUP_BRDF_ATLAS( &
             &  ERR,          &! OUT
             &  OPTS,         &! IN
             &  MONTH,        &! IN
             &  CONTEXT%RTH%BRDF_ATLAS, &
             &  ATLAS_ID = BRDF_ATLAS_ID, &
             &  PATH = TRIM(PATH_ATLAS)//"/brdf_data", &
             &  COEFS = CONTEXT%RTH%COEFS)
    THROW(ERR.NE.0)

    CONTEXT%RTH%BRDF_ATLAS_MONTH = MONTH
    CONTEXT%RTH%BRDF_ATLAS_ID    = BRDF_ATLAS_ID
    INFO("LOAD BRDF ATLAS END")
  ENDIF

  IF( DOTHERMAL ) THEN
    CALL RTTOV_HDF_LOAD( ERR, FILENAME_SURFACE, "/EMISSIVITY", EMISSIVITY=EMISSIVITY)
    THROWM(ERR.NE.0,"CANNOT LOAD EMISSIVITY FROM "//TRIM(FILENAME_SURFACE))
    CALL RTTOV_HDF_LOAD( ERR, FILENAME_SURFACE, "/EMISSIVITY", SNAME="CALCEMIS", PL1=CALCEMIS)
    THROWM(ERR.NE.0,"CANNOT LOAD CALCEMIS FROM "//TRIM(FILENAME_SURFACE))
    
    CALL RTTOV_GET_EMIS(    &
              & ERR,        &
              & OPTS,       &
              & CHANPROF,   &
              & CONTEXT%PROFILES,  &
              & CONTEXT%RTH%COEFS, &
              & CONTEXT%RTH%EMIS_ATLAS, &
              & EMISSIVITY(:)%EMIS_IN)
    THROW(ERR.NE.0)
    WHERE (EMISSIVITY(:)%EMIS_IN > 0.0_JPRB)
      CALCEMIS(:) = .FALSE.
    ELSEWHERE
      CALCEMIS(:) = .TRUE.
      EMISSIVITY(:)%EMIS_IN = 0.0_JPRB
    ENDWHERE
    EMISSIVITY(:)%EMIS_OUT = 0.0_JPRB

  ENDIF
  
  IF( DOSOLAR ) THEN
    ! in any case of solar calculation need to load reflectances, so arrays are allocated
    CALL RTTOV_HDF_LOAD( ERR, FILENAME_SURFACE, "/REFLECTANCE", REFLECTANCE=REFLECTANCE)
    THROWM(ERR.NE.0,"CANNOT LOAD REFLECTANCE FROM "//TRIM(FILENAME_SURFACE))
    CALL RTTOV_HDF_LOAD( ERR, FILENAME_SURFACE, "/REFLECTANCE", SNAME = "CALCREFL", PL1=CALCREFL)
    THROWM(ERR.NE.0,"CANNOT LOAD CALCREFL FROM "//TRIM(FILENAME_SURFACE))
    ! init values in case of channels in the range [2.5, 3.0]mu
    CALCREFL(:) = .TRUE.
    REFLECTANCE(:)%REFL_IN  = 0.0_JPRB
    REFLECTANCE(:)%REFL_OUT = 0.0_JPRB
  ENDIF

  IF( DOPURESOLAR ) THEN
  
    CALL RTTOV_GET_BRDF(    &
              & ERR,        &
              & OPTS,       &
              & CHANPROF,   &
              & CONTEXT%PROFILES,  &
              & CONTEXT%RTH%COEFS, &
              & CONTEXT%RTH%BRDF_ATLAS, &
              & REFLECTANCE(:)%REFL_IN)
    THROW(ERR.NE.0)
    WHERE (REFLECTANCE(:)%REFL_IN > 0.0_JPRB)
      CALCREFL(:) = .FALSE.
    ELSEWHERE
      CALCREFL(:) = .TRUE.
      REFLECTANCE(:)%REFL_IN = 0.0_JPRB
    ENDWHERE
    REFLECTANCE(:)%REFL_OUT = 0.0_JPRB

  ENDIF
  
  IF( DOTHERMAL ) THEN
      ! consider IR atlas is still valid up to 3300cm-1      
      ! Surface emissivity set to 
      !      1     in the range [2.5, 3.0]mu where we do not have any BRDF or EMIS atlas values
      !  1-refl*pi below 2.5mu where we have BRDF atlas values
    WHERE ( CONTEXT%RTH%COEFS%COEF%FF_CWN .GT. 3300.0_JPRB .AND. CONTEXT%RTH%COEFS%COEF%FF_CWN .LT. 4000.0_JPRB)
        EMISSIVITY%EMIS_IN = 1.0_JPRB
    ENDWHERE
    IF( DOPURESOLAR ) THEN
      WHERE ( (CONTEXT%RTH%COEFS%COEF%FF_CWN .GE. 4000.0_JPRB) .AND. .NOT. CALCREFL )
        EMISSIVITY%EMIS_IN = 1.0_JPRB - REFLECTANCE%REFL_IN*PI
      ENDWHERE
    ENDIF
    
    CALL RTTOV_HDF_SAVE( ERR, FILENAME_SURFACE, '/EMISSIVITY', CREATE=.FALSE., &
       &  EMISSIVITY = EMISSIVITY)
    THROW(ERR.NE.0)
    CALL RTTOV_HDF_SAVE( ERR, FILENAME_SURFACE, '/EMISSIVITY', CREATE=.FALSE., &
       &  L1 = CALCEMIS, SNAME='CALCEMIS')
    THROW(ERR.NE.0)

  ENDIF

  IF( DOSOLAR ) THEN
    ! Surface reflectance set to:
    !       0      in the range [2.5, 3.0]mu where we do not have any BRDF or EMIS atlas values
    !  (1-emis)/pi above 3.0mu where we have emissivity atlas values

    WHERE ( CONTEXT%RTH%COEFS%COEF%FF_CWN .GT. 3300.0_JPRB .AND. CONTEXT%RTH%COEFS%COEF%FF_CWN .LT. 4000.0_JPRB) 
        REFLECTANCE%REFL_IN = 0.0_JPRB
    ENDWHERE
 
    IF( DOTHERMAL ) THEN
      WHERE ( (CONTEXT%RTH%COEFS%COEF%FF_CWN .LT. 3300.0_JPRB) .AND. .NOT. CALCEMIS )
        REFLECTANCE%REFL_IN = (1.0_JPRB-EMISSIVITY%EMIS_IN)/PI
      ENDWHERE
    ENDIF
    
    CALL RTTOV_HDF_SAVE( ERR, FILENAME_SURFACE, '/REFLECTANCE', CREATE=.FALSE., &
        &  REFLECTANCE = REFLECTANCE)
    THROW(ERR.NE.0)
    CALL RTTOV_HDF_SAVE( ERR, FILENAME_SURFACE, '/REFLECTANCE', CREATE=.FALSE., &
        &  L1 = CALCREFL, SNAME='CALCREFL')
    THROW(ERR.NE.0)

  ENDIF
  
  CALL ADD_MISC_INFO( CONTEXT, FILENAME_SURFACE, ERR )
  THROW(ERR.NE.0)

  DEALLOCATE( CHANPROF, STAT= ERR)
  THROW(ERR.NE.0)
  
  IF( ASSOCIATED( CALCEMIS ) ) THEN
    DEALLOCATE( CALCEMIS, EMISSIVITY, STAT= ERR)
    THROW(ERR.NE.0)
  ENDIF
  
  IF( ASSOCIATED( CALCREFL ) ) THEN
    DEALLOCATE( CALCREFL, REFLECTANCE, STAT= ERR)
    THROW(ERR.NE.0)
  ENDIF

  ! remember profiles are allocated/deallocated by rttovcontextload/drop 
  CALL RTTOVCONTEXTDROP( CONTEXT, ERR )
  THROW(ERR.NE.0)


  CALL CLOSE_HDF( ERR )
  THROW(ERR.NE.0)

  CATCH
  
  CALL CLOSE_HDF( ERR )
  ERR = 1
  
  RETURN

END SUBROUTINE

SUBROUTINE RTTOV_GUI_RUN( FILENAME_PROF, FILENAME_SURFACE,  &
                     FILENAME_RADR, FILENAME_KMAT,          &
                     FILENAME_TRNS, MODE, NTHREADS, ERR )
!
#include "throw.h" 
USE parkind1, ONLY : jpim

!INTF_OFF
  USE RTTOV_TYPES, ONLY : &
      RTTOV_CHANPROF,     &
      RTTOV_TRANSMISSION, &
      RTTOV_RADIANCE,     &
      RTTOV_EMISSIVITY,   &
      RTTOV_REFLECTANCE,  &
      RTTOV_PROFILE,      &
      RTTOV_OPT_PARAM
  USE RTTOV_CONST
  USE RTTOV_HDF_MOD
  USE RTTOV_GUI_HANDLE
  USE RTTOV_GUI_CONTEXT
!INTF_ON
  
  IMPLICIT NONE
  CHARACTER(LEN=*),   INTENT(IN)  :: FILENAME_PROF, FILENAME_SURFACE,  &
                                     FILENAME_RADR, FILENAME_KMAT,   &
                                     FILENAME_TRNS
  CHARACTER(LEN=*),   INTENT(IN)  :: MODE
  INTEGER(KIND=JPIM), INTENT(IN)  :: NTHREADS
  INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
!INTF_END

!following file names should be in argument when user cld/aer opt param will be enabled
  CHARACTER(LEN=256) :: FILENAME_AER_OPT_PARAM, FILENAME_CLD_OPT_PARAM
 
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"

#include "rttov_hdf_load.interface"
#include "rttov_hdf_save.interface"

#include "rttov_parallel_direct.interface"
#include "rttov_parallel_k.interface"
#include "rttov_direct.interface"
#include "rttov_k.interface"

#include "rttov_print_opts.interface"
#include "rttov_user_options_checkinput.interface"

!
  TYPE(RTTOV_CHANPROF), ALLOCATABLE :: CHANPROF(:)      ! Profile/channel indices (nchanprof)
  TYPE(RTTOV_TRANSMISSION)          :: TRANSMISSION     ! transmittances and singlelayer optical depths (on User levels)
  TYPE(RTTOV_RADIANCE)              :: RADIANCE         ! radiances (mw/cm-1/ster/sq.m) and degK
  LOGICAL(KIND=JPLM),      POINTER  :: CALCEMIS(:)    => NULL()  ! switches for emis calcs
  TYPE(RTTOV_EMISSIVITY),  POINTER  :: EMISSIVITY(:)  => NULL()  ! surface emis
  LOGICAL(KIND=JPLM),      POINTER  :: CALCREFL(:)    => NULL()  ! switches for refl calcs
  TYPE(RTTOV_REFLECTANCE), POINTER  :: REFLECTANCE(:) => NULL()  ! surface refl

  TYPE(RTTOV_PROFILE),     POINTER  :: PROFILES_K(:)    => NULL()
  TYPE(RTTOV_TRANSMISSION)          :: TRANSMISSION_K
  TYPE(RTTOV_EMISSIVITY),  POINTER  :: EMISSIVITY_K(:)  => NULL()
  TYPE(RTTOV_REFLECTANCE), POINTER  :: REFLECTANCE_K(:) => NULL()
  TYPE(RTTOV_RADIANCE)              :: RADIANCE_K
  TYPE(RTTOV_OPT_PARAM)             :: AER_OPT_PARAM
  TYPE(RTTOV_OPT_PARAM)             :: CLD_OPT_PARAM
  TYPE(RTTOVGUICONTEXT_TYPE)        :: CONTEXT

  INTEGER(KIND=JPIM) :: NCHANNELS, NLEVELS, NLAYERS
  INTEGER(KIND=JPIM) :: ICHAN, I
  INTEGER(KIND=JPIM) :: ASW
  LOGICAL(KIND=JPLM) :: DOTHERMAL
  LOGICAL(KIND=JPLM) :: DOSOLAR

  TRY

  IF( .NOT. ASSOCIATED( THE_RTH ) ) THEN
    ERR = 1
    THROWM(ERR.NE.0,"No coefficients file loaded")
  ENDIF
  
  CALL OPEN_HDF( .FALSE., ERR ) ! 32 bits
  THROW(ERR.NE.0)


! Read profile and options from HDF
  CALL RTTOVCONTEXTINITFROMPROFILE( CONTEXT, FILENAME_PROF, THE_RTH, ERR )
  THROW(ERR.NE.0)
  
  CALL rttov_user_options_checkinput( &
       & err,                         &
       & CONTEXT%opts,            &
       & CONTEXT%RTH%coefs           )
  THROWM(ERR.NE.0, "ERROR COMPATIBILITY OPTIONS AND COEFS")

  CONTEXT%RTTOV_RUN_MODE = TRIM( MODE )
  WRITE ( 0, *) "contextmode =", CONTEXT%RTTOV_RUN_MODE

  NCHANNELS = CONTEXT%rth%COEFS%COEF%FMV_CHN
  NLAYERS   = CONTEXT%PROFILES(1)%NLAYERS
  NLEVELS   = CONTEXT%PROFILES(1)%NLEVELS

  DOTHERMAL = .TRUE.
  DOSOLAR   = .FALSE.
    ! ss_val_chn = 0 => thermal-only
    ! ss_val_chn = 1 => thermal + solar
    ! ss_val_chn = 2 => solar-only
  IF (ASSOCIATED(CONTEXT%RTH%COEFS%COEF%SS_VAL_CHN)) THEN
    DOTHERMAL = ANY(CONTEXT%RTH%COEFS%COEF%SS_VAL_CHN < 2_JPIM)
    DOSOLAR   = ANY(CONTEXT%RTH%COEFS%COEF%SS_VAL_CHN > 0_JPIM) .AND. &
                & CONTEXT%OPTS%RT_IR%ADDSOLAR
  ENDIF

  ! Check user aerosol/cloud optical parameters
  IF(CONTEXT%OPTS%RT_IR%USER_AER_OPT_PARAM) THEN
    ERR = 1_JPIM
    THROWM(ERR.NE.0,"USER AEROSOL OPTICAL PARAMETERS NOT IMPLEMENTED IN THIS VERSION OF GUI")
    CALL RTTOV_HDF_LOAD( ERR, FILENAME_AER_OPT_PARAM,"/USER_AEROSOL_OPT_PARAM", &
      & OPT_PARAM = AER_OPT_PARAM)
    THROWM(ERR.NE.0,"Cannot load user aerosol optical parameters from "//TRIM(FILENAME_AER_OPT_PARAM))
    ! check layers and channels
    IF ( SIZE(AER_OPT_PARAM%ABS,1_JPIM) .NE. NLAYERS ) THEN
      ERR = 1_JPIM
      THROWM(ERR.NE.0,"AEROSOL OPTICAL PARAMETERS DO NOT HAVE THE CORRECT NUMBER OF LAYERS")
    ENDIF
    IF ( SIZE(AER_OPT_PARAM%ABS,2_JPIM) .NE. NCHANNELS ) THEN
      ERR = 1_JPIM
      THROWM(ERR.NE.0,"AEROSOL OPTICAL PARAMETERS DO NOT HAVE THE CORRECT NUMBER OF CHANNELS")
    ENDIF
  ENDIF

  IF(CONTEXT%OPTS%RT_IR%USER_CLD_OPT_PARAM) THEN
    ERR = 1_JPIM
    THROWM(ERR.NE.0,"USER CLOUD OPTICAL PARAMETERS NOT IMPLEMENTED IN THIS VERSION OF GUI")
    CALL RTTOV_HDF_LOAD( ERR, FILENAME_CLD_OPT_PARAM,"/USER_CLOUD_OPT_PARAM", &
      & OPT_PARAM = CLD_OPT_PARAM)
    THROWM(ERR.NE.0,"Cannot load user cloud optical parameters from "//TRIM(FILENAME_CLD_OPT_PARAM))
    IF ( SIZE(CLD_OPT_PARAM%ABS,1_JPIM) .NE. NLAYERS ) THEN
      ERR = 1_JPIM
      THROWM(ERR.NE.0,"CLOUD OPTICAL PARAMETERS DO NOT HAVE THE CORRECT NUMBER OF LAYERS")
    ENDIF
    IF ( SIZE(CLD_OPT_PARAM%ABS,2_JPIM) .NE. NCHANNELS ) THEN
      ERR = 1_JPIM
      THROWM(ERR.NE.0,"CLOUD OPTICAL PARAMETERS DO NOT HAVE THE CORRECT NUMBER OF CHANNELS")
    ENDIF
  ENDIF

  ! allocate radiance results arrays with number of channels
  ASW = 1 ! ALLOCATE
  CALL RTTOV_ALLOC_RAD( ERR, NCHANNELS, RADIANCE, NLEVELS, ASW, INIT=.TRUE.)
  THROW(ERR.NE.0)
  
  ALLOCATE( CHANPROF ( NCHANNELS ), STAT= ERR)
  THROW(ERR.NE.0)
  
  IF( DOTHERMAL .AND. FILENAME_SURFACE.NE."" ) THEN
    CALL RTTOV_HDF_LOAD( ERR, FILENAME_SURFACE, "/EMISSIVITY", EMISSIVITY=emissivity)
    THROWM(ERR.NE.0,"Cannot load emissivity from "//TRIM(FILENAME_SURFACE))
    CALL RTTOV_HDF_LOAD( ERR, FILENAME_SURFACE, "/EMISSIVITY", sname="CALCEMIS", pl1=calcemis)
    THROWM(ERR.NE.0,"Cannot load calcemis from "//TRIM(FILENAME_SURFACE))
  ELSE
    ALLOCATE( CALCEMIS  ( NCHANNELS ) ,&
          &   EMISSIVITY( NCHANNELS ) ,&
              STAT= ERR)
    THROW(ERR.NE.0)
    CALCEMIS               = .TRUE.
    EMISSIVITY(:)%EMIS_IN  = 0._JPRB
    EMISSIVITY(:)%EMIS_OUT = 0._JPRB
  ENDIF
  
  If( DOSOLAR  .AND. FILENAME_SURFACE.NE."" ) Then
    CALL RTTOV_HDF_LOAD( ERR, FILENAME_SURFACE, "/REFLECTANCE", REFLECTANCE=reflectance)
    THROWM(ERR.NE.0,"Cannot load reflectance from "//TRIM(FILENAME_SURFACE))
    CALL RTTOV_HDF_LOAD( ERR, FILENAME_SURFACE, "/REFLECTANCE", sname = "CALCREFL", pl1=calcrefl)
    THROWM(ERR.NE.0,"Cannot load calcrefl from "//TRIM(FILENAME_SURFACE))
  Else
    ALLOCATE( CALCREFL   ( NCHANNELS ) ,&
          &   REFLECTANCE( NCHANNELS ) ,&
              STAT= ERR)
    THROW(ERR.NE.0)
    CALCREFL                = .TRUE.
    REFLECTANCE(:)%REFL_IN  = 0._JPRB
    REFLECTANCE(:)%REFL_OUT = 0._JPRB
  ENDIF


  ! allocate transmittance structure
  CALL RTTOV_ALLOC_TRANSMISSION( ERR, TRANSMISSION, NLEVELS, NCHANNELS,  &
      & ASW, INIT = .TRUE._JPLM)  
  THROW(ERR.NE.0)

  CHANPROF(:)%CHAN = (/ (I, I=1, NCHANNELS) /)
  CHANPROF(:)%PROF = 1_JPIM


  IF( CONTEXT%RTTOV_RUN_MODE(1:1) .EQ. 'K' ) THEN

    ALLOCATE( PROFILES_K( NCHANNELS ), EMISSIVITY_K( NCHANNELS ) , &
          &   REFLECTANCE_K( NCHANNELS ) ,&
              STAT= ERR)
    THROW(ERR.NE.0)

    CALL RTTOV_ALLOC_PROF(                                  &
      ERR, NCHANNELS, PROFILES_K, NLEVELS,  CONTEXT%OPTS  , &
      ASW, COEFS=CONTEXT%RTH%COEFS, INIT = .TRUE. )
    THROWM(ERR.NE.0,"ALLOCATION OF PROFILES")

    CALL RTTOV_ALLOC_TRANSMISSION(             &
      ERR, TRANSMISSION_K, NLEVELS, NCHANNELS, &
      ASW, INIT = .TRUE. )
    THROW(ERR.NE.0)

    DO ICHAN = 1, NCHANNELS
      PROFILES_K(ICHAN)%P = 0_JPRB
    ENDDO


    CALL RTTOV_ALLOC_RAD( ERR, NCHANNELS, RADIANCE_K, NLEVELS, ASW, INIT=.TRUE.)
    THROW(ERR.NE.0)

      EMISSIVITY_K(:)%EMIS_IN  = 0._JPRB
      EMISSIVITY_K(:)%EMIS_OUT = 0._JPRB
      REFLECTANCE_K(:)%REFL_IN  = 0._JPRB
      REFLECTANCE_K(:)%REFL_OUT = 0._JPRB

      IF( CONTEXT%OPTS%RT_ALL%SWITCHRAD ) THEN
        RADIANCE_K%BT    = 1._JPRB
      ELSE
        RADIANCE_K%TOTAL = 1._JPRB
      ENDIF
  ENDIF


  INFO("-- RTTOV CONFIGURATION --")
  CALL RTTOV_PRINT_OPTS (CONTEXT%OPTS, 0, "USER OPTIONS")
  INFO("-------------------------")


! Run RTTOV
  IF( CONTEXT%RTTOV_RUN_MODE .EQ. 'DIRECT' ) THEN
  INFO("Start Direct")
    IF( NTHREADS .GT. 1 ) THEN
    CALL RTTOV_PARALLEL_DIRECT( &
            & ERR,            &
            & CHANPROF,       &
            & CONTEXT%OPTS,   &
            & CONTEXT%PROFILES, &
            & CONTEXT%RTH%COEFS,&
            & TRANSMISSION,   &
            & RADIANCE,       &
            & CALCEMIS    = CALCEMIS,    &
            & EMISSIVITY  = EMISSIVITY,  &
            & CALCREFL    = CALCREFL,    &
            & REFLECTANCE = REFLECTANCE, &
            & NTHREADS=NTHREADS, &
            & DEBUG = .FALSE.  )

    ELSE
    CALL RTTOV_DIRECT(                  &
            & ERR,                      &
            & CHANPROF,                 &
            & CONTEXT%OPTS,             &
            & CONTEXT%PROFILES,         &
            & CONTEXT%RTH%COEFS,        &
            & TRANSMISSION,             &
            & RADIANCE,                 &
            & CALCEMIS    = CALCEMIS,   &
            & EMISSIVITY  = EMISSIVITY, &
            & CALCREFL    = CALCREFL,   &
            & REFLECTANCE = REFLECTANCE )
    ENDIF
    THROWM(ERR.NE.0,"An error occured while running RTTOV Direct")
    INFO("End Direct")

    CALL RTTOV_HDF_SAVE( ERR, FILENAME_RADR, '/RADIANCE', CREATE=.TRUE., &
       &  RADIANCE = RADIANCE )
    THROW(ERR.NE.0)
    CALL ADD_MISC_INFO( CONTEXT, FILENAME_RADR, ERR )
    THROW(ERR.NE.0)

    CALL RTTOV_HDF_SAVE( ERR, FILENAME_TRNS, '/TRANSMISSION', CREATE=.TRUE., &
       &  TRANSMISSION = TRANSMISSION )
    THROW(ERR.NE.0)
    CALL ADD_MISC_INFO( CONTEXT, FILENAME_TRNS, ERR )
    THROW(ERR.NE.0)

    IF ( DOTHERMAL  .AND. FILENAME_SURFACE.NE."" ) THEN
      CALL RTTOV_HDF_SAVE( ERR, FILENAME_SURFACE, '/EMISSIVITY', CREATE=.FALSE., &
         &  EMISSIVITY = EMISSIVITY)
      THROW(ERR.NE.0)
      CALL RTTOV_HDF_SAVE( ERR, FILENAME_SURFACE, '/EMISSIVITY', CREATE=.FALSE., &
         &  L1 = CALCEMIS, SNAME='CALCEMIS')
      THROW(ERR.NE.0)
    ENDIF
    
    IF ( DOSOLAR  .AND. FILENAME_SURFACE.NE."" ) THEN
      CALL RTTOV_HDF_SAVE( ERR, FILENAME_SURFACE, '/REFLECTANCE', CREATE=.FALSE., &
         &  REFLECTANCE = REFLECTANCE)
      THROW(ERR.NE.0)
      CALL RTTOV_HDF_SAVE( ERR, FILENAME_SURFACE, '/REFLECTANCE', CREATE=.FALSE., &
         &  L1 = CALCREFL, SNAME='CALCREFL')
      THROW(ERR.NE.0)
    ENDIF
        
    IF ( FILENAME_SURFACE.NE."" ) THEN
      CALL ADD_MISC_INFO( CONTEXT, FILENAME_SURFACE, ERR )
      THROW(ERR.NE.0)
    ENDIF

  ELSE IF( CONTEXT%RTTOV_RUN_MODE(1:1) .EQ. 'K' ) THEN
    INFO("Start K")

    IF( NTHREADS .GT. 1 ) THEN

      CALL RTTOV_PARALLEL_K( &
            & ERR,                           &
            & CHANPROF,                      &
            & CONTEXT%OPTS,                  &
            & CONTEXT%PROFILES,              &
            & PROFILES_K,                    &
            & CONTEXT%RTH%COEFS,             &
            & TRANSMISSION,                  &
            & TRANSMISSION_K,                &
            & RADIANCE,                      &
            & RADIANCE_K,                    &
            & CALCEMIS      = CALCEMIS,      &
            & EMISSIVITY    = EMISSIVITY,    &
            & EMISSIVITY_K  = EMISSIVITY_K,  &
            & CALCREFL      = CALCREFL,      &
            & REFLECTANCE   = REFLECTANCE,   &
            & REFLECTANCE_K = REFLECTANCE_K, &
            & NTHREADS      = NTHREADS)
    ELSE
    CALL RTTOV_K(  &
            & ERR,                           &
            & CHANPROF,                      &
            & CONTEXT%OPTS,                  &
            & CONTEXT%PROFILES,              &
            & PROFILES_K,                    &
            & CONTEXT%RTH%COEFS,             &
            & TRANSMISSION,                  &
            & TRANSMISSION_K,                &
            & RADIANCE,                      &
            & RADIANCE_K,                    &
            & CALCEMIS      = CALCEMIS,      &
            & EMISSIVITY    = EMISSIVITY,    &
            & EMISSIVITY_K  = EMISSIVITY_K,  &
            & CALCREFL      = CALCREFL,      &
            & REFLECTANCE   = REFLECTANCE,   &
            & REFLECTANCE_K = REFLECTANCE_K)
    ENDIF
    THROWM(ERR.NE.0,"An error occured while running RTTOV K")
    INFO("End K")

    
       CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/KMATRIX', CREATE=.TRUE., &
       &  KMATRIX = PROFILES_K, OPTS= CONTEXT%OPTS)
       THROW(ERR.NE.0)
                        
       CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/PROFILES', CREATE=.FALSE., &
       &  PROFILES = CONTEXT%PROFILES)
       THROW(ERR.NE.0)
                 
       CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/CHANPROF', CREATE=.FALSE., &
       &  CHANPROF = CHANPROF)
       THROW(ERR.NE.0)
                 
       CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/EMISSIVITY_K', CREATE=.FALSE., &
       &  EMISSIVITY = EMISSIVITY_K)
       THROW(ERR.NE.0)
       CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/EMISSIVITY_K', CREATE=.FALSE., &
       &  L1 = CALCEMIS, SNAME='CALCEMIS')
       THROW(ERR.NE.0)
         
       CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/OPTIONS', CREATE=.FALSE., &
       &  OPTIONS = CONTEXT%OPTS)
       THROW(ERR.NE.0)
       
       IF (CONTEXT%OPTS%RT_IR%ADDSOLAR) THEN
         CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/REFLECTANCE_K', CREATE=.FALSE., &
       &    REFLECTANCE = REFLECTANCE_K)
         THROW(ERR.NE.0)
         CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/REFLECTANCE_K', CREATE=.FALSE., &
       &  L1 = CALCREFL, SNAME='CALCREFL')
       THROW(ERR.NE.0)
       ENDIF

       CALL ADD_MISC_INFO( CONTEXT, FILENAME_KMAT, ERR )
       THROW(ERR.NE.0)

  ELSE
        ERR = 1
        THROWM(ERR.NE.0,"CONTEXT%RTTOV_RUN_MODE unknown "//CONTEXT%RTTOV_RUN_MODE)
  ENDIF


! Deallocate

  ASW = 0_JPIM
  CALL RTTOV_ALLOC_RAD( ERR, NCHANNELS, RADIANCE, NLEVELS, ASW)
  THROW(ERR.NE.0)
  
  DEALLOCATE( CHANPROF, STAT= ERR)
  THROW(ERR.NE.0)
  
  DEALLOCATE( CALCEMIS, EMISSIVITY, STAT= ERR)
  THROW(ERR.NE.0)
  
  DEALLOCATE( CALCREFL, REFLECTANCE, STAT= ERR)
  THROW(ERR.NE.0)

  CALL RTTOV_ALLOC_TRANSMISSION( ERR, TRANSMISSION, NLEVELS, NCHANNELS,  &
      & ASW)  
  THROW(ERR.NE.0)


  IF( CONTEXT%RTTOV_RUN_MODE(1:1) .EQ. 'K' ) THEN

    CALL RTTOV_ALLOC_PROF(                                  &
      ERR, NCHANNELS, PROFILES_K, NLEVELS,  CONTEXT%OPTS  , &
      ASW, COEFS=CONTEXT%RTH%COEFS)
    THROWM(ERR.NE.0,"DEALLOCATION OF PROFILES")

    CALL RTTOV_ALLOC_TRANSMISSION(             &
      ERR, TRANSMISSION_K, NLEVELS, NCHANNELS, &
      ASW )
    THROW(ERR.NE.0)

    CALL RTTOV_ALLOC_RAD( ERR, NCHANNELS, RADIANCE_K, NLEVELS, ASW)
    THROW(ERR.NE.0)
    DEALLOCATE( PROFILES_K, EMISSIVITY_K , &
          &   REFLECTANCE_K ,&
              STAT= ERR)
    THROW(ERR.NE.0)

  ENDIF

  ! remember profiles are allocated/deallocated by rttovcontextload/drop 
  CALL RTTOVCONTEXTDROP( CONTEXT, ERR )
  THROW(ERR.NE.0)

  CALL CLOSE_HDF( ERR )
  THROW(ERR.NE.0)

  CATCH
  
  CALL CLOSE_HDF( ERR )
  ERR = 1
  
  RETURN


END SUBROUTINE

SUBROUTINE RTTOV_GUI_PCRUN( FILENAME_PROF, FILENAME_SURFACE,  &
                     FILENAME_PC, FILENAME_KMAT,          &
                     MODE, NTHREADS, NPCSCORES, ERR )
!
#include "throw.h" 
USE parkind1, ONLY : jpim

!INTF_OFF
  USE RTTOV_TYPES, ONLY : &
      RTTOV_CHANPROF,     &
      RTTOV_TRANSMISSION, &
      RTTOV_RADIANCE,     &
      RTTOV_EMISSIVITY,   &
      RTTOV_REFLECTANCE,  &
      RTTOV_PROFILE,      &
!       RTTOV_OPT_PARAM,    &
      RTTOV_PCCOMP
  USE RTTOV_CONST
  USE RTTOV_HDF_MOD
  USE RTTOV_GUI_HANDLE
  USE RTTOV_GUI_CONTEXT
!INTF_ON
  
  IMPLICIT NONE
  CHARACTER(LEN=*),   INTENT(IN)  :: FILENAME_PROF, FILENAME_SURFACE,  &
                                     FILENAME_PC, FILENAME_KMAT
  CHARACTER(LEN=*),   INTENT(IN)  :: MODE
  INTEGER(KIND=JPIM), INTENT(IN)  :: NTHREADS
  INTEGER(KIND=JPIM), INTENT(IN)  :: NPCSCORES
  INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
!INTF_END

!following file names should be in argument when user cld/aer opt param will be enabled
!   CHARACTER(LEN=256) :: FILENAME_AER_OPT_PARAM, FILENAME_CLD_OPT_PARAM
 
#include "rttov_alloc_prof.interface"
#include "rttov_alloc_rad.interface"
#include "rttov_alloc_transmission.interface"
#include "rttov_alloc_pccomp.interface"

#include "rttov_hdf_save.interface"

#include "rttov_direct.interface"
#include "rttov_k.interface"

#include "rttov_print_opts.interface"
#include "rttov_user_options_checkinput.interface"

!
  TYPE(RTTOV_CHANPROF), ALLOCATABLE :: CHANPROF(:)      ! Profile/channel indices (nchanprof)
  TYPE(RTTOV_TRANSMISSION)          :: TRANSMISSION     ! transmittances and singlelayer optical depths (on User levels)
  TYPE(RTTOV_RADIANCE)              :: RADIANCE         ! radiances (mw/cm-1/ster/sq.m) and degK
  LOGICAL(KIND=JPLM),      POINTER  :: CALCEMIS(:)   => NULL()   ! switches for emis calcs
  TYPE(RTTOV_EMISSIVITY),  POINTER  :: EMISSIVITY(:) => NULL()   ! surface emis
!   LOGICAL(KIND=JPLM),      POINTER  :: CALCREFL(:)   => NULL()   ! switches for refl calcs
!   TYPE(RTTOV_REFLECTANCE), POINTER  :: REFLECTANCE(:) => NULL()  ! surface refl

  TYPE(RTTOV_PROFILE),     POINTER  :: PROFILES_K(:) => NULL()
  TYPE(RTTOV_TRANSMISSION)          :: TRANSMISSION_K
  TYPE(RTTOV_EMISSIVITY),  POINTER  :: EMISSIVITY_K(:)  => NULL()
  TYPE(RTTOV_REFLECTANCE), POINTER  :: REFLECTANCE_K(:) => NULL()
  TYPE(RTTOV_RADIANCE)              :: RADIANCE_K
!   TYPE(RTTOV_OPT_PARAM)             :: AER_OPT_PARAM
!   TYPE(RTTOV_OPT_PARAM)             :: CLD_OPT_PARAM
  TYPE(RTTOVGUICONTEXT_TYPE)        :: CONTEXT
  TYPE(RTTOV_PCCOMP)                :: PCCOMP
  
  TYPE(RTTOV_PROFILE),     POINTER  :: PROFILES_K_REC(:) => NULL()
  TYPE(RTTOV_PROFILE),     POINTER  :: PROFILES_K_PC(:)  => NULL()
  TYPE(RTTOV_PCCOMP)                :: PCCOMP_K

  INTEGER(KIND=JPIM) :: NCHANNELS, NLEVELS
  INTEGER(KIND=JPIM) :: ICHAN, J
  INTEGER(KIND=JPIM) :: ASW
  LOGICAL(KIND=JPLM) :: DOTHERMAL
  LOGICAL(KIND=JPLM) :: DOSOLAR

  INTEGER(KIND=jpim),      ALLOCATABLE :: channelsrec(:) ! Reconstructed radiance channel list
  INTEGER(KIND=jpim) :: nchannelsrec

  TRY

  IF( .NOT. ASSOCIATED( THE_RTH ) ) THEN
    ERR = 1
    THROWM(ERR.NE.0,"No coefficients file loaded")
  ENDIF
  
  CALL OPEN_HDF( .FALSE., ERR ) ! 32 bits
  THROW(ERR.NE.0)


! Read profile and options from HDF
  CALL RTTOVCONTEXTINITFROMPROFILE( CONTEXT, FILENAME_PROF, THE_RTH, ERR )
  THROW(ERR.NE.0)


  CALL rttov_user_options_checkinput( &
       & err,                         &
       & CONTEXT%opts,            &
       & CONTEXT%RTH%coefs           )
  THROWM(ERR.NE.0, "ERROR COMPATIBILITY OPTIONS AND COEFS")



  CONTEXT%RTTOV_RUN_MODE = TRIM( MODE )
  WRITE ( 0, *) "contextmode =", CONTEXT%RTTOV_RUN_MODE
  
  !NULLIFY(predictindex)
  !CALL rttov_get_pc_predictindex(ERR, CONTEXT%RTH%opts, predictindex, file_pccoef=pccoef_filename)
  !THROWM(ERR.NE.0,'rttov_get_pc_predictindex fatal error')
  
 
  nchannels = SIZE( CONTEXT%RTH%COEFS%coef_pccomp%pcreg ( &
                  & CONTEXT%opts%rt_ir%pc%ipcbnd,     &
                  & CONTEXT%opts%rt_ir%pc%ipcreg ) % predictindex)

  ALLOCATE( CHANPROF ( NCHANNELS ), STAT= ERR)
  THROW(ERR.NE.0)

  CHANPROF(:)%CHAN = CONTEXT%RTH%COEFS%coef_pccomp%pcreg ( &
                  & CONTEXT%opts%rt_ir%pc%ipcbnd,     &
                  & CONTEXT%opts%rt_ir%pc%ipcreg ) % predictindex(:)

  CHANPROF(:)%PROF = 1_JPIM

  NLEVELS   = CONTEXT%PROFILES(1)%NLEVELS


  IF (CONTEXT%opts%rt_ir%pc% addradrec) THEN
    nchannelsrec = CONTEXT%RTH%coefs % coef % fmv_chn
    ALLOCATE(channelsrec(nchannelsrec))
    channelsrec(:) = (/ (j, j = 1, nchannelsrec) /)
  ENDIF


  DOTHERMAL = .TRUE.
  DOSOLAR   = .FALSE.
    ! ss_val_chn = 0 => thermal-only
    ! ss_val_chn = 1 => thermal + solar
    ! ss_val_chn = 2 => solar-only
  IF (ASSOCIATED(CONTEXT%RTH%COEFS%COEF%SS_VAL_CHN)) THEN
    DOTHERMAL = ANY(CONTEXT%RTH%COEFS%COEF%SS_VAL_CHN < 2_JPIM)
    DOSOLAR   = ANY(CONTEXT%RTH%COEFS%COEF%SS_VAL_CHN > 0_JPIM) .AND. &
                & CONTEXT%OPTS%RT_IR%ADDSOLAR
  ENDIF

  ! Check user aerosol/cloud optical parameters
  IF(CONTEXT%OPTS%RT_IR%USER_AER_OPT_PARAM) THEN
    ERR = 1_JPIM
    THROWM(ERR.NE.0,"USER AEROSOL OPTICAL PARAMETERS NOT COMPATIBLE WITH PC CALCs")
  ENDIF

  IF(CONTEXT%OPTS%RT_IR%USER_CLD_OPT_PARAM) THEN
    ERR = 1_JPIM
    THROWM(ERR.NE.0,"USER CLOUD OPTICAL PARAMETERS NOT IMPLEMENTED IN THIS VERSION OF GUI")
  ENDIF

  ! allocate radiance results arrays with number of channels
  ASW = 1 ! ALLOCATE
  CALL RTTOV_ALLOC_RAD( ERR, NCHANNELS, RADIANCE, NLEVELS, ASW, INIT=.TRUE.)
  THROW(ERR.NE.0)
  
  
  ! SURFACE EMISSIVITY ALWAYS CALCULATED IN CASE OF PC
  ALLOCATE( CALCEMIS  ( NCHANNELS ) ,&
        &   EMISSIVITY( NCHANNELS ) ,&
            STAT= ERR)
  THROW(ERR.NE.0)
  CALCEMIS               = .TRUE.
  EMISSIVITY(:)%EMIS_IN  = 0._JPRB
  EMISSIVITY(:)%EMIS_OUT = 0._JPRB
  
  If( DOSOLAR ) Then  
    ERR = 1_JPIM
    THROWM(ERR.NE.0,"DOSOALR  NOT IMPLEMENTED IN THIS VERSION OF GUI")
  ENDIF



  ! allocate transmittance structure
  CALL RTTOV_ALLOC_TRANSMISSION( ERR, TRANSMISSION, NLEVELS, NCHANNELS,  &
      & ASW, INIT = .TRUE._JPLM)  
  THROW(ERR.NE.0)

  IF (CONTEXT%opts % rt_ir % pc % addradrec) THEN
    CALL rttov_alloc_pccomp(  &
        & ERR,        &
        & pccomp,             &
        & npcscores  ,  &
        & asw,                &
        & init = .TRUE._jplm, &
        & nchannels_rec = nchannelsrec)
  ELSE
    CALL rttov_alloc_pccomp( &
        & ERR,       &
        & pccomp,            &
        & npcscores , &
        & asw,               &
        & init = .TRUE._jplm)
  ENDIF
  THROWM(ERR.NE.0, 'allocation error for pccomp arrays')


  IF( CONTEXT%RTTOV_RUN_MODE(1:3) .EQ. 'PCK' ) THEN

    IF(CONTEXT%OPTS%RT_IR%PC%ADDRADREC) THEN
      ALLOCATE( PROFILES_K_REC( NCHANNELSREC ), EMISSIVITY_K( NCHANNELSREC ) , &
        &   REFLECTANCE_K( NCHANNELSREC ) ,&
            STAT= ERR)
      THROW(ERR.NE.0)
      CALL RTTOV_ALLOC_PROF(                                  &
          ERR, NCHANNELSREC, PROFILES_K_REC, NLEVELS,  CONTEXT%OPTS  , &
          ASW, COEFS=CONTEXT%RTH%COEFS, INIT = .TRUE. )
      THROWM(ERR.NE.0,"ALLOCATION OF PROFILES_K_REC")
    ELSE
      ALLOCATE( PROFILES_K_PC( NPCSCORES ), EMISSIVITY_K( NPCSCORES ) , &
        &   REFLECTANCE_K( NPCSCORES ) ,&
            STAT= ERR)
      THROW(ERR.NE.0)
      CALL RTTOV_ALLOC_PROF(                                  &
          ERR, NPCSCORES, PROFILES_K_PC, NLEVELS,  CONTEXT%OPTS  , &
          ASW, COEFS=CONTEXT%RTH%COEFS, INIT = .TRUE. )
      THROWM(ERR.NE.0,"ALLOCATION OF PROFILES_K_PC")
    ENDIF

    IF (CONTEXT%OPTS % RT_IR % PC % ADDRADREC) THEN
      CALL RTTOV_ALLOC_PCCOMP(  &
          & ERR,        &
          & PCCOMP_K,             &
          & NPCSCORES  ,  &
          & ASW,                &
          & INIT = .TRUE._JPLM, &
          & NCHANNELS_REC = NCHANNELSREC)
    ELSE
      CALL RTTOV_ALLOC_PCCOMP( &
          & ERR,       &
          & PCCOMP_K,            &
          & NPCSCORES , &
          & ASW,               &
          & INIT = .TRUE._JPLM)
    ENDIF
    THROWM(ERR.NE.0, 'allocation error for pccomp_k arrays')

    ALLOCATE( PROFILES_K( NCHANNELS ), EMISSIVITY_K( NCHANNELS ) , &
          &   REFLECTANCE_K( NCHANNELS ) ,&
              STAT= ERR)
    THROW(ERR.NE.0) 
    CALL RTTOV_ALLOC_PROF(                                  &
      ERR, NCHANNELS, PROFILES_K, NLEVELS,  CONTEXT%OPTS  , &
      ASW, COEFS=CONTEXT%RTH%COEFS, INIT = .TRUE. )
    THROWM(ERR.NE.0,"ALLOCATION OF PROFILES")

    CALL RTTOV_ALLOC_TRANSMISSION(             &
      ERR, TRANSMISSION_K, NLEVELS, NCHANNELS, &
      ASW, INIT = .TRUE. )
    THROW(ERR.NE.0)

    DO ICHAN = 1, NCHANNELS
      PROFILES_K(ICHAN)%P = 0_JPRB
    ENDDO


    CALL RTTOV_ALLOC_RAD( ERR, NCHANNELS, RADIANCE_K, NLEVELS, ASW, INIT=.TRUE.)
    THROW(ERR.NE.0)

    EMISSIVITY_K(:)%EMIS_IN  = 0._JPRB
    EMISSIVITY_K(:)%EMIS_OUT = 0._JPRB
    REFLECTANCE_K(:)%REFL_IN  = 0._JPRB
    REFLECTANCE_K(:)%REFL_OUT = 0._JPRB

    IF( CONTEXT%OPTS%RT_IR%PC%ADDRADREC )THEN
      IF ( CONTEXT%OPTS%RT_ALL%SWITCHRAD )THEN
        PCCOMP_K%BT_PCCOMP = 1._JPRB
      ELSE
        PCCOMP_K%TOTAL_PCCOMP = 1._JPRB
      ENDIF
    ELSE
      PCCOMP_K%TOTAL_PCSCORES = 1._JPRB
    ENDIF
      
      
  ENDIF


  INFO("-- RTTOV CONFIGURATION --")
  CALL RTTOV_PRINT_OPTS (CONTEXT%OPTS, 0, "USER OPTIONS")
  INFO("-------------------------")


! Run RTTOV
  IF( CONTEXT%RTTOV_RUN_MODE .EQ. 'PCDIRECT' ) THEN
  INFO("Start Direct")
    IF (CONTEXT%opts % rt_ir % pc % addradrec) THEN
      CALL RTTOV_DIRECT(                &
            & ERR,                      &
            & CHANPROF,                 &
            & CONTEXT%OPTS,             &
            & CONTEXT%PROFILES,         &
            & CONTEXT%RTH%COEFS,        &
            & TRANSMISSION,             &
            & RADIANCE,                 &
            & CALCEMIS    = CALCEMIS,   &
            & EMISSIVITY  = EMISSIVITY, &
            & pccomp       = pccomp,     &! inout computed PC scores
            & channels_rec = channelsrec) ! in    reconstructed channel list
    ELSE
      CALL RTTOV_DIRECT(                &
            & ERR,                      &
            & CHANPROF,                 &
            & CONTEXT%OPTS,             &
            & CONTEXT%PROFILES,         &
            & CONTEXT%RTH%COEFS,        &
            & TRANSMISSION,             &
            & RADIANCE,                 &
            & CALCEMIS    = CALCEMIS,   &
            & EMISSIVITY  = EMISSIVITY, &
            & pccomp       = pccomp   ) ! inout computed PC scores
    ENDIF
    THROWM(ERR.NE.0,"An error occured while running RTTOV Direct")
    INFO("End Direct")

    CALL RTTOV_HDF_SAVE( ERR, FILENAME_PC, '/PCCOMP', CREATE=.TRUE., &
       &  PCCOMP = PCCOMP )
    THROW(ERR.NE.0)
    CALL ADD_MISC_INFO( CONTEXT, FILENAME_PC, ERR )
    THROW(ERR.NE.0)

    IF ( DOTHERMAL  .AND. FILENAME_SURFACE.NE."" ) THEN
      CALL RTTOV_HDF_SAVE( ERR, FILENAME_SURFACE, '/EMISSIVITY', CREATE=.TRUE., &
         &  EMISSIVITY = EMISSIVITY)
      THROW(ERR.NE.0)
      CALL RTTOV_HDF_SAVE( ERR, FILENAME_SURFACE, '/EMISSIVITY', CREATE=.FALSE., &
         &  L1 = CALCEMIS, SNAME='CALCEMIS')
      THROW(ERR.NE.0)
    ENDIF
            
    IF ( FILENAME_SURFACE.NE."" ) THEN
      CALL ADD_MISC_INFO( CONTEXT, FILENAME_SURFACE, ERR )
      THROW(ERR.NE.0)
    ENDIF

  ELSE IF( CONTEXT%RTTOV_RUN_MODE(1:3) .EQ. 'PCK' ) THEN


  INFO("Start K")
      CALL RTTOV_K(  &
            & ERR,                             &
            & CHANPROF,                        &
            & CONTEXT%OPTS,                    &
            & CONTEXT%PROFILES,                &
            & PROFILES_K,                      &
            & CONTEXT%RTH%COEFS,               &
            & TRANSMISSION,                    &
            & TRANSMISSION_K,                  &
            & RADIANCE,                        &
            & RADIANCE_K,                      &
            & CALCEMIS       = CALCEMIS,       &
            & EMISSIVITY     = EMISSIVITY,     &
            & EMISSIVITY_K   = EMISSIVITY_K,   &
            & PCCOMP         = PCCOMP,         &
            & PCCOMP_K       = PCCOMP_K,       &
            & PROFILES_K_REC = PROFILES_K_REC, &
            & PROFILES_K_PC  = PROFILES_K_PC,  &
            & CHANNELS_REC   = CHANNELSREC )

    THROWM(ERR.NE.0,"An error occured while running RTTOV K")
    INFO("End K")

       CONTEXT%OPTS%RT_IR%PC%ADDPC = .FALSE.  ! force KMATRIX UNITS in HDF5
       CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/KMATRIX_PRED', CREATE=.TRUE., &
       &  KMATRIX = PROFILES_K, OPTS= CONTEXT%OPTS)
       THROWM(ERR.NE.0,"CANNOT CREATE FILE "//TRIM(FILENAME_KMAT))
       CONTEXT%OPTS%RT_IR%PC%ADDPC = .TRUE.

       IF( CONTEXT%OPTS%RT_IR%PC%ADDRADREC )THEN
         CONTEXT%OPTS%RT_IR%PC%ADDPC = .FALSE.  ! force KMATRIX UNITS in HDF5
         CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/KMATRIX', CREATE=.FALSE., &
         &  KMATRIX = PROFILES_K_REC, OPTS= CONTEXT%OPTS)
         THROW(ERR.NE.0)
         CONTEXT%OPTS%RT_IR%PC%ADDPC = .TRUE.
       ELSE
         CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/KMATRIX_K_PC', CREATE=.FALSE., &
         &  KMATRIX = PROFILES_K_PC, OPTS= CONTEXT%OPTS)
         THROW(ERR.NE.0)
       ENDIF
       
       CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/PCCOMP_K', CREATE=.FALSE., &
       &  PCCOMP = PCCOMP_K)
       THROW(ERR.NE.0)

       CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/PROFILES', CREATE=.FALSE., &
       &  PROFILES = CONTEXT%PROFILES)
       THROW(ERR.NE.0)
                 
       CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/CHANPROF', CREATE=.FALSE., &
       &  CHANPROF = CHANPROF)
       THROW(ERR.NE.0)
                 
       CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/EMISSIVITY_K', CREATE=.FALSE., &
       &  EMISSIVITY = EMISSIVITY_K)
       THROW(ERR.NE.0)
       CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/EMISSIVITY_K', CREATE=.FALSE., &
       &  L1 = CALCEMIS, SNAME='CALCEMIS')
       THROW(ERR.NE.0)
         
       CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/OPTIONS', CREATE=.FALSE., &
       &  OPTIONS = CONTEXT%OPTS)
       THROW(ERR.NE.0)
       
       CALL RTTOV_HDF_SAVE( ERR, FILENAME_KMAT, '/PCCOMP', CREATE=.FALSE., &
       &  PCCOMP = PCCOMP )
       THROW(ERR.NE.0)

       CALL ADD_MISC_INFO( CONTEXT, FILENAME_KMAT, ERR )
       THROW(ERR.NE.0)

  ELSE
        ERR = 1
        THROWM(ERR.NE.0,"CONTEXT%RTTOV_RUN_MODE unknown "//CONTEXT%RTTOV_RUN_MODE)
  ENDIF

 
! Deallocate

  ASW = 0_JPIM
  CALL RTTOV_ALLOC_RAD( ERR, NCHANNELS, RADIANCE, NLEVELS, ASW)
  THROW(ERR.NE.0)
  
  DEALLOCATE( CHANPROF, STAT= ERR)
  THROW(ERR.NE.0)
  
  DEALLOCATE( CALCEMIS, EMISSIVITY, STAT= ERR)
  THROW(ERR.NE.0)
  
 ! DEALLOCATE( CALCREFL, REFLECTANCE, STAT= ERR)
 ! THROW(ERR.NE.0)

  CALL RTTOV_ALLOC_TRANSMISSION( ERR, TRANSMISSION, NLEVELS, NCHANNELS,  &
      & ASW)  
  THROW(ERR.NE.0)
  
  
  IF(CONTEXT%OPTS%RT_IR%PC%ADDRADREC) THEN
    DEALLOCATE( CHANNELSREC, STAT= ERR)
    THROW(ERR.NE.0)
  ENDIF
  
  ! Deallocate pccomp arrays
  CALL rttov_alloc_pccomp(ERR, pccomp, npcscores, asw)
  THROW(ERR.NE.0)

  IF( CONTEXT%RTTOV_RUN_MODE(1:3) .EQ. 'KPC' ) THEN

    CALL RTTOV_ALLOC_PROF(                                  &
      ERR, NCHANNELS, PROFILES_K, NLEVELS,  CONTEXT%OPTS  , &
      ASW, COEFS=CONTEXT%RTH%COEFS)
    THROWM(ERR.NE.0,"DEALLOCATION OF PROFILES")

    CALL RTTOV_ALLOC_TRANSMISSION(             &
      ERR, TRANSMISSION_K, NLEVELS, NCHANNELS, &
      ASW )
    THROW(ERR.NE.0)

    CALL RTTOV_ALLOC_RAD( ERR, NCHANNELS, RADIANCE_K, NLEVELS, ASW)
    THROW(ERR.NE.0)
    DEALLOCATE( PROFILES_K, STAT= ERR)
    THROW(ERR.NE.0)
    
    DEALLOCATE( EMISSIVITY_K, REFLECTANCE_K, STAT= ERR)
    THROW(ERR.NE.0) 


    IF(CONTEXT%OPTS%RT_IR%PC%ADDRADREC) THEN
      CALL RTTOV_ALLOC_PROF(                                  &
          ERR, NCHANNELSREC, PROFILES_K_REC, NLEVELS,  CONTEXT%OPTS  , &
          ASW, COEFS=CONTEXT%RTH%COEFS, INIT = .TRUE. )
      THROWM(ERR.NE.0,"DE-ALLOCATION OF PROFILES_K_REC")

      DEALLOCATE( PROFILES_K_REC, EMISSIVITY_K , &
        &   REFLECTANCE_K ,&
            STAT= ERR)
      THROW(ERR.NE.0)
    ELSE
      CALL RTTOV_ALLOC_PROF(                                  &
          ERR, NPCSCORES, PROFILES_K_PC, NLEVELS,  CONTEXT%OPTS  , &
          ASW, COEFS=CONTEXT%RTH%COEFS, INIT = .TRUE. )
      THROWM(ERR.NE.0,"DEALLOCATION OF PROFILES_K_PC")

      DEALLOCATE( PROFILES_K_PC, EMISSIVITY_K, &
        &   REFLECTANCE_K ,&
            STAT= ERR)
      THROW(ERR.NE.0)
    ENDIF

    IF (CONTEXT%OPTS % RT_IR % PC % ADDRADREC) THEN
      CALL RTTOV_ALLOC_PCCOMP(  &
          & ERR,        &
          & PCCOMP_K,             &
          & NPCSCORES  ,  &
          & ASW,                &
          & INIT = .TRUE._JPLM, &
          & NCHANNELS_REC = NCHANNELSREC)
    ELSE
      CALL RTTOV_ALLOC_PCCOMP( &
          & ERR,       &
          & PCCOMP_K,            &
          & NPCSCORES , &
          & ASW,               &
          & INIT = .TRUE._JPLM)
    ENDIF
    THROWM(ERR.NE.0, 'deallocation error for pccomp_k arrays')

  ENDIF

  ! remember profiles are allocated/deallocated by rttovcontextload/drop 
  CALL RTTOVCONTEXTDROP( CONTEXT, ERR )
  THROW(ERR.NE.0)

  CALL CLOSE_HDF( ERR )
  THROW(ERR.NE.0)

  CATCH
  
  CALL CLOSE_HDF( ERR )
  ERR = 1
  
  RETURN


END SUBROUTINE


SUBROUTINE RTTOV_GUI_GET_COEF_VAL_I0( VARCH, I0, ERR)
!
#include "throw.h"

USE PARKIND1
!INTF_OFF
USE RTTOV_GUI_HANDLE
USE RTTOV_CONST
!INTF_ON
  IMPLICIT NONE
  CHARACTER(LEN=*),   INTENT(IN)  :: VARCH
  INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
  INTEGER(KIND=JPIM), INTENT(OUT) :: I0
!INTF_END

!
TRY

  ERR = 0

  IF( .NOT. ASSOCIATED( THE_RTH ) ) THEN
    ERR = 1_jpim
    THROWM(ERR.NE.0,"No coefficients loaded")
  ENDIF

  SELECT CASE (TRIM(VARCH))
    CASE ("ID_PLATFORM")
      I0 = THE_RTH%COEFS%COEF%ID_PLATFORM
    CASE ("ID_SAT")
      I0 = THE_RTH%COEFS%COEF%ID_SAT
    CASE ("ID_INST")
      I0 = THE_RTH%COEFS%COEF%ID_INST
    CASE ("ID_SENSOR")
      I0 = THE_RTH%COEFS%COEF%ID_SENSOR
    CASE ("ID_COMP_LVL")
      I0 = THE_RTH%COEFS%COEF%ID_COMP_LVL
    CASE ("ID_COMP_PC")
      I0 = THE_RTH%COEFS%COEF%ID_COMP_PC
    CASE ("FMV_MODEL_VER")
      I0 = THE_RTH%COEFS%COEF%FMV_MODEL_VER
    CASE ("FMV_CHN")
      I0 = THE_RTH%COEFS%COEF%FMV_CHN
    CASE ("FMV_GAS")
      I0 = THE_RTH%COEFS%COEF%FMV_GAS
    CASE ("NLEVELS")
      I0 = THE_RTH%COEFS%COEF%NLEVELS
    CASE ("NLAYERS")
      I0 = THE_RTH%COEFS%COEF%NLAYERS
    CASE ("NMIXED")
      I0 = THE_RTH%COEFS%COEF%NMIXED
    CASE ("NWATER")
      I0 = THE_RTH%COEFS%COEF%NWATER
    CASE ("NOZONE")
      I0 = THE_RTH%COEFS%COEF%NOZONE
    CASE ("NWVCONT")
      I0 = THE_RTH%COEFS%COEF%NWVCONT
    CASE ("NCO2")
      I0 = THE_RTH%COEFS%COEF%NCO2
    CASE ("NN2O")
      I0 = THE_RTH%COEFS%COEF%NN2O
    CASE ("NCO")
      I0 = THE_RTH%COEFS%COEF%NCO
    CASE ("NCH4")
      I0 = THE_RTH%COEFS%COEF%NCH4
    CASE ("NSO2")
      I0 = THE_RTH%COEFS%COEF%NSO2
    CASE ("INCZEEMAN")
      IF( THE_RTH%COEFS%COEF%INCZEEMAN ) THEN
        I0 = 1
      ELSE
        I0 = 0
      ENDIF
    CASE ("FMV_PC_BANDS")
      I0 = THE_RTH%COEFS%COEF_PCCOMP%FMV_PC_BANDS
    CASE ("FMV_PC_CLD")
      I0 = THE_RTH%COEFS%COEF_PCCOMP%FMV_PC_CLD
    CASE ("FMV_PC_MNUM")
      I0 = THE_RTH%COEFS%COEF_PCCOMP%FMV_PC_MNUM
    CASE ("SOLARCOEF")
      IF( THE_RTH%COEFS%COEF%SOLARCOEF ) THEN
        I0 = 1
      ELSE
        I0 = 0
      ENDIF
    CASE ("NLTECOEF")
      IF( THE_RTH%COEFS%COEF%NLTECOEF) THEN
        I0 = 1
      ELSE
        I0 = 0
      ENDIF
    CASE ("PMC_SHIFT")
      IF( THE_RTH%COEFS%COEF%PMC_SHIFT ) THEN
        I0 = 1
      ELSE
        I0 = 0
      ENDIF

! water optical constant 
! wave spectrum 
! fastem
! ssirem                 all omited
      
! aliases
    CASE ("NCHANNELS")
      I0 = THE_RTH%COEFS%COEF%FMV_CHN

! array sizes
    CASE ("SIZE_FMV_GAS_ID" ,&
        & "SIZE_FMV_VAR"    ,&
        & "SIZE_FMV_COE"    ,&
        & "SIZE_FMV_LVL"    )
      I0 = THE_RTH%COEFS%COEF%FMV_GAS
      
    CASE ("SIZE_FMV_GAS_POS")
      I0 = ngases_max
      
    CASE ("SIZE_FF_ORI_CHN" ,&
        & "SIZE_FF_VAL_CHN" ,&
        & "SIZE_FF_CWN" ,&
        & "SIZE_FF_BCO" ,&
        & "SIZE_FF_BCS" ,&
        & "SIZE_FF_GAM" ,&
        & "SIZE_TT_VAL_CHN" ,&
        & "SIZE_TT_A0"  ,&
        & "SIZE_TT_A1"  ,&
        & "SIZE_PW_VAL_CHN" ,&
        & "SIZE_SS_VAL_CHN" ,&
        & "SIZE_SS_SOLAR_SPECTRUM" )
      I0 = THE_RTH%COEFS%COEF%FMV_CHN
       
    CASE ("SIZE_REF_PRFL_P" ,&
        & "SIZE_LIM_PRFL_P" ,&
        & "SIZE_LIM_PRFL_TMAX" ,&
        & "SIZE_LIM_PRFL_TMIN" )
      I0 = THE_RTH%COEFS%COEF%NLEVELS
      
    CASE ("SIZE_NOISE_IN")
      I0=THE_RTH%COEFS%COEF_PCCOMP%FMV_PC_NCHN
    CASE ("SIZE_FMV_PC_SETS")
      I0 = THE_RTH%COEFS%COEF_PCCOMP%FMV_PC_BANDS

!    CASE ("SIZE_LINE_BY_LINE")
!      I0 = 20

! 2 dimensions arrays, should be specified in I0 and I1
    CASE ("SIZE_FMV_PC_NPRED")
      I0 = 2

! array size for aliases
    CASE ("SIZE_WAVENUMBERS")
      I0 = THE_RTH%COEFS%COEF%FMV_CHN
              
    CASE ("SIZE_REF_PRESSURE" ,&      
        & "SIZE_REF_TEMPERATURE" ,&      
        & "SIZE_REF_WATERVAPOR" ,&
        & "SIZE_REF_OZONE" ,&
        & "SIZE_REF_CO2" ,&
        & "SIZE_REF_N2O" ,&
        & "SIZE_REF_CO" ,&
        & "SIZE_REF_CH4" ,&
        & "SIZE_REF_SO2" ,&
        & "SIZE_MIN_TEMPERATURE" ,&
        & "SIZE_MIN_WATERVAPOR" ,&
        & "SIZE_MIN_OZONE" ,&
        & "SIZE_MIN_CO2" ,&
        & "SIZE_MIN_N2O" ,&
        & "SIZE_MIN_CO" ,&
        & "SIZE_MIN_CH4" ,&
        & "SIZE_MIN_SO2" ,&
        & "SIZE_MAX_TEMPERATURE" ,&
        & "SIZE_MAX_WATERVAPOR" ,&
        & "SIZE_MAX_OZONE" ,&
        & "SIZE_MAX_CO2" ,&
        & "SIZE_MAX_N2O" ,&
        & "SIZE_MAX_CO" ,&
        & "SIZE_MAX_CH4" ,&
        & "SIZE_MAX_SO2")
      I0 = THE_RTH%COEFS%COEF%NLEVELS

    CASE DEFAULT
      ERR = 1_JPIM
      THROWM(ERR.NE.0,"INTERGER SCALAR VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH
  
END SUBROUTINE

SUBROUTINE RTTOV_GUI_GET_COEF_VAL_R0( VARCH, R0, ERR)
!
#include "throw.h"

USE PARKIND1
!INTF_OFF
USE RTTOV_GUI_HANDLE
USE RTTOV_CONST
!INTF_ON
  IMPLICIT NONE
  CHARACTER(LEN=*),   INTENT(IN)  :: VARCH
  INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
  REAL(KIND=JPRB),    INTENT(OUT) :: R0
!INTF_END

!
TRY

  ERR = 0

  IF( .NOT. ASSOCIATED( THE_RTH ) ) THEN
    ERR = 1_jpim
    THROWM(ERR.NE.0,"No coefficients loaded")
  ENDIF

  SELECT CASE (TRIM(VARCH))
  
    CASE ("FC_SPEEDL")
      R0 = SPEEDL
    CASE ("FC_PLANCK_C1")
      R0 = THE_RTH%COEFS%COEF%FC_PLANCK_C1
    CASE ("FC_PLANCK_C2")
      R0 = THE_RTH%COEFS%COEF%FC_PLANCK_C2
    CASE ("FC_SAT_HEIGHT")
      R0 = THE_RTH%COEFS%COEF%FC_SAT_HEIGHT

    CASE DEFAULT
      ERR = 1_JPIM
      THROWM(ERR.NE.0,"REAL SCALAR VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH
  
END SUBROUTINE

SUBROUTINE RTTOV_GUI_GET_COEF_VAL_C0( VARCH, C0, ERR)
!
#include "throw.h"

USE PARKIND1
!INTF_OFF
USE RTTOV_GUI_HANDLE
USE RTTOV_CONST
!INTF_ON
  IMPLICIT NONE
  CHARACTER(LEN=*),   INTENT(IN)  :: VARCH
  INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
  CHARACTER(LEN=80),  INTENT(OUT) :: C0
!INTF_END

!
TRY

  ERR = 0

  IF( .NOT. ASSOCIATED( THE_RTH ) ) THEN
    ERR = 1_jpim
    THROWM(ERR.NE.0,"No coefficients loaded")
  ENDIF

  SELECT CASE (TRIM(VARCH))
  
    CASE ("ID_CREATION")
      C0 = THE_RTH%COEFS%COEF%ID_CREATION
    CASE ("ID_COMMON_NAME")
      C0 = THE_RTH%COEFS%COEF%ID_COMMON_NAME
    CASE ("FMV_MODEL_DEF")
      C0 = THE_RTH%COEFS%COEF%FMV_MODEL_DEF

    CASE DEFAULT
      ERR = 1_JPIM
      THROWM(ERR.NE.0,"CHARACTER STRING VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH
  
END SUBROUTINE

SUBROUTINE RTTOV_GUI_GET_COEF_VAL_I1( VARCH, M, I1, ERR)
!
#include "throw.h"

USE PARKIND1
!INTF_OFF
USE RTTOV_GUI_HANDLE
USE RTTOV_CONST
!INTF_ON
  IMPLICIT NONE
  CHARACTER(LEN=*),   INTENT(IN)  :: VARCH
  INTEGER(KIND=JPIM), INTENT(IN)  :: M

  INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
  INTEGER(KIND=JPIM), INTENT(OUT) :: I1(M)

!INTF_END

!
TRY

  ERR = 0

  IF( .NOT. ASSOCIATED( THE_RTH ) ) THEN
    ERR = 1_jpim
    THROWM(ERR.NE.0,"No coefficients loaded")
  ENDIF

  SELECT CASE (TRIM(VARCH))
    CASE ("FMV_GAS_ID")
      I1 = THE_RTH%COEFS%COEF%FMV_GAS_ID
    CASE ("FMV_GAS_POS")
      I1 = THE_RTH%COEFS%COEF%FMV_GAS_POS
    CASE ("FMV_VAR")
      I1 = THE_RTH%COEFS%COEF%FMV_VAR
    CASE ("FMV_COE")
      I1 = THE_RTH%COEFS%COEF%FMV_COE
    CASE ("FMV_LVL")
      I1 = THE_RTH%COEFS%COEF%FMV_LVL
    CASE ("FF_ORI_CHN")
      I1 = THE_RTH%COEFS%COEF%FF_ORI_CHN
    CASE ("FF_VAL_CHN")
      I1 = THE_RTH%COEFS%COEF%FF_VAL_CHN
    CASE ("TT_VAL_CHN")
      I1 = THE_RTH%COEFS%COEF%TT_VAL_CHN
    CASE ("PW_VAL_CHN")
      I1 = THE_RTH%COEFS%COEF%PW_VAL_CHN
    CASE ("SS_VAL_CHN")
      I1 = THE_RTH%COEFS%COEF%SS_VAL_CHN
    CASE ("FMV_PC_SETS")
      I1 = THE_RTH%COEFS%COEF_PCCOMP%FMV_PC_SETS

    CASE ("SIZE_FMV_PC_NPRED")
      I1 = (/THE_RTH%COEFS%COEF_PCCOMP%FMV_PC_BANDS ,THE_RTH%COEFS%COEF_PCCOMP%FMV_PC_SETS(1)/)

    CASE DEFAULT
      ERR = 1_JPIM
      THROWM(ERR.NE.0,"INTEGER ARRAY VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH
  
END SUBROUTINE

SUBROUTINE RTTOV_GUI_GET_COEF_VAL_I2( VARCH, M, N, I2, ERR)
!
#include "throw.h"

USE PARKIND1
!INTF_OFF
USE RTTOV_GUI_HANDLE
USE RTTOV_CONST
!INTF_ON
  IMPLICIT NONE
  CHARACTER(LEN=*),   INTENT(IN)  :: VARCH
  INTEGER(KIND=JPIM), INTENT(IN)  :: M
  INTEGER(KIND=JPIM), INTENT(IN)  :: N

  INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
  INTEGER(KIND=JPIM), INTENT(OUT) :: I2(M, N)

!INTF_END

!
TRY

  ERR = 0

  IF( .NOT. ASSOCIATED( THE_RTH ) ) THEN
    ERR = 1_jpim
    THROWM(ERR.NE.0,"No coefficients loaded")
  ENDIF

  SELECT CASE (TRIM(VARCH))

    CASE ("FMV_PC_NPRED")
      I2(:,:) = THE_RTH%COEFS%COEF_PCCOMP%PCREG(:,:)%FMV_PC_NPRED

    CASE DEFAULT
      ERR = 1_JPIM
      THROWM(ERR.NE.0,"INTEGER ARRAY VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH
  
END SUBROUTINE

SUBROUTINE RTTOV_GUI_GET_COEF_VAL_R1( VARCH, M, R1, ERR)
!
#include "throw.h"

USE PARKIND1
!INTF_OFF
USE RTTOV_GUI_HANDLE
USE RTTOV_CONST
!INTF_ON
  IMPLICIT NONE
  CHARACTER(LEN=*),   INTENT(IN)  :: VARCH
  INTEGER(KIND=JPIM), INTENT(IN)  :: M
  INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
  REAL(KIND=JPRB),    INTENT(OUT) :: R1(M)
!INTF_END

  INTEGER(KIND=JPIM) :: N

!
TRY

  ERR = 0


  IF( .NOT. ASSOCIATED( THE_RTH ) ) THEN
    ERR = 1_jpim
    THROWM(ERR.NE.0,"No coefficients loaded")
  ENDIF

  SELECT CASE (TRIM(VARCH))
  
    CASE ("FF_CWN")
      R1 = THE_RTH%COEFS%COEF%FF_CWN
    CASE ("FF_BCO")
      R1 = THE_RTH%COEFS%COEF%FF_BCO
    CASE ("FF_BCS")
      R1 = THE_RTH%COEFS%COEF%FF_BCS
    CASE ("FF_GAM")
      R1 = THE_RTH%COEFS%COEF%FF_GAM
    CASE ("TT_A0")
      R1 = THE_RTH%COEFS%COEF%TT_A0
    CASE ("TT_A1")
      R1 = THE_RTH%COEFS%COEF%TT_A1
    CASE ("SS_SOLAR_SPECTRUM")
      R1 = THE_RTH%COEFS%COEF%SS_SOLAR_SPECTRUM
    CASE ("REF_PRFL_P")
      R1 = THE_RTH%COEFS%COEF%REF_PRFL_P
    CASE ("LIM_PRFL_P")
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_P
    CASE ("LIM_PRFL_TMAX")
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_TMAX
    CASE ("LIM_PRFL_TMIN")
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_TMIN
    CASE ("NOISE_IN")
      R1 = THE_RTH%COEFS%COEF_PCCOMP%NOISE_IN

! Aliases
    CASE ("WAVENUMBERS")
      R1 = THE_RTH%COEFS%COEF%FF_CWN
      
    CASE ("REF_PRESSURE")
      R1 = THE_RTH%COEFS%COEF%REF_PRFL_P
      
    CASE ("REF_TEMPERATURE")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_MIXED)
      R1 = THE_RTH%COEFS%COEF%REF_PRFL_T(:,N)
      
    CASE ("REF_WATERVAPOR")
      IF( SIZE(R1) .NE. THE_RTH%COEFS%COEF%NLEVELS) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"REAL ARRAY ARGUMENT DOES NOT HAVE THE CORRECT DIMENSION")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_WATERVAPOUR)
      R1 = THE_RTH%COEFS%COEF%REF_PRFL_MR(:,N)
      
     CASE ("REF_OZONE")
      IF( THE_RTH%COEFS%COEF%NOZONE .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN OZONE")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_OZONE)
      R1 = THE_RTH%COEFS%COEF%REF_PRFL_MR(:,N)
        
     CASE ("REF_CO2")
      IF( THE_RTH%COEFS%COEF%NCO2 .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO2")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_CO2)
      R1 = THE_RTH%COEFS%COEF%REF_PRFL_MR(:,N)
      
     CASE ("REF_N2O")
      IF( THE_RTH%COEFS%COEF%NN2O .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN N2O")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_N2O)
      R1 = THE_RTH%COEFS%COEF%REF_PRFL_MR(:,N)
      
     CASE ("REF_CO")
      IF( THE_RTH%COEFS%COEF%NCO .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_CO)
      R1 = THE_RTH%COEFS%COEF%REF_PRFL_MR(:,N)
      
     CASE ("REF_CH4")
      IF( THE_RTH%COEFS%COEF%NCH4 .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CH4")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_CH4)
      R1 = THE_RTH%COEFS%COEF%REF_PRFL_MR(:,N)

     CASE ("REF_SO2")
      IF( THE_RTH%COEFS%COEF%NSO2 .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN SO2")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_SO2)
      R1 = THE_RTH%COEFS%COEF%REF_PRFL_MR(:,N)


    CASE ("MIN_TEMPERATURE")
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_TMIN(:)
      
    CASE ("MIN_WATERVAPOR")
      IF( SIZE(R1) .NE. THE_RTH%COEFS%COEF%NLEVELS) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"REAL ARRAY ARGUMENT DOES NOT HAVE THE CORRECT DIMENSION")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_WATERVAPOUR)
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_GMIN(:,N)
      
     CASE ("MIN_OZONE")
      IF( THE_RTH%COEFS%COEF%NOZONE .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN OZONE")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_OZONE)
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_GMIN(:,N)
        
     CASE ("MIN_CO2")
      IF( THE_RTH%COEFS%COEF%NCO2 .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO2")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_CO2)
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_GMIN(:,N)
      
     CASE ("MIN_N2O")
      IF( THE_RTH%COEFS%COEF%NN2O .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN N2O")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_N2O)
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_GMIN(:,N)
      
     CASE ("MIN_CO")
      IF( THE_RTH%COEFS%COEF%NCO .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_CO)
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_GMIN(:,N)
      
     CASE ("MIN_CH4")
      IF( THE_RTH%COEFS%COEF%NCH4 .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CH4")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_CH4)
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_GMIN(:,N)

     CASE ("MIN_SO2")
      IF( THE_RTH%COEFS%COEF%NSO2 .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN SO2")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_SO2)
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_GMIN(:,N)


    CASE ("MAX_TEMPERATURE")
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_TMAX(:)
      
    CASE ("MAX_WATERVAPOR")
      IF( SIZE(R1) .NE. THE_RTH%COEFS%COEF%NLEVELS) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"REAL ARRAY ARGUMENT DOES NOT HAVE THE CORRECT DIMENSION")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_WATERVAPOUR)
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_GMAX(:,N)
      
     CASE ("MAX_OZONE")
      IF( THE_RTH%COEFS%COEF%NOZONE .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN OZONE")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_OZONE)
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_GMAX(:,N)
        
     CASE ("MAX_CO2")
      IF( THE_RTH%COEFS%COEF%NCO2 .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO2")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_CO2)
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_GMAX(:,N)
      
     CASE ("MAX_N2O")
      IF( THE_RTH%COEFS%COEF%NN2O .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN N2O")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_N2O)
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_GMAX(:,N)
      
     CASE ("MAX_CO")
      IF( THE_RTH%COEFS%COEF%NCO .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CO")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_CO)
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_GMAX(:,N)
      
     CASE ("MAX_CH4")
      IF( THE_RTH%COEFS%COEF%NCH4 .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN CH4")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_CH4)
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_GMAX(:,N)

     CASE ("MAX_SO2")
      IF( THE_RTH%COEFS%COEF%NSO2 .LE. 0 ) THEN
        ERR = 1_JPIM
      ENDIF
      THROWM(ERR.NE.0,"DOES NOT CONTAIN SO2")
      N  = THE_RTH%COEFS%COEF%FMV_GAS_POS(GAS_ID_SO2)
      R1 = THE_RTH%COEFS%COEF%LIM_PRFL_GMAX(:,N)

    CASE DEFAULT
      ERR = 1_JPIM
      THROWM(ERR.NE.0,"REAL ARRAY VARIABLE "//TRIM(VARCH)//" NOT FOUND")

  END SELECT

CATCH
  
END SUBROUTINE

SUBROUTINE RTTOV_GUI_GET_SIZE_EIGENVECTOR( IBAND, I0, ERR)
!
#include "throw.h"

USE PARKIND1
!INTF_OFF
USE RTTOV_GUI_HANDLE
USE RTTOV_CONST
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=JPIM), INTENT(IN)  :: IBAND
  INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
  INTEGER(KIND=JPIM), INTENT(OUT) :: I0
!INTF_END


!
TRY

  ERR = 0


  IF( .NOT. ASSOCIATED( THE_RTH ) ) THEN
    ERR = 1_jpim
    THROWM(ERR.NE.0,"No coefficients loaded")
  ENDIF
  
  I0 = SIZE(THE_RTH%COEFS%COEF_PCCOMP%EIGEN(IBAND)%EIGENVECTORS, 1_JPIM)

CATCH
  
END SUBROUTINE

!SUBROUTINE RTTOV_GUI_GET_EIGENVECTOR( IBAND, IVECT, ICHAN, M, R1, ERR)
SUBROUTINE RTTOV_GUI_GET_EIGENVECTOR( IBAND, IVECT, M, R1, ERR)

#include "throw.h"

USE PARKIND1
!INTF_OFF
USE RTTOV_GUI_HANDLE
USE RTTOV_CONST
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=JPIM), INTENT(IN)  :: IBAND
  INTEGER(KIND=JPIM), INTENT(IN)  :: IVECT
  !INTEGER(KIND=JPIM), INTENT(IN)  :: ICHAN
  INTEGER(KIND=JPIM), INTENT(IN)  :: M
  INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
  REAL(KIND=JPRB),    INTENT(OUT) :: R1(M)
!INTF_END


!
TRY

  ERR = 0


  IF( .NOT. ASSOCIATED( THE_RTH ) ) THEN
    ERR = 1_jpim
    THROWM(ERR.NE.0,"No coefficients loaded")
  ENDIF
  
  !IF( ICHAN .GT. 0 ) THEN
  !  R1 = THE_RTH%COEFS%COEF_PCCOMP%EIGEN(IBAND)%EIGENVECTORS(ICHAN,:)
  !ELSE
    R1 = THE_RTH%COEFS%COEF_PCCOMP%EIGEN(IBAND)%EIGENVECTORS(:,IVECT)
  !ENDIF
  
CATCH
  
END SUBROUTINE


SUBROUTINE RTTOV_GUI_AER_CLIM_PROF( ML, P,T,Q,GAS_UNITS,MMR_AER,LEVSURF, &
           & LATITUDE,ELEVATION,SCALEFACTOR,AERPROF,ERR)
!
#include "throw.h"

USE PARKIND1
!INTF_OFF
USE RTTOV_GUI_HANDLE
USE RTTOV_CONST, ONLY : NAER_OPAC
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=JPIM), INTENT(IN)  :: ML
  INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
  REAL(KIND=JPRB),    INTENT(IN)  :: P(ML)
  REAL(KIND=JPRB),    INTENT(IN)  :: T(ML)
  REAL(KIND=jprb),    INTENT(IN)  :: Q(ML)               ! Input water vapour profile
  INTEGER(KIND=jpim), INTENT(IN)  :: GAS_UNITS           ! Input water vapour units
  INTEGER(KIND=jpim), INTENT(IN)  :: MMR_AER             ! Output aerosol units: 0=>cm-3; otherwise=>kg/kg
  INTEGER(KIND=jpim), INTENT(IN)  :: LEVSURF             ! Surface level
  REAL(KIND=jprb),    INTENT(IN)  :: LATITUDE            ! Latitude for profiles
  REAL(KIND=jprb),    INTENT(IN)  :: ELEVATION           ! Surface elevation (km)
  REAL(KIND=jprb),    INTENT(IN)  :: SCALEFACTOR         ! Factor by which to scale profiles
  REAL(KIND=jprb),    INTENT(OUT) :: AERPROF(ML-1,10,13) ! Output aerosol profiles (nlayer, nclim, naertype)

!INTF_END

  LOGICAL(KIND=jplm) :: LMMR_AER
#include "rttov_aer_clim_prof.interface"

!
TRY

  ERR = 0

  IF( naer_opac .NE. 13_JPIM ) THEN
    ERR = 1_JPIM
    THROWM(ERR.NE.0,"INCONSISTENCY  NAER_OPAC NOT EQUAL TO 13")
  ENDIF
  LMMR_AER = (MMR_AER .NE. 0)
  CALL RTTOV_AER_CLIM_PROF( P,T,Q,GAS_UNITS,LMMR_AER,LEVSURF,LATITUDE,ELEVATION,SCALEFACTOR,AERPROF)

CATCH
  
END SUBROUTINE


