! Description:
!> @file
!!   Subroutines for the GUI
!
!> @brief
!!   Subroutines for the GUI
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
MODULE RTTOV_GUI_CONTEXT

#include "throw.h"

USE PARKIND1
USE RTTOV_GUI_HANDLE
USE RTTOV_HDF_MOD
USE RTTOV_CONST, ONLY : &
    PLATFORM_NAME,      &
    INST_NAME
USE RTTOV_TYPES

IMPLICIT NONE

TYPE RTTOVGUICONTEXT_TYPE
  TYPE(RTTOVGUIHANDLE_TYPE), POINTER :: RTH
  TYPE(RTTOV_OPTIONS)             :: opts
  TYPE(RTTOV_PROFILE), POINTER    :: profiles(:)
  CHARACTER(LEN=16)               :: RTTOV_RUN_MODE
END TYPE

CONTAINS

SUBROUTINE RTTOVCONTEXTDROP( SELF, ERR )
  TYPE(RTTOVGUICONTEXT_TYPE), INTENT(INOUT) :: SELF
  INTEGER,                 INTENT(OUT)   :: ERR
  
#include "rttov_alloc_prof.interface"

  TYPE(RTTOV_OPTIONS)             :: default_opts
TRY
  ERR = 0
  ! do not drop Rttov Handle RTH
  !
  IF( ASSOCIATED(SELF%PROFILES)) THEN
    CALL rttov_alloc_prof( &
            & ERR,      &
            & 1_JPIM,    &
            & SELF%profiles, &
            & SELF%PROFILES(1_JPIM)%nlevels,  &
            & SELF%opts,     &
            & 0_JPIM)
    THROWM(ERR.NE.0,"DEALLOCATION OF PROFILES")
    DEALLOCATE( SELF%PROFILES )
    
    SELF%opts = default_opts
    SELF%RTTOV_RUN_MODE = ""
    
  ENDIF
CATCH
END SUBROUTINE


SUBROUTINE RTTOVCONTEXTINITFROMPROFILE( SELF, FILENAME, RTH, ERR )
  TYPE(RTTOVGUICONTEXT_TYPE),          INTENT(OUT) :: SELF
  CHARACTER(LEN=*),                 INTENT(IN)  :: FILENAME
  TYPE(RTTOVGUIHANDLE_TYPE), TARGET,   INTENT(IN)  :: RTH
  INTEGER(KIND=JPIM),                          INTENT(OUT) :: ERR
!

#include "rttov_hdf_load.interface"

TRY
  ERR = 0_JPIM
  
  SELF%RTH => RTH
  SELF%RTTOV_RUN_MODE = ' '

  ALLOCATE(SELF%profiles(1_JPIM) , STAT = ERR )
  THROW(ERR.NE.0)
  
  CALL RTTOV_HDF_LOAD( ERR, FILENAME, "/PROFILES/0001", PROFILE=SELF%profiles(1_JPIM) )
  THROWM(ERR.NE.0,"Cannot load profile from "//TRIM(FILENAME))

  CALL RTTOV_HDF_LOAD( ERR, FILENAME, "/OPTIONS", OPTIONS=SELF%opts )
  THROWM(ERR.NE.0,"Cannot load options from "//TRIM(FILENAME))

  CALL COMPARELEVELS( SELF )
  
CATCH    

END SUBROUTINE




SUBROUTINE COMPARELEVELS( SELF )
  TYPE(RTTOVGUICONTEXT_TYPE), INTENT(INOUT) :: SELF
!
  INTEGER :: I
 
  IF( ASSOCIATED( SELF%RTH ) ) THEN

    IF( SELF%PROFILES(1_JPIM)%NLEVELS .EQ. SELF%RTH%COEFS%COEF%NLEVELS ) THEN
      SELF%OPTS%INTERPOLATION%ADDINTERP = .FALSE.
      DO I = 1, SELF%PROFILES(1_JPIM)%NLEVELS
        IF( ABS( SELF%PROFILES(1_JPIM)%P( I ) - SELF%RTH%COEFS%COEF%REF_PRFL_P( I ) ) &
        &    / ( SELF%PROFILES(1_JPIM)%P( I ) + SELF%RTH%COEFS%COEF%REF_PRFL_P( I ) ) &
        &  .GT. 1e-4 ) THEN
          SELF%OPTS%INTERPOLATION%ADDINTERP = .TRUE.
        ENDIF
      ENDDO
    ELSE
      SELF%OPTS%INTERPOLATION%ADDINTERP = .TRUE.
    ENDIF
  
  ELSE 
  
    SELF%OPTS%INTERPOLATION%ADDINTERP = .FALSE.
   
  ENDIF
  
END SUBROUTINE


SUBROUTINE ADD_MISC_INFO( SELF, H5FILE, ERR )
  TYPE(RTTOVGUICONTEXT_TYPE), INTENT(IN)  :: SELF
  CHARACTER(LEN=*),           INTENT(IN)  :: H5FILE
  INTEGER,                    INTENT(OUT) :: ERR
#include "rttov_hdf_save.interface" 
 CHARACTER(LEN=128) :: SAT
!
TRY
  
  ERR = 0
       WRITE( SAT, "(A,I2)" ) TRIM(PLATFORM_NAME(SELF%RTH%COEFS%COEF%ID_PLATFORM))//"-", &
          & SELF%RTH%COEFS%COEF%ID_SAT
          
       CALL RTTOV_HDF_SAVE( ERR, H5FILE,     '/MISC', CREATE=.false., &
       &  I0 = SELF%RTH%coefs%coef%fmv_chn, SNAME='NCHANNELS' )
       THROW(err.ne.0)

       CALL RTTOV_HDF_SAVE( ERR, H5FILE,     '/MISC', CREATE=.false., &
       &  R1 = SELF%RTH%coefs%coef%ff_cwn, SNAME='WAVENUMBERS', UNITS='cm-1' )
       THROW(err.ne.0)

       CALL RTTOV_HDF_SAVE( ERR, H5FILE,     '/MISC', CREATE=.false., &
       &  c0 = TRIM(SAT), SNAME='SATELLITE')
       THROW(err.ne.0)

       CALL RTTOV_HDF_SAVE( ERR, H5FILE,     '/MISC', CREATE=.false., &
       &  c0 = TRIM(INST_NAME(SELF%RTH%COEFS%COEF%ID_INST)), SNAME='INSTRUMENT')
       THROW(err.ne.0)

       CALL RTTOV_HDF_SAVE( ERR, H5FILE,     '/MISC', CREATE=.false., &
       &  c0 = SELF%RTH%coefs%coef%ID_COMMON_NAME, SNAME='ID_COMMON_NAME')
       THROW(err.ne.0)
       
       CALL RTTOV_HDF_SAVE( ERR, H5FILE,     '/MISC', CREATE=.false., &
       &  c0 = SELF%RTH%file_coef, SNAME='COEF_FILENAME')
       THROW(err.ne.0)

       CALL RTTOV_HDF_SAVE( ERR, H5FILE,     '/MISC', CREATE=.false., &
       &  c0 = SELF%RTH%file_scaer, SNAME='SCAT_COEF_FILENAME')
       THROW(err.ne.0)

       CALL RTTOV_HDF_SAVE( ERR, H5FILE,     '/MISC', CREATE=.false., &
       &  c0 = SELF%RTH%file_sccld, SNAME='CLOUD_COEF_FILENAME')
       THROW(err.ne.0)

       CALL RTTOV_HDF_SAVE( ERR, H5FILE,     '/MISC', CREATE=.false., &
       &  c0 = SELF%RTH%file_pccoef, SNAME='PC_COEF_FILENAME')
       THROW(err.ne.0)

CATCH

END SUBROUTINE


END
