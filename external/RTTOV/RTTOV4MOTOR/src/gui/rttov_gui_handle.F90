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
MODULE RTTOV_GUI_HANDLE

#include "throw.h"

USE PARKIND1
USE RTTOV_TYPES
USE MOD_RTTOV_EMIS_ATLAS, ONLY : RTTOV_EMIS_ATLAS_DATA
USE MOD_RTTOV_BRDF_ATLAS, ONLY : RTTOV_BRDF_ATLAS_DATA

IMPLICIT NONE

#include "rttov_errorreport.interface"

TYPE RTTOVGUIHANDLE_TYPE
  TYPE(RTTOV_COEFS)  :: COEFS
  Character(len=256) :: file_coef
  Character(len=256) :: file_scaer
  Character(len=256) :: file_sccld
  Character(len=256) :: file_mfasis_cld
  Character(len=256) :: file_pccoef
  TYPE(RTTOV_OPTIONS) :: opts
  TYPE(RTTOV_EMIS_ATLAS_DATA) :: emis_atlas
  TYPE(RTTOV_BRDF_ATLAS_DATA) :: brdf_atlas
  Integer(Kind=jpim)  :: ems_atlas_month = -1_jpim
  Integer(Kind=jpim)  :: brdf_atlas_month = -1_jpim
  Integer(Kind=jpim)  :: emis_atlas_id = -1_jpim
  Integer(Kind=jpim)  :: brdf_atlas_id = -1_jpim
END TYPE RTTOVGUIHANDLE_TYPE


TYPE(RTTOVGUIHANDLE_TYPE), POINTER :: THE_RTH => NULL()

CONTAINS

SUBROUTINE RTTOV_GUI_HANDLELOAD( RTH, channels, file_coef, file_scaer, file_sccld, file_mfasis_cld, file_pccoef, ERR )

  IMPLICIT NONE


  TYPE(RTTOVGUIHANDLE_TYPE), POINTER :: RTH
  INTEGER(KIND=JPIM), intent(out) :: ERR
  Character(len=*),   Intent(in)  :: file_coef
  Character(len=*),   Intent(in)  :: file_scaer
  Character(len=*),   Intent(in)  :: file_sccld
  Character(len=*),   Intent(in)  :: file_mfasis_cld
  Character(len=*),   Intent(in)  :: file_pccoef
  Integer(Kind=jpim), Intent(in)  :: channels(:)

  
#include "rttov_read_coefs.interface"
! #include "rttov_user_options_checkinput.interface"

  !Character(len=32)  :: form_coef
  LOGICAL(KIND=JPLM) :: allchannels
TRY

  ALLOCATE( RTH, STAT = ERR )
  THROW(ERR.NE.0)
  
  RTH%file_coef   = file_coef
  RTH%file_scaer  = file_scaer
  RTH%file_sccld  = file_sccld
  RTH%file_pccoef = file_pccoef
  RTH%file_mfasis_cld = file_mfasis_cld
  
  RTH%opts%rt_ir%addaerosl = file_scaer  .NE. ""
  RTH%opts%rt_ir%addclouds = file_sccld  .NE. ""
  RTH%opts%rt_ir%pc%addpc  = file_pccoef .NE. ""

  IF (file_mfasis_cld .NE. "") THEN
     RTH%opts%rt_ir%vis_scatt_model = 3
  ENDIF
  
  allchannels = SIZE(channels) .EQ. 1 .AND. channels(1) .LE. 0
  
  IF( allchannels ) THEN
      CALL rttov_read_coefs( err, RTH%coefs, RTH%opts, &
                             file_coef =file_coef,     &
                             file_scaer=file_scaer,    &
                             file_sccld=file_sccld,    &
                             file_mfasis_cld=file_mfasis_cld, &
                             file_pccoef=file_pccoef    )
  ELSE
    CALL rttov_read_coefs( err, RTH%coefs, RTH%opts,   &
                             file_coef =file_coef,     &
                             file_scaer=file_scaer,    &
                             file_sccld=file_sccld,    &
                             file_mfasis_cld=file_mfasis_cld, &
                             file_pccoef=file_pccoef,  &
                             channels = channels        )

  ENDIF

  THROWM(ERR.NE.0, "ERROR READING COEFFICIENT FILES")
  
  ! ne pas teser ici car les options ne sont pas toujours a jour
  ! en particulier si on commence a charger des coefs PC avant
  ! avoir positionne ops%addpc=true et ops%ipcbnd
  !

  !CALL rttov_user_options_checkinput( &
  !     & err,                         &
  !     & RTH%opts,                    &
  !     & RTH%coefs               )
  !THROWM(ERR.NE.0, "ERROR COMPATIBILITY OPTIONS AND COEFS")

CATCH
    DEALLOCATE( RTH )
    NULLIFY( RTH )

END SUBROUTINE


SUBROUTINE RTTOV_GUI_HANDLEDROP( RTH, ERR )

  IMPLICIT NONE

  TYPE(RTTOVGUIHANDLE_TYPE), POINTER :: RTH
  INTEGER(KIND=JPIM) :: ERR
#include "rttov_deallocate_emis_atlas.interface"
#include "rttov_deallocate_brdf_atlas.interface"
#include "rttov_dealloc_coefs.interface"

TRY

  ERR = 0_JPIM

  ! In case emissivity or brdf atlas have been loaded
  If( RTH%ems_atlas_month .GT. 0_jpim) Then
    CALL rttov_deallocate_emis_atlas(RTH%emis_atlas)
    RTH%ems_atlas_month = -1_jpim
    RTH%emis_atlas_id   = -1_jpim
  EndIf
  
  If( RTH%brdf_atlas_month .GT. 0_jpim) Then
    CALL rttov_deallocate_brdf_atlas(RTH%brdf_atlas)
    RTH%brdf_atlas_month = -1_jpim 
    RTH%brdf_atlas_id    = -1_jpim
  EndIf
  
  CALL RTTOV_DEALLOC_COEFS( ERR, RTH%COEFS )
  THROWM(ERR.NE.0, "ERROR DE-ALLOC COEFFICIENTS")

  DEALLOCATE( RTH, STAT = ERR )
  THROWM(ERR.NE.0, "ERROR DE-ALLOC RTTOV HANDLE TYPE")



CATCH

END SUBROUTINE


END
