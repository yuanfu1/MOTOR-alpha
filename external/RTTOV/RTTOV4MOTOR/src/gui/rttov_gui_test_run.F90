! Description:
!> @file
!!   Test program for the GUI.
!
!> @brief
!!   Test program for the GUI.
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
PROGRAM MAIN

#include "throw.h"

  USE PARKIND1
  USE rttov_unix_env


  IMPLICIT NONE

#include "rttov_errorreport.interface"
!#include "rttov_gui_load.interface"
!#include "rttov_gui_run.interface"
!#include "rttov_gui_get_emisbrdf.interface"
!#include "rttov_gui_drop.interface"

  CHARACTER(LEN=128)      :: MODE = 'DIRECT'
  INTEGER(KIND=JPIM)            :: ERR

  CHARACTER(LEN=128) ::                                   &
    FILE_COEF = "rtcoef_noaa_16_hirs.dat",       &
    FILE_SCAER = "",                              &
    FILE_SCCLD = "",                              &
    FILE_PCCOEF = "",&
    FILE_PCCOMP = "",&
    FILE_SURFACE ="",&
    FILE_PROF = "profile.rtp",                &
    FILE_RADR = "radiance.rtr",                  &
    FILE_KMAT = "kmatrix.rtk",                   &
    FILE_TRNS = "transmission.rtt",             &
    PATH_ATLAS= "./"
  CHARACTER(LEN=8) CHNPCSCORES
  
  INTEGER  :: Riargc
  INTEGER(KIND=JPIM) :: NTHREADS = 1_JPIM
  INTEGER(KIND=JPIM) :: NCHANNELS
  INTEGER(KIND=JPIM) :: NPCSCORES
TRY
    
    ERR = 0
   Riargc= rttov_iargc()

    IF( riargc .LE. 0 ) THEN
          PRINT *, "Usage:   mode DIRECT K or GET_EMISBRDF"
          PRINT *, "Usage: DIRECT       FILE_PROF FILE_RADR FILE_TRNS FILE_COEF FILE_SCAER FILE_SCCLD"
          PRINT *, "Usage: K            FILE_PROF FILE_KMAT FILE_COEF FILE_SCAER FILE_SCCLD"     
          PRINT *, "Usage: GET_EMISBRDF FILE_PROF FILE_ISURFACE   FILE_OSURFACE FILE_COEF PATH_ATLAS"
          PRINT *, "Usage: PCDIRECT     FILE_PROF FILE_COEF FILE_PCCOEF FILE_PCCOMP FILE_SURFACE NPCSCORES"
          PRINT *, "Usage: PCK          FILE_PROF FILE_KMAT FILE_COEF FILE_PCCOEF NPCSCORES"

      STOP
    ENDIF
    CALL rttov_getarg( 1, MODE )

    IF( riargc .GE. 2 ) CALL rttov_getarg( 2, FILE_PROF )
    SELECT CASE( MODE )
      CASE( "DIRECT" ) 
        IF( riargc .GE. 3 ) CALL rttov_getarg( 3, FILE_RADR )
        IF( riargc .GE. 4 ) CALL rttov_getarg( 4, FILE_TRNS )
        IF( riargc .GE. 5 ) CALL rttov_getarg( 5, FILE_COEF )
        IF( riargc .GE. 6 ) CALL rttov_getarg( 6, FILE_SCAER )
        IF( riargc .GE. 7 ) CALL rttov_getarg( 7, FILE_SCCLD )
      CASE( "PCDIRECT" ) 
        IF( riargc .GE. 3 ) CALL rttov_getarg( 3, FILE_COEF )
        IF( riargc .GE. 4 ) CALL rttov_getarg( 4, FILE_PCCOEF )
        IF( riargc .GE. 5 ) CALL rttov_getarg( 5, FILE_PCCOMP )
        IF( riargc .GE. 6 ) CALL rttov_getarg( 6, FILE_SURFACE )
        IF( riargc .GE. 7 ) CALL rttov_getarg( 7, CHNPCSCORES )
        write(*,*) "CHNPCSCORES ",CHNPCSCORES
        read(CHNPCSCORES,*) NPCSCORES
        write(*,*) "NPCSCORES ",NPCSCORES

      CASE( "K" )
        IF( riargc .GE. 3 ) CALL rttov_getarg( 3, FILE_KMAT )
        IF( riargc .GE. 4 ) CALL rttov_getarg( 4, FILE_COEF )
        IF( riargc .GE. 5 ) CALL rttov_getarg( 5, FILE_SCAER )
        IF( riargc .GE. 6 ) CALL rttov_getarg( 6, FILE_SCCLD )
      CASE( "PCK" )
        IF( riargc .GE. 3 ) CALL rttov_getarg( 3, FILE_KMAT )
        IF( riargc .GE. 4 ) CALL rttov_getarg( 4, FILE_COEF )
        IF( riargc .GE. 5 ) CALL rttov_getarg( 5, FILE_PCCOEF )
        IF( riargc .GE. 6 ) CALL rttov_getarg( 6, CHNPCSCORES )
        write(*,*) "CHNPCSCORES ",CHNPCSCORES
        read(CHNPCSCORES,*) NPCSCORES
        write(*,*) "NPCSCORES ",NPCSCORES
        
        
      CASE( "GET_EMISBRDF" )
        IF( riargc .GE. 3 ) CALL rttov_getarg( 3, FILE_SURFACE )
        IF( riargc .GE. 4 ) CALL rttov_getarg( 5, FILE_COEF )
        IF( riargc .GE. 6 ) CALL rttov_getarg( 6, PATH_ATLAS )
      CASE DEFAULT
    END SELECT
    
    CALL RTTOV_GUI_LOAD( (/0/), file_coef, file_scaer, file_sccld, file_pccoef, NCHANNELS, ERR )
    THROW(ERR.NE.0)

    WRITE(*,*) "Number of channels ",NCHANNELS


    SELECT CASE( MODE )
      CASE( "DIRECT" , "K" ) 
        CALL RTTOV_GUI_RUN( FILE_PROF, FILE_SURFACE,  &
                     FILE_RADR, FILE_KMAT,     &
                     FILE_TRNS, MODE, NTHREADS, ERR )
      CASE( "PCDIRECT", "PCK" ) 
        CALL RTTOV_GUI_PCRUN( FILE_PROF, FILE_SURFACE,  &
                     FILE_PCCOMP, FILE_KMAT,          &
                     MODE, NTHREADS, NPCSCORES, ERR )

      CASE ( "GET_EMISBRDF" )
        CALL RTTOV_GUI_GET_EMISBRDF( PATH_ATLAS, FILE_PROF, FILE_SURFACE, ERR )
        
      CASE DEFAULT
    END SELECT
    THROW(ERR.NE.0)


    CALL RTTOV_GUI_DROP( ERR )
    THROW(ERR.NE.0)

PCATCH

END PROGRAM
