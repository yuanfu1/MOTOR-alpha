! Description:
!> @file
!!   Generic subroutine to read RTTOV structures or various kinds of scalar or
!!   vector data from an HDF5 file.
!
!> @brief
!!   Generic subroutine to read RTTOV structures or various kinds of scalar or
!!   vector data from an HDF5 file.
!!
!! @param[out]    err           return status
!! @param[in]     f             name of input HDF5 file
!! @param[in]     path          path of group within HDF5 file containing dataset to read
!! @param[out]    profile       single RTTOV profile, optional
!! @param         profiles      pointer to array of RTTOV profiles, optional
!! @param[out]    options       RTTOV options structure, optional
!! @param         chanprof      pointer to RTTOV chanprof structure, optional
!! @param         emissivity    pointer to RTTOV surface emissivity structure, optional
!! @param         reflectance   pointer to RTTOV surface reflectance structure, optional
!! @param         coef          RTTOV optical depth coefficient structure, optional
!! @param         sccldcoef     RTTOV cloud coefficient structure, optional
!! @param         scaercoef     RTTOV aerosol coefficient structure, optional
!! @param         mfasislutcld  MFASIS cloud LUT structure, optional
!! @param         mfasislutaer  MFASIS aerosol LUT structure, optional
!! @param         pccoef        RTTOV PC coefficients structure, optional
!! @param         opt_param     RTTOV aerosol/cloud explicit optical property structure, optional
!! @param         i0            integer scalar (kind jpim), optional
!! @param         r0            real scalar (kind jprb), optional
!! @param         rm0           real scalar (kind jprm), optional
!! @param         c0            character string, optional
!! @param         pl1           pointer to logical 1D array, optional
!! @param         pr1           pointer to real 1D array (kind jprb), optional
!! @param         prm1          pointer to real 1D array (kind jprm), optional
!! @param         pr2           pointer to real 2D array (kind jprb), optional
!! @param         prm2          pointer to real 2D array (kind jprm), optional
!! @param         pr3           pointer to real 3D array (kind jprb), optional
!! @param         prm3          pointer to real 3D array (kind jprm), optional
!! @param         pr4           pointer to real 4D array (kind jprb), optional
!! @param         prm4          pointer to real 4D array (kind jprm), optional
!! @param         pit1          pointer to integer 1D array (kind jpit), optional
!! @param         pit2          pointer to integer 2D array (kind jpit), optional
!! @param         pit3          pointer to integer 3D array (kind jpit), optional
!! @param         pis1          pointer to integer 1D array (kind jpis), optional
!! @param         pis2          pointer to integer 2D array (kind jpis), optional
!! @param         pis3          pointer to integer 3D array (kind jpis), optional
!! @param         pis4          pointer to integer 4D array (kind jpis), optional
!! @param         pi1           pointer to integer 1D array (kind jpim), optional
!! @param         pi2           pointer to integer 2D array (kind jpim), optional
!! @param         pi3           pointer to integer 3D array (kind jpim), optional
!! @param         pil1          pointer to integer 1D array (kind jpib), optional
!! @param         sname         dataset name for scalar/array arguments, optional
!! @param         lbl           set this to true only if reading optical depth coefs 
!!                                from LBL code (default false), optional
!!
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
      SUBROUTINE RTTOV_HDF_LOAD( ERR, F, PATH, PROFILE, PROFILES, OPTIONS, CHANPROF, &
        & EMISSIVITY, REFLECTANCE, COEF, SCCLDCOEF, SCAERCOEF, MFASISLUTCLD, MFASISLUTAER, &
        & PCCOEF, OPT_PARAM, I0, R0, RM0, C0, PL1, PR1, PRM1, PR2, PRM2, PR3, PRM3, PR4, PRM4,&
        & PIT1, PIT2, PIT3, PIS1, PIS2, PIS3, PIS4, PI1, PI2, PI3, PIL1,&
        & SNAME, LBL )
!INTF_OFF
#include "throw.h"
!INTF_ON
      USE PARKIND1
      USE RTTOV_TYPES
!INTF_OFF
      USE RTTOV_HDF_MOD
      USE HDF5
      USE H5LT
! 
      USE RTTOV_HDF_PROFILES
      USE RTTOV_HDF_COEFS
      USE RTTOV_HDF_PROFILE_IO
      USE RTTOV_HDF_S2M_IO
      USE RTTOV_HDF_SKIN_IO
      USE RTTOV_HDF_OPTIONS_IO
      USE RTTOV_HDF_OPTIONS_CONFIG_IO
      USE RTTOV_HDF_OPTIONS_RT_ALL_IO
      USE RTTOV_HDF_OPTIONS_RT_IR_IO
      USE RTTOV_HDF_OPTIONS_RT_MW_IO
      USE RTTOV_HDF_OPTIONS_PC_IO
      USE RTTOV_HDF_OPTIONS_INTERP_IO
      USE RTTOV_HDF_CHANPROF_IO
      USE RTTOV_HDF_EMISSIVITY_IO
      USE RTTOV_HDF_REFLECTANCE_IO
      USE RTTOV_HDF_RTTOV_COEF_IO
      USE RTTOV_HDF_RTTOV_FAST_COEF_IO 
      USE RTTOV_HDF_OPT_PARAM_IO
     
!INTF_ON
!
      IMPLICIT NONE
      CHARACTER(LEN=*),    INTENT(IN)  :: F
      CHARACTER(LEN=*),    INTENT(IN)  :: PATH
      INTEGER(KIND=JPIM),  INTENT(OUT) :: ERR
      TYPE(RTTOV_PROFILE), INTENT(OUT), OPTIONAL :: PROFILE
      TYPE(RTTOV_PROFILE), POINTER, OPTIONAL     :: PROFILES(:)
      TYPE(RTTOV_OPTIONS), INTENT(OUT), OPTIONAL :: OPTIONS
      TYPE(RTTOV_CHANPROF), POINTER, OPTIONAL    :: CHANPROF(:)
      TYPE(RTTOV_EMISSIVITY),  POINTER, OPTIONAL :: EMISSIVITY(:)
      TYPE(RTTOV_REFLECTANCE), POINTER, OPTIONAL :: REFLECTANCE(:)
      TYPE(RTTOV_COEF), OPTIONAL                 :: COEF
      TYPE(RTTOV_COEF_SCATT),    OPTIONAL        :: SCCLDCOEF
      TYPE(RTTOV_COEF_SCATT),    OPTIONAL        :: SCAERCOEF
      TYPE(RTTOV_COEF_MFASIS),   OPTIONAL        :: MFASISLUTCLD
      TYPE(RTTOV_COEF_MFASIS),   OPTIONAL        :: MFASISLUTAER
      TYPE(RTTOV_COEF_PCCOMP),   OPTIONAL        :: PCCOEF
      TYPE(RTTOV_OPT_PARAM),     OPTIONAL        :: OPT_PARAM

      INTEGER(KIND=JPIM),       OPTIONAL :: I0
      REAL(KIND=JPRB),          OPTIONAL :: R0
      REAL(KIND=JPRM),          OPTIONAL :: RM0
      CHARACTER(LEN=*),         OPTIONAL :: C0
      INTEGER(KIND=JPIT),     POINTER , OPTIONAL :: PIT1(:)
      INTEGER(KIND=JPIT),     POINTER , OPTIONAL :: PIT2(:,:)
      INTEGER(KIND=JPIT),     POINTER , OPTIONAL :: PIT3(:,:,:)
      INTEGER(KIND=JPIS),     POINTER , OPTIONAL :: PIS1(:)
      INTEGER(KIND=JPIS),     POINTER , OPTIONAL :: PIS2(:,:)
      INTEGER(KIND=JPIS),     POINTER , OPTIONAL :: PIS3(:,:,:)
      INTEGER(KIND=JPIS),     POINTER , OPTIONAL :: PIS4(:,:,:,:)
      INTEGER(KIND=JPIM),     POINTER , OPTIONAL :: PI1(:)
      INTEGER(KIND=JPIM),     POINTER , OPTIONAL :: PI2(:,:)
      INTEGER(KIND=JPIM),     POINTER , OPTIONAL :: PI3(:,:,:)
      INTEGER(KIND=JPIB),     POINTER , OPTIONAL :: PIL1(:)
      REAL(KIND=JPRB),        POINTER , OPTIONAL :: PR1(:)
      REAL(KIND=JPRB),        POINTER , OPTIONAL :: PR2(:,:)
      REAL(KIND=JPRB),        POINTER , OPTIONAL :: PR3(:,:,:)
      REAL(KIND=JPRB),        POINTER , OPTIONAL :: PR4(:,:,:,:)
      REAL(KIND=JPRM),        POINTER , OPTIONAL :: PRM1(:)
      REAL(KIND=JPRM),        POINTER , OPTIONAL :: PRM2(:,:)
      REAL(KIND=JPRM),        POINTER , OPTIONAL :: PRM3(:,:,:)
      REAL(KIND=JPRM),        POINTER , OPTIONAL :: PRM4(:,:,:,:)
      LOGICAL,                POINTER , OPTIONAL :: PL1(:)
      CHARACTER(LEN=*),       OPTIONAL :: SNAME
      LOGICAL,                OPTIONAL :: LBL
!INTF_END
#include "rttov_errorreport.interface"

      INTEGER(KIND=JPIM) :: I, IPROF, NPROF
      LOGICAL            :: LEXT
      CHARACTER(LEN=4)   :: CHPROF
!
      INTEGER(HID_T) :: FILE_ID, G_ID, G_ID_SUB
      LOGICAL        :: LBL1
!
TRY

      LBL1 = .FALSE.
      IF(PRESENT(LBL)) LBL1 = LBL

      CALL H5FOPEN_F( F, H5F_ACC_RDONLY_F, FILE_ID, ERR )
      THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(F))

      CALL H5GOPEN_F( FILE_ID, PATH, G_ID, ERR )
      THROWM(ERR.NE.0,"CANNOT OPEN GROUP "//TRIM(PATH)//" IN FILE "//TRIM(F))

      IF( PRESENT(PROFILE) ) THEN
      
        CALL RTTOV_HDF_one_profile_RH( PROFILE, G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT READ PROFILE IN FILE "//TRIM(F))

      ENDIF
            
      IF( PRESENT(PROFILES) ) THEN
      
        IPROF = 0
        DO
          IF( IPROF .LT. 9000_JPIM ) THEN
            WRITE( CHPROF, '(I4.4)' ) IPROF+1
          ENDIF
          CALL H5LEXISTS_F( G_ID, TRIM(CHPROF), LEXT, ERR )
          IF( .NOT. LEXT ) EXIT
          IPROF = IPROF+1
        ENDDO
        NPROF = IPROF

        ALLOCATE( PROFILES(NPROF), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE PROFILES")
        
        DO I = 1, NPROF
        
          WRITE(CHPROF,"(I4.4)") I

          CALL H5GOPEN_F( G_ID, CHPROF, G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT OPEN GROUP FOR PROFILE "//CHPROF)
      
          CALL RTTOV_HDF_one_profile_RH( PROFILES(I), G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT READ PROFILE"//CHPROF//" IN FILE "//TRIM(F))
          
          CALL H5GCLOSE_F( G_ID_SUB, ERR )
          THROWM(ERR.NE.0,"CANNOT CLOSE GROUP FOR PROFILE "//CHPROF)

        END DO
        
      ENDIF
      IF( PRESENT(OPTIONS) ) THEN
      
        CALL RTTOV_HDF_OPTIONS_RH( OPTIONS, G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT READ OPTIONS IN FILE "//TRIM(F))

        CALL RTTOV_HDF_OPTIONS_CONFIG_RH( OPTIONS%CONFIG, G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT READ OPTIONS%CONFIG IN FILE "//TRIM(F))

        CALL RTTOV_HDF_OPTIONS_RT_ALL_RH( OPTIONS%RT_ALL, G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT READ OPTIONS%RT_ALL IN FILE "//TRIM(F))

        CALL RTTOV_HDF_OPTIONS_RT_IR_RH( OPTIONS%RT_IR, G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT READ OPTIONS%RT_IR IN FILE "//TRIM(F))

        CALL RTTOV_HDF_OPTIONS_RT_MW_RH( OPTIONS%RT_MW, G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT READ OPTIONS%RT_MW IN FILE "//TRIM(F))

        CALL RTTOV_HDF_OPTIONS_PC_RH( OPTIONS%RT_IR%PC, G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT READ OPTIONS%RT_IR%PC IN FILE "//TRIM(F))

        CALL RTTOV_HDF_OPTIONS_INTERP_RH( OPTIONS%INTERPOLATION, G_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT READ OPTIONS%INTERPOLATION IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(CHANPROF) ) THEN
      
        CALL RTTOV_HDF_CHANPROF_RH(CHANPROF,G_ID,ERR)
        THROWM(ERR.NE.0,"CANNOT READ CHANPROF IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(EMISSIVITY) ) THEN
      
        CALL RTTOV_HDF_EMISSIVITY_RH(EMISSIVITY,G_ID,ERR)
        THROWM(ERR.NE.0,"CANNOT READ EMISSIVITY IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(REFLECTANCE) ) THEN
      
        CALL RTTOV_HDF_REFLECTANCE_RH(REFLECTANCE,G_ID,ERR)
        THROWM(ERR.NE.0,"CANNOT READ REFLECTANCE IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(SCAERCOEF) .AND. PRESENT(COEF)) THEN
      
        CALL RTTOV_HDF_SCAERCOEF_RH(SCAERCOEF, COEF, G_ID,ERR)
        THROWM(ERR.NE.0,"CANNOT READ AEROSOLS COEF IN FILE "//TRIM(F))
           
      ELSEIF( PRESENT(SCCLDCOEF) .AND. PRESENT(COEF)) THEN
      
        CALL RTTOV_HDF_SCCLDCOEF_RH(SCCLDCOEF, COEF, G_ID,ERR)
        THROWM(ERR.NE.0,"CANNOT READ CLOUDS COEF IN FILE "//TRIM(F))
            
      ELSEIF( PRESENT(MFASISLUTCLD)) THEN

        CALL RTTOV_HDF_MFASISCOEF_RH(MFASISLUTCLD, 'CLD', G_ID,ERR)
        THROWM(ERR.NE.0,"CANNOT READ MFASIS CLD LUT IN FILE "//TRIM(F))

      ELSEIF( PRESENT(MFASISLUTAER)) THEN

        CALL RTTOV_HDF_MFASISCOEF_RH(MFASISLUTAER, 'AER', G_ID,ERR)
        THROWM(ERR.NE.0,"CANNOT READ MFASIS AER LUT IN FILE "//TRIM(F))

      ELSEIF( PRESENT(PCCOEF) .AND. PRESENT(COEF)) THEN
      
        CALL RTTOV_HDF_PCCOEF_RH(PCCOEF, COEF, G_ID,ERR)
        THROWM(ERR.NE.0,"CANNOT READ PC COEFS COEF IN FILE "//TRIM(F))

      ELSEIF( PRESENT(COEF) ) THEN

        CALL RTTOV_HDF_COEF_RH(COEF,G_ID,ERR, LBL = LBL1)
        THROWM(ERR.NE.0,"CANNOT READ COEF IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(OPT_PARAM) ) THEN
      
        CALL RTTOV_HDF_OPT_PARAM_RH(OPT_PARAM, G_ID,ERR)
        THROWM(ERR.NE.0,"CANNOT READ OPT_PARAM IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(PR1) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,pr1=pr1)
        THROWM(ERR.NE.0,"CANNOT READ REAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(PRM1) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,prm1=prm1)
        THROWM(ERR.NE.0,"CANNOT READ REAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF


      IF( PRESENT(PR2) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,pr2=pr2)
        THROWM(ERR.NE.0,"CANNOT READ REAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(PRM2) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,prm2=prm2)
        THROWM(ERR.NE.0,"CANNOT READ REAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(PR3) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,pr3=pr3)
        THROWM(ERR.NE.0,"CANNOT READ REAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(PRM3) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,prm3=prm3)
        THROWM(ERR.NE.0,"CANNOT READ REAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(PR4) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,pr4=pr4)
        THROWM(ERR.NE.0,"CANNOT READ REAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(PRM4) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,prm4=prm4)
        THROWM(ERR.NE.0,"CANNOT READ REAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(PIT1) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,pit1=pit1)
        THROWM(ERR.NE.0,"CANNOT READ INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF 
      IF( PRESENT(PIT2) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,pit2=pit2)
        THROWM(ERR.NE.0,"CANNOT READ INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF 
      IF( PRESENT(PIT3) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,pit3=pit3)
        THROWM(ERR.NE.0,"CANNOT READ INTEGER  ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(PIS1) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,pis1=pis1)
        THROWM(ERR.NE.0,"CANNOT READ INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF 
      IF( PRESENT(PIS2) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,pis2=pis2)
        THROWM(ERR.NE.0,"CANNOT READ INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF 
      IF( PRESENT(PIS3) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,pis3=pis3)
        THROWM(ERR.NE.0,"CANNOT READ INTEGER  ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF 
      IF( PRESENT(PIS4) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,pis4=pis4)
        THROWM(ERR.NE.0,"CANNOT READ INTEGER  ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF 

      IF( PRESENT(PI1) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,pi1=pi1)
        THROWM(ERR.NE.0,"CANNOT READ INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF 
      IF( PRESENT(PI2) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,pi2=pi2)
        THROWM(ERR.NE.0,"CANNOT READ INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF 
      IF( PRESENT(PI3) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,pi3=pi3)
        THROWM(ERR.NE.0,"CANNOT READ INTEGER  ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(PIL1) .AND. PRESENT(SNAME) ) THEN

        call read_array_hdf(G_ID,SNAME,err,pil1=pil1)
        THROWM(ERR.NE.0,"CANNOT READ INTEGER ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF 

      IF( PRESENT(PL1) .AND. PRESENT(SNAME) ) THEN
      
        call read_array_hdf(G_ID,SNAME,err,pl1=pl1)
        THROWM(ERR.NE.0,"CANNOT READ LOGICAL ARRAY "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(C0) .AND. PRESENT(SNAME) ) THEN
      
        call read_array_hdf(G_ID,SNAME,err,c0=c0)
        THROWM(ERR.NE.0,"CANNOT READ STRING "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(I0) .AND. PRESENT(SNAME) ) THEN
                          
        call read_array_hdf(G_ID,SNAME,err,i0=i0)
        THROWM(ERR.NE.0,"CANNOT READ INTEGER "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF
      
      IF( PRESENT(R0) .AND. PRESENT(SNAME) ) THEN
                          
        call read_array_hdf(G_ID,SNAME,err,r0=r0)
        THROWM(ERR.NE.0,"CANNOT READ REAL SCALAR "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      IF( PRESENT(RM0) .AND. PRESENT(SNAME) ) THEN
                          
        call read_array_hdf(G_ID,SNAME,err,rm0=rm0)
        THROWM(ERR.NE.0,"CANNOT READ REAL SCALAR "//TRIM(SNAME)//" IN FILE "//TRIM(F))

      ENDIF

      CALL H5GCLOSE_F( G_ID, ERR )
      THROWM(ERR.NE.0,"CANNOT CLOSE GROUP "//TRIM(PATH)//" IN FILE "//TRIM(F))
          
      
      CALL H5FCLOSE_F( FILE_ID, ERR )
      THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(F))

CATCH
      END SUBROUTINE
