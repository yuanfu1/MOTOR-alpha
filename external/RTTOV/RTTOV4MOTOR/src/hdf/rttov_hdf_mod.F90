! Description:
!> @file
!!   Subroutines for HDF5 input/output
!
!> @brief
!!   Subroutines for HDF5 input/output
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
MODULE RTTOV_HDF_MOD
!
!     History:
!
!     Version   Date      Author   Comment
!     -------   ----      ------   -------
!
!     r1241     12/2012   DAR      added 4D real support for nlte_coefs
!     r1517+    04/2014   PB       added 8bits and 16bits integers for EMIS/BRDF Atlases
!     r1543+    05/2014   PB       added FORCE_DOUBLE in write_array_hdf for real variables
!     r1799     09/2015   PB       added JPRM reals, new macro allowing the use of C pointers 
!                                    in read routines if HDF5 library compiled with 
!                                    option --enable-fortran2003  (this reduce amount of memory)
!                                    macro name is _RTTOV_HDF_F2003
!
!  Note: H5DREAD_F normal return code is any positive value
!     In H5Dio.c which is the C routine that reads the data and returns the error code it is written:
!      * Function: H5D_read
!      *
!      * Purpose:  Reads (part of) a DATASET into application memory BUF. See
!      *           H5Dread() for complete details.
!      *
!      * Return:   Non-negative on success/Negative on failure
!      So failure read execution should be tested by ERR.LT.0
!
#include "throw.h"

USE PARKIND1, ONLY : JPIT, JPIS, JPIM, JPIB, JPRB, JPRM, JPLM
USE RTTOV_CONST, ONLY : ERRORSTATUS_FATAL
USE HDF5
USE H5LT

#ifdef _RTTOV_HDF_F2003
USE ISO_C_BINDING
#endif

IMPLICIT NONE

#include "rttov_errorreport.interface"

CHARACTER(LEN=*), PARAMETER :: H5_EXT = ".H5"
INTEGER :: LENSH = 64 ! length for dataset names
INTEGER(HID_T) :: H5T_JPRB_REAL, H5T_FILE_REAL, H5T_JPIM_INTEGER, H5T_FILE_INTEGER, H5T_FILE_CHARACTER
INTEGER(HID_T) :: H5T_JPRM_REAL, H5T_STD_I8, H5T_STD_I16, H5T_STD_I64

LOGICAL(JPLM) :: IS_HDF_OPEN = .FALSE.
LOGICAL(JPLM) :: IS_HDF_64BIT_REALS = .TRUE.

Public :: LENSH
PUBLIC :: H5_EXT, H5T_JPRB_REAL, H5T_FILE_REAL, H5T_JPIM_INTEGER, H5T_FILE_INTEGER, H5T_FILE_CHARACTER
PUBLIC :: H5T_STD_I8, H5T_STD_I16, H5T_STD_I64

PUBLIC :: OPEN_HDF
PUBLIC :: CLOSE_HDF
PUBLIC :: IS_HDF_OPEN, IS_HDF_64BIT_REALS
PUBLIC :: WRITE_ARRAY_HDF, WRITE_ARRAY_HDF_CMPLX
PUBLIC :: READ_ARRAY_HDF, READ_ARRAY_HDF_CMPLX
PUBLIC :: EXGRP, MKDSR, MKPAR
PUBLIC :: RTTOV_HDF_OPENW, RTTOV_HDF_OPENR, RTTOV_HDF_OPEN_TRUNC, RTTOV_HDF_OPEN_RDWR
PUBLIC :: RTTOV_HDF_LISTOBJ
PUBLIC :: RTTOV_HDF_INIT_GZIP

! Minimum array chunk size for gzip compression (default value)
INTEGER(KIND=JPIM), PARAMETER :: MINSAVALP = 4096_JPIM

LOGICAL(JPLM) :: ENABLE_SHUFFLE = .TRUE.
PUBLIC :: ENABLE_SHUFFLE

PRIVATE

CONTAINS

!> @brief Open HDF5 library
!! @param[in]  DBLE  if true reals are stored as H5T_NATIVE_DOUBLE, otherwise as H5T_NATIVE_REAL
!! @param[out] ERR   status on exit
      SUBROUTINE OPEN_HDF( DBLE, ERR )
!
        LOGICAL,            INTENT(IN)  :: DBLE
        INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
!
        INTEGER :: HDFERR
        LOGICAL, PARAMETER :: BIGEND = ICHAR(TRANSFER(1,'A')) == 0

TRY
        CALL H5OPEN_F( HDFERR )
        ERR = HDFERR
        THROW(ERR.NE.0)

        H5T_JPRM_REAL      = H5T_NATIVE_REAL

        SELECT CASE( JPRB )
          CASE( 4 )
            H5T_JPRB_REAL    = H5T_NATIVE_REAL
          CASE( 8 )
            H5T_JPRB_REAL    = H5T_NATIVE_DOUBLE
          CASE DEFAULT
            ERR = 1
            RETURN
        END SELECT

        IF( DBLE ) THEN
          H5T_FILE_REAL    = H5T_NATIVE_DOUBLE
        ELSE
          H5T_FILE_REAL    = H5T_NATIVE_REAL
        ENDIF

        H5T_JPIM_INTEGER = H5T_NATIVE_INTEGER
        H5T_FILE_INTEGER = H5T_NATIVE_INTEGER
        !H5T_FILE_CHARACTER = H5T_NATIVE_CHARACTER
        H5T_FILE_CHARACTER = H5T_STRING

        IF( BIGEND ) THEN
          H5T_STD_I8  = H5T_STD_I8BE
          H5T_STD_I16 = H5T_STD_I16BE
          H5T_STD_I64 = H5T_STD_I64BE
        ELSE
          H5T_STD_I8  = H5T_STD_I8LE
          H5T_STD_I16 = H5T_STD_I16LE
          H5T_STD_I64 = H5T_STD_I64LE
        ENDIF

        IS_HDF_OPEN = .TRUE.
        IS_HDF_64BIT_REALS = DBLE

CATCH     
      END SUBROUTINE
      
!> @brief Close HDF library
!! @param[out] ERR   status on exit
      SUBROUTINE CLOSE_HDF( ERR )
!
        INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
!
        INTEGER :: HDFERR
TRY
        IF (IS_HDF_OPEN) THEN
          CALL H5CLOSE_F( HDFERR )
          ERR = HDFERR
          THROW(ERR.NE.0)
          IS_HDF_OPEN = .FALSE.
        ENDIF
CATCH
      END SUBROUTINE

!> @brief Write an array in an opened HDF5 file
!! @param[in]  LUN  HDF5 identification
!! @param[in]  NAME datatset name
!! @param[in]  COMMENT comment for dataset
!! @param[out] ERR status on exit
!! @param[in]  R0 scalar of type real JPRB, optional
!! @param[in]  R1 array of type real JPRB, dimension 1, optional
!! @param[in]  R2 array of type real JPRB, dimension 2, optional
!! @param[in]  R3 array of type real JPRB, dimension 3, optional
!! @param[in]  R4 array of type real JPRB, dimension 4, optional
!! @param[in]  RM0 scalar of type real JPRM, optional
!! @param[in]  RM1 array of type real JPRM, dimension 1, optional
!! @param[in]  RM2 array of type real JPRM, dimension 2, optional
!! @param[in]  RM3 array of type real JPRM, dimension 3, optional
!! @param[in]  RM4 array of type real JPRM, dimension 4, optional
!! @param[in]  I0 scalar of type integer JPIM, optional
!! @param[in]  I1 array of type integer JPIM, dimension 1, optional
!! @param[in]  I2 array of type integer JPIM, dimension 2, optional
!! @param[in]  I3 array of type integer JPIM, dimension 3, optional
!! @param[in]  IT1 array of type integer JPIT, dimension 1, optional
!! @param[in]  IT2 array of type integer JPIT, dimension 2, optional
!! @param[in]  IT3 array of type integer JPIT, dimension 3, optional
!! @param[in]  IS1 array of type integer JPIS, dimension 1, optional
!! @param[in]  IS2 array of type integer JPIS, dimension 2, optional
!! @param[in]  IS3 array of type integer JPIS, dimension 3, optional
!! @param[in]  IL1 array of type integer JPIB, dimension 1, optional
!! @param[in]  C0 string, optional
!! @param[in]  C1 array of strings, dimension 1, optional
!! @param[in]  L0 logical, optional
!! @param[in]  L1 array of logicals, dimension 1, optional
!! @param[in]  UNITS dataset units attribute, optional
!! @param[in]  COMPRESS internal compression, optional
!! @param[in]  COMPDIMS dimensions to be addressed by compression, optional
!!   is used for definition of compression chunks, optional
!! @param[in]  FORCE_DOUBLE store real as H5T_NATIVE_DOUBLE, optional
      SUBROUTINE WRITE_ARRAY_HDF( LUN, NAME, COMMENT, ERR, I0, I1, &
        I2, I3, IT1, IT2, IT3, IS1, IS2, IS3, IL1, C0, C1, L0, L1, &
        R0, R1, R2, R3, R4, &
        RM0, RM1, RM2, RM3, RM4, &
        UNITS, COMPRESS, COMPDIMS, FORCE_DOUBLE )
!
      INTEGER(HID_T), INTENT(IN) :: LUN
      CHARACTER(LEN=*),  INTENT(IN) :: NAME
      CHARACTER(LEN=*),  INTENT(IN) :: COMMENT
      INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
      REAL(KIND=JPRB), INTENT(IN), OPTIONAL :: R0
      REAL(KIND=JPRB), INTENT(IN), OPTIONAL :: R1(:)
      REAL(KIND=JPRB), INTENT(IN), OPTIONAL :: R2(:,:)
      REAL(KIND=JPRB), INTENT(IN), OPTIONAL :: R3(:,:,:)
      REAL(KIND=JPRB), INTENT(IN), OPTIONAL :: R4(:,:,:,:)
      REAL(KIND=JPRM), INTENT(IN), OPTIONAL :: RM0
      REAL(KIND=JPRM), INTENT(IN), OPTIONAL :: RM1(:)
      REAL(KIND=JPRM), INTENT(IN), OPTIONAL :: RM2(:,:)
      REAL(KIND=JPRM), INTENT(IN), OPTIONAL :: RM3(:,:,:)
      REAL(KIND=JPRM), INTENT(IN), OPTIONAL :: RM4(:,:,:,:)

      INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: I0
      INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: I1(:)
      INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: I2(:,:)
      INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: I3(:,:,:)
      INTEGER(KIND=JPIT), INTENT(IN), OPTIONAL :: IT1(:)
      INTEGER(KIND=JPIT), INTENT(IN), OPTIONAL :: IT2(:,:)
      INTEGER(KIND=JPIT), INTENT(IN), OPTIONAL :: IT3(:,:,:)
      INTEGER(KIND=JPIS), INTENT(IN), OPTIONAL :: IS1(:)
      INTEGER(KIND=JPIS), INTENT(IN), OPTIONAL :: IS2(:,:)
      INTEGER(KIND=JPIS), INTENT(IN), OPTIONAL :: IS3(:,:,:)
      INTEGER(KIND=JPIB), INTENT(IN), OPTIONAL :: IL1(:)

      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: C0
      CHARACTER(LEN=*), INTENT(IN), target, OPTIONAL :: C1(:)
!       CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: C2(:,:)
!       CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: C3(:,:,:)

      LOGICAL, INTENT(IN), OPTIONAL :: L0
      LOGICAL, INTENT(IN), OPTIONAL :: L1(:)
      CHARACTER(LEN=*),  INTENT(IN) , OPTIONAL :: UNITS
      LOGICAL, INTENT(IN), OPTIONAL :: COMPRESS
      INTEGER(KIND=JPIM), INTENT(IN), OPTIONAL :: COMPDIMS(:)
      LOGICAL, INTENT(IN), OPTIONAL :: FORCE_DOUBLE

      INTEGER(HSIZE_T) :: DIMS(7), CHUNK_DIMS(7)
      INTEGER(SIZE_T) ::  RANK8
      INTEGER(KIND=JPIM) :: RANK, CHUNK_RANK
      INTEGER(HID_T) :: TYPE_F_ID, SPCE_ID, DSET_ID, CRP_LIST, FILE_REAL

      INTEGER(KIND=JPIM) :: UB(7), LB(7)
      INTEGER(KIND=JPIM) :: I
      INTEGER(KIND=JPIM) :: IL0
      INTEGER(KIND=JPIM), ALLOCATABLE:: LL1(:)

      INTEGER(KIND=JPIM) :: WI4

      LOGICAL :: SCAL, LEXT

      INTEGER(SIZE_T) :: SDIM

      LOGICAL :: LCOMPRESS 

      INTEGER(KIND=JPIM) :: CLEVEL
      !INTEGER(KIND=JPIM) :: SA

!
TRY
      CLEVEL = 6 ! COMPRESSION LEVEL
      DIMS = 0
      LB = 1
      UB = -1
      SCAL = .FALSE.
      TYPE_F_ID = -1

      ! USE FORCE_DOUBLE IN ORDER TO OVERIDE THE VALUE OF H5T_FILE_REAL
      ! THIS IS USEFULL IF USER WANTS TO STORE 64b DATA IN 32BITS
      ! AND VICE-VERSA
      IF( PRESENT( FORCE_DOUBLE ) ) THEN
        IF( FORCE_DOUBLE ) THEN
          FILE_REAL    = H5T_NATIVE_DOUBLE
        ELSE
          FILE_REAL    = H5T_NATIVE_REAL
        ENDIF
      ELSE
        FILE_REAL = H5T_FILE_REAL
      ENDIF

      IF( PRESENT( C0 ) ) THEN
        SCAL = .TRUE.
      ELSE IF( PRESENT( L0 ) ) THEN
        LB(1) = 1
        UB(1) = 1
        TYPE_F_ID = H5T_FILE_INTEGER
        SCAL = .TRUE.
      ELSE IF( PRESENT( R0 ) ) THEN
        LB(1) = 1
        UB(1) = 1
        TYPE_F_ID = FILE_REAL
        SCAL = .TRUE.
      ELSE IF( PRESENT( RM0 ) ) THEN
        LB(1) = 1
        UB(1) = 1
        TYPE_F_ID = FILE_REAL
        SCAL = .TRUE.
      ELSE IF( PRESENT( I0 ) ) THEN
        LB(1) = 1
        UB(1) = 1
        TYPE_F_ID = H5T_FILE_INTEGER
        SCAL = .TRUE.
      ELSE IF( PRESENT( C1 ) ) THEN
        LB(1:1) = LBOUND( C1 )
        UB(1:1) = UBOUND( C1 )
        SDIM = LEN( C1(1) )
        CALL H5Tcopy_f(H5T_FORTRAN_S1, TYPE_F_ID, ERR)
        CALL H5Tset_size_f(TYPE_F_ID, sdim, ERR)

      ELSE IF( PRESENT( R1 ) ) THEN
        LB(1:1) = LBOUND( R1 )
        UB(1:1) = UBOUND( R1 )
        TYPE_F_ID = FILE_REAL
      ELSE IF( PRESENT( RM1 ) ) THEN
        LB(1:1) = LBOUND( RM1 )
        UB(1:1) = UBOUND( RM1 )
        TYPE_F_ID = FILE_REAL
      ELSE IF( PRESENT( I1 ) ) THEN
        LB(1:1) = LBOUND( I1 )
        UB(1:1) = UBOUND( I1 )
        TYPE_F_ID = H5T_FILE_INTEGER
      ELSE IF( PRESENT( L1 ) ) THEN
        LB(1:1) = LBOUND( L1 )
        UB(1:1) = UBOUND( L1 )
        TYPE_F_ID = H5T_FILE_INTEGER
        DIMS(1) = SIZE( L1, 1 )
        ALLOCATE( LL1(DIMS(1)), STAT=ERR)
        THROWM(ERR.NE.0,"COULD NOT ALLOCATE LL1 ")
        LL1(:) = 0_JPIM
        WHERE ( L1 ) 
          LL1 = 1_JPIM
        END WHERE
      ELSE IF( PRESENT( R2 ) ) THEN
        LB(1:2) = LBOUND( R2 )
        UB(1:2) = UBOUND( R2 )
        TYPE_F_ID = FILE_REAL
      ELSE IF( PRESENT( RM2 ) ) THEN
        LB(1:2) = LBOUND( RM2 )
        UB(1:2) = UBOUND( RM2 )
        TYPE_F_ID = FILE_REAL
      ELSE IF( PRESENT( I2 ) ) THEN
        LB(1:2) = LBOUND( I2 )
        UB(1:2) = UBOUND( I2 )
        TYPE_F_ID = H5T_FILE_INTEGER
      ELSE IF( PRESENT( R3 ) ) THEN
        LB(1:3) = LBOUND( R3 )
        UB(1:3) = UBOUND( R3 )
        TYPE_F_ID = FILE_REAL
      ELSE IF( PRESENT( RM3 ) ) THEN
        LB(1:3) = LBOUND( RM3 )
        UB(1:3) = UBOUND( RM3 )
        TYPE_F_ID = FILE_REAL
      ELSE IF( PRESENT( I3 ) ) THEN
        LB(1:3) = LBOUND( I3 )
        UB(1:3) = UBOUND( I3 )
        TYPE_F_ID = H5T_FILE_INTEGER
      ELSE IF( PRESENT( R4 ) ) THEN
        LB(1:4) = LBOUND( R4 )
        UB(1:4) = UBOUND( R4 )
        TYPE_F_ID = FILE_REAL
       ELSE IF( PRESENT( RM4 ) ) THEN
        LB(1:4) = LBOUND( RM4 )
        UB(1:4) = UBOUND( RM4 )
        TYPE_F_ID = FILE_REAL
     ELSE IF( PRESENT( IT1 ) ) THEN
        LB(1:1) = LBOUND( IT1 )
        UB(1:1) = UBOUND( IT1 )
        TYPE_F_ID = H5T_STD_I8
      ELSE IF( PRESENT( IT2 ) ) THEN
        LB(1:2) = LBOUND( IT2 )
        UB(1:2) = UBOUND( IT2 )
        TYPE_F_ID = H5T_STD_I8
      ELSE IF( PRESENT( IT3 ) ) THEN
        LB(1:3) = LBOUND( IT3 )
        UB(1:3) = UBOUND( IT3 )
        TYPE_F_ID = H5T_STD_I8
      ELSE IF( PRESENT( IS1 ) ) THEN
        LB(1:1) = LBOUND( IS1 )
        UB(1:1) = UBOUND( IS1 )
        TYPE_F_ID = H5T_STD_I16
      ELSE IF( PRESENT( IS2 ) ) THEN
        LB(1:2) = LBOUND( IS2 )
        UB(1:2) = UBOUND( IS2 )
        TYPE_F_ID = H5T_STD_I16
      ELSE IF( PRESENT( IS3 ) ) THEN
        LB(1:3) = LBOUND( IS3 )
        UB(1:3) = UBOUND( IS3 )
        TYPE_F_ID = H5T_STD_I16
      ELSE IF( PRESENT( IL1 ) ) THEN
        LB(1:1) = LBOUND( IL1 )
        UB(1:1) = UBOUND( IL1 )
        TYPE_F_ID = H5T_STD_I64

      ELSE
        ERR = 1
      ENDIF
      THROWM(ERR.NE.0,"COULD NOT CREATE "//TRIM(NAME))

      DIMS = UB - LB + 1
      DO I = 7, 1, -1
        IF( DIMS(I) .GT. 0 ) EXIT
      ENDDO
      RANK = I


      IF( PRESENT(COMPRESS) ) THEN
         LCOMPRESS = COMPRESS
         CALL RTTOV_HDF_INIT_GZIP( LCOMPRESS, DIMS, CHUNK_DIMS, &
          &      CHUNK_RANK, ERR, COMPDIMS=COMPDIMS)
         THROW(ERR.NE.0)
      ELSE
         LCOMPRESS = .FALSE.
      ENDIF   


      IF( TYPE_F_ID .GE. 0 ) THEN
        CALL H5LEXISTS_F( LUN, NAME, LEXT, ERR )
        THROWM(ERR.NE.0,"COULD NOT TEST EXISTANCE "//TRIM(NAME))
        IF( LEXT ) THEN
          CALL H5LDELETE_F( LUN, NAME, ERR )
          THROWM(ERR.NE.0,"COULD NOT DELETE "//TRIM(NAME))
        ENDIF
        CALL H5SCREATE_SIMPLE_F(RANK, DIMS, SPCE_ID, ERR)
        THROWM(ERR.NE.0,"COULD NOT CREATE "//TRIM(NAME))


        IF( LCOMPRESS ) THEN
          CALL h5pcreate_f(H5P_DATASET_CREATE_F, crp_list, err)
          THROWM(ERR.NE.0,"COULD NOT CREATE PROPERTY LIST FOR "//TRIM(NAME))
          CALL h5pset_chunk_f(crp_list, chunk_rank, chunk_dims, err)
          THROWM(ERR.NE.0,"COULD NOT CREATE CHUNK FOR "//TRIM(NAME))
          IF( ENABLE_SHUFFLE ) THEN
            CALL h5pset_shuffle_f(crp_list, err)
            THROWM(ERR.NE.0,"COULD NOT SET SHUFFLE FOR "//TRIM(NAME))
          ENDIF
          CALL h5pset_deflate_f(crp_list, clevel, err)
          THROWM(ERR.NE.0,"COULD NOT SET DEFLATE FOR "//TRIM(NAME))
          CALL H5DCREATE_F(LUN, NAME, TYPE_F_ID, SPCE_ID, DSET_ID, ERR, CRP_LIST)
          THROWM(ERR.NE.0,"COULD NOT CREATE "//TRIM(NAME))
          CALL H5LTSET_ATTRIBUTE_STRING_F( LUN, NAME, "COMPRESSED", "DEFLATE FILTER"//CHAR(0), ERR)
          THROWM(ERR.NE.0,"COULD NOT WRITE COMPRESSED ATTRIBUTE IN "//TRIM(NAME))
        ELSE
          CALL H5DCREATE_F(LUN, NAME, TYPE_F_ID, SPCE_ID, DSET_ID, ERR)
          THROWM(ERR.NE.0,"COULD NOT CREATE "//TRIM(NAME))
        ENDIF

        CALL H5SCLOSE_F(SPCE_ID, ERR)
        THROWM(ERR.NE.0,"COULD NOT CLOSE "//TRIM(NAME))
        
        IF( LCOMPRESS ) THEN
          CALL H5PCLOSE_F(CRP_LIST, ERR)
          THROWM(ERR.NE.0,"COULD NOT CLOSE PROPERTY LIST FOR "//TRIM(NAME))
        ENDIF
      ENDIF

      IF( PRESENT( C0 ) ) THEN
        CALL H5LEXISTS_F( LUN, NAME, LEXT, ERR )
        THROWM(ERR.NE.0,"COULD NOT TEST EXISTANCE "//TRIM(NAME))
        IF( LEXT ) THEN
          CALL H5LDELETE_F( LUN, NAME, ERR )
          THROWM(ERR.NE.0,"COULD NOT DELETE "//TRIM(NAME))
        ENDIF
        CALL H5LTMAKE_DATASET_STRING_F( LUN, NAME, TRIM(C0)//CHAR(0), ERR)
      ELSE IF( PRESENT( R0 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_JPRB_REAL,    R0, DIMS, ERR)
      ELSE IF( PRESENT( RM0 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_JPRM_REAL,    RM0, DIMS, ERR)
      ELSE IF( PRESENT( I0 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_JPIM_INTEGER, I0, DIMS, ERR)
      ELSE IF( PRESENT( L0 ) ) THEN
        IF( L0 ) THEN 
          IL0 = 1_JPIM
        ELSE
          IL0 = 0_JPIM
        ENDIF
        CALL H5DWRITE_F(DSET_ID, H5T_JPIM_INTEGER, IL0, DIMS, ERR)
      ELSE IF( PRESENT( C1 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, TYPE_F_ID, C1(1)(1:1), DIMS, err)
      ELSE IF( PRESENT( R1 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_JPRB_REAL,    R1, DIMS, ERR)
      ELSE IF( PRESENT( RM1 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_JPRM_REAL,    RM1, DIMS, ERR)
      ELSE IF( PRESENT( I1 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_JPIM_INTEGER, I1, DIMS, ERR)
      ELSE IF( PRESENT( L1 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_JPIM_INTEGER, LL1, DIMS, ERR)
      ELSE IF( PRESENT( R2 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_JPRB_REAL,    R2, DIMS, ERR)
      ELSE IF( PRESENT( RM2 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_JPRM_REAL,    RM2, DIMS, ERR)
      ELSE IF( PRESENT( I2 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_JPIM_INTEGER, I2, DIMS, ERR)
      ELSE IF( PRESENT( R3 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_JPRB_REAL,    R3, DIMS, ERR)
      ELSE IF( PRESENT( RM3 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_JPRM_REAL,    RM3, DIMS, ERR)
      ELSE IF( PRESENT( I3 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_JPIM_INTEGER, I3, DIMS, ERR)
      ELSE IF( PRESENT( R4 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_JPRB_REAL,    R4, DIMS, ERR)
      ELSE IF( PRESENT( RM4 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_JPRM_REAL,    RM4, DIMS, ERR)
      ELSE IF( PRESENT( IT1) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_STD_I8, TRANSFER(IT1, WI4, DIMS(1)), DIMS, ERR)
      ELSE IF( PRESENT( IT2 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_STD_I8, TRANSFER(IT2, WI4, DIMS(1)*DIMS(2)), DIMS, ERR)
      ELSE IF( PRESENT( IT3 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_STD_I8, TRANSFER(IT3, WI4, DIMS(1)*DIMS(2)*DIMS(3)), DIMS, ERR)
      ELSE IF( PRESENT( IS1 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_STD_I16, TRANSFER(IS1, WI4, DIMS(1)), DIMS, ERR)
      ELSE IF( PRESENT( IS2 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_STD_I16, TRANSFER(IS2, WI4, DIMS(1)*DIMS(2)), DIMS, ERR)
       ELSE IF( PRESENT( IS3 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_STD_I16, TRANSFER(IS3, WI4, DIMS(1)*DIMS(2)*DIMS(3)), DIMS, ERR)
      ELSE IF( PRESENT( IL1 ) ) THEN
        CALL H5DWRITE_F(DSET_ID, H5T_STD_I64, TRANSFER(IL1, WI4, DIMS(1)*2), DIMS, ERR)
      ELSE
        ERR = 1
      ENDIF
      THROWM(ERR.NE.0,"COULD NOT WRITE "//TRIM(NAME))

      IF( PRESENT( L1 ) ) THEN
        DEALLOCATE(LL1, STAT=ERR)
        THROWM(ERR.NE.0,"COULD NOTDEALLOCATE LL1")
      ENDIF

      CALL H5LTSET_ATTRIBUTE_STRING_F( LUN, NAME, "COMMENT", COMMENT//CHAR(0), ERR)
      THROWM(ERR.NE.0,"COULD NOT WRITE COMMENT IN "//TRIM(NAME))

      IF( PRESENT (UNITS) ) THEN
        CALL H5LTSET_ATTRIBUTE_STRING_F( LUN, NAME, "UNITS", UNITS//CHAR(0), ERR)
        THROWM(ERR.NE.0,"COULD NOT WRITE UNITS IN "//TRIM(NAME))
      ENDIF

      IF( .NOT. SCAL ) THEN
        RANK8 = RANK
        CALL H5LTSET_ATTRIBUTE_INT_F( LUN, NAME, "LBOUND", LB, RANK8, ERR)
        CALL H5LTSET_ATTRIBUTE_INT_F( LUN, NAME, "UBOUND", UB, RANK8, ERR)
      ENDIF

      IF( TYPE_F_ID .GE. 0 ) THEN
        CALL H5DCLOSE_F(DSET_ID, ERR)
        THROWM(ERR.NE.0,"COULD NOT CLOSE "//TRIM(NAME))
      ENDIF

CATCH

      END SUBROUTINE

!> @brief Read an array in an opened HDF5 file
!! @param[in]  LUN  HDF5 identification
!! @param[in]  NAME datatset name
!! @param[out] ERR status on exit
!! @param[out] R0 scalar of type real JPRB, optional
!! @param[out] R1 array of type real JPRB, dimension 1, optional
!! @param[out] R2 array of type real JPRB, dimension 2, optional
!! @param[out] R3 array of type real JPRB, dimension 3, optional
!! @param[out] R4 array of type real JPRB, dimension 4, optional
!! @param[out] RM0 scalar of type real JPRM, optional
!! @param[out] RM1 array of type real JPRM, dimension 1, optional
!! @param[out] RM2 array of type real JPRM, dimension 2, optional
!! @param[out] RM3 array of type real JPRM, dimension 3, optional
!! @param[out] RM4 array of type real JPRM, dimension 4, optional
!! @param[out] I0 scalar of type integer JPIM, optional
!! @param[out] I1 array of type integer JPIM, dimension 1, optional
!! @param[out] I2 array of type integer JPIM, dimension 2, optional
!! @param[out] I3 array of type integer JPIM, dimension 3, optional
!! @param[out] C0 string, optional
!! @param[out] C1 array of strings, dimension 1, optional
!! @param[out] L0 logical, optional
!! @param[out] L1 array of logicals, dimension 1, optional
!! @param PR1 array of type real pointer JPRB, dimension 1, optional
!! @param PR2 array of type real pointer JPRB, dimension 2, optional
!! @param PR3 array of type real pointer JPRB, dimension 3, optional
!! @param PR4 array of type real pointer JPRB, dimension 4, optional
!! @param PRM1 array of type real pointer JPRM, dimension 1, optional
!! @param PRM2 array of type real pointer JPRM, dimension 2, optional
!! @param PRM3 array of type real pointer JPRM, dimension 3, optional
!! @param PRM4 array of type real pointer JPRM, dimension 4, optional
!! @param PZ1 array of type real pointer complex JPRB, dimension 1, optional
!! @param PZ2 array of type real pointer complex JPRB, dimension 2, optional
!! @param PZ3 array of type real pointer complex JPRB, dimension 3, optional
!! @param PIT1 array of type integer pointer JPIT, dimension 1, optional
!! @param PIT2 array of type integer pointer JPIT, dimension 2, optional
!! @param PIT3 array of type integer pointer JPIT, dimension 3, optional
!! @param PIS1 array of type integer pointer JPIS, dimension 1, optional
!! @param PIS2 array of type integer pointer JPIS, dimension 2, optional
!! @param PIS3 array of type integer pointer JPIS, dimension 3, optional
!! @param PIS4 array of type integer pointer JPIS, dimension 4, optional
!! @param PI1 array of type integer pointer JPIM, dimension 1, optional
!! @param PI2 array of type integer pointer JPIM, dimension 2, optional
!! @param PI3 array of type integer pointer JPIM, dimension 3, optional
!! @param PIL1 array of type integer pointer JPIB, dimension 1, optional
!! @param PL1 array of logicals pointer, dimension 1, optional
!! @param PC1 array of strings pointer, dimension 1, optional
      SUBROUTINE READ_ARRAY_HDF( LUN, NAME, ERR, R0, R1, R2, R3, R4, &
       & RM0, RM1, RM2, RM3, RM4, &
       & I0, I1, I2, I3, &
       & C0, C1, L0, L1, &
       & PR1, PR2, PR3, PR4, PZ1, PZ2, PZ3, &
       & PRM1, PRM2, PRM3, PRM4, &
       & PIT1, PIT2, PIT3, PIS1, PIS2, PIS3, PIS4, &
       & PI1, PI2, PI3, PIL1,  &
       & PL1, PC1)
!
      INTEGER(HID_T), INTENT(IN) :: LUN
      CHARACTER(LEN=*),  INTENT(IN) :: NAME
      INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
      REAL(KIND=JPRB),    INTENT(OUT), OPTIONAL :: R0
      REAL(KIND=JPRB),    INTENT(OUT), OPTIONAL :: R1(:)
      REAL(KIND=JPRB),    INTENT(OUT), OPTIONAL :: R2(:,:)
      REAL(KIND=JPRB),    INTENT(OUT), OPTIONAL :: R3(:,:,:)
      REAL(KIND=JPRB),    INTENT(OUT), OPTIONAL :: R4(:,:,:,:)
      REAL(KIND=JPRM),    INTENT(OUT), OPTIONAL :: RM0
      REAL(KIND=JPRM),    INTENT(OUT), OPTIONAL :: RM1(:)
      REAL(KIND=JPRM),    INTENT(OUT), OPTIONAL :: RM2(:,:)
      REAL(KIND=JPRM),    INTENT(OUT), OPTIONAL :: RM3(:,:,:)
      REAL(KIND=JPRM),    INTENT(OUT), OPTIONAL :: RM4(:,:,:,:)
      INTEGER(KIND=JPIM), INTENT(OUT), OPTIONAL :: I0
      INTEGER(KIND=JPIM), INTENT(OUT), OPTIONAL :: I1(:)
      INTEGER(KIND=JPIM), INTENT(OUT), OPTIONAL :: I2(:,:)
      INTEGER(KIND=JPIM), INTENT(OUT), OPTIONAL :: I3(:,:,:)
      CHARACTER(LEN=*),   INTENT(OUT), OPTIONAL :: C0
      CHARACTER(LEN=*),   INTENT(OUT), OPTIONAL :: C1(:)
      LOGICAL,            INTENT(OUT), OPTIONAL :: L0
      LOGICAL,            INTENT(OUT), OPTIONAL :: L1(:)

      REAL(KIND=JPRB),    POINTER, OPTIONAL :: PR1(:)
      REAL(KIND=JPRB),    POINTER, OPTIONAL :: PR2(:,:)
      REAL(KIND=JPRB),    POINTER, OPTIONAL :: PR3(:,:,:)
      REAL(KIND=JPRB),    POINTER, OPTIONAL :: PR4(:,:,:,:)
      REAL(KIND=JPRM),    POINTER, OPTIONAL :: PRM1(:)
      REAL(KIND=JPRM),    POINTER, OPTIONAL :: PRM2(:,:)
      REAL(KIND=JPRM),    POINTER, OPTIONAL :: PRM3(:,:,:)
      REAL(KIND=JPRM),    POINTER, OPTIONAL :: PRM4(:,:,:,:)
      COMPLEX(KIND=JPRB), POINTER, OPTIONAL :: PZ1(:)
      COMPLEX(KIND=JPRB), POINTER, OPTIONAL :: PZ2(:,:)
      COMPLEX(KIND=JPRB), POINTER, OPTIONAL :: PZ3(:,:,:)
      INTEGER(KIND=JPIT), POINTER, OPTIONAL :: PIT1(:)
      INTEGER(KIND=JPIT), POINTER, OPTIONAL :: PIT2(:,:)
      INTEGER(KIND=JPIT), POINTER, OPTIONAL :: PIT3(:,:,:)
      INTEGER(KIND=JPIS), POINTER, OPTIONAL :: PIS1(:)
      INTEGER(KIND=JPIS), POINTER, OPTIONAL :: PIS2(:,:)
      INTEGER(KIND=JPIS), POINTER, OPTIONAL :: PIS3(:,:,:)
      INTEGER(KIND=JPIS), POINTER, OPTIONAL :: PIS4(:,:,:,:)
      INTEGER(KIND=JPIM), POINTER, OPTIONAL :: PI1(:)
      INTEGER(KIND=JPIM), POINTER, OPTIONAL :: PI2(:,:)
      INTEGER(KIND=JPIM), POINTER, OPTIONAL :: PI3(:,:,:)
      INTEGER(KIND=JPIB), POINTER, OPTIONAL :: PIL1(:)
      LOGICAL,            POINTER, OPTIONAL :: PL1(:)
      CHARACTER(LEN=*),   POINTER, OPTIONAL :: PC1(:)

      INTEGER(HSIZE_T) :: DIMS(7), MAXDIMS(7)
      INTEGER(KIND=JPIM) :: UB(7), LB(7)
      INTEGER(KIND=JPIM) :: I0X(1)
      REAL(KIND=JPRB) :: R0X(1)
      REAL(KIND=JPRM) :: RM0X(1)
      CHARACTER(LEN=256) :: FF
      INTEGER(HID_T) :: DSET_ID, SPCE_ID, MEMTYPE
      LOGICAL :: HASATT
      INTEGER(KIND=JPIM) :: I
      INTEGER(KIND=JPIM), ALLOCATABLE:: LL1(:)
        INTEGER(SIZE_T)  :: SDIM
#ifdef _RTTOV_HDF_F2003
      TYPE(C_PTR)      :: PTR
#else
      INTEGER(KIND=JPIM), ALLOCATABLE :: WI4(:)
#endif

TRY
      ERR = 0

      DIMS = 0

      LB = 0
      UB = 0

      IF( PRESENT( C0 ) ) THEN
        FF = ""
        CALL H5LTREAD_DATASET_STRING_F( LUN, NAME, FF, ERR)
        C0 = FF
        DO I = 1, LEN( C0 )
          IF( C0(I:I) .EQ. CHAR(0)) C0(I:I) = ' '
        ENDDO
      ELSE IF( PRESENT( R0 ) ) THEN
        DIMS(1) = 1
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRB_REAL, R0X, DIMS(1:1), ERR)
        R0 = R0X(1)
      ELSE IF( PRESENT( RM0 ) ) THEN
        DIMS(1) = 1
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRM_REAL, RM0X, DIMS(1:1), ERR)
        RM0 = RM0X(1)
      ELSE IF( PRESENT( I0 ) ) THEN
        DIMS(1) = 1
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_FILE_INTEGER, I0X, DIMS(1:1), ERR)
        I0 = I0X(1) 
      ELSE IF( PRESENT( L0 ) ) THEN
        DIMS(1) = 1
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_FILE_INTEGER, I0X, DIMS(1:1), ERR)
        L0 = I0X(1) .EQ. 1_JPIM
      ELSE IF( PRESENT( C1 ) ) THEN
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        DIMS(1) = SIZE( C1, 1 )
        SDIM = LEN( C1(1) )
        CALL H5TCOPY_F(H5T_FORTRAN_S1, MEMTYPE, ERR)
        THROWM(ERR.NE.0,"CANNOT H5TCOPY "//TRIM(NAME))
        CALL H5TSET_SIZE_F(MEMTYPE, SDIM, ERR)
        THROWM(ERR.NE.0,"CANNOT H5TSET_SIZE "//TRIM(NAME))
        CALL H5DREAD_F(DSET_ID, MEMTYPE, C1(1)(1:1), DIMS, ERR, SPCE_ID)
        THROWM(ERR.LT.0,"CANNOT H5DREAD_F "//TRIM(NAME))
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
      ELSE IF( PRESENT( R1 ) ) THEN
        DIMS(1) = SIZE( R1, 1 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRB_REAL, R1, DIMS(1:1), ERR)
      ELSE IF( PRESENT( RM1 ) ) THEN
        DIMS(1) = SIZE( RM1, 1 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRM_REAL, RM1, DIMS(1:1), ERR)
      ELSE IF( PRESENT( I1 ) ) THEN
        DIMS(1) = SIZE( I1, 1 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_FILE_INTEGER, I1, DIMS(1:1), ERR)
      ELSE IF( PRESENT( L1 ) ) THEN
        DIMS(1) = SIZE( L1, 1 )
        ALLOCATE( LL1(DIMS(1)), STAT=ERR)
        THROWM(ERR.NE.0,"COULD NOT ALLOCATE LL1 ")
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_FILE_INTEGER, LL1, DIMS(1:1), ERR)
        WHERE( LL1 .EQ. 1_JPIM )
          L1 = .TRUE.
        ELSEWHERE
          L1 = .FALSE.
        END WHERE
        DEALLOCATE( LL1, STAT=ERR)
        THROWM(ERR.NE.0,"COULD NOT DEALLOCATE LL1 ")
      ELSE IF( PRESENT( R2 ) ) THEN
        DIMS(1) = SIZE( R2, 1 )
        DIMS(2) = SIZE( R2, 2 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRB_REAL, R2, DIMS(1:2), ERR)
      ELSE IF( PRESENT( RM2 ) ) THEN
        DIMS(1) = SIZE( RM2, 1 )
        DIMS(2) = SIZE( RM2, 2 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRM_REAL, RM2, DIMS(1:2), ERR)
      ELSE IF( PRESENT( I2 ) ) THEN
        DIMS(1) = SIZE( I2, 1 )
        DIMS(2) = SIZE( I2, 2 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_FILE_INTEGER, I2, DIMS(1:2), ERR)
      ELSE IF( PRESENT( R3 ) ) THEN
        DIMS(1) = SIZE( R3, 1 )
        DIMS(2) = SIZE( R3, 2 )
        DIMS(3) = SIZE( R3, 3 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRB_REAL, R3, DIMS(1:3), ERR)
      ELSE IF( PRESENT( RM3 ) ) THEN
        DIMS(1) = SIZE( RM3, 1 )
        DIMS(2) = SIZE( RM3, 2 )
        DIMS(3) = SIZE( RM3, 3 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRM_REAL, RM3, DIMS(1:3), ERR)
      ELSE IF( PRESENT( R4 ) ) THEN
        DIMS(1) = SIZE( R4, 1 )
        DIMS(2) = SIZE( R4, 2 )
        DIMS(3) = SIZE( R4, 3 )
        DIMS(4) = SIZE( R4, 4 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRB_REAL, R4, DIMS(1:4), ERR)
      ELSE IF( PRESENT( RM4 ) ) THEN
        DIMS(1) = SIZE( RM4, 1 )
        DIMS(2) = SIZE( RM4, 2 )
        DIMS(3) = SIZE( RM4, 3 )
        DIMS(4) = SIZE( RM4, 4 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRM_REAL, RM4, DIMS(1:4), ERR)
      ELSE IF( PRESENT( I3 ) ) THEN
        DIMS(1) = SIZE( I3, 1 )
        DIMS(2) = SIZE( I3, 2 )
        DIMS(3) = SIZE( I3, 3 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_FILE_INTEGER, I3, DIMS(1:3), ERR)
      ELSE IF( PRESENT( PR1 ) ) THEN
      ELSE IF( PRESENT( PRM1) ) THEN
      ELSE IF( PRESENT( PZ1 ) ) THEN
      ELSE IF( PRESENT( PI1 ) ) THEN
      ELSE IF( PRESENT( PL1 ) ) THEN
      ELSE IF( PRESENT( PC1 ) ) THEN
      ELSE IF( PRESENT( PR2 ) ) THEN
      ELSE IF( PRESENT( PRM2) ) THEN
      ELSE IF( PRESENT( PZ2 ) ) THEN
      ELSE IF( PRESENT( PI2 ) ) THEN
      ELSE IF( PRESENT( PR3 ) ) THEN
      ELSE IF( PRESENT( PRM3) ) THEN
      ELSE IF( PRESENT( PZ3 ) ) THEN
      ELSE IF( PRESENT( PI3 ) ) THEN
      ELSE IF( PRESENT( PR4 ) ) THEN
      ELSE IF( PRESENT( PRM4) ) THEN
      ELSE IF( PRESENT( PIT1) ) THEN
      ELSE IF( PRESENT( PIT2) ) THEN
      ELSE IF( PRESENT( PIT3) ) THEN
      ELSE IF( PRESENT( PIS1) ) THEN
      ELSE IF( PRESENT( PIS2) ) THEN
      ELSE IF( PRESENT( PIS3) ) THEN
      ELSE IF( PRESENT( PIS4) ) THEN
      ELSE IF( PRESENT( PIL1) ) THEN
      ELSE
        ERR = -1
      ENDIF

      IF( ERR .GE. 0 ) ERR = 0
      THROWM(ERR.NE.0,"COULD NOT READ "//TRIM(NAME))

!     CALL H5ESET_AUTO_F( 0, SAERR )

      CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
      THROWM(ERR.NE.0,"COULD NOT OPEN "//TRIM(NAME))
      CALL H5AEXISTS_F( DSET_ID, "LBOUND", HASATT, ERR )
      THROWM(ERR.NE.0,"COULD NOT READ "//TRIM(NAME))
      CALL H5DCLOSE_F( DSET_ID, ERR )
      THROWM(ERR.NE.0,"COULD NOT CLOSE "//TRIM(NAME))

      IF( HASATT ) THEN
        CALL H5LTGET_ATTRIBUTE_INT_F( LUN, NAME, "LBOUND", LB, ERR )
        THROWM(ERR.NE.0,"CANNOT READ LBOUND/"//TRIM(NAME))
        CALL H5LTGET_ATTRIBUTE_INT_F( LUN, NAME, "UBOUND", UB, ERR )
        THROWM(ERR.NE.0,"CANNOT READ UBOUND/"//TRIM(NAME))
      ELSE
        ERR = 0
        LB = 1
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        CALL H5SGET_SIMPLE_EXTENT_DIMS_F( SPCE_ID, DIMS, MAXDIMS, ERR ) 
        IF( ERR .GE. 0 ) ERR = 0
        THROWM(ERR.NE.0,"CANNOT GET DIMS OF "//TRIM(NAME))
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        UB = INT( DIMS, JPIM )
      ENDIF

      IF( PRESENT( PR1 ) ) THEN
        ALLOCATE( PR1(LB(1):UB(1)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        DIMS(1) = SIZE( PR1, 1 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRB_REAL, PR1, DIMS(1:1), ERR)
      ELSE IF( PRESENT( PRM1 ) ) THEN
        ALLOCATE( PRM1(LB(1):UB(1)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        DIMS(1) = SIZE( PRM1, 1 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRM_REAL, PRM1, DIMS(1:1), ERR)
      ELSE IF( PRESENT( PI1 ) ) THEN
        ALLOCATE( PI1(LB(1):UB(1)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        DIMS(1) = SIZE( PI1, 1 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_FILE_INTEGER, PI1, DIMS(1:1), ERR)
      ELSE IF( PRESENT( PL1 ) ) THEN
        ALLOCATE( PL1(LB(1):UB(1)), LL1(LB(1):UB(1)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        DIMS(1) = SIZE( PL1, 1 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_FILE_INTEGER, LL1, DIMS(1:1), ERR)
        WHERE( LL1 .EQ. 1_JPIM )
          PL1 = .TRUE.
        ELSEWHERE
          PL1 = .FALSE.
        END WHERE
        DEALLOCATE( LL1, STAT=ERR)
        THROWM(ERR.NE.0,"COULD NOT DEALLOCATE LL1 ")
      ELSE IF( PRESENT( PC1 ) ) THEN
        ALLOCATE( PC1(LB(1):UB(1)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        DIMS(1) = SIZE( PC1, 1 )
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        SDIM = LEN( PC1(1) )
        CALL H5TCOPY_F(H5T_FORTRAN_S1, MEMTYPE, ERR)
        THROWM(ERR.NE.0,"CANNOT H5TCOPY "//TRIM(NAME))
        CALL H5TSET_SIZE_F(MEMTYPE, SDIM, ERR)
        THROWM(ERR.NE.0,"CANNOT H5TSET_SIZE "//TRIM(NAME))
        CALL H5DREAD_F(DSET_ID, MEMTYPE, PC1(1)(1:1), DIMS, ERR, SPCE_ID)
        THROWM(ERR.LT.0,"CANNOT H5DREAD_F "//TRIM(NAME))
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
      ELSE IF( PRESENT( PR2 ) ) THEN
        ALLOCATE( PR2(LB(1):UB(1),LB(2):UB(2)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        DIMS(1) = SIZE( PR2, 1 )
        DIMS(2) = SIZE( PR2, 2 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRB_REAL, PR2, DIMS(1:2), ERR)
      ELSE IF( PRESENT( PRM2 ) ) THEN
        ALLOCATE( PRM2(LB(1):UB(1),LB(2):UB(2)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        DIMS(1) = SIZE( PRM2, 1 )
        DIMS(2) = SIZE( PRM2, 2 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRM_REAL, PRM2, DIMS(1:2), ERR)
      ELSE IF( PRESENT( PI2 ) ) THEN
        ALLOCATE( PI2(LB(1):UB(1),LB(2):UB(2)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        DIMS(1) = SIZE( PI2, 1 )
        DIMS(2) = SIZE( PI2, 2 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_FILE_INTEGER, PI2, DIMS(1:2), ERR)
      ELSE IF( PRESENT( PI3 ) ) THEN
        ALLOCATE( PI3(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        DIMS(1) = SIZE( PI3, 1 )
        DIMS(2) = SIZE( PI3, 2 )
        DIMS(3) = SIZE( PI3, 3 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_FILE_INTEGER, PI3, DIMS(1:3), ERR)
      ELSE IF( PRESENT( PR3 ) ) THEN
        ALLOCATE( PR3(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        DIMS(1) = SIZE( PR3, 1 )
        DIMS(2) = SIZE( PR3, 2 )
        DIMS(3) = SIZE( PR3, 3 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRB_REAL, PR3, DIMS(1:3), ERR)
      ELSE IF( PRESENT( PRM3 ) ) THEN
        ALLOCATE( PRM3(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        DIMS(1) = SIZE( PRM3, 1 )
        DIMS(2) = SIZE( PRM3, 2 )
        DIMS(3) = SIZE( PRM3, 3 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRM_REAL, PRM3, DIMS(1:3), ERR)
      ELSE IF( PRESENT( PR4 ) ) THEN
        ALLOCATE( PR4(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        DIMS(1) = SIZE( PR4, 1 )
        DIMS(2) = SIZE( PR4, 2 )
        DIMS(3) = SIZE( PR4, 3 )
        DIMS(4) = SIZE( PR4, 4 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRB_REAL, PR4, DIMS(1:4), ERR)
      ELSE IF( PRESENT( PRM4 ) ) THEN
        ALLOCATE( PRM4(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        DIMS(1) = SIZE( PRM4, 1 )
        DIMS(2) = SIZE( PRM4, 2 )
        DIMS(3) = SIZE( PRM4, 3 )
        DIMS(4) = SIZE( PRM4, 4 )
        CALL H5LTREAD_DATASET_F( LUN, NAME, H5T_JPRM_REAL, PRM4, DIMS(1:4), ERR)

#ifdef _RTTOV_HDF_F2003
      ELSE IF( PRESENT( PIT1 ) ) THEN
        DIMS(1) = UB(1)-LB(1)+1
        ALLOCATE( PIT1(LB(1):UB(1)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        PTR = C_LOC(PIT1(1))
        CALL H5DREAD_F(DSET_ID, H5T_STD_I8, PTR, err)
        THROWM(ERR.LT.0,"CANNOT H5DREAD_F "//TRIM(NAME))
        CALL H5ECLEAR_F(err)
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
      ELSE IF( PRESENT( PIT2 ) ) THEN
        DIMS(1) = UB(1)-LB(1)+1
        DIMS(2) = UB(2)-LB(2)+1
        ALLOCATE( PIT2(LB(1):UB(1),LB(2):UB(2)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        PTR = C_LOC(PIT2(1,1))
        CALL H5DREAD_F(DSET_ID, H5T_STD_I8, PTR, err)
        THROWM(ERR.LT.0,"CANNOT H5DREAD_F "//TRIM(NAME))
        CALL H5ECLEAR_F(err)
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
      ELSE IF( PRESENT( PIT3 ) ) THEN
        DIMS(1) = UB(1)-LB(1)+1
        DIMS(2) = UB(2)-LB(2)+1
        DIMS(3) = UB(3)-LB(3)+1
        ALLOCATE( PIT3(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        PTR = C_LOC(PIT3(1,1,1))
        CALL H5DREAD_F(DSET_ID, H5T_STD_I8, PTR, err)
        THROWM(ERR.LT.0,"CANNOT H5DREAD_F "//TRIM(NAME))
        CALL H5ECLEAR_F(err)
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
      ELSE IF( PRESENT( PIS1 ) ) THEN
        DIMS(1) = UB(1)-LB(1)+1
        ALLOCATE( PIS1(LB(1):UB(1)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        PTR = C_LOC(PIS1(1))
        CALL H5DREAD_F(DSET_ID, H5T_STD_I16, PTR, err)
        THROWM(ERR.LT.0,"CANNOT H5DREAD_F "//TRIM(NAME))
        CALL H5ECLEAR_F(err)
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
      ELSE IF( PRESENT( PIS2 ) ) THEN
        DIMS(1) = UB(1)-LB(1)+1
        DIMS(2) = UB(2)-LB(2)+1
        ALLOCATE( PIS2(LB(1):UB(1),LB(2):UB(2)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        PTR = C_LOC(PIS2(1,1))
        CALL H5DREAD_F(DSET_ID, H5T_STD_I16, PTR, err)
        THROWM(ERR.LT.0,"CANNOT H5DREAD_F "//TRIM(NAME))
        CALL H5ECLEAR_F(err)
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
      ELSE IF( PRESENT( PIS3 ) ) THEN
        DIMS(1) = UB(1)-LB(1)+1
        DIMS(2) = UB(2)-LB(2)+1
        DIMS(3) = UB(3)-LB(3)+1
        ALLOCATE( PIS3(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        PTR = C_LOC(PIS3(1,1,1))
        CALL H5DREAD_F(DSET_ID, H5T_STD_I16, PTR, err)
        THROWM(ERR.LT.0,"CANNOT H5DREAD_F "//TRIM(NAME))
        CALL H5ECLEAR_F(err)
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
      ELSE IF( PRESENT( PIS4 ) ) THEN
        DIMS(1) = UB(1)-LB(1)+1
        DIMS(2) = UB(2)-LB(2)+1
        DIMS(3) = UB(3)-LB(3)+1
        DIMS(4) = UB(4)-LB(4)+1
        ALLOCATE( PIS4(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        PTR = C_LOC(PIS4(1,1,1,1))
        CALL H5DREAD_F(DSET_ID, H5T_STD_I16, PTR, err)
        THROWM(ERR.LT.0,"CANNOT H5DREAD_F "//TRIM(NAME))
        CALL H5ECLEAR_F(err)
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
      ELSE IF( PRESENT( PIL1 ) ) THEN
        ALLOCATE( PIL1(LB(1):UB(1)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        DIMS(1) = SIZE( PIL1, 1 )
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        PTR = C_LOC(PIL1(1))
        CALL H5DREAD_F(DSET_ID, H5T_STD_I64, PTR, err)
        THROWM(ERR.LT.0,"CANNOT H5DREAD_F "//TRIM(NAME))
        CALL H5ECLEAR_F(err)
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))

#else
     ELSE IF( PRESENT( PIT1 ) ) THEN
        DIMS(1) = UB(1)-LB(1)+1
        ALLOCATE( PIT1(LB(1):UB(1)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        ALLOCATE( WI4(1:DIMS(1)/4+1), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        CALL h5dread_f(dset_id, H5T_STD_I8, WI4, dims, err)
        THROWM(ERR.LT.0,"CANNOT h5dread_f "//TRIM(NAME))
        PIT1  = TRANSFER(WI4, PIT1, DIMS(1) )
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        DEALLOCATE( WI4, STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT DEALLOCATE "//TRIM(NAME))
      ELSE IF( PRESENT( PIT2 ) ) THEN
        DIMS(1) = UB(1)-LB(1)+1
        DIMS(2) = UB(2)-LB(2)+1
        ALLOCATE( PIT2(LB(1):UB(1),LB(2):UB(2)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        ALLOCATE( WI4(1:DIMS(1)*DIMS(2)/4+1) , STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        CALL h5dread_f(dset_id, H5T_STD_I8, WI4, dims, err)
        THROWM(ERR.LT.0,"CANNOT h5dread_f "//TRIM(NAME))
        PIT2  = RESHAPE(TRANSFER(WI4, PIT2, DIMS(1)*DIMS(2) ), DIMS(1:2))
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        DEALLOCATE( WI4, STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT DEALLOCATE "//TRIM(NAME))
      ELSE IF( PRESENT( PIT3 ) ) THEN
        DIMS(1) = UB(1)-LB(1)+1
        DIMS(2) = UB(2)-LB(2)+1
        DIMS(3) = UB(3)-LB(3)+1
        ALLOCATE( PIT3(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        ALLOCATE( WI4(1:DIMS(1)*DIMS(2)*DIMS(3)/4+1) , STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        CALL h5dread_f(dset_id, H5T_STD_I8, WI4, dims, err)
        THROWM(ERR.LT.0,"CANNOT h5dread_f "//TRIM(NAME))
        PIT3  = RESHAPE(TRANSFER(WI4, PIT3, DIMS(1)*DIMS(2)*DIMS(3) ), DIMS(1:3))
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        DEALLOCATE( WI4, STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT DEALLOCATE "//TRIM(NAME))
      ELSE IF( PRESENT( PIS1 ) ) THEN
        DIMS(1) = UB(1)-LB(1)+1
        ALLOCATE( PIS1(LB(1):UB(1)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        ALLOCATE( WI4(1:DIMS(1)/2+1), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        CALL h5dread_f(dset_id, H5T_STD_I16, WI4, dims, err)
        THROWM(ERR.LT.0,"CANNOT h5dread_f "//TRIM(NAME))
        PIS1  = TRANSFER(WI4, PIS1, DIMS(1) )
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        DEALLOCATE( WI4, STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT DEALLOCATE "//TRIM(NAME))
      ELSE IF( PRESENT( PIS2 ) ) THEN
        DIMS(1) = UB(1)-LB(1)+1
        DIMS(2) = UB(2)-LB(2)+1
        ALLOCATE( PIS2(LB(1):UB(1),LB(2):UB(2)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        ALLOCATE( WI4(1:DIMS(1)*DIMS(2)/2+1) , STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        CALL h5dread_f(dset_id, H5T_STD_I16, WI4, dims, err)
        THROWM(ERR.LT.0,"CANNOT h5dread_f "//TRIM(NAME))
        PIS2  = RESHAPE(TRANSFER(WI4, PIS2, DIMS(1)*DIMS(2) ), DIMS(1:2))
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        DEALLOCATE( WI4, STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT DEALLOCATE "//TRIM(NAME))
      ELSE IF( PRESENT( PIS3 ) ) THEN
        DIMS(1) = UB(1)-LB(1)+1
        DIMS(2) = UB(2)-LB(2)+1
        DIMS(3) = UB(3)-LB(3)+1
        ALLOCATE( PIS3(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        ALLOCATE( WI4(1:DIMS(1)*DIMS(2)*DIMS(3)/2+1) , STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        CALL h5dread_f(dset_id, H5T_STD_I16, WI4, dims, err)
        THROWM(ERR.LT.0,"CANNOT h5dread_f "//TRIM(NAME))
        PIS3  = RESHAPE(TRANSFER(WI4, PIS3, DIMS(1)*DIMS(2)*DIMS(3) ), DIMS(1:3))
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        DEALLOCATE( WI4, STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT DEALLOCATE "//TRIM(NAME))
      ELSE IF( PRESENT( PIS4 ) ) THEN
        DIMS(1) = UB(1)-LB(1)+1
        DIMS(2) = UB(2)-LB(2)+1
        DIMS(3) = UB(3)-LB(3)+1
        DIMS(4) = UB(4)-LB(4)+1
        ALLOCATE( PIS4(LB(1):UB(1),LB(2):UB(2),LB(3):UB(3),LB(4):UB(4)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        ALLOCATE( WI4(1:DIMS(1)*DIMS(2)*DIMS(3)*DIMS(4)/2+1) , STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        CALL h5dread_f(dset_id, H5T_STD_I16, WI4, dims, err)
        THROWM(ERR.LT.0,"CANNOT h5dread_f "//TRIM(NAME))
        PIS4  = RESHAPE(TRANSFER(WI4, PIS4, DIMS(1)*DIMS(2)*DIMS(3)*DIMS(4) ), DIMS(1:4))
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        DEALLOCATE( WI4, STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT DEALLOCATE "//TRIM(NAME))
      ELSE IF( PRESENT( PIL1 ) ) THEN
        ALLOCATE( PIL1(LB(1):UB(1)), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        ALLOCATE( WI4(1:(UB(1)-LB(1)+1)*2), STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT ALLOCATE "//TRIM(NAME))
        DIMS(1) = SIZE( PIL1, 1 )
        CALL H5DOPEN_F( LUN, NAME, DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT OPEN "//TRIM(NAME))
        CALL H5DGET_SPACE_F( DSET_ID, SPCE_ID, ERR ) 
        THROWM(ERR.NE.0,"CANNOT GET SPACE OF "//TRIM(NAME))
        CALL h5dread_f(dset_id, H5T_STD_I64, WI4, dims, err)
        THROWM(ERR.LT.0,"CANNOT h5dread_f "//TRIM(NAME))
        PIL1  = TRANSFER(WI4, PIL1, DIMS(1) )
        CALL H5SCLOSE_F( SPCE_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        CALL H5DCLOSE_F( DSET_ID, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE "//TRIM(NAME))
        DEALLOCATE( WI4, STAT = ERR )
        THROWM(ERR.NE.0,"CANNOT DEALLOCATE "//TRIM(NAME))
#endif

      ELSE
        ERR = 1
      ENDIF
      IF( ERR .GE. 0 ) ERR = 0
      THROWM(ERR.NE.0,"COULD NOT READ "//TRIM(NAME))

CATCH

      END SUBROUTINE

!> @brief Create an HDF5 file
!! @param[out] ERR status on exit
!! @param[in]  F filename
!! @param[out] FILE_ID HDF5 file identification
      SUBROUTINE RTTOV_HDF_OPEN_TRUNC( ERR, F, FILE_ID)
!
! 
      CHARACTER(LEN=*),    INTENT(IN)  :: F
      INTEGER(KIND=JPIM),  INTENT(OUT) :: ERR
      INTEGER(HID_T),      INTENT(OUT) :: FILE_ID

      CHARACTER(LEN=8)   :: date
      CHARACTER(LEN=10)  :: time
      CHARACTER(LEN=21)  :: datetime
!
TRY

      CALL DATE_AND_TIME(date, time)

      WRITE(datetime,"(1X,a4,2('/',a2),2x,2(a2,':'),a2)")   &
              & date(1:4), date(5:6), date(7:8),        &
              & time(1:2), time(3:4), time(5:6)
              
      CALL H5FCREATE_F( F, H5F_ACC_TRUNC_F, FILE_ID, ERR )
      THROWM(ERR.NE.0,"CANNOT CREATE "//TRIM(F))

      CALL H5LTSET_ATTRIBUTE_STRING_F(FILE_ID, '.', "Creation date",   &
           TRIM(datetime)// &
           CHAR(0), ERR )
      THROWM(ERR.NE.0,"CANNOT WRITE ATTRIBUTE")

CATCH
      END SUBROUTINE

!> @brief Open an HDF5 file for read/write
!! @param[out] ERR status on exit
!! @param[in]  F filename
!! @param[out] FILE_ID HDF5 file identification
!! @param[in]  LWARN control the error messages, default true, optional
      SUBROUTINE RTTOV_HDF_OPEN_RDWR( ERR, F, FILE_ID, LWARN)
!
! 
      CHARACTER(LEN=*),    INTENT(IN)  :: F
      INTEGER(KIND=JPIM),  INTENT(OUT) :: ERR
      INTEGER(HID_T),      INTENT(OUT) :: FILE_ID
      LOGICAL(KIND=JPLM), OPTIONAL, INTENT(IN) :: LWARN
      
      CHARACTER(LEN=8)   :: date
      CHARACTER(LEN=10)  :: time
      CHARACTER(LEN=21)  :: datetime
      LOGICAL(KIND=JPLM) :: LWARN1

!
TRY

      IF( PRESENT(LWARN) ) THEN
        LWARN1 = LWARN
      ELSE
        LWARN1 = .TRUE.
      ENDIF

      CALL DATE_AND_TIME(date, time)

      WRITE(datetime,"(1X,a4,2('/',a2),2x,2(a2,':'),a2)")   &
              & date(1:4), date(5:6), date(7:8),        &
              & time(1:2), time(3:4), time(5:6)

      CALL H5FOPEN_F( F, H5F_ACC_RDWR_F, FILE_ID, ERR )
      THROWM(ERR.NE.0 .AND. LWARN1,"CANNOT OPEN (READ) "//TRIM(F))

      CALL H5LTSET_ATTRIBUTE_STRING_F(FILE_ID, '.', "Last update date",   &
           TRIM(datetime)// &
           CHAR(0), ERR )
      THROWM(ERR.NE.0 .AND. LWARN1,"CANNOT WRITE ATTRIBUTE")
   
CATCH
      END SUBROUTINE

!> @brief Open an HDF5 file for write: if file exists it is not overwritten, otherwise it is created
!! @param[out] ERR status on exit
!! @param[in]  F filename
!! @param[out] FILE_ID HDF5 file identification
      SUBROUTINE RTTOV_HDF_OPENW( ERR, F, FILE_ID)
! 
      CHARACTER(LEN=*),    INTENT(IN)  :: F
      INTEGER(KIND=JPIM),  INTENT(OUT) :: ERR
      INTEGER(HID_T),      INTENT(OUT) :: FILE_ID
      
      INTEGER :: SAERR
!
TRY

      CALL H5ESET_AUTO_F( 0, SAERR )
      CALL RTTOV_HDF_OPEN_RDWR( ERR, F, FILE_ID, LWARN=.FALSE.)
      CALL H5ESET_AUTO_F( 1, SAERR )

      IF( ERR .NE. 0 ) THEN
            CALL H5ESET_AUTO_F( 0, SAERR )
            CALL H5FCLOSE_F( FILE_ID, ERR )
            CALL H5ESET_AUTO_F( 1, SAERR )
        CALL RTTOV_HDF_OPEN_TRUNC( ERR, F, FILE_ID)
      ENDIF
      THROWM(ERR.NE.0,"CANNOT OPEN (WRITE) "//TRIM(F))

CATCH
      END SUBROUTINE

!> @brief Open an HDF5 file for read only
!! @param[out] ERR status on exit
!! @param[in]  F filename
!! @param[out] FILE_ID HDF5 file identification
      SUBROUTINE RTTOV_HDF_OPENR( ERR, F, FILE_ID)
! 
      CHARACTER(LEN=*),    INTENT(IN)  :: F
      INTEGER(KIND=JPIM),  INTENT(OUT) :: ERR
      INTEGER(HID_T),      INTENT(OUT) :: FILE_ID
!
TRY
      CALL H5FOPEN_F( F, H5F_ACC_RDONLY_F, FILE_ID, ERR )
      THROWM(ERR.NE.0,"CANNOT OPEN (READ) "//TRIM(F))
CATCH
      END SUBROUTINE

!> @brief Write an array of complex in an opened HDF5 file
!! @param[in]  LUN  HDF5 identification
!! @param[in]  NAME datatset name
!! @param[in]  COMMENT comment for dataset
!! @param[out] ERR status on exit
!! @param[in]  K1 array of type real complex JPRB, dimension 1
!! @param[in]  UNITS dataset units attribute, optional
      SUBROUTINE WRITE_ARRAY_HDF_CMPLX( LUN, NAME, COMMENT, ERR, K1, UNITS )
      
      INTEGER(HID_T),     INTENT(IN)  :: LUN
      CHARACTER(LEN=*),   INTENT(IN)  :: NAME
      CHARACTER(LEN=*),   INTENT(IN)  :: COMMENT
      INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
      COMPLEX(KIND=JPRB), INTENT(IN)  :: K1(:)
      CHARACTER(LEN=*),   INTENT(IN) , OPTIONAL :: UNITS

      INTEGER(HID_T) :: G_ID_SUB
!
TRY
     CALL MKPAR( LUN, TRIM(NAME), G_ID_SUB, ERR )
     THROWM(ERR.NE.0,"CANNOT CREATE PATH FOR COMPLEX "//TRIM(NAME))

     CALL H5LTSET_ATTRIBUTE_STRING_F( LUN, NAME, "COMMENT", COMMENT//CHAR(0), ERR)
      THROWM(ERR.NE.0,"COULD NOT WRITE COMMENT IN "//TRIM(NAME))

     IF( PRESENT (UNITS) ) THEN
       CALL H5LTSET_ATTRIBUTE_STRING_F( LUN, NAME, "UNITS", UNITS//CHAR(0), ERR)
       THROWM(ERR.NE.0,"COULD NOT WRITE UNITS IN "//TRIM(NAME))
     ENDIF

     call write_array_hdf(G_ID_SUB,'REAL','real part of complex',err,r1=real(K1), units=units )
     THROWM(ERR.NE.0,"CANNOT WRITE REAL PATH FOR COMPLEX "//TRIM(NAME))

     call write_array_hdf(G_ID_SUB,'IMAG','imaginary part of complex',err,r1=aimag(K1), units=units )
     THROWM(ERR.NE.0,"CANNOT WRITE IMAGINARY PATH FOR COMPLEX "//TRIM(NAME))

     CALL H5GCLOSE_F( G_ID_SUB, ERR )
     THROWM(ERR.NE.0,"CANNOT CLOSE PATH FOR COMPLEX "//TRIM(NAME))

CATCH
      END SUBROUTINE

!> @brief Read an array of complex in an opened HDF5 file
!! @param[in] LUN  HDF5 identification
!! @param[in]  NAME datatset name
!! @param[out] ERR status on exit
!! @param[out] K1 array of type real complex JPRB, dimension 1, optional
!! @param      PK1 array of type real pointer complex JPRB, dimension 1, optional
      SUBROUTINE READ_ARRAY_HDF_CMPLX( LUN, NAME, ERR, K1, PK1 )
!
      INTEGER(HID_T),    INTENT(IN) :: LUN
      CHARACTER(LEN=*),  INTENT(IN) :: NAME
      INTEGER(KIND=JPIM), INTENT(OUT) :: ERR
      COMPLEX(KIND=JPRB), INTENT(OUT), OPTIONAL :: K1(:)
      COMPLEX(KIND=JPRB), POINTER, OPTIONAL :: PK1(:)
!
      INTEGER(HID_T) :: G_ID_SUB
      REAL(KIND=JPRB), pointer :: R1(:), I1(:)

TRY

     CALL H5GOPEN_F( LUN, TRIM(NAME), G_ID_SUB, ERR )
     THROWM(ERR.NE.0,"CANNOT OPEN PATH FOR COMPLEX "//TRIM(NAME))

     call read_array_hdf(G_ID_SUB,'REAL',err,pr1=R1)
     THROWM(ERR.NE.0,"CANNOT READ REAL PATH FOR COMPLEX "//TRIM(NAME))

     call read_array_hdf(G_ID_SUB,'IMAG',err,pr1=I1)
     THROWM(ERR.NE.0,"CANNOT READ IMAGINARY PATH FOR COMPLEX "//TRIM(NAME))

     CALL H5GCLOSE_F( G_ID_SUB, ERR )
     THROWM(ERR.NE.0,"CANNOT CLOSE PATH FOR COMPLEX "//TRIM(NAME))

     IF( PRESENT(K1) ) THEN

       IF( size(K1) .EQ. size(R1) ) THEN
         K1 = cmplx( r1, i1, JPRB )
       ELSE
         ERR = 1_JPIM
         THROWM(ERR.NE.0,"WRONG SIZE FOR COMPLEX ARRAY "//TRIM(NAME))
      ENDIF

     ELSEIF( PRESENT(PK1) ) THEN

       Allocate(PK1(size(r1)), stat=err)
       THROWM(ERR.NE.0,"CANNOT ALLOCATE COMPLEX "//TRIM(NAME))

       PK1 = cmplx( r1, i1, JPRB )

       DeAllocate( r1, i1)

     ENDIF

CATCH
      END SUBROUTINE

!> @brief Check if an HDF5 group exists in opened file, group_id is closed
!! @param[in]  L_ID HDF5 file identification
!! @param[in]  NAME group name to check
!! @param[out] HDFERR return code (0 = group found)
      SUBROUTINE EXGRP( L_ID, NAME, HDFERR )
!
      INTEGER(HID_T),   INTENT(IN)  :: L_ID
      CHARACTER(LEN=*), INTENT(IN)  :: NAME
      INTEGER,          INTENT(OUT) :: HDFERR
!
      INTEGER :: I1, I2
      INTEGER :: SAERR, HDFERR1
      INTEGER(HID_T) :: L1_ID, L2_ID
      INTEGER(KIND=JPIM) :: N
!
      HDFERR = 0
      L1_ID = L_ID
      N = LEN( TRIM(NAME) )

      I1 = 1
      DO1: DO 
        DO2: DO WHILE( NAME(I1:I1) .EQ. '/' ) 
           I1 = I1 + 1
           IF( I1 .GT. N ) EXIT DO1
        ENDDO DO2
        I2 = I1 + 1
        DO3: DO WHILE( NAME(I2:I2) .NE. '/' ) 
           I2 = I2 + 1
           IF( I2 .GT. N ) EXIT DO3
        ENDDO DO3
        I2 = I2 - 1
        CALL H5ESET_AUTO_F( 0, SAERR )

        CALL H5GOPEN_F( L1_ID, NAME(I1:I2), L2_ID, HDFERR )
        CALL H5GCLOSE_F( L1_ID, HDFERR1 )

        IF( HDFERR .NE. 0 ) THEN 
          CALL H5ESET_AUTO_F( 1, SAERR )
          RETURN
        ENDIF
          
        I1    = I2 + 1
        L1_ID = L2_ID
        IF( I1 .GT. N ) EXIT DO1
      ENDDO DO1
      
      CALL H5ESET_AUTO_F( 0, SAERR )
      CALL H5GCLOSE_F( L1_ID, HDFERR1 )
      CALL H5ESET_AUTO_F( 1, SAERR )

      END SUBROUTINE

!> @brief Returns a datagroup ID and creates it if does not exist
!! @param[in]  L_ID top HDF5 identification group
!! @param[in]  NAME group or datatset name
!! @param[out] G_ID output group id
!! @param[out] HDFERR return code (0 = OK)
      SUBROUTINE MKPAR( L_ID, NAME, G_ID, HDFERR )
!
      INTEGER(HID_T),   INTENT(IN)  :: L_ID
      CHARACTER(LEN=*), INTENT(IN)  :: NAME
      INTEGER(HID_T),   INTENT(OUT) :: G_ID  
      INTEGER,          INTENT(OUT) :: HDFERR
!
      INTEGER :: I1, I2
      INTEGER :: SAERR
      INTEGER(HID_T) :: L1_ID, L2_ID
      INTEGER(KIND=JPIM) :: N
!
      L1_ID = L_ID
      N = LEN( TRIM(NAME) )

!     PRINT *, " NAME = ", TRIM(NAME)

      I1 = 1
      DO1: DO 
        DO2: DO WHILE( NAME(I1:I1) .EQ. '/' ) 
           I1 = I1 + 1
           IF( I1 .GT. N ) EXIT DO1
        ENDDO DO2
        I2 = I1 + 1
        DO3: DO WHILE( NAME(I2:I2) .NE. '/' ) 
           I2 = I2 + 1
           IF( I2 .GT. N ) EXIT DO3
        ENDDO DO3
        I2 = I2 - 1
        CALL H5ESET_AUTO_F( 0, SAERR )

!       PRINT *, " P = ", NAME(I1:I2)
        CALL H5GOPEN_F( L1_ID, NAME(I1:I2), L2_ID, HDFERR )
!       PRINT *, " HDFERR = ", HDFERR
        IF( HDFERR .NE. 0 ) &
          CALL H5GCREATE_F( L1_ID, NAME(I1:I2), L2_ID, HDFERR )
!       PRINT *, " HDFERR = ", HDFERR
        CALL H5ESET_AUTO_F( 1, SAERR )
        IF( HDFERR .NE. 0 ) RETURN
        IF( L1_ID .NE. L_ID ) &
          CALL H5GCLOSE_F( L1_ID, HDFERR )
        I1    = I2 + 1
        L1_ID = L2_ID
        IF( I1 .GT. N ) EXIT DO1
      ENDDO DO1

      G_ID = L2_ID


      END SUBROUTINE

!> @brief Create an HDF5 group and optionally a dataset
!! @param[in]  LOC_ID top HDF5 identification group
!! @param[in]  NAME group or datatset name
!! @param[out] DSET_ID output dataset id
!! @param[out] HDFERR return code (0 = OK)
!! @param[in]  DIMS dataset dimensions (so name is a dataset), optional
!! @param[in]  CRP_LIST dataset property, optional
      RECURSIVE SUBROUTINE MKDSR( LOC_ID, NAME, DSET_ID, HDFERR, DIMS, CRP_LIST)
!
      INTEGER(HID_T),   INTENT(IN)  :: LOC_ID
      CHARACTER(LEN=*), INTENT(IN)  :: NAME
      INTEGER(HID_T),   INTENT(OUT) :: DSET_ID  
      INTEGER,          INTENT(OUT) :: HDFERR
      INTEGER(HSIZE_T), INTENT(IN), OPTIONAL :: DIMS(:)
      INTEGER(HID_T),   INTENT(IN), OPTIONAL :: CRP_LIST

!

      INTEGER :: I1
      INTEGER :: SAERR
      INTEGER(HID_T) :: SPCE_ID
      INTEGER(HID_T) :: LOC1_ID

      IF( NAME(1:1) .EQ. '/' ) THEN
        CALL H5GOPEN_F( LOC_ID, '/', LOC1_ID, HDFERR )
        IF( HDFERR .NE. 0 ) RETURN
        CALL MKDSR( LOC1_ID, TRIM(NAME(2:)), DSET_ID, HDFERR, DIMS, CRP_LIST )
        RETURN
      ENDIF

      DO I1 = 1, LEN( TRIM( NAME ) )
        IF( NAME(I1:I1) .EQ. '/' ) EXIT
      ENDDO

      IF( I1 .LT. LEN( TRIM( NAME ) ) ) THEN
        CALL H5ESET_AUTO_F( 0, SAERR )
        CALL H5GOPEN_F( LOC_ID, NAME(1:I1-1), LOC1_ID, HDFERR )
        IF( HDFERR .NE. 0 ) &
          CALL H5GCREATE_F( LOC_ID, NAME(1:I1-1), LOC1_ID, HDFERR )
        CALL H5ESET_AUTO_F( 1, SAERR )
        IF( HDFERR .NE. 0 ) RETURN
        CALL MKDSR( LOC1_ID, TRIM(NAME(I1+1:)), DSET_ID, HDFERR, DIMS, CRP_LIST )
      ELSE
        CALL H5ESET_AUTO_F( 0, SAERR )
        CALL H5DOPEN_F( LOC_ID, NAME(1:I1-1), LOC1_ID, HDFERR )
        CALL H5ESET_AUTO_F( 1, SAERR )
        IF( HDFERR .NE. 0 ) THEN
          IF( PRESENT( DIMS ) ) THEN
            CALL H5SCREATE_SIMPLE_F( SIZE(DIMS), DIMS, SPCE_ID, HDFERR )
            IF( HDFERR .NE. 0 ) RETURN
            CALL H5DCREATE_F( LOC_ID, NAME(1:I1-1), H5T_FILE_REAL, SPCE_ID, LOC1_ID, HDFERR, CRP_LIST)
          ELSE
            HDFERR = 1
          ENDIF
        ENDIF
        IF( HDFERR .NE. 0 ) RETURN
        DSET_ID = LOC1_ID
      ENDIF

      END SUBROUTINE

!> @brief List to standard output the opened objects connected to FILE_ID
!! @param[out] ERR status on exit
!! @param[in]  FILE_ID HDF5 file identification
      SUBROUTINE RTTOV_HDF_LISTOBJ( ERR, FILE_ID)

      IMPLICIT NONE
      INTEGER(HID_T), INTENT(IN)       :: FILE_ID
      INTEGER(KIND=JPIM),  INTENT(OUT) :: ERR

      INTEGER(HID_T), ALLOCATABLE ::obj_ids(:)
      integer(SIZE_T)   :: maxobjs
      INTEGER(SIZE_T)   :: obj_count ! Number of opened objects
      INTEGER           ::obj_id
      character(len=64) :: obj_name
      INTEGER(SIZE_T)   :: buf_size =64  ! Buffer size
      INTEGER(SIZE_T)   :: name_size  ! Name size
!
TRY

      call h5fget_obj_count_f(file_id, H5F_OBJ_ALL_F, obj_count, err)
      THROWM(ERR.NE.0,"h5fget_obj_count_f")

      write(0,*) "Number of opened HDF5 objects found ",obj_count

      ALLOCATE(obj_ids(int(obj_count)), stat=err)
      THROWM(ERR.NE.0,"CANNOT ALLOCATE OBJ_IDS ")
     
      maxobjs = int(obj_count)

      obj_ids(:) = 0
      call h5fget_obj_ids_f(file_id, H5F_OBJ_ALL_F, maxobjs, obj_ids, err)
      THROWM(ERR.NE.0,"fget_obj_ids_f")

      do obj_id = 1, int(maxobjs)
        CALL h5iget_name_f(obj_ids(obj_id), obj_name, buf_size, name_size, err) 
        THROWM(ERR.NE.0,"h5iget_name_f")
        write(0,*) "     object id ", obj_ids(obj_id), "  is :", TRIM(obj_name)
      end do
      
      DEALLOCATE(obj_ids)
CATCH

      END SUBROUTINE

!> @brief Prepare HDF5 compression
!! @param[inout] LCOMPRESS internal compression
!!               this is input/output. True output means valid compression
!! @param[in]    DIMS dimensions of array
!! @param[out]   CHUNK_DIMS compression chunks
!! @param[out]   CHUNK_RANK rank of compression chunks
!! @param[out]   ERR status on exit
!! @param[in]    COMPDIMS dimensions to be addressed by compression
!!               is used for definition of compression chunks, optional
!! @param[in]    MINSA minimum number of element to allow compression, optional
      SUBROUTINE RTTOV_HDF_INIT_GZIP( LCOMPRESS, DIMS, &
         &  CHUNK_DIMS, CHUNK_RANK, ERR, COMPDIMS, MINSA)
      
      LOGICAL,           INTENT(INOUT) :: LCOMPRESS
      INTEGER(HSIZE_T),  INTENT(IN)    :: DIMS(7)
      INTEGER(HSIZE_T),  INTENT(OUT)   :: CHUNK_DIMS(7)
      INTEGER(KIND=JPIM),INTENT(OUT)   :: CHUNK_RANK
      INTEGER(KIND=JPIM),INTENT(OUT)   :: ERR

      INTEGER(KIND=JPIM),INTENT(IN), OPTIONAL :: COMPDIMS(:)
      INTEGER(KIND=JPIM),INTENT(IN), OPTIONAL :: MINSA

      INTEGER :: RANK

      INTEGER(KIND=JPIM) :: I

      LOGICAL :: GZIP_FLAG, SHUFFLE_FLAG
      INTEGER(KIND=JPIM) :: CONFIG_FLAG, CONFIG_FLAG_S, CONFIG_FLAG_BOTH
      INTEGER(KIND=JPIM) :: SA, MINSAVAL
TRY
      IF( .NOT. LCOMPRESS ) RETURN

      CHUNK_DIMS = DIMS

      IF(  PRESENT(MINSA) ) THEN
        IF( MINSA .GT. 0_JPIM ) THEN
          MINSAVAL = MINSA
        ELSE
          MINSAVAL = HUGE(MINSAVAL)
        ENDIF
      ELSE
        MINSAVAL = MINSAVALP
      ENDIF

      DO I = 7, 1, -1
        IF( DIMS(I) .GT. 0 ) EXIT
      ENDDO
      RANK = I

      IF(  PRESENT(COMPDIMS) ) THEN
        DO I = 1, SIZE(COMPDIMS)
          IF( COMPDIMS(I) .GT. 0) THEN
            CHUNK_DIMS(I) = COMPDIMS(I)
          ENDIF
        ENDDO
      ENDIF

      CHUNK_RANK = RANK

      SA = 1_JPIM
      DO I=1,RANK
        IF( CHUNK_DIMS(I) .GT. DIMS(I) ) THEN
          CHUNK_DIMS(I) = DIMS(I)
        ENDIF
        SA = SA * INT(CHUNK_DIMS(I))
      ENDDO

      CALL H5ZFILTER_AVAIL_F(H5Z_FILTER_DEFLATE_F, GZIP_FLAG, ERR)
      THROWM(ERR.NE.0,"h5zfilter_avail_f gzip")
      
      IF (.NOT. GZIP_FLAG) THEN
        LCOMPRESS = .FALSE.
      ENDIF

      CALL H5ZGET_FILTER_INFO_F(H5Z_FILTER_DEFLATE_F, CONFIG_FLAG, ERR)
      THROWM(ERR.NE.0,"h5zget_filter_info_f gzip")

      ! SHUFFLE Filter
      IF( ENABLE_SHUFFLE ) THEN
        CALL H5ZFILTER_AVAIL_F(H5Z_FILTER_SHUFFLE_F, SHUFFLE_FLAG, ERR)
        THROWM(ERR.NE.0,"h5zfilter_avail_f shuffle unavailable")

        CALL H5ZGET_FILTER_INFO_F(H5Z_FILTER_SHUFFLE_F, CONFIG_FLAG_S, ERR)
        THROWM(ERR.NE.0,"h5zget_filter_info_f shuffle")
      ENDIF
      !
      ! Make sure h5zget_filter_info_f returns the right flag
      !
      CONFIG_FLAG_BOTH=IOR(H5Z_FILTER_ENCODE_ENABLED_F,H5Z_FILTER_DECODE_ENABLED_F)
      IF( GZIP_FLAG ) THEN
        IF (CONFIG_FLAG .NE. CONFIG_FLAG_BOTH) THEN
          IF(CONFIG_FLAG .NE. H5Z_FILTER_DECODE_ENABLED_F)  THEN
            LCOMPRESS = .FALSE.
          ENDIF
        ENDIF
      ENDIF

      ! Continue only when encoder is available
      IF ( IAND(CONFIG_FLAG,  H5Z_FILTER_ENCODE_ENABLED_F) .EQ. 0 ) THEN
        LCOMPRESS = .FALSE.
      ENDIF

      ! Do not compress arrays less than MINSAVal elements
      IF(SA.LT.MINSAVAL) LCOMPRESS = .FALSE.
CATCH
      END SUBROUTINE

END MODULE
