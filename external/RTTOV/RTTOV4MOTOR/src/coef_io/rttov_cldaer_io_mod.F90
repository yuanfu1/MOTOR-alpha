! Description:
!> @file
!!   Subroutines for cloud/aerosol coefficient file I/O.
!
!> @brief
!!   Subroutines for cloud/aerosol coefficient file I/O.
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
!    Copyright 2018, EUMETSAT, All Rights Reserved.
!
MODULE rttov_cldaer_io_mod

  USE parkind1, ONLY : jprb, jpim, jplm

#include "throw.h"

  IMPLICIT NONE

#include "rttov_errorreport.interface"
#include "rttov_skipcommentline.interface"
#include "rttov_nullify_optp_data.interface"
#include "rttov_alloc_phfn_int.interface"
#include "rttov_channel_extract_sublist.interface"

  PRIVATE
  PUBLIC :: read_ascii_optp, write_ascii_optp, &
            read_binary_optp, write_binary_optp, &
            channel_extract_optp

CONTAINS

  !> Read ASCII cld/aer optical property header and data
  SUBROUTINE read_ascii_optp(err, file_id, section, coef, optp, channels)

    USE rttov_types, ONLY : rttov_coef, rttov_optp
    USE rttov_const, ONLY : scaercld_version_compatible_min, &
                            scaercld_version_compatible_max

    INTEGER(jpim),    INTENT(OUT)             :: err
    INTEGER(jpim),    INTENT(IN)              :: file_id
    CHARACTER(LEN=*), INTENT(IN)              :: section
    TYPE(rttov_coef), INTENT(IN)              :: coef
    TYPE(rttov_optp), INTENT(INOUT)           :: optp
    INTEGER(jpim),    INTENT(IN),    OPTIONAL :: channels(:)

    LOGICAL(jplm)              :: all_channels
    INTEGER(jpim)              :: i, j, k, n, irh, ideff
    INTEGER(jpim)              :: nchan_file, nchan_pha_file
    CHARACTER(LEN=32)          :: comp_name
    INTEGER(jpim), ALLOCATABLE :: list_of_channels(:)
    INTEGER(jpim), ALLOCATABLE :: phase_channels(:)  ! Solar channel numbers in original file
    INTEGER(jpim), ALLOCATABLE :: phase_ext_index(:) ! Indexes of phase data to extract
    REAL(jprb),    POINTER     :: abs_array(:,:,:)
    REAL(jprb),    POINTER     :: sca_array(:,:,:)
    REAL(jprb),    POINTER     :: bpr_array(:,:,:)
    INTEGER(jpim), POINTER     :: nmom_array(:,:)
    REAL(jprb),    POINTER     :: legcoef_array(:,:,:,:)
    REAL(jprb),    POINTER     :: pha_array(:,:,:,:)

    TRY

    all_channels = .NOT. PRESENT(channels)

    ! Header section

    CALL rttov_skipcommentline(file_id, err)
    THROWM(err.NE.0, "io status while reading section "//section)

    ! File format version
    READ(file_id, *, iostat=err) optp%version
    THROWM(err.NE.0, "io status while reading section "//section)

    IF (optp%version < scaercld_version_compatible_min .OR. &
        optp%version > scaercld_version_compatible_max) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, "Version of sccld/scaer file ("//TRIM(section)//") is incompatible with RTTOV library")
    ENDIF

    ! Contents ID
    READ(file_id, *, iostat=err) optp%id
    THROWM(err.NE.0, "io status while reading section "//section)

    ! Number of channels for which optical parameters are stored
    READ(file_id, *, iostat=err) optp%nchan
    THROWM(err.NE.0, "io status while reading section "//section)

    IF (coef%fmv_ori_nchn /= optp%nchan) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, "Incompatible channels between rtcoef and optical property files")
    ENDIF

    IF (.NOT. all_channels) THEN
      ALLOCATE (list_of_channels(SIZE(channels)))
      list_of_channels = channels
    ELSE
      ALLOCATE (list_of_channels(coef%fmv_chn))
      list_of_channels = (/(i, i = 1, coef%fmv_chn)/)
    ENDIF

    ! Take care of the user list of channels
    ! nchan_file is the number of channels in the file
    ! optp%nchan is the number of channels that the user requests
    nchan_file = optp%nchan
    IF (.NOT. all_channels) optp%nchan = SIZE(channels)

    ! Number of channels for which phase functions are stored
    READ(file_id, *, iostat=err) nchan_pha_file
    THROWM(err.NE.0, "io status while reading section "//section)

    ! Max number of Leg. coefs stored for each phase function
    READ(file_id, *, iostat=err) optp%maxnmom
    THROWM(err.NE.0, "io status while reading section "//section)

    ! Sort out the solar channels/phase functions
    IF (nchan_pha_file > 0) THEN

      ! Read the solar channel numbers
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0, "io status while reading section "//section)
      ALLOCATE(phase_channels(nchan_pha_file))
      READ(file_id, *, iostat=err) phase_channels(:)
      THROWM(err.NE.0, "io status while reading section "//section)

      ! Determine solar channels/phase functions to be extracted
      ALLOCATE (phase_ext_index(nchan_pha_file))
      CALL rttov_channel_extract_sublist( &
            err,                 &
            phase_channels,      &
            list_of_channels,    &
            optp%nchan_pha,      &
            optp%chan_pha,       &
            optp%chan_pha_index, &
            phase_ext_index)
      THROW(err.NE.0)

      ! At this point:
      !   nchan_pha_file         = total number of solar channels in file
      !   phase_channels(:)      = list of solar channel numbers
      !   optp%nchan_pha         = number of solar chans being extracted
      !   optp%chan_pha(:)       = list of solar channels to extract from file (indexes into extracted
      !                            chan list, NOT original channel numbers)
      !   phase_ext_index(:)     = list of indexes into the phase fns which are to be extracted to pha
      !   optp%chan_pha_index(:) = indexes for each extracted channel into the pha array

      ! Read phase function angle data
      READ(file_id, *, iostat=err) optp%nphangle
      THROWM(err.NE.0, "io status while reading section "//section)
      ALLOCATE(optp%phangle(optp%nphangle), STAT = err)
      THROWM(err.NE.0, "allocation of optp%phangle")
      READ(file_id, *, iostat=err) optp%phangle(:)
      THROWM(err.NE.0, "io status while reading section "//section)

      CALL rttov_alloc_phfn_int(err, optp%phangle, optp%phfn_int, 1_jpim)
      THROWM(err.NE.0, "initialisation of optp%phfn_int")
    ELSE
      optp%nchan_pha = 0
      optp%nphangle = 0
    ENDIF

    ! Number of particle types
    READ(file_id, *, iostat=err) optp%ntypes
    THROWM(err.NE.0, "io status while reading section "//section)

    ALLOCATE (optp%data(optp%ntypes), STAT = err)
    THROWM(err.NE.0, "allocation of optp%data")

    DO n = 1, optp%ntypes
      CALL rttov_nullify_optp_data(optp%data(n))

      IF (optp%ntypes > 1) THEN
        ! Particle name
        READ (file_id, '(a5)', iostat=err) comp_name
        THROWM(err.NE.0, "io status while reading section "//section)
        optp%data(n)%name = TRIM(ADJUSTL(comp_name))
      ENDIF

      ! Relative humidities
      READ(file_id, *, iostat=err) optp%data(n)%nrelhum
      THROWM(err.NE.0, "io status while reading section "//section)
      ALLOCATE (optp%data(n)%relhum(optp%data(n)%nrelhum), STAT = err)
      THROWM(err.NE.0, "allocation of optp%data(n)%relhum")
      READ(file_id, *, iostat=err) (optp%data(n)%relhum(i), i = 1, optp%data(n)%nrelhum)
      THROWM(err.NE.0, "io status while reading section "//section)

      ! Conversion factors
      READ(file_id, *, iostat=err) optp%data(n)%confac
      THROWM(err.NE.0, "io status while reading section "//section)

      ! Effective diameters
      READ(file_id, *, iostat=err) optp%data(n)%ndeff
      THROWM(err.NE.0, "io status while reading section "//section)
      ALLOCATE(optp%data(n)%deff(optp%data(n)%ndeff), STAT=err)
      THROWM(err.NE.0, "allocation of optp%data(n)%deff")
      READ(file_id, *, iostat=err) optp%data(n)%deff(:)
      THROWM(err.NE.0, "io status while reading section "//section)
    ENDDO

    DEALLOCATE (list_of_channels)

    ! Data section

    CALL rttov_skipcommentline(file_id, err)
    THROWM(err.NE.0, "io status while reading section "//section)

    DO n = 1, optp%ntypes
      ALLOCATE (optp%data(n)%abs(optp%data(n)%nrelhum, optp%data(n)%ndeff, optp%nchan), STAT = err)
      THROWM(err.NE.0, "allocation of optp%data(n)%abs)")

      ALLOCATE (optp%data(n)%sca(optp%data(n)%nrelhum, optp%data(n)%ndeff, optp%nchan), STAT = err)
      THROWM(err.NE.0, "allocation of optp%data(n)%sca")

      ALLOCATE (optp%data(n)%bpr(optp%data(n)%nrelhum, optp%data(n)%ndeff, optp%nchan), STAT = err)
      THROWM(err.NE.0, "allocation of optp%data(n)%bpr")

      IF (optp%maxnmom > 0) THEN
        ALLOCATE (optp%data(n)%nmom(optp%data(n)%nrelhum, optp%nchan), STAT = err)
        THROWM(err.NE.0, "allocation of optp%data(n)%nmom")

        ALLOCATE (optp%data(n)%legcoef(1:optp%maxnmom+1, optp%data(n)%nrelhum, optp%data(n)%ndeff, optp%nchan), STAT = err)
        THROWM(err.NE.0, "allocation of optp%data(n)%legcoef")
      ENDIF

      IF (optp%nchan_pha > 0) THEN
        ALLOCATE (optp%data(n)%pha(optp%nphangle, optp%data(n)%nrelhum, optp%data(n)%ndeff, optp%nchan_pha), STAT = err)
        THROWM(err.NE.0, "allocation of optp%data(n)%pha")
      ELSE
        NULLIFY(optp%data(n)%pha)
      ENDIF

      IF (all_channels) THEN
        abs_array     => optp%data(n)%abs
        sca_array     => optp%data(n)%sca
        bpr_array     => optp%data(n)%bpr
        nmom_array    => optp%data(n)%nmom
        legcoef_array => optp%data(n)%legcoef
        pha_array     => optp%data(n)%pha
      ELSE
        ALLOCATE (abs_array(optp%data(n)%nrelhum, optp%data(n)%ndeff, nchan_file), STAT = err)
        THROWM(err.NE.0, "allocation of abs_array")

        ALLOCATE (sca_array(optp%data(n)%nrelhum, optp%data(n)%ndeff, nchan_file), STAT = err)
        THROWM(err.NE.0, "allocation of sca_array")

        ALLOCATE (bpr_array(optp%data(n)%nrelhum, optp%data(n)%ndeff, nchan_file), STAT = err)
        THROWM(err.NE.0, "allocation of bpr_array")

        IF (optp%maxnmom > 0) THEN
          ALLOCATE (nmom_array(optp%data(n)%nrelhum, nchan_file), STAT = err)
          THROWM(err.NE.0, "allocation of nmom_array")

          ALLOCATE (legcoef_array(1:optp%maxnmom+1, optp%data(n)%nrelhum, optp%data(n)%ndeff, nchan_file), STAT = err)
          THROWM(err.NE.0, "allocation of legcoef_array")
        ENDIF

        IF (nchan_pha_file > 0) THEN
          ALLOCATE (pha_array(optp%nphangle, optp%data(n)%nrelhum, optp%data(n)%ndeff, nchan_pha_file), STAT = err)
          THROWM(err.NE.0, "allocation of pha_array")
        ENDIF
      ENDIF

      DO irh = 1, optp%data(n)%nrelhum
        CALL rttov_skipcommentline(file_id, err)
        THROWM(err.NE.0, "io status while reading section "//section)

        IF (optp%ntypes > 1) THEN
          READ(file_id, *, iostat=err) comp_name
          THROWM(err.NE.0, "io status while reading section "//section)
        ENDIF

        READ(file_id, *, iostat=err) abs_array(irh,:,:)
        THROWM(err.NE.0, "io status while reading section "//section)

        READ(file_id, *, iostat=err) sca_array(irh,:,:)
        THROWM(err.NE.0, "io status while reading section "//section)

        READ(file_id, *, iostat=err) bpr_array(irh,:,:)
        THROWM(err.NE.0, "io status while reading section "//section)

        k = 1
        DO i = 1, nchan_file
          IF (optp%maxnmom > 0) THEN
            READ(file_id, *, iostat=err) nmom_array(irh,i)
            THROWM(err.NE.0, "io status while reading nmom_array section "//section)

            legcoef_array(:,irh,:,i) = 0
            DO ideff = 1, optp%data(n)%ndeff
              READ(file_id, *, iostat=err) (legcoef_array(j,irh,ideff,i), j = 1, nmom_array(irh,i) + 1)
              THROWM(err.NE.0, "io status while reading legcoef_array section "//section)
            ENDDO
          ENDIF

          IF (nchan_pha_file > 0 .AND. k <= nchan_pha_file) THEN
            IF (phase_channels(k) == i) THEN
              DO ideff = 1, optp%data(n)%ndeff
                READ(file_id, *, iostat=err) (pha_array(j,irh,ideff,k), j = 1, optp%nphangle)
                THROWM(err.NE.0, "io status while reading pha_array section "//section)
              ENDDO
              k = k + 1
            ENDIF
          ENDIF
        ENDDO
      ENDDO

      IF (.NOT. all_channels) THEN
        optp%data(n)%abs(:,:,:)       = abs_array(:,:,channels(:))
        optp%data(n)%sca(:,:,:)       = sca_array(:,:,channels(:))
        optp%data(n)%bpr(:,:,:)       = bpr_array(:,:,channels(:))
        IF (optp%maxnmom > 0) THEN
          optp%data(n)%nmom(:,:)        = nmom_array(:,channels(:))
          optp%data(n)%legcoef(:,:,:,:) = legcoef_array(:,:,:,channels(:))
        ENDIF
        IF (optp%nchan_pha > 0) THEN
          optp%data(n)%pha(:,:,:,:) = pha_array(:,:,:,phase_ext_index(1:optp%nchan_pha))
        ENDIF

        DEALLOCATE (abs_array, sca_array, bpr_array, STAT = err)
        THROWM(err.NE.0, "deallocation of temporary arrays")

        IF (optp%maxnmom > 0) THEN
          DEALLOCATE (nmom_array, legcoef_array, STAT = err)
          THROWM(err.NE.0, "deallocation of nmom, legcoef arrays")
        ENDIF

        IF (nchan_pha_file > 0) THEN
          DEALLOCATE (pha_array, STAT = err)
          THROWM(err.NE.0, "deallocation of pha_array")
        ENDIF
      ENDIF
    ENDDO

    IF (ALLOCATED(phase_ext_index)) DEALLOCATE(phase_ext_index)
    IF (ALLOCATED(phase_channels))  DEALLOCATE(phase_channels)

    CATCH
  END SUBROUTINE read_ascii_optp

  !> Write ASCII cld/aer optical property header and data
  SUBROUTINE write_ascii_optp(err, file_id, coef, optp, ptype_str, routine_name)

    USE rttov_types, ONLY : rttov_coef, rttov_optp

    INTEGER(jpim),    INTENT(OUT) :: err
    INTEGER(jpim),    INTENT(IN)  :: file_id
    TYPE(rttov_coef), INTENT(IN)  :: coef
    TYPE(rttov_optp), INTENT(IN)  :: optp
    CHARACTER(LEN=*), INTENT(IN)  :: ptype_str
    CHARACTER(LEN=*), INTENT(IN)  :: routine_name

    INTEGER(jpim) :: i, n, irh, ideff, nchan_pha

    TRY

    ! Ensure we don't write out phase functions unnecessarily
    nchan_pha = 0
    IF (ASSOCIATED(coef%ss_val_chn)) THEN
      IF (ANY(coef%ss_val_chn(:) /= 0)) THEN
        nchan_pha = optp%nchan_pha
      ENDIF
    ENDIF

    ! Header section

    WRITE(file_id,'(1x,i8,t20,a)', iostat=err) &
      optp%version, '! File format version'
    THROW(err.NE.0)

    WRITE(file_id,'(1x,i8,t20,a)', iostat=err) &
      optp%id, '! Contents ID'
    THROW(err.NE.0)

    WRITE(file_id,'(1x,i8,t20,a)', iostat=err) &
      optp%nchan, '! Number of channels for which optical parameters are stored'
    THROW(err.NE.0)

    WRITE(file_id,'(1x,i8,t20,a)', iostat=err) &
      nchan_pha, '! Number of channels for which phase functions are stored'
    THROW(err.NE.0)

    WRITE(file_id,'(1x,i8,t20,a)', iostat=err) &
      optp%maxnmom, '! Maximum number of Legendre coefficients'
    THROW(err.NE.0)

    IF (nchan_pha > 0) THEN
      WRITE(file_id,'(1x,a)', iostat=err) '! Channel list for which phase functions are stored'
      THROW(err.NE.0)
      WRITE(file_id,'(10i6)', iostat=err) optp%chan_pha(1:nchan_pha)
      THROW(err.NE.0)

      WRITE(file_id,'(1x,i8,t20,a)', iostat=err) optp%nphangle, '! Phase function angles'
      THROW(err.NE.0)
      WRITE (file_id, '(10f8.3)', iostat=err) optp%phangle
      THROW(err.NE.0)
    ENDIF

    ! Number of particle types
    WRITE(file_id,'(1x,i8,t20,a)') optp%ntypes, '! Number of '//TRIM(ptype_str)//' types'

    DO n = 1, optp%ntypes

      IF (optp%ntypes > 1) THEN
        ! Particle name
        IF (LEN_TRIM(optp%data(n)%name) > 0) THEN
          WRITE (file_id, '(a5)', iostat=err) TRIM(optp%data(n)%name)
        ELSE
          WRITE (file_id, '(a)', iostat=err)' '//TRIM(ptype_str)//' ! default name for '//TRIM(routine_name)
        ENDIF
        THROW(err.NE.0)
      ENDIF

      ! Relative humidities
      WRITE (file_id, '(1x,i8,t20,a)', iostat=err) &
        optp%data(n)%nrelhum, '! Relative humidity values (%)'
      THROW(err.NE.0)
      WRITE (file_id, '(10f7.2)', iostat=err) optp%data(n)%relhum
      THROW(err.NE.0)

      ! Conversion factors
      WRITE (file_id, '(1x,e16.8,t20,a)', iostat=err) &
        optp%data(n)%confac, '! Conversion factor to particle density [part.cm^-3]/[g.m^-3]'
      THROW(err.NE.0)

      ! Effective diameters
      WRITE (file_id, '(1x,i8,t20,a)', iostat=err) &
        optp%data(n)%ndeff, '! Effective diameters (microns)'
      THROW(err.NE.0)
      WRITE (file_id, '(6e14.6)', iostat=err) optp%data(n)%deff
      THROW(err.NE.0)
    ENDDO

    ! Data section

    WRITE (file_id, '(a)', iostat=err) ' ! ------------------------------------------------------'
    THROW(err.NE.0)
    WRITE (file_id, '(a)', iostat=err) ' ! OPTICAL PROPERTY DATA'
    THROW(err.NE.0)

    DO n = 1, optp%ntypes
      DO irh = 1, optp%data(n)%nrelhum

        IF (optp%ntypes > 1) THEN
          WRITE (file_id, '(a)', iostat=err) ' ! ---------------------'
          THROW(err.NE.0)

          IF (LEN_TRIM(optp%data(n)%name) > 0) THEN
            WRITE (file_id, '(a5,i2.2)', iostat=err) TRIM(optp%data(n)%name), &
                                                     INT(optp%data(n)%relhum(irh))
          ELSE
            WRITE (file_id, '(a)', iostat=err) ' '//TRIM(ptype_str)//' ! default name for '//TRIM(routine_name)
          ENDIF
          THROW(err.NE.0)
        ENDIF

        WRITE (file_id, '(5e16.8)', iostat=err) optp%data(n)%abs(irh,:,:)
        THROW(err.NE.0)

        WRITE (file_id, '(5e16.8)', iostat=err) optp%data(n)%sca(irh,:,:)
        THROW(err.NE.0)

        WRITE (file_id, '(5e16.8)', iostat=err) optp%data(n)%bpr(irh,:,:)
        THROW(err.NE.0)

        DO i = 1, coef%fmv_chn
          IF (optp%maxnmom > 0) THEN
            WRITE (file_id, '(i6)', iostat=err) optp%data(n)%nmom(irh,i)
            THROW(err.NE.0)

            DO ideff = 1, optp%data(n)%ndeff
              WRITE (file_id, '(5e16.8)', iostat=err) &
                optp%data(n)%legcoef(1:optp%data(n)%nmom(irh,i)+1,irh,ideff,i)
              THROW(err.NE.0)
            ENDDO
          ENDIF

          IF (nchan_pha > 0) THEN
            IF (optp%chan_pha_index(i) > 0) THEN
              DO ideff = 1, optp%data(n)%ndeff
                WRITE (file_id, '(5e16.8)', iostat=err) optp%data(n)%pha(:,irh,ideff,optp%chan_pha_index(i))
                THROW(err.NE.0)
              ENDDO
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    CATCH
  END SUBROUTINE write_ascii_optp


  !> Read binary cld/aer optical property header and data
  SUBROUTINE read_binary_optp(err, file_id, coef, optp, channels)

    USE rttov_types, ONLY : rttov_coef, rttov_optp
    USE rttov_const, ONLY : scaercld_version_compatible_min, &
                            scaercld_version_compatible_max

    INTEGER(jpim),    INTENT(OUT)                :: err
    INTEGER(jpim),    INTENT(IN)                 :: file_id
    TYPE(rttov_coef), INTENT(IN)                 :: coef
    TYPE(rttov_optp), INTENT(INOUT)              :: optp
    INTEGER(jpim),    INTENT(IN),    OPTIONAL    :: channels(:)

    LOGICAL(jplm)              :: all_channels
    INTEGER(jpim)              :: i, n, irh, nchan_file, nchan_pha_file
    INTEGER(jpim), ALLOCATABLE :: list_of_channels(:)
    INTEGER(jpim), ALLOCATABLE :: phase_channels(:)  ! Solar channel numbers in original file
    INTEGER(jpim), ALLOCATABLE :: phase_ext_index(:) ! Indexes of phase data to extract
    REAL(jprb),    ALLOCATABLE :: dvalues0(:,:,:)
    REAL(jprb),    ALLOCATABLE :: tvalues0(:,:,:,:)
    INTEGER(jpim), ALLOCATABLE :: ivalues0(:,:)

    TRY

    all_channels = .NOT. PRESENT(channels)

    ! Header section

    READ (file_id, iostat=err) optp%version
    THROWM(err.NE.0,"reading file version")

    IF (optp%version < scaercld_version_compatible_min .OR. &
        optp%version > scaercld_version_compatible_max) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0, "Version of sccld/scaer file is incompatible with RTTOV library")
    ENDIF

    READ (file_id, iostat=err) optp%id
    THROWM(err.NE.0,"reading file id")

    READ (file_id, iostat=err) optp%nchan, nchan_pha_file, optp%ntypes, optp%maxnmom
    THROWM(err.NE.0,"reading dimensions")

    IF (coef%fmv_ori_nchn /= optp%nchan) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,"Incompatible channels between rtcoef and optical property files")
    ENDIF

    IF (.NOT. all_channels) THEN
      ALLOCATE (list_of_channels(SIZE(channels)))
      list_of_channels = channels
    ELSE
      ALLOCATE (list_of_channels(coef%fmv_chn))
      list_of_channels = (/(i, i = 1, coef%fmv_chn)/)
    ENDIF

    ! Take care of the user list of channels
    ! nchan_file is the number of channels in the file
    ! optp%nchan is the number of channels that the user requests
    nchan_file = optp%nchan
    IF (.NOT. all_channels) optp%nchan = SIZE(channels)

    ! Sort out the solar channels/phase functions
    IF (nchan_pha_file > 0) THEN

      ! Read the solar channel numbers
      ALLOCATE(phase_channels(nchan_pha_file))
      READ (file_id, iostat=err)phase_channels(:)
      THROWM(err.NE.0,"reading solar channel numbers")

      ! Determine solar channels/phase functions to be extracted
      ALLOCATE (phase_ext_index(nchan_pha_file))
      CALL rttov_channel_extract_sublist( &
            err,                 &
            phase_channels,      &
            list_of_channels,    &
            optp%nchan_pha,      &
            optp%chan_pha,       &
            optp%chan_pha_index, &
            phase_ext_index)
      THROW(err.NE.0)
      DEALLOCATE(phase_channels)

      ! See read_ascii_optp_hdr above for a description of what the various arrays contain.

      ! Read phase function angle data
      READ (file_id, iostat=err) optp%nphangle
      THROWM(err.NE.0,"reading nphangle")
      ALLOCATE(optp%phangle(optp%nphangle), STAT=err)
      THROWM(err.NE.0,"allocation of phangle")
      READ (file_id, iostat=err) optp%phangle(:)
      THROWM(err.NE.0,"reading phangle")

      CALL rttov_alloc_phfn_int(err, optp%phangle, optp%phfn_int, 1_jpim)
      THROWM(err.NE.0, "initialisation of optp%phfn_int")

    ELSE
      optp%nchan_pha = 0
      optp%nphangle = 0
    ENDIF

    ALLOCATE (optp%data(optp%ntypes), STAT = err)
    THROWM(err.NE.0, "allocation of optp%data")

    DO n = 1, optp%ntypes
      CALL rttov_nullify_optp_data(optp%data(n))
    ENDDO

    READ (file_id, iostat=err) optp%data(:)%name
    THROWM(err.NE.0,"reading optp%data(:)%name")

    READ (file_id, iostat=err) optp%data(:)%nrelhum
    THROWM(err.NE.0,"reading optp%data(:)%nrelhum")

    READ (file_id, iostat=err) optp%data(:)%ndeff
    THROWM(err.NE.0,"reading optp%data(:)%ndeff")

    READ (file_id, iostat=err) optp%data(:)%confac
    THROWM(err.NE.0,"reading optp%data(:)%confac")

    DO n = 1, optp%ntypes
      ALLOCATE (optp%data(n)%relhum(optp%data(n)%nrelhum), STAT=err)
      THROWM(err.NE.0,"allocation of optp%data(n)%relhum")
      READ (file_id, iostat=err) optp%data(n)%relhum
      THROWM(err.NE.0,"reading optp%data(n)%relhum")

      ALLOCATE (optp%data(n)%deff(optp%data(n)%ndeff), STAT=err)
      THROWM(err.NE.0,"allocation of optp%data(n)%deff")
      READ (file_id, iostat=err) optp%data(n)%deff
      THROWM(err.NE.0,"reading optp%data(n)%deff")
    ENDDO

    ! Data section

    DO n = 1, optp%ntypes
      ALLOCATE (optp%data(n)%abs(optp%data(n)%nrelhum, optp%data(n)%ndeff, optp%nchan), STAT=err)
      THROWM(err.NE.0,"allocation of optp%data(n)%abs")

      IF (all_channels) THEN
        READ (file_id, iostat=err)optp%data(n)%abs
        THROWM(err.NE.0,"reading optp%data(n)%abs")
      ELSE
        ALLOCATE (dvalues0(optp%data(n)%nrelhum, optp%data(n)%ndeff, nchan_file), STAT=err)
        THROWM(err.NE.0,"allocation of dvalues0")

        READ (file_id, iostat=err)dvalues0
        THROWM(err.NE.0,"reading dvalues0")

        optp%data(n)%abs(:,:,:) = dvalues0(:,:,channels(:))

        DEALLOCATE (dvalues0, STAT=err)
        THROWM(err.NE.0,"deallocation of dvalues0")
      ENDIF
    ENDDO

    DO n = 1, optp%ntypes
      ALLOCATE (optp%data(n)%sca(optp%data(n)%nrelhum, optp%data(n)%ndeff, optp%nchan), STAT=err)
      THROWM(err.NE.0,"allocation of optp%data(n)%sca")

      IF (all_channels) THEN
        READ (file_id, iostat=err)optp%data(n)%sca
        THROWM(err.NE.0,"reading optp%data(n)%sca")
      ELSE
        ALLOCATE (dvalues0(optp%data(n)%nrelhum, optp%data(n)%ndeff, nchan_file), STAT=err)
        THROWM(err.NE.0,"allocation of dvalues0")

        READ (file_id, iostat=err)dvalues0
        THROWM(err.NE.0,"reading dvalues0")

        optp%data(n)%sca(:,:,:) = dvalues0(:,:,channels(:))

        DEALLOCATE (dvalues0, STAT=err)
        THROWM(err.NE.0,"deallocation of dvalues0")
      ENDIF
    ENDDO

    DO n = 1, optp%ntypes
      ALLOCATE (optp%data(n)%bpr(optp%data(n)%nrelhum, optp%data(n)%ndeff, optp%nchan), STAT=err)
      THROWM(err.NE.0,"allocation of optp%data(n)%bpr")

      IF (all_channels) THEN
        READ (file_id, iostat=err)optp%data(n)%bpr
        THROWM(err.NE.0,"reading optp%data(n)%bpr")
      ELSE
        ALLOCATE (dvalues0(optp%data(n)%nrelhum, optp%data(n)%ndeff, nchan_file), STAT=err)
        THROWM(err.NE.0,"allocation of dvalues0")

        READ (file_id, iostat=err)dvalues0
        THROWM(err.NE.0,"reading dvalues0")

        optp%data(n)%bpr(:,:,:) = dvalues0(:,:,channels(:))

        DEALLOCATE (dvalues0, STAT=err)
        THROWM(err.NE.0,"deallocation of dvalues0")
      ENDIF
    ENDDO

    IF (optp%maxnmom > 0) THEN
      DO n = 1, optp%ntypes
        ALLOCATE (optp%data(n)%nmom(optp%data(n)%nrelhum,optp%nchan),STAT=err)
        THROWM(err.NE.0,"allocation of optp%data(n)%nmom")

        ALLOCATE (optp%data(n)%legcoef(1:optp%maxnmom+1,optp%data(n)%nrelhum,optp%data(n)%ndeff,optp%nchan), STAT=err)
        THROWM(err.NE.0,"allocation of optp%data(n)%legcoef")

        IF (all_channels) THEN
          DO irh = 1, optp%data(n)%nrelhum
            DO i = 1, optp%nchan
              READ (file_id, iostat=err) optp%data(n)%nmom(irh,i)
              THROWM(err.NE.0,"reading optp%data(n)%nmom")

              READ (file_id, iostat=err) optp%data(n)%legcoef(1:optp%data(n)%nmom(irh,i)+1,irh,:,i)
              THROWM(err.NE.0,"reading optp%data(n)%legcoef")
            ENDDO
          ENDDO
        ELSE
          ALLOCATE (ivalues0(optp%data(n)%nrelhum,nchan_file), &
                    tvalues0(0:optp%maxnmom,optp%data(n)%nrelhum,optp%data(n)%ndeff,nchan_file), STAT=err)
          THROWM(err.NE.0,"allocation of ivalues0, tvalues0")

          DO irh = 1, optp%data(n)%nrelhum
            DO i = 1, nchan_file
              READ (file_id, iostat=err) ivalues0(irh,i)
              THROWM(err.NE.0,"reading optp%data(n)%nmom")

              READ (file_id, iostat=err) tvalues0(0:ivalues0(irh,i),irh,:,i)
              THROWM(err.NE.0,"reading optp%data(n)%legcoef")
            ENDDO
          ENDDO

          optp%data(n)%nmom(:,:) = ivalues0(:,channels(:))
          optp%data(n)%legcoef(:,:,:,:) = tvalues0(:,:,:,channels(:))

          DEALLOCATE (ivalues0, tvalues0, STAT=err)
          THROWM(err.NE.0,"deallocation of ivalues0, tvalues0")
        ENDIF
      ENDDO
    ENDIF

    IF (optp%nchan_pha > 0) THEN
      DO n = 1, optp%ntypes
        ALLOCATE (optp%data(n)%pha(optp%nphangle,optp%data(n)%nrelhum,optp%data(n)%ndeff,optp%nchan_pha), STAT=err)
        THROWM(err.NE.0,"allocation of optp%data(n)%pha")

        IF (all_channels) THEN
          READ (file_id, iostat=err) optp%data(n)%pha(:,:,:,:)
          THROWM(err.NE.0,"reading optp%data(n)%pha")
        ELSE
          ALLOCATE (tvalues0(optp%nphangle,optp%data(n)%nrelhum,optp%data(n)%ndeff,nchan_pha_file), STAT=err)
          THROWM(err.NE.0,"allocation of tvalues0")

          READ (file_id, iostat=err) tvalues0(:,:,:,:)
          THROWM(err.NE.0,"reading optp%data(n)%pha")

          IF (optp%nchan_pha > 0) THEN
            optp%data(n)%pha(:,:,:,:) = tvalues0(:,:,:,phase_ext_index(1:optp%nchan_pha))
          ENDIF

          DEALLOCATE (tvalues0, STAT=err)
          THROWM(err.NE.0,"deallocation of tvalues0")
        ENDIF
      ENDDO
    ENDIF

    IF (ALLOCATED(phase_ext_index)) DEALLOCATE(phase_ext_index)
    DEALLOCATE (list_of_channels)

    CATCH
  END SUBROUTINE read_binary_optp

  !> Write binary cld/aer optical property header and data
  SUBROUTINE write_binary_optp(err, file_id, coef, optp)

    USE rttov_types, ONLY : rttov_coef, rttov_optp

    INTEGER(jpim),    INTENT(OUT)   :: err
    INTEGER(jpim),    INTENT(IN)    :: file_id
    TYPE(rttov_coef), INTENT(IN)    :: coef
    TYPE(rttov_optp), INTENT(IN)    :: optp

    INTEGER(jpim) :: i, n, irh, nchan_pha

    TRY

    ! Ensure we don't write out phase functions unnecessarily
    nchan_pha = 0
    IF (ASSOCIATED(coef%ss_val_chn)) THEN
      IF (ANY(coef%ss_val_chn(:) /= 0)) THEN
        nchan_pha = optp%nchan_pha
      ENDIF
    ENDIF

    ! Header section

    WRITE (file_id, iostat=err) optp%version
    THROW(err.NE.0)

    WRITE (file_id, iostat=err) optp%id
    THROW(err.NE.0)

    WRITE (file_id, iostat=err) optp%nchan, nchan_pha, optp%ntypes, optp%maxnmom
    THROW(err.NE.0)

    IF (nchan_pha > 0) THEN
      WRITE (file_id, iostat=err) optp%chan_pha(:)
      THROW(err.NE.0)

      WRITE (file_id, iostat=err) optp%nphangle
      THROW(err.NE.0)

      WRITE (file_id, iostat=err) optp%phangle(:)
      THROW(err.NE.0)
    ENDIF

    WRITE (file_id, iostat=err) optp%data(:)%name
    THROW(err.NE.0)

    WRITE (file_id, iostat=err) optp%data(:)%nrelhum
    THROW(err.NE.0)

    WRITE (file_id, iostat=err) optp%data(:)%ndeff
    THROW(err.NE.0)

    WRITE (file_id, iostat=err) optp%data(:)%confac
    THROW(err.NE.0)

    DO n = 1, optp%ntypes
      WRITE (file_id, iostat=err) optp%data(n)%relhum
      THROW(err.NE.0)

      WRITE (file_id, iostat=err) optp%data(n)%deff
      THROW(err.NE.0)
    ENDDO

    ! Data section

    DO n = 1, optp%ntypes
      WRITE (file_id, iostat=err) optp%data(n)%abs
      THROW(err.NE.0)
    ENDDO

    DO n = 1, optp%ntypes
      WRITE (file_id, iostat=err) optp%data(n)%sca
      THROW(err.NE.0)
    ENDDO

    DO n = 1, optp%ntypes
      WRITE (file_id, iostat=err) optp%data(n)%bpr
      THROW(err.NE.0)
    ENDDO

    IF (optp%maxnmom > 0) THEN
      DO n = 1, optp%ntypes
        DO irh = 1, optp%data(n)%nrelhum
          DO i = 1, optp%nchan
            WRITE (file_id, iostat=err) optp%data(n)%nmom(irh,i)
            THROW(err.NE.0)

            WRITE (file_id, iostat=err) optp%data(n)%legcoef(1:optp%data(n)%nmom(irh,i)+1,irh,:,i)
            THROW(err.NE.0)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    IF (nchan_pha > 0) THEN
      DO n = 1, optp%ntypes
        WRITE (file_id, iostat=err) optp%data(n)%pha
        THROW(err.NE.0)
      ENDDO
    ENDIF

    CATCH
  END SUBROUTINE write_binary_optp

  !> Extract a subset of channels from an aerosol/cloud optical property structure
  SUBROUTINE channel_extract_optp(err, optp1, optp2, channels)
    USE rttov_types, ONLY : rttov_optp

    INTEGER(jpim),    INTENT(OUT)   :: err
    TYPE(rttov_optp), INTENT(IN)    :: optp1
    TYPE(rttov_optp), INTENT(INOUT) :: optp2
    INTEGER(jpim),    INTENT(IN)    :: channels(:)

    INTEGER(jpim)              :: n
    INTEGER(jpim), ALLOCATABLE :: phase_channels(:)
    INTEGER(jpim), ALLOCATABLE :: phase_ext_index(:)

    TRY

    optp2%version = optp1%version
    optp2%id      = optp1%id
    optp2%ntypes  = optp1%ntypes
    optp2%maxnmom = optp1%maxnmom

    ! Work out the channel numbers including phase function (solar) channels
    optp2%nchan = SIZE(channels)

    IF (optp1%nchan_pha > 0) THEN
      ALLOCATE(phase_channels(optp1%nchan_pha))

      ! Determine the solar channel numbers
      phase_channels(:) = optp1%chan_pha(:)

      ! Determine solar channels/phase functions to be extracted
      ALLOCATE(phase_ext_index(optp1%nchan_pha))
      CALL rttov_channel_extract_sublist( &
            err,                  &
            phase_channels,       &
            channels,             &
            optp2%nchan_pha,      &
            optp2%chan_pha,       &
            optp2%chan_pha_index, &
            phase_ext_index)
      THROW(err.NE.0)

      ! See read_ascii_optp_hdr above for a description of what the various arrays contain.

      IF (optp2%nchan_pha > 0) THEN
        optp2%nphangle = optp1%nphangle
        ALLOCATE(optp2%phangle(optp2%nphangle))
        optp2%phangle = optp1%phangle

        CALL rttov_alloc_phfn_int(err, optp2%phangle, optp2%phfn_int, 1_jpim)
        THROWM(err.NE.0, "initialisation of optp2%phfn_int")
      ENDIF

      IF (ALLOCATED(phase_channels)) DEALLOCATE(phase_channels)
    ELSE
      ! What do we need to do if there are no phase fns at all?
      optp2%nchan_pha = 0
      optp2%nphangle = 0
    ENDIF

    ALLOCATE(optp2%data(optp2%ntypes), stat=err)
    THROWM(err.NE.0, 'allocation of optp2%data')

    DO n = 1, optp2%ntypes

      CALL rttov_nullify_optp_data(optp2%data(n))

      optp2%data(n)%name    = optp1%data(n)%name
      optp2%data(n)%confac  = optp1%data(n)%confac
      optp2%data(n)%nrelhum = optp1%data(n)%nrelhum
      optp2%data(n)%ndeff   = optp1%data(n)%ndeff

      IF (ASSOCIATED(optp1%data(n)%relhum)) THEN
        ALLOCATE(optp2%data(n)%relhum(optp2%data(n)%nrelhum), stat=err)
        THROWM(err.NE.0, 'allocation of optp2%data(n)%relhum')
        optp2%data(n)%relhum = optp1%data(n)%relhum
      ENDIF

      IF (ASSOCIATED(optp1%data(n)%deff)) THEN
        ALLOCATE(optp2%data(n)%deff(optp2%data(n)%ndeff), stat=err)
        THROWM(err.NE.0, 'allocation of optp2%data(n)%deff')
        optp2%data(n)%deff = optp1%data(n)%deff
      ENDIF

      ALLOCATE(optp2%data(n)%abs(optp2%data(n)%nrelhum,optp2%data(n)%ndeff,optp2%nchan), stat=err)
      THROWM(err.NE.0, 'allocation of data(n)%abs')
      optp2%data(n)%abs(:,:,:) = optp1%data(n)%abs(:,:,channels)

      ALLOCATE(optp2%data(n)%sca(optp2%data(n)%nrelhum,optp2%data(n)%ndeff,optp2%nchan), stat=err)
      THROWM(err.NE.0, 'allocation of data(n)%sca')
      optp2%data(n)%sca(:,:,:) = optp1%data(n)%sca(:,:,channels)

      ALLOCATE(optp2%data(n)%bpr(optp2%data(n)%nrelhum,optp2%data(n)%ndeff,optp2%nchan), stat=err)
      THROWM(err.NE.0, 'allocation of data(n)%bpr')
      optp2%data(n)%bpr(:,:,:) = optp1%data(n)%bpr(:,:,channels)

      IF (optp2%maxnmom > 0) THEN

        ALLOCATE(optp2%data(n)%nmom(optp2%data(n)%nrelhum,optp2%nchan), stat=err)
        THROWM(err.NE.0, 'allocation of data(n)%nmom')
        optp2%data(n)%nmom(:,:) = optp1%data(n)%nmom(:,channels)

        ALLOCATE(optp2%data(n)%legcoef(1:optp2%maxnmom+1, &
                                       optp2%data(n)%nrelhum, optp2%data(n)%ndeff, &
                                       optp2%nchan), stat=err)
        THROWM(err.NE.0, 'allocation of data(n)%legcoef')
        optp2%data(n)%legcoef(:,:,:,:) = optp1%data(n)%legcoef(1:optp2%maxnmom+1,:,:,channels)

      ENDIF

      IF (optp2%nchan_pha > 0) THEN

        ALLOCATE(optp2%data(n)%pha(optp2%nphangle, &
                                   optp2%data(n)%nrelhum, optp2%data(n)%ndeff, &
                                   optp2%nchan_pha), stat=err)
        THROWM(err.NE.0, 'allocation of data(n)%pha')
        optp2%data(n)%pha(:,:,:,:) = optp1%data(n)%pha(:,:,:,phase_ext_index(1:optp2%nchan_pha))

      ENDIF

    ENDDO

    IF (ALLOCATED(phase_ext_index)) DEALLOCATE (phase_ext_index)

    CATCH
  END SUBROUTINE channel_extract_optp

END MODULE rttov_cldaer_io_mod
