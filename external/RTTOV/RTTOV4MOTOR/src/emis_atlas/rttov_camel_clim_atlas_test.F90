! Description:
!> @file
!!   Executable for basic test of the CAMEL climatology IR emissivity atlas.
!
!> @brief
!!   Executable for basic test of the CAMEL climatology IR emissivity atlas.
!!
!! @details
!!   This should be run using the test_camel_clim_atlas.sh script in rttov_test/
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
program rttov_camel_clim_atlas_test

#include "throw.h"
  use parkind1, only : jpim, jprb, jplm

  use rttov_unix_env, only : rttov_exit

  use rttov_const, only : errorstatus_success

  use rttov_types, only : &
        rttov_coefs,      &
        rttov_options,    &
        rttov_profile,    &
        rttov_chanprof

  use mod_rttov_emis_atlas, only : atlas_type_ir, rttov_emis_atlas_data

  implicit none

  real(kind=jprb),    parameter :: difftol = 1.E-9_jprb

  integer(kind=jpim), parameter :: ioin  = 50 ! Input file unit
  integer(kind=jpim), parameter :: ioout = 51 ! Output file unit

  character(len=256) :: coef_filename ! Coefficient filename
  character(len=256) :: atlas_path    ! Path to atlas data
  integer(kind=jpim) :: imonth        ! Month for which to load data

  integer(kind=jpim) :: lo,hi,j,k
  integer(kind=jpim) :: err
  integer(kind=jpim) :: angcorr_run, init_run

  integer(kind=jpim) :: nprof  ! Number of profiles
  integer(kind=jpim) :: nchan  ! Number of channels per profile

  logical(kind=jplm) :: do_angcorr, single_instrument
  logical(kind=jplm) :: diffflag

  integer(kind=jpim), allocatable :: channels(:)

  type(rttov_coefs)                 :: coefs
  type(rttov_options)               :: opts
  type(rttov_profile),  allocatable :: profiles(:)
  type(rttov_chanprof), allocatable :: chanprof(:)

  type(rttov_emis_atlas_data)       :: atlas

  real(kind=jprb),    allocatable :: emissivity(:,:)
  real(kind=jprb),    allocatable :: emis_std(:,:)
  integer(kind=jpim), allocatable :: emis_flag(:,:)

#include "rttov_errorreport.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_skipcommentline.interface"

#include "rttov_setup_emis_atlas.interface"
#include "rttov_get_emis.interface"
#include "rttov_deallocate_emis_atlas.interface"

!-----------------------------------------------------------------------------------
TRY

! This test uses a simple input file for 'profiles' consisting of lat, lon, surface
! type and snow fraction. The controlling shell script passes the month, instrument,
! and channels for which to retrieve emissivities.

! Emissivity values, errors and the associated atlas flags are returned for all
! profiles and are written to an output file for comparison with reference output.

! The test is repeated for two cases, where the atlas is initialised for use with
! multiple instruments, and for use with a single instrument.

! Read profiles input file
  open(ioin, file='profiles_ir', status='old')
  call rttov_skipcommentline(ioin, err)
  THROWM(err.NE.0, "Error reading input profiles")

  read(ioin,*) nprof
  call rttov_skipcommentline(ioin, err)
  THROWM(err.NE.0, "Error reading input profiles")

  allocate(profiles(nprof))
  do k = 1, nprof
    read(ioin,*) profiles(k)%latitude,      &
                 profiles(k)%longitude,     &
                 profiles(k)%zenangle,      &
                 profiles(k)%sunzenangle,   &
                 profiles(k)%skin%surftype, &
                 profiles(k)%skin%snow_fraction
  enddo

  close(ioin)

! Read data from control script
  read(*,*) imonth
  read(*,'(a)') coef_filename
  read(*,*) nchan
  allocate(channels(nchan))
  read(*,*) channels(1:nchan)
  read(*,'(a)') atlas_path

! Generate the chanprof array
  allocate(chanprof(nchan*nprof))
  do k = 1, nprof
    lo = (k-1)*nchan+1
    hi = lo+nchan-1
    chanprof(lo:hi)%prof = k
    chanprof(lo:hi)%chan = (/ (j, j=1,nchan) /)
  enddo

  allocate(emissivity(nchan*nprof,2))
  allocate(emis_std(nchan*nprof,2))
  allocate(emis_flag(nchan*nprof,2))

! Set up RTTOV coefficients
  call rttov_read_coefs(                &
              err,                      &
              coefs,                    &
              opts,                     &
              channels = channels(:),   &
              file_coef = coef_filename)
  THROWM(err.NE.0, "Failure reading coefficients" )

! Open file for output
  open(ioout, file='output_camel_clim_atlas.ascii', form='formatted', status='replace', action='write')

! Test both ways of initialising the atlas, each with/without angular correction
  do angcorr_run = 1, 2

    if (angcorr_run == 1) then
      print *, 'Initialising atlas without angular correction'
      do_angcorr = .FALSE.   ! initialise without angular correction
    else
      print *, 'Initialising atlas with angular correction'
      do_angcorr = .TRUE.    ! initialise with angular correction
    endif

    do init_run = 1, 2

      if (init_run == 1) then
        print *, 'Initialising atlas for multiple instruments'
        single_instrument = .FALSE.   ! initialise for use with multiple instruments
      else
        print *, 'Initialising atlas for a single instrument'
        single_instrument = .TRUE.    ! initialise for use with a single instrument
      endif

!-------------------------
! Set up emissivity atlas
!-------------------------
      if (single_instrument) then
        call rttov_setup_emis_atlas(   &
                    err,               &
                    opts,              &
                    imonth,            &
                    atlas_type_ir,     &  ! IR atlas
                    atlas,             &
                    atlas_id = 3_jpim, &  ! To select CAMEL climatology atlas
                    path = atlas_path, &
                    coefs = coefs,     &
                    ir_atlas_read_std = .TRUE., &
                    ir_atlas_ang_corr = do_angcorr)
      else
        call rttov_setup_emis_atlas(   &
                    err,               &
                    opts,              &
                    imonth,            &
                    atlas_type_ir,     &  ! IR atlas
                    atlas,             &
                    atlas_id = 3_jpim, &  ! To select CAMEL climatology atlas
                    path = atlas_path, &
                    ir_atlas_read_std = .TRUE., &
                    ir_atlas_ang_corr = do_angcorr)
      endif
      THROWM(err.NE.0, "Failure in setting up emissivity atlas" )

!----------------------------
! Retrieve values from atlas
!----------------------------
      call rttov_get_emis(                           &
                  err,                               &
                  opts,                              &
                  chanprof,                          &
                  profiles,                          &
                  coefs,                             &
                  atlas,                             &
                  emissivity(:,init_run),            &
                  emis_std = emis_std(:,init_run),   &
                  emis_flag = emis_flag(:,init_run))
      THROWM(err.NE.0, "Failure in retrieving emissivity values" )

! Write out emissivity data
      write(ioout,'(a,l1)') 'Init for angular correction? ', do_angcorr
      write(ioout,'(a,l1)') 'Init for single inst? ', single_instrument
      do k = 1, nprof
        write(ioout,'(a,i2.2)') 'Profile ',k
        write(ioout,'(a)') ' Chan  Emissivity  Standard Dev  Flag'
        lo = (k-1)*nchan
        do j = 1, nchan
          write(ioout,'(i5,f11.4,f13.6,i6)') &
            channels(j),               &
            emissivity(lo+j,init_run), &
            emis_std(lo+j,init_run),   &
            emis_flag(lo+j,init_run)
        enddo
      enddo

!-----------------------
! Deallocate atlas data
!-----------------------
      call rttov_deallocate_emis_atlas(atlas)

    enddo ! init_run

! Ensure single and multi instrument runs gave same results
    diffflag = .false.
    if (any(abs(emissivity(:,1) - emissivity(:,2)) > difftol)) then
      print *, emissivity(:,1) - emissivity(:,2)
      write(ioout,'(a)') 'WARNING: difference in emissivity between single and multi-instrument results'
      diffflag = .true.
    endif
    if (any(abs(emis_std(:,1) - emis_std(:,2)) > difftol)) then
      print *,emis_std(:,1) - emis_std(:,2)
      write(ioout,'(a)') 'WARNING: difference in emissivity std between single and multi-instrument results'
      diffflag = .true.
    endif
    if (any(abs(emis_flag(:,1) - emis_flag(:,2)) > difftol)) then
      write(ioout,'(a)') 'WARNING: difference in emissivity flag between single and multi-instrument results'
      diffflag = .true.
    endif
    if (diffflag) then
      if (do_angcorr) then
        write(*,'(a)') 'WARNING: difference in results between single- and multi-instrument results (with angcorr)'
      else
        write(*,'(a)') 'WARNING: difference in results between single- and multi-instrument results (without angcorr)'
      endif
    endif

  enddo ! angcorr_run

! Close the output file
  close(ioout)

! Tidy up the remaining data
  call rttov_dealloc_coefs(err, coefs)

  deallocate(emissivity)
  deallocate(emis_std)
  deallocate(emis_flag)
  deallocate(chanprof)
  deallocate(channels)
  deallocate(profiles)

! End of test
PCATCH

end program rttov_camel_clim_atlas_test
