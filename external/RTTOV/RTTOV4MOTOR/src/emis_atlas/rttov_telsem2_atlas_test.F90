! Description:
!> @file
!!   Executable for basic test of the TELSEM2 emissivity atlas.
!
!> @brief
!!   Executable for basic test of the TELSEM2 emissivity atlas.
!!
!! @details
!!   This should be run using the test_telsem2_atlas.sh script in rttov_test/
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
program rttov_telsem2_atlas_test

#include "throw.h"
  use parkind1, only : jpim, jprb

  use rttov_unix_env, only : rttov_exit

  use rttov_const, only : errorstatus_success

  use rttov_types, only : &
        rttov_coefs,      &
        rttov_options,    &
        rttov_profile,    &
        rttov_chanprof

  use mod_rttov_emis_atlas, only : atlas_type_mw, rttov_emis_atlas_data

  implicit none

  integer(kind=jpim), parameter :: ioin  = 50 ! Input file unit
  integer(kind=jpim), parameter :: ioout = 51 ! Output file unit

  character(len=256)              :: coef_filename ! Coefficient filename
  character(len=256)              :: atlas_path    ! Path to atlas data
  real(kind=jprb)                 :: resol = 0.5   ! 'User-defined' resolution for testing
  integer(kind=jpim)              :: nmonth        ! Number of months for which to load data
  integer(kind=jpim), allocatable :: imonth(:)     ! Months for which to load data

  integer(kind=jpim) :: lo,hi,i,j,k
  integer(kind=jpim) :: err
  character(len=12)  :: formatstr

  real(kind=jprb)    :: zenangle
  integer(kind=jpim) :: nprof  ! Number of profiles
  integer(kind=jpim) :: nchan  ! Number of channels per profile

  type(rttov_coefs)                 :: coefs
  type(rttov_options)               :: opts
  type(rttov_profile),  allocatable :: profiles(:)
  type(rttov_chanprof), allocatable :: chanprof(:)

  type(rttov_emis_atlas_data), pointer :: atlas(:)

  real(kind=jprb), allocatable :: emissivity(:)
  real(kind=jprb), allocatable :: emis_std(:)
  real(kind=jprb), allocatable :: emis_cov(:,:,:)

#include "rttov_errorreport.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_skipcommentline.interface"

#include "rttov_setup_emis_atlas.interface"
#include "rttov_get_emis.interface"
#include "rttov_deallocate_emis_atlas.interface"

!-----------------------------------------------------------------------------------
TRY

! This test uses a simple input file for 'profiles' consisting of lat, lon and
! surface type. The controlling shell script passes the instrument and zenith angle
! for which to retrieve emissivities.

! Emissivity values, errors and covariances for all channels are returned for all
! profiles and are written to an output file for comparison with reference output.

! Atlas data for multiple months are loaded simultaneously

! Read profiles input file
  open(ioin, file='profiles_mw', status='old')
  call rttov_skipcommentline(ioin, err)
  THROWM(err.NE.0, "Error reading input profiles")

  read(ioin,*) nprof
  call rttov_skipcommentline(ioin, err)
  THROWM(err.NE.0, "Error reading input profiles")

  allocate(profiles(nprof))
  do k = 1, nprof
    read(ioin,*) profiles(k)%latitude,      &
                 profiles(k)%longitude,     &
                 profiles(k)%skin%surftype
  enddo

  close(ioin)

! Read data from control script
  read(*,*) nmonth
  ALLOCATE(imonth(nmonth), atlas(nmonth))
  read(*,*) imonth(:)
  read(*,'(a)') coef_filename
  read(*,*) zenangle
  read(*,'(a)') atlas_path

  profiles(:)%zenangle = zenangle

! Set up RTTOV coefficients
  call rttov_read_coefs(                &
              err,                      &
              coefs,                    &
              opts,                     &
              file_coef = coef_filename)

  nchan = coefs%coef%fmv_chn

! Generate the chanprof array
  allocate(chanprof(nchan*nprof))
  do k = 1, nprof
    lo = (k-1)*nchan+1
    hi = lo+nchan-1
    chanprof(lo:hi)%prof = k
    chanprof(lo:hi)%chan = (/ (j, j=1,nchan) /)
  enddo

  allocate(emissivity(nchan*nprof))
  allocate(emis_std(nchan*nprof))
  allocate(emis_cov(nprof,nchan,nchan))

!-------------------------
! Set up emissivity atlas
!-------------------------
  do i = 1, nmonth
    call rttov_setup_emis_atlas( &
                err,             &
                opts,            &
                imonth(i),       &
                atlas_type_mw,   &  ! MW atlas
                atlas(i),        &  ! TELSEM2 is the default MW atlas
                path = atlas_path)
    THROWM(err.NE.0, "Failure in setting up emissivity atlas" )
  enddo
!----------------------------
! Retrieve values from atlas
!----------------------------

! Open file for output
  open(ioout, file='output_telsem2_atlas.ascii', form='formatted', status='replace', action='write')

  do i = 1, nmonth
    write(ioout,'(a,i2.2)') 'Month ', imonth(i)
! Nominal atlas resolution
    call rttov_get_emis(             &
                err,                 &
                opts,                &
                chanprof,            &
                profiles,            &
                coefs,               &
                atlas(i),            &
                emissivity,          &
                emis_std = emis_std, &
                emis_cov = emis_cov)
    THROWM(err.NE.0, "Failure in retrieving emissivity values" )

    write(formatstr,'(a,i2.2,a)') '(',nchan,'f11.7)'
    do k = 1, nprof
      write(ioout,'(a,i2.2)') 'Nominal atlas resol. Profile ',k
      write(ioout,'(a)') ' Chan  Emissivity  Standard Dev'
      do j = 1, nchan
        lo = (k-1)*nchan
        write(ioout,'(i5,f11.4,f13.6)') j, emissivity(lo+j), emis_std(lo+j)
      enddo
      write(ioout,'(a)') 'Emissivity covariance matrix'
      do j = 1, nchan
        write(ioout,formatstr) emis_cov(k,j,1:nchan)
      enddo
    enddo

! User-defined resolution
    call rttov_get_emis(             &
                err,                 &
                opts,                &
                chanprof,            &
                profiles,            &
                coefs,               &
                atlas(i),            &
                emissivity,          &
                emis_std = emis_std, &
                emis_cov = emis_cov, &
                resolution = resol)
    THROWM(err.NE.0, "Failure in retrieving emissivity values" )

    write(formatstr,'(a,i2.2,a)') '(',nchan,'f11.7)'
    do k = 1, nprof
      write(ioout,'(a,i2.2)') 'User-defined resol. Profile ',k
      write(ioout,'(a)') ' Chan  Emissivity  Standard Dev'
      do j = 1, nchan
        lo = (k-1)*nchan
        write(ioout,'(i5,f11.4,f13.6)') j, emissivity(lo+j), emis_std(lo+j)
      enddo
      write(ioout,'(a)') 'Emissivity covariance matrix'
      do j = 1, nchan
        write(ioout,formatstr) emis_cov(k,j,1:nchan)
      enddo
    enddo

  enddo ! month loop

  close(ioout)

!-----------------------
! Deallocate atlas data
!-----------------------
  do i = 1, nmonth
    call rttov_deallocate_emis_atlas(atlas(i))
  enddo
  deallocate(imonth, atlas)

! Tidy up the remaining data
  call rttov_dealloc_coefs(err, coefs)

  deallocate(emissivity)
  deallocate(emis_std)
  deallocate(emis_cov)
  deallocate(chanprof)
  deallocate(profiles)

! End of test
PCATCH

end program rttov_telsem2_atlas_test
