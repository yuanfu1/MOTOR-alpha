! Description:
!> @file
!!   Subroutines for HDF5 I/O of chanprof structure
!
!> @brief
!!   Subroutines for HDF5 I/O of chanprof structure
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
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
module rttov_hdf_chanprof_io
#include "throw.h"
use parkind1
!use rttov_print_mod
use rttov_hdf_mod
use rttov_types

use hdf5
implicit none
#include "rttov_errorreport.interface"
private
public :: rttov_hdf_chanprof_init
public :: rttov_hdf_chanprof_rh
public :: rttov_hdf_chanprof_wh


contains

!> Nullify pointer to array of chanprof structures
!! param    x     pointer to array of chanprof structures
subroutine rttov_hdf_chanprof_init(x)
type(rttov_chanprof),pointer::x(:)

nullify(x)

end subroutine

!> Read chanprof structure from HDF5 file
!! param      x     unassociated pointer to array of chanprof structures
!! param[in]  lun   file ID of HDF5 file
!! param[out] err   return status
subroutine rttov_hdf_chanprof_rh(x,lun,err)
type(rttov_chanprof),pointer::x(:)
integer(hid_t),intent(in)::lun
integer(kind=jpim),intent(out)::err
character(len=lensh)::sname
integer(kind=jpim), pointer::i1(:)

TRY

!call rttov_hdf_chanprof_init(x)

sname='CHANNELS'
call read_array_hdf(lun,sname,err,pi1=i1)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

ALLOCATE( x(size(i1)), STAT = ERR )
THROWM(ERR.NE.0,"CANNOT ALLOCATE CHANPROF")

x(:)%chan = i1(:)

DEALLOCATE( i1, STAT = ERR )
THROWM(ERR.NE.0,"CANNOT DEALLOCATE I1")

nullify(i1)

sname='PROFILES'
call read_array_hdf(lun,sname,err,pi1=i1)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

x(:)%prof = i1(:)

DEALLOCATE( i1, STAT = ERR )
THROWM(ERR.NE.0,"CANNOT DEALLOCATE I1")

nullify(i1)

err=0_jpim
CATCH
end subroutine

!> Write chanprof structure to HDF5 file
!! param      x         pointer to array of chanprof structures
!! param[in]  lun       file ID of HDF5 file
!! param[out] err       return status
!! param[in]  compress  if true will apply internal HDF5 compression, optional
subroutine rttov_hdf_chanprof_wh(x,lun,err,compress)
type(rttov_chanprof),intent(in)::x(:)
integer(hid_t),intent(in)::lun
integer(kind=jpim),intent(out)::err
logical,intent(in),optional::compress
character(len=lensh)::sname

TRY

sname='CHANNELS'
call write_array_hdf(lun,sname,'CHANPROF_CHANNELS',err,i1=x(:)%chan, compress=compress)
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
sname='PROFILES'
call write_array_hdf(lun,sname,'CHANPROF_PROFILES',err,i1=x(:)%prof, compress=compress)
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

err=0_jpim
CATCH
end subroutine

END MODULE
