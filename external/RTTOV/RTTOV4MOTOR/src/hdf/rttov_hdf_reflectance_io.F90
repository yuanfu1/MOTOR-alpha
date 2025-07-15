! Description:
!> @file
!!   Subroutines for HDF5 I/O of reflectance structure
!
!> @brief
!!   Subroutines for HDF5 I/O of reflectance structure
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
module rttov_hdf_reflectance_io

#include "throw.h"
use parkind1
use rttov_hdf_mod
use rttov_types

use hdf5
implicit none
#include "rttov_errorreport.interface"
private
public :: rttov_hdf_reflectance_init
public :: rttov_hdf_reflectance_rh
public :: rttov_hdf_reflectance_wh


contains


!> Nullify pointer to array of reflectance structures
!! param    x     pointer to array of reflectance structures
subroutine rttov_hdf_reflectance_init(x)
type(rttov_reflectance),pointer::x(:)

nullify(x)

end subroutine

!> Read reflectance structure from HDF5 file
!! param      x     unassociated pointer to array of reflectance structures
!! param[in]  lun   file ID of HDF5 file
!! param[out] err   return status
subroutine rttov_hdf_reflectance_rh(x,lun,err)
type(rttov_reflectance),pointer::x(:)
integer(hid_t),intent(in)::lun
integer(kind=jpim),intent(out)::err
character(len=lensh)::sname
real(kind=jprb), pointer::r1(:)

TRY

!call rttov_hdf_reflectance_init(x)

sname='REFL_IN'
call read_array_hdf(lun,sname,err,pr1=r1)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

ALLOCATE( x(size(r1)), STAT = ERR )
THROWM(ERR.NE.0,"CANNOT ALLOCATE reflectance")

x(:)%refl_in = r1(:)

DEALLOCATE( r1, STAT = ERR )
THROWM(ERR.NE.0,"CANNOT DEALLOCATE R1")

nullify(r1)

sname='REFL_OUT'
call read_array_hdf(lun,sname,err,pr1=r1)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

x(:)%refl_out= r1(:)

DEALLOCATE( r1, STAT = ERR )
THROWM(ERR.NE.0,"CANNOT DEALLOCATE R1")

nullify(r1)

sname='DIFFUSE_REFL_IN'
call read_array_hdf(lun,sname,err,pr1=r1)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

ALLOCATE( x(size(r1)), STAT = ERR )
THROWM(ERR.NE.0,"CANNOT ALLOCATE reflectance")

x(:)%diffuse_refl_in = r1(:)

DEALLOCATE( r1, STAT = ERR )
THROWM(ERR.NE.0,"CANNOT DEALLOCATE R1")

nullify(r1)

sname='DIFFUSE_REFL_OUT'
call read_array_hdf(lun,sname,err,pr1=r1)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

x(:)%diffuse_refl_out= r1(:)

DEALLOCATE( r1, STAT = ERR )
THROWM(ERR.NE.0,"CANNOT DEALLOCATE R1")

nullify(r1)

sname='REFL_CLOUD_TOP'
call read_array_hdf(lun,sname,err,pr1=r1)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

x(:)%refl_cloud_top= r1(:)

DEALLOCATE( r1, STAT = ERR )
THROWM(ERR.NE.0,"CANNOT DEALLOCATE R1")

nullify(r1)

err=0_jpim
CATCH
end subroutine

!> Write reflectance structure to HDF5 file
!! param      x         pointer to array of reflectance structures
!! param[in]  lun       file ID of HDF5 file
!! param[out] err       return status
!! param[in]  compress  if true will apply internal HDF5 compression, optional
subroutine rttov_hdf_reflectance_wh(x,lun,err,compress)
type(rttov_reflectance),intent(in)::x(:)
integer(hid_t),intent(in)::lun
integer(kind=jpim),intent(out)::err
logical,intent(in),optional::compress
character(len=lensh)::sname

TRY

sname='REFL_IN'
call write_array_hdf(lun,sname,'input_reflectance',err,r1=x(:)%refl_in,compress=compress)
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
sname='REFL_OUT'
call write_array_hdf(lun,sname,'output_reflectance',err,r1=x(:)%refl_out,compress=compress)
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
sname='DIFFUSE_REFL_IN'
call write_array_hdf(lun,sname,'input_diffuse_reflectance',err,r1=x(:)%diffuse_refl_in,compress=compress)
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
sname='DIFFUSE_REFL_OUT'
call write_array_hdf(lun,sname,'output_diffuse_reflectance',err,r1=x(:)%diffuse_refl_out,compress=compress)
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
sname='REFL_CLOUD_TOP'
call write_array_hdf(lun,sname,'input_cloud_top_reflectance',err,r1=x(:)%refl_cloud_top,compress=compress)
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

err=0_jpim
CATCH
end subroutine

END MODULE
