! Description:
!> @file
!!   Subroutines for HDF5 I/O of emissivity structure
!
!> @brief
!!   Subroutines for HDF5 I/O of emissivity structure
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
module rttov_hdf_emissivity_io
#include "throw.h"
use parkind1
!use rttov_print_mod
use rttov_hdf_mod
use rttov_types

use hdf5
implicit none
#include "rttov_errorreport.interface"
private
public :: rttov_hdf_emissivity_init
public :: rttov_hdf_emissivity_rh
public :: rttov_hdf_emissivity_wh


contains


!> Nullify pointer to array of emissivity structures
!! param    x     pointer to array of emissivity structures
subroutine rttov_hdf_emissivity_init(x)
type(rttov_emissivity),pointer::x(:)

nullify(x)

end subroutine

!> Read emissivity structure from HDF5 file
!! param      x     unassociated pointer to array of emissivity structures
!! param[in]  lun   file ID of HDF5 file
!! param[out] err   return status
subroutine rttov_hdf_emissivity_rh(x,lun,err)
type(rttov_emissivity),pointer::x(:)
integer(hid_t),intent(in)::lun
integer(kind=jpim),intent(out)::err
character(len=lensh)::sname
real(kind=jprb), pointer::r1(:)

TRY

!call rttov_hdf_emissivity_init(x)

sname='EMIS_IN'
call read_array_hdf(lun,sname,err,pr1=r1)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

ALLOCATE( x(size(r1)), STAT = ERR )
THROWM(ERR.NE.0,"CANNOT ALLOCATE emissivity")

x(:)%emis_in = r1(:)

DEALLOCATE( r1, STAT = ERR )
THROWM(ERR.NE.0,"CANNOT DEALLOCATE R1")

nullify(r1)

sname='EMIS_OUT'
call read_array_hdf(lun,sname,err,pr1=r1)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

x(:)%emis_out= r1(:)

DEALLOCATE( r1, STAT = ERR )
THROWM(ERR.NE.0,"CANNOT DEALLOCATE R1")

nullify(r1)

sname='SPECULARITY'
call read_array_hdf(lun,sname,err,pr1=r1)
THROWM(err.ne.0,"CANNOT READ "//trim(sname))

ALLOCATE( x(size(r1)), STAT = ERR )
THROWM(ERR.NE.0,"CANNOT ALLOCATE emissivity")

x(:)%specularity = r1(:)

DEALLOCATE( r1, STAT = ERR )
THROWM(ERR.NE.0,"CANNOT DEALLOCATE R1")

nullify(r1)

err=0_jpim
CATCH
end subroutine

!> Write emissivity structure to HDF5 file
!! param      x         pointer to array of emissivity structures
!! param[in]  lun       file ID of HDF5 file
!! param[out] err       return status
!! param[in]  compress  if true will apply internal HDF5 compression, optional
subroutine rttov_hdf_emissivity_wh(x,lun,err,compress)
type(rttov_emissivity),intent(in)::x(:)
integer(hid_t),intent(in)::lun
integer(kind=jpim),intent(out)::err
logical,intent(in),optional::compress
character(len=lensh)::sname

TRY

sname='EMIS_IN'
call write_array_hdf(lun,sname,'input_emissivity',err,r1=x(:)%emis_in, compress=compress)
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
sname='EMIS_OUT'
call write_array_hdf(lun,sname,'output_emissivity',err,r1=x(:)%emis_out, compress=compress)
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
sname='SPECULARITY'
call write_array_hdf(lun,sname,'surface specularity',err,r1=x(:)%specularity, compress=compress)
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

err=0_jpim
CATCH
end subroutine

END MODULE
