! Description:
!> @file
!!   Subroutines for HDF5 I/O of profiles structure
!
!> @brief
!!   Subroutines for HDF5 I/O of profiles structure
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
     MODULE RTTOV_HDF_PROFILES
#include "throw.h"
      USE RTTOV_TYPES
      USE RTTOV_HDF_MOD
      USE HDF5
      USE H5LT
! 
      USE RTTOV_HDF_PROFILE_IO
      USE RTTOV_HDF_S2M_IO
      USE RTTOV_HDF_SKIN_IO
      !USE RTTOV_HDF_OPTIONS_IO
      !USE RTTOV_HDF_CHANPROF_IO
      !USE RTTOV_HDF_EMISSIVITY_IO
      !USE RTTOV_HDF_REFLECTANCE_IO
       
      IMPLICIT NONE
#include "rttov_errorreport.interface"

      PRIVATE
      PUBLIC :: RTTOV_HDF_ONE_PROFILE_WH, RTTOV_HDF_ONE_PROFILE_RH
      PUBLIC :: RTTOV_HDF_KMATRIX_WH
      
CONTAINS

!> Write a single profile structure to HDF5 file
!! param[in]  x         profile structure
!! param[in]  lun       file ID of HDF5 file
!! param[out] err       return status
!! param[in]  compress  if true will apply internal HDF5 compression, optional
      SUBROUTINE RTTOV_HDF_ONE_PROFILE_WH(X,LUN,ERR,COMPRESS)
      TYPE(RTTOV_PROFILE),INTENT(IN) ::X
      INTEGER(HID_T),INTENT(IN)      ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT) ::ERR
      LOGICAL,INTENT(IN),OPTIONAL    ::COMPRESS
!
      INTEGER(HID_T) :: G_ID_SUB
!
TRY

        CALL RTTOV_HDF_profile_WH( x, LUN, ERR, COMPRESS=COMPRESS )
        THROWM(ERR.NE.0,"CANNOT WRITE PROFILE")
        
        CALL H5LTSET_ATTRIBUTE_STRING_F(LUN, '.', "Description",   &
        "This is a RTTOV profile structure" // &
        CHAR(0), ERR )
        THROWM(ERR.NE.0,"CANNOT WRITE ATTRIBUTE")

        CALL MKPAR( LUN, "S2M", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CREATE GROUP S2M")
      
        CALL RTTOV_HDF_s2m_WH( x%S2M, G_ID_SUB, ERR, COMPRESS=COMPRESS )
        THROWM(ERR.NE.0,"CANNOT WRITE PROFILE %S2M ")

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP S2M")

        CALL MKPAR( LUN, "SKIN", G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CREATE GROUP SKIN")
      
        CALL RTTOV_HDF_skin_WH( x%SKIN, G_ID_SUB, ERR, COMPRESS=COMPRESS )
        THROWM(ERR.NE.0,"CANNOT WRITE PROFILE % SKIN")

        CALL H5GCLOSE_F( G_ID_SUB, ERR )
        THROWM(ERR.NE.0,"CANNOT CLOSE GROUP SKIN")

CATCH
      END SUBROUTINE

!> Read a single profile structure from HDF5 file
!! param[out] x         profile structure
!! param[in]  lun       file ID of HDF5 file
!! param[out] err       return status
      SUBROUTINE RTTOV_HDF_ONE_PROFILE_RH(X,LUN,ERR)
      TYPE(RTTOV_PROFILE),INTENT(OUT) ::X
      INTEGER(HID_T),INTENT(IN)       ::LUN
      INTEGER(KIND=JPIM),INTENT(OUT)  ::ERR

      INTEGER(HID_T) :: G_ID_SUB

TRY

      CALL RTTOV_HDF_profile_RH( x, LUN, ERR )
      THROWM(ERR.NE.0,"CANNOT READ PROFILE ")

      CALL H5GOPEN_F( LUN, "S2M", G_ID_SUB, ERR )
      THROWM(ERR.NE.0,"CANNOT OPEN GROUP S2M")

      CALL RTTOV_HDF_s2m_RH( x%s2m, G_ID_SUB, ERR )
      THROWM(ERR.NE.0,"CANNOT READ GROUP S2M")

      CALL H5GCLOSE_F( G_ID_SUB, ERR )
      THROWM(ERR.NE.0,"CANNOT CLOSE GROUP S2M")

      CALL H5GOPEN_F( LUN, "SKIN", G_ID_SUB, ERR )
      THROWM(ERR.NE.0,"CANNOT OPEN GROUP SKIN")

      CALL RTTOV_HDF_skin_RH( x%SKIN, G_ID_SUB, ERR )
      THROW(ERR.NE.0)

      CALL H5GCLOSE_F( G_ID_SUB, ERR )
      THROWM(ERR.NE.0,"CANNOT CLOSE GROUP SKIN")

CATCH
      END SUBROUTINE





!> Write a profile Jacobian to HDF5 file
!! param[in]  x         array of profile structures containing Jacobians
!! param[in]  opts      options used when simulation was run
!! param[in]  lun       file ID of HDF5 file
!! param[out] err       return status
!! param[in]  compress  if true will apply internal HDF5 compression, optional
subroutine rttov_hdf_kmatrix_wh(x,opts,lun,err,compress)
type(rttov_profile),  intent(in)  ::x(:)
type(rttov_options),  intent(in)  ::opts
integer(hid_t),       intent(in)  ::lun
integer(kind=jpim),   intent(out) ::err
logical,intent(in),optional       ::compress

real(kind=jprb),allocatable :: r1_chn(:)
real(kind=jprb),allocatable :: r2_lev(:,:)
real(kind=jprb),allocatable :: r2_lay(:,:)
real(kind=jprb),allocatable :: r3_lay(:,:,:)

integer(hid_t) ::g_id_sub

character(len=lensh) :: sname
character(len=32)    :: sunit
integer(kind=jpim)   :: nchn, i
integer(kind=jpim)   :: nlevels, nlayers, nc

#define R1C(N)  do i = 1, nchn; r1_chn(i)     = x(i)%N;      enddo
#define R2LV(N) do i = 1, nchn; r2_lev(i,:)   = x(i)%N(:);   enddo
#define R2LA(N) do i = 1, nchn; r2_lay(i,:)   = x(i)%N(:);   enddo
#define R3LA(N) do i = 1, nchn; r3_lay(i,:,:) = x(i)%N(:,:); enddo

TRY

if( opts%rt_all%switchrad ) then
   sunit = 'K / '
else
   sunit = '(mw/cm-1/ster/sq.m) / '
   if( opts%rt_ir%pc%addpc ) then
     sunit = '(pcscore) / '
   endif
endif

nchn    = size(x)
nlevels = x(1)%nlevels
nlayers = x(1)%nlayers


allocate(r1_chn(nchn),  STAT=ERR)
THROWM(ERR.NE.0,"COULD NOT ALLOCATE R1 ")
allocate(r2_lev(nchn,nlevels),  STAT=ERR)
THROWM(ERR.NE.0,"COULD NOT ALLOCATE R2 ")


sname='ID'
call write_array_hdf(lun,sname,'',err,c0=x(1)%ID )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='DATE'
call write_array_hdf(lun,sname,'Year Month Day',err,i1=x(1)%date)
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='TIME'
call write_array_hdf(lun,sname,'Hour Minute Second',err,i1=x(1)%TIME )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='NLEVELS'
call write_array_hdf(lun,sname,'number of atmospheric levels',err,i0=x(1)%NLEVELS )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='NLAYERS'
call write_array_hdf(lun,sname,'number of atmospheric layers',err,i0=x(1)%NLAYERS )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='GAS_UNITS'
call write_array_hdf(lun,sname,'units for gas profiles',err,i0=x(1)%GAS_UNITS )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='MMR_CLDAER'
call write_array_hdf(lun,sname,'units for cloud and aerosol profiles',err,l0=x(1)%MMR_CLDAER )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='P'
if(associated(x(1)%P))then
  R2LV(P)
  call write_array_hdf(lun,sname,'pressure',err,r2=r2_lev , units = TRIM(sunit)//'hPa',&
  & compress=compress)
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
endif

sname='T'
if(associated(x(1)%T))then
  R2LV(T)
  call write_array_hdf(lun,sname,'temperature',err,r2=r2_lev , units = TRIM(sunit)//'K',&
  & compress=compress)
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
endif

sname='Q'
if(associated(x(1)%Q))then
  R2LV(Q)
  call write_array_hdf(lun,sname,'water vapour',err,r2=r2_lev , units = TRIM(sunit)//'ppmv or kg/kg',&
  & compress=compress)
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
endif

sname='O3'
if(associated(x(1)%O3))then
  R2LV(O3)
  call write_array_hdf(lun,sname,'ozone',err,r2=r2_lev , units = TRIM(sunit)//'ppmv or kg/kg',&
  & compress=compress)
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
endif

sname='CO2'
if(associated(x(1)%CO2))then
  R2LV(CO2)
  call write_array_hdf(lun,sname,'carbon dioxide',err,r2=r2_lev , units = TRIM(sunit)//'ppmv or kg/kg',&
  & compress=compress)
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
endif

sname='N2O'
if(associated(x(1)%N2O))then
  R2LV(N2O)
  call write_array_hdf(lun,sname,'N2O',err,r2=r2_lev , units = TRIM(sunit)//'ppmv or kg/kg',&
  & compress=compress)
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
endif

sname='CO'
if(associated(x(1)%CO))then
  R2LV(CO)
  call write_array_hdf(lun,sname,'CO',err,r2=r2_lev , units = TRIM(sunit)//'ppmv or kg/kg',&
  & compress=compress)
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
endif

sname='CH4'
if(associated(x(1)%CH4))then
  R2LV(CH4)
  call write_array_hdf(lun,sname,'CH4',err,r2=r2_lev , units = TRIM(sunit)//'ppmv or kg/kg',&
  & compress=compress)
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
endif

sname='SO2'
if(associated(x(1)%SO2))then
  R2LV(SO2)
  call write_array_hdf(lun,sname,'SO2',err,r2=r2_lev , units = TRIM(sunit)//'ppmv or kg/kg',&
  & compress=compress)
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
endif

sname='CLW'
if(associated(x(1)%CLW))then
  R2LV(CLW)
  call write_array_hdf(lun,sname,'cloud liquid water (kg/kg)',err,r2=r2_lev, units = TRIM(sunit)//'kg/kg',&
  & compress=compress)
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
endif

sname='AEROSOLS'
if(associated(x(1)%AEROSOLS))then
  nc = size(x(1)%AEROSOLS, 1_jpim)
  allocate(r3_lay(nchn,nc,nlayers))
  R3LA(AEROSOLS)
  call write_array_hdf(lun,sname,'Aerosols',err,r3=r3_lay , units = TRIM(sunit)//'kg/kg or cm^-3',&
  & compress=compress)
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
  deallocate(r3_lay)
endif

sname='CLOUD'
if(associated(x(1)%CLOUD))then
  nc = size(x(1)%CLOUD, 1_jpim)
  allocate(r3_lay(nchn,nc,nlayers))
  R3LA(CLOUD)
  call write_array_hdf(lun,sname,'Clouds',err,r3=r3_lay , units = TRIM(sunit)//'kg/kg or g/m^3',&
  & compress=compress)
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
  deallocate(r3_lay)
endif

sname='CFRAC'
if(associated(x(1)%CFRAC))then
  allocate(r2_lay(nchn,nlayers))
  R2LA(CFRAC)
  call write_array_hdf(lun,sname,'Cloud fraction (0-1)',err,r2=r2_lay, units = TRIM(sunit)//'0:1',&
  & compress=compress)
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
  deallocate(r2_lay)  
endif

sname='CLWDE'
if(associated(x(1)%CLWDE))then
  allocate(r2_lay(nchn,nlayers),  STAT=ERR)
  THROWM(ERR.NE.0,"COULD NOT ALLOCATE R2 ")
  R2LA(CLWDE)
  call write_array_hdf(lun,sname,'Cloud liquid water diameter',err,r2=r2_lay , units = TRIM(sunit)//'microns',&
  & compress=compress)
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
  deallocate(r2_lay)
endif

sname='ICEDE'
if(associated(x(1)%ICEDE))then
  allocate(r2_lay(nchn,nlayers),  STAT=ERR)
  THROWM(ERR.NE.0,"COULD NOT ALLOCATE R2 ")
  R2LA(ICEDE)
  call write_array_hdf(lun,sname,'Ice crystals diameter',err,r2=r2_lay , units = TRIM(sunit)//'microns',&
  & compress=compress)
  THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))
  deallocate(r2_lay)
endif

sname='CLWDE_PARAM'
call write_array_hdf(lun,sname,'CLWDE_PARAM',err,i0=x(1)%CLWDE_PARAM )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='CLW_SCHEME'
call write_array_hdf(lun,sname,'CLW_SCHEME',err,i0=x(1)%CLW_SCHEME )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='ICEDE_PARAM'
call write_array_hdf(lun,sname,'ICEDE_PARAM',err,i0=x(1)%ICEDE_PARAM )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='ICE_SCHEME'
call write_array_hdf(lun,sname,'ICE_SCHEME',err,i0=x(1)%ICE_SCHEME )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='ZENANGLE'
call write_array_hdf(lun,sname,'Satellite zenith angle',err,r0=x(1)%ZENANGLE , units = 'degree')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='AZANGLE'
call write_array_hdf(lun,sname,'Satellite azimuth angle',err,r0=x(1)%AZANGLE , units = 'degree')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='SUNZENANGLE'
call write_array_hdf(lun,sname,'Sun zenith angle',err,r0=x(1)%SUNZENANGLE , units = 'degree')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='SUNAZANGLE'
call write_array_hdf(lun,sname,'Sun azimuth angle',err,r0=x(1)%SUNAZANGLE , units = 'degree')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='ELEVATION'
call write_array_hdf(lun,sname,'Surface elevaton',err,r0=x(1)%ELEVATION , units = 'km')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='LATITUDE'
call write_array_hdf(lun,sname,'Latitude (deg)',err,r0=x(1)%LATITUDE , units = 'degree')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='BE'
call write_array_hdf(lun,sname,'Earth magnetic field strength',err,r0=x(1)%BE , units = 'Gauss')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='COSBK'
call write_array_hdf(lun,sname,'cosine of the angle between the Earth magnetic field and wave propagation direction',&
                   err,r0=x(1)%COSBK )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='CTP'
call write_array_hdf(lun,sname,'Black body cloud top pressure',err,r1=x(:)%CTP , units = TRIM(sunit)//'hPa')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='CFRACTION'
call write_array_hdf(lun,sname,'Black body cloud fraction (0-1), 1 for 100% cloud cover',err,r1=x(:)%CFRACTION , &
                      units = TRIM(sunit)//'0:1')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))



CALL MKPAR( LUN, "S2M", G_ID_SUB, ERR )
THROWM(ERR.NE.0,"CANNOT CREATE GROUP S2M")

sname='T'
call write_array_hdf(g_id_sub,sname,'temperature',err,r1=x(:)%s2m%T , units = TRIM(sunit)//'K')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='Q'
call write_array_hdf(g_id_sub,sname,'water vapour',err,r1=x(:)%s2m%Q , units = TRIM(sunit)//'ppmv or kg/kg')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='O'
call write_array_hdf(g_id_sub,sname,'ozone',err,r1=x(:)%s2m%O , units = TRIM(sunit)//'ppmv or kg/kg')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='P'
call write_array_hdf(g_id_sub,sname,'surface pressure',err,r1=x(:)%s2m%P , units = TRIM(sunit)//'hPa')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='U'
call write_array_hdf(g_id_sub,sname,'U 10m wind component',err,r1=x(:)%s2m%U , units = TRIM(sunit)//'m/s')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='V'
call write_array_hdf(g_id_sub,sname,'V 10m wind component',err,r1=x(:)%s2m%V , units = TRIM(sunit)//'m/s')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='WFETC'
call write_array_hdf(g_id_sub,sname,'Wind fetch',err,r1=x(:)%s2m%WFETC , units = TRIM(sunit)//'m')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

CALL H5GCLOSE_F( G_ID_SUB, ERR )
THROWM(ERR.NE.0,"CANNOT CLOSE GROUP S2M")



CALL MKPAR( LUN, "SKIN", G_ID_SUB, ERR )
THROWM(ERR.NE.0,"CANNOT CREATE GROUP SKIN")

sname='SURFTYPE'
call write_array_hdf(g_id_sub,sname,'0=land, 1=sea, 2=sea-ice',err,i0=x(1)%skin%SURFTYPE )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='WATERTYPE'
call write_array_hdf(g_id_sub,sname,'0=fresh water, 1=ocean water',err,i0=x(1)%skin%WATERTYPE )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='T'
call write_array_hdf(g_id_sub,sname,'radiative skin temperature (K)',err,r1=x(:)%skin%T , units = TRIM(sunit)//'K')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='SNOW_FRACTION'
call write_array_hdf(lun,sname,'snow coverage fraction for IR emissivity atlas (0-1)',err,r0=x(1)%skin%SNOW_FRACTION,units='0:1')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='SOIL_MOISTURE'
call write_array_hdf(lun,sname,'soil moisture',err,r0=x(1)%skin%SOIL_MOISTURE , units = 'm^3/m^3')
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='SALINITY'
call write_array_hdf(g_id_sub,sname,'practical salinity unit %o',err,r1=x(:)%skin%SALINITY )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='FOAM_FRACTION'
call write_array_hdf(g_id_sub,sname,'ocean foam coverage fraction passed to FASTEM (0-1)',err,r1=x(:)%skin%FOAM_FRACTION )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))

sname='FASTEM'
call write_array_hdf(g_id_sub,sname,'land/sea-ice surface parameters for FASTEM',err,r1=x(1)%skin%FASTEM )
THROWM(err.ne.0,"CANNOT WRITE "//trim(sname))


CALL H5GCLOSE_F( G_ID_SUB, ERR )
THROWM(ERR.NE.0,"CANNOT CLOSE GROUP SKIN")


deallocate(r1_chn, r2_lev)

err=0_jpim
CATCH
end subroutine


END MODULE
