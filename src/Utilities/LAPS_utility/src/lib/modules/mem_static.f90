!dis
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis
!dis

MODULE mem_static

  TYPE static_type
    INTEGER                     :: nx             ! X grid points
    INTEGER                     :: ny             ! Y grid points
    INTEGER                     :: nz             ! Z grid points
    REAL                        :: swlat          ! SW Corner Lat
    REAL                        :: swlon          ! SE Corner Lon
    REAL                        :: nelat          ! NE Corner Lat
    REAL                        :: nelon          ! NE Corner Lon
    REAL                        :: dx             ! X-direction grid spacing
    REAL                        :: dy             ! Y-direction grid spacing
    REAL                        :: cenlon         ! center Longitude
    REAL                        :: cenlat         ! center Latitude
    REAL                        :: stdlon         ! Orientation Longitude (polelon)
    REAL                        :: stdlat1        ! Standard Lat 1
    REAL                        :: stdlat2        ! Standard Lat 2 (polelat)
    REAL, POINTER               :: topo(:, :)      ! LAPS Topographic height
    REAL, POINTER               :: glat(:, :)      ! LAPS Latitude Array
    REAL, POINTER               :: glon(:, :)      ! LAPS Longitude Array
    REAL, POINTER               :: ldf(:, :)       ! Land fraction
    REAL, POINTER               :: pbl(:, :)       ! PBL height
    REAL, POINTER               :: akk(:, :)       ! Drag coefficient
    CHARACTER(len=132)         :: grid_type      ! Map projection type
  END TYPE

  TYPE(static_type) :: stgrid

  REAL, POINTER, DIMENSION(:, :) :: lat, lon, topo, ldf

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE alloc_static_arrays(nx, ny)

    IMPLICIT NONE
    INTEGER :: nx, ny

    ALLOCATE (stgrid%glat(nx, ny))
    ALLOCATE (stgrid%glon(nx, ny))
    ALLOCATE (stgrid%topo(nx, ny))
    ALLOCATE (stgrid%ldf(nx, ny))
    ALLOCATE (stgrid%pbl(nx, ny))
    ALLOCATE (stgrid%akk(nx, ny))

    RETURN
  END SUBROUTINE alloc_static_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE point_static_arrays()

    IMPLICIT NONE

    lat => stgrid%glat
    lon => stgrid%glon
    topo => stgrid%topo
    ldf => stgrid%ldf

    RETURN
  END SUBROUTINE point_static_arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE dealloc_static_arrays()

    IMPLICIT NONE

    DEALLOCATE (stgrid%glat)
    DEALLOCATE (stgrid%glon)
    DEALLOCATE (stgrid%topo)
    DEALLOCATE (stgrid%ldf)
    DEALLOCATE (stgrid%pbl)
    DEALLOCATE (stgrid%akk)

    RETURN
  END SUBROUTINE dealloc_static_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE fill_static(imax, jmax &
                         , n_flds, dx, dy, lov, latin1 &
                         , latin2, origin, var, comment &
                         , DATA, model, grid_spacing, map_proj &
                         , la1, lo1, la2, lo2 &
                         , center_lat, center_lon, lli, llj, uri &
                         , urj, parent_id, ratio_2_parent &
                         , status)

    IMPLICIT NONE

    INTEGER     ::  imax, jmax, n_flds, status &
                   , lli, llj, uri, urj, parent_id, ratio_2_parent
    REAL        ::  dx, dy, lov, latin1, latin2
    CHARACTER(len=*) :: origin, var(n_flds), comment(n_flds), model, map_proj
    REAL        :: DATA(imax, jmax, n_flds)
    REAL        :: grid_spacing, la1, lo1, la2, lo2, center_lat, center_lon

    stgrid%nx = imax
    stgrid%ny = jmax
    stgrid%swlat = la1
    stgrid%swlon = lo1
    stgrid%nelat = la2
    stgrid%nelon = lo2
    stgrid%dx = grid_spacing
    stgrid%dy = grid_spacing
    stgrid%cenlon = center_lon
    stgrid%cenlat = center_lat
    stgrid%stdlon = lov
    stgrid%stdlat1 = latin1
    stgrid%stdlat2 = latin2
    stgrid%topo(:, :) = DATA(:, :, 3)
    stgrid%glat(:, :) = DATA(:, :, 1)
    stgrid%glon(:, :) = DATA(:, :, 2)
    stgrid%ldf(:, :) = DATA(:, :, 4)
    stgrid%grid_type = map_proj

!print*,'-----------:',imax,jmax,la1,lo1,la2,lo2,grid_spacing,grid_spacing  &
!,center_lon,center_lat,lov,latin1,latin2,data(1,1,3),data(1,1,1),data(1,1,2) &
!,map_proj

!         var(1)    = 'LAT'
!         var(2)    = 'LON'
!         var(3)    = 'AVG'
!         var(4)    = 'LDF'
!         var(5)    = 'USE'
!         var(6)    = 'ALB'  !now used for max snow alb 2-20-03 JS.
!         var(7)    = 'STD'
!         var(8)    = 'SLN'
!         var(9)    = 'SLT'
!         var(10)   = 'STL'
!         var(11)   = 'SBL'
!         var(12)   = 'LND'  ! Land-Water Mask based on USGS landuse
!         i=12
!         do j=1,12
!            write(cat,'(i2.2)')j
!            if(cat(1:1).eq.' ')cat(1:1)='0'
!            if(cat(2:2).eq.' ')cat(2:2)='0'
!            var(i+j)= 'G'//cat   ! vegetation greenness fraction
!         enddo
!
!         var(25)='TMP'
!         i=25
!         do j=1,12
!            write(cat,'(i2.2)')j
!            if(cat(1:1).eq.' ')cat(1:1)='0'
!            if(cat(2:2).eq.' ')cat(2:2)='0'
!            var(i+j)= 'A'//cat   ! monthly albedo
!         enddo
!
!         var(ngrids)   = 'ZIN'
!
!         comment(1) = 'Lat: From MODEL by J. Smart/ S. Albers 2-03\0'
!         comment(2) = 'Lon: From MODEL by J. Smart/ S. Albers 2-03\0'
!         comment(3) = 'Average terrain elevation (m) \0'
!         comment(4) = 'Land Fraction: derived from USGS land use \0'
!         comment(5) = 'Land Use dominant category (USGS 24 Category) \0'
!         comment(6) = 'Maximum Snow Albedo; defined over land only \0'
!         comment(7) = 'Standard Deviation of Elevation data (m)\0'
!         comment(8) = 'Mean longitudinal terrain slope (m/m)\0'
!         comment(9) = 'Mean latitudinal terrain slope (m/m)\0'
!         comment(10)= 'Top layer (0-30cm) dominant category soiltype\0'
!         comment(11)= 'Bot layer (30-90cm) dominant category soiltype\0'
!         comment(12)= 'Land-Water Mask (0=water; 1 otherwise) \0'

    RETURN
  END SUBROUTINE fill_static

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE mem_static
