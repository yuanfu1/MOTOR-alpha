
MODULE mem_allsky

!     Input arrays on model grid
  REAL, ALLOCATABLE, DIMENSION(:, :, :) :: pres_3d
  REAL, ALLOCATABLE, DIMENSION(:, :, :) :: heights_3d
  REAL, ALLOCATABLE, DIMENSION(:, :, :) :: clwc_3d
  REAL, ALLOCATABLE, DIMENSION(:, :, :) :: cice_3d
  REAL, ALLOCATABLE, DIMENSION(:, :, :) :: rain_3d
  REAL, ALLOCATABLE, DIMENSION(:, :, :) :: snow_3d
  REAL, ALLOCATABLE, DIMENSION(:, :, :) :: aod_3d

!     Local arrays on model grid
  REAL, ALLOCATABLE, DIMENSION(:, :, :) :: transm_3d
  REAL, ALLOCATABLE, DIMENSION(:, :, :, :) :: transm_4d
  REAL, ALLOCATABLE, DIMENSION(:, :, :, :) :: uprad_4d
  REAL, ALLOCATABLE, DIMENSION(:, :, :) :: upxrad_3d
  REAL, ALLOCATABLE, DIMENSION(:, :, :) :: upyrad_3d

!     2D arrays on sky grid
  REAL, ALLOCATABLE, DIMENSION(:, :) :: aod_ill_opac
  REAL, ALLOCATABLE, DIMENSION(:, :) :: aod_ill_opac_potl

!     Various non-gridded variables
  REAL ghi_sim
  INTEGER mode_aero_cld / 1/ ! treat aerosols more as clouds [1,2,3]

  PUBLIC alloc_allsky, dealloc_allsky

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE alloc_allsky(ni, nj, nk, nc, istatus)   ! I/O

!       Allocate some though not all arrays mentioned above

    ALLOCATE (pres_3d(ni, nj, nk))
    ALLOCATE (heights_3d(ni, nj, nk))
    ALLOCATE (clwc_3d(ni, nj, nk))
    ALLOCATE (cice_3d(ni, nj, nk))
    ALLOCATE (rain_3d(ni, nj, nk))
    ALLOCATE (snow_3d(ni, nj, nk))
    ALLOCATE (aod_3d(ni, nj, nk))
    ALLOCATE (transm_3d(ni, nj, nk))
    ALLOCATE (transm_4d(ni, nj, nk, nc))
    ALLOCATE (uprad_4d(ni, nj, nk, nc))
    ALLOCATE (upxrad_3d(ni, nj, nk))
    ALLOCATE (upyrad_3d(ni, nj, nk))

    WRITE (6, *) ' allsky successfully allocated'

    RETURN

  END SUBROUTINE alloc_allsky

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE dealloc_allsky()

!       Deallocate some though not all arrays mentioned above

    DEALLOCATE (pres_3d)
    DEALLOCATE (heights_3d)
    DEALLOCATE (clwc_3d)
    DEALLOCATE (cice_3d)
    DEALLOCATE (rain_3d)
    DEALLOCATE (snow_3d)
    DEALLOCATE (aod_3d)
    DEALLOCATE (transm_3d)
    DEALLOCATE (transm_4d)
    DEALLOCATE (uprad_4d)
    DEALLOCATE (upxrad_3d)
    DEALLOCATE (upyrad_3d)

    WRITE (6, *) ' allsky successfully deallocated'

    RETURN

  END SUBROUTINE dealloc_allsky

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE mem_allsky
