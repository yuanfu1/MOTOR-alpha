
MODULE mem_wind

  TYPE wind_fields
    REAL, POINTER, DIMENSION(:, :, :)  :: uanl, vanl, wanl
    REAL, POINTER, DIMENSION(:, :)    :: uanl_sfcitrp, vanl_sfcitrp
  END TYPE

  TYPE(wind_fields) :: wind

  INTEGER num_wind_obs

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE alloc_wind_arrays(nxl, nyl, nzl)

    IMPLICIT NONE
    INTEGER :: nxl, nyl, nzl

    INTEGER :: nt

    ALLOCATE (wind%uanl(nxl, nyl, nzl))
    ALLOCATE (wind%vanl(nxl, nyl, nzl))
    ALLOCATE (wind%wanl(nxl, nyl, nzl))
    ALLOCATE (wind%uanl_sfcitrp(nxl, nyl))
    ALLOCATE (wind%vanl_sfcitrp(nxl, nyl))

    RETURN
  END SUBROUTINE alloc_wind_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE deallocate_wind_arrays()

    IMPLICIT NONE

    INTEGER :: nt

    IF (ASSOCIATED(wind%uanl)) DEALLOCATE (wind%uanl)
    IF (ASSOCIATED(wind%vanl)) DEALLOCATE (wind%vanl)
    IF (ASSOCIATED(wind%wanl)) DEALLOCATE (wind%wanl)
    IF (ASSOCIATED(wind%uanl_sfcitrp)) DEALLOCATE (wind%uanl_sfcitrp)
    IF (ASSOCIATED(wind%vanl_sfcitrp)) DEALLOCATE (wind%vanl_sfcitrp)

    RETURN
  END SUBROUTINE deallocate_wind_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE nullify_wind_arrays()

    IMPLICIT NONE

    INTEGER :: nt

    IF (ASSOCIATED(wind%uanl)) NULLIFY (wind%uanl)
    IF (ASSOCIATED(wind%vanl)) NULLIFY (wind%vanl)
    IF (ASSOCIATED(wind%wanl)) NULLIFY (wind%wanl)
    IF (ASSOCIATED(wind%uanl_sfcitrp)) NULLIFY (wind%uanl_sfcitrp)
    IF (ASSOCIATED(wind%vanl_sfcitrp)) NULLIFY (wind%vanl_sfcitrp)

    RETURN
  END SUBROUTINE nullify_wind_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE
