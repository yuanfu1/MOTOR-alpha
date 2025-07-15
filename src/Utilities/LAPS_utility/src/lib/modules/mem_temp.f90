!
! Copyright (C) 2007  ; All Rights Reserved ; ATMET, LLC
!
! This file is free software; you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later version.
!
! This software is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with this
! program; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!======================================================================================

MODULE mem_temp

  TYPE temp_fields
    REAL, POINTER, DIMENSION(:, :, :)  :: t3, ht, p3
  END TYPE

  TYPE(temp_fields) :: temp

! Pointers for renaming arrays
  REAL, POINTER, DIMENSION(:, :, :) :: &
    temp_3d, heights_3d, pres_3d_pa

  INTEGER num_temp_obs

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE alloc_temp_arrays(nxl, nyl, nzl)

    IMPLICIT NONE
    INTEGER :: nxl, nyl, nzl

    INTEGER :: nt

    ALLOCATE (temp%t3(nxl, nyl, nzl))
    ALLOCATE (temp%ht(nxl, nyl, nzl))
    ALLOCATE (temp%p3(nxl, nyl, nzl))

    RETURN
  END SUBROUTINE alloc_temp_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE point_temp_arrays()

    IMPLICIT NONE

    temp_3d => temp%t3
    heights_3d => temp%ht
    pres_3d_pa => temp%p3

    RETURN
  END SUBROUTINE point_temp_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE deallocate_temp_arrays()

    IMPLICIT NONE

    INTEGER :: nt

    IF (ASSOCIATED(temp%t3)) DEALLOCATE (temp%t3)
    IF (ASSOCIATED(temp%ht)) DEALLOCATE (temp%ht)
    IF (ASSOCIATED(temp%p3)) DEALLOCATE (temp%p3)

    RETURN
  END SUBROUTINE deallocate_temp_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE nullify_temp_arrays()

    IMPLICIT NONE

    INTEGER :: nt

    IF (ASSOCIATED(temp%t3)) NULLIFY (temp%t3)
    IF (ASSOCIATED(temp%ht)) NULLIFY (temp%ht)
    IF (ASSOCIATED(temp%p3)) NULLIFY (temp%p3)

    RETURN
  END SUBROUTINE nullify_temp_arrays

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE
