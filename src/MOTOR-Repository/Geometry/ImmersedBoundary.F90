!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.MOTOR-Repository.Geometry
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.1
! HISTORY           :
!   Created by Yuanfu Xie, 2024/08, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!> @brief
!! Surface geometry for normal, tangential (2 components) and projections. This module is implemented
!! following immersed boundary method but modified to define the interface on the closest grids next
!! next to the actual
!! @author Yuanfu Xie
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @warning
! @attention
MODULE immersedBoundary_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  !USE linkedList_m, ONLY: linkedList_t,insert,free_list
  USE linked_list_module, ONLY: node, insert, print_list, free_list
  USE intArray_m, ONLY: intArray_t

  PUBLIC :: immersedBoundary_t
  TYPE :: immersedBoundary_t
    INTEGER(i_kind), ALLOCATABLE :: idxSurface(:) ! cell indices of singleGrid closest to the surface
    TYPE(SingleGrid_t), POINTER :: sg

  CONTAINS
    PROCEDURE :: initialization                 ! Initialize the boundary
    PROCEDURE :: gridProjection                 ! Project grid vector on to the interior grid cells
  END TYPE immersedBoundary_t

CONTAINS

  SUBROUTINE initialization(this, sg)
    IMPLICIT NONE

    CLASS(immersedBoundary_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg

    ! Local variables:
    INTEGER(i_kind) :: numSurfPts, i, j, k, ct
    INTEGER(i_kind), ALLOCATABLE :: values(:)
    ! TYPE(linkedList_t), pointer :: head => null()

    TYPE :: pair_t
      INTEGER(i_kind) :: iv, ih
    END TYPE
    TYPE(node), POINTER :: head => NULL(), temp
    TYPE(intArray_t), ALLOCATABLE :: e, surf(:)

    this%sg => sg

    ! CALL insert(head,10)
    ! CALL insert(head,20)

    ! CALL to1DArray(head,numSurfPts,values)
    ! PRINT*,'Number of grid cells: ',numSurfPts
    ! DO i=1,numSurfPts
    !   PRINT*,'Immersed initialization: ',i,values(i)
    ! END DO

    PRINT *, 'SG topo: ', MAXVAL(sg%topo), sg%ztop

    PRINT *, 'zHght 1: ', sg%zHght(1:5, 1000)
    PRINT *, 'topog 1: ', sg%topo(1000)
    PRINT *, 'above 1: ', sg%aboveTopo(:, 1000)
    PRINT *, 'zHght 2: ', sg%zHght(1:5, 2000)
    PRINT *, 'topog 2: ', sg%topo(2000)
    PRINT *, 'above 2: ', sg%aboveTopo(:, 2000)
    ! Search the surface or interface cells:
    ct = 0
    DO k = 1, this%sg%vLevel
      DO i = 1, this%sg%num_cell
        IF (this%sg%aboveTopo(k, i)) THEN
          DO j = 1, UBOUND(this%sg%cell_adnb, 1)
            IF (this%sg%cell_adnb(j, i) .GT. 0) THEN
              IF (.NOT. this%sg%aboveTopo(k, this%sg%cell_adnb(j, i))) THEN
                !ALLOCATE(e)
                e = intArray_t(2)
                e%iarray(1) = k; e%iarray(2) = i
                ct = ct + 1
                PRINT *, 'Surf: ', ct, i, k
                CALL insert(head, e)
                DEALLOCATE (e)
              END IF
            END IF
          END DO
        END IF
      END DO
    END DO
    PRINT *, 'How many elements?', ct
    numSurfPts = head%count
    PRINT *, 'Total number of interface points: ', numSurfPts, sg%num_cell
    CALL print_list(head)
    ! ALLOCATE(surf(numSurfPts))
    ! temp = head
    ! i = 0
    ! DO WHILE (ASSOCIATED(temp))
    !   i = i+1
    !   surf(i) = temp%value
    !   temp => temp%next
    ! END DO

    CALL free_list(head)

  END SUBROUTINE initialization

  SUBROUTINE gridProjection(this)
    IMPLICIT NONE

    CLASS(immersedBoundary_t) :: this
  END SUBROUTINE gridProjection

END MODULE immersedBoundary_m
