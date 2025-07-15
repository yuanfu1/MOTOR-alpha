!!--------------------------------------------------------------------------------------------------
! PROJECT           : Utilities.Utility
! AFFILIATION       : GuangDOng-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.1
! HISTORY           :
!   Created by Yuanfu Xie, 2024/08, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
!> @brief
!! A set of data structures
!! @author Yuanfu Xie
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @warning
! @attention
! ModIFied from a chatGPS version with a predefined data type for the list by Yuanfu Xie 2024-08-26

MODULE linked_list_module
  USE kinds_m, ONLY: i_kind
  USE intArray_m, ONLY: intArray_t

  IMPLICIT NONE
  PRIVATE

  TYPE :: node
    INTEGER(i_kind) :: count = 0
    CLASS(*), ALLOCATABLE :: VALUE
    TYPE(node), POINTER :: next => NULL()
  END TYPE node

  PUBLIC :: node, insert, print_list, free_list, get_type

CONTAINS

  SUBROUTINE insert(head, val)
    TYPE(node), POINTER :: head, temp, new_node
    CLASS(*), TARGET, INTENT(in) :: val

    ! Allocate a new node
    ALLOCATE (new_node)
    SELECT TYPE (val)
    TYPE is (intArray_t)
      ALLOCATE (intArray_t :: new_node%VALUE)
      new_node%VALUE = val
      PRINT *, 'Array size: ', val%num_elements
      ! type is (character(len=*))
      !   allocate(character(len=len(elememt)) :: new_node%value)
      ! type is (integer)
      !   allocate(integer :: new_node%value)
      ! type is (real)
      !   allocate(real :: new_node%value)
      ! type is (logical)
      !   allocate(logical :: new_node%value)
    CLASS default
      PRINT *, "Unsupported data type!"
      STOP
    END SELECT

    ! Assign the value to the new node
    new_node%VALUE = val
    new_node%next => NULL()

    ! Insert the new node into the list
    IF (.NOT. ASSOCIATED(head)) THEN
      head => new_node
      head%count = 1
    ELSE
      temp => head
      DO WHILE (ASSOCIATED(temp%next))
        temp => temp%next
      END DO
      temp%next => new_node
      head%count = head%count + 1
    END IF
    PRINT *, 'HEAD count: ', head%count
  END SUBROUTINE insert

  SUBROUTINE print_list(head)
    TYPE(node), POINTER :: head
    TYPE(node), POINTER :: temp

    ! Local variables:
    INTEGER(i_kind) :: i, j

    temp => head
    i = 0
    DO WHILE (ASSOCIATED(temp))
      !call get_type(temp%value)
      ASSOCIATE (val => temp%VALUE)
        SELECT TYPE (val)
        TYPE is (intArray_t)
          i = i + 1
          PRINT *, 'List: ', i, val%num_elements, (val%iarray(j), j=1, val%num_elements)
        END SELECT
      END ASSOCIATE
      temp => temp%next
    END DO
  END SUBROUTINE print_list

  SUBROUTINE get_type(val)
    CLASS(*), INTENT(in) :: val

    SELECT TYPE (val)
    TYPE is (intArray_t)
      PRINT *, "Person: Name =", val%num_elements
    TYPE is (CHARACTER(len=*))
      PRINT *, "Character:", val
    TYPE is (INTEGER)
      PRINT *, "Integer:", val
    TYPE is (REAL)
      PRINT *, "Real:", val
    TYPE is (LOGICAL)
      PRINT *, "Logical:", val
    CLASS default
      PRINT *, "Unknown or unsupported data type!"
    END SELECT
  END SUBROUTINE get_type

  SUBROUTINE free_list(head)
    TYPE(node), POINTER :: head, temp

    DO WHILE (ASSOCIATED(head))
      temp => head%next
      DEALLOCATE (head%VALUE)
      DEALLOCATE (head)
      head => temp
    END DO
  END SUBROUTINE free_list

END MODULE linked_list_module
