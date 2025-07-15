MODULE NameValueMap_m
  TYPE NameValue_t
    CHARACTER*256    :: name_
    INTEGER*8        :: value_
  END TYPE NameValue_t

  TYPE NameValueMap_t
    TYPE(NameValue_t), ALLOCATABLE    :: map_list(:)
    INTEGER     :: maxsize_
    INTEGER     :: size_

  CONTAINS
    PROCEDURE, PUBLIC    :: add_int4, add_int8
    GENERIC :: add => add_int4, add_int8
    PROCEDURE, PUBLIC    :: get_value_by_name
  END TYPE NameValueMap_t

CONTAINS
  FUNCTION create_map_list(m_maxsize)
    IMPLICIT NONE
    TYPE(NameValueMap_t)  :: create_map_list
    INTEGER     :: m_maxsize
    ALLOCATE (create_map_list%map_list(m_maxsize))
    create_map_list%maxsize_ = m_maxsize
    create_map_list%size_ = 0
  END FUNCTION create_map_list

  SUBROUTINE add_int4(this, m_name, m_value)
    IMPLICIT NONE
    CLASS(NameValueMap_t)     :: this
    CHARACTER(*)       :: m_name
    INTEGER*4           :: m_value
    this%size_ = this%size_ + 1
    this%map_list(this%size_)%name_ = m_name
    this%map_list(this%size_)%value_ = m_value

  END SUBROUTINE add_int4

  SUBROUTINE add_int8(this, m_name, m_value)
    IMPLICIT NONE
    CLASS(NameValueMap_t)     :: this
    CHARACTER(*)       :: m_name
    INTEGER*8           :: m_value
    this%size_ = this%size_ + 1
    this%map_list(this%size_)%name_ = m_name
    this%map_list(this%size_)%value_ = m_value

  END SUBROUTINE add_int8

  FUNCTION get_value_by_name(this, m_name)
    IMPLICIT NONE
    CLASS(NameValueMap_t)     :: this
    CHARACTER(*)          :: m_name
    INTEGER*8             :: get_value_by_name
    INTEGER             :: i
    DO i = 1, this%size_
      IF (this%map_list(i)%name_ == m_name) THEN
        get_value_by_name = this%map_list(i)%value_
      END IF
    END DO
  END FUNCTION get_value_by_name

  SUBROUTINE deallocate_map(this_map)
    TYPE(NameValueMap_t)    :: this_map

    DEALLOCATE (this_map%map_list)
  END SUBROUTINE deallocate_map
END MODULE NameValueMap_m

