!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.Read_Hdf5_interfaces
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2021/01/28, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief

MODULE Read_HDF5_interfaces_m
  USE kinds_m
  USE mpObs_m, ONLY: mpObs_t
  USE HDF5

  IMPLICIT NONE

  TYPE Read_HDF5_interfaces_t

  CONTAINS
    PROCEDURE, PUBLIC  :: open_hdf
    PROCEDURE, PUBLIC  :: close_hdf
    PROCEDURE, PUBLIC  :: get_hdf_array_dims
    PROCEDURE, PUBLIC  :: get_hdf_array_rank
    GENERIC, PUBLIC    :: read_hdf => &
                          read_hdf_int4_1D, &
                          read_hdf_int4_2D, &
                          read_hdf_int4_3D, &
                          read_hdf_int4_4D, &
                          read_hdf_real4_1D, &
                          read_hdf_real4_2D, &
                          read_hdf_real4_3D, &
                          read_hdf_real4_4D, &
                          read_hdf_real8_1D, & 
                          read_hdf_real8_2D, &
                          read_hdf_real8_3D, &
                          read_hdf_real8_4D, &
                          read_hdf_int8_1D, &
                          read_hdf_int8_2D, &
                          read_hdf_int8_3D, &
                          read_hdf_int8_4D
    PROCEDURE, PUBLIC :: read_hdf_int4_1D, read_hdf_int4_2D, read_hdf_int4_3D, read_hdf_int4_4D
    PROCEDURE, PUBLIC :: read_hdf_int8_1D, read_hdf_int8_2D, read_hdf_int8_3D, read_hdf_int8_4D
    PROCEDURE, PUBLIC :: read_hdf_real4_1D, read_hdf_real4_2D, read_hdf_real4_3D, read_hdf_real4_4D 
    PROCEDURE, PUBLIC :: read_hdf_real8_1D, read_hdf_real8_2D, read_hdf_real8_3D, read_hdf_real8_4D 
    PROCEDURE, PUBLIC :: read_hdf_attribute_real4, read_hdf_global_attribute

  END TYPE Read_HDF5_interfaces_t

CONTAINS

  SUBROUTINE open_hdf(this, filename, file_id)
    IMPLICIT NONE
    CLASS(Read_HDF5_interfaces_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER(HID_T), INTENT(OUT) :: file_id         ! Dataset identifier
    ! Local variables
    INTEGER(4) :: error

    CALL h5open_f(error)
    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "filename=", filename
      RETURN
    END IF

    CALL h5fopen_f(TRIM(filename), H5F_ACC_RDONLY_F, file_id, error)

  END SUBROUTINE open_hdf

  SUBROUTINE close_hdf(this)
    CLASS(Read_HDF5_interfaces_t) :: this
    INTEGER(4) :: error

    CALL h5close_f(error)

  END SUBROUTINE close_hdf

  SUBROUTINE get_hdf_array_dims(this, file_id, varname, dims)
    IMPLICIT NONE
    INTEGER(HID_T), INTENT(IN) :: file_id
    CLASS(Read_HDF5_interfaces_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(4), INTENT(INOUT) :: dims(:)

    ! Local variables
    INTEGER(HID_T) :: var_id
    INTEGER(HID_T) :: dspace_id
    INTEGER(4) :: error
    INTEGER(4) :: rank
    INTEGER(HSIZE_T) :: dset_dims(6), max_dims(6)

    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)
    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      RETURN
    END IF

    ! Get the dataspace ID
    CALL h5dget_space_f(var_id, dspace_id, error)

    ! Getting dims from dataspace
    CALL h5sget_simple_extent_ndims_f(dspace_id, rank, error)
    CALL h5sget_simple_extent_dims_f(dspace_id, dset_dims(1:rank), max_dims(1:rank), error)
    dims = dset_dims(1:rank)

    CALL h5sclose_f(dspace_id, error)
    CALL h5dclose_f(var_id, error)

  END SUBROUTINE get_hdf_array_dims

  SUBROUTINE get_hdf_array_rank(this, file_id, varname, rank)
    IMPLICIT NONE
    INTEGER(HID_T), INTENT(IN) :: file_id
    CLASS(Read_HDF5_interfaces_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(4), INTENT(OUT) :: rank

    ! Local variables
    INTEGER(HID_T) :: var_id
    INTEGER(HID_T) :: dspace_id
    INTEGER(4) :: error

    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)
    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      RETURN
    END IF

    ! Get the dataspace ID
    CALL h5dget_space_f(var_id, dspace_id, error)

    ! Getting dims from dataspace
    CALL h5sget_simple_extent_ndims_f(dspace_id, rank, error)

    CALL h5sclose_f(dspace_id, error)
    CALL h5dclose_f(var_id, error)

  END SUBROUTINE get_hdf_array_rank

  SUBROUTINE read_hdf_int4_1D(this, file_id, varname, var_out)
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndims = 1
    CLASS(Read_HDF5_interfaces_t) :: this
    INTEGER(4), INTENT(INOUT) :: var_out(:)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(HSIZE_T) :: var_shape(ndims)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id
    INTEGER(HID_T) :: var_id
    INTEGER(4) :: error

    var_shape = SHAPE(var_out)
    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      RETURN
    END IF

    CALL h5dread_f(var_id, H5T_NATIVE_INTEGER, var_out, var_shape, error)

    CALL h5dclose_f(var_id, error)

  END SUBROUTINE read_hdf_int4_1D

  SUBROUTINE read_hdf_int4_2D(this, file_id, varname, var_out)
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndims = 2
    CLASS(Read_HDF5_interfaces_t) :: this
    INTEGER(4), INTENT(INOUT) :: var_out(:, :)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(HSIZE_T) :: var_shape(ndims)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id
    INTEGER(HID_T) :: var_id
    INTEGER(4) :: error

    var_shape = SHAPE(var_out)
    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      RETURN
    END IF

    CALL h5dread_f(var_id, H5T_NATIVE_INTEGER, var_out, var_shape, error)

    CALL h5dclose_f(var_id, error)

  END SUBROUTINE read_hdf_int4_2D

  SUBROUTINE read_hdf_int4_3D(this, file_id, varname, var_out)
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndims = 3
    CLASS(Read_HDF5_interfaces_t) :: this
    INTEGER(4), INTENT(INOUT) :: var_out(:, :, :)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(HSIZE_T) :: var_shape(ndims)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id
    INTEGER(HID_T) :: var_id
    INTEGER(4) :: error

    var_shape = SHAPE(var_out)
    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      RETURN
    END IF

    CALL h5dread_f(var_id, H5T_NATIVE_INTEGER, var_out, var_shape, error)

    CALL h5dclose_f(var_id, error)

  END SUBROUTINE read_hdf_int4_3D

  SUBROUTINE read_hdf_int4_4D(this, file_id, varname, var_out)
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndims = 4
    CLASS(Read_HDF5_interfaces_t) :: this
    INTEGER(4), INTENT(INOUT) :: var_out(:,:,:,:)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(HSIZE_T) :: var_shape(ndims)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id         
    INTEGER(HID_T) :: var_id
    INTEGER(4) :: error

    var_shape=shape(var_out)
    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      return
    END IF

    CALL h5dread_f(var_id, H5T_NATIVE_INTEGER, var_out, var_shape, error)

    CALL h5dclose_f(var_id, error)

  END SUBROUTINE read_hdf_int4_4D

  SUBROUTINE read_hdf_int8_1D(this, file_id, varname, var_out)
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndims = 1
    CLASS(Read_HDF5_interfaces_t) :: this
    INTEGER(8), INTENT(INOUT) :: var_out(:)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(HSIZE_T) :: var_shape(ndims)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id
    INTEGER(HID_T) :: var_id
    INTEGER(4) :: error

    var_shape = SHAPE(var_out)
    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      RETURN
    END IF

    BLOCK
      INTEGER(4), ALLOCATABLE :: temp(:)
      CALL h5dread_f(var_id, H5T_NATIVE_INTEGER, temp, var_shape, error)
      var_out = temp
    END BLOCK

    CALL h5dclose_f(var_id, error)

  END SUBROUTINE read_hdf_int8_1D

  SUBROUTINE read_hdf_int8_2D(this, file_id, varname, var_out)
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndims = 2
    CLASS(Read_HDF5_interfaces_t) :: this
    INTEGER(8), INTENT(INOUT) :: var_out(:, :)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(HSIZE_T) :: var_shape(ndims)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id
    INTEGER(HID_T) :: var_id
    INTEGER(4) :: error

    var_shape = SHAPE(var_out)
    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      RETURN
    END IF

    BLOCK
      INTEGER(4), ALLOCATABLE :: temp(:, :)
      CALL h5dread_f(var_id, H5T_NATIVE_INTEGER, temp, var_shape, error)
      var_out = temp
    END BLOCK

    CALL h5dclose_f(var_id, error)

  END SUBROUTINE read_hdf_int8_2D

  SUBROUTINE read_hdf_int8_3D(this, file_id, varname, var_out)
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndims = 3
    CLASS(Read_HDF5_interfaces_t) :: this
    INTEGER(8), INTENT(INOUT) :: var_out(:, :, :)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(HSIZE_T) :: var_shape(ndims)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id
    INTEGER(HID_T) :: var_id
    INTEGER(4) :: error

    var_shape = SHAPE(var_out)
    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      RETURN
    END IF

    BLOCK
      INTEGER(4), ALLOCATABLE :: temp(:, :, :)
      CALL h5dread_f(var_id, H5T_NATIVE_INTEGER, temp, var_shape, error)
      var_out = temp
    END BLOCK

    CALL h5dclose_f(var_id, error)

  END SUBROUTINE read_hdf_int8_3D

  SUBROUTINE read_hdf_int8_4D(this, file_id, varname, var_out)
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndims = 4
    CLASS(Read_HDF5_interfaces_t) :: this
    INTEGER(8), INTENT(INOUT) :: var_out(:, :, :, :)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(HSIZE_T) :: var_shape(ndims)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id         
    INTEGER(HID_T) :: var_id
    INTEGER(4) :: error

    var_shape=shape(var_out)
    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      return
    END IF

    BLOCK
    INTEGER(4), ALLOCATABLE :: temp(:, :, :, :)
    CALL h5dread_f(var_id, H5T_NATIVE_INTEGER, temp, var_shape, error)
    var_out = temp
    END BLOCK

    CALL h5dclose_f(var_id, error)

  END SUBROUTINE read_hdf_int8_4D

  SUBROUTINE read_hdf_real4_1D(this, file_id, varname, var_out)
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndims = 1
    CLASS(Read_HDF5_interfaces_t) :: this
    REAL(4), INTENT(INOUT) :: var_out(:)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(HSIZE_T) :: var_shape(ndims)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id
    INTEGER(HID_T) :: var_id
    INTEGER(4) :: error

    var_shape = SHAPE(var_out)
    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      RETURN
    END IF

    CALL h5dread_f(var_id, H5T_NATIVE_REAL, var_out, var_shape, error)

    CALL h5dclose_f(var_id, error)

  END SUBROUTINE read_hdf_real4_1D

  SUBROUTINE read_hdf_real4_2D(this, file_id, varname, var_out)
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndims = 2
    CLASS(Read_HDF5_interfaces_t) :: this
    REAL(4), INTENT(INOUT) :: var_out(:, :)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(HSIZE_T) :: var_shape(ndims)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id
    INTEGER(HID_T) :: var_id
    INTEGER(4) :: error

    var_shape = SHAPE(var_out)
    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      RETURN
    END IF

    CALL h5dread_f(var_id, H5T_NATIVE_REAL, var_out, var_shape, error)

    CALL h5dclose_f(var_id, error)

  END SUBROUTINE read_hdf_real4_2D

  SUBROUTINE read_hdf_real4_3D(this, file_id, varname, var_out)
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndims = 3
    CLASS(Read_HDF5_interfaces_t) :: this
    REAL(4), INTENT(INOUT) :: var_out(:, :, :)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(HSIZE_T) :: var_shape(ndims)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id
    INTEGER(HID_T) :: var_id
    INTEGER(4) :: error

    var_shape = SHAPE(var_out)
    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      RETURN
    END IF

    CALL h5dread_f(var_id, H5T_NATIVE_REAL, var_out, var_shape, error)

    CALL h5dclose_f(var_id, error)

  END SUBROUTINE read_hdf_real4_3D

  SUBROUTINE read_hdf_real4_4D(this, file_id, varname, var_out)
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndims = 4
    CLASS(Read_HDF5_interfaces_t) :: this
    REAL(4), INTENT(INOUT) :: var_out(:, :, :, :)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(HSIZE_T) :: var_shape(ndims)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id         
    INTEGER(HID_T) :: var_id
    INTEGER(4) :: error

    var_shape=shape(var_out)
    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      return
    END IF

    CALL h5dread_f(var_id, H5T_NATIVE_REAL, var_out, var_shape, error)

    CALL h5dclose_f(var_id, error)

  END SUBROUTINE read_hdf_real4_4D

  SUBROUTINE read_hdf_real8_1D(this, file_id, varname, var_out)
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndims = 1
    CLASS(Read_HDF5_interfaces_t) :: this
    REAL(8), INTENT(INOUT) :: var_out(:)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(HSIZE_T) :: var_shape(ndims)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id
    INTEGER(HID_T) :: var_id
    INTEGER(4) :: error

    var_shape = SHAPE(var_out)
    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      RETURN
    END IF

    CALL h5dread_f(var_id, H5T_NATIVE_DOUBLE, var_out, var_shape, error)

    CALL h5dclose_f(var_id, error)

  END SUBROUTINE read_hdf_real8_1D

  SUBROUTINE read_hdf_real8_2D(this, file_id, varname, var_out)
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndims = 2
    CLASS(Read_HDF5_interfaces_t) :: this
    REAL(8), INTENT(INOUT) :: var_out(:, :)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(HSIZE_T) :: var_shape(ndims)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id
    INTEGER(HID_T) :: var_id
    INTEGER(4) :: error

    var_shape = SHAPE(var_out)
    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      RETURN
    END IF

    CALL h5dread_f(var_id, H5T_NATIVE_DOUBLE, var_out, var_shape, error)

    CALL h5dclose_f(var_id, error)

  END SUBROUTINE read_hdf_real8_2D

  SUBROUTINE read_hdf_real8_3D(this, file_id, varname, var_out)
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndims = 3
    CLASS(Read_HDF5_interfaces_t) :: this
    REAL(8), INTENT(INOUT) :: var_out(:, :, :)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(HSIZE_T) :: var_shape(ndims)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id
    INTEGER(HID_T) :: var_id
    INTEGER(4) :: error

    var_shape = SHAPE(var_out)
    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      RETURN
    END IF

    CALL h5dread_f(var_id, H5T_NATIVE_DOUBLE, var_out, var_shape, error)

    CALL h5dclose_f(var_id, error)

  END SUBROUTINE read_hdf_real8_3D

  SUBROUTINE read_hdf_real8_4D(this, file_id, varname, var_out)
    IMPLICIT NONE
    INTEGER, PARAMETER :: ndims = 4
    CLASS(Read_HDF5_interfaces_t) :: this
    REAL(8), INTENT(INOUT) :: var_out(:,:,:,:)
    CHARACTER(LEN=*), INTENT(IN) :: varname
    INTEGER(HSIZE_T) :: var_shape(ndims)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id         
    INTEGER(HID_T) :: var_id
    INTEGER(4) :: error

    var_shape=shape(var_out)
    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      return
    END IF

    CALL h5dread_f(var_id, H5T_NATIVE_DOUBLE, var_out, var_shape, error)

    CALL h5dclose_f(var_id, error)

  END SUBROUTINE read_hdf_real8_4D

  SUBROUTINE read_hdf_attribute_real4(this, file_id, varname, attrname, attr_shape, attr_out)
    IMPLICIT NONE
    CLASS(Read_HDF5_interfaces_t) :: this
    REAL(4), INTENT(INOUT) :: attr_out
    CHARACTER(LEN=*), INTENT(IN) :: varname, attrname
    INTEGER(HSIZE_T), INTENT(IN) :: attr_shape(:)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id         
    INTEGER(HID_T) :: var_id, attr_id
    INTEGER(4) :: error

    CALL h5dopen_f(file_id, TRIM(varname), var_id, error)
    CALL h5aopen_f(var_id, attrname, attr_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing HDF5 Lib, can't read HDF5 data."
      PRINT *, "varname=", varname
      return
    END IF

    CALL h5aread_f(attr_id, H5T_NATIVE_REAL, attr_out, attr_shape, error)

    CALL h5dclose_f(var_id, error)
    CALL h5aclose_f(attr_id, error)

  END SUBROUTINE read_hdf_attribute_real4

  SUBROUTINE read_hdf_global_attribute(this, file_id, attrname, attr_shape, attr_out)
    IMPLICIT NONE
    CLASS(Read_HDF5_interfaces_t) :: this
    REAL(4), INTENT(INOUT) :: attr_out
    CHARACTER(LEN=*), INTENT(IN) :: attrname
    INTEGER(HSIZE_T), INTENT(IN) :: attr_shape(:)

    ! Local variables
    INTEGER(HID_T), INTENT(IN) :: file_id         
    INTEGER(HID_T) :: var_id, attr_id
    integer(HID_T) :: dtype_id
    integer(HID_T) :: space_id, attr_type
    integer(HSIZE_T) :: attr_size
    INTEGER(4) :: error

    ! Global attribute
    ! attrname = 'your_global_attribute_name'
    call h5aopen_name_f(file_id, attrname, attr_id, error)

    IF (error < 0) THEN
      PRINT *, "Problems initializing global attribute."
      PRINT *, "attrname=", attrname
      return
    END IF

    ! 获取全局属性的数据类型和大小
    call h5aget_type_f(attr_id, dtype_id, error)
    call h5aget_space_f(attr_id, space_id, error)
    call h5aget_type_f(space_id, attr_type, error)
    ! call h5aget_simple_extent_npoints_f(space_id, attr_size, error)

    ! 读取全局属性数据
    ! allocate(attr_out(attr_size))
    call h5aread_f(attr_id, dtype_id, attr_out, attr_shape, error)

    ! 输出全局属性数据
    print *, 'Global Attribute Name: ', trim(attrname)
    print *, 'Global Attribute Data: ', attr_out

    CALL h5aclose_f(attr_id, error)

  END SUBROUTINE read_hdf_global_attribute

END MODULE Read_HDF5_interfaces_m
