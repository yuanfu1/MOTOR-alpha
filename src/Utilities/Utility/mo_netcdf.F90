!> \file mo_netcdf.f90

!> \brief NetCDF Fortran 90 interface wrapper

!> \details A wrapper around the NetCDF Fortran 90 interface.
!
!> \authors David Schaefer
!> \date Jun 2015

MODULE mo_netcdf

  ! This module provides a thin wrapper around the NetCDF Fortran 90 interface,
  ! following a object-oriented approach.

  ! Written  David Schaefer, Jun 2015
  ! Modified Matthias Cuntz, Jan 2016 - compiled with PGI Fortran rev 15.9 - no automatic allocation of left-hand-side
  ! Modified Ricardo Torres, Feb 2017 - add derived type NcGroup and NcAttributable. NcAttributable is the base derived type,
  !                                     NcGroup and NcVariable are extended from it. NcDataset is extended from NcGroup. No more
  !                                     duplicated routines to set attributes.

  ! License
  ! -------
  ! GNU Lesser General Public License http://www.gnu.org/licenses/

  USE mo_types
  USE netcdf, ONLY: &
    nf90_open, nf90_close, nf90_strerror, nf90_def_dim, nf90_def_var, &
    nf90_put_var, nf90_get_var, nf90_put_att, nf90_get_att, &
    nf90_inquire, nf90_inq_dimid, nf90_inquire_dimension, &
    nf90_inq_varid, nf90_inq_varids, nf90_inquire_variable, nf90_inquire_attribute, &
    nf90_inq_ncid, nf90_inq_grp_parent, nf90_inq_grpname, nf90_def_grp, &
    nf90_rename_dim, nf90_rename_var, nf90_rename_att, nf90_sync, &
    NF90_OPEN, NF90_NETCDF4, NF90_CREATE, NF90_WRITE, NF90_NOWRITE, &
    NF90_BYTE, NF90_SHORT, NF90_INT, NF90_FLOAT, NF90_DOUBLE, &
    NF90_FILL_BYTE, NF90_FILL_SHORT, NF90_FILL_INT, NF90_FILL_FLOAT, NF90_FILL_DOUBLE, &
    NF90_NOERR, NF90_UNLIMITED, NF90_GLOBAL, NF90_SHARE, NF90_HDF5, &
    NF90_64BIT_OFFSET, NF90_CLASSIC_MODEL

  IMPLICIT NONE

  ! --------------------------------------------------------------------------------------

  TYPE, ABSTRACT :: NcBase

    INTEGER(i32) :: id

  CONTAINS

    PROCEDURE(getNameInterface), DEFERRED   :: getName
    PROCEDURE(getParentInterface), DEFERRED :: getParent

  END TYPE NcBase

  TYPE, ABSTRACT, EXTENDS(NcBase) :: NcAttributable

  CONTAINS

    PROCEDURE, PUBLIC  :: hasAttribute
    PROCEDURE, PUBLIC  :: renameAttribute

    PROCEDURE, PRIVATE :: getAttributableIds
    PROCEDURE, PRIVATE :: setAttributeChar
    PROCEDURE, PRIVATE :: setAttributeI8
    PROCEDURE, PRIVATE :: setAttributeI16
    PROCEDURE, PRIVATE :: setAttributeI32
    PROCEDURE, PRIVATE :: setAttributeF32
    PROCEDURE, PRIVATE :: setAttributeF64

    PROCEDURE, PRIVATE :: getAttributeChar
    PROCEDURE, PRIVATE :: getAttributeI8
    PROCEDURE, PRIVATE :: getAttributeI16
    PROCEDURE, PRIVATE :: getAttributeI32
    PROCEDURE, PRIVATE :: getAttributeF32
    PROCEDURE, PRIVATE :: getAttributeF64

    GENERIC, PUBLIC :: setAttribute => &
      setAttributeChar, &
      setAttributeI8, &
      setAttributeI16, &
      setAttributeI32, &
      setAttributeF32, &
      setAttributeF64

    GENERIC, PUBLIC :: getAttribute => &
      getAttributeChar, &
      getAttributeI8, &
      getAttributeI16, &
      getAttributeI32, &
      getAttributeF32, &
      getAttributeF64
  END TYPE NcAttributable

  ! --------------------------------------------------------------------------------------

  TYPE, EXTENDS(NcAttributable) :: NcGroup

  CONTAINS

    ! getter
    PROCEDURE, PRIVATE :: getVariableIds
    PROCEDURE, PUBLIC  :: getVariables
    PROCEDURE, PUBLIC  :: getUnlimitedDimension
    PROCEDURE, PUBLIC  :: getNoVariables

    PROCEDURE, PRIVATE :: getDimensionByName
    PROCEDURE, PRIVATE :: getDimensionById

    PROCEDURE, PUBLIC  :: getParent => getGroupParent
    PROCEDURE, PUBLIC  :: getName => getGroupName
    PROCEDURE, PUBLIC  :: getGroup => getGroupByName
    PROCEDURE, PUBLIC  :: getVariable => getVariableByName
    GENERIC, PUBLIC  :: getDimension => &
      getDimensionById, &
      getDimensionByName

    ! checker
    PROCEDURE, PUBLIC  :: hasVariable
    PROCEDURE, PUBLIC  :: hasDimension
    PROCEDURE, PUBLIC  :: hasGroup
    PROCEDURE, PUBLIC  :: isUnlimited => isDatasetUnlimited

    ! setter
    PROCEDURE, PUBLIC  :: setGroup
    PROCEDURE, PRIVATE :: setLimitedDimension
    PROCEDURE, PRIVATE :: setUnlimitedDimension
    GENERIC, PUBLIC    :: setDimension => &
      setLimitedDimension, &
      setUnlimitedDimension
    PROCEDURE, PRIVATE :: setVariableWithTypes
    PROCEDURE, PRIVATE :: setVariableWithNames
    PROCEDURE, PRIVATE :: setVariableWithIds
    PROCEDURE, PRIVATE :: setVariableScalar

    GENERIC, PUBLIC  :: setVariable => &
      setVariableWithNames, &
      setVariableWithTypes, &
      setVariableWithIds, &
      setVariableScalar

  END TYPE NcGroup

  INTERFACE NcGroup
    PROCEDURE newNcGroup
  END INTERFACE NcGroup

  ! --------------------------------------------------------------------------------------

  TYPE, EXTENDS(NcGroup) :: NcDataset

    CHARACTER(256) :: fname !> Filename of the opened dataset
    CHARACTER(1)   :: mode  !> File open mode

  CONTAINS

    PROCEDURE, PUBLIC :: sync
    PROCEDURE, PUBLIC :: CLOSE

  END TYPE NcDataset

  INTERFACE NcDataset
    PROCEDURE newNcDataset
  END INTERFACE NcDataset

! --------------------------------------------------------------------------------------

  TYPE, EXTENDS(NcBase) :: NcDimension

    TYPE(NcGroup)     :: parent  !> The dimension's parent

  CONTAINS

    PROCEDURE, PUBLIC :: renameDimension
    PROCEDURE, PUBLIC :: getParent => getDimensionParent
    PROCEDURE, PUBLIC :: getName => getDimensionName
    PROCEDURE, PUBLIC :: getLength => getDimensionLength

    PROCEDURE, PUBLIC :: isUnlimited => isUnlimitedDimension

  END TYPE NcDimension

  INTERFACE NcDimension
    PROCEDURE newNcDimension
  END INTERFACE NcDimension

! --------------------------------------------------------------------------------------

  TYPE, EXTENDS(NcAttributable) :: NcVariable

    TYPE(NcGroup)      :: parent   !> The variables's parent

  CONTAINS

    PROCEDURE, PUBLIC  :: renameVariable
    PROCEDURE, PUBLIC  :: getParent => getVariableParent
    PROCEDURE, PUBLIC  :: getName => getVariableName
    PROCEDURE, PRIVATE :: getSlicingShape

    PROCEDURE, PRIVATE :: setDataScalarI8
    PROCEDURE, PRIVATE :: setData1dI8
    PROCEDURE, PRIVATE :: setData2dI8
    PROCEDURE, PRIVATE :: setData3dI8
    PROCEDURE, PRIVATE :: setData4dI8
    PROCEDURE, PRIVATE :: setData5dI8
    PROCEDURE, PRIVATE :: setDataScalarI16
    PROCEDURE, PRIVATE :: setData1dI16
    PROCEDURE, PRIVATE :: setData2dI16
    PROCEDURE, PRIVATE :: setData3dI16
    PROCEDURE, PRIVATE :: setData4dI16
    PROCEDURE, PRIVATE :: setData5dI16
    PROCEDURE, PRIVATE :: setDataScalarI32
    PROCEDURE, PRIVATE :: setData1dI32
    PROCEDURE, PRIVATE :: setData2dI32
    PROCEDURE, PRIVATE :: setData3dI32
    PROCEDURE, PRIVATE :: setData4dI32
    PROCEDURE, PRIVATE :: setData5dI32
    PROCEDURE, PRIVATE :: setDataScalarF32
    PROCEDURE, PRIVATE :: setData1dF32
    PROCEDURE, PRIVATE :: setData2dF32
    PROCEDURE, PRIVATE :: setData3dF32
    PROCEDURE, PRIVATE :: setData4dF32
    PROCEDURE, PRIVATE :: setData5dF32
    PROCEDURE, PRIVATE :: setDataScalarF64
    PROCEDURE, PRIVATE :: setData1dF64
    PROCEDURE, PRIVATE :: setData2dF64
    PROCEDURE, PRIVATE :: setData3dF64
    PROCEDURE, PRIVATE :: setData4dF64
    PROCEDURE, PRIVATE :: setData5dF64

    PROCEDURE, PRIVATE :: getDataScalarI8
    PROCEDURE, PRIVATE :: getData1dI8
    PROCEDURE, PRIVATE :: getData2dI8
    PROCEDURE, PRIVATE :: getData3dI8
    PROCEDURE, PRIVATE :: getData4dI8
    PROCEDURE, PRIVATE :: getData5dI8
    PROCEDURE, PRIVATE :: getDataScalarI16
    PROCEDURE, PRIVATE :: getData1dI16
    PROCEDURE, PRIVATE :: getData2dI16
    PROCEDURE, PRIVATE :: getData3dI16
    PROCEDURE, PRIVATE :: getData4dI16
    PROCEDURE, PRIVATE :: getData5dI16
    PROCEDURE, PRIVATE :: getDataScalarI32
    PROCEDURE, PRIVATE :: getData1dI32
    PROCEDURE, PRIVATE :: getData2dI32
    PROCEDURE, PRIVATE :: getData3dI32
    PROCEDURE, PRIVATE :: getData4dI32
    PROCEDURE, PRIVATE :: getData5dI32
    PROCEDURE, PRIVATE :: getDataScalarF32
    PROCEDURE, PRIVATE :: getData1dF32
    PROCEDURE, PRIVATE :: getData2dF32
    PROCEDURE, PRIVATE :: getData3dF32
    PROCEDURE, PRIVATE :: getData4dF32
    PROCEDURE, PRIVATE :: getData5dF32
    PROCEDURE, PRIVATE :: getDataScalarF64
    PROCEDURE, PRIVATE :: getData1dF64
    PROCEDURE, PRIVATE :: getData2dF64
    PROCEDURE, PRIVATE :: getData3dF64
    PROCEDURE, PRIVATE :: getData4dF64
    PROCEDURE, PRIVATE :: getData5dF64
    PROCEDURE, PRIVATE :: getData1dStr  ! get string vector variable, TS@0207

    PROCEDURE, PRIVATE :: setVariableFillValueI8
    PROCEDURE, PRIVATE :: setVariableFillValueI16
    PROCEDURE, PRIVATE :: setVariableFillValueI32
    PROCEDURE, PRIVATE :: setVariableFillValueF32
    PROCEDURE, PRIVATE :: setVariableFillValueF64

    PROCEDURE, PRIVATE :: getVariableFillValueI8
    PROCEDURE, PRIVATE :: getVariableFillValueI16
    PROCEDURE, PRIVATE :: getVariableFillValueI32
    PROCEDURE, PRIVATE :: getVariableFillValueF32
    PROCEDURE, PRIVATE :: getVariableFillValueF64

    PROCEDURE, PUBLIC  :: getNoDimensions

    PROCEDURE, PUBLIC  :: getDimensions => getVariableDimensions

    PROCEDURE, PUBLIC  :: getRank => getVariableRank

    PROCEDURE, PUBLIC  :: getShape => getVariableShape

    PROCEDURE, PUBLIC  :: getDtype => getVariableDtype

    PROCEDURE, PUBLIC  :: isUnlimited => isUnlimitedVariable

    GENERIC, PUBLIC :: setData => &
      setDataScalarI8, &
      setData1dI8, &
      setData2dI8, &
      setData3dI8, &
      setData4dI8, &
      setData5dI8, &
      setDataScalarI16, &
      setData1dI16, &
      setData2dI16, &
      setData3dI16, &
      setData4dI16, &
      setData5dI16, &
      setDataScalarI32, &
      setData1dI32, &
      setData2dI32, &
      setData3dI32, &
      setData4dI32, &
      setData5dI32, &
      setDataScalarF32, &
      setData1dF32, &
      setData2dF32, &
      setData3dF32, &
      setData4dF32, &
      setData5dF32, &
      setDataScalarF64, &
      setData1dF64, &
      setData2dF64, &
      setData3dF64, &
      setData4dF64, &
      setData5dF64

    GENERIC, PUBLIC :: getData => &
      getDataScalarI8, &
      getData1dI8, &
      getData2dI8, &
      getData3dI8, &
      getData4dI8, &
      getData5dI8, &
      getDataScalarI16, &
      getData1dI16, &
      getData2dI16, &
      getData3dI16, &
      getData4dI16, &
      getData5dI16, &
      getDataScalarI32, &
      getData1dI32, &
      getData2dI32, &
      getData3dI32, &
      getData4dI32, &
      getData5dI32, &
      getDataScalarF32, &
      getData1dF32, &
      getData2dF32, &
      getData3dF32, &
      getData4dF32, &
      getData5dF32, &
      getDataScalarF64, &
      getData1dF64, &
      getData2dF64, &
      getData3dF64, &
      getData4dF64, &
      getData5dF64, &
      getData1dStr  ! get string verctor, TS@0207

    GENERIC, PUBLIC :: setFillValue => &
      setVariableFillValueI8, &
      setVariableFillValueI16, &
      setVariableFillValueI32, &
      setVariableFillValueF32, &
      setVariableFillValueF64

    GENERIC, PUBLIC :: getFillValue => &
      getVariableFillValueI8, &
      getVariableFillValueI16, &
      getVariableFillValueI32, &
      getVariableFillValueF32, &
      getVariableFillValueF64

  END TYPE NcVariable

  INTERFACE NcVariable
    PROCEDURE newNcVariable
  END INTERFACE NcVariable
  ! --------------------------------------------------------------------------------------

  ! abstract interfaces
  INTERFACE
    FUNCTION getNameInterface(self)
      IMPORT NcBase
      CLASS(NcBase), INTENT(in) :: self
      CHARACTER(len=256)        :: getNameInterface
    END FUNCTION getNameInterface

    FUNCTION getParentInterface(self)
      IMPORT NcBase, NcGroup
      CLASS(NcBase), INTENT(in) :: self
      TYPE(NcGroup)             :: getParentInterface
    END FUNCTION getParentInterface
  END INTERFACE

  INTERFACE OPERATOR(==)
    PROCEDURE equalNcBases
  END INTERFACE OPERATOR(==)

CONTAINS

  FUNCTION newNcDataset(fname, fmode, cmode) RESULT(out)
    CHARACTER(*), INTENT(in)              :: fname
    CHARACTER(1), INTENT(in)              :: fmode
    CHARACTER(*), INTENT(inout), OPTIONAL :: cmode
    INTEGER(i32)                          :: status
    TYPE(NcDataset)                       :: out

    SELECT CASE (fmode)
    CASE ("w")
      status = nf90_create(TRIM(fname), getCreationMode(cmode), out%id)
    CASE ("r")
      status = nf90_open(TRIM(fname), NF90_NOWRITE, out%id)
    CASE ("a")
      status = nf90_open(TRIM(fname), NF90_WRITE, out%id)
    CASE default
      WRITE (*, *) "Mode argument must be in 'w','r','a' ! "
      STOP 1
    END SELECT
    CALL check(status, "Failed to open file: "//fname)

    out%fname = fname
    out%mode = fmode
  END FUNCTION newNcDataset

  FUNCTION newNcVariable(id, parent) RESULT(out)
    INTEGER(i32), INTENT(in) :: id
    TYPE(NcGroup), INTENT(in) :: parent
    TYPE(NcVariable)          :: out

    out%id = id
    out%parent = parent
  END FUNCTION newNcVariable

  FUNCTION newNcDimension(id, parent) RESULT(out)
    INTEGER(i32), INTENT(in) :: id
    TYPE(NcGroup), INTENT(in) :: parent
    TYPE(NcDimension)         :: out

    out%id = id
    out%parent = parent
  END FUNCTION newNcDimension

  FUNCTION newNcGroup(id) RESULT(out)
    INTEGER(i32), INTENT(in) :: id
    TYPE(NcGroup)                :: out

    out%id = id
  END FUNCTION newNcGroup

  SUBROUTINE sync(self)
    CLASS(NcDataset) :: self

    CALL check(nf90_sync(self%id), "Failed to sync file: "//self%fname)
  END SUBROUTINE sync

  SUBROUTINE CLOSE (self)
    CLASS(NcDataset) :: self

    CALL check(nf90_close(self%id), "Failed to close file: "//self%fname)
  END SUBROUTINE CLOSE

  FUNCTION setGroup(self, name)
    CLASS(NcGroup), INTENT(inout) :: self
    CHARACTER(*), INTENT(in)    :: name
    INTEGER(i32)                  :: id
    TYPE(NcGroup)                 :: setGroup

    CALL check(nf90_def_grp(self%id, name, id), "Failed to create new group: "//name)
    setGroup = NcGroup(id)
  END FUNCTION setGroup

  FUNCTION getGroupParent(self)
    CLASS(NcGroup), INTENT(in) :: self
    INTEGER(i32)               :: id
    TYPE(NcGroup)              :: getGroupParent

    CALL check(nf90_inq_grp_parent(self%id, id), "Failed to get parent group of: "//self%getName())
    getGroupParent = NcGroup(id)
  END FUNCTION getGroupParent

  FUNCTION getGroupName(self)
    CLASS(NcGroup), INTENT(in) :: self
    CHARACTER(256)             :: getGroupName

    CALL check(nf90_inq_grpname(self%id, getGroupName), "Failed to inquire group name")
  END FUNCTION getGroupName

  FUNCTION getNoVariables(self)
    CLASS(NcGroup), INTENT(in) :: self
    INTEGER(i32)               :: getNoVariables

    CALL check(nf90_inquire(self%id, nvariables=getNoVariables), "Failed inquire number of variables")
  END FUNCTION getNoVariables

  FUNCTION getDimensionParent(self)
    CLASS(NcDimension), INTENT(in) :: self
    TYPE(NcGroup)                  :: getDimensionParent

    getDimensionParent = self%parent
  END FUNCTION getDimensionParent

  FUNCTION getVariableParent(self)
    CLASS(NcVariable), INTENT(in) :: self
    TYPE(NcGroup)                  :: getVariableParent

    getVariableParent = self%parent
  END FUNCTION getVariableParent

  FUNCTION getVariableIds(self)
    CLASS(NcGroup), INTENT(in)              :: self
    INTEGER(i32), DIMENSION(:), ALLOCATABLE :: getVariableIds
    INTEGER(i32)                            :: tmp

    ALLOCATE (getVariableIds(self%getNoVariables()))
    CALL check(nf90_inq_varids(self%id, tmp, getVariableIds), "Failed to inquire variable ids")
  END FUNCTION getVariableIds

  FUNCTION getVariables(self)
    CLASS(NcGroup), INTENT(in)                  :: self
    TYPE(NcVariable), DIMENSION(:), ALLOCATABLE :: getVariables
    INTEGER(i32), DIMENSION(:), ALLOCATABLE     :: varids
    INTEGER(i32)                                :: i, nvars

    nvars = self%getNoVariables()
    ALLOCATE (getVariables(nvars), varids(nvars))

    varids = self%getVariableIds()
    DO i = 1, SIZE(varids)
      getVariables(i) = NcVariable(varids(i), self)
    END DO

  END FUNCTION getVariables

  FUNCTION getDimensionName(self)
    CLASS(NcDimension), INTENT(in) :: self
    CHARACTER(len=256)             :: getDimensionName

    CALL check(nf90_inquire_dimension(self%parent%id, self%id, name=getDimensionName), &
               "Failed to inquire dimension name")
  END FUNCTION getDimensionName

  FUNCTION getDimensionLength(self)
    CLASS(NcDimension), INTENT(in) :: self
    INTEGER(i32)                    :: getDimensionLength

    CALL check(nf90_inquire_dimension(self%parent%id, self%id, len=getDimensionLength), &
               "Failed to inquire dimension: "//self%getName())
  END FUNCTION getDimensionLength

  FUNCTION isDatasetUnlimited(self)
    CLASS(NcGroup), INTENT(in)   :: self
    LOGICAL                      :: isDatasetUnlimited
    INTEGER(i32)                  :: dimid

    CALL check(nf90_inquire(self%id, unlimitedDimId=dimid), &
               "Failed to inquire group "//self%getName())
    isDatasetUnlimited = (dimid .NE. -1)
  END FUNCTION isDatasetUnlimited

  FUNCTION getUnlimitedDimension(self)
    CLASS(NcGroup), INTENT(in)   :: self
    TYPE(NcDimension)            :: getUnlimitedDimension
    INTEGER(i32)                  :: dimid

    CALL check(nf90_inquire(self%id, unlimitedDimId=dimid), &
               "Failed to inquire group "//self%getName())

    IF (dimid .EQ. -1) THEN
      WRITE (*, *) "Dataset has no unlimited dimension"
      STOP 1
    END IF

    getUnlimitedDimension = self%getDimension(dimid)
  END FUNCTION getUnlimitedDimension

  FUNCTION equalNcBases(left, right) RESULT(out)
    CLASS(NcBase), INTENT(in) :: left, right
    LOGICAL                   :: out

    out = (left%id .EQ. right%id)
  END FUNCTION equalNcBases

  FUNCTION isUnlimitedDimension(self)
    CLASS(NcDimension), INTENT(in) :: self
    LOGICAL                        :: isUnlimitedDimension

    isUnlimitedDimension = .FALSE.
    IF (self%parent%isUnlimited()) THEN
      isUnlimitedDimension = (self == self%parent%getUnlimitedDimension())
    END IF
  END FUNCTION isUnlimitedDimension

  FUNCTION setUnlimitedDimension(self, name)
    CLASS(NcGroup), INTENT(in)           :: self
    CHARACTER(*), INTENT(in)           :: name
    TYPE(NcDimension)                    :: setUnlimitedDimension

    setUnlimitedDimension = self%setLimitedDimension(name, -1)
  END FUNCTION setUnlimitedDimension

  FUNCTION setLimitedDimension(self, name, length)
    CLASS(NcGroup), INTENT(in) :: self
    CHARACTER(*), INTENT(in) :: name
    INTEGER(i32), INTENT(in) :: length
    TYPE(NcDimension)          :: setLimitedDimension
    INTEGER(i32)               :: id, dimlength

    IF (length .LE. 0) THEN
      dimlength = NF90_UNLIMITED
    ELSE
      dimlength = length
    END IF

    CALL check(nf90_def_dim(self%id, name, dimlength, id), &
               "Failed to create dimension: "//name)

    setLimitedDimension = NcDimension(id, self)
  END FUNCTION setLimitedDimension

  FUNCTION hasVariable(self, name)
    CLASS(NcGroup), INTENT(in) :: self
    CHARACTER(*), INTENT(in) :: name
    LOGICAL                      :: hasVariable
    INTEGER(i32)                  :: tmpid

    hasVariable = (nf90_inq_varid(self%id, name, tmpid) .EQ. NF90_NOERR)
  END FUNCTION hasVariable

  FUNCTION hasDimension(self, name)
    CLASS(NcGroup), INTENT(in) :: self
    CHARACTER(*), INTENT(in) :: name
    LOGICAL                      :: hasDimension
    INTEGER(i32)                  :: tmpid

    hasDimension = (nf90_inq_dimid(self%id, name, tmpid) .EQ. NF90_NOERR)
  END FUNCTION hasDimension

  FUNCTION hasGroup(self, name)
    CLASS(NcGroup), INTENT(in) :: self
    CHARACTER(*), INTENT(in) :: name
    LOGICAL                      :: hasGroup
    INTEGER(i32)                  :: tmpid

    hasGroup = (nf90_inq_ncid(self%id, name, tmpid) .EQ. NF90_NOERR)
  END FUNCTION hasGroup

  FUNCTION setVariableWithIds(self, name, dtype, dimensions, CONTIGUOUS, &
                              chunksizes, deflate_level, shuffle, fletcher32, endianness, &
                              cache_size, cache_nelems, cache_preemption)
    CLASS(NcGroup), INTENT(in)           :: self
    CHARACTER(*), INTENT(in)           :: name
    CHARACTER(*), INTENT(in)           :: dtype
    INTEGER(i32), INTENT(in)           :: dimensions(:)
    LOGICAL, INTENT(in), OPTIONAL :: CONTIGUOUS, shuffle, fletcher32
    INTEGER(i32), INTENT(in), OPTIONAL :: endianness, deflate_level, cache_size, &
                                          cache_nelems, cache_preemption, chunksizes(:)
    TYPE(NcVariable)                       :: setVariableWithIds
    INTEGER(i32)                            :: varid, status

    status = nf90_def_var(self%id, name, getDtypeFromString(dtype), dimensions, varid, CONTIGUOUS, &
                          chunksizes, deflate_level, shuffle, fletcher32, endianness, &
                          cache_size, cache_nelems, cache_preemption)
    CALL check(status, "setVariableWithIds - Failed to create variable: "//name)
    setVariableWithIds = NcVariable(varid, self)
  END FUNCTION setVariableWithIds

  FUNCTION setVariableWithNames(self, name, dtype, dimensions, CONTIGUOUS, &
                                chunksizes, deflate_level, shuffle, fletcher32, endianness, &
                                cache_size, cache_nelems, cache_preemption)

    CLASS(NcGroup), INTENT(in)              :: self
    CHARACTER(*), INTENT(in)              :: name
    CHARACTER(*), INTENT(in)              :: dtype
    CHARACTER(*), INTENT(in)              :: dimensions(:)
    LOGICAL, INTENT(in), OPTIONAL    :: CONTIGUOUS, shuffle, fletcher32
    INTEGER(i32), INTENT(in), OPTIONAL    :: endianness, deflate_level, cache_size, &
                                             cache_nelems, cache_preemption, chunksizes(:)
    TYPE(NcVariable)                          :: setVariableWithNames
    TYPE(NcDimension)                         :: dim
    INTEGER(i32)                               :: i, dimids(SIZE(dimensions))

    DO i = 1, SIZE(dimensions)
      dim = self%getDimension(dimensions(i))
      dimids(i) = dim%id
    END DO

    setVariableWithNames = setVariableWithIds(self, name, dtype, dimids, CONTIGUOUS, &
                                              chunksizes, deflate_level, shuffle, fletcher32, endianness, &
                                              cache_size, cache_nelems, cache_preemption)
  END FUNCTION setVariableWithNames

  FUNCTION setVariableWithTypes(self, name, dtype, dimensions, CONTIGUOUS, &
                                chunksizes, deflate_level, shuffle, fletcher32, endianness, &
                                cache_size, cache_nelems, cache_preemption)
    CLASS(NcGroup), INTENT(in)              :: self
    CHARACTER(*), INTENT(in)              :: name
    CHARACTER(*), INTENT(in)              :: dtype
    TYPE(NcDimension), INTENT(in)              :: dimensions(:)
    LOGICAL, INTENT(in), OPTIONAL    :: CONTIGUOUS, shuffle, fletcher32
    INTEGER(i32), INTENT(in), OPTIONAL    :: endianness, deflate_level, cache_size, &
                                             cache_nelems, cache_preemption, chunksizes(:)
    TYPE(NcVariable)                           :: setVariableWithTypes
    TYPE(NcDimension)                          :: dim
    INTEGER(i32)                                :: i, dimids(SIZE(dimensions))

    DO i = 1, SIZE(dimensions)
      dim = dimensions(i)
      dimids(i) = dim%id
    END DO

    setVariableWithTypes = setVariableWithIds(self, name, dtype, dimids, CONTIGUOUS, &
                                              chunksizes, deflate_level, shuffle, fletcher32, endianness, &
                                              cache_size, cache_nelems, cache_preemption)
  END FUNCTION setVariableWithTypes

  FUNCTION setVariableScalar(self, name, dtype)
    CLASS(NcGroup), INTENT(in)           :: self
    CHARACTER(*), INTENT(in)           :: name
    CHARACTER(*), INTENT(in)           :: dtype
    TYPE(NcVariable)                       :: setVariableScalar
    INTEGER(i32)                           :: varid, status

    status = nf90_def_var(self%id, name, getDtypeFromString(dtype), varid=varid)
    CALL check(status, "setVariableScalar - Failed to create variable: "//name)
    setVariableScalar = NcVariable(varid, self)
  END FUNCTION setVariableScalar

  FUNCTION getDimensionById(self, id)
    CLASS(NcGroup), INTENT(in) :: self
    INTEGER(i32)                  :: id
    TYPE(NcDimension)            :: getDimensionById
    CHARACTER(32)                :: msg, name

    WRITE (msg, *) id
    CALL check(nf90_inquire_dimension(self%id, id, name), &
               "Could not inquire dimension: "//msg)
    getDimensionById = NcDimension(id, self)
  END FUNCTION getDimensionById

  FUNCTION getDimensionByName(self, name)
    CLASS(NcGroup), INTENT(in) :: self
    CHARACTER(*)                 :: name
    TYPE(NcDimension)            :: getDimensionByName
    INTEGER(i32)                  :: id

    CALL check(nf90_inq_dimid(self%id, name, id), &
               "Could not inquire dimension: "//name)
    getDimensionByName = self%getDimensionById(id)
  END FUNCTION getDimensionByName

  FUNCTION getGroupByName(self, name)
    CLASS(NcGroup), INTENT(in) :: self
    CHARACTER(*), INTENT(in) :: name
    TYPE(NcGroup)              :: getGroupByName
    INTEGER(i32)                :: id

    CALL check(nf90_inq_ncid(self%id, name, id), &
               "Could not inquire variable: "//name)
    getGroupByName = NcGroup(id)
  END FUNCTION getGroupByName

  FUNCTION getVariableByName(self, name)
    CLASS(NcGroup), INTENT(in) :: self
    CHARACTER(*), INTENT(in) :: name
    TYPE(NcVariable)             :: getVariableByName
    INTEGER(i32)                  :: id

    CALL check(nf90_inq_varid(self%id, name, id), &
               "Could not inquire variable: "//name)
    getVariableByName = NcVariable(id, self)

  END FUNCTION getVariableByName

  FUNCTION getVariableName(self)
    CLASS(NcVariable), INTENT(in) :: self
    CHARACTER(len=256)            :: getVariableName

    CALL check(nf90_inquire_variable(self%parent%id, self%id, name=getVariableName), &
               "Could not inquire variable name")
  END FUNCTION getVariableName

  FUNCTION getNoDimensions(self)
    CLASS(NcVariable), INTENT(in) :: self
    INTEGER(i32)                   :: getNoDimensions

    CALL check(nf90_inquire_variable(self%parent%id, self%id, ndims=getNoDimensions), &
               "Could not inquire variable: "//self%getName())
  END FUNCTION getNoDimensions

  FUNCTION getVariableDimensions(self)
    CLASS(NcVariable), INTENT(in)  :: self
    TYPE(NcDimension), ALLOCATABLE :: getVariableDimensions(:)
    INTEGER(i32), ALLOCATABLE :: dimids(:)
    INTEGER(i32)                    :: ii, ndims

    ndims = self%getNoDimensions()
    ALLOCATE (dimids(ndims), getVariableDimensions(ndims))
    CALL check(nf90_inquire_variable(self%parent%id, self%id, dimids=dimids), &
               "Could not inquire variable: "//self%getName())

    DO ii = 1, ndims
      getVariableDimensions(ii) = self%parent%getDimension(dimids(ii))
    END DO
  END FUNCTION getVariableDimensions

  FUNCTION getVariableShape(self)
    CLASS(NcVariable), INTENT(in)  :: self
    INTEGER(i32), ALLOCATABLE :: getVariableShape(:)
    TYPE(NcDimension), ALLOCATABLE :: dims(:)
    INTEGER(i32)                   :: ii, ndims

    ndims = self%getNoDimensions()
    ALLOCATE (getVariableShape(ndims), dims(ndims))

    dims = self%getDimensions()
    DO ii = 1, SIZE(dims)
      getVariableShape(ii) = dims(ii)%getLength()
    END DO
  END FUNCTION getVariableShape

  FUNCTION getVariableRank(self)
    CLASS(NcVariable), INTENT(in)  :: self
    INTEGER(i32)                   :: getVariableRank

    getVariableRank = SIZE(self%getDimensions())
  END FUNCTION getVariableRank

  FUNCTION getVariableDtype(self)
    CLASS(NcVariable), INTENT(in) :: self
    INTEGER(i32)                   :: dtype
    CHARACTER(3)                  :: getVariableDtype

    CALL check(nf90_inquire_variable(self%parent%id, self%id, xtype=dtype), &
               "Could not inquire variable: "//self%getName())
    getVariableDtype = getDtypeFromInteger(dtype)
  END FUNCTION getVariableDtype

  FUNCTION isUnlimitedVariable(self)
    CLASS(NcVariable), INTENT(in)  :: self
    LOGICAL                        :: isUnlimitedVariable
    TYPE(NcDimension), ALLOCATABLE :: dims(:)
    TYPE(NcDimension)              :: dim
    INTEGER(i32)                    :: ii

    ALLOCATE (dims(self%getNoDimensions()))

    isUnlimitedVariable = .FALSE.
    dims = self%getDimensions()

    DO ii = 1, SIZE(dims)
      dim = dims(ii)
      IF (dim%isUnlimited()) THEN
        isUnlimitedVariable = .TRUE.
      END IF
    END DO
  END FUNCTION isUnlimitedVariable

  LOGICAL FUNCTION hasAttribute(self, name)
    CLASS(NcAttributable), INTENT(inout) :: self
    CHARACTER(*), INTENT(in) :: name
    INTEGER(i32)                   :: status

    SELECT TYPE (self)
    CLASS is (NcGroup)
      status = nf90_inquire_attribute(self%id, NF90_GLOBAL, name)
    CLASS is (NcVariable)
      status = nf90_inquire_attribute(self%parent%id, self%id, name)
    END SELECT

    hasAttribute = (status .EQ. NF90_NOERR)
  END FUNCTION hasAttribute

  SUBROUTINE setAttributeChar(self, name, DATA)
    CLASS(NcAttributable), INTENT(in) :: self
    CHARACTER(*), INTENT(in) :: name
    CHARACTER(*), INTENT(in) :: DATA
    INTEGER(i32)                      :: ids(2)

    ids = self%getAttributableIds()
    CALL check(nf90_put_att(ids(1), ids(2), name, DATA), &
               "Failed to write attribute: "//name)

  END SUBROUTINE setAttributeChar

  SUBROUTINE setAttributeI8(self, name, DATA)
    CLASS(NcAttributable), INTENT(in) :: self
    CHARACTER(*), INTENT(in) :: name
    INTEGER(i8), INTENT(in) :: DATA
    INTEGER(i32)                      :: ids(2)

    ids = self%getAttributableIds()
    CALL check(nf90_put_att(ids(1), ids(2), name, DATA), &
               "Failed to write attribute: "//name)

  END SUBROUTINE setAttributeI8

  SUBROUTINE setAttributeI16(self, name, DATA)
    CLASS(NcAttributable), INTENT(in) :: self
    CHARACTER(*), INTENT(in) :: name
    INTEGER(i16), INTENT(in) :: DATA
    INTEGER(i32)                      :: ids(2)

    ids = self%getAttributableIds()
    CALL check(nf90_put_att(ids(1), ids(2), name, DATA), &
               "Failed to write attribute: "//name)

  END SUBROUTINE setAttributeI16

  SUBROUTINE setAttributeI32(self, name, DATA)
    CLASS(NcAttributable), INTENT(in) :: self
    CHARACTER(*), INTENT(in) :: name
    INTEGER(i32), INTENT(in) :: DATA
    INTEGER(i32)                      :: ids(2)

    ids = self%getAttributableIds()
    CALL check(nf90_put_att(ids(1), ids(2), name, DATA), &
               "Failed to write attribute: "//name)

  END SUBROUTINE setAttributeI32

  SUBROUTINE setAttributeF32(self, name, DATA)
    CLASS(NcAttributable), INTENT(in) :: self
    CHARACTER(*), INTENT(in) :: name
    REAL(f32), INTENT(in) :: DATA
    INTEGER(i32)                      :: ids(2)

    ids = self%getAttributableIds()
    CALL check(nf90_put_att(ids(1), ids(2), name, DATA), &
               "Failed to write attribute: "//name)

  END SUBROUTINE setAttributeF32

  SUBROUTINE setAttributeF64(self, name, DATA)
    CLASS(NcAttributable), INTENT(in) :: self
    CHARACTER(*), INTENT(in) :: name
    REAL(f64), INTENT(in) :: DATA
    INTEGER(i32)                      :: ids(2)

    ids = self%getAttributableIds()
    CALL check(nf90_put_att(ids(1), ids(2), name, DATA), &
               "Failed to write attribute: "//name)

  END SUBROUTINE setAttributeF64

  SUBROUTINE getAttributeChar(self, name, avalue)
    CLASS(NcAttributable), INTENT(in)  :: self
    CHARACTER(*), INTENT(in)  :: name
    CHARACTER(*), INTENT(out) :: avalue
    INTEGER(i32)                       :: length, ids(2)

    ids = self%getAttributableIds()
    CALL check(nf90_inquire_attribute(ids(1), ids(2), name, len=length), &
               "Could not inquire attribute "//name)
    CALL check(nf90_get_att(ids(1), ids(2), name, avalue), &
               "Could not read attribute "//name)
  END SUBROUTINE getAttributeChar

  SUBROUTINE getAttributeI8(self, name, avalue)
    CLASS(NcAttributable), INTENT(in)  :: self
    CHARACTER(*), INTENT(in)  :: name
    INTEGER(i8), INTENT(out) :: avalue
    INTEGER(i32)                        :: length, ids(2)

    ids = self%getAttributableIds()
    CALL check(nf90_inquire_attribute(ids(1), ids(2), name, len=length), &
               "Could not inquire attribute "//name)
    CALL check(nf90_get_att(ids(1), ids(2), name, avalue), &
               "Could not read attribute "//name)

  END SUBROUTINE getAttributeI8

  SUBROUTINE getAttributeI16(self, name, avalue)
    CLASS(NcAttributable), INTENT(in)  :: self
    CHARACTER(*), INTENT(in)  :: name
    INTEGER(i16), INTENT(out) :: avalue
    INTEGER(i32)                        :: length, ids(2)

    ids = self%getAttributableIds()
    CALL check(nf90_inquire_attribute(ids(1), ids(2), name, len=length), &
               "Could not inquire attribute "//name)
    CALL check(nf90_get_att(ids(1), ids(2), name, avalue), &
               "Could not read attribute "//name)

  END SUBROUTINE getAttributeI16

  SUBROUTINE getAttributeI32(self, name, avalue)
    CLASS(NcAttributable), INTENT(in)  :: self
    CHARACTER(*), INTENT(in)  :: name
    INTEGER(i32), INTENT(out) :: avalue
    INTEGER(i32)                        :: length, ids(2)

    ids = self%getAttributableIds()
    CALL check(nf90_inquire_attribute(ids(1), ids(2), name, len=length), &
               "Could not inquire attribute "//name)
    CALL check(nf90_get_att(ids(1), ids(2), name, avalue), &
               "Could not read attribute "//name)

  END SUBROUTINE getAttributeI32

  SUBROUTINE getAttributeF32(self, name, avalue)
    CLASS(NcAttributable), INTENT(in)  :: self
    CHARACTER(*), INTENT(in)  :: name
    REAL(f32), INTENT(out) :: avalue
    INTEGER(i32)                        :: length, ids(2)

    ids = self%getAttributableIds()
    CALL check(nf90_inquire_attribute(ids(1), ids(2), name, len=length), &
               "Could not inquire attribute "//name)
    CALL check(nf90_get_att(ids(1), ids(2), name, avalue), &
               "Could not read attribute "//name)

  END SUBROUTINE getAttributeF32

  SUBROUTINE getAttributeF64(self, name, avalue)
    CLASS(NcAttributable), INTENT(in)  :: self
    CHARACTER(*), INTENT(in)  :: name
    REAL(f64), INTENT(out) :: avalue
    INTEGER(i32)                        :: length, ids(2)

    ids = self%getAttributableIds()
    CALL check(nf90_inquire_attribute(ids(1), ids(2), name, len=length), &
               "Could not inquire attribute "//name)
    CALL check(nf90_get_att(ids(1), ids(2), name, avalue), &
               "Could not read attribute "//name)

  END SUBROUTINE getAttributeF64

  FUNCTION getAttributableIds(self)
    CLASS(NcAttributable), INTENT(in) :: self
    INTEGER(i32)                      :: getAttributableIds(2)
    SELECT TYPE (self)
    CLASS is (NcGroup)
      getAttributableIds(1) = self%id
      getAttributableIds(2) = NF90_GLOBAL
    CLASS is (NcVariable)
      getAttributableIds(1) = self%parent%id
      getAttributableIds(2) = self%id
    END SELECT
  END FUNCTION getAttributableIds

  SUBROUTINE renameAttribute(self, oldname, newname)
    CLASS(NcAttributable), INTENT(inout) :: self
    CHARACTER(len=*), INTENT(in)         :: oldname, newname
    INTEGER(i32)                         :: ids(2)
    ids = self%getAttributableIds()
    CALL check(nf90_rename_att(ids(1), ids(2), oldname, newname), "Failed to rename attribute: "//oldname)
  END SUBROUTINE renameAttribute

  SUBROUTINE renameVariable(self, name)
    CLASS(NcVariable), INTENT(inout) :: self
    CHARACTER(len=*), INTENT(in)    :: name
    CALL check(nf90_rename_var(self%parent%id, self%id, name), "Failed to rename variable: "//self%getName())
  END SUBROUTINE renameVariable

  SUBROUTINE renameDimension(self, name)
    CLASS(NcDimension), INTENT(inout) :: self
    CHARACTER(len=*), INTENT(in)     :: name
    CALL check(nf90_rename_dim(self%parent%id, self%id, name), "Failed to rename dimension: "//self%getName())
  END SUBROUTINE renameDimension

  SUBROUTINE setVariableFillValueI8(self, fvalue)
    CLASS(NcVariable), INTENT(inout)  :: self
    INTEGER(i8), INTENT(in)  :: fvalue

    IF (.NOT. self%hasAttribute("_FillValue")) THEN
      CALL self%setAttribute("_FillValue", fvalue)
    END IF
  END SUBROUTINE setVariableFillValueI8

  SUBROUTINE setVariableFillValueI16(self, fvalue)
    CLASS(NcVariable), INTENT(inout)  :: self
    INTEGER(i16), INTENT(in)  :: fvalue

    IF (.NOT. self%hasAttribute("_FillValue")) THEN
      CALL self%setAttribute("_FillValue", fvalue)
    END IF
  END SUBROUTINE setVariableFillValueI16

  SUBROUTINE setVariableFillValueI32(self, fvalue)
    CLASS(NcVariable), INTENT(inout)  :: self
    INTEGER(i32), INTENT(in)  :: fvalue

    IF (.NOT. self%hasAttribute("_FillValue")) THEN
      CALL self%setAttribute("_FillValue", fvalue)
    END IF
  END SUBROUTINE setVariableFillValueI32

  SUBROUTINE setVariableFillValueF32(self, fvalue)
    CLASS(NcVariable), INTENT(inout)  :: self
    REAL(f32), INTENT(in)  :: fvalue

    IF (.NOT. self%hasAttribute("_FillValue")) THEN
      CALL self%setAttribute("_FillValue", fvalue)
    END IF
  END SUBROUTINE setVariableFillValueF32

  SUBROUTINE setVariableFillValueF64(self, fvalue)
    CLASS(NcVariable), INTENT(inout)  :: self
    REAL(f64), INTENT(in)  :: fvalue

    IF (.NOT. self%hasAttribute("_FillValue")) THEN
      CALL self%setAttribute("_FillValue", fvalue)
    END IF
  END SUBROUTINE setVariableFillValueF64

  SUBROUTINE getVariableFillValueI8(self, fvalue)
    CLASS(NcVariable), INTENT(inout)  :: self
    INTEGER(i8), INTENT(out) :: fvalue

    IF (self%hasAttribute("_FillValue")) THEN
      CALL self%getAttribute("_FillValue", fvalue)
    ELSE
      fvalue = NF90_FILL_BYTE
    END IF
  END SUBROUTINE getVariableFillValueI8

  SUBROUTINE getVariableFillValueI16(self, fvalue)
    CLASS(NcVariable), INTENT(inout)  :: self
    INTEGER(i16), INTENT(out) :: fvalue

    IF (self%hasAttribute("_FillValue")) THEN
      CALL self%getAttribute("_FillValue", fvalue)
    ELSE
      fvalue = NF90_FILL_SHORT
    END IF
  END SUBROUTINE getVariableFillValueI16

  SUBROUTINE getVariableFillValueI32(self, fvalue)
    CLASS(NcVariable), INTENT(inout)  :: self
    INTEGER(i32), INTENT(out) :: fvalue

    IF (self%hasAttribute("_FillValue")) THEN
      CALL self%getAttribute("_FillValue", fvalue)
    ELSE
      fvalue = NF90_FILL_INT
    END IF
  END SUBROUTINE getVariableFillValueI32

  SUBROUTINE getVariableFillValueF32(self, fvalue)
    CLASS(NcVariable), INTENT(inout)  :: self
    REAL(f32), INTENT(out) :: fvalue

    IF (self%hasAttribute("_FillValue")) THEN
      CALL self%getAttribute("_FillValue", fvalue)
    ELSE
      fvalue = NF90_FILL_FLOAT
    END IF
  END SUBROUTINE getVariableFillValueF32

  SUBROUTINE getVariableFillValueF64(self, fvalue)
    CLASS(NcVariable), INTENT(inout)  :: self
    REAL(f64), INTENT(out) :: fvalue

    IF (self%hasAttribute("_FillValue")) THEN
      CALL self%getAttribute("_FillValue", fvalue)
    ELSE
      fvalue = NF90_FILL_DOUBLE
    END IF
  END SUBROUTINE getVariableFillValueF64

  SUBROUTINE setDataScalarI8(self, values, start)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i8), INTENT(in)           :: values
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setDataScalarI8

  SUBROUTINE setData1dI8(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i8), INTENT(in)           :: values(:)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData1dI8

  SUBROUTINE setData2dI8(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i8), INTENT(in)           :: values(:, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData2dI8

  SUBROUTINE setData3dI8(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i8), INTENT(in)           :: values(:, :, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData3dI8

  SUBROUTINE setData4dI8(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i8), INTENT(in)           :: values(:, :, :, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData4dI8

  SUBROUTINE setData5dI8(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i8), INTENT(in)           :: values(:, :, :, :, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData5dI8

  SUBROUTINE setDataScalarI16(self, values, start)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i16), INTENT(in)           :: values
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setDataScalarI16

  SUBROUTINE setData1dI16(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i16), INTENT(in)           :: values(:)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData1dI16

  SUBROUTINE setData2dI16(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i16), INTENT(in)           :: values(:, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData2dI16

  SUBROUTINE setData3dI16(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i16), INTENT(in)           :: values(:, :, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData3dI16

  SUBROUTINE setData4dI16(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i16), INTENT(in)           :: values(:, :, :, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData4dI16

  SUBROUTINE setData5dI16(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i16), INTENT(in)           :: values(:, :, :, :, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData5dI16

  SUBROUTINE setDataScalarI32(self, values, start)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i32), INTENT(in)           :: values
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setDataScalarI32

  SUBROUTINE setData1dI32(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i32), INTENT(in)           :: values(:)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData1dI32

  SUBROUTINE setData2dI32(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i32), INTENT(in)           :: values(:, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData2dI32

  SUBROUTINE setData3dI32(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i32), INTENT(in)           :: values(:, :, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData3dI32

  SUBROUTINE setData4dI32(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i32), INTENT(in)           :: values(:, :, :, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData4dI32

  SUBROUTINE setData5dI32(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i32), INTENT(in)           :: values(:, :, :, :, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData5dI32

  SUBROUTINE setDataScalarF32(self, values, start)
    CLASS(NcVariable), INTENT(in)           :: self
    REAL(f32), INTENT(in)           :: values
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setDataScalarF32

  SUBROUTINE setData1dF32(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    REAL(f32), INTENT(in)           :: values(:)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData1dF32

  SUBROUTINE setData2dF32(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    REAL(f32), INTENT(in)           :: values(:, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData2dF32

  SUBROUTINE setData3dF32(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    REAL(f32), INTENT(in)           :: values(:, :, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData3dF32

  SUBROUTINE setData4dF32(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    REAL(f32), INTENT(in)           :: values(:, :, :, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData4dF32

  SUBROUTINE setData5dF32(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    REAL(f32), INTENT(in)           :: values(:, :, :, :, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData5dF32

  SUBROUTINE setDataScalarF64(self, values, start)
    CLASS(NcVariable), INTENT(in)           :: self
    REAL(f64), INTENT(in)           :: values
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setDataScalarF64

  SUBROUTINE setData1dF64(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    REAL(f64), INTENT(in)           :: values(:)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData1dF64

  SUBROUTINE setData2dF64(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    REAL(f64), INTENT(in)           :: values(:, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData2dF64

  SUBROUTINE setData3dF64(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    REAL(f64), INTENT(in)           :: values(:, :, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData3dF64

  SUBROUTINE setData4dF64(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    REAL(f64), INTENT(in)           :: values(:, :, :, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData4dF64

  SUBROUTINE setData5dF64(self, values, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)           :: self
    REAL(f64), INTENT(in)           :: values(:, :, :, :, :)
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)

    CALL check(nf90_put_var(self%parent%id, self%id, values, start, cnt, stride, map), &
               "Failed to write data into variable: "//TRIM(self%getName()))
  END SUBROUTINE setData5dF64

  SUBROUTINE getDataScalarI8(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)            :: self
    INTEGER(i32), INTENT(in), OPTIONAL :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i8), INTENT(out)           :: DATA
    INTEGER(i8)                              :: tmp(1)

    CALL check(nf90_get_var(self%parent%id, self%id, tmp, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
    DATA = tmp(1)
  END SUBROUTINE getDataScalarI8

  SUBROUTINE getData1dI8(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i8), INTENT(out), ALLOCATABLE :: DATA(:)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData1dI8

  SUBROUTINE getData2dI8(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i8), INTENT(out), ALLOCATABLE :: DATA(:, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData2dI8

  SUBROUTINE getData3dI8(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i8), INTENT(out), ALLOCATABLE :: DATA(:, :, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2), datashape(3)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData3dI8

  SUBROUTINE getData4dI8(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i8), INTENT(out), ALLOCATABLE :: DATA(:, :, :, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2), datashape(3), datashape(4)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData4dI8

  SUBROUTINE getData5dI8(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i8), INTENT(out), ALLOCATABLE :: DATA(:, :, :, :, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2), datashape(3), datashape(4), datashape(5)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData5dI8

  SUBROUTINE getDataScalarI16(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i16), INTENT(out)              :: DATA
    INTEGER(i16)                                 :: tmp(1)

    CALL check(nf90_get_var(self%parent%id, self%id, tmp, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
    DATA = tmp(1)
  END SUBROUTINE getDataScalarI16

  SUBROUTINE getData1dI16(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i16), INTENT(out), ALLOCATABLE :: DATA(:)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData1dI16

  SUBROUTINE getData2dI16(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i16), INTENT(out), ALLOCATABLE :: DATA(:, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData2dI16

  SUBROUTINE getData3dI16(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i16), INTENT(out), ALLOCATABLE :: DATA(:, :, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2), datashape(3)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData3dI16

  SUBROUTINE getData4dI16(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i16), INTENT(out), ALLOCATABLE :: DATA(:, :, :, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2), datashape(3), datashape(4)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData4dI16

  SUBROUTINE getData5dI16(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i16), INTENT(out), ALLOCATABLE :: DATA(:, :, :, :, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2), datashape(3), datashape(4), datashape(5)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData5dI16

  SUBROUTINE getDataScalarI32(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i32), INTENT(out)              :: DATA
    INTEGER(i32)                                :: tmp(1)

    CALL check(nf90_get_var(self%parent%id, self%id, tmp, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
    DATA = tmp(1)
  END SUBROUTINE getDataScalarI32

  SUBROUTINE getData1dI32(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i32), INTENT(out), ALLOCATABLE :: DATA(:)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData1dI32

  SUBROUTINE getData2dI32(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i32), INTENT(out), ALLOCATABLE :: DATA(:, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData2dI32

  SUBROUTINE getData3dI32(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i32), INTENT(out), ALLOCATABLE :: DATA(:, :, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2), datashape(3)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData3dI32

  SUBROUTINE getData4dI32(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i32), INTENT(out), ALLOCATABLE :: DATA(:, :, :, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2), datashape(3), datashape(4)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData4dI32

  SUBROUTINE getData5dI32(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    INTEGER(i32), INTENT(out), ALLOCATABLE :: DATA(:, :, :, :, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2), datashape(3), datashape(4), datashape(5)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData5dI32

  SUBROUTINE getDataScalarF32(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)             :: self
    INTEGER(i32), INTENT(in), OPTIONAL  :: start(:), cnt(:), stride(:), map(:)
    REAL(f32), INTENT(out)            :: DATA
    REAL(f32)                                 :: tmp(1)

    CALL check(nf90_get_var(self%parent%id, self%id, tmp, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
    DATA = tmp(1)
  END SUBROUTINE getDataScalarF32

  SUBROUTINE getData1dF32(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    REAL(f32), INTENT(out), ALLOCATABLE :: DATA(:)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData1dF32

  SUBROUTINE getData2dF32(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    REAL(f32), INTENT(out), ALLOCATABLE :: DATA(:, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData2dF32

  SUBROUTINE getData3dF32(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    REAL(f32), INTENT(out), ALLOCATABLE :: DATA(:, :, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2), datashape(3)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData3dF32

  SUBROUTINE getData4dF32(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    REAL(f32), INTENT(out), ALLOCATABLE :: DATA(:, :, :, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2), datashape(3), datashape(4)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData4dF32

  SUBROUTINE getData5dF32(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    REAL(f32), INTENT(out), ALLOCATABLE :: DATA(:, :, :, :, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2), datashape(3), datashape(4), datashape(5)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData5dF32

  SUBROUTINE getDataScalarF64(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)             :: self
    INTEGER(i32), INTENT(in), OPTIONAL  :: start(:), cnt(:), stride(:), map(:)
    REAL(f64), INTENT(out)            :: DATA
    REAL(f64)                                 :: tmp(1)

    CALL check(nf90_get_var(self%parent%id, self%id, tmp, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
    DATA = tmp(1)
  END SUBROUTINE getDataScalarF64

  SUBROUTINE getData1dF64(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    REAL(f64), INTENT(out), ALLOCATABLE :: DATA(:)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData1dF64

  SUBROUTINE getData2dF64(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    REAL(f64), INTENT(out), ALLOCATABLE :: DATA(:, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData2dF64

  SUBROUTINE getData3dF64(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    REAL(f64), INTENT(out), ALLOCATABLE :: DATA(:, :, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2), datashape(3)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData3dF64

  SUBROUTINE getData4dF64(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    REAL(f64), INTENT(out), ALLOCATABLE :: DATA(:, :, :, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2), datashape(3), datashape(4)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData4dF64

  SUBROUTINE getData5dF64(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)               :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    REAL(f64), INTENT(out), ALLOCATABLE :: DATA(:, :, :, :, :)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    slcshape = self%getSlicingShape(start, cnt, stride)
    datashape = getReadShape(slcshape, SIZE(SHAPE(DATA)))

    ALLOCATE (DATA(datashape(1), datashape(2), datashape(3), datashape(4), datashape(5)))
    CALL check(nf90_get_var(self%parent%id, self%id, DATA, start, cnt, stride, map), &
               "Could not read data from variable: "//TRIM(self%getName()))
  END SUBROUTINE getData5dF64

  ! get string vector subroutine, TS@0207
  SUBROUTINE getData1dStr(self, DATA, start, cnt, stride, map)
    CLASS(NcVariable), INTENT(in)              :: self
    INTEGER(i32), INTENT(in), OPTIONAL    :: start(:), cnt(:), stride(:), map(:)
    CHARACTER(len=20), INTENT(out), ALLOCATABLE :: DATA(:)
    INTEGER(i32), ALLOCATABLE :: slcshape(:), datashape(:)

    ! local variables
    INTEGER(i32) :: i, status, data_ndim
    INTEGER(i32), ALLOCATABLE :: data_dimids(:), data_dimlen(:)

    ! get number of data dimensions
    status = nf90_inquire_variable(self%parent%id, self%id, ndims=data_ndim)
    ! allocate dimension related vectors
    ALLOCATE (data_dimids(data_ndim))
    ALLOCATE (data_dimlen(data_ndim))
    ! get dimension id
    status = nf90_inquire_variable(self%parent%id, self%id, dimids=data_dimids)
    ! get length of each dimension
    DO i = 1, data_ndim
      status = nf90_inquire_dimension(self%parent%id, data_dimids(i), len=data_dimlen(i))
    END DO

    ! slcshape = self%getSlicingShape(start, cnt, stride)
    ! datashape = getReadShape(slcshape, size(shape(data)))

    ALLOCATE (DATA(data_dimlen(1)))

    DO i = 1, data_dimlen(1)
      CALL check(nf90_get_var(self%parent%id, self%id, DATA(i), start=(/i, 1/), count=(/1, data_dimlen(2)/)), &
                 "Could not read data from variable: "//TRIM(self%getName()))
    END DO
  END SUBROUTINE getData1dStr

  FUNCTION getSlicingShape(self, instart, incnt, instride) RESULT(out)
    CLASS(NcVariable), INTENT(in)           :: self
    INTEGER(i32), INTENT(in), OPTIONAL :: instart(:), incnt(:), instride(:)
    INTEGER(i32), ALLOCATABLE          :: out(:)

    out = self%getShape()

    IF (PRESENT(incnt)) THEN
      out(:SIZE(incnt)) = incnt
      ! out = incnt
    ELSE
      IF (PRESENT(instart)) THEN
        out(:SIZE(instart)) = out(:SIZE(instart)) - (instart - 1)
      END IF
      IF (PRESENT(instride)) THEN
        out(:SIZE(instride)) = out(:SIZE(instride)) / instride
      END IF
    END IF

  END FUNCTION getSlicingShape

  FUNCTION getReadShape(slcshape, outrank) RESULT(out)
    INTEGER(i32), INTENT(in)   :: slcshape(:)
    INTEGER(i32), INTENT(in)   :: outrank
    INTEGER(i32)               :: naxis
    INTEGER(i32), ALLOCATABLE  :: out(:)

    naxis = COUNT(slcshape .GT. 1)

    IF (ALL(slcshape .EQ. 1)) THEN
      ! return 1-element array
      out(:) = 1
    ELSE IF (SIZE(slcshape) .EQ. outrank) THEN
      ! sizes fit
      out = slcshape
    ELSE IF (naxis .EQ. outrank) THEN
      out = PACK(slcshape, slcshape .GT. 1)
      ! else if (naxis .lt. outrank) then
      ! would be nice...
    ELSE
      WRITE (*, *) "Given indices do not match output variable rank!"
      STOP 1
    END IF
  END FUNCTION getReadShape

  FUNCTION getDtypeFromString(dtype)
    INTEGER(i32)          :: getDtypeFromString
    CHARACTER(*)         :: dtype

    SELECT CASE (dtype)
    CASE ("f32")
      getDtypeFromString = NF90_FLOAT
    CASE ("f64")
      getDtypeFromString = NF90_DOUBLE
    CASE ("i8")
      getDtypeFromString = NF90_BYTE
    CASE ("i16")
      getDtypeFromString = NF90_SHORT
    CASE ("i32")
      getDtypeFromString = NF90_INT
    CASE default
      WRITE (*, *) "Datatype not understood: ", dtype
      STOP 1
    END SELECT
  END FUNCTION getDtypeFromString

  FUNCTION getDtypeFromInteger(dtype)
    CHARACTER(3) :: getDtypeFromInteger
    INTEGER(i32)  :: dtype

    SELECT CASE (dtype)
    CASE (NF90_FLOAT)
      getDtypeFromInteger = "f32"
    CASE (NF90_DOUBLE)
      getDtypeFromInteger = "f64"
    CASE (NF90_BYTE)
      getDtypeFromInteger = "i8"
    CASE (NF90_SHORT)
      getDtypeFromInteger = "i16"
    CASE (NF90_INT)
      getDtypeFromInteger = "i32"
    CASE default
      WRITE (*, *) "Datatype not understood: ", dtype
      STOP 1
    END SELECT
  END FUNCTION getDtypeFromInteger

  FUNCTION getCreationMode(cmode)
    CHARACTER(*), INTENT(in), OPTIONAL :: cmode
    INTEGER(i32)                       :: getCreationMode
    CHARACTER(256)                     :: mode

    IF (.NOT. (PRESENT(cmode))) THEN
      mode = "NETCDF4"
    ELSE
      mode = cmode
    END IF

    SELECT CASE (TRIM(mode))
    CASE ("NETCDF4")
      getCreationMode = NF90_NETCDF4
    CASE ("SHARE")
      getCreationMode = NF90_SHARE
    CASE ("CLASSIC")
      getCreationMode = NF90_CLASSIC_MODEL
    CASE ("HDF5")
      getCreationMode = NF90_HDF5
    CASE ("64BIT_OFFSET")
      getCreationMode = NF90_64BIT_OFFSET
    CASE default
      PRINT *, "Creation mode not understood: "//TRIM(mode)
      STOP 1
    END SELECT

  END FUNCTION getCreationMode

  SUBROUTINE check(status, msg)
    INTEGER(i32), INTENT(in) :: status
    CHARACTER(*), INTENT(in) :: msg

    IF (status .NE. NF90_NOERR) THEN
      WRITE (*, *) msg
      WRITE (*, *) nf90_strerror(status)
      STOP 1
    END IF
  END SUBROUTINE check

END MODULE mo_netcdf

