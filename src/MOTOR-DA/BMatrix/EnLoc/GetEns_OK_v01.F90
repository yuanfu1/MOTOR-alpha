MODULE GetEns_m
  USE kinds_m
  USE NMLRead_m
  USE State_m, ONLY: State_t
  USE WRFIO_m, ONLY: WRFIO_t
  USE Widgets_m

  TYPE EnsData_t
    REAL(r_kind), ALLOCATABLE :: DATA(:, :, :)     ! vertical X horizontal X numEns
    CHARACTER(LEN=20) ::  varName                ! the name of one variable
    INTEGER(i_kind)   ::  numEns                 ! the number of ensembles
    ! INTEGER(i_kind)   ::  nx, ny, nz             ! dimensions of one variable from model domain
    INTEGER(i_kind)   ::  numLevel, numHorz      ! dimensions of one variable from DA domain
  END TYPE EnsData_t

  TYPE Ens_t
    TYPE(EnsData_t), ALLOCATABLE :: EnsData(:)
    INTEGER(i_kind)              :: numVar       ! the number of variable
    INTEGER(i_kind)              :: ensBKDiag    ! the flag of BEC temporary output file format
    LOGICAL :: dataReady
  END TYPE Ens_t

  INTERFACE Ens_t
    PROCEDURE :: constructor
  END INTERFACE

CONTAINS

  ! SUBROUTINE constructor(this, configFile)
  FUNCTION constructor(configFile, Xm) RESULT(this)
    TYPE(Ens_t) :: this

    CHARACTER(LEN=1024), INTENT(IN)  ::  configFile
    TYPE(State_t), INTENT(IN)  ::  Xm(:)

    TYPE(WRFIO_t), ALLOCATABLE  ::  WRFIO
    CHARACTER(LEN=1024)               ::  inFileDir
    CHARACTER(LEN=1024), ALLOCATABLE  ::  inFileNames(:)
    CHARACTER(LEN=20), ALLOCATABLE  ::  varList(:), varWRFList(:)
    CHARACTER(LEN=20)  ::  ensForm
    CHARACTER(LEN=3)   ::  memid
    INTEGER(i_kind)      ::  numEns, i, j, z, shape_in(3), nx, ny, nz, nt, varIndx, iens

    ! get the number of ensembles and the names
    CALL namelist_read(TRIM(configFile), "NumEns", numEns)
    PRINT *, 'number of ensembles:', numEns
    ALLOCATE (inFileNames(numEns))
    CALL GET_ENVIRONMENT_VARIABLE("INPUT_DIR", inFileDir)
    DO i = 1, numEns
      WRITE (memid, '(A1,I2.2)') "m", i
      inFileNames(i) = TRIM(inFileDir)//"/EnsData/wrfout_"//TRIM(memid)
    END DO

    CALL namelist_read(TRIM(configFile), "varList", varList)
    numVars = SIZE(varList)
    this%numVar = numVars
    PRINT *, 'numVars, varList:', numVars, varList

    CALL namelist_read(TRIM(configFile), "varWRFList", varWRFList)
    PRINT *, 'varWRFList:', varWRFList

    CALL namelist_read(TRIM(configFile), "ensForm", ensForm)
    PRINT *, 'Ensemble data for BEC calculation were from ', TRIM(ensForm)

    ! this%ensBKDiag = 0
    CALL namelist_read(TRIM(configFile), "ensBKDiag", this%ensBKDiag)
    PRINT *, 'ensBKDiag ', this%ensBKDiag

    ! initialization
    shape_in = 1
    this%dataReady = .FALSE.
    ALLOCATE (this%EnsData(numVars))

    IF (TRIM(ensForm) == 'model') THEN
      !!! This form was used for test, not applied in minimization.

      ! get the variables with original dimensions through 'WRFIO_t'
      ALLOCATE (WRFIO)
      WRFIO = WRFIO_t(inFileNames, varWRFList)
      nx = WRFIO%nx
      ny = WRFIO%ny
      nz = WRFIO%nz
      nt = WRFIO%nt

      ! t2m
      varIndx = getStrIndex(varList, 't2m')
      IF (varIndx > 0) THEN
        IF (ALLOCATED(WRFIO%t2m)) THEN
          ALLOCATE (this%EnsData(varIndx)%DATA(1, nx * ny, nt))
          this%EnsData(varIndx)%DATA(1, :, :) = RESHAPE(WRFIO%t2m, (/nx * ny, nt/))
          ! this%EnsData(varIndx)%nx = nx
          ! this%EnsData(varIndx)%ny = ny
          ! this%EnsData(varIndx)%nz = 1
          this%EnsData(varIndx)%numHorz = nx * ny
          this%EnsData(varIndx)%numLevel = 1
          this%EnsData(varIndx)%numEns = nt
          this%EnsData(varIndx)%varName = 't2m'
        ELSE
          PRINT *, "ERROR! 't2m' is not in WRFIO, inconsistent with the setting. Please check the setting."
          STOP
        END IF
      END IF

      ! u10m
      varIndx = getStrIndex(varList, 'u10m')
      IF (varIndx > 0) THEN
        IF (ALLOCATED(WRFIO%u10m)) THEN
          ALLOCATE (this%EnsData(varIndx)%DATA(1, nx * ny, nt))
          this%EnsData(varIndx)%DATA(1, :, :) = RESHAPE(WRFIO%u10m, (/nx * ny, nt/))
          ! this%EnsData(varIndx)%nx = nx
          ! this%EnsData(varIndx)%ny = ny
          ! this%EnsData(varIndx)%nz = 1
          this%EnsData(varIndx)%numHorz = nx * ny
          this%EnsData(varIndx)%numLevel = 1
          this%EnsData(varIndx)%numEns = nt
          this%EnsData(varIndx)%varName = 'u10m'
        ELSE
          PRINT *, "ERROR! 'u10m' is not in WRFIO, inconsistent with the setting. Please check the setting."
          STOP
        END IF
      END IF

      ! v10m
      varIndx = getStrIndex(varList, 'v10m')
      IF (varIndx > 0) THEN
        IF (ALLOCATED(WRFIO%v10m)) THEN
          ALLOCATE (this%EnsData(varIndx)%DATA(1, nx * ny, nt))
          this%EnsData(varIndx)%DATA(1, :, :) = RESHAPE(WRFIO%v10m, (/nx * ny, nt/))
          ! this%EnsData(varIndx)%nx = nx
          ! this%EnsData(varIndx)%ny = ny
          ! this%EnsData(varIndx)%nz = 1
          this%EnsData(varIndx)%numHorz = nx * ny
          this%EnsData(varIndx)%numLevel = 1
          this%EnsData(varIndx)%numEns = nt
          this%EnsData(varIndx)%varName = 'v10m'
        ELSE
          PRINT *, "ERROR! 'v10m' is not in WRFIO, inconsistent with the setting. Please check the setting."
          STOP
        END IF
      END IF

      ! temp
      varIndx = getStrIndex(varList, 'temp')
      IF (varIndx > 0) THEN
        IF (ALLOCATED(WRFIO%temp)) THEN
          ALLOCATE (this%EnsData(varIndx)%DATA(nz, nx * ny, nt))
          DO z = 1, nz
            this%EnsData(varIndx)%DATA(z, :, :) = RESHAPE(WRFIO%temp(:, :, z, :), (/nx * ny, nt/))
          END DO
          ! this%EnsData(varIndx)%nx = nx
          ! this%EnsData(varIndx)%ny = ny
          ! this%EnsData(varIndx)%nz = nz
          this%EnsData(varIndx)%numHorz = nx * ny
          this%EnsData(varIndx)%numLevel = nz
          this%EnsData(varIndx)%numEns = nt
          this%EnsData(varIndx)%varName = 'temp'
        ELSE
          PRINT *, "ERROR! 'temp' is not in WRFIO, inconsistent with the setting. Please check the setting."
          STOP
        END IF
      END IF

      ! uwnd
      varIndx = getStrIndex(varList, 'uwnd')
      IF (varIndx > 0) THEN
        IF (ALLOCATED(WRFIO%uwnd)) THEN
          ALLOCATE (this%EnsData(varIndx)%DATA(nz, nx * ny, nt))
          DO z = 1, nz
            this%EnsData(varIndx)%DATA(z, :, :) = RESHAPE(WRFIO%uwnd(:, :, z, :), (/nx * ny, nt/))
          END DO
          ! this%EnsData(varIndx)%nx = nx
          ! this%EnsData(varIndx)%ny = ny
          ! this%EnsData(varIndx)%nz = nz
          this%EnsData(varIndx)%numHorz = nx * ny
          this%EnsData(varIndx)%numLevel = nz
          this%EnsData(varIndx)%numEns = nt
          this%EnsData(varIndx)%varName = 'uwnd'
        ELSE
          PRINT *, "ERROR! 'uwnd' is not in WRFIO, inconsistent with the setting. Please check the setting."
          STOP
        END IF
      END IF

      ! vwnd
      varIndx = getStrIndex(varList, 'vwnd')
      IF (varIndx > 0) THEN
        IF (ALLOCATED(WRFIO%vwnd)) THEN
          ALLOCATE (this%EnsData(varIndx)%DATA(nz, nx * ny, nt))
          DO z = 1, nz
            this%EnsData(varIndx)%DATA(z, :, :) = RESHAPE(WRFIO%vwnd(:, :, z, :), (/nx * ny, nt/))
          END DO
          ! this%EnsData(varIndx)%nx = nx
          ! this%EnsData(varIndx)%ny = ny
          ! this%EnsData(varIndx)%nz = nz
          this%EnsData(varIndx)%numHorz = nx * ny
          this%EnsData(varIndx)%numLevel = nz
          this%EnsData(varIndx)%numEns = nt
          this%EnsData(varIndx)%varName = 'vwnd'
        ELSE
          PRINT *, "ERROR! 'vwnd' is not in WRFIO, inconsistent with the setting. Please check the setting."
          STOP
        END IF
      END IF

      ! pres
      varIndx = getStrIndex(varList, 'pres')
      IF (varIndx > 0) THEN
        IF (ALLOCATED(WRFIO%pres)) THEN
          ALLOCATE (this%EnsData(varIndx)%DATA(nz, nx * ny, nt))
          DO z = 1, nz
            this%EnsData(varIndx)%DATA(z, :, :) = RESHAPE(WRFIO%pres(:, :, z, :), (/nx * ny, nt/))
          END DO
          ! this%EnsData(varIndx)%nx = nx
          ! this%EnsData(varIndx)%ny = ny
          ! this%EnsData(varIndx)%nz = nz
          this%EnsData(varIndx)%numHorz = nx * ny
          this%EnsData(varIndx)%numLevel = nz
          this%EnsData(varIndx)%numEns = nt
          this%EnsData(varIndx)%varName = 'pres'
        ELSE
          PRINT *, "ERROR! 'pres' is not in WRFIO, inconsistent with the setting. Please check the setting."
          STOP
        END IF
      END IF

      ! psfc
      varIndx = getStrIndex(varList, 'psfc')
      IF (varIndx > 0) THEN
        IF (ALLOCATED(WRFIO%psfc)) THEN
          ALLOCATE (this%EnsData(varIndx)%DATA(1, nx * ny, nt))
          this%EnsData(varIndx)%DATA(1, :, :) = RESHAPE(WRFIO%psfc, (/nx * ny, nt/))
          ! this%EnsData(varIndx)%nx = nx
          ! this%EnsData(varIndx)%ny = ny
          ! this%EnsData(varIndx)%nz = 1
          this%EnsData(varIndx)%numHorz = nx * ny
          this%EnsData(varIndx)%numLevel = 1
          this%EnsData(varIndx)%numEns = nt
          this%EnsData(varIndx)%varName = 'psfc'
        ELSE
          PRINT *, "ERROR! 'psfc' is not in WRFIO, inconsistent with the setting. Please check the setting."
          STOP
        END IF
      END IF

      ! q2m
      varIndx = getStrIndex(varList, 'q2m')
      IF (varIndx > 0) THEN
        IF (ALLOCATED(WRFIO%q2m)) THEN
          ALLOCATE (this%EnsData(varIndx)%DATA(1, nx * ny, nt))
          this%EnsData(varIndx)%DATA(1, :, :) = RESHAPE(WRFIO%q2m, (/nx * ny, nt/))
          ! this%EnsData(varIndx)%nx = nx
          ! this%EnsData(varIndx)%ny = ny
          ! this%EnsData(varIndx)%nz = 1
          this%EnsData(varIndx)%numHorz = nx * ny
          this%EnsData(varIndx)%numLevel = 1
          this%EnsData(varIndx)%numEns = nt
          this%EnsData(varIndx)%varName = 'q2m'
        ELSE
          PRINT *, "ERROR! 'q2m' is not in WRFIO, inconsistent with the setting. Please check the setting."
          STOP
        END IF
      END IF

      ! qvapor
      varIndx = getStrIndex(varList, 'qvapor')
      IF (varIndx > 0) THEN
        IF (ALLOCATED(WRFIO%qvapor)) THEN
          ALLOCATE (this%EnsData(varIndx)%DATA(nz, nx * ny, nt))
          DO z = 1, nz
            this%EnsData(varIndx)%DATA(z, :, :) = RESHAPE(WRFIO%qvapor(:, :, z, :), (/nx * ny, nt/))
          END DO
          ! this%EnsData(varIndx)%nx = nx
          ! this%EnsData(varIndx)%ny = ny
          ! this%EnsData(varIndx)%nz = nz
          this%EnsData(varIndx)%numHorz = nx * ny
          this%EnsData(varIndx)%numLevel = nz
          this%EnsData(varIndx)%numEns = nt
          this%EnsData(varIndx)%varName = 'qvapor'
        ELSE
          PRINT *, "ERROR! 'qvapor' is not in WRFIO, inconsistent with the setting. Please check the setting."
          STOP
        END IF
      END IF

      this%dataReady = .TRUE.
      DEALLOCATE (WRFIO)

    ELSE IF (TRIM(ensForm) == 'state') THEN
      !!! This form was used in minimization, get the data from State_t.

      numEns = UBOUND(Xm, 1)
      DO i = 1, numVars
        varIndx = Xm(1)%getVarIdx(TRIM(varList(i)))
        IF (varIndx .NE. 0) THEN
          shape_in = SHAPE(Xm(1)%Fields(varIndx)%DATA)
          ALLOCATE (this%EnsData(i)%DATA(shape_in(1), shape_in(2), numEns))
          this%EnsData(i)%numLevel = shape_in(1)
          this%EnsData(i)%numHorz = shape_in(2)
          this%EnsData(i)%numEns = numEns
          this%EnsData(i)%varName = TRIM(varList(i))
          DO iens = 1, numEns
            this%EnsData(i)%DATA(:, :, iens) = Xm(iens)%Fields(varIndx)%DATA(:, :, 1)
          END DO
        ELSE
          PRINT *, "ERROR! ", TRIM(varList(i)), "is not in Xm. Please check the setting or IO code."
          STOP
        END IF
      END DO

      this%dataReady = .TRUE.

    ELSE
      PRINT *, "ERROR! MOTOR-DA doesn't know where the sample data from, please check 'ensForm' in setting file."
      STOP
    END IF

  END FUNCTION constructor

END MODULE GetEns_m
