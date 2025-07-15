MODULE GetEns_m
  USE kinds_m
  USE YAMLRead_m
  USE State_m, ONLY: State_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE Widgets_m

  TYPE EnsData_t
    REAL(r_kind), ALLOCATABLE :: DATA(:, :, :, :)             ! vertical X horizontal X tSlots X ensNum
    CHARACTER(LEN=20)         ::  varName                     ! the name of one variable
    INTEGER(i_kind)           ::  ensNum                      ! the number of ensembles
    INTEGER(i_kind)           ::  numLevel, numHorz, tSlots   ! dimensions of one variable from DA domain
    INTEGER(i_kind)           ::  numC                        ! the number of control variables
  END TYPE EnsData_t

  TYPE Ens_t
    TYPE(EnsData_t), ALLOCATABLE :: EnsData(:)
    !TYPE(SingleGrid_t)           :: sg
    INTEGER(i_kind)              :: numVar       ! the number of variable
    LOGICAL                      :: ensBKDiag, & ! the flag of BEC temporary output file format
                                    dataReady
    REAL(r_kind)                 :: lozRadius
  END TYPE Ens_t

  INTERFACE Ens_t
    PROCEDURE :: constructor
  END INTERFACE

CONTAINS


  FUNCTION constructor(configFile, Xm) RESULT(this)
    TYPE(Ens_t) :: this

    CHARACTER(LEN=1024), INTENT(IN)     ::  configFile
    TYPE(State_t), INTENT(IN), OPTIONAL ::  Xm(:)
    CHARACTER(LEN=20), ALLOCATABLE      ::  varList(:)
    CHARACTER(LEN=3)                    ::  memid
    INTEGER(i_kind)                     ::  ensNum, i, j, z, nx, ny, nz, nt, varIndx, iens, ifile
    INTEGER(i_kind)                     ::  numLevel, numHorz, numTime, numVars
    REAL(r_kind), ALLOCATABLE           ::  valueGlob(:, :, :)

    ! get the number of ensembles and the names
    ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ensNum', ensNum)
    PRINT *, 'number of ensembles:', ensNum

    ifile = yaml_get_var(TRIM(configFile), 'modelState', 'varList', varList)
    numVars = SIZE(varList)
    this%numVar = numVars
    PRINT *, 'numVars, varList:', numVars, varList


    ! initialization
    print *, 'this%dataReady:', this%dataReady
    this%dataReady = .FALSE.
    print *, 'this%dataReady:', this%dataReady
    !this%sg = Xm(1)%sg

    ALLOCATE (this%EnsData(numVars))

    
    !this%sg = Xm(1)%sg
    ensNum = UBOUND(Xm, 1)
    numHorz = Xm(1)%sg%num_icell_global
    numTime = Xm(1)%sg%tSlots
    numLevel = Xm(1)%sg%vLevel

    DO i = 1, numVars
      varIndx = Xm(1)%getVarIdx(TRIM(varList(i)))
      print *, 'varIndx: ', varIndx

      IF (varIndx .EQ. 0) THEN
        PRINT *, "ERROR! ", TRIM(varList(i)), "is not in Xm. Please check the setting or IO code."
        STOP
      ELSE
        this%EnsData(i)%numLevel = numLevel
        this%EnsData(i)%numHorz  = numHorz
        this%EnsData(i)%tSlots   = numTime
        this%EnsData(i)%ensNum   = ensNum
        this%EnsData(i)%numC     = numLevel * numHorz * numTime
        this%EnsData(i)%varName  = TRIM(varList(i))

        !IF (Xm(1)%sg%isBaseProc())  THEN
        !  IF (ALLOCATED(this%EnsData(i)%DATA)) DEALLOCATE(this%EnsData(i)%DATA)
        !  IF (ALLOCATED(valueGlob))            DEALLOCATE (valueGlob)
        !  ALLOCATE( this%EnsData(i)%DATA(numLevel, numHorz, numTime, ensNum), &
        !             valueGlob(numLevel, numHorz, numTime) )
        !END IF

        !DO iens = 1, ensNum
        !  CALL Xm(1)%sg%aggrGridRealForFieldGrid(Xm(iens)%fields(varIndx)%DATA, valueGlob, &
        !                                        [numLevel, numHorz, numTime])
        !  IF (Xm(1)%sg%isBaseProc()) THEN
        !    this%EnsData(i)%DATA(:, :, :, iens) = valueGlob
        !  !this%EnsData(i)%DATA(:, :, :, iens) = Xm(iens)%fields(varIndx)%DATA
        !    valueGlob = 0.0D0
        !  END IF
        !END DO

      END IF
    END DO

    this%dataReady = .TRUE.

  END FUNCTION constructor

END MODULE GetEns_m
