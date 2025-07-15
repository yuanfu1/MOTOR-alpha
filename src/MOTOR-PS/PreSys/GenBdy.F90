MODULE GenBdy_m
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE YAMLRead_m
  IMPLICIT NONE

  TYPE :: GenBdy_t
    INTEGER(i_kind) :: buffer
  CONTAINS
    PROCEDURE :: GenBdyIdx
    PROCEDURE :: GenBdyBufferWeight
    PROCEDURE :: GenBdyTendency
  END TYPE

  INTERFACE GenBdy_t
    PROCEDURE :: constructor
  END INTERFACE GenBdy_t
CONTAINS

  FUNCTION constructor(configFile) RESULT(this)
    IMPLICIT NONE

    TYPE(GenBdy_t) :: this
    CHARACTER(*), INTENT(IN) :: configFile
    INTEGER(i_kind) :: istatus

    istatus = yaml_get_var(TRIM(configFile), 'TimeIntegral', 'buffer', this%buffer)

    IF (istatus .NE. 0) this%buffer = 5

    this%buffer = this%buffer * (-1)

  END FUNCTION

  SUBROUTINE GenBdyIdx(this, sg)
    IMPLICIT NONE

    CLASS(GenBdy_t) :: this
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    INTEGER(i_kind) :: i, idx1, idx2
    INTEGER(i_kind), ALLOCATABLE :: BdyIdxTemp(:)

    ALLOCATE (BdyIdxTemp(sg%num_icell_global), sg%bdy_type(sg%num_cell))

    !   CALL sg%aggrGridInt(sg%cell_type, CellTypeTemp, sg%num_icell_global)
    !   CALL sg%aggrGridInt(sg%cell_stcl, CellStclTemp, [2, sg%num_icell_global])

    !   BdyIdxTemp = CellTypeTemp

    PRINT *, 'The global Latlon gridpoints are ',sg%dimCell_global,sg%gLevel

    DO i = 1, sg%num_icell_global
      idx1 = i / sg%dimCell_global(2) + 1
      idx2 = MOD(i, sg%dimCell_global(2))
      IF (idx2 .EQ. 0) THEN
        idx2 = sg%dimCell_global(2)
        idx1 = idx1 - 1
      END IF
      IF (idx1 .GT. sg%dimCell_global(1) / 2) THEN
        idx1 = sg%dimCell_global(1) - idx1 + 1
      END IF
      IF (idx2 .GT. sg%dimCell_global(2) / 2) THEN
        idx2 = sg%dimCell_global(2) - idx2 + 1
      END IF

      IF (idx1 .LE. idx2) THEN
        BdyIdxTemp(i) = -(idx1 - 2)
      ELSE
        BdyIdxTemp(i) = -(idx2 - 2)
      END IF
    END DO

    CALL sg%DistGridInt1D(BdyIdxTemp, sg%bdy_type, (/sg%num_icell_global/))
    ! CALL sg%ExchangeMatOnHalo2D(1, sg%bdy_type)
  END SUBROUTINE

  SUBROUTINE GenBdyBufferWeight(this, sg)
    IMPLICIT NONE

    CLASS(GenBdy_t) :: this
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    INTEGER(i_kind) :: i, idx1, idx2
    INTEGER(i_kind), ALLOCATABLE :: temp(:)
    ! REAL(r_kind), ALLOCATABLE :: BdyAlphaTemp(:), BdyBetaTemp(:)

    IF (.NOT. ALLOCATED(sg%bdy_type)) THEN
      PRINT *, '#==============================='
      PRINT *, 'ERROR, sg%bdy_type should be set before calculating BufferWeight'
      PRINT *, 'Use GenBdYIdx in GenBdy module to set sg%bdy'
      PRINT *, '#==============================='
      STOP
    END IF

    ALLOCATE (sg%BdyAlpha(sg%num_cell), &
              sg%BdyBeta(sg%num_cell))

    temp = PACK([(i, i=1, SIZE(sg%bdy_type))], sg%bdy_type >= this%buffer)
    ! temp = PACK(sg%bdy_type, sg%bdy_type >= this%buffer)
    IF (SIZE(temp) .GE. 1) THEN
      ALLOCATE (sg%BufferHalo(SIZE(temp)))
      sg%BufferHalo = temp
    END IF

    DO i = 1, sg%num_cell
      IF (sg%bdy_type(i) .GE. 0) THEN
        sg%BdyAlpha(i) = 1.0D0
        sg%BdyBeta(i) = 0.0D0
      ELSEIF (sg%bdy_type(i) .GE. this%buffer) THEN
        sg%BdyAlpha(i) = 1.0D0 - DBLE(sg%bdy_type(i)) / DBLE(this%buffer)
        sg%BdyBeta(i) = 1.0D0 - sg%BdyAlpha(i)
      ELSE
        sg%BdyAlpha(i) = 0.0D0
        sg%BdyBeta(i) = 1.0D0
      END IF
    END DO

    DEALLOCATE (temp)
  END SUBROUTINE

  SUBROUTINE GenBdyTendency(this, XFull, XHalf)
    CLASS(GenBdy_t) :: this
    TYPE(State_t), INTENT(INOUT) :: XFull, XHalf
    INTEGER(i_kind) :: i, st, sHalf, sFull
    REAL(r_kind) :: dt
    CHARACTER(6) :: varNamesHalf(4) = ['vor   ', 'div   ', 'rho   ', 'qvapor'], &
                    varNamesFull(2) = ['w     ', 'theta ']
    CHARACTER(11) :: TenVarNamesHalf(4) = ['BdyTenVor  ', 'BdyTenDiv  ', 'BdyTenRho  ', &
                                           'BdyTenq    '], &
                     TenVarNamesFull(2) = ['BdyTenW    ', 'BdyTenTheta']

    sHalf = SIZE(varNamesHalf)
    sFull = SIZE(varNamesFull)

    DO i = 1, sHalf
      CALL XHalf%addVar(TRIM(TenVarNamesHalf(i)))
    END DO

    DO i = 1, sFull
      CALL XFull%addVar(TRIM(TenVarNamesFull(i)))
    END DO

    dt = XHalf%sg%tt(2) - XHalf%sg%tt(1)
    st = SIZE(XHalf%fields(XHalf%getVarIdx(TRIM(varNamesHalf(1))))%DATA, 3)

    DO i = 1, sHalf
      XHalf%fields(XHalf%getVarIdx(TRIM(TenVarNamesHalf(i))))%DATA(:, :, 1:st - 1) = &
        (XHalf%fields(XHalf%getVarIdx(TRIM(varNamesHalf(i))))%DATA(:, :, 2:st) - &
         XHalf%fields(XHalf%getVarIdx(TRIM(varNamesHalf(i))))%DATA(:, :, 1:st - 1)) / dt
    END DO

    DO i = 1, sFull
      XFull%fields(XFull%getVarIdx(TRIM(TenVarNamesFull(i))))%DATA(:, :, 1:st - 1) = &
        (XFull%fields(XFull%getVarIdx(TRIM(varNamesFull(i))))%DATA(:, :, 2:st) - &
         XFull%fields(XFull%getVarIdx(TRIM(varNamesFull(i))))%DATA(:, :, 1:st - 1)) / dt
    END DO

  END SUBROUTINE

END MODULE
