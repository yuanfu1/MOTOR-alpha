MODULE nwp_m
  USE TenDiv_m, ONLY: TenDiv_t
  USE TenVor_m, ONLY: TenVor_t
  USE TenRho_m, ONLY: TenRho_t
  USE TenTheta_m, ONLY: TenTheta_t
  USE TenW_m, ONLY: Tenw_t
  USE Tenqvapor_m, ONLY: Tenqvapor_t
  USE CalPres_m, ONLY: CalPres_t
  USE State_m, ONLY: State_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE PreCal_m, ONLY: PreCal_t
  USE GenBdy_m, ONLY: GenBdy_t

  IMPLICIT NONE

  TYPE :: nwp_t

    TYPE(TenDiv_t) :: TenDiv
    TYPE(TenVor_t) :: TenVor
    TYPE(TenTheta_t) :: TenTheta
    TYPE(TenRho_t) :: TenRho
    TYPE(TenW_t) :: Tenw
    TYPE(Tenqvapor_t) :: Tenqvapor
    TYPE(CalPres_t) :: CalPres
    TYPE(SingleGrid_t), POINTER :: sgHalf, sgFull
    TYPE(precal_t) :: precal
    TYPE(GenBdy_t) :: GenBdy

  CONTAINS
    FINAL :: destructor
    PROCEDURE :: righthands
    PROCEDURE :: UpdateBdy
    PROCEDURE :: destroy

  END TYPE nwp_t

  INTERFACE nwp_t
    PROCEDURE :: constructor
  END INTERFACE

CONTAINS

  FUNCTION constructor(sgFull, sgHalf) RESULT(this)

    IMPLICIT NONE
    TYPE(nwp_t) :: this
    TYPE(SingleGrid_t), TARGET :: sgFull, sgHalf

    this%sgFull => sgFull
    this%sgHalf => sgHalf

    this%PreCal = PreCal_t(sgHalf) ! div, vor, rho and q are on half layers

    this%TenDiv = TenDiv_t(this%sgHalf)
    this%TenVor = TenVor_t(this%sgHalf)
    this%TenQvapor = TenQvapor_t(this%sgHalf)
    this%TenRho = TenRho_t(this%sgHalf)
    this%CalPres = CalPres_t(this%sgHalf)
    this%TenW = TenW_t(this%sgFull)
    this%TenTheta = TenTheta_t(this%sgFull)

  END FUNCTION constructor
  SUBROUTINE righthands(this, XFull, XHalf)
    IMPLICIT NONE

    CLASS(nwp_t) :: this
    TYPE(State_t), INTENT(INOUT) :: XFull, XHalf
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: div, vor, psi, chi, &
                                                  w, theta, qvapor, pres, rho, &
                                                  ten_div, ten_vor, &
                                                  ten_w, ten_theta, ten_q, ten_rho
    REAL(r_kind) :: start_time_t, end_time_t, start_div, end_div, start_vor, end_vor, &
                    start_time_al, end_time_al

    ASSOCIATE (vLevel1 => this%sgHalf%vLevel, &
               vLevel2 => this%sgFull%vLevel, &
               num_cell => this%sgFull%num_cell)

      CALL CPU_TIME(start_time_al)

      ALLOCATE (div(vLevel1, num_cell), vor(vLevel1, num_cell), &
                psi(vLevel1, num_cell), chi(vLevel1, num_cell), &
                qvapor(vLevel1, num_cell), pres(vLevel1, num_cell), &
                rho(vLevel1, num_cell), &
                w(vLevel2, num_cell), theta(vLevel1, num_cell), &
                ten_div(vLevel1, num_cell), ten_vor(vLevel1, num_cell), &
                ten_q(vLevel1, num_cell), ten_rho(vLevel1, num_cell), &
                ten_w(vLevel2, num_cell), ten_theta(vLevel2, num_cell))
      CALL CPU_TIME(end_time_al)

      PRINT *, 'allocate in nwp takes time: ', end_time_al - start_time_al, 'S'

      div = XHalf%Fields(XHalf%getVarIdx('div'))%DATA(:, :, 1)
      vor = XHalf%Fields(XHalf%getVarIdx('vor'))%DATA(:, :, 1)
      psi = XHalf%Fields(XHalf%getVarIdx('psi'))%DATA(:, :, 1)
      chi = XHalf%Fields(XHalf%getVarIdx('chi'))%DATA(:, :, 1)
      rho = XHalf%Fields(XHalf%getVarIdx('rho'))%DATA(:, :, 1)
      w = XFull%Fields(XHalf%getVarIdx('w'))%DATA(:, :, 1)
      theta = XFull%Fields(XHalf%getVarIdx('theta'))%DATA(:, :, 1)
      qvapor = XHalf%Fields(XHalf%getVarIdx('qvapor'))%DATA(:, :, 1)
      pres = XHalf%Fields(XHalf%getVarIdx('pres'))%DATA(:, :, 1)
      ! pres = 0.0D0

      IF (XHalf%getVarIdx('tendiv') .EQ. 0) THEN
        CALL XHalf%addVar('tendiv')
        CALL XHalf%addVar('tenvor')
        CALL XHalf%addVar('tenqvapor')
        CALL XHalf%addVar('tenrho')
        CALL XFull%addVar('tenw')
        CALL XFull%addVar('tentheta')
      END IF

      CALL this%precal%PreCals(XHalf)

      ten_div = 0.0D0
      ten_vor = 0.0D0
      ten_rho = 0.0D0
      ten_w = 0.0D0
      ten_theta = 0.0D0
      ten_q = 0.0D0

      ! CALL this%CalPres%CalPres(rho, theta, pres)
      ! CALL this%sgHalf%ExchangeMatOnHalo2D(this%sgHalf%vLevel, pres)
      CALL CPU_TIME(start_div)
      CALL this%TenDiv%Tendency(div, vor, psi, chi, w, rho, pres, this%precal, ten_div)
      CALL CPU_TIME(end_div)
      CALL CPU_TIME(start_vor)
      CALL this%TenVor%Tendency(vor, psi, chi, w, rho, pres, this%precal, ten_vor)
      CALL CPU_TIME(end_vor)
      PRINT *, 'div in nwp takes time: ', end_div - start_div, 'S'
      PRINT *, 'vor in nwp takes time: ', end_vor - start_vor, 'S'
      PRINT *, 'div and vor in nwp takes time: ', start_vor - end_div, 'S'

      CALL this%TenRho%Tendency(rho, div, psi, chi, w, this%precal, ten_rho)
      CALL this%TenW%Tendency(w, div, psi, chi, rho, pres, ten_w)
      CALL this%TenTheta%Tendency(theta, div, psi, chi, w, ten_theta)
      CALL this%Tenqvapor%Tendency(qvapor, div, psi, chi, w, ten_q)

      CALL this%sgHalf%ExchangeMatOnHalo2D(this%sgHalf%vLevel, ten_div)
      CALL this%sgHalf%ExchangeMatOnHalo2D(this%sgHalf%vLevel, ten_vor)
      CALL this%sgHalf%ExchangeMatOnHalo2D(this%sgHalf%vLevel, ten_rho)
      CALL this%sgHalf%ExchangeMatOnHalo2D(this%sgHalf%vLevel, ten_q)
      CALL this%sgFull%ExchangeMatOnHalo2D(this%sgFull%vLevel, ten_w)
      CALL this%sgFull%ExchangeMatOnHalo2D(this%sgFull%vLevel, ten_theta)

      XHalf%Fields(XHalf%getVarIdx('tendiv'))%DATA(:, :, 1) = ten_div
      XHalf%Fields(XHalf%getVarIdx('tenvor'))%DATA(:, :, 1) = ten_vor
      XHalf%Fields(XHalf%getVarIdx('tenrho'))%DATA(:, :, 1) = ten_rho
      XHalf%Fields(XHalf%getVarIdx('tenqvapor'))%DATA(:, :, 1) = ten_q
      XHalf%Fields(XHalf%getVarIdx('pres'))%DATA(:, :, 1) = pres
      XFull%Fields(XFull%getVarIdx('tenw'))%DATA(:, :, 1) = ten_w
      XFull%Fields(XFull%getVarIdx('tentheta'))%DATA(:, :, 1) = ten_theta

    END ASSOCIATE

    DEALLOCATE (div, vor, psi, chi, &
                w, theta, qvapor, pres, rho, &
                ten_div, ten_vor, &
                ten_w, ten_theta, ten_q, ten_rho)

  END SUBROUTINE

  SUBROUTINE UpdateBdy(this, Xin, XBdy, k, layertype)
    CLASS(nwp_t) :: this
    TYPE(State_t), INTENT(IN) :: XBdy
    TYPE(State_t), INTENT(INOUT) :: Xin
    INTEGER(i_kind), INTENT(IN) :: k
    CHARACTER(4), INTENT(IN) :: layertype
    INTEGER(i_kind) :: i, j, sHalf, sFull
    CHARACTER(9) :: varNamesHalf(4) = ['tenvor   ', 'tendiv   ', 'tenrho   ', 'tenqvapor'], &
                    varNamesFull(2) = ['tenw     ', 'tentheta ']
    CHARACTER(11) :: TenVarNamesHalf(4) = ['BdyTenVor  ', 'BdyTenDiv  ', 'BdyTenRho  ', &
                                           'BdyTenq    '], &
                     TenVarNamesFull(2) = ['BdyTenW    ', 'BdyTenTheta']

    sHalf = SIZE(varNamesHalf)
    sFull = SIZE(varNamesFull)

    ASSOCIATE (sg => Xin%sg)

      IF ((.NOT. ALLOCATED(sg%bdy_type)) .OR. (.NOT. ALLOCATED(sg%BdyAlpha))) THEN
        PRINT *, '#==============================='
        PRINT *, 'ERROR, sg%bdy_type should be set and BufferWeight should be calculated before UpdateBdy'
        PRINT *, 'Use GenBdYIdx in GenBdy module to set sg%bdy_type'
        PRINT *, 'Use GenBdyBufferWeight in GenBdy module to set sg%BdyAlpha and sg%BdyBeta'
        PRINT *, '#==============================='
        STOP
      END IF

      IF (.NOT. ALLOCATED(sg%BufferHalo)) RETURN

      IF (layertype .EQ. 'HALF') THEN
        DO i = 1, sHalf
          DO j = 1, sg%vLevel
            Xin%Fields(Xin%getVarIdx(TRIM(varNamesHalf(i))))%DATA(j, sg%BufferHalo, 1) = &
              Xin%Fields(Xin%getVarIdx(TRIM(varNamesHalf(i))))%DATA(j, sg%BufferHalo, 1) * sg%BdyBeta(sg%BufferHalo) &
              + XBdy%Fields(XBdy%getVarIdx(TRIM(TenVarNamesHalf(i))))%DATA(j, sg%BufferHalo, k) * sg%BdyAlpha(sg%BufferHalo)
          END DO
        END DO
      ELSEIF (layertype .EQ. 'FULL') THEN
        DO i = 1, sFull
          DO j = 1, sg%vLevel
            Xin%Fields(Xin%getVarIdx(TRIM(varNamesFull(i))))%DATA(j, sg%BufferHalo, 1) = &
              Xin%Fields(Xin%getVarIdx(TRIM(varNamesFull(i))))%DATA(j, sg%BufferHalo, 1) * sg%BdyBeta(sg%BufferHalo) &
              + XBdy%Fields(XBdy%getVarIdx(TRIM(TenVarNamesFull(i))))%DATA(j, sg%BufferHalo, k) * sg%BdyAlpha(sg%BufferHalo)
          END DO
        END DO
      ELSE
        PRINT *, 'ERROR: Layer type should be FULL or HALF'
        STOP
      END IF
    END ASSOCIATE

  END SUBROUTINE

  SUBROUTINE destroy(this)
    IMPLICIT NONE
    CLASS(nwp_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sgFull)) NULLIFY (this%sgFull)
    IF (ASSOCIATED(this%sgHalf)) NULLIFY (this%sgHalf)

    CALL this%PreCal%destroy()
    CALL this%TenDiv%destroy()
    CALL this%TenVor%destroy()
    CALL this%TenQvapor%destroy()
    CALL this%TenRho%destroy()
    CALL this%CalPres%destroy()
    CALL this%TenW%destroy()
    CALL this%TenTheta%destroy()

  END SUBROUTINE destroy

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(nwp_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sgFull)) NULLIFY (this%sgFull)
    IF (ASSOCIATED(this%sgHalf)) NULLIFY (this%sgHalf)

  END SUBROUTINE destructor

END MODULE nwp_m
