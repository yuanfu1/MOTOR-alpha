!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ObsBase
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2021-09-01   Created by Yuanfu Xie
!                     2022-01-19   Added sub SuperObs to the ABSTRACT Data Type by Jiongming Pang
!
!   Created by Yuanfu Xie (xieyf@gbamwf.com), 2021/9/01, @GBA-MWF, Shenzhen
!   Modified by Jiongming Pang (pang.j.m@hotmail.com), 2022/01/19, @GBA-MWF, Shenzhen
!   Modified by Zilong Qin (zilong.qin@gmail.com), 2022/2/27, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/6/3, for changing thinning weighting func
!               using the sum of Gaussian weight to rescale thinning instead of simple average.
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/6/21, for adding a new retriction of thinned
!               observations from a finer grid to a coarser grid, griddedObsMap.
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/9/16, for modifying GetForwardValue to
!     return 0 background if the obs variables are not model states, since this operator defined
!     in this module is for convention obs only. All other obs will have their own GetForwardValue.
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2023/10/27, for adding landmask, topo and bkgd
!     flow dependence to thinning.
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/05/30, for extending M2O, the observation
!     forward operators. M is predetermined model states, that Control and background are required
!     to convert to, in C2M and background ingest.
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides an abstract data structure of general observation data.
MODULE ObsBase_m
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE mpObs_m, ONLY: mpObs_t
  USE geoTools_m, ONLY: GeoBox_t
  USE slint, ONLY: slint_init, tgt_grid
  USE ObsField_m, ONLY: ObsField_t
  USE AdvanceTime_m
  USE domainCheck_m, ONLY: domainCheck_t
  USE parameters_m, ONLY: degree2radian, radian2degree, missing, invalid, machineEps
  USE parameters_m, ONLY: spec_heat_const_pres, dry_air_gas_const, surface_ref_pres
  USE obsCommon_m, ONLY: obsCommon_t
  USE ObsRaw_m, ONLY: ObsRaw_t

  IMPLICIT NONE

  TYPE, ABSTRACT :: ObsBase_t
    CHARACTER(LEN=40) :: obsType
    REAL(r_kind) :: correlation_threshold
    REAL(r_kind), ALLOCATABLE :: obsData(:, :), &      ! First: data; second: obs elements e.g. (u,v)
                                 obsErrs(:, :)         ! observation errors
    REAL(r_kind), ALLOCATABLE :: olatlon(:, :), &      ! obsSlot: observation latlon location
                                 obsHght(:)
    INTEGER(i_kind), ALLOCATABLE :: obsTime(:)
    CHARACTER(LEN=1024) :: configFile

    INTEGER(i_kind) :: numObs, numVars, numFiles       ! Number of obs, vars and obs files
    INTEGER(i_kind) :: interpolation                   ! 1: using slint for unstructed grid; 2: latlon grid interpolation
    CHARACTER(LEN=1024), ALLOCATABLE :: fileNames(:)   ! Obs file names
    CHARACTER(LEN=10), ALLOCATABLE :: stNames(:)       ! Obs station names
    CHARACTER(LEN=20), ALLOCATABLE :: varNames(:)      ! Obs variable names
    CHARACTER(LEN=20), ALLOCATABLE :: types(:)         ! obsType_varName
    REAL(r_kind), ALLOCATABLE :: radius(:, :)          ! Influence radius in 4 directions of each variable
    REAL(r_kind), ALLOCATABLE :: land(:)               ! Topography and land at obs site
    REAL(r_kind), ALLOCATABLE :: sizeInc(:)            ! Size of increment in 4 directions of each
    REAL(r_kind), ALLOCATABLE :: qcThreshold(:)        ! Threshold values of QC for each variable
    INTEGER(r_kind) :: presIdx = 0
    INTEGER(r_kind) :: qvaporIdx = 0
    INTEGER(r_kind) :: tempIdx = 0                     ! The pressure and temperature index of vars in obs
    INTEGER(r_kind) :: lnpIdx = 0                      ! The pressure and temperature index of vars in obs

    TYPE(obsCommon_t) :: obsCommon

    ! SoloSite information for holding off obs from analysis:
    CHARACTER(LEN=5), ALLOCATABLE :: soloNames(:)
    REAL(r_kind), ALLOCATABLE :: soloLat(:), soloLon(:)

  CONTAINS
    ! PROCEDURE, PUBLIC :: griddedObsMap! Added by Xie 2022-06-21
    PROCEDURE, PUBLIC :: getVarIdx! Added by Xie 2022-08-02
    ! PROCEDURE, PUBLIC :: ObsMinusState! Added by Xie 2022-06-11
    PROCEDURE, PUBLIC :: ObsThinning  ! Modified ObsSuper by Xie
    PROCEDURE, PUBLIC :: ObsDestroy   ! destroy arrays of observations
    PROCEDURE(oneAug), DEFERRED :: ObsInitial
    PROCEDURE(ingest), DEFERRED :: ObsIngest
    PROCEDURE(forward), DEFERRED :: ObsForward
    PROCEDURE(tangent), DEFERRED :: ObsTangent
    PROCEDURE(Adjoint), DEFERRED :: ObsAdjoint
    ! PROCEDURE(noAugs), DEFERRED :: ObsQC
    !PROCEDURE(allBC),  DEFERRED :: ObsBC

    PROCEDURE, PUBLIC :: GetForwardValue  ! Used by superObs to calculate the increment of Obs
    PROCEDURE, PUBLIC :: ObsPrepareForSg  ! Used for observations which has to be prepared before thinning, like radar and satellite
    PROCEDURE, NOPASS, PUBLIC :: Field2ObsIsExisted

    ! A function saving gridded (thinned) obs to raw obs:
    PROCEDURE, PUBLIC :: GriddedObs2ObsBase
    PROCEDURE, PUBLIC :: GetBcgAtObs
  END TYPE ObsBase_t

  ABSTRACT INTERFACE
    SUBROUTINE oneAug(this, configFile, idxFile)
      IMPORT :: ObsBase_t
      CLASS(ObsBase_t) :: this
      CHARACTER(LEN=1024), INTENT(IN) :: configFile
      INTEGER, OPTIONAL :: idxFile
    END SUBROUTINE oneAug

    SUBROUTINE ingest(this, X)
      IMPORT :: ObsBase_t, State_t
      CLASS(ObsBase_t) :: this
      TYPE(State_t) :: X
    END SUBROUTINE ingest
    FUNCTION forward(this, X, O) RESULT(Y)
      IMPORT :: ObsBase_t, State_t, ObsSet_t
      CLASS(ObsBase_t) :: this
      TYPE(State_t) :: X
      TYPE(ObsSet_t), INTENT(IN) :: O ! Telling the operator to map X on to this obs structure
      TYPE(ObsSet_t) :: Y
    END FUNCTION forward
    FUNCTION tangent(this, dX, X) RESULT(Y)
      IMPORT :: ObsBase_t, State_t, ObsSet_t
      CLASS(ObsBase_t) :: this
      TYPE(State_t) :: dX, X
      TYPE(ObsSet_t) :: Y
    END FUNCTION tangent
    FUNCTION adjoint(this, dY, X) RESULT(dX)
      IMPORT :: ObsBase_t, State_t, ObsSet_t
      CLASS(ObsBase_t) :: this
      TYPE(State_t) :: dX, X
      TYPE(ObsSet_t) :: dY
    END FUNCTION adjoint

    SUBROUTINE noAugs(this)
      IMPORT :: ObsBase_t
      CLASS(ObsBase_t) :: this
    END SUBROUTINE noAugs

    SUBROUTINE allBC(this)
      IMPORT :: ObsBase_t
      CLASS(ObsBase_t) :: this
    END SUBROUTINE allBC

  END INTERFACE

CONTAINS

  SUBROUTINE ObsPrepareForSg(this, X)
    CLASS(ObsBase_t) :: this
    TYPE(State_t) :: X

    ! This is an optional holder for data prepare polymorphic subroutine
  END SUBROUTINE

  FUNCTION GetForwardValue(this, X, varnametmp, vIdx, hIdx, tIdx, iv) RESULT(val)
    CLASS(ObsBase_t) :: this
    TYPE(State_t) :: X
    INTEGER(i_kind), INTENT(IN) :: vIdx, hIdx, tIdx
    INTEGER(i_kind), INTENT(IN), OPTIONAL :: iv
    REAL(r_kind) :: val
    CHARACTER(LEN=*) :: varnametmp
    REAL(r_kind) :: k1, k2, k3
    REAL(r_kind) :: pres, p_dry, p_vap, temp, shum, R_dry, R_vap, epsilon_water

    ! val = 0.0D0
    ! Only for those obs match the model state; otherwise returns 0.0D0
    ! IF (X%getVarIdx(TRIM(varnametmp)) .GT. 0) &

    !! TODO **HYJ**: 可能需要换成调用观测算子
    IF (TRIM(varnametmp) .EQ. 'refractivity') THEN
      pres = X%Fields(X%getVarIdx('pres'))%DATA(vIdx, hIdx, tIdx)
      temp = X%Fields(X%getVarIdx('temp'))%DATA(vIdx, hIdx, tIdx)
      shum = X%Fields(X%getVarIdx('qvapor'))%DATA(vIdx, hIdx, tIdx)
      k1 = 77.60E-2
      k2 = 3.73E3
      k3 = 77.60E-2
      R_dry = 287.0597
      R_vap = 461.5250
      epsilon_water = R_dry / R_vap

      ! dry pres and partial pres of vapor
      p_vap = pres * shum / (epsilon_water + (1.0 - epsilon_water) * shum)
      p_dry = pres - p_vap

      val = k1 * (p_dry / temp) + k2 * (p_vap / (temp * temp)) + k3 * (p_vap / temp)

      IF (val .GT. 500) THEN
        PRINT *, 'WARNNING: GNSSRO GetForwardValue greater than 500, ', p_dry, p_vap, temp, shum, val
      END IF
    ELSE
      val = X%Fields(X%getVarIdx(TRIM(varnametmp)))%DATA(vIdx, hIdx, tIdx)
    END IF

  END FUNCTION

  FUNCTION Field2ObsIsExisted(X, varnametmp)
    TYPE(State_t) :: X
    CHARACTER(LEN=*) :: varnametmp
    LOGICAL :: Field2ObsIsExisted

    IF (X%getVarIdx(TRIM(varnametmp)) .EQ. 0) THEN
      Field2ObsIsExisted = .FALSE.
    ELSE
      Field2ObsIsExisted = .TRUE.
      RETURN
    END IF

  END FUNCTION

  ! @brief Yuanfu Xie added this function to calculate variable index in obs:
  FUNCTION getVarIdx(this, varName) RESULT(idx)
    CLASS(ObsBase_t) :: this
    CHARACTER(*) :: varName
    INTEGER(i_kind) :: idx
    INTEGER(i_kind) :: i

    idx = 0
    DO i = LBOUND(this%varNames, 1), UBOUND(this%varNames, 1)
      IF (TRIM(this%varNames(i)) .EQ. TRIM(varName)) THEN
        idx = i
        RETURN
      END IF
    END DO
  END FUNCTION

  !> @brief
  !==================================================================
  !  This routine saves a gridded obs type to raw obs type. It aims
  !  at using thinned obs at the finest level as raw obs in analysis.
  !  For observation datasets with much more obs than analysis grids,
  !  it saves the computing time.
  !
  !  Author: Yuanfu Xie 2022-05-11
  SUBROUTINE GriddedObs2ObsBase(this, state, thinObs)
    CLASS(ObsBase_t) :: this
    TYPE(State_t), INTENT(IN)  :: state
    TYPE(ObsSet_t), INTENT(IN) :: thinObs

    ! Local variables:
    INTEGER(i_kind) :: iv, io, im, nv

    ! Allocate raw obs variable name array:
    ALLOCATE (this%varNames(UBOUND(thinObs%ObsFields, 1)))

    im = 0
    nv = UBOUND(thinObs%ObsFields, 1)
    DO iv = 1, nv
      IF (UBOUND(thinObs%ObsFields(iv)%values, 1) .GT. im) im = iv
    END DO
    this%numobs = UBOUND(thinObs%ObsFields(im)%values, 1)

    ! For all variables:
    DO iv = 1, nv
      this%varNames(iv) = thinObs%ObsFields(iv)%Get_Name()
    END DO

    DO io = 1, this%numObs
      this%olatlon(:, io) = state%sg%cell_cntr(:, thinObs%ObsFields(im)%idx(io)%hIdx)
      this%obsHght(io) = state%sg%zHght(thinObs%ObsFields(im)%idx(io)%vIdx, io)
      this%obsTime(io) = state%sg%tt(thinObs%ObsFields(im)%idx(io)%tIdx)
    END DO
  END SUBROUTINE GriddedObs2ObsBase

  !> @brief
  SUBROUTINE GetBcgAtObs(this, state, rawObs, bcgObs)
    IMPLICIT NONE

    CLASS(ObsBase_t) :: this
    TYPE(State_t), INTENT(IN)  :: state
    TYPE(ObsRaw_t), INTENT(OUT) :: rawObs, bcgObs

    ! Local variables:
    CHARACTER(LEN=20) :: varname
    INTEGER(i_kind) :: numObs(this%numVars)
    INTEGER(i_kind) :: nhstencil
    TYPE(domainCheck_t) :: domain

    ! Loop variables:
    INTEGER(i_kind) :: io, iu, iv, i, j, k, istatus

    ! Consider to save these variables in ObsBase:
    INTEGER(i_kind) :: numValided, maskValided(this%numobs)
    INTEGER(i_kind), ALLOCATABLE :: &
      idxGrdValided(:, :), idxHgtValided(:, :, :), idxTimValided(:, :)
    REAL(r_kind), ALLOCATABLE :: &
      coeGrdValided(:, :), coeHgtValided(:, :, :), coeTimValided(:, :), &
      forward(:, :, :) ! vertical, horizontal and time

!#define TRACE_OBS
    ! Skip if the current processor is not active:
    IF (.NOT. state%sg%isActiveProc()) RETURN

    ! Output obs data information:
    DO iu = 1, this%numVars
      ! Check field to Obs exist:
      varname = TRIM(this%varNames(iu))
      PRINT *, varname, this%numObs
      IF (.NOT. this%Field2ObsIsExisted(state, varname)) CYCLE
    END DO

    ! Validating the obs and saving their interpolation coefficients:
    IF (TRIM(this%obsType) == 'SOUND') THEN
      CALL domain%validation(state%sg, this%numobs, this%olatlon, this%obsHght, &
                             this%obsTime, numValided, maskValided, &
                             idxGrdValided, idxHgtValided, idxTimValided, &
                             coeGrdValided, coeHgtValided, coeTimValided, &
                             this%interpolation, nhstencil, "SOUND_d")

    ELSE
      CALL domain%validation(state%sg, this%numobs, this%olatlon, this%obsHght, &
                             this%obsTime, numValided, maskValided, &
                             idxGrdValided, idxHgtValided, idxTimValided, &
                             coeGrdValided, coeHgtValided, coeTimValided, &
                             this%interpolation, nhstencil, this%obsType)
    END IF

    ! Allocate valided parameter arrays:
    ALLOCATE (rawObs%obsData(numValided, this%numVars), &
              rawObs%obsErrs(numValided, this%numVars), &
              rawObs%olatlon(2, numValided), &
              rawObs%height(numValided), &
              rawObs%tempor(numValided), &
              rawObs%varNames(this%numVars))

    ALLOCATE (forward(2, nhstencil, 2))  ! vertical, horizontal, temporal

    ! Pass the original arrays to the valided arrays:
    i = 0
    DO io = 1, this%numObs
      IF (maskValided(io) .GT. 0) THEN
        i = i + 1
        rawObs%obsData(i, :) = this%obsData(io, 1:this%numVars)
        rawObs%obsErrs(i, :) = this%obsErrs(io, 1:this%numVars)
        rawObs%olatlon(:, i) = this%olatlon(:, io)
        rawObs%height(i) = this%obsHght(io)
        rawObs%tempor(i) = this%obsTime(io)

      END IF
    END DO
    rawObs%varNames = this%varNames
    rawObs%obsType = this%obsType

    bcgObs = rawObs
    bcgObs%obsType = this%obsType

    ! Thinning valid data only:
    DO io = 1, numValided
      DO iv = 1, this%numVars
        varname = TRIM(this%varNames(iv))
        IF (.NOT. this%Field2ObsIsExisted(state, varname)) CYCLE

        ! rawObs%varNames(iv) = TRIM(this%varNames(iv))//TRIM(this%obsType)
        ! bcgObs%varNames(iv) = TRIM(this%varNames(iv))//TRIM(this%obsType)

        ! Background value at obs location if required:
        bcgObs%obsData(io, iv) = 0.0D0
        DO i = 1, 2 ! 2 time frames of interpolation
          DO j = 1, nhstencil ! 3 horizontal interpolation points
            DO k = 1, 2 ! 2 vertical levels
              forward(k, j, i) = this%GetForwardValue(state, varname, &
                                                      idxHgtValided(k, j, io), idxGrdValided(j, io), idxTimValided(i, io), iv)

              ! Interpolate the background at observation:
              bcgObs%obsData(io, iv) = bcgObs%obsData(io, iv) + forward(k, j, i) * &
                                       coeHgtValided(k, j, io) * coeGrdValided(j, io) * coeTimValided(i, io)
            END DO
          END DO
        END DO

      END DO
    END DO

  END SUBROUTINE GetBcgAtObs

  !> @brief
  !==================================================================
  !  This is an improved superobing with unified 4D interpolation
  !   Input:
  !     state:    state_t with geometry and background fields
  !     mpObs:    mpObs_t Obs multiprocessers
  !   Output:
  !     thinObs:  ObsSet_t thinned or super obs
  !
  !   Modified by Yuanfu Xie changing the thinning weighting. 2022-06-03
  !   Modified by Yuanfu Xie, adding mask check for thinned obs, 2024-07-2
  !
  SUBROUTINE ObsThinning(this, state, thinObs, mpObs, useBkg, vector)
    IMPLICIT NONE

    CLASS(ObsBase_t) :: this
    TYPE(State_t), INTENT(IN)  :: state
    TYPE(ObsSet_t), INTENT(INOUT) :: thinObs
    TYPE(mpObs_t), INTENT(IN)  :: mpObs
    LOGICAL, INTENT(IN) :: useBkg, vector

    ! Local variables:
    CHARACTER*20, PARAMETER :: header = 'ObsBase>ObsThinning: '
    INTEGER(i_kind), PARAMETER :: header_len = 22

    ! Thinning parameters:
    REAL(r_kind), ALLOCATABLE :: weights(:, :, :, :)       ! Observation weight
    REAL(r_kind), ALLOCATABLE :: wghtobs(:, :, :, :)       ! Weighted obs: sum w*obs
    REAL(r_kind), ALLOCATABLE :: wghterr(:, :, :, :)       ! Weighted err: sum w*err

    ! Local variables:
    CHARACTER(LEN=20) :: varname
    INTEGER(i_kind) :: numObs(this%numVars)
    INTEGER(i_kind) :: nhstencil
    TYPE(domainCheck_t) :: domain

    REAL(r_kind), PARAMETER :: sigma2 = 4.0D0 !8.0

    ! Loop variables:
    INTEGER(i_kind) :: io, iu, iv, i, j, k, istatus
    REAL(r_kind) :: bkgdAtObs, gaussian, landmaskAtObs, topographyAtObs, t1, t2, t11, t22, t111, t222, bkgd

    ! Consider to save these variables in ObsBase:
    INTEGER(i_kind) :: numValided, maskValided(this%numobs), numQC(this%numVars), numMiss(this%numVars)
    INTEGER(i_kind), ALLOCATABLE :: &
      idxGrdValided(:, :), idxHgtValided(:, :, :), idxTimValided(:, :)
    REAL(r_kind), ALLOCATABLE :: &
      coeGrdValided(:, :), coeHgtValided(:, :, :), coeTimValided(:, :), &
      forward(:, :, :) ! vertical, horizontal and time

    ! All parameters of valided obs:
    CHARACTER(LEN=1024), ALLOCATABLE :: valided_stName(:)
    CHARACTER(LEN=5) OneStationName
    REAL(r_kind), ALLOCATABLE :: valided_obsData(:, :), valided_obsErrs(:, :), &
                                 valided_olatlon(:, :), valided_height(:), valided_tempor(:)
    LOGICAL :: thresQC

!#define TRACE_OBS

#ifdef TRACE_OBS
    INTEGER(i_kind), PARAMETER :: io_debug = 104, iovalid_debug = 639, myrank_debug = 2, glevel_debug = 5, &
                                  icell_debug = 72, ihght_debug = 43, iv_debug = 1 ! iv changes with obsType!!!
#endif

    ! Skip if the current processor is not active:
    IF (.NOT. state%sg%isActiveProc()) RETURN

! #ifdef TRACE_OBS
!     IF (TRIM(this%obsType) .EQ. 'SYNOP') THEN
!       ! debugging 2022-08-26: for synop data only for now
!       ! Check if an obs with given station name is in the original obs set and its LL/H/T info:
!       !+++++++++++++++++++++++++++
!       ! This is the THIRD step to find an obs with given station found in the SECOND step:
!       !+++++++++++++++++++++++++++
!       OneStationName = 'G3502' !'56964' !'57494' !'98233' !'98328' !'56985'
!       DO io = 1, this%numObs
!         IF (this%stNames(io) (1:5) .EQ. OneStationName .AND. state%sg%glevel .EQ. glevel_debug .AND. &
!             state%sg%mpddInfo_sg%myrank .EQ. myrank_debug) & ! .AND. io .EQ. io_debug) &
!           WRITE (*, 702) this%stNames(io) (1:5), io, this%olatlon(:, io) / degree2radian, &
!           this%obsData(io, 2), this%obsHght(io), this%obsTime(io) - state%sg%tt(1), &
!           state%sg%tt(state%sg%tSlots) - state%sg%tt(1), state%sg%gLevel
! 702     FORMAT('OBS at a station: ', A5, I5, ' LL:', 2D12.5, ' q:', D12.4, ' H/T: ', 2D13.5, &
!                ' WDW: ', D12.5, ' pc:', I2, ' glvl: ', I2)

!         IF (this%stNames(io) (1:5) .EQ. 'G3502' ) THEN !& .OR. this%stNames(io) (1:5) .EQ. '59486' .OR. &
!           ! this%stNames(io) (1:5) .EQ. 'G3672') &
!           WRITE (*, 777) io, &
!             this%stNames(io) (1:5), this%olatlon(:, io) / degree2radian, &
!             this%obsData(io, 6), this%obsTime(io) - state%sg%tt(1), state%sg%gLevel
! 777       FORMAT('ThreeStation !: ', I6, ' Name', A, ' LL:', 2E16.8, ' Obs', E12.4, ' Tim: ', E12.4, ' G', I2)
!         END IF
!       END DO
!     END IF
! #endif

    ! Output obs data information:
    DO iu = 1, this%numVars
      ! Check field to Obs exist:
      varname = TRIM(this%varNames(iu))
      PRINT *, 'check check ', varname, this%Field2ObsIsExisted(state, varname)
      WRITE (*, 1111) TRIM(this%varNames(iu)), TRIM(this%obsType), state%sg%mpddInfo_sg%myrank
1111  FORMAT('Thinning check var name: ', A, ' TYPE: ', A, ' pc:', I2)
      IF (.NOT. this%Field2ObsIsExisted(state, varname)) CYCLE

      WRITE (*, 1) TRIM(varname), MAXVAL(this%obsData(:, iu), this%obsData(:, iu) .LT. 1.0D8), &
        MINVAL(this%obsData(:, iu), this%obsData(:, iu) .GT. -1.0D8), &
        state%sg%mpddInfo_sg%myrank, TRIM(this%obsType), state%sg%gLevel
1     FORMAT('MAX/MIN obs: ', A, ' before thinning: ', 2D12.4, ' at proc: ', I1, ' Type: ', A, ' G: ', I2)
    END DO

    ! Initialize thinObs
    thinObs = ObsSet_t(this%configFile, mpObs)

    ! Allocate memory for thinning weights:
    ALLOCATE ( &
      weights(state%sg%vLevel, state%sg%num_cell, state%sg%tSlots, this%numVars), &
      wghtobs(state%sg%vLevel, state%sg%num_cell, state%sg%tSlots, this%numVars), &
      wghterr(state%sg%vLevel, state%sg%num_cell, state%sg%tSlots, this%numVars), &
      STAT=istatus)
    IF (istatus .NE. 0) THEN
      PRINT *, header, 'Cannot allocate memory for weights'
      STOP
    END IF

    ! Initializing:
    weights = 0.0D0
    wghtobs = 0.0D0
    wghterr = 0.0D0

    ! Validating the obs and saving their interpolation coefficients:
    IF (TRIM(this%obsType) == 'SOUND') THEN
      CALL domain%validation(state%sg, this%numobs, this%olatlon, this%obsHght, &
                             this%obsTime, numValided, maskValided, &
                             idxGrdValided, idxHgtValided, idxTimValided, &
                             coeGrdValided, coeHgtValided, coeTimValided, &
                             this%interpolation, nhstencil, "SOUND_d")

    ELSE
      CALL domain%validation(state%sg, this%numobs, this%olatlon, this%obsHght, &
                             this%obsTime, numValided, maskValided, &
                             idxGrdValided, idxHgtValided, idxTimValided, &
                             coeGrdValided, coeHgtValided, coeTimValided, &
                             this%interpolation, nhstencil, this%obsType)
    END IF

    WRITE (*, 21) numValided, this%numobs, state%mpddGlob%myrank, state%sg%gLevel, TRIM(this%obsType)
21  FORMAT('Valided data in thinning: ', I10, ' out of raw: ', I10, ' at proc: ', I2, ' at Glevel: ', I2, ' DataType: ', A)

    ! Check if there is any valid data:
    ! IF (numValided .EQ. 0) THEN
    !   ALLOCATE(thinObs%ObsFields(0))
    !   RETURN
    ! ELSE
    ALLOCATE (thinObs%ObsFields(this%numVars))
    ! END IF

    ! Allocate valided parameter arrays:
    ALLOCATE (valided_obsData(numValided, this%numVars), &
              valided_obsErrs(numValided, this%numVars), &
              valided_olatlon(2, numValided), &
              valided_height(numValided), &
              valided_tempor(numValided))

    !IF (TRIM(this%obsType) .EQ. 'SYNOP') &
    ALLOCATE (valided_stName(numValided))

    ALLOCATE (forward(2, nhstencil, 2))  ! vertical, horizontal, temporal

    CALL CPU_TIME(t1)

    ! Pass the original arrays to the valided arrays:
    i = 0
    DO io = 1, this%numObs
      IF (maskValided(io) .GT. 0) THEN
        i = i + 1
        valided_obsData(i, :) = this%obsData(io, 1:this%numVars)
        valided_obsErrs(i, :) = this%obsErrs(io, 1:this%numVars)
        valided_olatlon(:, i) = this%olatlon(:, io)
        valided_height(i) = this%obsHght(io)
        valided_tempor(i) = this%obsTime(io)

        ! Track the station names for debugging:
        IF (TRIM(this%obsType) .EQ. 'SYNOP' .OR. TRIM(this%obsType) .EQ. 'SOUND') &
          valided_stName(i) = this%stNames(io)

      END IF

! #ifdef TRACE_OBS
!       ! Debugging:
!       IF (this%stNames(io) (1:5) .EQ. OneStationName .AND. &
!           io .GE. io_debug .AND. io .LE. io_debug .AND. &
!           state%sg%glevel .EQ. glevel_debug .AND. &
!           state%sg%mpddInfo_sg%myrank .EQ. myrank_debug) &
!         WRITE (*, 123) i, maskValided(io), this%obsData(io, 2), &
!         this%obsHght(io), state%sg%mpddGlob%myrank, state%sg%gLevel
! 123   FORMAT('Check the mask of station: ', 2I4, ' Q/H: ', 2D12.4, ' pc:', I2, ' glvl: ', I2)

!       IF (this%stNames(io) (1:5) .EQ. 'G3502') THEN
!         WRITE (*, 701) io, this%obsData(io, 6), &
!           this%obsTime(io) - state%sg%tt(1), &
!           state%sg%tt(state%sg%tSlots) - state%sg%tt(1), state%sg%tSlots, maskValided(io), &
!           state%sg%gLevel, state%sg%mpddInfo_sg%myrank
! 701     FORMAT('G3502 Valid check: ', I5, ' obs:', E12.4, ' obT-T1:', E12.4, &
!                ' T2-T1:', E12.4, ' Slot:', I3, ' Mask: ', I6, ' G', I2, ' pc:', I2)
!       END IF
! #endif

    END DO

    ! Correct the pressure and temperature of the surface observations
    IF (TRIM(this%obsType) == 'SYNOP' .AND. state%sg%vLevel > 1) THEN
      BLOCK
        REAL(r_kind) :: diffTopo
        REAL(r_kind), ALLOCATABLE :: topo(:)

        ALLOCATE (topo(numValided))
        topo = 0.0D0

        DO i = 1, numValided
          DO j = 1, nhstencil ! 3 horizontal interpolation points if slint, 4 for regular
            topo(i) = topo(i) + state%sg%topo(idxGrdValided(j, i)) * coeGrdValided(j, i)
          END DO
        END DO

        IF (this%tempIdx /= 0) THEN
          diffTopo = 0.0D0
          DO i = 1, numValided
            IF (ABS(topo(i) - valided_height(i)) < 2000.0D0) THEN
              valided_obsData(i, this%tempIdx) = valided_obsData(i, this%tempIdx) - 0.65 * (topo(i) - valided_height(i)) / 100.0D0
            ELSE
              IF (diffTopo .LT. ABS(topo(i) - valided_height(i))) &
                diffTopo = ABS(topo(i) - valided_height(i))
              valided_obsData(i, this%tempIdx) = 1E9
            END IF
          END DO
          WRITE (*, 51) diffTopo
51        FORMAT('Synop - Max Topo and obs height diff: ', E14.6)
        END IF

        IF (this%presIdx /= 0 .AND. this%tempIdx /= 0) THEN
          DO i = 1, numValided
            IF (ABS(topo(i) - valided_height(i)) < 2000.0D0) THEN
              bkgd = valided_obsData(i, this%presIdx)
              valided_obsData(i, this%presIdx) = &
                valided_obsData(i, this%presIdx) * EXP(-(topo(i) - valided_height(i)) * 9.80655D0 / (dry_air_gas_const * &
                                                                                                     ((valided_obsData(i, this%tempIdx) - 0.65 * (topo(i) - valided_height(i)) / 200.0D0) &
                                                                                                      * (1.0D0 + 0.608 * valided_obsData(i, this%qvaporIdx)))))
              ! write(*,77) i,numValided,valided_obsData(i, this%presIdx)-bkgd, &
              !  valided_obsData(i, this%presIdx),bkgd, &
              !  valided_obsData(i, this%tempIdx), &
              !  valided_height(i), topo(i),state%sg%gLevel
! 77            format('Hydro P at obs: ', 2I6, ' q/p/pold/t:', 4D14.6, ' h/topo', 2D12.4, ' G:', I2)
            ELSE
              PRINT *, 'The diff is ', topo(i) - valided_height(i)
              valided_obsData(i, this%tempIdx) = 1E9
              valided_obsData(i, this%presIdx) = 1E9
            END IF
          END DO
        END IF

        IF (this%lnpIdx /= 0 .AND. this%tempIdx /= 0) THEN
          DO i = 1, numValided
            IF (ABS(topo(i) - valided_height(i)) < 2000.0D0) THEN
              bkgd = valided_obsData(i, this%lnpIdx)
              valided_obsData(i, this%lnpIdx) = &
                valided_obsData(i, this%lnpIdx) + (-(topo(i) - valided_height(i)) * 9.80655D0 / (dry_air_gas_const * &
                                                                                                 ((valided_obsData(i, this%tempIdx) - 0.65 * (topo(i) - valided_height(i)) / 200.0D0) &
                                                                                                  * (1.0D0 + 0.608 * valided_obsData(i, this%qvaporIdx)))))
            ELSE
              PRINT *, 'The diff is ', topo(i) - valided_height(i)
              valided_obsData(i, this%tempIdx) = 1E9
              valided_obsData(i, this%lnpIdx) = 1E9
            END IF
          END DO
        END IF

        DEALLOCATE (topo)
      END BLOCK
    END IF

    CALL CPU_TIME(t11)

    ! Thinning valid data only:
    numQC = 0
    numMiss = 0
    DO io = 1, numValided

! #ifdef TRACE_OBS
!       ! debugging 2022-08-26:
!       ! Following the given obs station in the valided obs set:
!       IF (valided_stName(io) (1:5) .EQ. OneStationName .AND. io .EQ. iovalid_debug) &
!         WRITE (*, 903) io, valided_stName(io) (1:5), valided_olatlon(:, io) / degree2radian, &
!         valided_tempor(io) - state%sg%tt(1), valided_obsData(io, 2), valided_height(io), &
!         state%sg%num_icell, state%sg%mpddInfo_sg%myrank, state%sg%gLevel
! 903   FORMAT('valided station - order: ', I5, ' SN: ', A, ' oLL:', 2D12.4, ' oT:', D12.5, &
!              ' Q/H: ', 2D12.5, ' NiC', I7, ' pc:', I2, ' gLvl:', I2)

!       ! Find out a given station horizontal indices:
!       IF (valided_stName(io) (1:5) .EQ. OneStationName) THEN
!         WRITE (*, 910) io, valided_stName(io) (1:5), valided_olatlon(:, io) / degree2radian, &
!           valided_tempor(io) - state%sg%tt(1), idxGrdValided(1, io), &
!           state%sg%num_icell, state%sg%mpddInfo_sg%myrank, state%sg%gLevel, this%varNames(6) (1:5), &
!           this%obsData(io, 6)
! 910     FORMAT('Find a station indices: ', I6, ' Station: ', A, 2E12.4, ' obsT: ', E12.4, &
!                ' Horizontal index: ', 2I8, ' pc: ', I2, ' G:', I2, A, ' PCP: ', D12.4)
!       END IF
! #endif

      DO iv = 1, this%numVars
        varname = TRIM(this%varNames(iv))
        !print*,'VARNAME: ',iv,io,TRIM(varname),this%getVarIdx(varname),state%mpddGlob%myrank
        IF (.NOT. this%Field2ObsIsExisted(state, varname) .OR. &
            ABS(valided_obsData(io, iv)) .GE. ABS(missing)) THEN
          numMiss(iv) = numMiss(iv) + 1
          CYCLE
        END IF

        ! Background value at obs location if required:
        bkgdAtObs = 0.0D0
        IF (useBkg) THEN
          DO i = 1, 2 ! 2 time frames of interpolation
            DO j = 1, nhstencil ! 3 horizontal interpolation points
              DO k = 1, 2 ! 2 vertical levels
                forward(k, j, i) = this%GetForwardValue(state, varname, &
                                                        idxHgtValided(k, j, io), idxGrdValided(j, io), idxTimValided(i, io), iv)
                ! Interpolate the background at observation:
                bkgdAtObs = bkgdAtObs + forward(k, j, i) * &
                            coeHgtValided(k, j, io) * coeGrdValided(j, io) * coeTimValided(i, io)
              END DO
            END DO
          END DO
        END IF

        ! Landmask and topography at obs location:
        landmaskAtObs = 0.0D0
        topographyAtObs = 0.0D0
        DO j = 1, nhstencil ! 3 horizontal interpolation points
          landmaskAtObs = landmaskAtObs + state%sg%landmask(idxGrdValided(j, io)) * coeGrdValided(j, io)
          topographyAtObs = topographyAtObs + state%sg%topo(idxGrdValided(j, io)) * coeGrdValided(j, io)
        END DO

! #ifdef TRACE_OBS
!         ! Debugging: Check if the obs passes a QC check:
!         IF (iv .EQ. iv_debug .AND. io .EQ. iovalid_debug .AND. &
!             state%sg%glevel .EQ. glevel_debug .AND. &
!             state%sg%mpddInfo_sg%myrank .EQ. myrank_debug) &
!           WRITE (*, 888) bkgdAtObs - valided_obsData(io, iv), this%qcThreshold(iv), &
!           bkgdAtObs, valided_obsData(io, iv), valided_height(io), valided_stName(io) (1:5), TRIM(this%obsType), &
!           TRIM(this%varNames(iv))
! 888     FORMAT('QCcheck: ', D12.4, ' threshold: ', D12.4, ' bkg-obs:', 2D12.4, &
!                ' oHght: ', D12.4, ' valiedStation: ', A, ' obsType: ', A, ' varName: ', A)

!         IF (valided_stName(io) (1:5) .EQ. 'G3502' .AND. iv .EQ. 6) THEN
!           WRITE (*, 101) io, bkgdAtObs, valided_obsData(io, iv), &
!             valided_tempor(io) - state%sg%tt(1), &
!             state%sg%tt(state%sg%tSlots) - state%sg%tt(1), &
!             state%sg%gLevel, state%sg%mpddInfo_sg%myrank
! 101       FORMAT('G3502 check: ', I5, ' bkg:', E12.4, ' obs:', E12.4, ' obT-T1:', E12.4, &
!                  ' T2-T1:', E12.4, ' G', I2, ' pc:', I2)
!           WRITE (*, 201) io, state%sg%gLevel, bkgdAtObs, valided_obsData(io, iv), &
!             idxHgtValided(1:2, j, io), idxTimValided(1:2, io), idxGrdValided(:, io), this%qcThreshold(iv)
! 201       FORMAT('G3502--:', I7, ' G', I1, ' bkg', E12.4, ' obs', E12.4, ' Hgh:', 2I2, ' Tim', 2I2, ' Hori: ', 3I7, ' thresh: ', E12.4)
!         END IF

!         ! Check if any other obs affecting the grid:
!         IF (iv .EQ. 6 .AND. state%sg%gLevel .EQ. 8 .AND. idxTimValided(2, io) .EQ. 9 .AND. &
!             (ANY(idxGrdValided(:, io) .EQ. 207788) .AND. ANY(idxGrdValided(:, io) .EQ. 207789) .AND. &
!              ANY(idxGrdValided(:, io) .EQ. 208365))) THEN
!           WRITE (*, 301) io, valided_stName(io) (1:5), idxGrdValided(:, io), &
!             valided_olatlon(:, io) / degree2radian
! 301       FORMAT('Who uses these grid: ', I6, ' StName: ', A, ' idx: ', 3I7, ' LL:', 2E12.4)
!         END IF
! #endif

        ! Check if the obs passes a QC check:
        thresQC = .FALSE.
        IF (this%lnpIdx == iv) THEN
          IF (ABS(EXP(bkgdAtObs) - EXP(valided_obsData(io, iv))) .LE. this%qcThreshold(iv)) thresQC = .TRUE.
        ELSE
          IF (ABS(bkgdAtObs - valided_obsData(io, iv)) .LE. this%qcThreshold(iv)) thresQC = .TRUE.
        END IF

        IF (thresQC) THEN
          numQC(iv) = numQC(iv) + 1
          ! Calculate the Gaussian weights:
          DO i = 1, 2 ! 2 time frames of interpolation
            ! Yuanfu Xie added to prevent weighting on grid with zero coefficient 2023-07-15
            IF (coeTimValided(i, io) .LE. machineEps) CYCLE
            DO j = 1, nhstencil ! 3 horizontal interpolation points if slint, 4 for regular
              ! Yuanfu Xie added to prevent weighting on grid with zero coefficient 2023-07-15
              IF (coeGrdValided(j, io) .LE. machineEps) CYCLE
              DO k = 1, 2 ! 2 vertical levels
                ! Yuanfu Xie added to prevent weighting on grid with zero coefficient 2023-07-15
                IF (coeHgtValided(k, j, io) .LE. machineEps) CYCLE
                gaussian = EXP(-(1.0D0 - coeHgtValided(k, j, io))**2 / sigma2 &
                               - (1.0D0 - coeGrdValided(j, io))**2 / sigma2 &
                               - (1.0D0 - coeTimValided(i, io))**2 / sigma2)

                ! Below is the new weighting scheme pending testing, committed by Zilong 2024-04-28 in merge -----
                ! Consider the background flow dependency: Yuanfu Xie reactivates this on 2023-10-12:
                ! IF (state%sg%flowDependentThinning) gaussian = gaussian * &
                !                                                EXP(-4.0D0 * (bkgdAtObs - forward(k, j, i))**2)

                ! Consider the landmask:
                ! gaussian = gaussian * ((1.0D0 - state%sg%vertcalSimilarity(idxHgtValided(k, j, io))) * &
                !                        EXP(-state%sg%scales4Similarity(1) * (landmaskAtObs - state%sg%landmask(idxGrdValided(j, io)))**2) + &
                !                        EXP(-state%sg%scales4Similarity(2) * (topographyAtObs - state%sg%topo(idxGrdValided(j, io)))**2) + &
                !                        state%sg%vertcalSimilarity(idxHgtValided(k, j, io)))
                ! ------------------------------------------------------------------------------

                ! If Gaussian shows good correlation:
                ! Yuanfu Xie adds mask check on 2024-07-02: mask >= 1 leaves the option for using 2 for dependent variables in a 4DVAR
                IF (gaussian .GE. this%correlation_threshold &
                    ! .AND. state%fields(state%getVarIdx(TRIM(varname)))%maskVertical(idxHgtValided(k, j, io)) .GE. 1 .AND. &
                    ! state%fields(state%getVarIdx(TRIM(varname)))%maskHorizontal .GE. 1 .AND. &
                    ! state%fields(state%getVarIdx(TRIM(varname)))%maskTemporal(idxTimValided(i, io)) .GE. 1 &
                    ) THEN

                  weights(idxHgtValided(k, j, io), idxGrdValided(j, io), idxTimValided(i, io), iv) = &
                    weights(idxHgtValided(k, j, io), idxGrdValided(j, io), idxTimValided(i, io), iv) + gaussian !1.0 ! Yuanfu Xie changed the weight 2022-06-3

                  ! Sum of weighted innovation: W*(obs-bkg)
                  wghtobs(idxHgtValided(k, j, io), idxGrdValided(j, io), idxTimValided(i, io), iv) = &
                    wghtobs(idxHgtValided(k, j, io), idxGrdValided(j, io), idxTimValided(i, io), iv) + &
                    gaussian * (valided_obsData(io, iv) - bkgdAtObs)

                  wghterr(idxHgtValided(k, j, io), idxGrdValided(j, io), idxTimValided(i, io), iv) = &
                    wghterr(idxHgtValided(k, j, io), idxGrdValided(j, io), idxTimValided(i, io), iv) + &
                    gaussian * valided_obsErrs(io, iv)

! #ifdef TRACE_OBS
!                   ! IF (valided_stName(io) (1:5) .EQ. 'G3502' .AND. iv .EQ. 6 .AND. state%sg%gLevel .EQ. 8 .AND. &
!                   IF (idxGrdValided(j, io) .EQ. 207788 .OR. idxGrdValided(j, io) .EQ. 207789 .OR. &
!                       idxGrdValided(j, io) .EQ. 208365) THEN
!                     IF (iv .EQ. 6 .AND. state%sg%gLevel .EQ. 8 .AND. &
!                         valided_tempor(io) - state%sg%tt(1) .GE. 2400.0D0) THEN
!                       WRITE (*, 444) io, idxHgtValided(k, j, io), idxGrdValided(j, io), idxTimValided(i, io), &
!                         wghtobs(idxHgtValided(k, j, io), idxGrdValided(j, io), idxTimValided(i, io), iv), &
!                         gaussian, valided_obsData(io, iv), bkgdAtObs
! 444                   FORMAT('WeightObs G3502:', I7, ' HXT: ', I2, I7, I3, ' wobs: ', E12.4, ' GOB', 3E12.4)
!                       WRITE (*, 666) k, idxHgtValided(k, j, io), i, idxTimValided(i, io), &
!                         (1.0D0 - coeHgtValided(k, j, io))**2, &
!                         (1.0D0 - coeGrdValided(j, io))**2, (1.0D0 - coeTimValided(i, io))**2, &
!                         valided_tempor(io) - state%sg%tt(1), valided_obsData(io, iv)
! 666                   FORMAT('G-P: hidx', 2I3, ' tidx:', 2I3, ' HXT:', 3E12.4, ' oT:', E12.4, ' o:', E12.4)
!                     END IF
!                   END IF
!                   IF (TRIM(this%obsType) .EQ. 'SOUND') THEN
!                     ! Debugging 2022-08-26:
!                     ! The following is to check the row and column number: 258 is G8 with single processor;
!                     ! 129 is G8 with 4 processors;
!                     ! This is usually the second step to trace an obs being thinned to a grid cell in the
!                     ! next grid weighting check.
!                     ! Check if horizontal interpolated index match a given grid number along with the weighting info:
!                     !+++++++++++++++++++++++++++
!                     ! This is the SECOND step to find the obs station information etc.:
!                     !+++++++++++++++++++++++++++
!                     IF (iv .EQ. iv_debug .AND. io .EQ. iovalid_debug .AND. &
!                         state%sg%glevel .EQ. glevel_debug .AND. &
!                         state%sg%mpddInfo_sg%myrank .EQ. myrank_debug) &
!                       WRITE (*, 505) idxGrdValided(j, io), idxHgtValided(k, j, io)
! 505                 FORMAT('The data passed QC threshold check: ', I6, ' hght idx: ', I3)
!                     IF (iv .EQ. iv_debug .AND. idxGrdValided(j, io) .EQ. icell_debug .AND. &
!                         idxHgtValided(k, j, io) .EQ. ihght_debug .AND. &
!                         state%sg%mpddInfo_sg%myrank .EQ. myrank_debug .AND. &
!                         state%sg%gLevel .EQ. glevel_debug .AND. &
!                         TRIM(this%obsType) .EQ. 'SOUND') &
!                       WRITE (*, 501) io, valided_stName(io) (1:5), i, j, k, valided_olatlon(:, io) / degree2radian, &
!                       wghtobs(idxHgtValided(k, j, io), idxGrdValided(j, io), idxTimValided(i, io), iv), &
!                       weights(idxHgtValided(k, j, io), idxGrdValided(j, io), idxTimValided(i, io), iv), &
!                       valided_obsData(io, iv), bkgdAtObs, &
!                       valided_tempor(io) - state%sg%tt(1), &
!                       valided_height(io), &
!                       coeHgtValided(k, j, io), coeGrdValided(j, io), coeTimValided(i, io)
! 501                 FORMAT('Obs to Find: ', I5, ' stNames: ', A5, ' T/H/V: ', 3I2, ' obsLL', 2D12.5, &
!                            ' W: ', 2D12.4, ' obs/bkg: ', 2D12.4 &
!                            ' oT: ', D12.5, ' Hgt: ', D12.4, ' coe:', 3D12.4)
!                   END IF
! #endif

                END IF
              END DO
            END DO
          END DO
        END IF
      END DO
    END DO

    CALL CPU_TIME(t22)
    WRITE (*, 56) t22 - t11, state%sg%mpddInfo_sg%myrank, state%sg%gLevel, TRIM(this%obsType), this%numObs
56  FORMAT('Time spent in Gaussian weights: ', D12.4, ' proc: ', I1, ' Glevel: ', I2, ' TYPE: ', A, ' NOB: ', I10)
    DO iv = 1, this%numVars
      WRITE (*, 57) numQC(iv), numMiss(iv), state%sg%gLevel, TRIM(this%obsType), TRIM(this%varNames(iv))
57    FORMAT('Number of QC obs: ', I6, ' number of missing value: ', I6, ' at G:', I2, ' obsType: ', A, ' of varname: ', A)
    END DO

    ! For each variables: this run counts the numObs
    DO iv = 1, this%numVars
      ! Count the thinned obs:
      numObs(iv) = COUNT(weights(:, 1:state%sg%num_icell, :, iv) .GT. 0.0)
      thinObs%ObsFields(iv) = ObsField_t(this%configFile, mpObs)
      ALLOCATE ( &
        thinObs%ObsFields(iv)%values(numObs(iv)), &
        thinObs%ObsFields(iv)%errors(numObs(iv)), &
        thinObs%ObsFields(iv)%idx(numObs(iv)))

      CALL thinObs%ObsFields(iv)%Set_Name(TRIM(this%varNames(iv)))
      CALL thinObs%ObsFields(iv)%Set_ObsType(TRIM(this%obsType))

      WRITE (*, 58) numObs(iv), state%sg%gLevel, TRIM(this%obsType), TRIM(this%varNames(iv)), &
        thinObs%ObsFields(iv)%Get_Id_Name()
58    FORMAT('Number of Weight>0: ', I6, ' G:', I2, ' Type: ', A, ' varName: ', 2(A, ' '))
    END DO

    CALL CPU_TIME(t111)

    ! Create the obsFields:
    numObs = 0
    DO iu = 1, this%numVars

      ! Check field to Obs exist:
      varname = TRIM(this%varNames(iu))
      IF (.NOT. this%Field2ObsIsExisted(state, varname)) THEN
        ! Changed to halt as all variables read in are matching with the obs: Yuanfu Xie 2022-08-05
        WRITE (*, 41) varname
41      FORMAT('ObsThinning - variable name is not supported: ', A, ' Check and rerun')
        CYCLE !STOP
      END IF

      DO i = 1, state%sg%tSlots
        DO j = 1, state%sg%num_icell
          DO k = 1, state%sg%vLevel

            IF (weights(k, j, i, iu) .GT. 0.0) THEN
              ! Found observation associated with this grid point:
              ! Weighted obs:
              iv = iu
              IF (vector) iv = 1

              numObs(iv) = numObs(iv) + 1
              thinObs%ObsFields(iv)%idx(numObs(iv))%hIdx = j
              thinObs%ObsFields(iv)%idx(numObs(iv))%tIdx = i
              thinObs%ObsFields(iv)%idx(numObs(iv))%vIdx = k

              bkgd = 0.0D0
              IF (useBkg) bkgd = this%GetForwardValue(state, varname, k, j, i, iv)
              IF (vector) THEN
                thinObs%ObsFields(iv)%valueArray(numObs(iv), iu) = &
                  bkgd + &
                  wghtobs(k, j, i, iv) / weights(k, j, i, iv)
              ELSE
                thinObs%ObsFields(iv)%values(numObs(iv)) = &
                  bkgd + &
                  wghtobs(k, j, i, iv) / weights(k, j, i, iv)

! #ifdef TRACE_OBS
!                 IF (iv .EQ. 6 .AND. k .EQ. 1 .AND. &
!                     (j .EQ. 208365 .OR. j .EQ. 207788 .OR. j .EQ. 207789) .AND. i .EQ. 33 .AND. &
!                     state%sg%gLevel .EQ. 8) THEN
!                   WRITE (*, 555) numObs(iv), j, thinObs%ObsFields(iv)%values(numObs(iv)), &
!                     wghtobs(k, j, i, iv), weights(k, j, i, iv), &
!                     MAXVAL(state%sg%edge_leng(:, :))
! 555               FORMAT('Thinned G3502: ', I6, ' J: ', I6, ' obsThin: ', E12.4, ' wob/w: ', 2E12.4, &
!                          ' EdgeLength: ', E14.6)
!                 END IF
! #endif
              END IF

! #ifdef TRACE_OBS
!               IF (.TRUE.) THEN

!                 BLOCK
!                   INTEGER(i_kind) :: ii, jj, iii, jjj, ir

!                   ir = 1

!                   ! Multi processors:
!                   ! Calculate iii jjj for interior point: the grid number is set based on
!                   ! G5 33 x 17 interior points on each processor of 4
!                   iii = MOD(j - 1, 33) + 1
!                   jjj = INT(j / 33.0D0) + 1

!                   ! Check if VWPW has data above vertical level 30:
!                   IF (TRIM(this%obsType) .EQ. 'VWPW' .AND. k .GE. 30 .AND. state%sg%mpddInfo_sg%myrank .EQ. myrank_debug) &
!                     WRITE (*, 432) iii, jjj, k, state%sg%mpddInfo_sg%myrank, iv
! 432               FORMAT('Found VWPW above 30 at rectangle grid: ', 2I3, ' vertical level: ', I2, ' pc: ', I2, ' vid: ', I1)

!                   ! Mark the obs with a special value so that to see it in ncview:
!                   !+++++++++++++++++++++++++++
!                   ! This is the FIRST step to find a point of an obs from ncview:
!                   !+++++++++++++++++++++++++++
!                   IF (TRIM(this%obsType) .EQ. 'VWPW' .AND. iii .EQ. 1 .AND. jjj .EQ. 9 .AND. &
!                       state%sg%gLevel .EQ. glevel_debug .AND. & !k .EQ. 43 .AND. &
!                       iv .EQ. iv_debug .AND. state%sg%mpddInfo_sg%myrank .EQ. myrank_debug) THEN
!                     WRITE (*, 402) iii, jjj, j, i, k, iv, numObs(iv), TRIM(this%obsType), state%sg%glevel, &
!                       state%sg%cell_cntr(:, j) / degree2radian, state%sg%zHght(k, j), state%sg%topo(j), TRIM(thinObs%ObsFields(iv)%Get_Name())
! 402                 FORMAT('Find the cell: ', 2I4, ' J:', I6, ' tm-vertcl:', 2I3, ' var: ', I2, ' i-obs: ', I6, ' type: ', A, &
!                            ' Glvl: ', I2, ' gridLL: ', 2D12.5, ' gridHgt/topo: ', 2D12.4, ' varname: ', A)
!                     ! Set the cell with special value to see it in a ncview plot:
!                     ! thinObs%ObsFields(iv)%values(numObs(iv)) = 500.0D0
!                   END IF
!                 END BLOCK
!               END IF
! #endif

              ! Thinned obs errors:
              thinObs%ObsFields(iv)%errors(numObs(iv)) = &
                1.0 / weights(k, j, i, iv) * MAX(wghterr(k, j, i, iv), 1.0E-6)  ! Yuanfu and Yali corrected this formula to * from /
            END IF
          END DO
        END DO
      END DO
    END DO

    CALL CPU_TIME(t222)
!     WRITE(*,67) t222 - t111,state%sg%mpddInfo_sg%myrank,state%sg%num_cell, &
!     TRIM(this%obsType),this%numObs,numObs(:)
! 67  FORMAT('Time spent on gridded Obs: ',D10.2,' proc: ',I1,' NCELL: ',I6, &
!       ' Type: ',A,' Num thinning: ',I10,' Num thinned: ',10I6)

    ! Check the maximum values of thinned obs:
    IF (vector) THEN
      DO iv = 1, this%numVars
        PRINT *, 'MAXVAL(thinObs%ObsFields(iv)%valueArray): ', &
          state%sg%mpddInfo_sg%myrank, iv, MAXVAL(thinObs%ObsFields(iv)%valueArray(:, iv))
      END DO
    ELSE
      DO iv = 1, this%numVars
        ! Check field to Obs exist:
        varname = TRIM(this%varNames(iv))
        IF (.NOT. this%Field2ObsIsExisted(state, varname)) CYCLE
        IF (numObs(iv) .GT. 0) WRITE (*, 11) state%sg%gLevel, state%sg%mpddInfo_sg%myrank, &
          TRIM(this%obsType), TRIM(varname), iv, MAXVAL(thinObs%ObsFields(iv)%values), numObs(iv)
11      FORMAT('Max thinned obs at Glevel: ', I2, ' proc: ', I2, &
               ' Type: ', A, ' var: ', A, ' iv: ', I1, ' is: ', D20.12, ' NumThinned: ', I6)
      END DO
    END IF

    CALL CPU_TIME(t2)
    WRITE (*, 45) t2 - t1, state%sg%mpddInfo_sg%myrank, this%numObs, TRIM(this%obsType), state%sg%gLevel, state%sg%tSlots
45  FORMAT('Time spent on Thinning: ', D12.4, ' proc: ', I1, ' NOB: ', I10, ' Type: ', A, ' Glevel: ', I2, ' time slots: ', I2)

    ! Deallocate memory:
    DEALLOCATE (idxGrdValided, idxHgtValided, idxTimValided, &
                coeGrdValided, coeHgtValided, coeTimValided)
    DEALLOCATE (weights, wghtobs, wghterr)
    DEALLOCATE (valided_obsData, valided_obsErrs, &
                valided_olatlon, valided_height, valided_tempor)
    IF (TRIM(this%obsType) .EQ. 'SYNOP' .OR. ALLOCATED(valided_stName)) DEALLOCATE (valided_stName)
    DEALLOCATE (forward)

    DO iv = 1, this%numVars
      IF (numObs(iv) .GT. 0) THEN
        WRITE (*, 2) iv, MAXVAL(thinObs%ObsFields(iv)%values(1:numObs(iv))), MINVAL(thinObs%ObsFields(iv)%values(1:numObs(iv))), &
          state%sg%mpddInfo_sg%myrank, TRIM(this%obsType)
2       FORMAT('obsThinning: MAX/MIN thinned: ', I1, 2D12.4, ' proc: ', I2, ' Type: ', A)
      END IF
    END DO

    ! thinObs%ObsFields(this%qvaporIdx)%values( thinObs%ObsFields(this%qvaporIdx)%values<0) = 0
  END SUBROUTINE ObsThinning

  SUBROUTINE ObsDestroy(this)
    CLASS(ObsBase_t) :: this

    IF (ALLOCATED(this%obsData)) DEALLOCATE (this%obsData)
    IF (ALLOCATED(this%types)) DEALLOCATE (this%types)
    IF (ALLOCATED(this%obsErrs)) DEALLOCATE (this%obsErrs)
    IF (ALLOCATED(this%olatlon)) DEALLOCATE (this%olatlon)

    IF (ALLOCATED(this%obsTime)) DEALLOCATE (this%obsTime)
    IF (ALLOCATED(this%obsHght)) DEALLOCATE (this%obsHght)
    IF (ALLOCATED(this%land)) DEALLOCATE (this%land)
    IF (ALLOCATED(this%radius)) DEALLOCATE (this%radius)
    IF (ALLOCATED(this%fileNames)) DEALLOCATE (this%fileNames)
    IF (ALLOCATED(this%StNames)) DEALLOCATE (this%StNames)
    IF (ALLOCATED(this%varNames)) DEALLOCATE (this%varNames)

    IF (ALLOCATED(this%sizeInc)) DEALLOCATE (this%sizeInc)
    IF (ALLOCATED(this%qcThreshold)) DEALLOCATE (this%qcThreshold)

    ! SoloSite if allocated:
    IF (ALLOCATED(this%soloNames)) DEALLOCATE (this%soloNames)
    IF (ALLOCATED(this%soloLat)) DEALLOCATE (this%soloLat)
    IF (ALLOCATED(this%soloLon)) DEALLOCATE (this%soloLon)

  END SUBROUTINE ObsDestroy

END MODULE ObsBase_m
