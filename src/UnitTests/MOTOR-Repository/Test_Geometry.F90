!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.geometry.Test
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : Beta 0.1
! HISTORY           :
!  Created by Zilong Qin (zilong.qin@gmail.com), 2020/11/19, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
PROGRAM test_geometry
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  ! USE mpi_f08
  INCLUDE "mpif.h"

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  REAL(r_kind), ALLOCATABLE :: Xm1(:, :), Xm2(:, :), Xm3(:, :), XmT(:, :)
  INTEGER(i_kind), PARAMETER :: vLevel = 100
  INTEGER(i_kind) :: i, gLevel1, gLevel2, gLevel3
  REAL(r_kind) :: sumRelt1, sumRelt2, sumTemp
  CHARACTER(LEN=1024) :: configFile

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testGM.yaml"

  gLevel1 = 1
  gLevel2 = 2
  gLevel3 = 3

  PRINT *, 'Start the program.'
  ! Initialize the mpdd
  CALL mpddGlob%initialize()

  ! Initialize the geometry
  CALL geometry%initialize(configFile, mpddGlob)

  CALL geometry%mg%sg(gLevel1)%allocateMat(vLevel, Xm1)
  CALL geometry%mg%sg(gLevel2)%allocateMat(vLevel, Xm2)
  CALL geometry%mg%sg(gLevel3)%allocateMat(vLevel, Xm3)

  IF (mpddGlob%isBaseProc()) ALLOCATE (XmT(vLevel, geometry%mg%sg(gLevel2)%num_icell_global))
  IF (mpddGlob%isBaseProc()) FORALL (i=1:geometry%mg%sg(gLevel2)%num_icell_global) XmT(:, i) = i

  IF (mpddGlob%isBaseProc()) PRINT *, XmT(1, :)
  CALL geometry%mg%sg(gLevel2)%distGridReal(XmT, Xm2, [vLevel, geometry%mg%sg(gLevel2)%num_icell_global])
  CALL geometry%mg%sg(gLevel2)%aggrGridReal(Xm2, XmT, [vLevel, geometry%mg%sg(gLevel2)%num_icell_global])
  IF (mpddGlob%isBaseProc()) PRINT *, XmT(1, :)

  ! Test restriction
  DO i = 1, geometry%mg%sg(gLevel3)%num_icell
    Xm3(:, i) = geometry%mg%sg(gLevel3)%sp_t_g_idx(i)
  END DO

  CALL geometry%mg%sg(gLevel3)%ExchangeMatOnHalo2D(vLevel, Xm3)
  CALL geometry%mg%restrictionAtGLevel(gLevel3, vLevel, Xm3, Xm2)

  sumTemp = SUM(Xm2(1, 1:geometry%mg%sg(gLevel2)%num_icell))
  sumRelt1 = 0.0D0

  CALL MPI_ALLREDUCE(sumTemp, sumRelt1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpddGlob%comm, mpddGlob%ierr)

  IF (mpddGlob%isBaseProc()) THEN
    PRINT *, 'Sum is', sumRelt1
  END IF

  ! Test prolongation
  DO i = 1, geometry%mg%sg(gLevel1)%num_icell
    Xm1(:, i) = geometry%mg%sg(gLevel1)%sp_t_g_idx(i)
  END DO

  CALL geometry%mg%sg(gLevel1)%ExchangeMatOnHalo2D(vLevel, Xm1)
  CALL geometry%mg%prolongationAtGLevel(gLevel1, vLevel, Xm2, Xm1)

  sumTemp = SUM(Xm2(1, 1:geometry%mg%sg(gLevel2)%num_icell))
  sumRelt2 = 0.0D0
  CALL MPI_ALLREDUCE(sumTemp, sumRelt2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpddGlob%comm, mpddGlob%ierr)

  IF (mpddGlob%isBaseProc()) THEN
    PRINT *, 'Sum is', sumRelt2
  END IF

  CALL geometry%mg%sg(gLevel3)%m_interp_points_on_bndy_linear(vLevel, Xm3)

  CALL mpddGlob%barrier
  ! Destroy the mpdd
  DEALLOCATE (Xm1, Xm2, Xm3)
  IF (mpddGlob%isBaseProc()) DEALLOCATE (XmT)
  IF (mpddGlob%isBaseProc()) THEN
    IF ((sumRelt1 .EQ. 1818) .AND. (sumRelt2 .EQ. 306)) PRINT *, 'Test passed!'
  END IF

  CALL mpddGlob%finalize
END PROGRAM test_geometry
