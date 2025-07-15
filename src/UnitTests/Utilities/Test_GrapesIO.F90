!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.IO.GrapesIO
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Sanshan Tu, Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Sanshan Tu (tss71618@163.com), 2020/12/31, @SZSC, Shenzhen
!   Modified by Sanshan Tu (tss71618@163.com), 2021/11/01, @SZSC, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! GrapesIO module contains the data type for model_states
!! method.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
program GRAPES_INPUT_TEST
  USE ncWriteGrid_m
  USE GrapesIO_m, ONLY: GrapesIO_t
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  ! USE Filter_m, ONLY: smoothField, guidedfilter

  IMPLICIT NONE
  REAL(r_kind)      :: bgn_time, end_time, value_tmp
  REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE :: value_sin, value_u, value_cmp
  INTEGER           :: ierr, i, j, k
  type(GrapesIO_t)  :: GrapesIO
  INTEGER*8         :: offset_1d_i, offset_1d_j, offset_1d_k
  CHARACTER(LEN=1024) :: ncOutDir, nlFileName, giFileName

  CALL GET_ENVIRONMENT_VARIABLE("INPUT_DIR", nlFileName)
  CALL GET_ENVIRONMENT_VARIABLE("OUTPUT_DIR", ncOutDir)
  CALL GET_ENVIRONMENT_VARIABLE("INPUT_DIR", giFileName)

  ! yaml_get_var(configFile, 'IO', 'input_dir', nlFileName)
  ! yaml_get_var(configFile, 'IO', 'output_dir', ncOutDir)
  ! yaml_get_var(configFile, 'IO', 'input_dir', giFileName)

  nlFileName = TRIM(nlFileName)//"/namelist.input"
  ! giFileName = TRIM(giFileName)//"/grapesinput-4dv-2021053006"
  giFileName = TRIM(ncOutDir)//"/grapesinput-MOTORDA"
  ! giFileName = "/public/home/simi/optest/3DVarVerification/srf/220527/220527_0000/output/grapesinput-MOTORDA"
  giFileName = "/public/home/simi/optest/3DVarVerification/srf/220527/220527_0000/output/grapesinput-MOTORDA"
  giFileName = "/public/home/simi/optest/3DVarVerification/srf/220527/220527_0000/input/model/grapesinput2022052700"
  giFileName = "/public/home/simi/optest/SITest/grapesinput"
  nlFileName = "/public/home/simi/optest/SITest/namelist.input"
  ncOutDir = "/public/home/simi/optest/SITest/"
  ! giFileName = "/Users/zilongqin/source/220527_0000/input/model/grapesinput2022052700"
  nlFileName = ''
  giFileName = ''
  ncOutDir = ''

  CALL getarg(1, giFileName)
  CALL getarg(2, nlFileName)
  CALL getarg(3, ncOutDir)

  PRINT *, 'giFileName: ', TRIM(giFileName)
  PRINT *, 'nlFileName: ', TRIM(nlFileName)
  PRINT *, 'ncOutDir: ', TRIM(ncOutDir)

  IF (TRIM(giFileName) .EQ. '' .OR. TRIM(nlFileName) .EQ. '' &
      .OR. TRIM(ncOutDir) .EQ. '') THEN
    PRINT *, 'Wrong input arguements, test failed!'
    STOP
  END IF

  print *, nlFileName

  CALL CPU_TIME(bgn_time)
  GrapesIO = GrapesIO_t(nlFileName, giFileName, isTest=.true., grapes_model_type="CMA-GD-V3")
  CALL CPU_TIME(end_time)

  ! ALLOCATE (value_sin(GrapesIO%hgrid%ids:GrapesIO%hgrid%ide, &
  !                     GrapesIO%hgrid%kds:GrapesIO%hgrid%kde, &
  !                     GrapesIO%hgrid%jds:GrapesIO%hgrid%jde))

  ! call GrapesIO%read_value("u", value_u)
  ! print *, "u(1, 2, 1), u(2, 2, 1), u(100, 2, 1) = ", &
  !   value_u(1, 2, 1), ", ", value_u(2, 2, 1), ", ", value_u(100, 2, 1)

  ! DO i = GrapesIO%hgrid%ids, GrapesIO%hgrid%ide
  !   DO k = GrapesIO%hgrid%kds, GrapesIO%hgrid%kde
  !     DO j = GrapesIO%hgrid%jds, GrapesIO%hgrid%jde
  !       value_sin(i, k, j) = sin(1.0*i)*sin(1.0*j) + k
  !     END DO
  !   END DO
  ! END DO

  ! call GrapesIO%write_value("u", value_sin)

  ! call GrapesIO%read_value("u", value_cmp)
  ! print *, "u(1, 2, 1), u(2, 2, 1), u(100, 2, 1) = ", &
  !   value_cmp(1, 2, 1), ", ", value_cmp(2, 2, 1), ", ", value_cmp(100, 2, 1)

  ! call GrapesIO%write_value("u", value_u) !write back the original data

  ! DO i = GrapesIO%hgrid%ids, GrapesIO%hgrid%ide
  !   DO k = GrapesIO%hgrid%kds, GrapesIO%hgrid%kde
  !     DO j = GrapesIO%hgrid%jds, GrapesIO%hgrid%jde
  !       if (abs(value_sin(i, k, j) - value_cmp(i, k, j)) > 0.00001) then
  !         print *, "Test Fail"
  !         stop
  !       end if
  !     END DO
  !   END DO
  ! END DO

  ! print *, "Test Pass"

  ! offset_1d_i = (GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1)*4
  ! offset_1d_j = (GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1)*4
  ! offset_1d_k = (GrapesIO%hgrid%kde - GrapesIO%hgrid%kds + 1)*4
  ! PRINT *, "offset_1d", offset_1d_i, offset_1d_j, offset_1d_k, GrapesIO%hgrid%config%num_soil_layers*4

  !here using writegrid3D, output "NetCDF : Invalide dimension ID or name"
  PRINT *, 'SHAPE', SHAPE(GrapesIO%hgrid%u)

  ! DO i = GrapesIO%hgrid%ids, GrapesIO%hgrid%ide
  !   DO j = GrapesIO%hgrid%jds, GrapesIO%hgrid%jde
  !     IF(GrapesIO%hgrid%pi(i, GrapesIO%kds_u+2, j) > GrapesIO%hgrid%pi(i, GrapesIO%kds_u+1, j)) THEN
  !       ! swap = GrapesIO%hgrid%pi(i, GrapesIO%kds_u+2, j)
  !       ! GrapesIO%hgrid%pi(i, GrapesIO%kds_u+2, j) = GrapesIO%hgrid%pi(i, GrapesIO%kds_u+1, j)
  !       ! GrapesIO%hgrid%pi(i, GrapesIO%kds_u+1, j) = swap

  !       ! swap = GrapesIO%hgrid%pip(i, GrapesIO%kds_u+2, j)
  !       ! GrapesIO%hgrid%pip(i, GrapesIO%kds_u+2, j) = GrapesIO%hgrid%pip(i, GrapesIO%kds_u+1, j)
  !       ! GrapesIO%hgrid%pip(i, GrapesIO%kds_u+1, j) = swap
  !       GrapesIO%hgrid%pi(i, GrapesIO%kds_u+1, j) = GrapesIO%hgrid%pi(i, GrapesIO%kds_u+2, j) + 1e-6
  !     ELSE IF (GrapesIO%hgrid%pi(i, GrapesIO%kds_u+2, j) == GrapesIO%hgrid%pi(i, GrapesIO%kds_u+1, j)) THEN
  !       GrapesIO%hgrid%pi(i, GrapesIO%kds_u+1, j) = GrapesIO%hgrid%pi(i, GrapesIO%kds_u+2, j) + 1e-6
  !     END IF

  !     IF(GrapesIO%hgrid%pi(i, GrapesIO%kds_u+1, j) > GrapesIO%hgrid%pi(i, GrapesIO%kds_u, j)) THEN
  !       ! swap = GrapesIO%hgrid%pi(i, GrapesIO%kds_u+1, j)
  !       ! GrapesIO%hgrid%pi(i, GrapesIO%kds_u+1, j) = GrapesIO%hgrid%pi(i, GrapesIO%kds_u, j)
  !       ! GrapesIO%hgrid%pi(i, GrapesIO%kds_u, j) = swap

  !       ! swap = GrapesIO%hgrid%pip(i, GrapesIO%kds_u+1, j)
  !       ! GrapesIO%hgrid%pip(i, GrapesIO%kds_u+1, j) = GrapesIO%hgrid%pip(i, GrapesIO%kds_u, j)
  !       ! GrapesIO%hgrid%pip(i, GrapesIO%kds_u, j) = swap
  !       GrapesIO%hgrid%pi(i, GrapesIO%kds_u, j) = GrapesIO%hgrid%pi(i, GrapesIO%kds_u+1, j) + 1e-6
  !     ELSE IF (GrapesIO%hgrid%pi(i, GrapesIO%kds_u+1, j) == GrapesIO%hgrid%pi(i, GrapesIO%kds_u, j)) THEN
  !       GrapesIO%hgrid%pi(i, GrapesIO%kds_u, j) = GrapesIO%hgrid%pi(i, GrapesIO%kds_u+1, j) + 1e-6
  !     END IF
  !   END DO
  ! END DO
  ! ! PRINT *, 'MINVAL(GrapesIO%hgrid%q): ', MINVAL(GrapesIO%hgrid%q)

  ! CALL GrapesIO%write_value_3d('pi', GrapesIO%hgrid%pi)
  ! PRINT *, 'pi is flushed to the file.'

  BLOCK
    REAL(r_single), ALLOCATABLE :: lat1D(:), lon1D(:)             !< Degree

    ALLOCATE (lat1D(GrapesIO%hgrid%jds:GrapesIO%hgrid%jde), lon1D(GrapesIO%hgrid%ids:GrapesIO%hgrid%ide))
    FORALL (i=GrapesIO%hgrid%ids:GrapesIO%hgrid%ide) lon1D(i) = &
      GrapesIO%hgrid%config%xs_we + (i - 1)*GrapesIO%hgrid%config%xd
    FORALL (j=GrapesIO%hgrid%jds:GrapesIO%hgrid%jde) lat1D(j) = &
      GrapesIO%hgrid%config%ys_sn + (j - 1)*GrapesIO%hgrid%config%yd

    PRINT *, 'Going here.'

    GrapesIO%zSigma_s(0) = GrapesIO%zSigma_s(1) - GrapesIO%zSigma_s(2)
    CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_u.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
                       GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, GrapesIO%hgrid%kde - GrapesIO%hgrid%kds + 1, &
                       lon1D, lat1d, GrapesIO%zSigma_u, "u", GrapesIO%hgrid%u)

    CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_v.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
                       GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, GrapesIO%hgrid%kde - GrapesIO%hgrid%kds + 1, &
                       lon1D, lat1d, GrapesIO%zSigma_u, "v", GrapesIO%hgrid%u)

    CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_q.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
                       GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, GrapesIO%hgrid%kde - GrapesIO%hgrid%kds + 1, &
                       lon1D, lat1d, GrapesIO%zSigma_s, "q", GrapesIO%hgrid%q)

    CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_thp.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
                       GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, GrapesIO%hgrid%kde - GrapesIO%hgrid%kds + 1, &
                       lon1D, lat1d, GrapesIO%zSigma_s, "thp", GrapesIO%hgrid%thp)

    CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_pip.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
                       GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, GrapesIO%hgrid%kde - GrapesIO%hgrid%kds + 1, &
                       lon1D, lat1d, GrapesIO%zSigma_u, "pip", GrapesIO%hgrid%pip)

    CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_th.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
                       GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, GrapesIO%hgrid%kde - GrapesIO%hgrid%kds + 1, &
                       lon1D, lat1d, GrapesIO%zSigma_s, "th", GrapesIO%hgrid%th)

    CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_pi.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
                       GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, GrapesIO%hgrid%kde - GrapesIO%hgrid%kds + 1, &
                       lon1D, lat1d, GrapesIO%zSigma_u, "pi", GrapesIO%hgrid%pi)

    CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_ps.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
                       GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, 1, &
                       lon1D, lat1d, GrapesIO%zSigma_s(1:1), "ps", GrapesIO%hgrid%ps)

    ! BLOCK
    !   REAL(r_kind) :: cc(size(GrapesIO%hgrid%ht, 1), size(GrapesIO%hgrid%ht, 2))

    !   cc = GrapesIO%hgrid%ht*1.0D0
    !   CALL guidedfilter(cc, cc, 3, 0.4D0, cc, GrapesIO, lat1D, lon1D)

    !   ! BLOCK
    !   !   INTEGER(i_kind) :: hei, wid, i, j
    !   !   INTEGER(i_kind) :: r
    !   !   REAL(r_kind) :: imDst(size(GrapesIO%hgrid%ht, 1), size(GrapesIO%hgrid%ht, 2))
    !   !   REAL(r_kind) :: imCum(size(GrapesIO%hgrid%ht, 1), size(GrapesIO%hgrid%ht, 2))
    !   !   !  REAL(r_kind) :: imDst(size( GrapesIO%hgrid%ht, 1), size( GrapesIO%hgrid%ht, 2))
    !   !   REAL(r_single) :: imCuma(size(GrapesIO%hgrid%ht, 1), size(GrapesIO%hgrid%ht, 2))

    !   !   hei = size(GrapesIO%hgrid%ht, 1)
    !   !   wid = size(GrapesIO%hgrid%ht, 2)

    !   !   imDst = 0.0D0
    !   !   r = 3

    !   !   ! cumulative sum over X axis
    !   !   imCum = GrapesIO%hgrid%ht*1.0D0
    !   !   imCuma = imCum*1.0D0
    !   !   CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_imCum0.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
    !   !                      GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, 1, &
    !   !                      lon1D, lat1d, GrapesIO%zSigma_s(1:1), "imCuma", imCuma)

    !   !   DO i = 2, wid
    !   !     imCum(:, i) = imCum(:, i) + imCum(:, i - 1)
    !   !   END DO
    !   !   imCuma = imCum
    !   !   CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_imCum1.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
    !   !                      GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, 1, &
    !   !                      lon1D, lat1d, GrapesIO%zSigma_s(1:1), "imCuma", imCuma)
    !   !   ! difference over X axis
    !   !   imDst(:, 1:r + 1) = imCum(:, 1 + r:2*r + 1); 
    !   !   imDst(:, r + 2:wid - r) = imCum(:, 2*r + 2:wid) - imCum(:, 1:wid - 2*r - 1); 
    !   !   DO j = wid - r + 1, wid
    !   !     imDst(:, j) = imCum(:, wid) - imCum(:, j - r - 1)
    !   !   END DO
    !   !   imCuma = imDst
    !   !   CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_imCum2.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
    !   !                      GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, 1, &
    !   !                      lon1D, lat1d, GrapesIO%zSigma_s(1:1), "imCuma", imCuma)
    !   !   ! cumulative sum over Y axis
    !   !   imCum = imDst
    !   !   DO i = 2, hei
    !   !     imCum(i, :) = imCum(i, :) + imCum(i - 1, :)
    !   !   END DO
    !   !   imCuma = imCum
    !   !   CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_imCum3.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
    !   !                      GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, 1, &
    !   !                      lon1D, lat1d, GrapesIO%zSigma_s(1:1), "imCuma", imCuma)
    !   !   ! difference over Y axis
    !   !   imDst(1:r + 1, :) = imCum(1 + r:2*r + 1, :)
    !   !   imDst(r + 2:hei - r, :) = imCum(2*r + 2:hei, :) - imCum(1:hei - 2*r - 1, :); 
    !   !   DO i = hei - r + 1, hei
    !   !     imDst(i, :) = imCum(hei, :) - imCum(i - r - 1, :)
    !   !   END DO

    !   !   imCuma = imDst
    !   !   CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_imCum4.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
    !   !                      GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, 1, &
    !   !                      lon1D, lat1d, GrapesIO%zSigma_s(1:1), "imCuma", imCuma)

    !   !   ! CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_mean_I.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
    !   !   !                    GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, 1, &
    !   !   !                    lon1D, lat1d, GrapesIO%zSigma_s(1:1), "mean_I", imDst)
    !   ! END BLOCK

    !   GrapesIO%hgrid%ht = cc

    ! END BLOCK

    CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_ht.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
                       GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, 1, &
                       lon1D, lat1d, GrapesIO%zSigma_s(1:1), "ht", GrapesIO%hgrid%ht)
    PRINT * , "GrapesIO%hgrid%zz", shape(GrapesIO%hgrid%zz)
    CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_zz.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
                       GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, GrapesIO%hgrid%kde - GrapesIO%hgrid%kds + 1, &
                       lon1D, lat1d, GrapesIO%zSigma_s, "zz", GrapesIO%hgrid%zz)

    ! CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_wzet.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
    !                    GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, GrapesIO%hgrid%kde - GrapesIO%hgrid%kds + 1, &
    !                    lon1D, lat1d, GrapesIO%zSigma_u, "wzet", GrapesIO%hgrid%wzet)

    ! CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_xland.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
    !                    GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, 1, &
    !                    lon1D, lat1d, GrapesIO%zSigma_u(1:1), "xland", GrapesIO%hgrid%xland)

    ! CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_xice.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
    !                    GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, 1, &
    !                    lon1D, lat1d, GrapesIO%zSigma_u(1:1), "xice", GrapesIO%hgrid%xice)

    ! CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_snowc.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
    !                    GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, 1, &
    !                    lon1D, lat1d, GrapesIO%zSigma_u(1:1), "snowc", GrapesIO%hgrid%snowc)

    DEALLOCATE (lat1D, lon1D)
  END BLOCK

! CONTAINS

!   SUBROUTINE guidedfilter(I, p, r, eps, q, GrapesIO, lat1D, lon1D)
!     INTEGER :: hei, wid
!     INTEGER(i_kind), INTENT(IN) :: r
!     REAL(r_kind), INTENT(IN) :: I(:, :), p(:, :), eps
!     REAL(r_kind), INTENT(OUT) :: q(size(I, 1), size(I, 2))
!     type(GrapesIO_t)  :: GrapesIO
!     REAL(r_single), ALLOCATABLE :: lat1D(:), lon1D(:)             !< Degree

!     REAL(r_kind) :: One(size(I, 1), size(I, 2)), &
!                     NMat(size(I, 1), size(I, 2)), &
!                     mean_I(size(I, 1), size(I, 2)), &
!                     mean_p(size(I, 1), size(I, 2)), &
!                     mean_Ip(size(I, 1), size(I, 2)), &
!                     cov_Ip(size(I, 1), size(I, 2)), &
!                     ! mean_II(size(I, 1), size(I, 2)), &
!                     var_I(size(I, 1), size(I, 2)), &
!                     a(size(I, 1), size(I, 2)), &
!                     b(size(I, 1), size(I, 2)), &
!                     mean_a(size(I, 1), size(I, 2)), &
!                     mean_b(size(I, 1), size(I, 2))

!     REAL(r_single) :: cc(size(I, 1), size(I, 2))

!     hei = size(I, 1)
!     wid = size(I, 2)

!     One = 1.0
!     NMat = boxfilter(One, r)
!     ! PRINT*, 'NMat:', NMat

!     ! PRINT*, 'NMat:', MAXVAL(NMat), MINVAL(NMat)

!     mean_I = boxfilter(I, r)/NMat
!     PRINT *, 'mean_I:', MAXVAL(mean_I), MINVAL(mean_I)

!     mean_p = boxfilter(p, r)/NMat
!     mean_Ip = boxfilter(I*p, r)/NMat
!     cov_Ip = mean_Ip - mean_I*mean_p
!     PRINT *, 'mean_I:', MAXVAL(mean_I), MINVAL(mean_I)
!     PRINT *, 'mean_p:', MAXVAL(mean_p), MINVAL(mean_p)
!     PRINT *, 'mean_Ip:', MAXVAL(mean_Ip), MINVAL(mean_Ip)
!     PRINT *, 'mean_I*mean_p:', MAXVAL(mean_I*mean_p), MINVAL(mean_I*mean_p)

!     cc = mean_I
!     CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_mean_I.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
!                        GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, 1, &
!                        lon1D, lat1d, GrapesIO%zSigma_s(1:1), "mean_I", cc)

!     cc = mean_Ip
!     CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_mean_Ip.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
!                        GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, 1, &
!                        lon1D, lat1d, GrapesIO%zSigma_s(1:1), "mean_Ip", cc)

!     cc = cov_Ip
!     CALL writegrid3D_s(TRIM(ncOutDir)//"/grapesinput_cov_Ip.nc", GrapesIO%hgrid%ide - GrapesIO%hgrid%ids + 1, &
!                        GrapesIO%hgrid%jde - GrapesIO%hgrid%jds + 1, 1, &
!                        lon1D, lat1d, GrapesIO%zSigma_s(1:1), "cov_Ip", cc)

!     PRINT *, 'cov_Ip:', MAXVAL(cov_Ip), MINVAL(cov_Ip)

!     ! mean_II = boxfilter(I*I, r)/N
!     var_I = boxfilter(I*I, r)/NMat - mean_I*mean_I
!     PRINT *, 'var_I:', MAXVAL(var_I), MINVAL(var_I)

!     a = cov_Ip/(var_I + eps)
!     b = mean_p - a*mean_I
!     PRINT *, 'a:', MAXVAL(a), MINVAL(a)
!     PRINT *, 'b:', MAXVAL(b), MINVAL(b)

!     mean_a = boxfilter(a, r)/NMat
!     mean_b = boxfilter(b, r)/NMat
!     PRINT *, 'mean_a:', MAXVAL(boxfilter(a, r)), MINVAL(boxfilter(a, r))
!     PRINT *, 'mean_b:', MAXVAL(boxfilter(b, r)), MINVAL(boxfilter(b, r))
!     PRINT *, 'NMat:', MAXVAL(NMat), MINVAL(NMat)

!     q = mean_a*I + mean_b
!   END SUBROUTINE

!   FUNCTION boxfilter(imSrc, r) RESULT(imDst)
!     INTEGER(i_kind) :: hei, wid, i, j
!     INTEGER(i_kind), INTENT(IN) :: r
!     REAL(r_kind), INTENT(IN) :: imSrc(:, :)
!     REAL(r_kind) :: imDst(size(imSrc, 1), size(imSrc, 2))
!     REAL(r_kind) :: imCum(size(imSrc, 1), size(imSrc, 2))

!     hei = size(imSrc, 1)
!     wid = size(imSrc, 2)

!     imDst = 0.0D0

!     ! cumulative sum over X axis
!     imCum = imSrc
!     DO i = 2, wid
!       imCum(:, i) = imCum(:, i) + imCum(:, i - 1)
!     END DO

!     ! difference over X axis
!     imDst(:, 1:r + 1) = imCum(:, 1 + r:2*r + 1); 
!     imDst(:, r + 2:wid - r) = imCum(:, 2*r + 2:wid) - imCum(:, 1:wid - 2*r - 1); 
!     DO j = wid - r + 1, wid
!       imDst(:, j) = imCum(:, wid) - imCum(:, j - r - 1)
!     END DO

!     ! cumulative sum over Y axis
!     imCum = imDst
!     DO i = 2, hei
!       imCum(i, :) = imCum(i, :) + imCum(i - 1, :)
!     END DO

!     ! difference over Y axis
!     imDst(1:r + 1, :) = imCum(1 + r:2*r + 1, :)
!     imDst(r + 2:hei - r, :) = imCum(2*r + 2:hei, :) - imCum(1:hei - 2*r - 1, :); 
!     DO i = hei - r + 1, hei
!       imDst(i, :) = imCum(hei, :) - imCum(i - r - 1, :)
!     END DO

!   END FUNCTION

end program GRAPES_INPUT_TEST
