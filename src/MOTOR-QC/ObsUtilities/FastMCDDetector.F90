!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.ObsUtilities
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yongjian Huang
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yongjian Huang, 2024/10/31, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

MODULE FastMCDDetector_m
  USE ObsOutlierDetector_m
  USE stdlib_sorting
  USE stdlib_linalg
  USE stdlib_stats

  TYPE, PUBLIC, EXTENDS(ObsOutlierDetector_t) :: FastMCDDetector_t
      
      CONTAINS
      PROCEDURE, PUBLIC :: detect => detect_FASTMCD
  END TYPE FastMCDDetector_t

  CONTAINS


  SUBROUTINE detect_FASTMCD(this, outlier_type, obs_data, bck_at_obs, upper_bound, lower_bound, outlier_mask)
    CLASS(FastMCDDetector_t)           :: this
    CHARACTER(len=*) , INTENT(IN)         :: outlier_type
    REAL(r_kind),    INTENT(IN)         :: obs_data(:,:)  !> obs_data(n_obs, n_dim)
    REAL(r_kind),    INTENT(IN)         :: bck_at_obs(:,:) !> bck_at_obs(n_obs, n_dim)
    REAL(r_kind) , INTENT(IN)         :: upper_bound
    REAL(r_kind) , INTENT(IN)         :: lower_bound
    INTEGER(i_kind), ALLOCATABLE, INTENT(INOUT)     :: outlier_mask(:) !> ture is valid
    
    INTEGER(i_kind) :: nObs, nDim 
    INTEGER(i_kind) :: h
    INTEGER(i_kind) :: i, j, it, idx_min, thr
    INTEGER(i_kind), ALLOCATABLE :: h2(:,:), hidx(:,:)
    INTEGER(i_kind) :: htop10(10)
    REAL(r_kind) :: det_s_top10(10)
    REAL(r_kind) :: det_fnl, md_median
    INTEGER(i_kind), ALLOCATABLE :: idx_fnl(:)
    REAL(r_kind), ALLOCATABLE :: dobs(:,:)
    REAL(r_kind), ALLOCATABLE :: t_fnl(:), s_fnl(:,:), md(:), s_mcd(:,:)
    REAL(r_kind) :: lb1, lb2, ub1, ub2
  

    nObs = size(obs_data, 1)
    nDim = size(obs_data, 2)

    ALLOCATE(dobs(nObs, nDim) )


    dobs = obs_data - bck_at_obs
    ! dobs = ABS(obs_data)

    h = NINT((nObs + nDim + 1) / 2.0D0)
   
    !! step1: generate 500 samples randomly
    CALL draw_h(dobs, h, 500, h2)
    
    !! step2: pick 10 samples with minimum covariance determinant
    CALL pick10(dobs, h2, htop10)
    
    !! step3: find the global minimum det iteratively
    ALLOCATE(hidx(10, nObs))
    it = 100
    DO i = 1, 10
      CALL step_it(dobs, h2(htop10(i),:), h, it, hidx(i,:), det_s_top10(i))
    END DO 

    !! step4: get mah dist using the smallest covariance determinant
    ALLOCATE(idx_fnl(nObs), t_fnl(nDim),  s_mcd(nDim, nDim), s_fnl(nDim, nDim), md(nObs))
    idx_min = MINLOC(det_s_top10, 1)
    det_fnl = det_s_top10(idx_min)
    idx_fnl = hidx(idx_min,:)
    t_fnl = MEAN(dobs(idx_fnl,: ), 1)
    s_fnl = COV(dobs(idx_fnl,: ), 1)
    md = mah_dist(dobs, t_fnl, s_fnl)
    CALL SORT_INDEX(md, idx_fnl)
 
    !! step5: Reweight and recal the mah dist, using condition X_0.975
    thr = NINT(nObs * 0.975)

    t_fnl = MEAN(dobs(idx_fnl(1:thr),: ), 1)
    s_fnl = COV(dobs(idx_fnl(1:thr),: ), 1)
    md = mah_dist(dobs, t_fnl, s_fnl)

    ! 1.39 chi dist p=0.5 and fd = 2
    s_mcd = md(NINT(nObs/2.0D0)) * s_fnl / 1.39
    ! print*, 'HYJ+++, tfnl', t_fnl, 's_fnl:', s_fnl,  's_mcd', s_mcd

    
    !! step6: set the outlier mask

    outlier_mask(:) = 1

    lb1 = (t_fnl(1) - SQRT(s_mcd(1,1)))
    ub1 = (t_fnl(1) + SQRT(s_mcd(1,1)))

    lb2 = (t_fnl(2) - SQRT(s_mcd(2,2)))
    ub2 = (t_fnl(2) + SQRT(s_mcd(2,2)))

    ! PRINT*, 'HYJ+++ var1, lb1, ub1, var2, lb2, lb2', lb1, ub1, lb2, ub2

    DO i = 1, nObs
      ! IF (md(i) .GT. md_median) THEN 
      !   print*, 'HYJ md', md(i)
      !   outlier_mask(idx_fnl(i)) = 0
      ! END IF

      IF ( ( dobs(i, 1) .GT. ub1) .OR. (dobs(i, 1) .LT. lb1) &
      .OR. (dobs(i, 2) .GT. ub2) .OR. (dobs(i, 2) .LT. lb2)) THEN
        outlier_mask(i) = 0
      END IF
    END DO

    
    ! DO i = thr, nObs
    !   outlier_mask(idx_fnl(i)) = 0
    ! END DO 

    DEALLOCATE(dobs, hidx, idx_fnl, t_fnl, s_fnl, md, s_mcd)

  END SUBROUTINE detect_FASTMCD


  SUBROUTINE draw_h(x, h, ss, h2)
    REAL(r_kind), DIMENSION(:,:), INTENT(IN)    :: x !> observations data
    INTEGER(i_kind), INTENT(IN)    :: h !> num obs of each samples
    INTEGER(i_kind), INTENT(IN)    :: ss !> size of samples
    INTEGER(i_kind), ALLOCATABLE, INTENT(OUT)  ::h2(:,:) !> samples index

    INTEGER(i_kind) :: i, j, rn
    INTEGER(i_kind) :: nobs, ndim

    nobs = size(x,1)
    ndim = size(x,2)
    ALLOCATE(h2(ss, h))
    xs = 0

    CALL init_random_seed()
    DO i = 1, ss
      DO j = 1, h
          h2(i,j) = gen_rand_int(nobs)
      END DO
    END DO
  END SUBROUTINE


  SUBROUTINE pick10(x, h2, h_top10)
    REAL(r_kind), INTENT(IN)  :: x(:,:)
    INTEGER(i_kind), INTENT(IN)  :: h2(:,:)
    INTEGER(i_kind), INTENT(OUT) :: h_top10(10)

    INTEGER(i_kind) :: h, nDim, nSamples, i, nObs
    REAL(r_kind), ALLOCATABLE :: s(:, :, :), dets(:)
    INTEGER(i_kind), ALLOCATABLE :: hidx(:)

    TYPE(linalg_state_type) :: status
    

    nSamples = size(h2, 1)
    nObs = size(x, 1)
    nDim = size(x, 2)

    ALLOCATE(s(nSamples, nDim, nDim), dets(nSamples))
    ALLOCATE(hidx(nSamples))

    DO i = 1, nSamples
      s(i,:, :) = COV(x(h2(i,:),:), 1)
      ! PRINT*, 'hyj+++ NSAMPLE', s(i,:,:)
      dets(i) = DET(s(i,:,:), err=status)

      IF (status%state .NE. 0) THEN
        dets(i) = invalid
      END IF
    END DO

    CALL sort_index(dets, hidx)
    h_top10(:) = hidx(:10)

    DEALLOCATE(s, dets, hidx)
  END SUBROUTINE

  FUNCTION gen_rand_int(n) RESULT(res)
    INTEGER(i_kind), INTENT(IN)    :: n
    INTEGER(i_kind) :: res

    REAL(r_kind) :: tmp

    CALL random_number(tmp)

    res = NINT(n * tmp + 1)

    IF (res .GT. n ) THEN
      res = n
    END IF

  END FUNCTION

  SUBROUTINE init_random_seed()
    INTEGER(i_kind) :: i, n, clock
    INTEGER(i_kind), ALLOCATABLE :: seed(:)

    CALL random_seed(size=n)
    ALLOCATE(seed(n))

    i = 0

    CALL system_clock(count=clock)

    seed = clock + 37*(/(i-1, i = 1, n)/)

    CALL random_seed(put = seed)

    DEALLOCATE(seed)

  END SUBROUTINE


  SUBROUTINE step_it(x, h2, h, it, hidx, det_s)
    REAL(r_kind), INTENT(IN)  ::x(:,:) 
    INTEGER(i_kind), INTENT(IN) :: h
    INTEGER(i_kind), INTENT(IN) :: it
    INTEGER(i_kind), INTENT(IN) :: h2(h)

    INTEGER(i_kind), INTENT(INOUT) :: hidx(:)
    REAL(r_kind), INTENT(INOUT) :: det_s

    REAL(r_kind) :: det_old, det_new
    INTEGER(i_kind) :: nt, nObs, nDim

    REAL(r_kind), ALLOCATABLE :: t(:), s(:,:)
    REAL(r_kind), ALLOCATABLE :: md(:)
    TYPE(linalg_state_type) :: status

    nObs = SIZE(x, 1)
    nDim = SIZE(x, 2)
    ALLOCATE(md(nObs),  t(nDim), s(nDim, nDim))
    
    nt = 1

    det_old = invalid

    t(:)  = MEAN(x(h2,:), 1)
    s(:,:) = COV(x(h2,:), 1)
    md = mah_dist(x, t, s)

    PRINT*, 'md size:', size(md), size(hidx)

    CALL SORT_INDEX(md, hidx)
    det_new = DET(s)

    DO WHILE ((det_new .LT. det_old) .and. (nt .LT. it))
      det_old = det_new
      
      t(:) = MEAN(x(hidx(1:h), :), 1)
      s(:,:) = COV(x(hidx(1:h), :), 1)

      ! PRINT*, 'hyJ+++ STEP IT ', s, 'h2:', h2
      md = mah_dist(x, t, s)
      CALL SORT_INDEX(md, hidx)
      det_new = DET(s, err=status)

      IF (status%state .NE. 0) THEN
        det_s = invaild
        EXIT
      END IF

      PRINT*, 'iter:', nt, 'det_new:', det_new, 'cov:', s

      nt = nt + 1
    END DO

    det_s = det_new

    DEALLOCATE(t, s, md)
    
  END SUBROUTINE


  ! mah distance
  FUNCTION mah_dist(x, mu, sigma2) RESULT(res)
    REAL(r_kind), INTENT(IN) :: x(:,:)
    REAL(r_kind), INTENT(IN) :: mu(:)
    REAL(r_kind), INTENT(IN) :: sigma2(:,:)

    REAL(r_kind), ALLOCATABLE:: res(:)
    REAL(r_kind), ALLOCATABLE :: sigma2_inv(:,:)
    REAL(r_kind), ALLOCATABLE :: tmp(:), tmp2(:)
    type(linalg_state_type) :: status

    n_elem = size(x,1)
    n_dim = size(x,2)
    
    ALLOCATE(res(n_elem), sigma2_inv(size(sigma2,1), size(sigma2, 2)), tmp(n_dim), tmp2(n_dim))

    sigma2_inv = INV(sigma2, status)

    IF ( status%state .NE. 0 ) THEN 
      DO i = 1, n_elem
        tmp = x(i, :) - mu
        res(i) = invalid
      END DO
    END IF

    DO i = 1, n_elem
      tmp = x(i, :) - mu
      DO j = 1, n_dim
        tmp2(j) = dot_product(sigma2_inv(:, j), tmp)
        
      END DO

      res(i) = dot_product(tmp, tmp2)
    END DO 

    DEALLOCATE(tmp, tmp2, sigma2_inv)

  END FUNCTION 


END MODULE 