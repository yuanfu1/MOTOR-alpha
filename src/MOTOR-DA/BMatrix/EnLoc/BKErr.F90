MODULE BKErr_m

  !>
   !!============================================================================
   !!  This module defines an object for background error covariance modeling.
   !!  It can be used for any ensemble perturbations, single/multi variables,
   !!  cross variables, since this is an abstract method.
   !!
   !!  Reference: see Yuanfu Xie's doc:
   !!    "A New Background Error Modeling Scheme.docx" under statistics dir.
   !!
   !!  Using pseudo inverse method, B matrix is saved through its "useful"
   !!  eigenvectors and eigenvalues.This object provides an interface for users
   !!  multiplying Bx or B^-1 x.
   !!
   !!  \author Yuanfu Xie
   !!  \b Creation April 2015
   !!  \modified by Jiongming Pang (pang.j.m@hotmail.com), April 2022
   !!============================================================================
  !

  ! USE kinds_m, ONLY : r_single,r_kind,i_kind,r_quad,i_long
  USE kinds_m

  TYPE BKErr_t

    !PRIVATE                     ! Default is always private

    CHARACTER(LEN=17) :: header = "BK error module: "
    INTEGER(i_kind)           :: nc        ! Number of control variables
    INTEGER(i_kind)           :: nz        ! Number of non-zero eigenvalues
    LOGICAL                   :: b_ready   ! B-matrix is ready to use
    REAL(r_kind)              :: threshold ! A value for detecting singularity
    REAL(r_kind), ALLOCATABLE :: eigenvalues(:), eigenvectors(:, :)
    REAL(r_kind), ALLOCATABLE :: Q(:, :), R(:, :)
    REAL(r_kind), ALLOCATABLE :: Rinv(:, :)
    REAL(r_kind), ALLOCATABLE :: qr_a(:, :), qr_aux(:)
    REAL(r_kind), ALLOCATABLE :: xe_dr(:)
    REAL(r_kind), ALLOCATABLE :: RinvQinvI(:, :), QRinvTI(:, :)

  CONTAINS
    PROCEDURE, PUBLIC :: init
    PROCEDURE, PUBLIC :: free
    PROCEDURE, PUBLIC :: b_qrDecomposition    ! QR Decomposition, return qr_a and qr_aux.
    PROCEDURE, PUBLIC :: b_qrReconstruction   ! Reconstruction of Q and R matrixes.
    PROCEDURE, PUBLIC :: b_qrDecomp_save      ! Save the qr_a and qr_aux to a file
    PROCEDURE, PUBLIC :: b_qrDecomp_read      ! Read the qr_a and qr_aux
    PROCEDURE, PUBLIC :: b_qrDecomp_Rinv      ! Calculate the inverse of R
    PROCEDURE, PUBLIC :: b_forward            ! calculate FWD with B
    PROCEDURE, PUBLIC :: b_adjoint            ! calculate ADJ with B
    PROCEDURE, PUBLIC :: b_forwardT_save      ! save FWD in terms of transpose
    PROCEDURE, PUBLIC :: b_adjoint_save       ! save ADJ
    PROCEDURE, PUBLIC :: b_adjoint_read_dims  ! read ADJ dims info
    PROCEDURE, PUBLIC :: b_adjoint_read_data  ! read ADJ data

  END TYPE BKErr_t

CONTAINS

  !>
   !!============================================================================
   !!  This routine initializes a background error data type, allocating memory,
   !!  assigning initial values of needed variables.
   !!  Input:
   !!    nv         : an integer specifies how many controls, i.e., the dimension
   !!                 of the background error covariance.
   !!    threshold  : a value to determine the rank of a QR decomposed matrix.
   !!============================================================================
  !
  SUBROUTINE init(this, nc, threshold)
    CLASS(BKErr_t) :: this
    INTEGER(i_kind), INTENT(IN) :: nc
    REAL(r_kind), INTENT(IN)    :: threshold

    this%nc = nc
    this%nz = 0                     ! No eigenvalues determined yet
    this%threshold = threshold
    this%b_ready = .FALSE.

    ! Note: at this point, we do not know how many useful eigenvectors to keep!

  END SUBROUTINE init

  !>
   !!============================================================================
   !!  This routine frees up memory allocated by this data type.
   !!============================================================================
  !
  SUBROUTINE free(this)
    CLASS(BKErr_t) :: this

    ! Deallocate memory:
    IF (ALLOCATED(this%eigenvalues)) DEALLOCATE (this%eigenvalues)
    IF (ALLOCATED(this%eigenvectors)) DEALLOCATE (this%eigenvectors)
    IF (ALLOCATED(this%qr_a)) DEALLOCATE (this%qr_a)
    IF (ALLOCATED(this%qr_aux)) DEALLOCATE (this%qr_aux)
    IF (ALLOCATED(this%Q)) DEALLOCATE (this%Q)
    IF (ALLOCATED(this%R)) DEALLOCATE (this%R)
    ! DEALLOCATE(this%eigenvalues,this%eigenvectors)

  END SUBROUTINE free

  !>
  !!============================================================================
  !!  This routine calculates QR decomposition by using dqrdc in Linpack,
  !!  return qr_a and qr_aux.
  !!
  !!  Created by Jiongming Pang (pang.j.m@hotmail.com), SIMI, Shenzhen
  !!
  !!  Input:
  !!    ne        : number of ensembles
  !!    xe        : ensemble perturbations
  !!    job       : LINPACK parameter requiring colume pivoting, 1/0 for use/no
  !!  Output:
  !!    status    :  0 normal;
  !!                -1 zero rank, which means there is a bug in the code.
  !!      qr_a    : contained the upper triangular matrix and the information
  !!                to recover the whole Q matrix
  !!    qr_aux    : contained the information to recover the whole Q matrix
  !!============================================================================
  !
  SUBROUTINE b_qrDecomposition(this, ne, xe, job, status)
    IMPLICIT NONE
    CLASS(BKErr_t) :: this
    INTEGER(i_kind), INTENT(IN)  :: ne
    REAL(r_kind), INTENT(IN)     :: xe(this%nc, ne)   ! Ensemble perturbations
    INTEGER, INTENT(IN)          :: job
    INTEGER(i_kind), INTENT(OUT) :: status

    ! Local variables:
    INTEGER                       :: i, j, k, num, info
    INTEGER, ALLOCATABLE          :: jpvt(:)
    DOUBLE PRECISION, ALLOCATABLE :: xe_d(:, :), qraux(:), work(:), y(:), qy(:)
    DOUBLE PRECISION, ALLOCATABLE :: v(:, :), IEA(:, :), IE(:, :), tmp(:, :), Rinvtmp(:, :)

    ! Initialize status:
    status = 0

    ! Allocate memory for ensembles:
    ALLOCATE (xe_d(this%nc, ne), qraux(ne), jpvt(ne), work(ne))

    xe_d = xe

    ! LINPACK QR factorization:
    PRINT *, 'Calling LINPACK dqrdc for QR factorization...'
    jpvt = 0
    CALL dqrdc(xe_d, this%nc, this%nc, ne, qraux, jpvt, work, job)
    PRINT *, 'End calling dqrdc with status: ', job
    PRINT *, 'Pivoting: ', jpvt
    PRINT *, 'qauxqr: ', qraux

    ! Determine the rank of R:
    this%nz = 0
    DO i = 1, ne
      PRINT *, "PJM DEBUG --- xe_d(i,i):", i, xe_d(i, i)
      IF (ABS(xe_d(i, i)) .GT. this%threshold) THEN   !this%threshold
        this%nz = this%nz + 1
      ELSE
        EXIT
      END IF
    END DO
    IF (this%nz .EQ. 0) THEN
      status = -1
      PRINT *, "ERROR! The rank of R is 0!"
    ELSE
      ! WRITE(*,"(1I4)") "The rank of R is", this%nz
      PRINT *, "The rank of R is", this%nz
      ALLOCATE (this%qr_a(this%nc, this%nz), this%qr_aux(this%nz))
      this%qr_a = xe_d(:, 1:this%nz)
      this%qr_aux = qraux(1:this%nz)
      this%b_ready = .TRUE.

      ! for log output
      ALLOCATE (this%xe_dr(this%nz))
      DO i = 1, this%nz
        this%xe_dr(i) = xe_d(i, i)
      END DO
    END IF

    ALLOCATE (this%R(this%nz, this%nz))

    this%R = 0.0D0
    DO i = 1, this%nz
      !DO j = i, this%nz
      !  this%R(i, j) = this%qr_a(i, j)
      !END DO
      this%R(i, i:) = this%qr_a(i, 1:this%nz-i+1)
    END DO

    ALLOCATE (Rinvtmp(this%nz, this%nz))
    Rinvtmp = this%R
    CALL DTRTRI('U', 'N', this%nz, Rinvtmp, this%nz, info)
    IF (info .EQ. 0) THEN
      ALLOCATE (this%Rinv(this%nz, this%nz))
      this%Rinv = Rinvtmp
    ELSE
      PRINT *, "ERROR! Failed in calculating the inverse of R matirx."
      STOP
    END IF

    DEALLOCATE (Rinvtmp)

  END SUBROUTINE b_qrDecomposition

  !>
  !!============================================================================
  !!  This routine recover Q and R matrixes by using the outputs for dqrdc, and
  !!  return them.
  !!
  !!  Created by Jiongming Pang (pang.j.m@hotmail.com), SIMI, Shenzhen
  !!
  !!  Input:
  !!    qr_a      : output from dqrdc in routine b_qrDecomposition, contained
  !!                the upper triangular matrix and the information to recover
  !!                the whole Q matrix
  !!    qr_aux    : output from dqrdc in routine b_qrDecomposition, contained
  !!                the information to recover the whole Q matrix
  !!  Output:
  !!    status    :  0 normal;
  !!                -1 zero rank, which means there is a bug in the code.
  !!         Q    : the whole Q matrix
  !!         R    : the whole R matrix
  !!============================================================================
  !
  SUBROUTINE b_qrReconstruction(this, status)
    IMPLICIT NONE
    CLASS(BKErr_t) :: this
    INTEGER(i_kind), INTENT(OUT) :: status

    ! Local variables:
    INTEGER                       :: i, j, k, num
    DOUBLE PRECISION, ALLOCATABLE :: v(:, :), IEA(:, :), IE(:, :), tmp(:, :)

    ! Initialize status:
    status = 0

    ! Revover Q matrix
    ALLOCATE (IEA(this%nc, this%nc), this%Q(this%nc, this%nc), this%R(this%nc, this%nz))

    ! initialize the whole I matrix
    IEA = 0.0D0
    DO i = 1, this%nc
      IEA(i, i) = 1.0D0
    END DO

    DO k = 1, this%nz

      num = this%nc - k + 1
      ALLOCATE (IE(num, num), v(num, 1), tmp(num, num))
      IE = IEA(k:this%nc, k:this%nc)

      ! construct v vector
      v(1, 1) = this%qr_aux(k)
      v(2:num, 1) = this%qr_a(k + 1:this%nc, k)

      tmp = IE - MATMUL(v, TRANSPOSE(v)) / v(1, 1)
      IF (k .EQ. 1) THEN
        this%Q = tmp
      ELSE
        this%Q(:, k:this%nc) = MATMUL(this%Q(:, k:this%nc), tmp)
      END IF

      DEALLOCATE (IE)
      DEALLOCATE (v)
      DEALLOCATE (tmp)

    END DO


    this%R = 0.0D0
    DO i = 1, this%nz
      DO j = i, this%nz
        this%R(i, j) = this%qr_a(i, j)
      END DO
    END DO

    this%b_ready = .TRUE.

  END SUBROUTINE b_qrReconstruction


  SUBROUTINE b_qrDecomp_save(this, fname, bkDiag, status)
    CLASS(BKErr_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: fname
    LOGICAL, INTENT(IN)          :: bkDiag
    INTEGER(i_kind), INTENT(OUT) :: status

    INTEGER(i_kind)              :: i

    status = 0

    IF (this%b_ready) THEN
      OPEN (813, FILE=TRIM(fname), FORM='UNFORMATTED', STATUS='UNKNOWN')
      WRITE (813) this%nc, this%nz
      WRITE (813) this%qr_aux, this%qr_a
      CLOSE (813)

      IF (bkDiag) THEN
        OPEN (813, FILE=TRIM(fname)//'_head', FORM='FORMATTED', STATUS='UNKNOWN')
        WRITE (813, *) this%nc, this%nz
        CLOSE (813)

        OPEN (813, FILE=TRIM(fname)//'_qraux', FORM='FORMATTED', STATUS='UNKNOWN')
        DO i = 1, this%nz
          WRITE (813, 6) this%qr_aux(i)
        END DO
6       FORMAT(e18.8)
        CLOSE (813)

        OPEN (813, FILE=TRIM(fname)//'_qra', FORM='FORMATTED', STATUS='UNKNOWN')
        DO i = 1, this%nc
          WRITE (813, 7) this%qr_a(i, :)
        END DO
7       FORMAT(100E18.8)   !100 is a temporary value
        CLOSE (813)
      END IF
    ELSE
      status = -1       ! Background error is not ready
    END IF

  END SUBROUTINE b_qrDecomp_save

  SUBROUTINE b_qrDecomp_read(this, fname, status)
    IMPLICIT NONE
    CLASS(BKErr_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER(i_kind), INTENT(OUT) :: status

    INTEGER(i_kind)              :: i

    status = 0

    OPEN (813, FILE=TRIM(fname), FORM='UNFORMATTED', STATUS='UNKNOWN')
    READ (813) this%nc, this%nz
    PRINT *, "this%nc,this%nz:", this%nc, this%nz
    IF (.NOT. ALLOCATED(this%qr_aux)) ALLOCATE (this%qr_aux(this%nz))
    IF (.NOT. ALLOCATED(this%qr_a)) ALLOCATE (this%qr_a(this%nc, this%nz))
    READ (813) this%qr_aux, this%qr_a
    CLOSE (813)

  END SUBROUTINE b_qrDecomp_read

  SUBROUTINE b_qrDecomp_Rinv(this)
    IMPLICIT NONE
    CLASS(BKErr_t)  :: this

    INTEGER(i_kind) :: info
    INTEGER(i_kind) :: i, j, status
    DOUBLE PRECISION, ALLOCATABLE :: Rinvtmp(:, :)

    status = 0

    IF (.NOT. ALLOCATED(this%R)) ALLOCATE (this%R(this%nz, this%nz))

    this%R = 0.0D0
    DO i = 1, this%nz
      DO j = i, this%nz
        this%R(i, j) = this%qr_a(i, j)
      END DO
    END DO

    ALLOCATE (Rinvtmp(this%nz, this%nz))
    Rinvtmp = this%R
    CALL DTRTRI('U', 'N', this%nz, Rinvtmp, this%nz, info)
    IF (info .EQ. 0) THEN
      IF (.NOT. ALLOCATED(this%Rinv)) ALLOCATE (this%Rinv(this%nz, this%nz))
      this%Rinv = Rinvtmp
    ELSE
      PRINT *, "ERROR! Failed in calculating the inverse of R matirx."
      STOP
    END IF

    DEALLOCATE (Rinvtmp)

  END SUBROUTINE b_qrDecomp_Rinv

  SUBROUTINE b_forward(this)
    IMPLICIT NONE
    CLASS(BKErr_t) :: this

    REAL(r_kind), ALLOCATABLE :: v(:, :), w(:, :), htw(:, :), RinvQinvI_tmp(:, :), RinvA(:, :)
    REAL(r_kind)              :: vTw
    INTEGER(i_kind)           :: i, j

    REAL(r_kind), ALLOCATABLE :: eye_1(:, :), eye_2(:, :), q1(:, :), q2(:, :)
    REAL(r_kind)              :: q1_q1, q2_q2, q1_q2

    IF (.NOT. ALLOCATED(this%QRinvTI)) ALLOCATE (this%RinvQinvI(this%nz, this%nc))
    ALLOCATE (v(this%nc, 1), w(this%nc, 1), htw(this%nc, 1))
    ALLOCATE (RinvQinvI_tmp(this%nc, 1), RinvA(this%nz, this%nc))

    RinvA = 0.0D0
    RinvA(1:this%nz, 1:this%nz) = this%Rinv

    ! for H_i check
    ALLOCATE (eye_1(this%nc, 1), eye_2(this%nc, 1), q1(this%nc, 1), q2(this%nc, 1))
    eye_1(:, 1) = 0.0D0
    eye_1(1, 1) = 1.0D0
    eye_2(:, 1) = 0.0D0
    eye_2(2, 1) = 1.0D0
    v(:, 1) = 0.0D0
    v(1, 1) = this%qr_aux(1)
    v(2:this%nc, 1) = this%qr_a(2:this%nc, 1)
    q1 = eye_1 - v(1, 1) / this%qr_aux(1) * v
    q2 = eye_2 - v(2, 1) / this%qr_aux(1) * v
    q1_q1 = DOT_PRODUCT(q1(:, 1), q1(:, 1))
    q2_q2 = DOT_PRODUCT(q2(:, 1), q2(:, 1))
    q1_q2 = DOT_PRODUCT(q1(:, 1), q2(:, 1))

    eye_1(:, 1) = 0.0D0
    eye_1(1, 1) = 1.0D0
    eye_2(:, 1) = 0.0D0
    eye_2(2, 1) = 1.0D0
    v(:, 1) = 0.0D0
    v(2, 1) = this%qr_aux(1)
    v(3:this%nc, 1) = this%qr_a(3:this%nc, 1)
    q1 = eye_1 - v(1, 1) / this%qr_aux(1) * v
    q2 = eye_2 - v(2, 1) / this%qr_aux(1) * v
    q1_q1 = DOT_PRODUCT(q1(:, 1), q1(:, 1))
    q2_q2 = DOT_PRODUCT(q2(:, 1), q2(:, 1))
    q1_q2 = DOT_PRODUCT(q1(:, 1), q2(:, 1))

    DO j = 1, this%nc

      w(:, 1) = 0.0D0
      w(j, 1) = 1.0D0

      DO i = 1, this%nz
        v(:, 1) = 0.0D0
        v(i, 1) = this%qr_aux(i)
        v(i + 1:this%nc, 1) = this%qr_a(i + 1:this%nc, i)

        IF (i .EQ. 1) THEN
          htw = w
        END IF

        vTw = DOT_PRODUCT(v(:, 1), htw(:, 1)) / this%qr_aux(i)
        htw = htw - vTw * v
      END DO

      RinvQinvI_tmp = 0.0D0
      RinvQinvI_tmp = MATMUL(RinvA, htw)
      this%RinvQinvI(:, j) = RinvQinvI_tmp(:, 1)

    END DO

    DEALLOCATE (v, w, htw, RinvQinvI_tmp, RinvA)

  END SUBROUTINE b_forward

  SUBROUTINE b_adjoint(this)
    IMPLICIT NONE
    CLASS(BKErr_t) :: this

    REAL(r_kind), ALLOCATABLE :: v(:), w(:), RinvTw(:), HRinvTw(:), RinvA(:, :), qr_aux(:), qr_a(:, :)
    REAL(r_kind)              :: vTw
    REAL(r_kind), ALLOCATABLE :: qr_aux_inv(:), QRinvTI(:,:)
    INTEGER(i_kind)           :: i, j, k, nz, nc
    REAL(r_kind), EXTERNAL    :: DDOT

    nz = this%nz
    nc = this%nc

    ALLOCATE (this%QRinvTI(nc, nz))
    ALLOCATE (v(nc), w(nz)) !, RinvTw(nc))
    ALLOCATE (HRinvTw(nc), RinvA(nc, nz))
    ALLOCATE (qr_aux_inv(nz))

    RinvA = 0.0D0
    RinvA(1:nz, 1:nz) = TRANSPOSE(this%Rinv)
    qr_aux_inv = 1.0D0 / this%qr_aux

    DO j = 1, nz

      w = 0.0D0
      w(j) = 1.0D0
      CALL DGEMM('N', 'N', nc, 1, nz, 1.0D0, RinvA, nc, w, nz, 0.0D0, HRinvTw, nc)

      DO k = 1, nz
        i = nz - k + 1
        v = 0.0D0
        v(i) = this%qr_aux(i)
        v(i+1:nc) = this%qr_a(i+1:nc, i)

        vTw = DOT_PRODUCT(v, HRinvTw) * qr_aux_inv(i)
        CALL DAXPY(nc, -vTw, v, 1, HRinvTw, 1)
      END DO

      this%QRinvTI(:, j) = HRinvTw

    END DO

    DEALLOCATE (v, w, HRinvTw, RinvA, qr_aux_inv)

  END SUBROUTINE b_adjoint

  SUBROUTINE b_forwardT_save(this, fname, status)
    CLASS(BKErr_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER(i_kind), INTENT(OUT) :: status

    INTEGER(i_kind) :: i

    status = 0

    IF (ALLOCATED(this%RinvQinvI)) THEN
      OPEN (813, FILE=TRIM(fname), STATUS='UNKNOWN')
      WRITE (UNIT=813, FMT="(2I10)") this%nz, this%nc
      DO i = 1, this%nc
        WRITE (UNIT=813, FMT="(200F18.8)") this%RinvQinvI(:, i)
      END DO
      CLOSE (813)
    ELSE
      status = -1       ! Background error is not ready
    END IF

  END SUBROUTINE b_forwardT_save

  SUBROUTINE b_adjoint_save(this, fname, status)
    CLASS(BKErr_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER(i_kind), INTENT(OUT) :: status
    INTEGER(i_kind) :: i

    CHARACTER(LEN=1024)          :: fnameT

    status = 0

    IF (ALLOCATED(this%QRinvTI)) THEN
      OPEN (813, FILE=TRIM(fname), ACCESS='STREAM', FORM='UNFORMATTED', STATUS='REPLACE')
      WRITE (813) this%nz, this%nc
      WRITE (813) REAL(this%QRinvTI, kind=4)
      CLOSE (813)
    ELSE
      status = -1       ! Background error is not ready
    END IF

  END SUBROUTINE b_adjoint_save


  SUBROUTINE b_adjoint_read_dims(this, fname)
    CLASS(BKErr_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: fname

    INTEGER(i_kind) :: i

    OPEN (813, FILE=TRIM(fname), ACCESS='STREAM', FORM='UNFORMATTED', STATUS='OLD', IOSTAT=ierr)
    IF (ierr /= 0) THEN
      PRINT*, "Error opening file"
      STOP
    END IF
    READ (UNIT=813) this%nz, this%nc
    CLOSE (813)

  END SUBROUTINE b_adjoint_read_dims

  SUBROUTINE b_adjoint_read_data(this, fname)
    CLASS(BKErr_t) :: this
    CHARACTER(LEN=*), INTENT(IN)   :: fname

    INTEGER(i_kind) :: i, nc, nz
    REAL(4), ALLOCATABLE :: tmp(:,:)

    OPEN (813, FILE=TRIM(fname), ACCESS='STREAM', FORM='UNFORMATTED', STATUS='OLD', IOSTAT=ierr)
    IF (ierr /= 0) THEN
      PRINT*, "Error opening file"
      STOP
    END IF
    READ (UNIT=813) nz, nc
    ALLOCATE (this%QRinvTI(nc, nz))
    ALLOCATE (tmp(nc, nz))
    READ (813) tmp
    this%QRinvTI = REAL(tmp, kind=r_kind)
    CLOSE (813)
    DEALLOCATE (tmp)

  END SUBROUTINE b_adjoint_read_data


END MODULE BKErr_m
