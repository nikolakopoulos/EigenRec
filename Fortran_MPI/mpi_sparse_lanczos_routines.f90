!--------------------------------------------------!
!   This file is part of the EIGENREC library      !                                                                                                                                                              
!--------------------------------------------------!                                                                                                                                                             
!Vassilis Kalantzis, Athanasios N. Nikolakopoulos  !                                                                                                                                                             
!          University of Minnesota                 !                                                                                                                                       
!--------------------------------------------------! 

SUBROUTINE EIGENREC(local_a, local_ia, local_ja, local_nnz, local_rows, local_cols, N, CONV_EIGVECS, Da, nev, tol, msteps, timings, SZ, RANK, COMM)

IMPLICIT NONE
INCLUDE 'mpif.h'

!---------------------!
!Arguments declaration!
!---------------------!
INTEGER, INTENT(IN)                                      :: local_nnz, local_rows, local_cols, N, RANK, SZ, COMM
INTEGER, INTENT(INOUT)                                   :: nev, msteps
INTEGER, DIMENSION(local_nnz),          INTENT(IN)       :: local_ja
DOUBLE PRECISION, DIMENSION(local_nnz), INTENT(IN)       :: local_a
INTEGER, DIMENSION(local_rows+1),       INTENT(IN)       :: local_ia
DOUBLE PRECISION, INTENT(IN)                             :: tol
DOUBLE PRECISION, DIMENSION(7), INTENT(OUT)              :: timings
DOUBLE PRECISION, DIMENSION(local_rows,nev), INTENT(OUT) :: CONV_EIGVECS
DOUBLE PRECISION, DIMENSION(msteps),         INTENT(OUT) :: Da

INTEGER, DIMENSION(local_nnz)                            :: local_jao
DOUBLE PRECISION, DIMENSION(local_nnz)                   :: local_ao
INTEGER, DIMENSION(local_cols+1)                         :: local_iao
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)            :: V, T, Z, Vn
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)              :: y, yold, Ay, t1, t2
DOUBLE PRECISION, DIMENSION(msteps)                      :: RES2
DOUBLE PRECISION, DIMENSION(msteps-1)                    :: Db
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)              :: WORK
INTEGER, ALLOCATABLE, DIMENSION(:)                       :: IWORK
INTEGER                                                  :: LDWORK
INTEGER                                                  :: LIWORK
DOUBLE PRECISION                                         :: normy, beta, alpha, time1, time2
DOUBLE PRECISION, DIMENSION(2)                           :: timings2
INTEGER                                                  :: W, I, IERR, INFO, CONV_EIGVECS_COUNT, CHUNK, nev2, conv_flag
INTEGER                                                  :: status(MPI_STATUS_SIZE)
INTEGER                                                  :: ISEED, TEN

nev2 = 0D0
conv_flag = 0D0
TEN = 10D0

LDWORK = 2D0*(1D0 + 4D0*msteps + msteps**2)
LIWORK = 6D0*msteps
CONV_EIGVECS_COUNT = 0D0
timings = 0.0D0
beta = 0.0D0
alpha = 0.0D0
RES2 = 0.0D0
CHUNK = local_rows
CONV_EIGVECS = 0.0D0
local_ao = 0.0D0
local_iao = 0D0
local_jao = 0D0

ALLOCATE(V(CHUNK,msteps+1), T(msteps+1,msteps+1), Z(msteps,msteps), Vn(CHUNK,msteps))
ALLOCATE(y(CHUNK), yold(CHUNK), Ay(CHUNK), t2(CHUNK))
ALLOCATE(WORK(LDWORK), IWORK(LIWORK))

!Z = 0.0D0, T = 0.0D0, V = 0.0D0, Vn = 0.0D0

CALL SRAND(ISEED)
!CALL RANDOM_NUMBER(y)
y = 1.0D0

!--- Normalize starting vector
CALL pddot(CHUNK, y, y, normy, timings2, RANK, COMM)
normy = sqrt(normy)
y = y / normy
V(:,1) = y
yold = 0.0D0


DO I = 1, msteps

   CALL MM_VVT_1D_new(local_a, local_ja, local_ia, y, Ay, timings2, CHUNK, local_cols, local_nnz, RANK, COMM)
   timings(1) = timings(1) + timings2(1); timings(4) = timings(4) + timings2(2)

   !--- Perform the 3-term recurrence
   Ay = Ay - beta*yold
   CALL pddot(CHUNK, Ay, y, alpha, timings2, RANK, COMM)
   timings(3) = timings(3) + timings2(1); timings(5) = timings(5) + timings2(2)
   T(I,I) = alpha
   Ay = Ay - alpha*y

   !--- Perform full reorthogonalization
   CALL MPI_BARRIER(COMM, IERR)
   time1 = MPI_WTIME()
   ALLOCATE(t1(I))
   CALL pbdot2(CHUNK, V(1:CHUNK,1:I), Ay, I, 1, t1, timings2, RANK, COMM)
   CALL DGEMM("N", "N", CHUNK, 1, I, 1.0d0, V(1:CHUNK,1:I), CHUNK, t1, I, 0.0D0, t2, CHUNK)
   Ay = Ay - t2
   DEALLOCATE(t1)
   CALL MPI_BARRIER(COMM, IERR)
   time2 = MPI_WTIME()
   timings(6) = timings(6) + (time2-time1)

   !--- Augment the tridiagonal matrix and continue
   CALL pddot(CHUNK, Ay, Ay, beta, timings2, RANK, COMM)
   timings(3) = timings(3) + timings2(1); timings(5) = timings(5) + timings2(2)
   beta = sqrt(beta)
   yold = y
   y = Ay / beta
   V(1:CHUNK, I+1) = y ! the msteps + 1 parts are not referenced in T, V
   T(I, I+1) = beta
   T(I+1, I) = beta

   IF (MOD(I,TEN) .EQ. 0D0) THEN
     CALL CHECK_CONVERGENCE(T, msteps, I, nev, nev2, conv_flag, tol, RANK)
     IF (conv_flag > 0D0) THEN
        EXIT
     ENDIF
   ENDIF

ENDDO

! # of Lanczos steps performed
msteps = I
IF (conv_flag == 0D0) THEN
  msteps = msteps - 1D0
ENDIF

!IF (RANK == 0) THEN; WRITE(*,*) "Main Lanczos loop is over"; ENDIF

!--- Extract off-diagonal entries of the tridiagonal matrix
DO I = 1, msteps
   Da(I) = T(I,I)
   IF (I < msteps) THEN
      Db(I) = T(I,I+1)
   ENDIF
ENDDO

!IF (RANK == 0) THEN; WRITE(*,*) "Tridiagonal matrix formed"; ENDIF

!--- Solve the tridiagonal eigenvalue problem and form approximate eigvectors
CALL MPI_BARRIER(COMM, IERR)
time1 = MPI_WTIME()

beta = Db(msteps-1)
CALL DSTEVD('V', msteps, Da, Db, Z, msteps, WORK, LDWORK, IWORK, LIWORK, INFO)
CALL DGEMM("N", "N", CHUNK, msteps, msteps, 1.0D0, V(1:CHUNK,1:msteps), CHUNK, Z, msteps, 0.0D0, Vn, CHUNK)
!IF (RANK == 0) THEN; WRITE(*,*) "Tridiagonal eigenvalue problem solved"; ENDIF

!--- Trick to check residuals
RES2 = 0.0D0
nev2 = 0D0
DO I = 1, msteps
   RES2(msteps-I+1) = ABS(beta)*ABS(Z(msteps,msteps-I+1))
   IF (RES2(msteps-I+1) .LE. tol) THEN
      nev2 = nev2 + 1D0
      CONV_EIGVECS(:,nev2) = Vn(:,msteps-I+1)
      !IF (RANK == 0) THEN
      !   WRITE(*,*) nev, RES2(msteps-I+1), Da(msteps-I+1)
      !ENDIF
   ENDIF
   IF (nev2 == nev) THEN
     nev = nev2
     EXIT
   ENDIF
ENDDO
!IF (RANK == 0) THEN; WRITE(*,*) "Residuals checked"; ENDIF

CALL MPI_BARRIER(COMM, IERR)
time2 = MPI_WTIME()
timings(7) = time2 - time1

! Return # of converged eigenpairs
nev = nev2

DEALLOCATE(V, T, Z, Vn, y, yold, Ay, t2, WORK, IWORK)

!IF (RANK == 0) THEN; WRITE(*,*) "Deallocation of resources completed -- returning to main"; ENDIF

ENDSUBROUTINE EIGENREC

!-----------------------------------------

SUBROUTINE CHECK_CONVERGENCE(T, msteps, cur_it, nev, nev2, conv_flag, tol, RANK)

IMPLICIT NONE

INTEGER, INTENT(IN)                                      :: nev, msteps, cur_it, RANK
INTEGER, INTENT(INOUT)                                   :: conv_flag
DOUBLE PRECISION, DIMENSION(msteps+1,msteps+1), INTENT(IN)   :: T
DOUBLE PRECISION, INTENT(IN)                             :: tol

DOUBLE PRECISION, DIMENSION(cur_it)                      :: Da
DOUBLE PRECISION, DIMENSION(cur_it-1)                    :: Db
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)            :: Z
DOUBLE PRECISION, DIMENSION(cur_it)                      :: RES2
INTEGER                                                  :: nev2, I
DOUBLE PRECISION                                         :: beta

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)              :: WORK
INTEGER, ALLOCATABLE, DIMENSION(:)                       :: IWORK
INTEGER                                                  :: LDWORK, LIWORK, INFO

conv_flag = 0D0
nev2      = 0D0
LDWORK    = 2D0*(1D0 + 4D0*cur_it + cur_it**2)
LIWORK    = 6D0*cur_it

ALLOCATE(WORK(LDWORK), IWORK(LIWORK), Z(cur_it,cur_it))

!--- Form vectors Da and Db
DO I = 1, cur_it
   Da(I) = T(I,I)
   IF (I < cur_it) THEN
      Db(I) = T(I,I+1)
   ENDIF
ENDDO

!--- Solve the tridiagonal eigenvalue problem and form approximate eigvectors
beta = Db(cur_it-1)
CALL DSTEVD('V', cur_it, Da, Db, Z, cur_it, WORK, LDWORK, IWORK, LIWORK, INFO)
!--- Trick to check residuals
RES2 = 0.0D0
DO I = 1, cur_it
   RES2(cur_it-I+1) = ABS(beta)*ABS(Z(cur_it,cur_it-I+1))
   IF (RES2(cur_it-I+1) .LE. tol) THEN
      nev2 = nev2 + 1D0
      IF (RANK == 0) THEN
         !WRITE(*,*) nev2, RES2(cur_it-I+1)
      ENDIF
   ENDIF
ENDDO

IF (RANK == 0) THEN; WRITE(*,*) "After:", cur_it, "iterations, # of eigenpairs converged:", nev2; ENDIF

IF (nev2 >= nev) THEN
   conv_flag = 1D0
ENDIF

ENDSUBROUTINE


