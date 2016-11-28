!-------------------------------------------!
!This file is part of the EIGENREC library  !
!-------------------------------------------!
!Vassilis Kalantzis, University of Minnesota!
!Thanos Nikolakopoulos, University of Patras!
!-------------------------------------------!

SUBROUTINE EIGENREC(local_a, local_ia, local_ja, local_nnz, local_rows, local_cols, N, CONV_EIGVECS, Da, nev, tol, msteps, timings, SZ, RANK, COMM)

IMPLICIT NONE
INCLUDE 'mpif.h'

!---------------------!
!Arguments declaration!
!---------------------!
INTEGER, INTENT(IN)                                      :: local_nnz, local_rows, local_cols, N, RANK, msteps, SZ, COMM
INTEGER, INTENT(INOUT)                                   :: nev
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
INTEGER                                                  :: W, I, IERR, INFO, CONV_EIGVECS_COUNT, CHUNK
INTEGER                                                  :: status(MPI_STATUS_SIZE)
INTEGER                                                  :: ISEED

nev = 0D0

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

ENDDO

!--- Extract off-diagonal entries of the tridiagonal matrix
DO I = 1, msteps
   Da(I) = T(I,I)
   IF (I < msteps) THEN
      Db(I) = T(I,I+1)
   ENDIF
ENDDO

!--- Solve the tridiagonal eigenvalue problem and form approximate eigvectors
CALL MPI_BARRIER(COMM, IERR)
time1 = MPI_WTIME()

beta = Db(msteps-1)
CALL DSTEVD('V', msteps, Da, Db, Z, msteps, WORK, LDWORK, IWORK, LIWORK, INFO)
CALL DGEMM("N", "N", CHUNK, msteps, msteps, 1.0D0, V(1:CHUNK,1:msteps), CHUNK, Z, msteps, 0.0D0, Vn, CHUNK)

!--- Trick to check residuals
RES2 = 0.0D0
DO I = 1, msteps
   RES2(msteps-I+1) = ABS(beta)*ABS(Z(msteps,msteps-I+1))
   IF (RES2(msteps-I+1) .LE. tol) THEN
      nev = nev + 1D0
      CONV_EIGVECS(:,nev) = Vn(:,msteps-I+1)
      !IF (RANK == 0) THEN
      !   WRITE(*,*) nev, abs(Z(msteps,msteps-I+1))
      !ENDIF
   ENDIF
ENDDO

CALL MPI_BARRIER(COMM, IERR)
time2 = MPI_WTIME()
timings(7) = time2 - time1

DEALLOCATE(V, T, Z, Vn, y, yold, Ay, t2, WORK, IWORK)

ENDSUBROUTINE EIGENREC

