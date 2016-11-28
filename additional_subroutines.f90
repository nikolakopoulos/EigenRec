!-------------------------------------------!
!This file is part of the EIGENREC library  !
!-------------------------------------------!
!Vassilis Kalantzis, University of Minnesota!
!Thanos Nikolakopoulos, University of Patras!
!-------------------------------------------!

SUBROUTINE pddot(CHUNK, x, y, glob, timings2, RANK, COMM)
!Element-wise dot product
IMPLICIT NONE
INCLUDE 'mpif.h'

INTEGER, INTENT(IN)                             :: CHUNK, RANK, COMM
DOUBLE PRECISION, DIMENSION(CHUNK), INTENT(IN)  :: x, y
DOUBLE PRECISION, INTENT(OUT)                   :: glob
DOUBLE PRECISION                                :: loc, ddot
INTEGER                                         :: IERR
DOUBLE PRECISION, DIMENSION(2), INTENT(OUT)     :: timings2
DOUBLE PRECISION                                :: t1, t2

loc = 0.0D0;    glob = 0.0D0

CALL MPI_BARRIER(COMM, IERR)
IF (RANK == 0) THEN; t1 = MPI_WTIME(); ENDIF
loc = ddot(CHUNK, x, 1, y, 1)
CALL MPI_BARRIER(COMM, IERR)
IF (RANK == 0) THEN; t2 = MPI_WTIME(); timings2(1) = t2-t1; ENDIF

CALL MPI_BARRIER(COMM, IERR)
IF (RANK == 0) THEN; t1 = MPI_WTIME(); ENDIF
CALL MPI_ALLREDUCE(loc, glob, 1, MPI_DOUBLE_PRECISION, MPI_SUM, COMM, IERR)
CALL MPI_BARRIER(COMM, IERR)
IF (RANK == 0) THEN; t2 = MPI_WTIME(); timings2(2) = t2-t1; ENDIF

ENDSUBROUTINE pddot

SUBROUTINE pbdot2(CHUNK, x, y, col1, col2, glob, timings2, RANK, COMM)
! Block-wise dot product -- here col1 might differ from col2
IMPLICIT NONE
INCLUDE 'mpif.h'

INTEGER, INTENT(IN)                                   :: CHUNK, RANK, col1, col2, COMM
DOUBLE PRECISION, DIMENSION(CHUNK, col1), INTENT(IN)  :: x
DOUBLE PRECISION, DIMENSION(CHUNK, col2), INTENT(IN)  :: y
DOUBLE PRECISION, DIMENSION(col1, col2), INTENT(OUT)  :: glob
DOUBLE PRECISION, DIMENSION(col1, col2)               :: loc
INTEGER                                               :: IERR
DOUBLE PRECISION, DIMENSION(2), INTENT(OUT)           :: timings2
DOUBLE PRECISION                                      :: t1, t2

CALL MPI_BARRIER(COMM, IERR)
IF (RANK == 0) THEN; t1 = MPI_WTIME(); ENDIF
CALL DGEMM("T", "N", col1, col2, CHUNK, 1.0D0, x, CHUNK, y, CHUNK, 0.0D0, loc, col1)
CALL MPI_BARRIER(COMM, IERR)
IF (RANK == 0) THEN; t2 = MPI_WTIME(); timings2(1) = t2-t1; ENDIF

CALL MPI_BARRIER(COMM, IERR)
IF (RANK == 0) THEN; t1 = MPI_WTIME(); ENDIF
CALL MPI_ALLREDUCE(loc, glob, col1*col2, MPI_DOUBLE_PRECISION, MPI_SUM, COMM, IERR)
CALL MPI_BARRIER(COMM, IERR)
IF (RANK == 0) THEN; t2 = MPI_WTIME(); timings2(1) = t2-t1; ENDIF

ENDSUBROUTINE pbdot2


SUBROUTINE MM_VVT_1D_new(av, jav, iav, x, y, timings, m, n, nnz, RANK, COMM)
! Implements the MV product with VV'
! m : number of LOCAL rows
! n : number of GLOBAL cols

IMPLICIT NONE
INCLUDE 'mpif.h'

INTEGER, INTENT(IN)                           :: m, n, nnz, RANK, COMM
INTEGER, INTENT(IN), DIMENSION(m+1)           :: iav
INTEGER, INTENT(IN), DIMENSION(nnz)           :: jav
DOUBLE PRECISION, INTENT(IN), DIMENSION(nnz)  :: av
DOUBLE PRECISION, INTENT(IN), DIMENSION(m)    :: x
DOUBLE PRECISION, INTENT(OUT), DIMENSION(m)   :: y
DOUBLE PRECISION, INTENT(OUT), DIMENSION(2)   :: timings
CHARACTER*1                                   :: transa

DOUBLE PRECISION                              :: t1, t2
DOUBLE PRECISION, DIMENSION(n)                :: ytemp, ytemp2
INTEGER                                       :: I, J, K, IERR
INTEGER, DIMENSION(m)                         :: pntrb, pntre

y      = 0.0D0
ytemp  = 0.0D0
ytemp2 = 0.0D0

! Fix pointers for MKL's MV routine
DO I = 1, m
   pntrb(I) = iav(I)
   pntre(I) = iav(I+1)
ENDDO

! Perform V'*x -- this is local in each processor
CALL MPI_BARRIER(COMM, IERR)
IF (RANK == 0) THEN; t1 = MPI_WTIME(); ENDIF
transa = "T"
CALL mkl_dcsrmv(transa, m, n, 1.D0, "G", av, jav, pntrb, pntre, x, 0.D0, ytemp)
CALL MPI_BARRIER(COMM, IERR)
IF (RANK == 0) THEN; t2 = MPI_WTIME(); timings(1) = timings(1) + t2-t1; ENDIF


! Perform global reduction of V'*x -- this is a "global" op (communication is needed)
IF (RANK == 0) THEN; t1 = MPI_WTIME(); ENDIF
CALL MPI_ALLREDUCE(ytemp, ytemp2, n, MPI_DOUBLE_PRECISION, MPI_SUM, COMM, IERR)
IF (RANK == 0) THEN; t2 = MPI_WTIME(); timings(2) = timings(2) + t2-t1; ENDIF

! Perform V*x  -- this is local in each processor
CALL MPI_BARRIER(COMM, IERR)
IF (RANK == 0) THEN; t1 = MPI_WTIME(); ENDIF
transa = "N"
CALL mkl_dcsrmv(transa, m, n, 1.D0, "G", av, jav, pntrb, pntre, ytemp2, 0.D0, y)
CALL MPI_BARRIER(COMM, IERR)
IF (RANK == 0) THEN; t2 = MPI_WTIME(); timings(1) = timings(1) + t2-t1; ENDIF

ENDSUBROUTINE MM_VVT_1D_new



SUBROUTINE csr2csc(n, n2, nnz, a, ja, ia, ao, jao, iao)
! Converts CSR to CSC sparse storage format -- 
! this routine was used in an earlier phase, 
! when we wrote our own MV routine -- after 
! switching to MKL, this routine is not used 
! in EIGENREC anymore.

  IMPLICIT NONE

  INTEGER, INTENT(IN)                              :: n, n2, nnz
  INTEGER, DIMENSION(nnz),          INTENT(INOUT)  :: ja, jao
  DOUBLE PRECISION, DIMENSION(nnz), INTENT(INOUT)  :: a, ao
  INTEGER, DIMENSION(n+1),          INTENT(INOUT)  :: ia
  INTEGER, DIMENSION(n2+1),         INTENT(INOUT)  :: iao

  INTEGER                                          :: i, j, k, next, ipos

  ipos = 1D0

  DO i = 1, n2+1
     iao(i) = 0
  ENDDO

  DO i = 1, n
     DO k = ia(i), ia(i+1) - 1
        j = ja(k)+1
        iao(j) = iao(j)+1
     ENDDO
  ENDDO
  
  iao(1) = ipos
  DO i=1,n2
     iao(i+1) = iao(i) + iao(i+1)
  ENDDO

  DO i=1,n
     DO k=ia(i),ia(i+1)-1
        j = ja(k)
        next = iao(j)
        ao(next) = a(k)
        jao(next) = i
        iao(j) = next+1
     ENDDO
  ENDDO

  DO i=n2,1,-1
     iao(i+1) = iao(i)
  ENDDO
  iao(1) = ipos

ENDSUBROUTINE csr2csc
