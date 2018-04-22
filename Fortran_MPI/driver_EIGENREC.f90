!--------------------------------------------------!
!   This file is part of the EIGENREC library      !                                                                                                                                                              
!--------------------------------------------------!                                                                                                                                                             
!Vassilis Kalantzis, Athanasios N. Nikolakopoulos  !                                                                                                                                                             
!          University of Minnesota                 !                                                                                                                                       
!--------------------------------------------------! 

PROGRAM driver_eigenrec

IMPLICIT NONE
INCLUDE 'mpif.h'

!-------------------------------------------!                                                                                                                                                             
!This file is part of the EIGENREC library  !                                                                                                                                                             
!-------------------------------------------!                                                                                                                                                             
!Vassilis Kalantzis, Thanos Nikolakopoulos  !                                                                                                                                                             
!          University of Minnesota          !                                                                                                                                                             
!-------------------------------------------! 

!-------------------------------------------!
!Date: 09 / 13 / 2017, Minneapolis, MN, USA !
!-------------------------------------------!

!----------------!
!Define variables!
!----------------!
INTEGER                                           :: N, len, IERR, atoi, status, I, J, RANK, SIZE, dummy, nev, nev2, lanczos_steps
INTEGER                                           :: choice, counter, lansteps_init
DOUBLE PRECISION                                  :: atof, tol, t1, t2, time_2_read_matrix
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE       :: timings_lanczos, Da
CHARACTER(len=32)                                 :: BUFFER, infile, name
INTEGER                                           :: nzmax, dest, job, sbw
INTEGER                                           :: info, m, damping, m2, chunk_nnz
CHARACTER(len=10)                                 :: rep
CHARACTER(len=7)                                  :: field
CHARACTER(len=19)                                 :: symm
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE       :: a, a_gen, a_gen2, IDummy, CDummy, iwork, iw, Av
INTEGER, DIMENSION(:), ALLOCATABLE                :: Ai, Aj, ia, ja, ia_gen, ja_gen, ia_gen2, ja_gen2
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE     :: W

INTEGER, DIMENSION(:), ALLOCATABLE                :: local_ia, local_ia_per_proc_last, local_ia_per_proc
INTEGER, DIMENSION(:), ALLOCATABLE                :: local_ja, local_ja_per_proc, local_ja_per_proc_last
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE       :: local_a, local_a_per_proc, local_a_per_proc_last
INTEGER                                           :: row_wise, ALL_CHUNKS_SAME, rows_to_send, root_proc_rows_to_send, rows_to_send_last
INTEGER                                           :: local_nnz, local_nnz_per_proc_last, local_rows
INTEGER, DIMENSION(:), ALLOCATABLE                :: local_nnz_per_proc

! For nnz-aware partitioning
INTEGER, DIMENSION(:), ALLOCATABLE                :: nnz_wise_start_index, nnz_wise_rows_size, nnz_wise_nnz
INTEGER                                           :: row_start, counter2, set, J2

!----------------------------------------------------!
!Choice to distribute the matrix row-wise or nnz-wise!
!----------------------------------------------------!
row_wise = 0D0

!--------------------------!
!Initialize MPI environment!
!--------------------------!
CALL MPI_INIT(IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE, IERR)

IF (RANK==0) THEN
WRITE(*,*) " "
WRITE(*,*) "/----------------------/"
WRITE(*,*) "EIGENREC library, Version 1.1"
WRITE(*,*) "For more details, see the paper: 'EIGENREC: An Efficient Latent Factor Framework for Top-N Recommendations'"
WRITE(*,*) "Send any any comments or bug reports to: anikolak@umn.edu, or kalan019@umn.edu"
WRITE(*,*) "/----------------------/"
WRITE(*,*) " "
ENDIF

ALLOCATE(local_nnz_per_proc(SIZE))
local_nnz_per_proc = 0D0

lanczos_steps = 0D0

!---------------------------!
!Take command line arguments!
!---------------------------!
CALL GETARG(1, name)

CALL GETARG(2, BUFFER)
READ (BUFFER, "(i)") nev

CALL GETARG(3, BUFFER)
READ (BUFFER, "(i)") lanczos_steps

IF (lanczos_steps < nev) THEN
   WRITE(*,*) "max # of lanczos steps is smaller than 'nev'. Aborting..."
   STOP
ENDIF

!-----------------------------------------!
!Root processor reads the iteration matrix!
!-----------------------------------------!
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
t1 = MPI_WTIME()
IF (RANK .EQ. 0) THEN
   open(100, file = name, status = 'old')
   CALL mminfo(100, rep, field, symm, N, m2, nzmax)

   WRITE(*,*) " "
   WRITE(*,*) "/----------------------/"
   WRITE(*,*) "Matrix: ", name
   WRITE(*,*) "Symmetry: ", symm
   WRITE(*,*) "rep: ", rep
   WRITE(*,*) "field: ", field
   WRITE(*,*) "ROWS: ", N
   WRITE(*,*) "COLUMNS: ", m2
   WRITE(*,*) "NNZ: ", nzmax
   WRITE(*,*) "/----------------------/"
   WRITE(*,*) " "

   ALLOCATE(Ai(nzmax), Aj(nzmax), Av(nzmax))
   ALLOCATE(IDummy(nzmax))
   ALLOCATE(CDummy(nzmax))
   Ai = 0.0D0; Aj = 0.0D0; Av = 0.0D0

   CALL mmread(100, rep, field, symm, N, m2, nzmax, nzmax, Ai, Aj, IDummy, Av, CDummy)
   ALLOCATE(ia(N+1), ja(nzmax), a(nzmax))
   CALL coocsr(N, nzmax, Av, Ai, Aj, a, ja, ia)
   ALLOCATE(iwork(max(N+1, 2*nzmax)))

   IF (symm .EQ. 'symmetric') THEN 

      ALLOCATE(a_gen(nzmax), ia_gen(N+1), ja_gen(nzmax))

      CALL coocsr(N, nzmax, Av, Aj, Ai, a_gen, ja_gen, ia_gen)
      ALLOCATE(a_gen2(2*nzmax), ia_gen2(N+1), ja_gen2(2*nzmax), iw(N))
      dest = 2*nzmax
      job = 1

      CALL aplb(N, m2, job, a, ja, ia, a_gen, ja_gen, ia_gen, a_gen2, ja_gen2, ia_gen2, dest, iw, info)

      DO i = 1, N
         DO j = ia_gen2(i), ia_gen2(i+1) - 1
            IF (ja_gen2(j) .eq. i) then
               a_gen2(j) = a_gen2(j)/2
            ENDIF
         ENDDO
      ENDDO

      DEALLOCATE(ja_gen, a_gen, ia_gen, ja, a)

      nzmax = ia_gen2(N + 1) - 1
      ALLOCATE(a(nzmax), ja(nzmax))

      ia(1:N+1) = ia_gen2(1:N+1)
      ja(1:nzmax) = ja_gen2(1:nzmax)
      a(1:nzmax) = a_gen2(1:nzmax)

      DEALLOCATE(ia_gen2, ja_gen2, a_gen2, iwork, iw)

   ENDIF

   DEALLOCATE(IDummy)
   DEALLOCATE(CDummy)

ENDIF

!--- Distribute #rows of matrix
CALL MPI_BCAST(N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)

!--- Distribute #cols of matrix
CALL MPI_BCAST(m2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, IERR)


!-----------------------------------------------!
!Root processor distributes the iteration matrix!
!-----------------------------------------------!
ALLOCATE( nnz_wise_start_index(SIZE), nnz_wise_rows_size(SIZE), nnz_wise_nnz(SIZE) );

IF (RANK == 0)THEN

   IF (row_wise .EQ. 1D0) THEN
      !---------------------------------------------!
      !Find how many rows to distribute to each proc!
      !---------------------------------------------!
      ALL_CHUNKS_SAME = 0D0

      IF (MOD(N,SIZE) .EQ. 0D0) THEN ! If N divides # PROCS (SIZE) evenly then all PROCS get the same # or rows
         ALL_CHUNKS_SAME = 1D0
      ENDIF

      rows_to_send = N / SIZE ! Rows to each PROC
      root_proc_rows_to_send = rows_to_send ! Rows for the root PROC

      DO I = 1, SIZE-1
         IF (I .LE. SIZE-2) THEN ! Send corresponding number of rows to each PROC
            CALL MPI_SEND(rows_to_send, 1, MPI_INTEGER, I, I, MPI_COMM_WORLD, IERR)
         ELSE
            IF (ALL_CHUNKS_SAME .EQ. 1D0) THEN ! The number of rows of the last PROC depends on ALL_CHUNKS_SAME
               rows_to_send_last = rows_to_send
               CALL MPI_SEND(rows_to_send_last, 1, MPI_INTEGER, I, I, MPI_COMM_WORLD, IERR) ! last proc takes equal number of rows
            ELSE
               rows_to_send_last = N - (SIZE-1)*rows_to_send ! Last PROC gets different number of rows
               CALL MPI_SEND(rows_to_send_last, 1, MPI_INTEGER, I, I, MPI_COMM_WORLD, IERR) ! last proc takes unequal number of rows
            ENDIF
         ENDIF
      ENDDO

      !----------------------!
      !Send local row indices!
      !----------------------!
      ALLOCATE(local_ia(root_proc_rows_to_send+1))
      ALLOCATE(local_ia_per_proc(rows_to_send+1))
      ALLOCATE(local_ia_per_proc_last(rows_to_send_last+1))

      local_ia = ia(1:root_proc_rows_to_send+1)
      DO I = 1, SIZE-1
         IF (I .LE. SIZE-2) THEN
            local_ia_per_proc = ia(I*rows_to_send + 1:(I+1)*rows_to_send + 1)
            CALL MPI_SEND(local_ia_per_proc, rows_to_send + 1, MPI_INTEGER, I, I+SIZE, MPI_COMM_WORLD, IERR)
         ELSE
            local_ia_per_proc_last = ia(I*rows_to_send+1:N+1)
            CALL MPI_SEND(local_ia_per_proc_last, rows_to_send_last + 1, MPI_INTEGER, I, I+SIZE, MPI_COMM_WORLD, IERR)
         ENDIF
      ENDDO

      !----------------------------------------------------!
      !Find how many nnz entries are allocated to each proc!
      !----------------------------------------------------!
      local_nnz = ia(root_proc_rows_to_send + 1) - ia(1)
      local_nnz_per_proc(1) = local_nnz
      ALLOCATE(local_ja(local_nnz), local_a(local_nnz))

      DO I = 1, SIZE-1
         IF (I .LE. SIZE-2) THEN
            local_nnz_per_proc(I+1) = ia((I+1)*rows_to_send + 1) - ia(I*rows_to_send + 1)
            CALL MPI_SEND(local_nnz_per_proc(I+1), 1, MPI_INTEGER, I, I+2*SIZE, MPI_COMM_WORLD, IERR)
         ELSE
            local_nnz_per_proc_last = ia(N+1) - ia(I*rows_to_send + 1)
            local_nnz_per_proc(SIZE) = local_nnz_per_proc_last
            CALL MPI_SEND(local_nnz_per_proc_last, 1, MPI_INTEGER, I, I+2*SIZE, MPI_COMM_WORLD, IERR)
         ENDIF
      ENDDO

      !---------------------!
      !Now send local values!
      !---------------------!
      local_a = a(1:local_nnz)

      DO I = 1, SIZE-1
         IF (I .LE. SIZE-2) THEN
            CALL MPI_SEND(a(ia(I*rows_to_send + 1):ia((I+1)*rows_to_send+1)-1), local_nnz_per_proc(I+1), MPI_DOUBLE_PRECISION, I, I+3*SIZE, MPI_COMM_WORLD, IERR)
         ELSE
            CALL MPI_SEND(a(ia(I*rows_to_send + 1):ia(N+1)-1), local_nnz_per_proc_last, MPI_DOUBLE_PRECISION, I, I+3*SIZE, MPI_COMM_WORLD, IERR)
         ENDIF
      ENDDO

      !-----------------------------!
      !Now send local column indices!
      !-----------------------------!
      local_ja = ja(1:local_nnz)

      DO I = 1, SIZE-1
         IF (I .LE. SIZE-2) THEN
            CALL MPI_SEND(ja(ia(I*rows_to_send + 1):ia((I+1)*rows_to_send+1)-1), local_nnz_per_proc(I+1), MPI_INTEGER, I, I+4*SIZE, MPI_COMM_WORLD, IERR)
         ELSE
            CALL MPI_SEND(ja(ia(I*rows_to_send + 1):ia(N+1)-1), local_nnz_per_proc_last, MPI_INTEGER, I, I+4*SIZE, MPI_COMM_WORLD, IERR)
         ENDIF
      ENDDO

      local_rows = root_proc_rows_to_send

   ELSE
      !---------------------------------------------!
      !Here I should treat the nnz-wise distribution!
      !---------------------------------------------!
      chunk_nnz = nzmax / SIZE
      row_start = 1D0
      counter2 = 0D0
      set = 0D0
      J2 = 0D0

      DO I = 1, SIZE-1
         nnz_wise_start_index(I) = row_start
         DO J = row_start, N
            IF (set == 0D0) THEN
            counter2 = counter2 + ia(J+1) - ia(J)
            IF ( counter2 >= chunk_nnz .AND. set == 0D0 ) THEN
               nnz_wise_nnz(I) = counter2
               nnz_wise_rows_size(I) = J - row_start + 1D0
               set = 1D0
               J2 = J
            ENDIF
            ENDIF
         ENDDO
         set = 0D0
         counter2 = 0D0
         row_start = J2+1
      ENDDO
      nnz_wise_start_index(SIZE) = row_start
      nnz_wise_rows_size(SIZE) = N - row_start + 1D0
      nnz_wise_nnz(SIZE) = ia(N+1) - ia(row_start)

      local_rows = nnz_wise_rows_size(1)
      DO I = 1, SIZE-1 ! Send corresponding number of rows to each PROC
         CALL MPI_SEND(nnz_wise_rows_size(I+1), 1, MPI_INTEGER, I, I, MPI_COMM_WORLD, IERR)
      ENDDO
      ALLOCATE( local_ia(local_rows+1) )
      local_ia = ia(1:local_rows+1)
      DO I = 1, SIZE-1 ! Send the rows
         CALL MPI_SEND(ia(nnz_wise_start_index(I+1):nnz_wise_start_index(I+1)+nnz_wise_rows_size(I+1)), nnz_wise_rows_size(I+1) + 1, MPI_INTEGER, I, I+SIZE, MPI_COMM_WORLD, IERR)
      ENDDO
      local_nnz = nnz_wise_nnz(1)
      DO I = 1, SIZE-1 ! Send the number of non-zeros
         CALL MPI_SEND(nnz_wise_nnz(I+1), 1, MPI_INTEGER, I, I+2*SIZE, MPI_COMM_WORLD, IERR)
      ENDDO
      ALLOCATE( local_a(local_nnz) )
      local_a = a(1:local_nnz)
      DO I = 1, SIZE-1 ! Send values
         CALL MPI_SEND(a(ia(nnz_wise_start_index(I+1)):ia(nnz_wise_start_index(I+1)+nnz_wise_rows_size(I+1))-1), nnz_wise_nnz(I+1), MPI_DOUBLE_PRECISION, I, I+3*SIZE, MPI_COMM_WORLD, IERR)
      ENDDO
      ALLOCATE( local_ja(local_nnz) )
      local_ja = ja(1:local_nnz)
      DO I = 1, SIZE-1 ! Send column indices
         CALL MPI_SEND(ja(ia(nnz_wise_start_index(I+1)):ia(nnz_wise_start_index(I+1)+nnz_wise_rows_size(I+1))-1), nnz_wise_nnz(I+1), MPI_INTEGER, I, I+4*SIZE, MPI_COMM_WORLD, IERR)
      ENDDO

   ENDIF


ELSE


   IF (row_wise .EQ. 1D0) THEN

      CALL MPI_RECV(local_rows, 1, MPI_INTEGER, 0, RANK, MPI_COMM_WORLD, status, IERR) ! number of rows
      ALLOCATE( local_ia(local_rows+1) )
      local_ia = 0D0
      CALL MPI_RECV(local_ia, local_rows + 1, MPI_INTEGER, 0, RANK+SIZE, MPI_COMM_WORLD, status, IERR)
      CALL MPI_RECV(local_nnz, 1, MPI_INTEGER, 0, RANK+2*SIZE, MPI_COMM_WORLD, status, IERR) ! number of nnz
      ALLOCATE( local_ja(local_nnz) )
      ALLOCATE( local_a(local_nnz) )
      local_ja = 0D0
      local_a = 0.0D0
      CALL MPI_RECV(local_a, local_nnz, MPI_DOUBLE_PRECISION, 0, RANK+3*SIZE, MPI_COMM_WORLD, status, IERR)
      CALL MPI_RECV(local_ja, local_nnz, MPI_INTEGER, 0, RANK+4*SIZE, MPI_COMM_WORLD, status, IERR)

      !-----------------------------------------------------------!
      !The non-root processors should refix their local_ia indices!
      !-----------------------------------------------------------!
      DO I = 2, local_rows + 1
         local_ia(I) = local_ia(I) - local_ia(1)
         local_ia(I) = local_ia(I) + 1D0
      ENDDO

      local_ia(1) = 1D0

   ELSE

      !---------------------------------------------!
      !Here I should treat the nnz-wise distribution!
      !---------------------------------------------!
      CALL MPI_RECV(local_rows, 1, MPI_INTEGER, 0, RANK, MPI_COMM_WORLD, status, IERR) ! number of rows
      ALLOCATE( local_ia(local_rows+1) )
      local_ia = 0D0
      CALL MPI_RECV(local_ia, local_rows + 1, MPI_INTEGER, 0, RANK+SIZE, MPI_COMM_WORLD, status, IERR)
      CALL MPI_RECV(local_nnz, 1, MPI_INTEGER, 0, RANK+2*SIZE, MPI_COMM_WORLD, status, IERR) ! number of nnz
      ALLOCATE( local_ja(local_nnz) )
      ALLOCATE( local_a(local_nnz) )
      local_ja = 0D0
      local_a = 0.0D0
      CALL MPI_RECV(local_a, local_nnz, MPI_DOUBLE_PRECISION, 0, RANK+3*SIZE, MPI_COMM_WORLD, status, IERR)
      CALL MPI_RECV(local_ja, local_nnz, MPI_INTEGER, 0, RANK+4*SIZE, MPI_COMM_WORLD, status, IERR)

      !-----------------------------------------------------------!
      !The non-root processors should refix their local_ia indices!
      !-----------------------------------------------------------!
      DO I = 2, local_rows + 1
         local_ia(I) = local_ia(I) - local_ia(1)
         local_ia(I) = local_ia(I) + 1D0
      ENDDO

      local_ia(1) = 1D0

   ENDIF

ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
t2 = MPI_WTIME()
time_2_read_matrix = t2-t1
!--------------------------!
!End of matrix distribution!
!--------------------------!

!------------------------!
!Distributed Lanczos test!
!------------------------!
tol = 1.0D-8                     ! tolerance for the sought eigenpairs

choice = 1D0
IF (choice == 1D0) THEN
   ALLOCATE(timings_lanczos(7))  ! a vector that keeps track of the amount of time spent in different tasks
   timings_lanczos = 0.0D0

   !nev = 50D0                    ! # of eigenvectors sought
   nev2 = nev
   IF (lanczos_steps .EQ. 0D0) THEN
      lanczos_steps = 50D0       ! # of Lanczos steps
   ENDIF
   ALLOCATE(W(local_rows, nev))  ! buffer to keep the approximate eigenvectors
   ALLOCATE(Da(lanczos_steps))   ! buffer to kep the Ritz values

   CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
   t1 = MPI_WTIME()
   CALL EIGENREC(local_a, local_ia, local_ja, local_nnz, local_rows, m2, N, W, Da, nev, tol, lanczos_steps, timings_lanczos, SIZE, RANK, MPI_COMM_WORLD)
   CALL MPI_BARRIER(MPI_COMM_WORLD, IERR)
   t2 = MPI_WTIME()

   IF (RANK == 0) THEN
      WRITE(*, *) " "
      WRITE(*, *) "/--------------------------------------------/"
      WRITE(*, *) "# of MPI processes: ", SIZE
      WRITE(*, *) "Matrix distribution (1: by rows, 0: by # of nonzeros): ", row_wise
      WRITE(*, *) "Number of MV products with V (2*Lanczos iterations): ", 2*lanczos_steps
      WRITE(*, *) "Number of eigenpairs requested: ", nev2
      WRITE(*, *) "Number of eigenpairs converged: ", nev
      WRITE(*, *) "Accuracy of each eigenpair sought: ", tol
      WRITE(*, *) "/--------------------------------------------/"
      WRITE(*, *) "Time spent on reading and distributing users-items matrix: ", time_2_read_matrix
      WRITE(*, *) "Time spent on Lanczos: ", t2 - t1, "seconds."
      WRITE(*, *) "Time spent on performing MV products: ", timings_lanczos(1)
      WRITE(*, *) "Time spent on communication for the MV products: ", timings_lanczos(4)
      WRITE(*, *) "Time spent on performing DOTS: ", timings_lanczos(3)
      WRITE(*, *) "Time spent on communication for the DOTS: ", timings_lanczos(5)
      WRITE(*, *) "Time spent on performing orthogonalization: ", timings_lanczos(6)
      WRITE(*, *) "Time spent on the Rayleigh-Ritz procedure: ", timings_lanczos(7)
      WRITE(*, *) "/--------------------------------------------/"
      WRITE(*, *) " "
   ENDIF

   DEALLOCATE(W, timings_lanczos, Da)
ENDIF


!--------------------!
!Deallocate resources!
!--------------------!
DEALLOCATE(local_a, local_ia, local_ja)


!------------!
!Finalize MPI!
!------------!
CALL MPI_FINALIZE(IERR)


END PROGRAM driver_eigenrec
