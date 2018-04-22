MPIF90 = mpiifort -heap-arrays -traceback -g -openmp -lpthread
MPIRUN = mpirun
F90 = ifort

MKL_ROOT = /soft/intel/mkl/10.3.6.233/mkl
MKL_LIBROOT = $(MKL_ROOT)/lib/intel64
MKL_LIB = -L$(MKL_LIBROOT) -lmkl_blas95_lp64 -lmkl_lapack95_lp64  -Wl,--start-group  $(MKL_LIBROOT)/libmkl_intel_lp64.a $(MKL_LIBROOT)/libmkl_intel_thread.a $(MKL_LIBROOT)/libmkl_core.a $(MKL_LIBROOT)/libmkl_sequential.a $(MKL_LIBROOT)/libmkl_scalapack_lp64.a $(MKL_LIBROOT)/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -openmp -lpthread
MKL_SP = -L/home/saady/kalantzi/SPARSKIT2 -lskit

SOURCE1 = additional_subroutines.f90
SOURCE  = driver_EIGENREC.f90
SOURCE2 = mpi_sparse_lanczos_routines.f90
EXE1 = eigenrec.exe

all: build_all

build_all: eigenrec

eigenrec:
	$(MPIF90) $(SOURCE) $(SOURCE1) $(SOURCE2) mmio.f -o $(EXE1) $(MKL_LIB) $(MKL_SP)

clean:
	rm -f *.exe
	
	
