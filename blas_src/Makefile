DIR = src.90/blas_src

FC = ifort

# opzioni per compilazione con ifort
ifeq ($(FC),ifort)
#FOPT =  -ipo -fp-model fast=2  -fp-stack-check -C
FOPT =  -fp-model fast=2
FCOMP = -c -r8 -diag-disable 8290,8291
FLINK =
endif

# opzioni per compilazione con gfortran
ifeq ($(FC),gfortran)
FOPT =  -O3 -ffast-math -funroll-loops 
FCOMP = -c -fno-automatic -fdefault-real-8 -fdefault-double-8 
FLINK = 
endif

FPROF = 
PARALLEL = 

LIB_TOOL = ar crs


BLAS_OBJS = \
   dgemm.o dtrsm.o xerbla.o lsame.o dger.o \
   dcopy.o dgemv.o dtbsv.o daxpy.o ddot.o dgbmv.o dasum.o \
   dtrsv.o dscal.o dswap.o idamax.o \
   strsm.o sgemm.o slamch.o isamax.o \
   sswap.o sscal.o sger.o \
   zaxpy.o zdotc.o zgerc.o zhemv.o zher2.o ztrsv.o

BLAS_LIB = libblas.a


.PHONY : clean all

all: $(BLAS_LIB)


$(BLAS_LIB) : $(BLAS_OBJS)
	$(LIB_TOOL) $(BLAS_LIB) $(BLAS_OBJS)

# regola per convertire .f in .o
.SUFFIXES: .f
.f.o:
	$(FC) $(FCOMP) $(FOPT) $(FPROF) $(PARALLEL) $< 

clean: 
	@rm -f *.o *.a
