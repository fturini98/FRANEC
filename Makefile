DIR = src.90

FC = ifort

# opzioni per compilazione con ifort
ifeq ($(FC),ifort)
#FOPT =  -ipo -fp-model fast=2
FOPT =  -fp-model fast=2
FCOMP = -c -r8 -diag-disable 8290,8291
#FLINK = -vec-report3 -opt-report -opt-report-file=report.txt
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



OBJ = moduli.o Cattura_DM.o mylapack.o mylapack2.o myblas.o  dvode_f90_m.o\
	atm_bh05.o atm_ck03.o atmos.o carbon.o  \
	condint2006.o cross.o cub.o elio.o epsi.o epsig.o\
	idreli.o idro.o \
	chimica.o chim_convettivo.o \
	eqlb.o evolut.o \
	ciacioLi_nstop_novisua.o \
	extkappa.o fato.o fitta.o hbrezo.o   \
	kappa.o io.o local.o  neutr.o optim.o \
	henyey.o innes.o pastem.o maslos.o mixing.o \
	plotta.o quatm.o resnuc.o rilassamento.o routinedoc.o simq.o \
	splin.o stampa.o  state.o stop_evolut.o\
	supera.o sk.o veiove.o xcotrin2006.o zone_convettive.o franec.o

.PHONY : clean all

all: franec

franec : $(OBJ) 
	$(FC) $(FPROF) $(OBJ) $(PARALLEL) -o $@ $(FLINK)
	@ln -fs $(DIR)/$@ -t ..
	@ln -fs $(DIR)/$@ ../pepper

# regola per convertire .f90 in .o
.SUFFIXES: .f90
.f90.o:
	$(FC) $(FCOMP) $(FOPT) $(FPROF) $(PARALLEL) $< 
.SUFFIXES: .f
.f.o:
	$(FC) $(FCOMP) $(FOPT) $(FPROF) $(PARALLEL) $<

clean: 
	@rm -f *.o *.mod franec pepper
	@rm -f ../franec ../pepper 
