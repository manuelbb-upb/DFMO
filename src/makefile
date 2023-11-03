FC = gfortran
RM = rm -f

FFLAGS = -g
FFLAGS = -O3

EXE = multiobj
OBJS = modules_DFMO.o main.o DFMO.o subroutines_DFMO.o halton.o sobol.o qsortd.o problem.o

all:  $(OBJS)
	$(FC) -o $(EXE) $(OBJS)

.SUFFIXES: .f90 .o
	.f   .o

.f90.o: $* ; $(FC) $(FFLAGS) -c $*.f90

.f.o:   $* ; $(FC) $(FFLAGS) -c $*.f

clean: 
	$(RM) *.o
	$(RM) *.mod
	$(RM) *~
	$(RM) multiobj
	$(RM) pareto_fobs.out
	$(RM) pareto_vars.out
	$(RM) fort.1
	$(RM) fort.2
	$(RM) fort.3
	$(RM) fort.33
	$(RM) fort.55

