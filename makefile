FC = gfortran
RM = rm -f

FFLAGS = -O3 -g
FFLAGS_SHARED = $(FFLAGS) -fPIC

EXE = multiobj
EXE_SHARED = $(EXE).so
OBJS =  main.o opt_multiobj.o DFMO.o subroutines_DFMO.o halton.o sobol.o qsortd.o
OBJS_SHARED = opt_multiobj.fpic.o DFMO.fpic.o subroutines_DFMO.fpic.o halton.fpic.o sobol.fpic.o qsortd.fpic.o 

all: $(EXE)

$(EXE): modules_DFMO.o $(OBJS) problem.o
	$(FC) -o $@ $^

modules_DFMO.o: problem.o modules_DFMO.f90
	$(FC) -c $(FFLAGS) modules_DFMO.f90 -o modules_DFMO.o

problem.o: problem_shim.f90
	$(FC) -c $(FFLAGS) $< -o $@

$(OBJS): %.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

shared: $(EXE_SHARED)

$(EXE_SHARED): modules_DFMO.fpic.o $(OBJS_SHARED) problem.fpic.o
	$(FC) -shared -o $@ $^

modules_DFMO.fpic.o: problem.fpic.o modules_DFMO.f90
	$(FC) -c $(FFLAGS_SHARED) modules_DFMO.f90 -o modules_DFMO.fpic.o

problem.fpic.o: problem_standalone.f03
	$(FC) -c $(FFLAGS_SHARED) $< -o $@

$(OBJS_SHARED): %.fpic.o: %.f90
	$(FC) -c $(FFLAGS_SHARED) $< -o $@

clean: 
	$(RM) *.o
	$(RM) *.so
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

