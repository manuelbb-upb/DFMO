FC = gfortran
RM = rm -f

FFLAGS = -O3
FFLAGS_SHARED = $(FFLAGS) -fPIC

EXE = multiobj
EXE_SHARED = $(EXE).so
OBJS_COMMON = modules_DFMO.o main.o DFMO.o subroutines_DFMO.o halton.o sobol.o qsortd.o
OBJS = problem_shim.o $(OBJS_COMMON)
OBJS_SHARED = $(OBJS_COMMON)

all: $(EXE)

$(EXE): $(OBJS)
	$(FC) -o $@ $^

%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

shared: $(EXE_SHARED)

$(EXE_SHARED): problem_standalone.o $(OBJS_SHARED)
	$(FC) -shared -o $@ $^

problem_standalone.o: problem_standalone.f03
	$(FC) -c $(FFLAGS_SHARED) $< -o $@

$(OBJS_SHARED): %.o: %.f90
	$(FC) -c $(FFLAGS_SHARED) $< -o $@
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

