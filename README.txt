-----------------------------------------------------------
 How to use the derivative-free optimizer DFMO 
-----------------------------------------------------------

1- Gunzip and untar the archive dfmo.tar.gz in a folder on your computer.

2- Edit file problem.f90 to define your own problem.
   In particular, you need to define the routines
   setdim		which sets the dimensions of the problem
				i.e. number of variables (n),
					 number of objectives (q),
					 number of inequality constraints (m)
				N.B. any equality constraint h(x) (if present) must be rewritten
				by the user as h(x) <= 0, -h(x) <= 0

   startp	 	which sets the starting point
   setbounds 	which sets upper and lower bounds on the variables
   functs 		which returns the obj. functions values
   fconstriq 	which returns the constraints values
				N.B. constraints are assumed to be in the form g(x) <= 0

3- Edit file main.f90 to define algorithm parameters.
   In particular, before the call to routine sd_box, the user can specify:
		alfa_stop	: tolerance for step_length termination. DFMO will terminate as soon as all of the 
					  step_lengths fall below alfa_stop
		nf_max		: maximum number of allowed functions evaluations
		iprint		: printing level. 0 - no console output, >0 different levels of printing
		hschoice	: which type of dense direction is used. 1 - HALTON-type, 2 - SOBOL-type
		dir_dense	: whether to use the dense direction or not
		dir_coord	: whether to use the coordinate directions or not
 
4- At command prompt execute 

     $> make
 
   which will create the executable 'multiobj'

5- execute

     $> ./multiobj

   The computed Pareto-front can be found in files pareto_fobs.out and pareto_vars.out

   pareto_vars.out : contains the coordinate of the non-dominated points (one per each line) 
				     along with the violation of the constraints
   pareto_fobs.out : contains the function values of the non-dominated points 
                    (one per each line) along with the violation of the constraints
   N.B. both the above files have one heading line 
