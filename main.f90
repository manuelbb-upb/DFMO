program main_multiobj
    use, intrinsic :: iso_c_binding
	use main_mod
	implicit none

	real*8 		:: alfa_stop
	integer 	:: nf_max, iprint, hschoice
	logical		:: dir_dense, dir_coord

	!====================================================================================================
	!
	! Choice of the parameters that define DFMO algorithm. They are:
	!
	!	alfa_stop	: tolerance for step_length termination. DFMO will terminate as soon as all of the 
	!				  step_lengths fall below alfa_stop
	!	nf_max		: maximum number of allowed functions evaluations
	!	iprint		: printing level. 0 - no console output, >0 different levels of printing
	!	hschoice	: which type of dense direction is used. 1 - HALTON-type, 2 - SOBOL-type
	!	dir_dense	: whether to use the dense direction or not
	!	dir_coord	: whether to use the coordinate directions or not
	!
	!====================================================================================================
	alfa_stop	= 1.d-9
	nf_max		= 2000 !-500*n !20000
	iprint		= 0
	hschoice	= 2
	dir_dense	=.true.
	dir_coord	=.true.

	call opt_multiobj(alfa_stop, nf_max, iprint, hschoice, dir_dense, dir_coord)
end program main_multiobj