subroutine opt_multiobj()

	use eps_mod
	use filtro
	use cache
	use vincoli
	use alfa_mod
	implicit none

	integer :: n,p,q
	integer			 :: n_int

	integer	:: icontr,numvar,j,istop,nobjn,icheck
	integer, allocatable	:: index_int(:)
	real		:: tbegin, tend
	character*30 :: nomefun
	integer ::i,k

	real*8, allocatable :: x(:),bl(:),bu(:),step(:)    
	real*8, allocatable :: punti(:,:)
	real*8, allocatable :: ciq(:)
	real*8, allocatable :: fob(:), finiz(:), alfaciniz(:),f(:)
	real*8, allocatable :: fob_media(:), ciq_media(:)

	integer ::            num_funct
	real*8             :: alfamax,delta 
	real*8			   :: fapp, alfainiz
	real*8			   :: violiniz, violint, fex
	real*8             :: alfa_stop
	integer            :: nf_max,iprint,hschoice
	logical			   :: flag, dir_dense, dir_coord
    integer 		   :: t

!------------------------------------------------------------------------------
  	  common /calfamax/alfamax
!------------------------------------------------------------------------------

	  num_funct = 0
	  call setdim(n,m,q)
	  p = 0
	  mm = m
	  write(*,*) 'n=',n,' q=',q,' m=',m,' mm=',mm,' p=',p
	  !pause
	  call setup_cache(n,q,m)
	  
	  allocate(Lkappa(ndim,n+q+2+n),Lnew(ndim,n+q+2+n),Ltilde(ndim,n+q+2+n))
	  Lkappa(:,n+q+1) = 0.d0
	  Lnew  (:,n+q+1) = 0.d0
	  Ltilde(:,n+q+1) = 0.d0
      t=n
      
	  allocate(index_int(n),alfaciniz(n))
	  allocate(x(n),bl(n),bu(n),step(n))
	  allocate(punti(n,t))
        
	  allocate(fob(q),finiz(q), f(q))
	  allocate(fob_media(q))
      allocate(fob_m(q),finiz_m(q))

      if(m.ge.1) then
	    allocate (ciq(m),epsiq(q,m),epsiq_in(q,m))
        allocate(eps(q,m),constr(m),epsiniz(q,m))
	    allocate(ciq_media(m))
        allocate(ciq_m(mm))
	  endif

	  allocate(alfaciniz_m(n))

	  fob = 1.d0
	  finiz = 0.d0

 	  call cpu_time(tbegin)

	  write(*,*) 'call setbounds...'
	  call setbounds(n,bl,bu)

	  write(*,*) 'bl = ',bl
	  write(*,*) 'bu = ',bu

	  alfainiz=1.0d1

	  alfainiz=min(10.d0,maxval(bu-bl)/10.d0)

	  do i = 1,n
        	alfaciniz(i) = min(10.d0,(bu(i)-bl(i))/10.d0)
	  enddo

        write(*,*) 'call startp...'
        call startp(n,x)

        call functs(n,x,q,fob)
        if(m.ge.1) call fconstriq(n,m,x,ciq)
		num_funct = num_funct+1

        write(*,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
        write(*,*) 'fob = ',fob
        write(*,*) 'ciq = ',ciq
        write(*,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

        do i=1,n

            if((x(i).lt.bl(i)).or.(x(i).gt.bu(i))) then
            write(*,*) ' starting point violates bound constraints'
            write(*,*) bl
            write(*,*) x
            write(*,*) bu
            !!!pause
            stop
            endif

        enddo

2002    format(2d20.10)

        !-----------------------------------------------------------------------
        !     print starting point info
        !-----------------------------------------------------------------------

        if(.not.find_in_cache(n,q,m,x,fob,ciq)) then   
            call functs(n,x,q,fob)
            if(m.ge.1) call fconstriq(n,m,x,ciq)
            call insert_in_cache(n,q,m,x,fob,ciq)
			num_funct = num_funct+1
        endif

        write(*,*) '--------------------------------------'
        write(*,*) '--------- EPSILON INIZIALI -----------'
        write(*,*) '--------------------------------------'
        do k = 1,q
            do i = 1,m
                if(max(0.d0,ciq(i)) < 1.d-0) then
	            	epsiq(k,i) = 1.d-3
                else
                	epsiq(k,i) = 1.d-1
                endif
            enddo
        enddo
             
    
    call startp(n,x)
	write(*,113) x
113	format(10(1x,es16.9))

	write(*,*) 'call diagonale...'
	call diagonale(n,t,bl,bu,punti)

	if(.not.find_in_cache(n,q,m,x,fob,ciq)) then   
        call functs(n,x,q,fob)
		if(m.ge.1) call fconstriq(n,m,x,ciq)
		call insert_in_cache(n,q,m,x,fob,ciq)
		num_funct = num_funct+1
	endif
	call functs_pen(n, q, m, x, finiz, fob, ciq)	
	call update_point(n,q,x,finiz,alfainiz, alfaciniz, flag)

	t = n
	do i = 1,t
		x = punti(:,i)
		write(*,113) x

		if(.not.find_in_cache(n,q,m,x,fob,ciq)) then   
			call functs(n,x,q,fob)
			if(m.ge.1) call fconstriq(n,m,x,ciq)
			call insert_in_cache(n,q,m,x,fob,ciq)
			num_funct = num_funct+1
		endif
		call functs_pen(n, q, m, x, finiz, fob, ciq)	
		if(intorno(n,q,x,finiz,1.d-2)) then
			call update_point(n,q,x,finiz,alfainiz, alfaciniz, flag)
		endif
	enddo

    write(*,*) ' ------------------------------------------------- '
    write(1,*) ' ------------------------------------------------- '

    write(*,*) '      f(xo) = ',fob
    write(1,*) '      f(xo) = ',fob

    viol=0.d0

    do i = 1,m
   	   viol=viol+max(0.d0,ciq(i))
    enddo

    write(*,*) '  cviol(xo) = ',viol
    write(1,*) '  cviol(xo) = ',viol



    write(*,*) ' ------------------------------------------------- '
    write(1,*) ' ------------------------------------------------- ' 		      
    do i=1,n
	   write(*,*) ' xo(',i,') =',x(i)
	   write(1,*) ' xo(',i,') =',x(i)
    enddo


    write(*,*) ' ------------------------------------------------- '
    write(1,*) ' ------------------------------------------------- '


!-----------------------------------------------------------------------
!	choice of starting penalty parameter values
!-----------------------------------------------------------------------

	if(m >= 1) epsiq_in     = epsiq

	finiz       = fob
	violiniz    = viol

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
   
    write(*,*)
    write(*,*) ' call the optimizer ...',n,q,m

    call sd_box(n,q,bl,bu,alfa_stop,nf_max,num_funct,iprint,hschoice,dir_dense,dir_coord,istop)

    write(*,*) ' ... done'
    write(*,*)

    call cpu_time(tend)

    write(*,*) ' ------------------------------------------------------------------------------'     
    write(1,*) ' ------------------------------------------------------------------------------'     
    if(istop.eq.1) then
       write(*,*)  ' STOP - stopping criterion satisfied. alfa_max <= ',alfa_stop
       write(1,*)  ' STOP - stopping criterion satisfied. alfa_max <= ',alfa_stop
    endif
    if(istop.eq.2) then
       write(*,*)  ' STOP - max number function evaluations exceeded. nf > ',nf_max
       write(1,*)  ' STOP - max number function evaluations exceeded. nf > ',nf_max
    endif

    if(istop.eq.3) then
       write(*,*)  ' STOP - max number iterations exceeded. ni > ',nf_max
       write(1,*)  ' STOP - max number iterations exceeded. ni > ',nf_max
    endif

    write(*,*) ' total time:',tend-tbegin
    write(1,*) ' total time:',tend-tbegin
    write(*,*) ' number of function evaluations = ',num_funct 
    write(1,*) ' number of function evaluations = ',num_funct     		    
    write(*,*) ' ------------------------------------------------------------------------------'  
    write(1,*) ' ------------------------------------------------------------------------------'  

	write(*,*)
	write(1,*)

    write(*,*) ' ------------------------------------------------- '
    write(1,*) ' ------------------------------------------------- '

    write(*,*) 
    write(1,*) 

133 format(a30,' & ',i6,' & ',es15.6,'\\\hline')

	  call destroy_cache()

	  deallocate(index_int,alfaciniz)
	  deallocate(x,bl,bu,step)
	  deallocate(punti)
	  deallocate(Lkappa,Lnew,Ltilde)
	  deallocate(fob,finiz, f, fob_media)
      deallocate(fob_m,finiz_m)
	  deallocate(alfaciniz_m)

      if(m.ge.1) then
	    deallocate (ciq,epsiq,epsiq_in,ciq_media)
        deallocate(eps,constr,epsiniz)
        deallocate(ciq_m)
      endif

end subroutine opt_multiobj

program main_multiobj
	call opt_multiobj()
end program main_multiobj