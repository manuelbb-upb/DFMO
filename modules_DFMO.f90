!============================================================================================
!    DFMO - Derivative-Free Linesearch program for constrained multiobjective 
!    nonsmooth optimization 
!    Copyright (C) 2016  G.Liuzzi, S.Lucidi, F.Rinaldi
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    G.Liuzzi, S.Lucidi, F.Rinaldi. A derivative-free approach to constrained multiobjective 
!    nonsmooth optimization, submitted for pubblication on SIAM J. on Optimization 
!============================================================================================
module alfa_mod
	integer			:: mm
	real*8, allocatable	:: fob_m(:), ciq_m(:), finiz_m(:)
	real*8, allocatable 	:: alfaciniz_m(:)
	real*8		    	:: alfainiz_m
end module alfa_mod

module eps_mod
     double precision, allocatable :: epseq(:,:), epsiq(:,:)
     double precision, allocatable :: epseq_in(:,:), epsiq_in(:,:)
end module eps_mod

module vincoli
     integer m
     double precision, allocatable :: eps(:,:),constr(:),epsiniz(:,:)
     double precision viol,violz
end module vincoli

module cache
	integer, parameter 	:: maxstore = 100000
	real*8, parameter	:: soglia = 1.d-6
	integer			:: ncache

	type tcache
		real*8, allocatable	:: x(:), fobj(:), fconstr(:)
	end type tcache
	type(tcache)		:: store(maxstore)
	integer			:: pos_free, cur_dim

	contains

	subroutine setup_cache(n,q,m)
		implicit none
		integer		:: n, q, m
		integer		:: i

		do i = 1,maxstore
			allocate(store(i)%x(n))
			allocate(store(i)%fobj(q))
			allocate(store(i)%fconstr(m))
		enddo
		pos_free = 1
		cur_dim = 0
 		ncache  = 0

		return

	end subroutine setup_cache
	subroutine destroy_cache()
		implicit none
		integer		:: i

		do i = 1,maxstore
			deallocate(store(i)%x)
			deallocate(store(i)%fobj)
			deallocate(store(i)%fconstr)
		enddo

		return

	end subroutine destroy_cache
	logical function find_in_cache(n,q,m,x,fob,fconstr)
		implicit none
		integer		:: n, q, m
		real*8		:: x(n), fob(q), fconstr(m)
		logical 	:: flag
		integer		:: i

		ncache = ncache+1
		flag = .false.
		!find_in_cache = .false.
		!return
		do i = 1,cur_dim
			if(maxval(abs(x-store(i)%x)) <= soglia) then
				!write(*,*) 'CACHE: x =',x
				!write(*,*) 'CACHE: xs=',store(i)%x
				
				flag = .true.
				fob  = store(i)%fobj
				fconstr = store(i)%fconstr
				exit
			endif			
		enddo

		find_in_cache = flag

		return
	end function find_in_cache
	subroutine insert_in_cache(n,q,m,x,fob,fconstr)
		implicit none
		integer		:: n, q, m
		real*8		:: x(n), fob(q), fconstr(m)
		logical 	:: flag
		integer		:: i

		!return
		if(.not.find_in_cache(n,q,m,x,fob,fconstr)) then
			store(pos_free)%x = x
			store(pos_free)%fobj = fob
			store(pos_free)%fconstr = fconstr
			pos_free = pos_free+1
			if(pos_free > maxstore) pos_free = 1
			if(cur_dim < maxstore) cur_dim = cur_dim+1
		else
			write(*,*) 'error in insert_in_cache: the passed point is already present'
			!pause
		endif

		return
	end subroutine insert_in_cache
end module cache

module filtro
	use problem_mod
	integer, parameter :: ndim = 100000
	!----------------------------------------------------------------
	! il filtro e' un array la cui generica riga di dimensione n+q+n+2
	! e' intesa come segue:
	!
	!   ------------------------------------------------------------------------------------
	!   | x(1)  ....   x(n) | f(1) .... f(q) | flag_busy | alfa | alfa_c(1) .... alfa_c(n) |
	!   ------------------------------------------------------------------------------------
	! 
	! ndim indica il numero di righe nella matrice
	! flag_busy = 100 se la riga corrispondente e' occupata da un non-dominato
	!           = 0   se la riga corrispondente e' vuota/libera
	!
	!----------------------------------------------------------------
	real*8, allocatable :: Lkappa(:,:), Lnew(:,:), Ltilde(:,:)

	contains 

	logical function intorno(n,q,x,f,delta)
		implicit none
		integer n,q,i
		real*8 x(n), f(q), delta
		real*8 theta

		theta = 1.d-0
		intorno = .true.
		do i = 1,ndim
			if(Lkappa(i,n+q+1) > 50.d0) then
				if(maxval(abs(Lkappa(i,n+1:n+q)-f)) < min(1.d-0,theta*delta)) then
					intorno = .false.
					exit
				endif
			endif
			if(Lnew(i,n+q+1) > 50.d0) then
				if(maxval(abs(Lnew(i,n+1:n+q)-f)) < min(1.d-0,theta*delta)) then
					intorno = .false.
					exit
				endif
			endif
			if(Ltilde(i,n+q+1) > 50.d0) then
				if(maxval(abs(Ltilde(i,n+1:n+q)-f)) < min(1.d-0,theta*delta)) then
					intorno = .false.
					exit
				endif
			endif
		enddo

		return
	end function intorno

	subroutine add_point(n,q,x,f,alfa,alfac)
		implicit none
		integer n,q,i,indfree
		real*8 x(n), f(q), alfa, alfac(n)
		integer, parameter :: iprint = 0
		indfree = -1
		do i = 1,ndim
			if(Lnew(i,n+q+1) < 50.d0) then
				indfree = i
				exit
			endif
		enddo

		if(indfree /= -1) then
			Lnew(indfree,1:n) = x
			Lnew(indfree,n+1:n+q) = f
			Lnew(indfree,n+q+1) = 100.d0
			Lnew(indfree,n+q+2) = alfa
			Lnew(indfree,n+q+2+1:n+q+2+n) = alfac
		else
			write(*,*) 'WARNING: Lnew pieno'
			!pause
		endif

		if(iprint > 0) then
			write(*,*) '---------- ADD_POINT START ----------'
			do i = 1,ndim
				if(Lnew(i,n+q+1) > 50.d0) then
					write(*,*) 'Lnw:',i,Lnew(i,n+1:n+q),Lnew(i,n+q+2:n+q+2+n)
				endif
			enddo
			write(*,*) '---------- ADD_POINT END ----------'
		endif

		return
	end subroutine add_point

	subroutine update_point(n,q,x,f,alfa,alfac,flag_change)
		implicit none
		integer n,q
		real*8 x(n), f(q), alfa, alfac(n)
		integer i, j, indfree
		integer flag
		logical flag_change, domina, minore_uguale
		!--------------------------------------------------------------
		! flag assume tre valori: 1,2,3
		! flag = 1 : Lkappa(i,:) < (x,f) 
		! flag = 2 : (x,f) < Lkappa(i,:)
		! flag = 3 : (x,f) non domina e non e' dominato Lkappa(i,:)
		!--------------------------------------------------------------
		flag_change = .false.
		indfree = -1
		flag = 0
		do i = 1,ndim
			flag = 0
			if((Lkappa(i,n+q+1) < 50.d0).and.(indfree == -1)) then
				indfree = i
			elseif(Lkappa(i,n+q+1) > 50.d0) then
			    if( (minore_uguale(q,Lkappa(i,n+1:n+q),f)) .and. &
		            	(minore_uguale(q,f,Lkappa(i,n+1:n+q))) ) then
				    flag = 1
				    exit
                	    endif                    
                	    if(domina(q,Lkappa(i,n+1:n+q),f)) then
					flag = 1
				elseif(domina(q,f,Lkappa(i,n+1:n+q))) then
					flag = 2
				else
					flag = 3
				endif
				!write(*,*) i,flag
				!flag = 0
				!do j = 1,q
				!	if(Lkappa(i,n+j) >= f(j)) then
				!		flag = .false.
				!		exit
				!	endif
				!enddomergese
				select case (flag)
				case(1)
					!------------------------------------------------
					!il punto passato e' dominato da uno nel filtro
					!------------------------------------------------
					exit
				case(2)
					!------------------------------------------------
					!il punto passato domina uno nel filtro
					!------------------------------------------------
					Lkappa(i,n+q+1) = 0.d0
					if(indfree == -1) indfree = i
				case(3)
				end select
			endif
		enddo
		if(flag /= 1) then
			flag_change = .true.
			if(indfree /= -1) then
				Lkappa(indfree,1:n) = x
				Lkappa(indfree,n+1:n+q) = f
				Lkappa(indfree,n+q+1) = 100.d0
				Lkappa(indfree,n+q+2) = alfa
				Lkappa(indfree,n+q+2+1:n+q+2+n) = alfac
			else
				write(*,*) 'WARNING: Filtro pieno'
				!pause
			endif
		endif
		return
	end subroutine update_point

	subroutine mergesets(n,q,flag_change)
		implicit none
		integer n,q,i
		logical flag_change, flag

		flag_change = .false.
		do i = 1,ndim
			if(Lnew(i,n+q+1) > 50.d0) then
				write(*,*) 'Lnew:',i,Lnew(i,n+1:n+q)
			endif
			if(Lkappa(i,n+q+1) > 50.d0) then
				write(*,*) 'Lkap:',i,Lkappa(i,n+1:n+q)
			endif
			if(Lnew(i,n+q+1) < 50.d0) cycle
			call update_point(n,q,Lnew(i,1:n),Lnew(i,n+1:n+q), Lnew(i,n+q+2), Lnew(i,n+q+2+1:n+q+2+n), flag)
			flag_change = flag_change.or.flag
		enddo
		!pause
		return
	end subroutine mergesets

	subroutine print_filter(n,q)
		use vincoli
		implicit none
		integer n,q,i,j
		real*8	x(n), ciq(m)
		integer totfeas

		totfeas = 0
		do i = 1,ndim
			if(Lkappa(i,n+q+1) > 50.d0) then
				write(*,100) i
				write(*,200)
				do j = 1,n
					write(*,110) Lkappa(i,j)
					write(33,110) Lkappa(i,j)
				enddo
				write(*,300)
				do j = 1,q
					write(*,110) Lkappa(i,n+j)
					write(33,110) Lkappa(i,n+j)
				enddo
				write(*,400)
				write(*,110) Lkappa(i,n+q+2)
				write(33,110) Lkappa(i,n+q+2)
				do j = 1,n
					write(*,110) Lkappa(i,n+q+2+j)
					write(33,110) Lkappa(i,n+q+2+j)
				enddo
				x = Lkappa(i,1:n)
				if(m.ge.1) then
					call fconstriq(n,m,x,ciq)
					viol = max(0.d0,maxval(ciq))
				else	
					viol = 0.d0
				endif
				if(viol < 1.d-3) totfeas = totfeas+1
				write(*,110) viol
				write(33,110) viol
				write(*,*)
				write(33,*)
			endif
		enddo
		write(*,*) 'totfeas = ',totfeas
		return
	100 format(1x,i7,$)
        200 format(1x,'x=',$)
        300 format(1x,'f=',$)
        400 format(1x,'alfa=',$)
	110 format(1x,es20.13,$)
	end subroutine print_filter
end module filtro
