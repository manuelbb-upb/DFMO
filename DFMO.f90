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
subroutine sd_box(n,q,bl,bu,alfa_stop,nf_max,nf,iprint,hschoice,dir_dense,dir_coord,istop)
	use vincoli
	use cache
	use filtro
	implicit none
	integer :: n,q,i,j,k,ii,nf,nf_max, pn, pnum
	real*8  :: alfas(ndim), Lktemp(ndim,n+q+2+n), deltaf(ndim)
	integer :: ipt(ndim)
	logical :: dir_dense, dir_coord
	integer :: istop, ncomp
	integer :: iprint
	integer :: index_halton
	integer :: ibest, ielle, ifront
	logical :: flag_spread
	integer*8 :: index_sobol
	integer*8 ::n8

	integer :: hschoice

	real*8	:: ciq(m)

	real*8 :: x(n),z(n),d(n), d_dense(n), H(n,n), d_coord(n,n), alfainiz, alfaciniz(n)
	real*8 :: alfa_d(n),alfa,alfa_max, alfa_d_old, alfa_dense,Lmaxalfa,Lmaxalfa2
	real*8 :: f(q),fz(q) , eta, fob(q), fpar(q), gamma, coef_delta
	real*8 :: bl(n),bu(n),alfa_stop, fbest
	logical:: flag_dense, flag_change, trovato, condizione, aspetta
	logical:: maggiore_stretto, minore_uguale, Lkappa_domina
	real*8 :: dummy(n+q+2+n)
	real*8 :: violmin, violmax

	alfainiz=min(10.d0,maxval(bu-bl)/10.d0)
	do i = 1,n
		alfaciniz(i)=min(10.d0,(bu(i)-bl(i))/10.d0)
	enddo

	coef_delta = 1.d0
	eta = 1.d-6
	gamma=1.d-6    
	
    n8 = n
    
	istop = 0
	alfa_dense = 0.0d0

	index_halton = 1000
	index_sobol  = 10000

	d_coord = 0.d0
        do i=1, n
           d_coord(i,i)=1.0d0
   	end do


	index_halton=index_halton+2*n
	if(n>1) then
		if (hschoice.eq.1) then
			call halton(n,index_halton,d_dense)
		else 
			call i8_sobol(n8,index_sobol,d_dense)
		endif
	else
		d_dense = 1.d0
	endif

	call gen_base(n,d_dense,H)
	call gram_schmidt(n,H)
    
    !=================================================
    !                     INIZIO PROCEDURA
    !=================================================
    do 
        
	    do i = 1,q
		    fbest = 1.d+100
		    do pn=1,ndim
			    if(Lkappa(pn,n+q+1) > 50.d0) then
				    if(fbest > Lkappa(pn,n+i)) then
					    fbest = Lkappa(pn,n+i)
					    ibest = pn
				    endif
			    endif
		    enddo

		    dummy = Lkappa(ibest,:)
		    Lkappa(ibest,n+q+1) = 0.d0
		    do pn = ibest,2,-1
			    Lkappa(pn,:) = Lkappa(pn-1,:)
		    enddo
		    Lkappa(1,:) = dummy
	    enddo

	
	!==========================================
	!   ORDINAMENTO RISP. A SPREAD(GAMMA)
	!==========================================
       
	call spread_ord(n,q,deltaf,flag_spread)

	Lmaxalfa = 0.d0
	Lmaxalfa2 = 0.d0
	do pn=1,ndim
		if (Lkappa(pn,n+q+1)>50.d0) then
			if(Lmaxalfa < Lkappa(pn,n+q+2)) then
				Lmaxalfa = Lkappa(pn,n+q+2)
			endif
			if(Lmaxalfa2 < Lkappa(pn,n+q+2)) then
				Lmaxalfa2 = Lkappa(pn,n+q+2)
			endif
			do i = 1,n
				if(Lmaxalfa2 < Lkappa(pn,n+q+2+i)) then
					Lmaxalfa2 = Lkappa(pn,n+q+2+i)
				endif
			enddo
			
		endif
	enddo
	if(iprint > 0) then
		write(*,*) 'Lmaxalfa = ',Lmaxalfa
		write(*,*) 'Lmaxalfa2 = ',Lmaxalfa2
	endif

	!if (Lmaxalfa2<alfa_stop) exit

	aspetta=.false.
        trovato=.false.
        
	!====================================================================
	! Refine Lkappa (the filter)
	!====================================================================
	do pnum = 1,ndim
		pn = pnum
		
		if(Lkappa(pn,n+q+1)>50.d0) then		
			!pause
			viol = 0.d0
			if( m >= 1) then
				x = Lkappa(pn,1:n)
				if(.not.find_in_cache(n,q,m,x,fob,ciq)) then   
					!call functs(n,x,q,fob)
					call fconstriq(n,m,x,ciq)
					!call insert_in_cache(n,q,m,x,fob,ciq)
				endif
				viol = max(0.d0,maxval(ciq))
			
			endif
		endif		
		if(Lkappa(pn,n+q+1)>50.d0) then		
			alfa_max = 0.d0
			do i = 1,n
				if(alfa_max < Lkappa(pn,n+q+2+i)) then
					alfa_max = Lkappa(pn,n+q+2+i)
				endif
			enddo
			if(alfa_max < Lkappa(pn,n+q+2)) then
				alfa_max = Lkappa(pn,n+q+2)
			endif
		
		endif
	        if((deltaf(pn).eq.0.0d0).and.flag_spread) then
	             if(iprint > 0) write(*,*) pn
			endif

		condizione = (Lkappa(pn,n+q+1)>50.d0).and.   &
			(Lkappa(pn,n+q+2) > alfa_stop).and. &
 			((deltaf(pn) >= coef_delta).or..not.flag_spread).and. &
			((Lkappa(pn,n+q+2) > 1.d-1*Lmaxalfa).or.(pn<=q))

		if ( condizione ) then
		        
			if(iprint > 0) write(*,*)'ANALIZZA PUNTO ', pn
			!============================
			!   ASSEGNA PUNTO e ALFA
			!============================
			x=Lkappa(pn,1:n)
			f=Lkappa(pn,n+1:n+q)

			do i=1,n
				alfa_d(i)=Lkappa(pn,n+q+2+i)
			end do
			alfa_dense=Lkappa(pn,n+q+2)
			trovato=.true.

			if(iprint >= 1) then
				if(n<4) write(*,*) '         x= ', x
			            write(*,*) '         f= ', f
			            write(*,*) 'alfa_dense= ', alfa_dense
				write(*,*) 'alfa_coord= ', alfa_d
			endif

			z=x

			!================================================
			!               ANALISI SINGOLO PUNTO
			!================================================
			if(iprint.ge.1) then
				write(*,*) 'PROVA', alfa_max
				write(*,*) '----------------------------------------------'
				write(1,*) '----------------------------------------------'
				write(*,100) nf,alfa_max
				write(1,100) nf,alfa_max

100        format('nf=',i5,'   alfamax=',d12.5)

				write(1,*) 'f=',f
			endif
			if(iprint.ge.2) then
				do i=1,n
					write(*,*) ' x(',i,')=',x(i)
					write(1,*) ' x(',i,')=',x(i)
				enddo
			endif

			
			if(dir_coord.or.(n==1)) then
			do ii = 1,n
				d = d_coord(ii,:)
				if(iprint > 0) write(*,*) 'uso la ',ii,'direzione coord'

			   ifront = 0
			   DO IELLE = 1,2
				if(ifront == 1) exit
				!ifront = 0

				if(iprint > 0) write(*,*) 'ENTRO NELLA LINESEARCH COORD CON alfa=',alfa_d(ii),nf

				 if(d(ii).gt.0.d0) then

				     if((alfa_d(ii)-(bu(ii)-x(ii))).lt.(-1.d-6)) then  
		   			    alfa=dmax1(1.d-24,alfa_d(ii))
					 else
					    alfa=bu(ii)-x(ii)
						ifront=1
						if(iprint.ge.1) then
	                                                   write(*,*) 'a', x(ii),  bl(ii), bu(ii), alfa_d(ii), d(ii), alfa
							   write(*,*) ' punto espan. sulla front. *'
							   write(1,*) ' punto espan. sulla front. *'
						endif
					 endif

				  else

					 if((alfa_d(ii)-(x(ii)-bl(ii))).lt.(-1.d-6)) then
					    alfa=dmax1(1.d-24,alfa_d(ii))
					 else
						alfa=x(ii)-bl(ii)
						ifront=1
						if(iprint.ge.1) then
							   write(*,*) 'b', x(ii), bl(ii), bu(ii), alfa_d(ii), d(ii), alfa
							   write(*,*) ' punto espan. sulla front. *'
							   write(1,*) ' punto espan. sulla front. *'
						endif
					 endif

				  endif

	                          if(iprint > 0) write(*,*)'alfa_max==================>',alfa_max

				  if(dabs(alfa).le.1.d-3*dmin1(1.d0,alfa_max)) then
					 d(ii)=-d(ii)
					 ifront=0

					 if(iprint.ge.1) then
						   write(*,*) ' direzione opposta per alfa piccolo'
						   write(1,*) ' direzione opposta per alfa piccolo'
						   write(*,*) ' j =',ii,'    d(j) =',d(ii)
						   write(1,*) ' j =',ii,'    d(j) =',d(ii)
						   write(*,*) ' alfa=',alfa,'    alfamax=',alfa_max
						   write(1,*) ' alfa=',alfa,'    alfamax=',alfa_max
					  endif
					  if ((alfa.eq.0.0d0).and.(IELLE.eq.2) ) then
						Lkappa(pn,n+q+2+ii) = Lkappa(pn,n+q+2+ii)/2.d0
					  endif
					  alfa=0.d0
					  cycle
				  endif

				z(ii) = x(ii)+alfa*d(ii)

				if(.not.find_in_cache(n,q,m,z,fob,ciq)) then   
					call functs(n,z,q,fob)
					if(m>=1) call fconstriq(n,m,z,ciq)
					call insert_in_cache(n,q,m,z,fob,ciq)
					nf=nf+1
				else
					if(iprint > 0) write(*,*) 'trovato punto in cache'
				endif
				call functs_pen(n, q, m, z, fz, fob, ciq)
				violz=viol

				if(iprint.ge.1) then
						write(*,*) ' fz =',fz,'   alfa =',alfa
						write(1,*) ' fz =',fz,'   alfa =',alfa
				endif
				if(iprint.ge.2) then
					  do i=1,n
						  write(*,*) ' z(',i,')=',z(i)
						  write(1,*) ' z(',i,')=',z(i)
					  enddo
				endif

				fpar= f-gamma*alfa*alfa

				!if(maggiore_stretto(q,fz,fpar)) then
				if(Lkappa_domina(n,q,fz,gamma*alfa*alfa)) then 
	                                if(iprint > 0) write(*,*) 'FAILURE STEP'
					!-----------------------------
					! (Failure step)
					!-----------------------------
				
					IF(IELLE == 2) THEN
						Lkappa(pn,n+q+2+ii) = Lkappa(pn,n+q+2+ii)/2.d0	

					ELSE
						D(ii) = -D(ii)
					ENDIF
					z(ii) = x(ii)
				else 
					!-----------------------------
					! (Do a projected strong expansion)
					!-----------------------------
					call linesearchbox_cont(n,q,x,f,d,alfa,alfa_d,alfa_dense,z,fz,ii,&
						   alfa_max,iprint,bl,bu,nf, pn, ifront)
						  ! pause
					if(iprint > 0) write(*,*) 'ESCO  DALLA LINESEARCH COORD CON alfa=',alfa_d(ii),nf

		
					if(dabs(alfa).ge.1.d-12) then
						if(iprint > 0) then
							write(*,*) 'coord con successo ============'
							write(1,*) 'coord con successo ============'
						endif
					
					end if
				endif

			    ENDDO  ! DO IELLE

			    !pause
			enddo
			end if

			if(dir_dense.and.(n>1)) then
			do ii = 1,n
	
				d_dense = H(:,ii)
				if(iprint > 0) write(*,*) 'uso la',ii,' direzione densa'

			    DO IELLE = 1,2
				if(iprint > 0) write(*,*) 'ENTRO NELLA LINESEARCH DENSE CON alfa=',alfa_dense,nf
				!pause

				z = x+alfa_dense*d_dense
				z = max(bl,min(bu,z))

				 if(dsqrt(dot_product(z-x,z-x)).le.1.d-8) then
					d_dense = -d_dense
					!pause
					cycle
				 endif
				if(.not.find_in_cache(n,q,m,z,fob,ciq)) then   
					call functs(n,z,q,fob)
					if(m>=1) call fconstriq(n,m,z,ciq)
					call insert_in_cache(n,q,m,z,fob,ciq)
					nf=nf+1
				endif
				call functs_pen(n, q, m, z, fz, fob, ciq)				 
				viol = max(0.d0,maxval(ciq))

				fpar = f - gamma*alfa_dense*alfa_dense

				!if((maggiore_stretto(q,fz,fpar)).or.(dsqrt(dot_product(z-x,z-x)).le.1.d-8)) then
				if(Lkappa_domina(n,q,fz,gamma*alfa_dense*alfa_dense).or.(dsqrt(dot_product(z-x,z-x)).le.1.d-8)) then 
					if(iprint.ge.1) write(*,*) '   fz > fpar'
					!-----------------------------
					! (Failure step)
					!-----------------------------
					IF(IELLE == 2) THEN
						if(iprint.ge.1) write(*,*) '   ielle = 2, riduco alfa'
						Lkappa(pn,n+q+2) = alfa_dense/2.d0
					ELSE
						if(iprint.ge.1) write(*,*) '   ielle = 1, giro d'
						D_DENSE = -D_DENSE
					ENDIF
				else
					!-----------------------------
					! (Do a projected strong expansion)
					!-----------------------------
					call linesearchbox_dense(n,q,x,f,d_dense,alfa,alfa_dense,alfa_d,z,fz,&
						   alfa_max,iprint,bl,bu,nf, pn)
						  ! pause
					if(iprint > 0) write(*,*) 'ESCO  DALLA LINESEARCH DENSE CON alfa=',alfa,nf

					if(dabs(alfa).ge.1.d-12) then
						if(iprint > 0) then
							write(*,*) 'denso con successo ============'
							write(1,*) 'denso con successo ============'
						endif
				
					end if
				endif

			    ENDDO  ! DO IELLE

			enddo
			end if

			if(iprint > 0) then
				write(*,*) 'index_halton =',index_halton,alfa
				write(1,*) 'index_halton =',index_halton,alfa
			endif

			if(.not.find_in_cache(n,q,m,x,fob,ciq)) then   
				call functs(n,x,q,fob)
				if(m>=1) call fconstriq(n,m,x,ciq)
				call insert_in_cache(n,q,m,x,fob,ciq)
			!	nf=nf+1
			endif
			call functs_pen(n, q, m, x, f, fob, ciq)				 
			viol = max(0.d0,maxval(ciq))



			if(maxval(Lnew(:,n+q+1)) > 50.d0) then
				if(iprint > 0) then
					write(*,*) 'merge Lnew Ltilde'
					write(1,*) 'merge Lnew Ltilde'
				endif
				call merge(n,q,ndim,Lnew,Ltilde)
				Lnew(:,n+q+1) = 0.d0
			endif	

			if(Lkappa(pn,n+q+2) > alfa_stop) Lkappa(pn,n+q+2) = 5.d-1*Lkappa(pn,n+q+2)
			do i = 1,n
				if(Lkappa(pn,n+q+2+i) > alfa_stop) Lkappa(pn,n+q+2+i) = 5.d-1*Lkappa(pn,n+q+2+i)
			enddo

			call stopp(n,q,alfa_d,alfa_dense,istop,alfa_max,nf,f,alfa_stop,nf_max)
		         
			if (istop.ge.2) exit

		end if !if ((Lkappa(pn,n+q+1)>50.d0).and.(Lkappa(pn,n+q+2) > alfa_stop)) then
         
	!=======================================
	!          FINE ANALISI SINGOLO PUNTO
	!=======================================
        
	enddo   ! finita l'analisi di tutti gli elementi delle liste

	!call print_filter(n,q)


	!=======================================
	!              MERGE
	!=======================================
	if(iprint > 0) then
		write(*,*) 'merge Ltilde Lkappa'
		write(1,*) 'merge Ltilde Lkappa'
	endif
	call merge(n,q,ndim,Ltilde,Lkappa)

	!=======================================
	!SVUOTA Lnew e Ltilde
	!=======================================
	Lnew(:,n+q+1) = 0.d0
	Ltilde(:,n+q+1) = 0.d0

	index_halton=index_halton+2*n
	if(n>1) then
		if (hschoice.eq.1) then
			call halton(n,index_halton,d_dense)
		else 
			call i8_sobol(n8,index_sobol,d_dense)
		endif
	else
		d_dense = 1.d0
	endif
	call gen_base(n,d_dense,H)
	call gram_schmidt(n,H)

	if(istop >= 2) then
		if(iprint > 0) write(*,*) 'istop = ',istop
		exit
	endif

	!========================================
	!           STAMPA PUNTI
	!========================================
	if (.not.trovato) then
		coef_delta = 0.95d0*coef_delta
		if(iprint > 0) write(*,*) 'not trovato', coef_delta
		if(coef_delta < 1.d-16) exit
		!exit
		!pause
	endif

      end do


	  open(2,file='pareto_fobs.out',status='replace')
	  open(3,file='pareto_vars.out',status='replace')

	  write(3,*) '     x(1) ... x(n)  viol  '
      write(2,*) '     f(1) ... f(q)  viol  '

      j = 0
     
      do i=1,ndim

	   if(Lkappa(i,n+q+1)>50.0d0) then
             j = j+1
	     if(m.ge.1) then
		call fconstriq(n,m,Lkappa(i,1:n),ciq)
		viol = max(0.d0,maxval(ciq))
	     else
		viol = 0.d0
	     endif
	     
	     write(3,*) Lkappa(i,  1:n  ), viol
	     write(2,*) Lkappa(i,n+1:n+q), viol
	   end if
      end do
	  if(iprint > 0) then
		  write(*,700) j, nf, ncache
		  write(*,*) 'Nondom. = ',j

		  write(*,*) 'n.f.evals(inc. cache hits) = ',ncache,' n.f.evals(eff) = ',nf
	  endif
700   format(I7,' & ',I7,' & ',I7,' & ',I7,' \\\hline')

	  close(2)
	  close(3)
      return

end

!     #######################################################

subroutine stopp(n,q, alfa_d,alfa_dense,istop,alfa_max,nf,f,alfa_stop,nf_max)
	implicit none

	integer :: n,istop,i,nf,nf_max, q


	real*8 :: alfa_d(n),alfa_max,f(q),alfa_stop
	real*8 :: alfa_dense

	logical :: test

	istop=0

	alfa_max=0.0d0

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(n>1) alfa_max = max(alfa_max,alfa_dense)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if(alfa_max.le.alfa_stop) then
	test=.true.
	if (test.eqv..true.) then
	   istop = 1
	end if

	end if

	if(nf.gt.nf_max) then
	istop = 2
	end if

	return

end subroutine stopp
 
subroutine linesearchbox_cont(n,q,x,f,d,alfa,alfa_d,alfa_dense,z,fz,i_corr,&
                         alfa_max,iprint,bl,bu,nf,pn, ifront)
      
	  use vincoli
	  use filtro
	  use cache
      implicit none

	real*8	:: ciq(m)

      integer :: n,i_corr,nf,q
      integer :: i,j,k
      integer :: iprint, pn
      integer :: ifront
      real*8 :: x(n),d(n),alfa_d(n),z(n),bl(n),bu(n), alfac(n),zdelta(n)
      real*8 :: fob(q), f(q),alfa,alfa_dense,alfa_max,alfaex, fz(q),gamma
      real*8 ::  gamma_int, fpar2(q)
      real*8 :: delta,delta1,fpar(q),fzdelta(q),violzdelta
      logical :: minore_uguale, maggiore_stretto, domina, lkappa_domina
	  
	  gamma=1.d-6    
 	  
        delta =0.5d0
        delta1 =0.5d0

     

	alfac = alfa_d
	j=i_corr

	zdelta = x

	if(iprint.ge.1) then
			write(*,*) 'variabile continua  j =',j,'    d(j) =',d(j),' alfa=',alfa_d(j)
			write(1,*) 'variabile continua  j =',j,'    d(j) =',d(j),' alfa=',alfa_d(j)
	endif



          !================================================================
	  !               PROJECTED STRONG EXPANSION 
	  !================================================================

	do

		if((ifront.eq.1)) then

			if(iprint.ge.1) then
				 write(*,*) ' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa
				 write(1,*) ' accetta punto sulla frontiera fz =',fz,'   alfa =',alfa
			endif


			alfa_d(j) = alfa
			call add_point(n,q,z,fz,alfa_dense,alfa_d)

			return

		end if

		if(d(j).gt.0.d0) then
					
			 if((alfa/delta1-(bu(j)-x(j))).lt.(-1.d-6)) then
				 alfaex=alfa/delta1
			 else
				 alfaex=bu(j)-x(j)
				 ifront=1
				 if(iprint.ge.1) then
					write(*,*) ' punto espan. sulla front.'
					write(1,*) ' punto espan. sulla front.'
				 endif
			 end if

		else

			 if((alfa/delta1-(x(j)-bl(j))).lt.(-1.d-6)) then
				 alfaex=alfa/delta1
			 else
				 alfaex=x(j)-bl(j)
				 ifront=1
				 if(iprint.ge.1) then
					write(*,*) ' punto espan. sulla front.'
					write(1,*) ' punto espan. sulla front.'
				 endif
			 end if

		endif

		zdelta(j) = x(j)+alfaex*d(j) 

		if(.not.find_in_cache(n,q,m,zdelta,fob,ciq)) then   
			call functs(n,zdelta,q,fob)
			if(m>=1) call fconstriq(n,m,zdelta,ciq)
			call insert_in_cache(n,q,m,zdelta,fob,ciq)
			nf=nf+1
		endif
		call functs_pen(n, q, m, zdelta, fzdelta, fob, ciq)				 
		violzdelta=viol

		if(iprint.ge.1) then
			  write(*,*) 'n,q,m=',n,q,m
			  write(*,*) '  fob=',fob
			  write(*,*) '  ciq=',ciq
			  write(*,*) ' fzex=',fzdelta,'  alfaex=',alfaex  
			  write(1,*) ' fzex=',fzdelta,'  alfaex=',alfaex
		endif
		if(iprint.ge.2) then
			  do i=1,n
				 write(*,*) ' z(',i,')=',z(i),zdelta(i)
				 write(1,*) ' z(',i,')=',z(i),zdelta(i)
			  enddo
		endif

		fpar = f-gamma*alfaex*alfaex
		fpar2= fz-gamma*(alfaex-alfa)*(alfaex-alfa)
		fpar2= fz-gamma*alfaex*alfaex+gamma*alfa*alfa
		if( .not.minore_uguale(q,fzdelta,fpar2) ) then
			if(iprint.ge.1) then
				write(*,*) ' AGGIUNGE PUNTo Lnew fz =',fz,'   alfa =',alfa
				write(1,*) ' AGGIUNGE PUNTo Lnew fz =',fz,'   alfa =',alfa
			endif

			alfa_d(j) = alfa
			call add_point(n,q,z,fz,alfa_dense,alfa_d)
			
		endif


		
		!if(maggiore_stretto(q,fzdelta,fpar)) then
		if( Lkappa_domina(n,q,fzdelta,gamma*alfaex*alfaex) ) then
			!write(*,*) 'fzdelta > fpar'
			if(iprint > 0) write(*,*) 'Lkappa_domina'
				
			return

		else    

			if(iprint > 0) write(*,*) 'fzdelta NOT > fpar'
			
			if(ifront.ne.1) then
				fz=fzdelta
				z = zdelta
	     			violz=violzdelta
				alfa=alfaex
			else
				if(iprint.ge.1) then
					write(*,*) ' AGGIUNGE PUNTo Lnew fz =',fz,'   alfa =',alfa
					write(1,*) ' AGGIUNGE PUNTo Lnew fz =',fz,'   alfa =',alfa
				endif

				alfa_d(j) = alfaex
				call add_point(n,q,zdelta,fzdelta,alfa_dense,alfa_d)
			endif

		endif
	enddo

end subroutine linesearchbox_cont

subroutine linesearchbox_dense(n,q,x,f,d,alfa,alfa_dense,alfac,z,fz,&
                                 alfa_max,iprint,bl,bu,nf, pn)
	  use vincoli
	  use filtro
	  use cache
      implicit none

!---------------------------------------------------------------------      
! QUESTA SUBROUTINE IMPLEMENTA: PROJECTED STRONG-EXPANSION
!---------------------------------------------------------------------      

	real*8	:: ciq(m)

      integer :: n,i_corr,nf,q
      integer :: i,j, pn
      integer :: iprint
	  integer :: ifront,ielle
      real*8 :: x(n),d(n),z(n),z1(n),zz(n),bl(n),bu(n),alfac(n),zdelta(n)
      real*8 :: fob(q),f(q),alfa,alfa_max,alfaex, fz(q), fz1(q),gamma
      real*8 :: gamma_int, alfa_dense, alfa_front, fpar2(q)
      real*8 :: delta,delta1,fpar(q),fzdelta(q),violzdelta
      logical :: minore_uguale, maggiore_stretto, domina, Lkappa_domina
	  
	  gamma=1.d-6  
	  

      	  delta =0.5d0
          delta1 =0.5d0
          


	  ifront=0

	  if(iprint.ge.1) then
			write(*,*) 'direzione halton, alfa=',alfa_dense
			write(1,*) 'direzione halton, alfa=',alfa_dense
	  endif
     
 

	  alfa=alfa_dense
     	  alfaex = alfa


	  !================================================================
	  !               PROJECTED STRONG EXPANSION 
	  !================================================================


	  do

		 alfaex= alfa/delta1
	     

		 zdelta = x+alfaex*d 
		 
		 zz= max(bl,min(bu,x+alfa*d))  
		 zdelta = max(bl,min(bu,zdelta))

		 if(dsqrt(dot_product(zdelta-zz,zdelta-zz)).le.1.d-8) then
			d = -d
			alfa = alfa/2.d0
			!pause
			exit
		 endif
	

		if(.not.find_in_cache(n,q,m,zdelta,fob,ciq)) then   
			call functs(n,zdelta,q,fob)
			if(m>=1) call fconstriq(n,m,zdelta,ciq)
			call insert_in_cache(n,q,m,zdelta,fob,ciq)
			nf=nf+1
		endif

		call functs_pen(n, q, m, zdelta, fzdelta, fob, ciq)				 
		  violzdelta=viol
	
		 if(iprint.ge.1) then
			  write(*,*) ' fzex=',fzdelta,'  alfaex=',alfaex  
			  write(1,*) ' fzex=',fzdelta,'  alfaex=',alfaex
		 endif
		 if(iprint.ge.2) then
			  do i=1,n
				 write(*,*) ' z(',i,')=',z(i)
				 write(1,*) ' z(',i,')=',z(i)
			  enddo
		 endif

		 
		fpar = f -gamma*alfaex*alfaex
		fpar2= fz-gamma*(alfaex-alfa)*(alfaex-alfa)
		fpar2= fz-gamma*alfaex*alfaex+gamma*alfa*alfa 
		!if( (.not.minore_uguale(q,fzdelta,fpar2)).and..not.Lkappa_domina(n,q,fz,gamma*alfa*alfa)) then
		if( .not.minore_uguale(q,fzdelta,fpar2) ) then
		
			if(iprint.ge.1) then
				write(*,*) ' AGGIUNGE PUNTo Lnew fz =',fz,'   alfa =',alfa
				write(1,*) ' AGGIUNGE PUNTo Lnew fz =',fz,'   alfa =',alfa
			endif

			call add_point(n,q,z,fz,alfa,alfac)
			
		endif
		
		if(Lkappa_domina(n,q,fzdelta,gamma*alfaex*alfaex)) then
		!if(maggiore_stretto(q,fzdelta,fpar)) then
			!write(*,*) 'fzdelta > fpar'
			if(iprint > 0) write(*,*) 'Lkappa_domina'
				
			return

		else    

			if(iprint > 0) write(*,*) 'fzdelta NOT > fpar'

			fz=fzdelta
			z = zdelta
     			violz=violzdelta
			alfa=alfaex

		endif

	  enddo

end subroutine linesearchbox_dense

subroutine spread_ord(n,q,app,flag)
	use filtro
	implicit none
	integer	:: n, q, i, j
	real*8		:: alfas(ndim), fobj(ndim,q), delta(ndim), app(ndim)	
	real*8		:: dummy(n+q+2+n), dummyr
	real*8		:: Lktemp(ndim,n+q+2+n)
	integer	:: dummyi
	integer	:: ipt(ndim), jlast, indf(ndim), ncomp, iapp(ndim)
	logical	:: primo,flag

	flag = .true.
	delta = 0.d0
	ncomp = 0
	do i = 1,ndim
		if(Lkappa(i,n+q+1) > 50.d0) then
			ncomp = ncomp+1
			fobj(ncomp,:) = Lkappa(i,n+1:n+q)
			indf(ncomp) = i
		endif
	enddo
	if(ncomp.le.1) then
		flag = .false.
		return
	endif
	
	do j = 1,q

		alfas = fobj(:,j)
		call qsortd(alfas(1:ncomp),ipt(1:ncomp),ncomp)
		do i=1, ncomp
			app(ncomp+1-i)=alfas(ipt(i))
			iapp(ncomp+1-i)=indf(ipt(i))
		end do
		alfas(1:ncomp) = app(1:ncomp)
		indf(1:ncomp) = iapp(1:ncomp)
		
		delta(indf(1)) = delta(indf(1))  + 1.d+10
		delta(indf(ncomp)) = delta(indf(ncomp)) + 1.d+10
		do i = 2,ncomp-1
			delta(indf(i)) = delta(indf(i)) + abs(alfas(i-1)-alfas(i+1))/(alfas(1)-alfas(ncomp))
		enddo
	enddo

	delta = delta/dble(ncomp)
	do i = 1,ndim
		if(Lkappa(i,n+q+1) < 50.d0) then
			delta(i) = -100.d0
		endif
	enddo

	call qsortd(delta,ipt,ndim)
	do i=1,ndim
		app(ndim-i+1) = delta(ipt(i))
		Lktemp(ndim-i+1,:) = Lkappa(ipt(i),:)
	end do
	Lkappa = Lktemp
	
	return
end subroutine spread_ord
