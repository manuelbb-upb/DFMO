!###############################################################################
!#
!#   As described by Huband et al. in "A review of multiobjective test problems
!#   and a scalable test problem toolkit", IEEE Transactions on Evolutionary 
!#   Computing 10(5): 477-506, 2006.
!#           
!#   Example WFG7.
!#
!#   This file is part of a collection of problems developed for
!#   derivative-free multiobjective optimization in
!#   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
!#   Direct Multisearch for Multiobjective Optimization, 2010.
!#
!#   Written by the authors in June 1, 2010.
!#
!#   with constraints (6,3) from:
!#
!*     Test problems for NonSmooth InEquality Constrained minimization
!*
!*     Napsu Karmitsa (2003, inequality constrained version 2006-2007)      
!*
!*     Haarala M., Miettinen K. and Mäkelä M.M.: New Limited Memory
!*     Bundle Method for Large-Scale Nonsmooth Optimization, Optimization
!*     Methods and Software, Vol. 19, No. 6, 2004, 673-692.
!*
!*     Karmitsa N.: Test Problems for Large-Scale Nonsmooth Minimization,
!*     Reports of the Department of Mathematical Information Technology,
!*     Series B, Scientific Computing, B 4/2007, University of Jyväskylä, 
!*     Jyväskylä, 2007.
!
!###############################################################################
subroutine setdim(n,m,q)
	implicit none
	integer	:: n,m,q

	n = 8
	m = n-1
	q = 3

	return
end subroutine setdim

subroutine startp(n,x)
	implicit none
	integer	:: n
	real*8		:: x(n), l(n), u(n)

	!call setbounds(n,l,u)

	x = 0.d0
	
	return
end subroutine startp

subroutine setbounds(n,lb,ub)
	implicit none
	integer	:: n, i
	real*8		:: lb(n), ub(n)
	real*8, parameter :: pi = 4.d0*atan(1.d0);

	lb = 0.d0
	do i = 1,n
		ub(i) = 2.d0*dble(i)
	enddo

	return
end subroutine setbounds

subroutine functs(n,z,M,f)
	implicit none
	integer	:: n, M, i, j
	real*8		:: z(n), f(M)
	integer, parameter :: k = 4
	integer, parameter :: l = 4
	real*8, parameter :: pi = 4.d0*atan(1.d0)
	real*8, parameter :: pi2 = 2.d0*atan(1.d0)
	real*8		:: S(3), zmax(8), A(2), y(8), h(3), rsum(4)
	real*8		:: t1(8), t2(8), t3(3), w(8), t4(3), x(3)
	real*8		:: gsum1, gsum2
	real*8, parameter :: AA = 0.98d0/49.98d0;
	real*8, parameter :: BB = 0.02d0;
	real*8, parameter :: CC = 50.d0;
	real*8, parameter :: AAA = 0.02d0;
	real*8, parameter :: alpha = 1.d0
	real*8, parameter :: AAAA = 5.d0

	do i = 1,M
		S(i) = 2.d0*dble(i)
	enddo
	do i = 1,n
		zmax(i) = 2.d0*dble(i)
	enddo
	do i = 1,M-1
		A(i) = 1.d0
	enddo

	do i = 1,n
		y(i) = z(i)/zmax(i)
	enddo

	do i = 1,n
		w(i) = 1.d0
	enddo
	do i = 1,k
		gsum1 = 0.d0
		gsum2 = 0.d0
		do j = i+1,n
			gsum1 = gsum1 + w(j)*y(j)
			gsum2 = gsum2 + w(j)
		enddo
		rsum(i) = gsum1 / gsum2
	enddo
	do i = 1,k
		t1(i) = y(i)**(BB+(CC-BB)*(AA-(1.d0-2.d0*rsum(i))*abs(floor(0.5d0-rsum(i))+AA)))
	enddo
	do i = k+1,n
		t1(i) = y(i)
	enddo

	do i = 1,k
		t2(i) = t1(i)
	enddo
	do i = k+1,n
		t2(i) = abs(t1(i)-0.35d0)/abs(floor(0.35d0-t1(i))+0.35d0)
	enddo

	do i = 1,M-1
		gsum1 = 0.d0
		gsum2 = 0.d0
		do j = ((i-1)*k/(M-1)+1),(i*k/(M-1))
			gsum1 = gsum1 + w(j)*t2(j)
			gsum2 = gsum2 + w(j)
		enddo
		t3(i) = gsum1 / gsum2
	enddo
	gsum1 = 0.d0
	gsum2 = 0.d0
	do j = k+1,n
		gsum1 = gsum1 + w(j)*t2(j)
		gsum2 = gsum2 + w(j)
	enddo
	t3(M) = gsum1 / gsum2

	do i = 1,M-1
		x(i) = max(t3(M),A(i))*(t3(i)-0.5d0)+0.5d0
	enddo
	x(M) = t3(M)

	h(1) = 1.d0
	do i = 1,M-1
		h(1) = h(1)*sin(x(i)*pi2)
	enddo
	do j = 2,M-1
		h(j) = 1.d0
		do i = 1,M-j
			h(j) = h(j)*sin(x(i)*pi2)
		enddo
		h(j) = h(j)*cos(x(M-j+1)*pi2)
	enddo
	h(M) = cos(x(1)*pi2)

	do i = 1,M
		f(i) = x(M)+S(i)*h(i)
	enddo

	return
end subroutine functs

subroutine fconstriq(n,m,x,ciq)
	implicit none
	integer n,m,i
	real*8 x(n), ciq(m)

	do i = 1,m
		ciq(i) = x(i)**2 + x(i+1)**2 + x(i)*x(i+1) - 1.0D+00
	enddo

	return
 
end subroutine fconstriq
