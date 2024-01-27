subroutine halton(n,index,xx)
	integer		:: n, index
	real*8		:: x(n), xx(n), f
	integer		:: base(n), j
	real*8, parameter :: pi = 3.1415926535897932384626433832795
	real*8		:: prodsin

	integer, save	:: primes(1000)

	include 'primes.inc'

	if(n==1) then
		xx(1) = 1.d0
		return
	endif

	do i = 1,n
		base(i) = primes(i)
		x(i)    = 0.d0
		f       = 1.d0 / dble(base(i))
		j       = index
		do while (j > 0)
			x(i) = x(i) + f * mod(j,base(i))
			!write(*,*) 'j = ',j
			j    = j / base(i)
			!write(*,*) 'j = ',j
			f    = f / dble(base(i))
		enddo
		!write(*,*) 'x(',i,')=',x(i)
	enddo
	xx = x
	return

	do i = 1,n-1
		base(i) = primes(i)
		x(i)    = 0.d0
		f       = 1.d0 / dble(base(i))
		j       = index
		do while (j > 0)
			x(i) = x(i) + f * mod(j,base(i))
			!write(*,*) 'j = ',j
			j    = j / base(i)
			!write(*,*) 'j = ',j
			f    = f / dble(base(i))
		enddo
		!write(*,*) 'x(',i,')=',x(i)
	enddo

	x(1:n-2) = x(1:n-2)*pi
	x(  n-1) = x(  n-1)*2.d0*pi

	xx = 1.d0
	xx(1) = cos(x(1))

	prodsin = sin(x(1))
	do i = 2,n-1
		xx(i) = prodsin*cos(x(i))
		prodsin = prodsin*sin(x(i))
	enddo
	xx(n) = prodsin
	!open(11,file='denso.txt',status='unknown',access='append')
	!write(11,100) xx
	!close(11)
100	format(10(1x,es20.13))
	
!   BEGIN
!       result = 0;
!       f = 1 / base;
!       i = index;
!       WHILE (i > 0) 
!       BEGIN
!           result = result + f * (i % base);
!           i = FLOOR(i / base);
!           f = f / base;
!       END
!       RETURN result;
END subroutine halton

subroutine gram_schmidt(n,H)
	implicit none
	integer			:: n
	real*8			:: proj(n)
	real*8			:: H(n,n)
	integer			:: i,j

	do i = 2,n
		proj = 0.d0
		do j = 1,i-1
			proj = proj + dot_product(H(:,i),H(:,j))/dot_product(H(:,j),H(:,j))*H(:,j)
		enddo
		H(:,i) = H(:,i) - proj
	enddo

	do i = 1,n
		H(:,i) = H(:,i)/dsqrt(dot_product(H(:,i),H(:,i)))
	enddo

	do i = 1,n
		do j = 1,n
			if (abs(H(j,i)) < 1.d-4) H(j,i) = 0.d0
		enddo
	enddo
	
	do i = 1,n
		H(:,i) = H(:,i)/dsqrt(dot_product(H(:,i),H(:,i)))
	enddo

end subroutine gram_schmidt

subroutine gen_base(n,d,H)
	implicit none
	integer			:: n
	real*8			:: d(n)
	real*8			:: H(n,n)
	integer			:: i,loc(1), ind, j
	
	loc = maxloc(abs(d))
	ind = loc(1)

	H = 0.d0
	H(:,1) = d

	do i = 2,ind
		H(i-1,i) = 1.d0
	enddo
	do i = ind+1,n
		H(i,i) = 1.d0
	enddo

	!write(*,*) 'ind = ',ind
	!do i = 1,n
	!	write(*,100) (H(i,j),j=1,n)
	!enddo

100 format(4(1x,es15.9))

	!pause

end subroutine gen_base
