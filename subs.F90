! ================================================================
!                     SUBROUTINES
! ================================================================

! ================================================================
SUBROUTINE initialconditions(q,x,dx,pi,fx,N,traceval)

! Subroutine to set up initial tracer

IMPLICIT NONE
DOUBLE PRECISION, INTENT(in) :: dx, x(N), pi
INTEGER, INTENT(in) :: N, fx(N), traceval
DOUBLE PRECISION, INTENT(out) :: q(N)
INTEGER :: i, j
DOUBLE PRECISION :: r

	q(fx)=0.0

	! tracerval==1 sine function
	! tracerval==2 step function
	! tracerval==3 one cloud
	! tracerval==4 many clouds

	if (traceval==1) then

		do i=1, N

			q(i)=0.5*(1.0+sin(2.0*pi*x(i)))
                        !q(i)= q(i)+1

		enddo

	elseif (traceval==2) then

		do i=1, N

			if (x(i) .gt. 0.25 .and. x(i) .lt. 0.75) then

				q(i)=1.0

			endif

		enddo

	elseif (traceval==3) then

		q(floor(N/2.0))=1.0

	elseif (traceval==4) then

		q(floor(4.0*N/32.0))=1.0
		q(floor(5.0*N/32.0))=1.0
		q(floor(10.0*N/32.0))=1.0
		q(floor(14.0*N/32.0))=1.0
		q(floor(16.0*N/32.0))=1.0
		q(floor(18.0*N/32.0))=1.0
		q(floor(26.0*N/32.0))=1.0
		q(floor(29.0*N/32.0))=1.0

	endif


END SUBROUTINE initialconditions
! ================================================================

! ================================================================
SUBROUTINE gridsetup(x,dx,fx,f1x,f2x,f3x,f4x,g1x,g2x,g3x,g4x,N)

! Subroutine to set up the grid

IMPLICIT NONE
DOUBLE PRECISION, INTENT(in) :: dx
INTEGER, INTENT(in) :: N
DOUBLE PRECISION, INTENT(out) :: x(N)
INTEGER, INTENT(out) :: fx(N), f1x(N), f2x(N),  & 
		f3x(N), f4x(N), g1x(N), & 
		g2x(N), g3x(N), g4x(N)
INTEGER :: m, mm
	DO m = 1, N
		x(m) = m*dx
	ENDDO

	DO mm=1, N
		fx(mm)=mm
	enddo
		f1x=mod(fx,N)+1
		f2x=mod(f1x,N)+1
		f3x=mod(f2x,N)+1
		f4x=mod(f3x,N)+1
		g1x=mod(fx-1,N)
		g1x(1)=N
		g2x=mod(g1x-1,N)
		g2x(2)=N
		g3x=mod(g2x-1,N)
		g3x(3)=N
		g4x=mod(g3x-1,N)
		g4x(4)=N

END SUBROUTINE gridsetup
! ================================================================

! ===============================================================
SUBROUTINE dump1(u,ytitle,N)

! Subroutine to output to matlab

IMPLICIT NONE

INTEGER, INTENT(IN) :: N
INTEGER :: ii
REAL*8, INTENT(IN) :: u(N)
REAL :: ustar4(N)
CHARACTER*(*) :: ytitle

! Convert to single precision to reduce size of output and improve readability

	ustar4 = u

	open (unit=23, file='output.m', status='unknown')
	WRITE(23,*) ytitle ,'=['
	WRITE(23,*) ustar4(:),';...'
	WRITE(23,*) '];'

END SUBROUTINE dump1
! ===============================================================

! ===============================================================
SUBROUTINE dump2(u,ytitle,N)

! Subroutine to output to matlab

IMPLICIT NONE

INTEGER, INTENT(IN) :: N
INTEGER :: ii
REAL*8, INTENT(IN) :: u(N)
REAL :: ustar4(N)
CHARACTER*(*) :: ytitle

! Convert to single precision to reduce size of output and improve readability

	ustar4 = u

	open (unit=23, file='evals.m', status='unknown')
	WRITE(23,*) ytitle ,'=['
	WRITE(23,*) ustar4(:),';...'
	WRITE(23,*) '];'

END SUBROUTINE dump2
! ===============================================================


! =========================================================
SUBROUTINE unlimitedppm1d(flux,z,v,dt,dx,N)

IMPLICIT NONE
INTEGER, INTENT(in):: N
DOUBLE PRECISION, INTENT(in) :: dx, dt, z(N), v
DOUBLE PRECISION, INTENT(out) :: flux(N)
DOUBLE PRECISION :: zL(N), zR(N), z6(N), zD(N), &
		zLL(N), zRR(N), zp(N), &
		D2s, D2sL, D2sR, s, ones, minpart, D2sLIM
INTEGER :: i, ai, bi, ci, di, ei, fi, j

! flux to left, use u from left cell edge

	ones=1.0

	DO i=1,N
		ci=i+2
		ai=i+1
		bi=i-1
		if (bi==0) then
			bi=N
		endif
		if (ai==N+1) then
			ai=1
		endif
		if (ci==N+1) then
			ci=1
		endif
		if (ci==N+2) then
			ci=2
		endif
	
		zp(i)=(-z(ci)+7.0*z(ai)+7.0*z(i)-z(bi))/12.0

	ENDDO

	DO i=1,N
		bi=i-1
		if (bi==0) then
			bi=N
		endif

		zL(i) = zp(bi)
		zR(i) = zp(i)

	ENDDO

	DO i=1,N
		ai=i+1
		bi=i-1
		if (bi==0) then
			bi=N
		endif
		if (ai==N+1) then
			ai=1
		endif

		zD(i) = zR(i) - zL(i)
		z6(i) = 6.0*(z(i)-0.5*(zL(i)+zR(i)) )

	ENDDO

	DO i=1,N
		ai=i+1
		bi=i-1
		if (bi==0) then
			bi=N
		endif
		if (ai==N+1) then
			ai=1
		endif

		if (v .ge. 0.0) then
    
        		flux(i) = zR(bi) - 0.5*(abs(v)*dt/dx)*( zD(bi) - z6(bi)*(1.0-(2.0/3.0)*(abs(v)*dt/dx)) )
		else    
        		flux(i) = zL(i) + 0.5*(abs(v)*dt/dx)*( zD(i) + z6(i)*(1.0-(2.0/3.0)*(abs(v)*dt/dx)) )
		endif
    
	ENDDO

END SUBROUTINE unlimitedppm1d
! =========================================================

! =========================================================
SUBROUTINE UNLIMITEDPPM1D_TLM(flux, fluxd, z, zd0, v, vd, dt, dx, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: dx, dt, z(n), v
  DOUBLE PRECISION, INTENT(IN) :: zd0(n), vd
  DOUBLE PRECISION, INTENT(OUT) :: flux(n)
  DOUBLE PRECISION, INTENT(OUT) :: fluxd(n)
  DOUBLE PRECISION :: zl(n), zr(n), z6(n), zd(n), zll(n), zrr(n), zp(n)&
& , d2s, d2sl, d2sr, s, ones, minpart, d2slim
  DOUBLE PRECISION :: zld(n), zrd(n), z6d(n), zdd(n), zpd(n)
  INTEGER :: i, ai, bi, ci, di, ei, fi, j
  INTRINSIC ABS
  DOUBLE PRECISION :: abs1d
  DOUBLE PRECISION :: abs0d
  DOUBLE PRECISION :: abs3d
  DOUBLE PRECISION :: abs3
  DOUBLE PRECISION :: abs2
  DOUBLE PRECISION :: abs2d
  DOUBLE PRECISION :: abs1
  DOUBLE PRECISION :: abs0
! flux to left, use u from left cell edge
  ones = 1.0
  zpd = 0.D0
  DO i=1,n
    ci = i + 2
    ai = i + 1
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    IF (ai .EQ. n + 1) ai = 1
    IF (ci .EQ. n + 1) ci = 1
    IF (ci .EQ. n + 2) ci = 2
    zpd(i) = (7.0*zd0(ai)-zd0(ci)+7.0*zd0(i)-zd0(bi))/12.0
    zp(i) = (-z(ci)+7.0*z(ai)+7.0*z(i)-z(bi))/12.0
  END DO
  zld = 0.D0
  zrd = 0.D0
  DO i=1,n
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    zld(i) = zpd(bi)
    zl(i) = zp(bi)
    zrd(i) = zpd(i)
    zr(i) = zp(i)
  END DO
  z6d = 0.D0
  zdd = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    IF (ai .EQ. n + 1) ai = 1
    zdd(i) = zrd(i) - zld(i)
    zd(i) = zr(i) - zl(i)
    z6d(i) = 6.0*(zd0(i)-0.5*(zld(i)+zrd(i)))
    z6(i) = 6.0*(z(i)-0.5*(zl(i)+zr(i)))
  END DO
  fluxd = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    IF (ai .EQ. n + 1) ai = 1
    IF (v .GE. 0.0) THEN
      IF (v .GE. 0.) THEN
        abs0d = vd
        abs0 = v
      ELSE
        abs0d = -vd
        abs0 = -v
      END IF
      IF (v .GE. 0.) THEN
        abs2d = vd
        abs2 = v
      ELSE
        abs2d = -vd
        abs2 = -v
      END IF
      fluxd(i) = zrd(bi) - 0.5*(dt*abs0d*(zd(bi)-z6(bi)*(1.0-2.0/3.0*(&
&       abs2*dt/dx)))/dx+abs0*dt*(zdd(bi)-z6d(bi)*(1.0-2.0/3.0*(abs2*dt/&
&       dx))+z6(bi)*2.0*dt*abs2d/(3.0*dx))/dx)
      flux(i) = zr(bi) - 0.5*(abs0*dt/dx)*(zd(bi)-z6(bi)*(1.0-2.0/3.0*(&
&       abs2*dt/dx)))
    ELSE
      IF (v .GE. 0.) THEN
        abs1d = vd
        abs1 = v
      ELSE
        abs1d = -vd
        abs1 = -v
      END IF
      IF (v .GE. 0.) THEN
        abs3d = vd
        abs3 = v
      ELSE
        abs3d = -vd
        abs3 = -v
      END IF
      fluxd(i) = zld(i) + 0.5*(dt*abs1d*(zd(i)+z6(i)*(1.0-2.0/3.0*(abs3*&
&       dt/dx)))/dx+abs1*dt*(zdd(i)+z6d(i)*(1.0-2.0/3.0*(abs3*dt/dx))-z6&
&       (i)*2.0*dt*abs3d/(3.0*dx))/dx)
      flux(i) = zl(i) + 0.5*(abs1*dt/dx)*(zd(i)+z6(i)*(1.0-2.0/3.0*(abs3&
&       *dt/dx)))
    END IF
  END DO
END SUBROUTINE UNLIMITEDPPM1D_TLM
! =========================================================


! =========================================================
SUBROUTINE ppmlin1d(flux,zlin,u,dt,dx,N)

! subroutine for Lin scheme in 1D

IMPLICIT NONE
INTEGER, INTENT(in):: N
DOUBLE PRECISION, INTENT(in) :: dx, dt, zlin(N), u
DOUBLE PRECISION, INTENT(out) :: flux(N)
DOUBLE PRECISION :: deltaz(N), deltazmax(N), deltazmin(N), &
	minpart, deltazmono(N), zplus(N), zminus(N), &
	dl(N), dr(N), zminusnew(N), zplusnew(N), z6(N)
	
INTEGER :: i, ai, bi

! flux to left 

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif
    
    		deltaz(i)=0.25*(zlin(ai)-zlin(bi))
    
    		deltazmax(i) = max(zlin(bi),zlin(i),zlin(ai))-zlin(i)
    		deltazmin(i) = zlin(i)-min(zlin(bi),zlin(i),zlin(ai))

	ENDDO

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif
    
    		minpart=min(abs(deltaz(i)),deltazmin(i),deltazmax(i));
    
    		if (deltaz(i) .ge. 0.0) then
       			deltazmono(i)=abs(minpart)
    		else
        		deltazmono(i)=-abs(minpart)
    		endif
    
	ENDDO

! Approx of z at the cell edge

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif
    
    		zminus(i)=0.5*(zlin(bi)+zlin(i))+(1.0/3.0)*(deltazmono(bi)-deltazmono(i))
    
	ENDDO

! continuous at cell edges

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif
    
    		zplus(i)=zminus(ai)
    
	ENDDO

! limit cell edges

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif
    
    		dl(i)=min(abs(2.0*deltazmono(i)),abs(zminus(i)-zlin(i)))
    		dr(i)=min(abs(2.0*deltazmono(i)),abs(zplus(i)-zlin(i)))
    
    		if (2.0*deltazmono(i) .ge. 0.0) then
        		dl(i) = abs(dl(i))
        		dr(i) = abs(dr(i))
    		else
        		dl(i) = -abs(dl(i))
        		dr(i) = -abs(dr(i))
    		endif
    
    		zminusnew(i)=zlin(i) - dl(i)
    		zplusnew(i)=zlin(i) + dr(i)
    
	ENDDO

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif
    
    		zminus(i)=zminusnew(i)
		zplus(i)=zplusnew(i)  
    
    		!z6(i) = 6.0*(zlin(i) - 0.5*(zminus(i)+zplus(i)))
    		z6(i) = 3.0*(dl(i)-dr(i))
    
	ENDDO

! calculate flux based on velocity direction

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif
    
    	if (u .ge. 0.0) then
        	flux(i) = zplus(bi) - 0.5*(u*dt/dx)*( -zminus(bi)+zplus(bi) - z6(bi)*(1.0-(2.0/3.0)*(u*dt/dx)) )
    	else
        	flux(i) = zminus(i) + 0.5*(abs(u)*dt/dx)*( zplus(i)-zminus(i) + z6(i)*(1.0-(2.0/3.0)*(abs(u)*dt/dx)) )
    	endif
    
	ENDDO

	

END SUBROUTINE ppmlin1d
! =========================================================

! =========================================================
SUBROUTINE PPMLIN1D_TLM(flux, fluxd, zlin, zlind, u, ud, dt, dx, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: dx, dt, zlin(n), u
  DOUBLE PRECISION, INTENT(IN) :: zlind(n), ud
  DOUBLE PRECISION, INTENT(OUT) :: flux(n)
  DOUBLE PRECISION, INTENT(OUT) :: fluxd(n)
  DOUBLE PRECISION :: deltaz(n), deltazmax(n), deltazmin(n), minpart, &
& deltazmono(n), zplus(n), zminus(n), dl(n), dr(n), zminusnew(n), &
& zplusnew(n), z6(n)
  DOUBLE PRECISION :: deltazd(n), deltazmaxd(n), deltazmind(n), minpartd&
& , deltazmonod(n), zplusd(n), zminusd(n), dld(n), drd(n), zminusnewd(n)&
& , zplusnewd(n), z6d(n)
  INTEGER :: i, ai, bi
  INTRINSIC MAX
  INTRINSIC MIN
  INTRINSIC ABS
  DOUBLE PRECISION :: abs1d
  DOUBLE PRECISION :: min1
  DOUBLE PRECISION :: abs4d
  DOUBLE PRECISION :: min1d
  DOUBLE PRECISION :: x3
  DOUBLE PRECISION :: x2
  DOUBLE PRECISION :: x2d
  DOUBLE PRECISION :: x1
  DOUBLE PRECISION :: abs0d
  DOUBLE PRECISION :: max1d
  DOUBLE PRECISION :: abs3d
  DOUBLE PRECISION :: x1d
  DOUBLE PRECISION :: y2d
  DOUBLE PRECISION :: abs4
  DOUBLE PRECISION :: abs3
  DOUBLE PRECISION :: abs2
  DOUBLE PRECISION :: abs2d
  DOUBLE PRECISION :: abs1
  DOUBLE PRECISION :: abs0
  DOUBLE PRECISION :: max1
  DOUBLE PRECISION :: y2
  DOUBLE PRECISION :: x3d
  DOUBLE PRECISION :: y1
  DOUBLE PRECISION :: y1d
  deltazmaxd = 0.D0
  deltazmind = 0.D0
  deltazd = 0.D0
! flux to left 
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    deltazd(i) = 0.25*(zlind(ai)-zlind(bi))
    deltaz(i) = 0.25*(zlin(ai)-zlin(bi))
    IF (zlin(bi) .LT. zlin(i)) THEN
      IF (zlin(i) .LT. zlin(ai)) THEN
        max1d = zlind(ai)
        max1 = zlin(ai)
      ELSE
        max1d = zlind(i)
        max1 = zlin(i)
      END IF
    ELSE IF (zlin(bi) .LT. zlin(ai)) THEN
      max1d = zlind(ai)
      max1 = zlin(ai)
    ELSE
      max1d = zlind(bi)
      max1 = zlin(bi)
    END IF
    deltazmaxd(i) = max1d - zlind(i)
    deltazmax(i) = max1 - zlin(i)
    IF (zlin(bi) .GT. zlin(i)) THEN
      IF (zlin(i) .GT. zlin(ai)) THEN
        min1d = zlind(ai)
        min1 = zlin(ai)
      ELSE
        min1d = zlind(i)
        min1 = zlin(i)
      END IF
    ELSE IF (zlin(bi) .GT. zlin(ai)) THEN
      min1d = zlind(ai)
      min1 = zlin(ai)
    ELSE
      min1d = zlind(bi)
      min1 = zlin(bi)
    END IF
    deltazmind(i) = zlind(i) - min1d
    deltazmin(i) = zlin(i) - min1
  END DO
  deltazmonod = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (deltaz(i) .GE. 0.) THEN
      x1d = deltazd(i)
      x1 = deltaz(i)
    ELSE
      x1d = -deltazd(i)
      x1 = -deltaz(i)
    END IF
    IF (x1 .GT. deltazmin(i)) THEN
      IF (deltazmin(i) .GT. deltazmax(i)) THEN
        minpartd = deltazmaxd(i)
        minpart = deltazmax(i)
      ELSE
        minpartd = deltazmind(i)
        minpart = deltazmin(i)
      END IF
    ELSE IF (x1 .GT. deltazmax(i)) THEN
      minpartd = deltazmaxd(i)
      minpart = deltazmax(i)
    ELSE
      minpartd = x1d
      minpart = x1
    END IF
    IF (deltaz(i) .GE. 0.0) THEN
      IF (minpart .GE. 0.) THEN
        deltazmonod(i) = minpartd
        deltazmono(i) = minpart
      ELSE
        deltazmonod(i) = -minpartd
        deltazmono(i) = -minpart
      END IF
    ELSE
      IF (minpart .GE. 0.) THEN
        abs0d = minpartd
        abs0 = minpart
      ELSE
        abs0d = -minpartd
        abs0 = -minpart
      END IF
      deltazmonod(i) = -abs0d
      deltazmono(i) = -abs0
    END IF
  END DO
  zminusd = 0.D0
! Approx of z at the cell edge
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    zminusd(i) = 0.5*(zlind(bi)+zlind(i)) + (deltazmonod(bi)-deltazmonod&
&     (i))/3.0
    zminus(i) = 0.5*(zlin(bi)+zlin(i)) + 1.0/3.0*(deltazmono(bi)-&
&     deltazmono(i))
  END DO
  zplusd = 0.D0
! continuous at cell edges
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    zplusd(i) = zminusd(ai)
    zplus(i) = zminus(ai)
  END DO
  dld = 0.D0
  drd = 0.D0
  zplusnewd = 0.D0
  zminusnewd = 0.D0
! limit cell edges
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (2.0*deltazmono(i) .GE. 0.) THEN
      x2d = 2.0*deltazmonod(i)
      x2 = 2.0*deltazmono(i)
    ELSE
      x2d = -(2.0*deltazmonod(i))
      x2 = -(2.0*deltazmono(i))
    END IF
    IF (zminus(i) - zlin(i) .GE. 0.) THEN
      y1d = zminusd(i) - zlind(i)
      y1 = zminus(i) - zlin(i)
    ELSE
      y1d = -(zminusd(i)-zlind(i))
      y1 = -(zminus(i)-zlin(i))
    END IF
    IF (x2 .GT. y1) THEN
      dld(i) = y1d
      dl(i) = y1
    ELSE
      dld(i) = x2d
      dl(i) = x2
    END IF
    IF (2.0*deltazmono(i) .GE. 0.) THEN
      x3d = 2.0*deltazmonod(i)
      x3 = 2.0*deltazmono(i)
    ELSE
      x3d = -(2.0*deltazmonod(i))
      x3 = -(2.0*deltazmono(i))
    END IF
    IF (zplus(i) - zlin(i) .GE. 0.) THEN
      y2d = zplusd(i) - zlind(i)
      y2 = zplus(i) - zlin(i)
    ELSE
      y2d = -(zplusd(i)-zlind(i))
      y2 = -(zplus(i)-zlin(i))
    END IF
    IF (x3 .GT. y2) THEN
      drd(i) = y2d
      dr(i) = y2
    ELSE
      drd(i) = x3d
      dr(i) = x3
    END IF
    IF (2.0*deltazmono(i) .GE. 0.0) THEN
      IF (dl(i) .GE. 0.) THEN
        dl(i) = dl(i)
      ELSE
        dld(i) = -dld(i)
        dl(i) = -dl(i)
      END IF
      IF (dr(i) .GE. 0.) THEN
        dr(i) = dr(i)
      ELSE
        drd(i) = -drd(i)
        dr(i) = -dr(i)
      END IF
    ELSE
      IF (dl(i) .GE. 0.) THEN
        abs1d = dld(i)
        abs1 = dl(i)
      ELSE
        abs1d = -dld(i)
        abs1 = -dl(i)
      END IF
      dld(i) = -abs1d
      dl(i) = -abs1
      IF (dr(i) .GE. 0.) THEN
        abs2d = drd(i)
        abs2 = dr(i)
      ELSE
        abs2d = -drd(i)
        abs2 = -dr(i)
      END IF
      drd(i) = -abs2d
      dr(i) = -abs2
    END IF
    zminusnewd(i) = zlind(i) - dld(i)
    zminusnew(i) = zlin(i) - dl(i)
    zplusnewd(i) = zlind(i) + drd(i)
    zplusnew(i) = zlin(i) + dr(i)
  END DO
  z6d = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    zminusd(i) = zminusnewd(i)
    zminus(i) = zminusnew(i)
    zplusd(i) = zplusnewd(i)
    zplus(i) = zplusnew(i)
!z6(i) = 6.0*(zlin(i) - 0.5*(zminus(i)+zplus(i)))
    z6d(i) = 3.0*(dld(i)-drd(i))
    z6(i) = 3.0*(dl(i)-dr(i))
  END DO
  fluxd = 0.D0
! calculate flux based on velocity direction
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (u .GE. 0.0) THEN
      fluxd(i) = zplusd(bi) - 0.5*(dt*ud*(-zminus(bi)+zplus(bi)-z6(bi)*(&
&       1.0-2.0/3.0*(u*dt/dx)))/dx+u*dt*(zplusd(bi)-zminusd(bi)-z6d(bi)*&
&       (1.0-2.0/3.0*(u*dt/dx))+z6(bi)*2.0*dt*ud/(3.0*dx))/dx)
      flux(i) = zplus(bi) - 0.5*(u*dt/dx)*(-zminus(bi)+zplus(bi)-z6(bi)*&
&       (1.0-2.0/3.0*(u*dt/dx)))
    ELSE
      IF (u .GE. 0.) THEN
        abs3d = ud
        abs3 = u
      ELSE
        abs3d = -ud
        abs3 = -u
      END IF
      IF (u .GE. 0.) THEN
        abs4d = ud
        abs4 = u
      ELSE
        abs4d = -ud
        abs4 = -u
      END IF
      fluxd(i) = zminusd(i) + 0.5*(dt*abs3d*(zplus(i)-zminus(i)+z6(i)*(&
&       1.0-2.0/3.0*(abs4*dt/dx)))/dx+abs3*dt*(zplusd(i)-zminusd(i)+z6d(&
&       i)*(1.0-2.0/3.0*(abs4*dt/dx))-z6(i)*2.0*dt*abs4d/(3.0*dx))/dx)
      flux(i) = zminus(i) + 0.5*(abs3*dt/dx)*(zplus(i)-zminus(i)+z6(i)*(&
&       1.0-2.0/3.0*(abs4*dt/dx)))
    END IF
  END DO
END SUBROUTINE PPMLIN1D_TLM
! =========================================================

! =========================================================
SUBROUTINE PPMLIN1D0_TLM(flux, fluxd, zlin, zlind, u, ud, dt, dx, n)

  IMPLICIT NONE

  !Inputs
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: dx, dt, zlin(n), u
  DOUBLE PRECISION, INTENT(IN) :: zlind(n), ud

  !Outputs
  DOUBLE PRECISION, INTENT(OUT) :: flux(n)
  DOUBLE PRECISION, INTENT(OUT) :: fluxd(n)

  !Locals
  INTEGER :: i, ai, bi
  DOUBLE PRECISION :: deltaz (n), deltazmax (n), deltazmin (n), minpart
  DOUBLE PRECISION :: deltazd(n), deltazmaxd(n), deltazmind(n), minpartd
  DOUBLE PRECISION :: deltazmono (n), zplus (n), zminus (n)
  DOUBLE PRECISION :: deltazmonod(n), zplusd(n), zminusd(n)
  DOUBLE PRECISION :: dl (n), dr (n)
  DOUBLE PRECISION :: dld(n), drd(n)
  DOUBLE PRECISION :: zminusnew (n), zplusnew (n)
  DOUBLE PRECISION :: zminusnewd(n), zplusnewd(n)
  DOUBLE PRECISION :: z6 (n)
  DOUBLE PRECISION :: z6d(n)

   INTRINSIC MAX
  INTRINSIC MIN
  INTRINSIC ABS
  DOUBLE PRECISION :: abs1d
  DOUBLE PRECISION :: min1
  DOUBLE PRECISION :: abs4d
  DOUBLE PRECISION :: min1d
  DOUBLE PRECISION :: x3
  DOUBLE PRECISION :: x2
  DOUBLE PRECISION :: x2d
  DOUBLE PRECISION :: x1
  DOUBLE PRECISION :: abs0d
  DOUBLE PRECISION :: max1d
  DOUBLE PRECISION :: abs3d
  DOUBLE PRECISION :: x1d
  DOUBLE PRECISION :: y2d
  DOUBLE PRECISION :: abs4
  DOUBLE PRECISION :: abs3
  DOUBLE PRECISION :: abs2
  DOUBLE PRECISION :: abs2d
  DOUBLE PRECISION :: abs1
  DOUBLE PRECISION :: abs0
  DOUBLE PRECISION :: max1
  DOUBLE PRECISION :: y2
  DOUBLE PRECISION :: x3d
  DOUBLE PRECISION :: y1
  DOUBLE PRECISION :: y1d

! flux to left 

 DO i=1,N

    ai=i+1
    bi=i-1

    if (ai==N+1) then
       ai=1
    endif

    if (bi==0) then
       bi=N
    endif
    
    deltaz (i)=0.25*(zlin(ai)-zlin(bi))
    deltazd(i)=0.25*(zlind(ai)-zlind(bi))
    
    deltazmax(i) = max(zlin(bi),zlin(i),zlin(ai))-zlin(i)
    if     ( deltazmax(i) == zlin(bi)-zlin(i) ) then
       deltazmaxd(i) = zlind(bi)-zlind(i)
    elseif ( deltazmax(i) == zlin(i) -zlin(i) ) then
       deltazmaxd(i) = zlind(i) -zlind(i)
    elseif ( deltazmax(i) == zlin(ai)-zlin(i) ) then
       deltazmaxd(i) = zlind(ai)-zlind(i)
    endif

    deltazmin(i) = zlin(i)-min(zlin(bi),zlin(i),zlin(ai))
    if     ( deltazmin(i) == zlin(i)-zlin(bi) ) then
       deltazmind(i) = zlind(i)-zlind(bi)
    elseif ( deltazmin(i) == zlin(i)-zlin(i)  ) then
       deltazmind(i) = zlind(i)-zlind(i)
    elseif ( deltazmin(i) == zlin(i)-zlin(ai) ) then
       deltazmind(i) = zlind(i)-zlind(ai)
    endif

 ENDDO


 DO i=1,n

    ai=i+1
    bi=i-1

    if (ai==N+1) then
       ai=1
    endif

    if (bi==0) then
       bi=N
    endif

    minpart=min(abs(deltaz(i)),deltazmin(i),deltazmax(i))

    if ( minpart == abs(deltaz(i)) ) then
       if (deltaz(i) .ge. 0.0) then
          minpartd = deltazd(i)
       else
          minpartd = -deltazd(i)
       endif
    elseif ( minpart == deltazmin(i) ) then
       minpartd = deltazmind(i)
    elseif ( minpart == deltazmax(i) ) then
       minpartd = deltazmaxd(i)
    endif

    if (deltaz(i) .ge. 0.0) then
       if (minpart .ge. 0) then
          deltazmono (i)= minpart
          deltazmonod(i)= minpartd
       else
          deltazmono (i)=-minpart
          deltazmonod(i)=-minpartd
       endif
    else
       if (minpart .ge. 0) then
          deltazmono (i)=-minpart
          deltazmonod(i)=-minpartd
       else
          deltazmono (i)= minpart
          deltazmonod(i)= minpartd
       endif
    endif


 END DO

! Approx of z at the cell edge

 DO i=1,N

    ai=i+1
    bi=i-1

    if (ai==N+1) then
       ai=1
    endif

    if (bi==0) then
       bi=N
    endif
    
    zminus(i)  = 0.5*(zlin (bi)+zlin (i))+(1.0/3.0)*(deltazmono (bi)-deltazmono (i))
    zminusd(i) = 0.5*(zlind(bi)+zlind(i))+(1.0/3.0)*(deltazmonod(bi)-deltazmonod(i))

 ENDDO

! continuous at cell edges

 DO i=1,N
    ai=i+1
    bi=i-1

    if (ai==N+1) then
       ai=1
    endif

    if (bi==0) then
       bi=N
    endif
    
    zplus (i)=zminus (ai)
    zplusd(i)=zminusd(ai)
     
 ENDDO

! limit cell edges

  DO i=1,n

    ai = i+1
    bi = i-1

    IF (ai .EQ. n + 1) THEN
       ai = 1
    ENDIF

    IF (bi .EQ. 0) THEN
       bi = n
    ENDIF

    IF (2.0*deltazmono(i) .GE. 0.) THEN
      x2d = 2.0*deltazmonod(i)
      x2 = 2.0*deltazmono(i)
    ELSE
      x2d = -(2.0*deltazmonod(i))
      x2 = -(2.0*deltazmono(i))
    END IF

    IF (zminus(i) - zlin(i) .GE. 0.) THEN
      y1d = zminusd(i) - zlind(i)
      y1 = zminus(i) - zlin(i)
    ELSE
      y1d = -(zminusd(i)-zlind(i))
      y1 = -(zminus(i)-zlin(i))
    END IF

    IF (x2 .GT. y1) THEN
      dld(i) = y1d
      dl(i) = y1
    ELSE
      dld(i) = x2d
      dl(i) = x2
    END IF

    IF (2.0*deltazmono(i) .GE. 0.) THEN
      x3d = 2.0*deltazmonod(i)
      x3 = 2.0*deltazmono(i)
    ELSE
      x3d = -(2.0*deltazmonod(i))
      x3 = -(2.0*deltazmono(i))
    END IF
    IF (zplus(i) - zlin(i) .GE. 0.) THEN
      y2d = zplusd(i) - zlind(i)
      y2 = zplus(i) - zlin(i)
    ELSE
      y2d = -(zplusd(i)-zlind(i))
      y2 = -(zplus(i)-zlin(i))
    END IF
    IF (x3 .GT. y2) THEN
      drd(i) = y2d
      dr(i) = y2
    ELSE
      drd(i) = x3d
      dr(i) = x3
    END IF
    IF (2.0*deltazmono(i) .GE. 0.0) THEN
      IF (dl(i) .GE. 0.) THEN
        dl(i) = dl(i)
      ELSE
        dld(i) = -dld(i)
        dl(i) = -dl(i)
      END IF
      IF (dr(i) .GE. 0.) THEN
        dr(i) = dr(i)
      ELSE
        drd(i) = -drd(i)
        dr(i) = -dr(i)
      END IF
    ELSE
      IF (dl(i) .GE. 0.) THEN
        abs1d = dld(i)
        abs1 = dl(i)
      ELSE
        abs1d = -dld(i)
        abs1 = -dl(i)
      END IF
      dld(i) = -abs1d
      dl(i) = -abs1
      IF (dr(i) .GE. 0.) THEN
        abs2d = drd(i)
        abs2 = dr(i)
      ELSE
        abs2d = -drd(i)
        abs2 = -dr(i)
      END IF
      drd(i) = -abs2d
      dr(i) = -abs2
    END IF
    zminusnewd(i) = zlind(i) - dld(i)
    zminusnew(i) = zlin(i) - dl(i)
    zplusnewd(i) = zlind(i) + drd(i)
    zplusnew(i) = zlin(i) + dr(i)
  END DO


 DO i=1,N
    ai=i+1
    bi=i-1

    if (ai==N+1) then
       ai=1
    endif

    if (bi==0) then
       bi=N
    endif
    
    zminus (i)=zminusnew (i)
    zminusd(i)=zminusnewd(i)
    zplus (i)=zplusnew (i)  
    zplusd(i)=zplusnewd(i)  
    
    z6 (i) = 3.0*(dl (i)-dr (i))
    z6d(i) = 3.0*(dld(i)-drd(i))
    
 ENDDO

! calculate flux based on velocity direction

 DO i=1,N
    ai=i+1
    bi=i-1

    if (ai==N+1) then
       ai=1
    endif

    if (bi==0) then
       bi=N
    endif
    
    if (u .ge. 0.0) then
       flux (i) = zplus (bi) - 0.5*(u *dt/dx)*( -zminus (bi)+zplus (bi) - z6 (bi)*(1.0-(2.0/3.0)*(u *dt/dx)) )
       fluxd(i) = zplusd(bi) - 0.5*(u *dt/dx)*( -zminusd(bi)+zplusd(bi) - z6d(bi)*(1.0-(2.0/3.0)*(u *dt/dx))   &
                                                                        - z6 (bi)*(   -(2.0/3.0)*(ud*dt/dx)) ) &
                             - 0.5*(ud*dt/dx)*( -zminus (bi)+zplus (bi) - z6 (bi)*(1.0-(2.0/3.0)*(u *dt/dx)) )

    else
       flux (i) = zminus(i) + 0.5*(abs(u)*dt/dx)*( zplus(i)-zminus(i) + z6(i)*(1.0-(2.0/3.0)*(abs(u)*dt/dx)) )
       if (u .ge. 0) then
          fluxd(i) = zminusd(i) + 0.5*(u *dt/dx)*( zplusd(i)-zminusd(i) + z6d(i)*(1.0-(2.0/3.0)*(u *dt/dx))   &
                                                                        + z6 (i)*(1.0-(2.0/3.0)*(ud*dt/dx)) ) &
                                + 0.5*(ud*dt/dx)*( zplus (i)-zminus (i) + z6 (i)*(1.0-(2.0/3.0)*(u *dt/dx)) )
       else
          fluxd(i) = zminusd(i) + 0.5*(-u *dt/dx)*( zplusd(i)-zminusd(i) + z6d(i)*(1.0-(2.0/3.0)*(-u *dt/dx))   &
                                                                         + z6 (i)*(   -(2.0/3.0)*(-ud*dt/dx)) ) &
                                + 0.5*(-ud*dt/dx)*( zplus (i)-zminus (i) + z6 (i)*(1.0-(2.0/3.0)*(-u *dt/dx)) )
       endif
    endif
    
 ENDDO

END SUBROUTINE PPMLIN1D0_TLM
! =========================================================


! ===============================================================
SUBROUTINE getqhalf(qhalf,qin,order,N,C,fx,f1x,g1x,g2x,g3x)

IMPLICIT NONE
INTEGER, INTENT(in) :: N, order
INTEGER, INTENT(in) :: fx(N),f1x(N),g1x(N),g2x(N),g3x(N)
DOUBLE PRECISION, INTENT(in)    :: C, qin(N)
DOUBLE PRECISION, INTENT(inout) :: qhalf(N)


 if (order == 1) then

    qhalf(fx) = qin(g1x)

 elseif (order == 2) then

    qhalf(fx)=0.5*(qin(fx)+qin(g1x)) - 0.5*C*(qin(fx)-qin(g1x))

 elseif (order == 3) then

    qhalf(fx)=(1.0/6.0)*(2.0*qin(fx)+5.0*qin(g1x)-qin(g2x))-0.5*C*(qin(fx)-qin(g1x))+(C**2/6)*(qin(fx)-2.0* qin(g1x)+ qin(g2x))

 elseif (order == 4) then

    qhalf(fx)=(1.0/12.0)*(-qin(f1x)+7.0*qin(fx)+7.0*qin(g1x)-qin(g2x)) & 
			- 0.5*C*(-1.0/12.0*qin(f1x)+5.0/4.0*qin(fx)-5.0/4.0*qin(g1x)+1.0/12.0*qin(g2x)) & 
			+ (C**2.0/6.0)*(1.0/2.0*qin(f1x)-1.0/2.0*qin(fx)-1.0/2.0*qin(g1x)+1.0/2.0*qin(g2x)) & 
			- (C**3.0/24.0)*(qin(f1x)-3.0*qin(fx)+3.0*qin(g1x)-qin(g2x))

 elseif (order == 5) then

    qhalf(fx)=(1.0/60.0)*(-3.0*qin(f1x)+27.0*qin(fx)+47.0*qin(g1x)-13.0*qin(g2x)&
                                +2.0*qin(g3x)) - 0.5*C*(-1.0/12.0*qin(f1x)+5.0/4.0*qin(fx)&
                                -5.0/4.0*qin(g1x)+1.0/12.0*qin(g2x)) + (C**2.0/6.0)*(1.0/4.0*qin(f1x)&
                                +1.0/2.0*qin(fx)-2.0*qin(g1x)+3.0/2.0*qin(g2x)-1.0/4.0*qin(g3x)) & 
				- (C**3.0/24.0)*(qin(f1x)-3.0*qin(fx)+3.0*qin(g1x)-qin(g2x)) & 
				+ (C**4.0/120.0)*(qin(f1x)-4.0*qin(fx)+6.0*qin(g1x)-4.0*qin(g2x)+qin(g3x))

 else

   print*, 'ORDER OF ADVECTION NOT SUPPORTED'

 endif

END SUBROUTINE getqhalf
! ===============================================================

! ===============================================================
SUBROUTINE GETQHALF_TLM(qhalf, qhalfd, qin, qind, order, n, c, fx, f1x, &
& g1x, g2x, g3x)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, order
  INTEGER, INTENT(IN) :: fx(n), f1x(n), g1x(n), g2x(n), g3x(n)
  DOUBLE PRECISION, INTENT(IN) :: c, qin(n)
  DOUBLE PRECISION, INTENT(IN) :: qind(n)
  DOUBLE PRECISION, INTENT(INOUT) :: qhalf(n)
  DOUBLE PRECISION, INTENT(INOUT) :: qhalfd(n)
  IF (order .EQ. 1) THEN
    qhalfd = 0.D0
    qhalfd(fx) = qind(g1x)
    qhalf(fx) = qin(g1x)
  ELSE IF (order .EQ. 2) THEN
    qhalfd = 0.D0
    qhalfd(fx) = 0.5*(qind(fx)+qind(g1x)) - 0.5*c*(qind(fx)-qind(g1x))
    qhalf(fx) = 0.5*(qin(fx)+qin(g1x)) - 0.5*c*(qin(fx)-qin(g1x))
  ELSE IF (order .EQ. 3) THEN
    qhalfd = 0.D0
    qhalfd(fx) = (2.0*qind(fx)+5.0*qind(g1x)-qind(g2x))/6.0 - 0.5*c*(&
&     qind(fx)-qind(g1x)) + c**2*(qind(fx)-2.0*qind(g1x)+qind(g2x))/6
    qhalf(fx) = 1.0/6.0*(2.0*qin(fx)+5.0*qin(g1x)-qin(g2x)) - 0.5*c*(qin&
&     (fx)-qin(g1x)) + c**2/6*(qin(fx)-2.0*qin(g1x)+qin(g2x))
  ELSE IF (order .EQ. 4) THEN
    qhalfd = 0.D0
    qhalfd(fx) = (7.0*qind(fx)-qind(f1x)+7.0*qind(g1x)-qind(g2x))/12.0 -&
&     0.5*c*(5.0*qind(fx)/4.0-qind(f1x)/12.0-5.0*qind(g1x)/4.0+qind(g2x)&
&     /12.0) + c**2.0*(qind(f1x)/2.0-qind(fx)/2.0-qind(g1x)/2.0+qind(g2x&
&     )/2.0)/6.0 - c**3.0*(qind(f1x)-3.0*qind(fx)+3.0*qind(g1x)-qind(g2x&
&     ))/24.0
    qhalf(fx) = 1.0/12.0*(-qin(f1x)+7.0*qin(fx)+7.0*qin(g1x)-qin(g2x)) -&
&     0.5*c*(-(1.0/12.0*qin(f1x))+5.0/4.0*qin(fx)-5.0/4.0*qin(g1x)+1.0/&
&     12.0*qin(g2x)) + c**2.0/6.0*(1.0/2.0*qin(f1x)-1.0/2.0*qin(fx)-1.0/&
&     2.0*qin(g1x)+1.0/2.0*qin(g2x)) - c**3.0/24.0*(qin(f1x)-3.0*qin(fx)&
&     +3.0*qin(g1x)-qin(g2x))
  ELSE IF (order .EQ. 5) THEN
    qhalfd = 0.D0
    qhalfd(fx) = (27.0*qind(fx)-3.0*qind(f1x)+47.0*qind(g1x)-13.0*qind(&
&     g2x)+2.0*qind(g3x))/60.0 - 0.5*c*(5.0*qind(fx)/4.0-qind(f1x)/12.0-&
&     5.0*qind(g1x)/4.0+qind(g2x)/12.0) + c**2.0*(qind(f1x)/4.0+qind(fx)&
&     /2.0-2.0*qind(g1x)+3.0*qind(g2x)/2.0-qind(g3x)/4.0)/6.0 - c**3.0*(&
&     qind(f1x)-3.0*qind(fx)+3.0*qind(g1x)-qind(g2x))/24.0 + c**4.0*(&
&     qind(f1x)-4.0*qind(fx)+6.0*qind(g1x)-4.0*qind(g2x)+qind(g3x))/&
&     120.0
    qhalf(fx) = 1.0/60.0*(-(3.0*qin(f1x))+27.0*qin(fx)+47.0*qin(g1x)-&
&     13.0*qin(g2x)+2.0*qin(g3x)) - 0.5*c*(-(1.0/12.0*qin(f1x))+5.0/4.0*&
&     qin(fx)-5.0/4.0*qin(g1x)+1.0/12.0*qin(g2x)) + c**2.0/6.0*(1.0/4.0*&
&     qin(f1x)+1.0/2.0*qin(fx)-2.0*qin(g1x)+3.0/2.0*qin(g2x)-1.0/4.0*qin&
&     (g3x)) - c**3.0/24.0*(qin(f1x)-3.0*qin(fx)+3.0*qin(g1x)-qin(g2x)) &
&     + c**4.0/120.0*(qin(f1x)-4.0*qin(fx)+6.0*qin(g1x)-4.0*qin(g2x)+qin&
&     (g3x))
  ELSE
    PRINT*, 'ORDER OF ADVECTION NOT SUPPORTED'
    qhalfd = 0.D0
  END IF
END SUBROUTINE GETQHALF_TLM
! ===============================================================

! ===============================================================
SUBROUTINE relaxedlimiter(flux,q,v,dt,dx,N)

IMPLICIT NONE
INTEGER, INTENT(in):: N
DOUBLE PRECISION, INTENT(in) :: dx, dt, q(N), v
DOUBLE PRECISION, INTENT(inout) :: flux(N)
DOUBLE PRECISION :: qLinmin(N), qLinmax(N), qLp(N), qLpp(N), &
	qnewmax(N), qnewmin(N), qoutmin(N), qoutmax(N), qRppp(N), qR(N), &
	lowb(N), upb(N), zeros, maxq, minq
INTEGER :: i, ai, bi, ci, di, fi

	zeros=0.0

	DO i=1,N
		fi=i-3
		di=i-2
		ci=i+2
		ai=i+1
		bi=i-1
		if (bi==0) then
			bi=N
		endif
		if (ai==N+1) then
			ai=1
		endif
		if (ci==N+1) then
			ci=1
		endif
		if (ci==N+2) then
			ci=2
		endif
		if (di==0) then
			di=N
		endif
		if (di==-1) then
			di=N-1
		endif
		if (fi==0) then
			fi=N
		endif
		if (fi==-1) then
			fi=N-1
		endif
		if (fi==-2) then
			fi=N-2
		endif
                
                maxq = max(q(bi),q(i),q(ai))
                minq = min(q(bi),q(i),q(ai))
                
                if ( (flux(i)-q(bi))*(q(i)-flux(i)) .lt. 0.0) then 
                    
                    if ( (q(bi)-q(di))*(q(ai)-q(i)) .ge. 0.0 &
                            .or. (q(bi)-q(di))*(q(di)-q(fi)) .le. 0.0 &
                            .or. (q(ai)-q(i))*(q(ci)-q(ai)) .le. 0.0 &
                            .or. (flux(i)-q(bi))*(q(bi)-q(di)) .le. 0.0 ) then
                        
                        upb(i) = 0.0
                        lowb(i) = 0.0
                        
                    else
                        
                        upb(i) = (maxq-minq)
                        lowb(i) = (maxq-minq)
                        
                    endif
                    
                else
                    
                    upb(i) = (maxq-minq)
                    lowb(i) = (maxq-minq)
                    
                endif
	
	ENDDO    

	DO i=1,N
		ai=i+1
		bi=i-1
		if (bi==0) then
			bi=N
		endif
		if (ai==N+1) then
			ai=1
		endif

		qLinmin(i) = max(min(q(bi),q(i))-lowb(i),zeros)
		qLinmax(i) = max(q(bi),q(i))+upb(i) 
	ENDDO
    


	DO i=1,N
		ai=i+1
		bi=i-1
		if (bi==0) then
			bi=N
		endif
		if (ai==N+1) then
			ai=1
		endif
		qLp(i) = min(flux(i),qLinmax(i)) 
		qLpp(i) = max(qLp(i),qLinmin(i)) 
	ENDDO
    
	qnewmax(1:N) = qLinmax(1:N) 
	qnewmin(1:N) = qLinmin(1:N)

	DO i=1,N
		ai=i+1
		bi=i-1
		if (bi==0) then
			bi=N
		endif
		if (ai==N+1) then
			ai=1
		endif
		qoutmin(i) = (q(i) + (v*dt/dx)*qLinmax(i) - qnewmax(i))/(v*dt/dx) 
		qoutmax(i) = (q(i) + (v*dt/dx)*qLinmin(i) - qnewmin(i))/(v*dt/dx)
	ENDDO
    


	DO i=1,N
		ai=i+1
		bi=i-1
		if (bi==0) then
			bi=N
		endif
		if (ai==N+1) then
			ai=1
		endif
		qRppp(i) = min(qLpp(ai),qoutmax(i)) 
		qR(i) = max(qRppp(i),qoutmin(i)) 
	ENDDO

	DO i=1,N
		ai=i+1
		bi=i-1
		if (bi==0) then
			bi=N
		endif
		if (ai==N+1) then
			ai=1
		endif
		flux(i) = qR(bi)
	ENDDO
    

END SUBROUTINE relaxedlimiter
! =========================================================

! =========================================================
SUBROUTINE RELAXEDLIMITER_TLM(flux, fluxd, q, qd, v, vd, dt, dx, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: dx, dt, q(n), v
  DOUBLE PRECISION, INTENT(IN) :: qd(n), vd
  DOUBLE PRECISION, INTENT(INOUT) :: flux(n)
  DOUBLE PRECISION, INTENT(INOUT) :: fluxd(n)
  DOUBLE PRECISION :: qlinmin(n), qlinmax(n), qlp(n), qlpp(n), qnewmax(n&
& ), qnewmin(n), qoutmin(n), qoutmax(n), qrppp(n), qr(n), lowb(n), upb(n&
& ), zeros, maxq, minq
  DOUBLE PRECISION :: qlinmind(n), qlinmaxd(n), qlpd(n), qlppd(n), &
& qnewmaxd(n), qnewmind(n), qoutmind(n), qoutmaxd(n), qrpppd(n), qrd(n)&
& , lowbd(n), upbd(n), maxqd, minqd
  INTEGER :: i, ai, bi, ci, di, fi
  INTRINSIC MAX
  INTRINSIC MIN
  DOUBLE PRECISION :: min1
  DOUBLE PRECISION :: min1d
  DOUBLE PRECISION :: x1
  DOUBLE PRECISION :: max1d
  DOUBLE PRECISION :: x1d
  DOUBLE PRECISION :: max1
  zeros = 0.0
  upbd = 0.D0
  lowbd = 0.D0
  DO i=1,n
    fi = i - 3
    di = i - 2
    ci = i + 2
    ai = i + 1
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    IF (ai .EQ. n + 1) ai = 1
    IF (ci .EQ. n + 1) ci = 1
    IF (ci .EQ. n + 2) ci = 2
    IF (di .EQ. 0) di = n
    IF (di .EQ. -1) di = n - 1
    IF (fi .EQ. 0) fi = n
    IF (fi .EQ. -1) fi = n - 1
    IF (fi .EQ. -2) fi = n - 2
    IF (q(bi) .LT. q(i)) THEN
      IF (q(i) .LT. q(ai)) THEN
        maxqd = qd(ai)
        maxq = q(ai)
      ELSE
        maxqd = qd(i)
        maxq = q(i)
      END IF
    ELSE IF (q(bi) .LT. q(ai)) THEN
      maxqd = qd(ai)
      maxq = q(ai)
    ELSE
      maxqd = qd(bi)
      maxq = q(bi)
    END IF
    IF (q(bi) .GT. q(i)) THEN
      IF (q(i) .GT. q(ai)) THEN
        minqd = qd(ai)
        minq = q(ai)
      ELSE
        minqd = qd(i)
        minq = q(i)
      END IF
    ELSE IF (q(bi) .GT. q(ai)) THEN
      minqd = qd(ai)
      minq = q(ai)
    ELSE
      minqd = qd(bi)
      minq = q(bi)
    END IF
    IF ((flux(i)-q(bi))*(q(i)-flux(i)) .LT. 0.0) THEN
      IF ((((q(bi)-q(di))*(q(ai)-q(i)) .GE. 0.0 .OR. (q(bi)-q(di))*(q(di&
&         )-q(fi)) .LE. 0.0) .OR. (q(ai)-q(i))*(q(ci)-q(ai)) .LE. 0.0) &
&         .OR. (flux(i)-q(bi))*(q(bi)-q(di)) .LE. 0.0) THEN
        upbd(i) = 0.D0
        upb(i) = 0.0
        lowbd(i) = 0.D0
        lowb(i) = 0.0
      ELSE
        upbd(i) = maxqd - minqd
        upb(i) = maxq - minq
        lowbd(i) = maxqd - minqd
        lowb(i) = maxq - minq
      END IF
    ELSE
      upbd(i) = maxqd - minqd
      upb(i) = maxq - minq
      lowbd(i) = maxqd - minqd
      lowb(i) = maxq - minq
    END IF
  END DO
  qlinmaxd = 0.D0
  qlinmind = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    IF (ai .EQ. n + 1) ai = 1
    IF (q(bi) .GT. q(i)) THEN
      min1d = qd(i)
      min1 = q(i)
    ELSE
      min1d = qd(bi)
      min1 = q(bi)
    END IF
    x1d = min1d - lowbd(i)
    x1 = min1 - lowb(i)
    IF (x1 .LT. zeros) THEN
      qlinmind(i) = 0.D0
      qlinmin(i) = zeros
    ELSE
      qlinmind(i) = x1d
      qlinmin(i) = x1
    END IF
    IF (q(bi) .LT. q(i)) THEN
      max1d = qd(i)
      max1 = q(i)
    ELSE
      max1d = qd(bi)
      max1 = q(bi)
    END IF
    qlinmaxd(i) = max1d + upbd(i)
    qlinmax(i) = max1 + upb(i)
  END DO
  qlpd = 0.D0
  qlppd = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    IF (ai .EQ. n + 1) ai = 1
    IF (flux(i) .GT. qlinmax(i)) THEN
      qlpd(i) = qlinmaxd(i)
      qlp(i) = qlinmax(i)
    ELSE
      qlpd(i) = fluxd(i)
      qlp(i) = flux(i)
    END IF
    IF (qlp(i) .LT. qlinmin(i)) THEN
      qlppd(i) = qlinmind(i)
      qlpp(i) = qlinmin(i)
    ELSE
      qlppd(i) = qlpd(i)
      qlpp(i) = qlp(i)
    END IF
  END DO
  qnewmaxd(1:n) = qlinmaxd(1:n)
  qnewmax(1:n) = qlinmax(1:n)
  qnewmind(1:n) = qlinmind(1:n)
  qnewmin(1:n) = qlinmin(1:n)
  qoutmind = 0.D0
  qoutmaxd = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    IF (ai .EQ. n + 1) ai = 1
    qoutmind(i) = ((qd(i)+dt*vd*qlinmax(i)/dx+v*dt*qlinmaxd(i)/dx-&
&     qnewmaxd(i))*v*dt/dx-(q(i)+v*dt/dx*qlinmax(i)-qnewmax(i))*dt*vd/dx&
&     )/(v*dt/dx)**2
    qoutmin(i) = (q(i)+v*dt/dx*qlinmax(i)-qnewmax(i))/(v*dt/dx)
    qoutmaxd(i) = ((qd(i)+dt*vd*qlinmin(i)/dx+v*dt*qlinmind(i)/dx-&
&     qnewmind(i))*v*dt/dx-(q(i)+v*dt/dx*qlinmin(i)-qnewmin(i))*dt*vd/dx&
&     )/(v*dt/dx)**2
    qoutmax(i) = (q(i)+v*dt/dx*qlinmin(i)-qnewmin(i))/(v*dt/dx)
  END DO
  qrd = 0.D0
  qrpppd = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    IF (ai .EQ. n + 1) ai = 1
    IF (qlpp(ai) .GT. qoutmax(i)) THEN
      qrpppd(i) = qoutmaxd(i)
      qrppp(i) = qoutmax(i)
    ELSE
      qrpppd(i) = qlppd(ai)
      qrppp(i) = qlpp(ai)
    END IF
    IF (qrppp(i) .LT. qoutmin(i)) THEN
      qrd(i) = qoutmind(i)
      qr(i) = qoutmin(i)
    ELSE
      qrd(i) = qrpppd(i)
      qr(i) = qrppp(i)
    END IF
  END DO
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    IF (ai .EQ. n + 1) ai = 1
    fluxd(i) = qrd(bi)
    flux(i) = qr(bi)
  END DO
END SUBROUTINE RELAXEDLIMITER_TLM
! =========================================================

! =========================================================
SUBROUTINE universallimiter(flux,q,v,dt,dx,N)

IMPLICIT NONE
INTEGER, INTENT(in):: N
DOUBLE PRECISION, INTENT(in) :: dx, dt, q(N), v
DOUBLE PRECISION, INTENT(inout) :: flux(N)
DOUBLE PRECISION :: qLinmin(N), qLinmax(N), qLp(N), qLpp(N), &
	qnewmax(N), qnewmin(N), qoutmin(N), qoutmax(N), qRppp(N), qR(N)
INTEGER :: i, ai, bi



	DO i=1,N
		ai=i+1
		bi=i-1
		if (bi==0) then
			bi=N
		endif
		if (ai==N+1) then
			ai=1
		endif

		qLinmin(i) = min(q(bi),q(i)) 
		qLinmax(i) = max(q(bi),q(i)) 
	ENDDO
    


	DO i=1,N
		ai=i+1
		bi=i-1
		if (bi==0) then
			bi=N
		endif
		if (ai==N+1) then
			ai=1
		endif
		qLp(i) = min(flux(i),qLinmax(i)) 
		qLpp(i) = max(qLp(i),qLinmin(i)) 
	ENDDO
    
	qnewmax(1:N) = qLinmax(1:N) 
	qnewmin(1:N) = qLinmin(1:N)

	DO i=1,N
		ai=i+1
		bi=i-1
		if (bi==0) then
			bi=N
		endif
		if (ai==N+1) then
			ai=1
		endif
		qoutmin(i) = (q(i) + (v*dt/dx)*qLinmax(i) - qnewmax(i))/(v*dt/dx) 
		qoutmax(i) = (q(i) + (v*dt/dx)*qLinmin(i) - qnewmin(i))/(v*dt/dx)
	ENDDO
    


	DO i=1,N
		ai=i+1
		bi=i-1
		if (bi==0) then
			bi=N
		endif
		if (ai==N+1) then
			ai=1
		endif
		qRppp(i) = min(qLpp(ai),qoutmax(i)) 
		qR(i) = max(qRppp(i),qoutmin(i)) 
	ENDDO

	DO i=1,N
		ai=i+1
		bi=i-1
		if (bi==0) then
			bi=N
		endif
		if (ai==N+1) then
			ai=1
		endif
		flux(i) = qR(bi)
	ENDDO
    

END SUBROUTINE universallimiter
! =========================================================

! =========================================================
SUBROUTINE UNIVERSALLIMITER_TLM(flux, fluxd, q, qd, v, vd, dt, dx, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: dx, dt, q(n), v
  DOUBLE PRECISION, INTENT(IN) :: qd(n), vd
  DOUBLE PRECISION, INTENT(INOUT) :: flux(n)
  DOUBLE PRECISION, INTENT(INOUT) :: fluxd(n)
  DOUBLE PRECISION :: qlinmin(n), qlinmax(n), qlp(n), qlpp(n), qnewmax(n&
& ), qnewmin(n), qoutmin(n), qoutmax(n), qrppp(n), qr(n)
  DOUBLE PRECISION :: qlinmind(n), qlinmaxd(n), qlpd(n), qlppd(n), &
& qnewmaxd(n), qnewmind(n), qoutmind(n), qoutmaxd(n), qrpppd(n), qrd(n)
  INTEGER :: i, ai, bi
  INTRINSIC MIN
  INTRINSIC MAX
  qlinmaxd = 0.D0
  qlinmind = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    IF (ai .EQ. n + 1) ai = 1
    IF (q(bi) .GT. q(i)) THEN
      qlinmind(i) = qd(i)
      qlinmin(i) = q(i)
    ELSE
      qlinmind(i) = qd(bi)
      qlinmin(i) = q(bi)
    END IF
    IF (q(bi) .LT. q(i)) THEN
      qlinmaxd(i) = qd(i)
      qlinmax(i) = q(i)
    ELSE
      qlinmaxd(i) = qd(bi)
      qlinmax(i) = q(bi)
    END IF
  END DO
  qlpd = 0.D0
  qlppd = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    IF (ai .EQ. n + 1) ai = 1
    IF (flux(i) .GT. qlinmax(i)) THEN
      qlpd(i) = qlinmaxd(i)
      qlp(i) = qlinmax(i)
    ELSE
      qlpd(i) = fluxd(i)
      qlp(i) = flux(i)
    END IF
    IF (qlp(i) .LT. qlinmin(i)) THEN
      qlppd(i) = qlinmind(i)
      qlpp(i) = qlinmin(i)
    ELSE
      qlppd(i) = qlpd(i)
      qlpp(i) = qlp(i)
    END IF
  END DO
  qnewmaxd(1:n) = qlinmaxd(1:n)
  qnewmax(1:n) = qlinmax(1:n)
  qnewmind(1:n) = qlinmind(1:n)
  qnewmin(1:n) = qlinmin(1:n)
  qoutmind = 0.D0
  qoutmaxd = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    IF (ai .EQ. n + 1) ai = 1
    qoutmind(i) = ((qd(i)+dt*vd*qlinmax(i)/dx+v*dt*qlinmaxd(i)/dx-&
&     qnewmaxd(i))*v*dt/dx-(q(i)+v*dt/dx*qlinmax(i)-qnewmax(i))*dt*vd/dx&
&     )/(v*dt/dx)**2
    qoutmin(i) = (q(i)+v*dt/dx*qlinmax(i)-qnewmax(i))/(v*dt/dx)
    qoutmaxd(i) = ((qd(i)+dt*vd*qlinmin(i)/dx+v*dt*qlinmind(i)/dx-&
&     qnewmind(i))*v*dt/dx-(q(i)+v*dt/dx*qlinmin(i)-qnewmin(i))*dt*vd/dx&
&     )/(v*dt/dx)**2
    qoutmax(i) = (q(i)+v*dt/dx*qlinmin(i)-qnewmin(i))/(v*dt/dx)
  END DO
  qrd = 0.D0
  qrpppd = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    IF (ai .EQ. n + 1) ai = 1
    IF (qlpp(ai) .GT. qoutmax(i)) THEN
      qrpppd(i) = qoutmaxd(i)
      qrppp(i) = qoutmax(i)
    ELSE
      qrpppd(i) = qlppd(ai)
      qrppp(i) = qlpp(ai)
    END IF
    IF (qrppp(i) .LT. qoutmin(i)) THEN
      qrd(i) = qoutmind(i)
      qr(i) = qoutmin(i)
    ELSE
      qrd(i) = qrpppd(i)
      qr(i) = qrppp(i)
    END IF
  END DO
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    IF (ai .EQ. n + 1) ai = 1
    fluxd(i) = qrd(bi)
    flux(i) = qr(bi)
  END DO
END SUBROUTINE UNIVERSALLIMITER_TLM
! =========================================================

! =========================================================
SUBROUTINE SPECTRALADVECTION(qspec,qspecold,u,N,x,twodt)

IMPLICIT NONE

INTEGER, INTENT(in) :: N
DOUBLE PRECISION, INTENT(in) :: x(N), twodt, u
DOUBLE PRECISION, INTENT(inout) :: qspec(N), qspecold(N)

INTEGER :: j, k
DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384

! Coefficients etc
COMPLEX(8) :: minus1, sqrtm1, dbdx(N)
DOUBLE PRECISION :: kwave(N)
COMPLEX(8) :: qspect(N), dqspectdx(N), dqspecdx(N), qspecfftshift(N), qspecnew(N)


 !Square root of minus 1
 minus1=(-1.0,0.0)
 sqrtm1=sqrt(minus1)

 !Save spectral coefficient
 DO j=1, N
    kwave(j)=j-1
    dbdx(j)=kwave(j)*2.0*pi*sqrtm1
 ENDDO

 !Transform to spectral space
 do k=1, N/2+1
    qspect(k) = 0.0
    DO j=1, N
       qspect(k) = qspect(k) + qspec(j)*exp(-2.0*pi*sqrtm1*(k-1.0)*x(j))
    ENDDO
 enddo    

!optional Ensure parts that are zero are zero!
!e.g. imag(qspec(1)) = 0.0
!     imag(qspec(N/2+1)) = 0.0 

 !Derivative
 dqspectdx = dbdx*qspect

 !Inverse Transform    
 DO j=1, N
    dqspecdx(j) = 0.0
    do k=1, N/2+1
       if (k .eq. 1) then
          dqspecdx(j) = dqspecdx(j) + 1.0*(dqspectdx(k)*exp(2.0*pi*sqrtm1*(k-1.0)*x(j))/(1.0*N))
       else
          dqspecdx(j) = dqspecdx(j) + 2.0*(dqspectdx(k)*exp(2.0*pi*sqrtm1*(k-1.0)*x(j))/(1.0*N))
       endif
    enddo
 ENDDO

 !New Value
 do j=1, N        
    qspecnew(j) = qspecold(j) - u*twodt*real(dqspecdx(j))
 enddo        

 !Update old values
 qspecold=qspec
 qspec=qspecnew

END SUBROUTINE SPECTRALADVECTION
! =========================================================


! =======================================================================================
SUBROUTINE SPECTRALADVECTION_TLM(qspec, qspecd, qspecold, qspecoldd, u, ud , n, x, twodt)

  IMPLICIT NONE

  !Inputs
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: x(n), twodt
  DOUBLE PRECISION, INTENT(IN) :: u
  DOUBLE PRECISION, INTENT(IN) :: ud

  !Prognostic
  DOUBLE PRECISION, INTENT(INOUT) :: qspec (n), qspecold (n)
  DOUBLE PRECISION, INTENT(INOUT) :: qspecd(n), qspecoldd(n)

  !Locals
  INTEGER :: j, k
  INTEGER :: result_alias

  DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384
  DOUBLE PRECISION :: kwave(n)

  COMPLEX(8) :: earg
  COMPLEX(8) :: minus1, sqrtm1, dbdx(n)
  COMPLEX(8) :: qspect (n), dqspectdx (n), dqspecdx (n), qspecfftshift (n), qspecnew (n)
  COMPLEX(8) :: qspectd(n), dqspectdxd(n), dqspecdxd(n), qspecfftshiftd(n), qspecnewd(n)

  !Square root of minus 1
  minus1 = (-1.0,0.0)
  sqrtm1 = SQRT(minus1)

  !Save spectral coefficient
  DO j=1,n
    kwave(j) = j - 1
    dbdx(j) = kwave(j)*2.0*pi*sqrtm1
  END DO

  !Transform to spectral space
  DO k=1,N/2+1
    qspectd(k) = 0.0
    qspect (k) = 0.0
    DO j=1,n
      earg = -(2.0*pi*sqrtm1*(k-1.0)*x(j))
      qspectd(k) = qspectd(k) + EXP(earg)*qspecd(j)
      qspect (k) = qspect (k) + EXP(earg)*qspec (j)
    END DO
  END DO

  !Derivative 
  dqspectdxd = dbdx*qspectd
  dqspectdx  = dbdx*qspect

  !Inverse Transform    
  DO j=1,n
    dqspecdxd(j) = 0.0
    dqspecdx (j) = 0.0
    DO k=1,N/2+1
      earg = 2.0*pi*sqrtm1*(k-1.0)*x(j)
      if (k .eq. 1) then
         dqspecdxd(j) = dqspecdxd(j) + 1.0*EXP(earg)*dqspectdxd(k)/(1.0*n)
         dqspecdx (j) = dqspecdx (j) + 1.0*EXP(earg)*dqspectdx (k)/(1.0*n)
      else
         dqspecdxd(j) = dqspecdxd(j) + 2.0*EXP(earg)*dqspectdxd(k)/(1.0*n)
         dqspecdx (j) = dqspecdx (j) + 2.0*EXP(earg)*dqspectdx (k)/(1.0*n)
      endif
    END DO
  END DO

  !New Value
  DO j=1,n
    qspecnewd(j) = qspecoldd(j) - twodt*(ud*REAL(dqspecdx(j))+u*REAL(dqspecdxd(j)))
    qspecnew (j) = qspecold (j) - twodt*(u *REAL(dqspecdx(j))                     )
  END DO

  !Update old values
  qspecoldd = qspecd
  qspecold  = qspec
  qspecd = qspecnewd
  qspec  = qspecnew

END SUBROUTINE SPECTRALADVECTION_TLM
! =========================================================



! =========================================================
SUBROUTINE ppmpart11d(flux,z,u,dt,dx,N)

! subroutine for original PPM scheme in 1D
! but only limiting the first part

IMPLICIT NONE
INTEGER, INTENT(in):: N
DOUBLE PRECISION, INTENT(in) :: dx, dt, z(N), u
DOUBLE PRECISION, INTENT(out) :: flux(N)
DOUBLE PRECISION :: zplus(N), zminus(N), &
	zp(N), zminusnew(N), zplusnew(N), z6(N), &
	D2s, D2sL, D2sR, s, D2sLIM, ones, zeros, D2sMIN
	
INTEGER :: i, ai, bi, ci, di

	ones=1.0
	zeros=0.0

! flux to left 

	DO i=1,N
    		ai=i+1
		bi=i-1
		ci=i+2

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif

		if (ci==N+1) then
			ci=1
		endif

		if (ci==N+2) then
			ci=2
		endif

		! Fourth-Order Edge-Reconstruction
	
		zp(i)=(-z(ci)+7.0*z(ai)+7.0*z(i)-z(bi))/12.0

	ENDDO

	DO i=1,N
    		ai=i+1
		bi=i-1
		ci=i+2

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif

		if (ci==N+1) then
			ci=1
		endif

		if (ci==N+2) then
			ci=2
		endif

		! Limit edges if needed
        
        	if ((zp(i) - z(i))*(z(ai) - zp(i)) .lt. 0.0) then
            
            		D2s = 3.0*(z(i) - 2.0*zp(i) + z(ai))
            		D2sL = (z(bi) - 2.0*z(i) + z(ai))
            		D2sR = (z(i) - 2.0*z(ai) + z(ci))
            
            		s = sign(ones,D2s)
            
			D2sMIN = min(s*abs(D2sL),s*abs(D2sR),s*abs(D2s))
            		D2sLIM = s*max(D2sMIN,0.0)
            
            		zp(i) = 0.5*(z(i)+z(ai)) - (1.0/6.0)*D2sLIM
            
        	endif
    
	ENDDO

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif

		! cell edges are therefore
    
    		zminus(i) = zp(bi)
    		zplus(i) = zp(i)
    
	ENDDO

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif 
    
    		z6(i) = 6.0*(z(i) - 0.5*(zminus(i)+zplus(i)))
    
	ENDDO

	! calculate flux based on velocity direction

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif
    
    	if (u .ge. 0.0) then
        	flux(i) = zplus(bi) - 0.5*(u*dt/dx)*( -zminus(bi)+zplus(bi) - z6(bi)*(1.0-(2.0/3.0)*(u*dt/dx)) )
    	else
        	flux(i) = zminus(i) + 0.5*(abs(u)*dt/dx)*( zplus(i)-zminus(i) + z6(i)*(1.0-(2.0/3.0)*(abs(u)*dt/dx)) )
    	endif
    
	ENDDO

	

END SUBROUTINE ppmpart11d
! =========================================================

! =========================================================
SUBROUTINE PPMPART11D_TLM(flux, fluxd, z, zd, u, ud, dt, dx, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: dx, dt, z(n), u
  DOUBLE PRECISION, INTENT(IN) :: zd(n), ud
  DOUBLE PRECISION, INTENT(OUT) :: flux(n)
  DOUBLE PRECISION, INTENT(OUT) :: fluxd(n)
  DOUBLE PRECISION :: zplus(n), zminus(n), zp(n), zminusnew(n), zplusnew&
& (n), z6(n), d2s, d2sl, d2sr, s, d2slim, ones, zeros, d2smin
  DOUBLE PRECISION :: zplusd(n), zminusd(n), zpd(n), z6d(n), d2sd, d2sld&
& , d2srd, d2slimd, d2smind
  INTEGER :: i, ai, bi, ci, di
  INTRINSIC SIGN
  INTRINSIC ABS
  INTRINSIC MIN
  INTRINSIC MAX
  DOUBLE PRECISION :: abs1d
  DOUBLE PRECISION :: abs4d
  DOUBLE PRECISION :: x1
  DOUBLE PRECISION :: z1d
  DOUBLE PRECISION :: abs0d
  DOUBLE PRECISION :: max1d
  DOUBLE PRECISION :: abs3d
  DOUBLE PRECISION :: x1d
  DOUBLE PRECISION :: z1
  DOUBLE PRECISION :: abs4
  DOUBLE PRECISION :: abs3
  DOUBLE PRECISION :: abs2
  DOUBLE PRECISION :: abs2d
  DOUBLE PRECISION :: abs1
  DOUBLE PRECISION :: abs0
  DOUBLE PRECISION :: max1
  DOUBLE PRECISION :: y1
  DOUBLE PRECISION :: y1d
  ones = 1.0
  zeros = 0.0
  zpd = 0.D0
! flux to left 
  DO i=1,n
    ai = i + 1
    bi = i - 1
    ci = i + 2
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (ci .EQ. n + 1) ci = 1
    IF (ci .EQ. n + 2) ci = 2
! Fourth-Order Edge-Reconstruction
    zpd(i) = (7.0*zd(ai)-zd(ci)+7.0*zd(i)-zd(bi))/12.0
    zp(i) = (-z(ci)+7.0*z(ai)+7.0*z(i)-z(bi))/12.0
  END DO
  DO i=1,n
    ai = i + 1
    bi = i - 1
    ci = i + 2
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (ci .EQ. n + 1) ci = 1
    IF (ci .EQ. n + 2) ci = 2
! Limit edges if needed
    IF ((zp(i)-z(i))*(z(ai)-zp(i)) .LT. 0.0) THEN
      d2sd = 3.0*(zd(i)-2.0*zpd(i)+zd(ai))
      d2s = 3.0*(z(i)-2.0*zp(i)+z(ai))
      d2sld = zd(bi) - 2.0*zd(i) + zd(ai)
      d2sl = z(bi) - 2.0*z(i) + z(ai)
      d2srd = zd(i) - 2.0*zd(ai) + zd(ci)
      d2sr = z(i) - 2.0*z(ai) + z(ci)
      s = SIGN(ones, d2s)
      IF (d2sl .GE. 0.) THEN
        abs1d = d2sld
        abs1 = d2sl
      ELSE
        abs1d = -d2sld
        abs1 = -d2sl
      END IF
      x1d = s*abs1d
      x1 = s*abs1
      IF (d2sr .GE. 0.) THEN
        abs2d = d2srd
        abs2 = d2sr
      ELSE
        abs2d = -d2srd
        abs2 = -d2sr
      END IF
      y1d = s*abs2d
      y1 = s*abs2
      IF (d2s .GE. 0.) THEN
        abs3d = d2sd
        abs3 = d2s
      ELSE
        abs3d = -d2sd
        abs3 = -d2s
      END IF
      z1d = s*abs3d
      z1 = s*abs3
      IF (x1 .GT. y1) THEN
        IF (y1 .GT. z1) THEN
          d2smind = z1d
          d2smin = z1
        ELSE
          d2smind = y1d
          d2smin = y1
        END IF
      ELSE IF (x1 .GT. z1) THEN
        d2smind = z1d
        d2smin = z1
      ELSE
        d2smind = x1d
        d2smin = x1
      END IF
      IF (d2smin .LT. 0.0) THEN
        max1 = 0.0
        max1d = 0.D0
      ELSE
        max1d = d2smind
        max1 = d2smin
      END IF
      d2slimd = s*max1d
      d2slim = s*max1
      zpd(i) = 0.5*(zd(i)+zd(ai)) - d2slimd/6.0
      zp(i) = 0.5*(z(i)+z(ai)) - 1.0/6.0*d2slim
    END IF
  END DO
  zminusd = 0.D0
  zplusd = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
! cell edges are therefore
    zminusd(i) = zpd(bi)
    zminus(i) = zp(bi)
    zplusd(i) = zpd(i)
    zplus(i) = zp(i)
  END DO
  z6d = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    z6d(i) = 6.0*(zd(i)-0.5*(zminusd(i)+zplusd(i)))
    z6(i) = 6.0*(z(i)-0.5*(zminus(i)+zplus(i)))
  END DO
  fluxd = 0.D0
! calculate flux based on velocity direction
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (u .GE. 0.0) THEN
      fluxd(i) = zplusd(bi) - 0.5*(dt*ud*(-zminus(bi)+zplus(bi)-z6(bi)*(&
&       1.0-2.0/3.0*(u*dt/dx)))/dx+u*dt*(zplusd(bi)-zminusd(bi)-z6d(bi)*&
&       (1.0-2.0/3.0*(u*dt/dx))+z6(bi)*2.0*dt*ud/(3.0*dx))/dx)
      flux(i) = zplus(bi) - 0.5*(u*dt/dx)*(-zminus(bi)+zplus(bi)-z6(bi)*&
&       (1.0-2.0/3.0*(u*dt/dx)))
    ELSE
      IF (u .GE. 0.) THEN
        abs0d = ud
        abs0 = u
      ELSE
        abs0d = -ud
        abs0 = -u
      END IF
      IF (u .GE. 0.) THEN
        abs4d = ud
        abs4 = u
      ELSE
        abs4d = -ud
        abs4 = -u
      END IF
      fluxd(i) = zminusd(i) + 0.5*(dt*abs0d*(zplus(i)-zminus(i)+z6(i)*(&
&       1.0-2.0/3.0*(abs4*dt/dx)))/dx+abs0*dt*(zplusd(i)-zminusd(i)+z6d(&
&       i)*(1.0-2.0/3.0*(abs4*dt/dx))-z6(i)*2.0*dt*abs4d/(3.0*dx))/dx)
      flux(i) = zminus(i) + 0.5*(abs0*dt/dx)*(zplus(i)-zminus(i)+z6(i)*(&
&       1.0-2.0/3.0*(abs4*dt/dx)))
    END IF
  END DO
END SUBROUTINE PPMPART11D_TLM
! =========================================================

! =========================================================
SUBROUTINE ppmpart21d(flux,z,u,dt,dx,N)

! subroutine for original PPM scheme in 1D
! but only limiting second part

IMPLICIT NONE
INTEGER, INTENT(in):: N
DOUBLE PRECISION, INTENT(in) :: dx, dt, z(N), u
DOUBLE PRECISION, INTENT(out) :: flux(N)
DOUBLE PRECISION :: zplus(N), zminus(N), &
	zp(N), zminusnew(N), zplusnew(N), z6(N), &
	D2s, D2sL, D2sR, s, D2sLIM, ones, zeros, D2sMIN
	
INTEGER :: i, ai, bi, ci, di

	ones=1.0
	zeros=0.0

! flux to left 

	DO i=1,N
    		ai=i+1
		bi=i-1
		ci=i+2

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif

		if (ci==N+1) then
			ci=1
		endif

		if (ci==N+2) then
			ci=2
		endif

		! Fourth-Order Edge-Reconstruction
	
		zp(i)=(-z(ci)+7.0*z(ai)+7.0*z(i)-z(bi))/12.0

	ENDDO

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif

		! cell edges are therefore
    
    		zminus(i) = zp(bi)
    		zminusnew(i) = zp(bi)
    		zplus(i) = zp(i)
    		zplusnew(i) = zp(i)
    
	ENDDO

	DO i=1,N
    		ai=i+1
		bi=i-1
		ci=i+2

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif

		if (ci==N+1) then
			ci=1
		endif

		if (ci==N+2) then
			ci=2
		endif

		! Limit subgrid reconstruction

      		if ( (zplus(i)-z(i))*(zminus(i)-z(i)) .gt. 0.0) then 
            		zplusnew(i)=z(i)
            		zminusnew(i)=z(i)
        	else
            
            		if ( (zplus(i)-zminus(i))*(z(i)-0.5*(zplus(i)+zminus(i))) & 
				.gt. ((zplus(i)-zminus(i))**2.0)/6.0 ) then
                		zminusnew(i) = 3.0*z(i) - 2.0*zplus(i)
            		endif
                        
            		if ((zplus(i)-zminus(i))*(z(i)-0.5*(zplus(i)+zminus(i))) & 
				.lt. -((zplus(i)-zminus(i))**2)/6.0 ) then
                		zplusnew(i) = 3.0*z(i) - 2.0*zminus(i)
            		endif
                        
        	endif
    
	ENDDO

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif
    
    		zminus(i)=zminusnew(i)
		zplus(i)=zplusnew(i)  
    
    		z6(i) = 6.0*(z(i) - 0.5*(zminus(i)+zplus(i)))
    
	ENDDO

	! calculate flux based on velocity direction

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif
    
    	if (u .ge. 0.0) then
        	flux(i) = zplus(bi) - 0.5*(u*dt/dx)*( -zminus(bi)+zplus(bi) - z6(bi)*(1.0-(2.0/3.0)*(u*dt/dx)) )
    	else
        	flux(i) = zminus(i) + 0.5*(abs(u)*dt/dx)*( zplus(i)-zminus(i) + z6(i)*(1.0-(2.0/3.0)*(abs(u)*dt/dx)) )
    	endif
    
	ENDDO

	

END SUBROUTINE ppmpart21d
! =========================================================

! =========================================================
SUBROUTINE PPMPART21D_TLM(flux, fluxd, z, zd, u, ud, dt, dx, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: dx, dt, z(n), u
  DOUBLE PRECISION, INTENT(IN) :: zd(n), ud
  DOUBLE PRECISION, INTENT(OUT) :: flux(n)
  DOUBLE PRECISION, INTENT(OUT) :: fluxd(n)
  DOUBLE PRECISION :: zplus(n), zminus(n), zp(n), zminusnew(n), zplusnew&
& (n), z6(n), d2s, d2sl, d2sr, s, d2slim, ones, zeros, d2smin
  DOUBLE PRECISION :: zplusd(n), zminusd(n), zpd(n), zminusnewd(n), &
& zplusnewd(n), z6d(n)
  INTEGER :: i, ai, bi, ci, di
  INTRINSIC ABS
  DOUBLE PRECISION :: abs1d
  DOUBLE PRECISION :: abs0d
  DOUBLE PRECISION :: abs1
  DOUBLE PRECISION :: abs0
  ones = 1.0
  zeros = 0.0
  zpd = 0.D0
! flux to left 
  DO i=1,n
    ai = i + 1
    bi = i - 1
    ci = i + 2
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (ci .EQ. n + 1) ci = 1
    IF (ci .EQ. n + 2) ci = 2
! Fourth-Order Edge-Reconstruction
    zpd(i) = (7.0*zd(ai)-zd(ci)+7.0*zd(i)-zd(bi))/12.0
    zp(i) = (-z(ci)+7.0*z(ai)+7.0*z(i)-z(bi))/12.0
  END DO
  zminusd = 0.D0
  zplusnewd = 0.D0
  zplusd = 0.D0
  zminusnewd = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
! cell edges are therefore
    zminusd(i) = zpd(bi)
    zminus(i) = zp(bi)
    zminusnewd(i) = zpd(bi)
    zminusnew(i) = zp(bi)
    zplusd(i) = zpd(i)
    zplus(i) = zp(i)
    zplusnewd(i) = zpd(i)
    zplusnew(i) = zp(i)
  END DO
  DO i=1,n
    ai = i + 1
    bi = i - 1
    ci = i + 2
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (ci .EQ. n + 1) ci = 1
    IF (ci .EQ. n + 2) ci = 2
! Limit subgrid reconstruction
    IF ((zplus(i)-z(i))*(zminus(i)-z(i)) .GT. 0.0) THEN
      zplusnewd(i) = zd(i)
      zplusnew(i) = z(i)
      zminusnewd(i) = zd(i)
      zminusnew(i) = z(i)
    ELSE
      IF ((zplus(i)-zminus(i))*(z(i)-0.5*(zplus(i)+zminus(i))) .GT. (&
&         zplus(i)-zminus(i))**2.0/6.0) THEN
        zminusnewd(i) = 3.0*zd(i) - 2.0*zplusd(i)
        zminusnew(i) = 3.0*z(i) - 2.0*zplus(i)
      END IF
      IF ((zplus(i)-zminus(i))*(z(i)-0.5*(zplus(i)+zminus(i))) .LT. -((&
&         zplus(i)-zminus(i))**2/6.0)) THEN
        zplusnewd(i) = 3.0*zd(i) - 2.0*zminusd(i)
        zplusnew(i) = 3.0*z(i) - 2.0*zminus(i)
      END IF
    END IF
  END DO
  z6d = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    zminusd(i) = zminusnewd(i)
    zminus(i) = zminusnew(i)
    zplusd(i) = zplusnewd(i)
    zplus(i) = zplusnew(i)
    z6d(i) = 6.0*(zd(i)-0.5*(zminusd(i)+zplusd(i)))
    z6(i) = 6.0*(z(i)-0.5*(zminus(i)+zplus(i)))
  END DO
  fluxd = 0.D0
! calculate flux based on velocity direction
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (u .GE. 0.0) THEN
      fluxd(i) = zplusd(bi) - 0.5*(dt*ud*(-zminus(bi)+zplus(bi)-z6(bi)*(&
&       1.0-2.0/3.0*(u*dt/dx)))/dx+u*dt*(zplusd(bi)-zminusd(bi)-z6d(bi)*&
&       (1.0-2.0/3.0*(u*dt/dx))+z6(bi)*2.0*dt*ud/(3.0*dx))/dx)
      flux(i) = zplus(bi) - 0.5*(u*dt/dx)*(-zminus(bi)+zplus(bi)-z6(bi)*&
&       (1.0-2.0/3.0*(u*dt/dx)))
    ELSE
      IF (u .GE. 0.) THEN
        abs0d = ud
        abs0 = u
      ELSE
        abs0d = -ud
        abs0 = -u
      END IF
      IF (u .GE. 0.) THEN
        abs1d = ud
        abs1 = u
      ELSE
        abs1d = -ud
        abs1 = -u
      END IF
      fluxd(i) = zminusd(i) + 0.5*(dt*abs0d*(zplus(i)-zminus(i)+z6(i)*(&
&       1.0-2.0/3.0*(abs1*dt/dx)))/dx+abs0*dt*(zplusd(i)-zminusd(i)+z6d(&
&       i)*(1.0-2.0/3.0*(abs1*dt/dx))-z6(i)*2.0*dt*abs1d/(3.0*dx))/dx)
      flux(i) = zminus(i) + 0.5*(abs0*dt/dx)*(zplus(i)-zminus(i)+z6(i)*(&
&       1.0-2.0/3.0*(abs1*dt/dx)))
    END IF
  END DO
END SUBROUTINE PPMPART21D_TLM
! =========================================================

! =========================================================
SUBROUTINE ppmoriginal1d(flux,z,u,dt,dx,N)

! subroutine for original PPM scheme in 1D

IMPLICIT NONE
INTEGER, INTENT(in):: N
DOUBLE PRECISION, INTENT(in) :: dx, dt, z(N), u
DOUBLE PRECISION, INTENT(out) :: flux(N)
DOUBLE PRECISION :: zplus(N), zminus(N), &
	zp(N), zminusnew(N), zplusnew(N), z6(N), &
	D2s, D2sL, D2sR, s, D2sLIM, ones, zeros, D2sMIN
	
INTEGER :: i, ai, bi, ci, di

	ones=1.0
	zeros=0.0

! flux to left 

	DO i=1,N
    		ai=i+1
		bi=i-1
		ci=i+2

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif

		if (ci==N+1) then
			ci=1
		endif

		if (ci==N+2) then
			ci=2
		endif

		! Fourth-Order Edge-Reconstruction
	
		zp(i)=(-z(ci)+7.0*z(ai)+7.0*z(i)-z(bi))/12.0

	ENDDO

	DO i=1,N
    		ai=i+1
		bi=i-1
		ci=i+2

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif

		if (ci==N+1) then
			ci=1
		endif

		if (ci==N+2) then
			ci=2
		endif

		! Limit edges if needed
        
        	if ((zp(i) - z(i))*(z(ai) - zp(i)) .lt. 0.0) then
            
            		D2s = 3.0*(z(i) - 2.0*zp(i) + z(ai))
            		D2sL = (z(bi) - 2.0*z(i) + z(ai))
            		D2sR = (z(i) - 2.0*z(ai) + z(ci))
            
            		s = sign(ones,D2s)
            
			D2sMIN = min(s*abs(D2sL),s*abs(D2sR),s*abs(D2s))
            		D2sLIM = s*max(D2sMIN,0.0)
            
            		zp(i) = 0.5*(z(i)+z(ai)) - (1.0/6.0)*D2sLIM
            
        	endif
    
	ENDDO

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif

		! cell edges are therefore
    
    		zminus(i) = zp(bi)
    		zminusnew(i) = zp(bi)
    		zplus(i) = zp(i)
    		zplusnew(i) = zp(i)
    
	ENDDO

	DO i=1,N
    		ai=i+1
		bi=i-1
		ci=i+2

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif

		if (ci==N+1) then
			ci=1
		endif

		if (ci==N+2) then
			ci=2
		endif

		! Limit subgrid reconstruction

      		if ( (zplus(i)-z(i))*(zminus(i)-z(i)) .gt. 0.0) then 
            		zplusnew(i)=z(i)
            		zminusnew(i)=z(i)
        	else
            
            		if ( (zplus(i)-zminus(i))*(z(i)-0.5*(zplus(i)+zminus(i))) & 
				.gt. ((zplus(i)-zminus(i))**2.0)/6.0 ) then
                		zminusnew(i) = 3.0*z(i) - 2.0*zplus(i)
            		endif
                        
            		if ((zplus(i)-zminus(i))*(z(i)-0.5*(zplus(i)+zminus(i))) & 
				.lt. -((zplus(i)-zminus(i))**2)/6.0 ) then
                		zplusnew(i) = 3.0*z(i) - 2.0*zminus(i)
            		endif
                        
        	endif
    
	ENDDO

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif
    
    		zminus(i)=zminusnew(i)
		zplus(i)=zplusnew(i)  
    
    		z6(i) = 6.0*(z(i) - 0.5*(zminus(i)+zplus(i)))
    
	ENDDO

	! calculate flux based on velocity direction

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif
    
    	if (u .ge. 0.0) then
        	flux(i) = zplus(bi) - 0.5*(u*dt/dx)*( -zminus(bi)+zplus(bi) - z6(bi)*(1.0-(2.0/3.0)*(u*dt/dx)) )
    	else
        	flux(i) = zminus(i) + 0.5*(abs(u)*dt/dx)*( zplus(i)-zminus(i) + z6(i)*(1.0-(2.0/3.0)*(abs(u)*dt/dx)) )
    	endif
    
	ENDDO

	

END SUBROUTINE ppmoriginal1d
! =========================================================

! =========================================================
SUBROUTINE PPMORIGINAL1D_TLM(flux, fluxd, z, zd, u, ud, dt, dx, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: dx, dt, z(n), u
  DOUBLE PRECISION, INTENT(IN) :: zd(n), ud
  DOUBLE PRECISION, INTENT(OUT) :: flux(n)
  DOUBLE PRECISION, INTENT(OUT) :: fluxd(n)
  DOUBLE PRECISION :: zplus(n), zminus(n), zp(n), zminusnew(n), zplusnew&
& (n), z6(n), d2s, d2sl, d2sr, s, d2slim, ones, zeros, d2smin
  DOUBLE PRECISION :: zplusd(n), zminusd(n), zpd(n), zminusnewd(n), &
& zplusnewd(n), z6d(n), d2sd, d2sld, d2srd, d2slimd, d2smind
  INTEGER :: i, ai, bi, ci, di
  INTRINSIC SIGN
  INTRINSIC ABS
  INTRINSIC MIN
  INTRINSIC MAX
  DOUBLE PRECISION :: abs1d
  DOUBLE PRECISION :: abs4d
  DOUBLE PRECISION :: x1
  DOUBLE PRECISION :: z1d
  DOUBLE PRECISION :: abs0d
  DOUBLE PRECISION :: max1d
  DOUBLE PRECISION :: abs3d
  DOUBLE PRECISION :: x1d
  DOUBLE PRECISION :: z1
  DOUBLE PRECISION :: abs4
  DOUBLE PRECISION :: abs3
  DOUBLE PRECISION :: abs2
  DOUBLE PRECISION :: abs2d
  DOUBLE PRECISION :: abs1
  DOUBLE PRECISION :: abs0
  DOUBLE PRECISION :: max1
  DOUBLE PRECISION :: y1
  DOUBLE PRECISION :: y1d
  ones = 1.0
  zeros = 0.0
  zpd = 0.D0
! flux to left 
  DO i=1,n
    ai = i + 1
    bi = i - 1
    ci = i + 2
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (ci .EQ. n + 1) ci = 1
    IF (ci .EQ. n + 2) ci = 2
! Fourth-Order Edge-Reconstruction
    zpd(i) = (7.0*zd(ai)-zd(ci)+7.0*zd(i)-zd(bi))/12.0
    zp(i) = (-z(ci)+7.0*z(ai)+7.0*z(i)-z(bi))/12.0
  END DO
  DO i=1,n
    ai = i + 1
    bi = i - 1
    ci = i + 2
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (ci .EQ. n + 1) ci = 1
    IF (ci .EQ. n + 2) ci = 2
! Limit edges if needed
    IF ((zp(i)-z(i))*(z(ai)-zp(i)) .LT. 0.0) THEN
      d2sd = 3.0*(zd(i)-2.0*zpd(i)+zd(ai))
      d2s = 3.0*(z(i)-2.0*zp(i)+z(ai))
      d2sld = zd(bi) - 2.0*zd(i) + zd(ai)
      d2sl = z(bi) - 2.0*z(i) + z(ai)
      d2srd = zd(i) - 2.0*zd(ai) + zd(ci)
      d2sr = z(i) - 2.0*z(ai) + z(ci)
      s = SIGN(ones, d2s)
      IF (d2sl .GE. 0.) THEN
        abs1d = d2sld
        abs1 = d2sl
      ELSE
        abs1d = -d2sld
        abs1 = -d2sl
      END IF
      x1d = s*abs1d
      x1 = s*abs1
      IF (d2sr .GE. 0.) THEN
        abs2d = d2srd
        abs2 = d2sr
      ELSE
        abs2d = -d2srd
        abs2 = -d2sr
      END IF
      y1d = s*abs2d
      y1 = s*abs2
      IF (d2s .GE. 0.) THEN
        abs3d = d2sd
        abs3 = d2s
      ELSE
        abs3d = -d2sd
        abs3 = -d2s
      END IF
      z1d = s*abs3d
      z1 = s*abs3
      IF (x1 .GT. y1) THEN
        IF (y1 .GT. z1) THEN
          d2smind = z1d
          d2smin = z1
        ELSE
          d2smind = y1d
          d2smin = y1
        END IF
      ELSE IF (x1 .GT. z1) THEN
        d2smind = z1d
        d2smin = z1
      ELSE
        d2smind = x1d
        d2smin = x1
      END IF
      IF (d2smin .LT. 0.0) THEN
        max1 = 0.0
        max1d = 0.D0
      ELSE
        max1d = d2smind
        max1 = d2smin
      END IF
      d2slimd = s*max1d
      d2slim = s*max1
      zpd(i) = 0.5*(zd(i)+zd(ai)) - d2slimd/6.0
      zp(i) = 0.5*(z(i)+z(ai)) - 1.0/6.0*d2slim
    END IF
  END DO
  zminusd = 0.D0
  zplusnewd = 0.D0
  zplusd = 0.D0
  zminusnewd = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
! cell edges are therefore
    zminusd(i) = zpd(bi)
    zminus(i) = zp(bi)
    zminusnewd(i) = zpd(bi)
    zminusnew(i) = zp(bi)
    zplusd(i) = zpd(i)
    zplus(i) = zp(i)
    zplusnewd(i) = zpd(i)
    zplusnew(i) = zp(i)
  END DO
  DO i=1,n
    ai = i + 1
    bi = i - 1
    ci = i + 2
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (ci .EQ. n + 1) ci = 1
    IF (ci .EQ. n + 2) ci = 2
! Limit subgrid reconstruction
    IF ((zplus(i)-z(i))*(zminus(i)-z(i)) .GT. 0.0) THEN
      zplusnewd(i) = zd(i)
      zplusnew(i) = z(i)
      zminusnewd(i) = zd(i)
      zminusnew(i) = z(i)
    ELSE
      IF ((zplus(i)-zminus(i))*(z(i)-0.5*(zplus(i)+zminus(i))) .GT. (&
&         zplus(i)-zminus(i))**2.0/6.0) THEN
        zminusnewd(i) = 3.0*zd(i) - 2.0*zplusd(i)
        zminusnew(i) = 3.0*z(i) - 2.0*zplus(i)
      END IF
      IF ((zplus(i)-zminus(i))*(z(i)-0.5*(zplus(i)+zminus(i))) .LT. -((&
&         zplus(i)-zminus(i))**2/6.0)) THEN
        zplusnewd(i) = 3.0*zd(i) - 2.0*zminusd(i)
        zplusnew(i) = 3.0*z(i) - 2.0*zminus(i)
      END IF
    END IF
  END DO
  z6d = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    zminusd(i) = zminusnewd(i)
    zminus(i) = zminusnew(i)
    zplusd(i) = zplusnewd(i)
    zplus(i) = zplusnew(i)
    z6d(i) = 6.0*(zd(i)-0.5*(zminusd(i)+zplusd(i)))
    z6(i) = 6.0*(z(i)-0.5*(zminus(i)+zplus(i)))
  END DO
  fluxd = 0.D0
! calculate flux based on velocity direction
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (u .GE. 0.0) THEN
      fluxd(i) = zplusd(bi) - 0.5*(dt*ud*(-zminus(bi)+zplus(bi)-z6(bi)*(&
&       1.0-2.0/3.0*(u*dt/dx)))/dx+u*dt*(zplusd(bi)-zminusd(bi)-z6d(bi)*&
&       (1.0-2.0/3.0*(u*dt/dx))+z6(bi)*2.0*dt*ud/(3.0*dx))/dx)
      flux(i) = zplus(bi) - 0.5*(u*dt/dx)*(-zminus(bi)+zplus(bi)-z6(bi)*&
&       (1.0-2.0/3.0*(u*dt/dx)))
    ELSE
      IF (u .GE. 0.) THEN
        abs0d = ud
        abs0 = u
      ELSE
        abs0d = -ud
        abs0 = -u
      END IF
      IF (u .GE. 0.) THEN
        abs4d = ud
        abs4 = u
      ELSE
        abs4d = -ud
        abs4 = -u
      END IF
      fluxd(i) = zminusd(i) + 0.5*(dt*abs0d*(zplus(i)-zminus(i)+z6(i)*(&
&       1.0-2.0/3.0*(abs4*dt/dx)))/dx+abs0*dt*(zplusd(i)-zminusd(i)+z6d(&
&       i)*(1.0-2.0/3.0*(abs4*dt/dx))-z6(i)*2.0*dt*abs4d/(3.0*dx))/dx)
      flux(i) = zminus(i) + 0.5*(abs0*dt/dx)*(zplus(i)-zminus(i)+z6(i)*(&
&       1.0-2.0/3.0*(abs4*dt/dx)))
    END IF
  END DO
END SUBROUTINE PPMORIGINAL1D_TLM
! =========================================================

! =========================================================
SUBROUTINE ppmcs1d(flux,z,u,dt,dx,N)

! subroutine for PPM with CS limiter in 1D

IMPLICIT NONE
INTEGER, INTENT(in):: N
DOUBLE PRECISION, INTENT(in) :: dx, dt, z(N), u
DOUBLE PRECISION, INTENT(out) :: flux(N)
DOUBLE PRECISION :: zplus(N), zminus(N), &
	zp(N), zminusnew(N), zplusnew(N), z6(N), &
	D2s, D2sL, D2sR, D2sC, s, D2sLIM, ones, zeros, D2sMIN, &
	alphaplus(N), alphaminus(N), alphaminusnew(N), alphaplusnew(N), &
	DsfaceM, DsfaceP, DsfaceMIN, DsccM, DsccP, DsP, DsM, &
	delIext, dels, DsccMIN
	
INTEGER :: i, ai, bi, ci, di, extreme(N)

	ones=1.0
	zeros=0.0

! flux to left 

	DO i=1,N
    		ai=i+1
		bi=i-1
		ci=i+2

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif

		if (ci==N+1) then
			ci=1
		endif

		if (ci==N+2) then
			ci=2
		endif

		! Fourth-Order Edge-Reconstruction
	
		zp(i)=(-z(ci)+7.0*z(ai)+7.0*z(i)-z(bi))/12.0

	ENDDO

	DO i=1,N
    		ai=i+1
		bi=i-1
		ci=i+2

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif

		if (ci==N+1) then
			ci=1
		endif

		if (ci==N+2) then
			ci=2
		endif

		! Limit edges if needed
        
        	if ((zp(i) - z(i))*(z(ai) - zp(i)) .lt. 0.0) then
            
            		D2s = 3.0*(z(i) - 2.0*zp(i) + z(ai))
            		D2sL = (z(bi) - 2.0*z(i) + z(ai))
            		D2sR = (z(i) - 2.0*z(ai) + z(ci))
            
            		s = sign(ones,D2s)
            
			D2sMIN = min(1.25*s*abs(D2sL),1.25*s*abs(D2sR),s*abs(D2s))
            		D2sLIM = s*max(D2sMIN,0.0)
            
            		zp(i) = 0.5*(z(i)+z(ai)) - (1.0/6.0)*D2sLIM
            
        	endif
    
	ENDDO

	! Is it an extremum

	extreme(1:N)=0

	DO i=1,N
    		ai=i+1
		bi=i-1
		ci=i+2
                di=i-2

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif

		if (ci==N+1) then
			ci=1
		endif

		if (ci==N+2) then
			ci=2
		endif

		if (di==0) then
			di=N
		endif

		if (di==-1) then
			di=N-1
		endif
        
        	alphaplus(i) = zp(i) - z(i)
        	alphaminus(i) = zp(bi) - z(i)
        	alphaplusnew(i) = alphaplus(i)
        	alphaminusnew(i) = alphaminus(i)
        
        	if ( alphaplus(i)*alphaminus(i) .ge. 0.0 ) then
            
            		extreme(i)=1
            
        	endif
        
        	if ( abs(alphaplus(i))>2.0*abs(alphaminus(i)) ) then
        
            		DsfaceM = zp(bi) - zp(di)
            		DsfaceP = zp(ai) - zp(bi)
            	
            		DsfaceMIN = min(abs(DsfaceM),abs(DsfaceP))
    
            		DsccM = z(i) - z(bi)
            		DsccP = z(ai) - z(i)
            
            		DsccMIN = min(abs(DsccM),abs(DsccP))
            
            		if (DsfaceMIN .ge. DsccMIN) then
                
                		DsP = DsfaceP
                		DsM = DsfaceM
                
            		else
                
                		DsP = DsccP
                		DsM = DsccM
                
            		endif
            
            		if ( DsP*DsM .le. 0.0 ) then
                
                		extreme(i)=1
                
            		endif
            
        	endif
        
        	if ( abs(alphaminus(i))>2.0*abs(alphaplus(i)) ) then
        
            		DsfaceM = zp(bi) - zp(di)
            		DsfaceP = zp(ai) - zp(bi)
            
            		DsfaceMIN = min(abs(DsfaceM),abs(DsfaceP))
    
            		DsccM = z(i) - z(bi)
            		DsccP = z(ai) - z(i)
            
            		DsccMIN = min(abs(DsccM),abs(DsccP))
            
            		if (DsfaceMIN .ge. DsccMIN) then
                
                		DsP = DsfaceP
                		DsM = DsfaceM
                
            		else
                
                		DsP = DsccP
                		DsM = DsccM
                
            		endif
            
            		if (DsP*DsM .le. 0.0) then
                
                		extreme(i)=1
                
            		endif
                        
        	endif		
    
	ENDDO

	! At Extremum...

	DO i=1,N
    		ai=i+1
		bi=i-1
		ci=i+2
		di=i-2		

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif

		if (ci==N+1) then
			ci=1
		endif

		if (ci==N+2) then
			ci=2
		endif

		if (di==0) then
			di=N
		endif

		if (di==-1) then
			di=N-1
		endif

		if (extreme(i) .eq. 1) then

			D2s = 6.0*(alphaplus(i) - alphaminus(i))           
            		D2sL = (z(di) - 2.0*z(bi) + z(i))
            		D2sR = (z(i) - 2.0*z(ai) + z(ci))
            		D2sC = (z(bi) - 2.0*z(i) + z(i))
            
            		s = sign(ones,D2s)
            
			D2sMIN = min(1.25*s*abs(D2sL),1.25*s*abs(D2sR),1.25*s*abs(D2sC),s*abs(D2s))
            		D2sLIM = max(D2sMIN,0.0)
                        
            		alphaplusnew(i) = (alphaplusnew(i)*D2sLIM)/(max(D2s,10.0**(-10.0)))
            		alphaminusnew(i) = (alphaminusnew(i)*D2sLIM)/(max(D2s,10.0**(-10.0)))

		else
            
            		if ( abs(alphaplusnew(i)) .gt. 2.0*abs(alphaminusnew(i)) ) then
                
                		s = sign(ones,alphaminusnew(i))
                
                		delIext = -( alphaplusnew(i)**2.0)/(4.0*(alphaplusnew(i) + alphaminusnew(i)))
                
                		dels = z(bi) - z(i)
                
                		if ( s*delIext .ge. s*dels ) then
                    
                    			if ( (s*dels - alphaminusnew(i)) .gt. 10.0**(-10.0) ) then
                        
                        		alphaplusnew(i) = -2.0*dels - 2.0*s*sqrt(abs(dels**2.0 - dels*alphaminusnew(i)))
                        
                    			else
                        
                        		alphaplusnew(i) = - 2.0*alphaminusnew(i)
                        
                    			endif
                    
                		endif                
                
            		endif    
                
            		if ( abs(alphaminusnew(i)) .gt. 2.0*abs(alphaplusnew(i)) ) then
                
                		s = sign(ones,alphaplusnew(i))
                
                		delIext = -( alphaminusnew(i)**2.0)/(4.0*(alphaplusnew(i) + alphaminusnew(i)))
                
                		dels = z(ai) - z(i)
                
                		if (s*delIext .ge. s*dels) then
                    
                    			if ( (s*dels - alphaplusnew(i)) .gt. 10.0**(-10.0) ) then
                        
                        		alphaminusnew(i) = -2.0*dels - 2.0*s*sqrt(abs(dels**2.0 - dels*alphaplusnew(i)))
                        
                    			else
                        
                        		alphaminusnew(i) = - 2.0*alphaplusnew(i)
                        
                    			endif
                    
                		endif               
                
            		endif

		endif

        	zplus(i) = z(i) + alphaplusnew(i)
        	zminus(i) = z(i) + alphaminusnew(i)

	ENDDO

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif 
    
    		z6(i) = 6.0*(z(i) - 0.5*(zminus(i)+zplus(i)))
    
	ENDDO

	! calculate flux based on velocity direction

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif
    
    	if (u .ge. 0.0) then
        	flux(i) = zplus(bi) - 0.5*(u*dt/dx)*( -zminus(bi)+zplus(bi) - z6(bi)*(1.0-(2.0/3.0)*(u*dt/dx)) )
    	else
        	flux(i) = zminus(i) + 0.5*(abs(u)*dt/dx)*( zplus(i)-zminus(i) + z6(i)*(1.0-(2.0/3.0)*(abs(u)*dt/dx)) )
    	endif
    
	ENDDO

	

END SUBROUTINE ppmcs1d
! =========================================================

! =========================================================
SUBROUTINE PPMCS1D_TLM(flux, fluxd, z, zd, u, ud, dt, dx, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: dx, dt, z(n), u
  DOUBLE PRECISION, INTENT(IN) :: zd(n), ud
  DOUBLE PRECISION, INTENT(OUT) :: flux(n)
  DOUBLE PRECISION, INTENT(OUT) :: fluxd(n)
  DOUBLE PRECISION :: zplus(n), zminus(n), zp(n), zminusnew(n), zplusnew&
& (n), z6(n), d2s, d2sl, d2sr, d2sc, s, d2slim, ones, zeros, d2smin, &
& alphaplus(n), alphaminus(n), alphaminusnew(n), alphaplusnew(n), &
& dsfacem, dsfacep, dsfacemin, dsccm, dsccp, dsp, dsm, deliext, dels, &
& dsccmin
  DOUBLE PRECISION :: zplusd(n), zminusd(n), zpd(n), z6d(n), d2sd, d2sld&
& , d2srd, d2scd, d2slimd, d2smind, alphaplusd(n), alphaminusd(n), &
& alphaminusnewd(n), alphaplusnewd(n), delsd
  INTEGER :: i, ai, bi, ci, di, extreme(n)
  INTRINSIC SIGN
  INTRINSIC ABS
  INTRINSIC MIN
  INTRINSIC MAX
  INTRINSIC SQRT
  REAL :: pwr1
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION :: result1d
  DOUBLE PRECISION :: abs18d
  DOUBLE PRECISION :: x6d
  DOUBLE PRECISION :: z2d
  DOUBLE PRECISION :: y7d
  DOUBLE PRECISION :: max2d
  DOUBLE PRECISION :: x7
  DOUBLE PRECISION :: abs7d
  DOUBLE PRECISION :: x6
  DOUBLE PRECISION :: x5
  DOUBLE PRECISION :: x4
  DOUBLE PRECISION :: x3
  DOUBLE PRECISION :: x2
  DOUBLE PRECISION :: x1
  DOUBLE PRECISION :: abs17d
  DOUBLE PRECISION :: z1d
  DOUBLE PRECISION :: max1d
  DOUBLE PRECISION :: y6d
  DOUBLE PRECISION :: abs3d
  DOUBLE PRECISION :: abs18
  DOUBLE PRECISION :: abs17
  DOUBLE PRECISION :: abs13d
  DOUBLE PRECISION :: abs6d
  DOUBLE PRECISION :: abs16
  DOUBLE PRECISION :: abs15
  DOUBLE PRECISION :: abs14
  DOUBLE PRECISION :: abs13
  DOUBLE PRECISION :: abs12
  DOUBLE PRECISION :: x1d
  DOUBLE PRECISION :: abs11
  DOUBLE PRECISION :: abs10
  DOUBLE PRECISION :: abs9d
  DOUBLE PRECISION :: abs16d
  DOUBLE PRECISION :: z2
  DOUBLE PRECISION :: z1
  DOUBLE PRECISION :: abs9
  DOUBLE PRECISION :: abs8
  DOUBLE PRECISION :: abs7
  DOUBLE PRECISION :: x7d
  DOUBLE PRECISION :: abs6
  DOUBLE PRECISION :: abs5
  DOUBLE PRECISION :: abs4
  DOUBLE PRECISION :: abs3
  DOUBLE PRECISION :: abs2
  DOUBLE PRECISION :: abs1
  DOUBLE PRECISION :: abs0
  DOUBLE PRECISION :: max3d
  DOUBLE PRECISION :: abs12d
  DOUBLE PRECISION :: abs5d
  DOUBLE PRECISION :: max3
  DOUBLE PRECISION :: y7
  DOUBLE PRECISION :: max2
  DOUBLE PRECISION :: abs8d
  DOUBLE PRECISION :: y6
  DOUBLE PRECISION :: max1
  DOUBLE PRECISION :: y5
  DOUBLE PRECISION :: y4
  DOUBLE PRECISION :: y3
  DOUBLE PRECISION :: y2
  DOUBLE PRECISION :: y1
  DOUBLE PRECISION :: y1d
  ones = 1.0
  zeros = 0.0
  zpd = 0.D0
! flux to left 
  DO i=1,n
    ai = i + 1
    bi = i - 1
    ci = i + 2
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (ci .EQ. n + 1) ci = 1
    IF (ci .EQ. n + 2) ci = 2
! Fourth-Order Edge-Reconstruction
    zpd(i) = (7.0*zd(ai)-zd(ci)+7.0*zd(i)-zd(bi))/12.0
    zp(i) = (-z(ci)+7.0*z(ai)+7.0*z(i)-z(bi))/12.0
  END DO
  DO i=1,n
    ai = i + 1
    bi = i - 1
    ci = i + 2
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (ci .EQ. n + 1) ci = 1
    IF (ci .EQ. n + 2) ci = 2
! Limit edges if needed
    IF ((zp(i)-z(i))*(z(ai)-zp(i)) .LT. 0.0) THEN
      d2sd = 3.0*(zd(i)-2.0*zpd(i)+zd(ai))
      d2s = 3.0*(z(i)-2.0*zp(i)+z(ai))
      d2sld = zd(bi) - 2.0*zd(i) + zd(ai)
      d2sl = z(bi) - 2.0*z(i) + z(ai)
      d2srd = zd(i) - 2.0*zd(ai) + zd(ci)
      d2sr = z(i) - 2.0*z(ai) + z(ci)
      s = SIGN(ones, d2s)
      IF (d2sl .GE. 0.) THEN
        abs7d = d2sld
        abs7 = d2sl
      ELSE
        abs7d = -d2sld
        abs7 = -d2sl
      END IF
      x1d = 1.25*s*abs7d
      x1 = 1.25*s*abs7
      IF (d2sr .GE. 0.) THEN
        abs8d = d2srd
        abs8 = d2sr
      ELSE
        abs8d = -d2srd
        abs8 = -d2sr
      END IF
      y1d = 1.25*s*abs8d
      y1 = 1.25*s*abs8
      IF (d2s .GE. 0.) THEN
        abs9d = d2sd
        abs9 = d2s
      ELSE
        abs9d = -d2sd
        abs9 = -d2s
      END IF
      z1d = s*abs9d
      z1 = s*abs9
      IF (x1 .GT. y1) THEN
        IF (y1 .GT. z1) THEN
          d2smind = z1d
          d2smin = z1
        ELSE
          d2smind = y1d
          d2smin = y1
        END IF
      ELSE IF (x1 .GT. z1) THEN
        d2smind = z1d
        d2smin = z1
      ELSE
        d2smind = x1d
        d2smin = x1
      END IF
      IF (d2smin .LT. 0.0) THEN
        max1 = 0.0
        max1d = 0.D0
      ELSE
        max1d = d2smind
        max1 = d2smin
      END IF
      d2slimd = s*max1d
      d2slim = s*max1
      zpd(i) = 0.5*(zd(i)+zd(ai)) - d2slimd/6.0
      zp(i) = 0.5*(z(i)+z(ai)) - 1.0/6.0*d2slim
    END IF
  END DO
! Is it an extremum
  extreme(1:n) = 0
  alphaplusnewd = 0.D0
  alphaplusd = 0.D0
  alphaminusnewd = 0.D0
  alphaminusd = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    ci = i + 2
    di = i - 2
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (ci .EQ. n + 1) ci = 1
    IF (ci .EQ. n + 2) ci = 2
    IF (di .EQ. 0) di = n
    IF (di .EQ. -1) di = n - 1
    alphaplusd(i) = zpd(i) - zd(i)
    alphaplus(i) = zp(i) - z(i)
    alphaminusd(i) = zpd(bi) - zd(i)
    alphaminus(i) = zp(bi) - z(i)
    alphaplusnewd(i) = alphaplusd(i)
    alphaplusnew(i) = alphaplus(i)
    alphaminusnewd(i) = alphaminusd(i)
    alphaminusnew(i) = alphaminus(i)
    IF (alphaplus(i)*alphaminus(i) .GE. 0.0) extreme(i) = 1
    IF (alphaplus(i) .GE. 0.) THEN
      abs0 = alphaplus(i)
    ELSE
      abs0 = -alphaplus(i)
    END IF
    IF (alphaminus(i) .GE. 0.) THEN
      abs10 = alphaminus(i)
    ELSE
      abs10 = -alphaminus(i)
    END IF
    IF (abs0 .GT. 2.0*abs10) THEN
      dsfacem = zp(bi) - zp(di)
      dsfacep = zp(ai) - zp(bi)
      IF (dsfacem .GE. 0.) THEN
        x2 = dsfacem
      ELSE
        x2 = -dsfacem
      END IF
      IF (dsfacep .GE. 0.) THEN
        y2 = dsfacep
      ELSE
        y2 = -dsfacep
      END IF
      IF (x2 .GT. y2) THEN
        dsfacemin = y2
      ELSE
        dsfacemin = x2
      END IF
      dsccm = z(i) - z(bi)
      dsccp = z(ai) - z(i)
      IF (dsccm .GE. 0.) THEN
        x3 = dsccm
      ELSE
        x3 = -dsccm
      END IF
      IF (dsccp .GE. 0.) THEN
        y3 = dsccp
      ELSE
        y3 = -dsccp
      END IF
      IF (x3 .GT. y3) THEN
        dsccmin = y3
      ELSE
        dsccmin = x3
      END IF
      IF (dsfacemin .GE. dsccmin) THEN
        dsp = dsfacep
        dsm = dsfacem
      ELSE
        dsp = dsccp
        dsm = dsccm
      END IF
      IF (dsp*dsm .LE. 0.0) extreme(i) = 1
    END IF
    IF (alphaminus(i) .GE. 0.) THEN
      abs1 = alphaminus(i)
    ELSE
      abs1 = -alphaminus(i)
    END IF
    IF (alphaplus(i) .GE. 0.) THEN
      abs11 = alphaplus(i)
    ELSE
      abs11 = -alphaplus(i)
    END IF
    IF (abs1 .GT. 2.0*abs11) THEN
      dsfacem = zp(bi) - zp(di)
      dsfacep = zp(ai) - zp(bi)
      IF (dsfacem .GE. 0.) THEN
        x4 = dsfacem
      ELSE
        x4 = -dsfacem
      END IF
      IF (dsfacep .GE. 0.) THEN
        y4 = dsfacep
      ELSE
        y4 = -dsfacep
      END IF
      IF (x4 .GT. y4) THEN
        dsfacemin = y4
      ELSE
        dsfacemin = x4
      END IF
      dsccm = z(i) - z(bi)
      dsccp = z(ai) - z(i)
      IF (dsccm .GE. 0.) THEN
        x5 = dsccm
      ELSE
        x5 = -dsccm
      END IF
      IF (dsccp .GE. 0.) THEN
        y5 = dsccp
      ELSE
        y5 = -dsccp
      END IF
      IF (x5 .GT. y5) THEN
        dsccmin = y5
      ELSE
        dsccmin = x5
      END IF
      IF (dsfacemin .GE. dsccmin) THEN
        dsp = dsfacep
        dsm = dsfacem
      ELSE
        dsp = dsccp
        dsm = dsccm
      END IF
      IF (dsp*dsm .LE. 0.0) extreme(i) = 1
    END IF
  END DO
  zminusd = 0.D0
  zplusd = 0.D0
! At Extremum...
  DO i=1,n
    ai = i + 1
    bi = i - 1
    ci = i + 2
    di = i - 2
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (ci .EQ. n + 1) ci = 1
    IF (ci .EQ. n + 2) ci = 2
    IF (di .EQ. 0) di = n
    IF (di .EQ. -1) di = n - 1
    IF (extreme(i) .EQ. 1) THEN
      d2sd = 6.0*(alphaplusd(i)-alphaminusd(i))
      d2s = 6.0*(alphaplus(i)-alphaminus(i))
      d2sld = zd(di) - 2.0*zd(bi) + zd(i)
      d2sl = z(di) - 2.0*z(bi) + z(i)
      d2srd = zd(i) - 2.0*zd(ai) + zd(ci)
      d2sr = z(i) - 2.0*z(ai) + z(ci)
      d2scd = zd(bi) - 2.0*zd(i) + zd(i)
      d2sc = z(bi) - 2.0*z(i) + z(i)
      s = SIGN(ones, d2s)
      IF (d2sl .GE. 0.) THEN
        abs12d = d2sld
        abs12 = d2sl
      ELSE
        abs12d = -d2sld
        abs12 = -d2sl
      END IF
      x6d = 1.25*s*abs12d
      x6 = 1.25*s*abs12
      IF (d2sr .GE. 0.) THEN
        abs13d = d2srd
        abs13 = d2sr
      ELSE
        abs13d = -d2srd
        abs13 = -d2sr
      END IF
      y6d = 1.25*s*abs13d
      y6 = 1.25*s*abs13
      IF (d2sc .GE. 0.) THEN
        abs17d = d2scd
        abs17 = d2sc
      ELSE
        abs17d = -d2scd
        abs17 = -d2sc
      END IF
      x7d = 1.25*s*abs17d
      x7 = 1.25*s*abs17
      IF (d2s .GE. 0.) THEN
        abs18d = d2sd
        abs18 = d2s
      ELSE
        abs18d = -d2sd
        abs18 = -d2s
      END IF
      y7d = s*abs18d
      y7 = s*abs18
      IF (x7 .GT. y7) THEN
        z2d = y7d
        z2 = y7
      ELSE
        z2d = x7d
        z2 = x7
      END IF
      IF (x6 .GT. y6) THEN
        IF (y6 .GT. z2) THEN
          d2smind = z2d
          d2smin = z2
        ELSE
          d2smind = y6d
          d2smin = y6
        END IF
      ELSE IF (x6 .GT. z2) THEN
        d2smind = z2d
        d2smin = z2
      ELSE
        d2smind = x6d
        d2smin = x6
      END IF
      IF (d2smin .LT. 0.0) THEN
        d2slim = 0.0
        d2slimd = 0.D0
      ELSE
        d2slimd = d2smind
        d2slim = d2smin
      END IF
      pwr1 = 10.0**(-10.0)
      IF (d2s .LT. pwr1) THEN
        max2d = 0.D0
        max2 = 10.0**(-10.0)
      ELSE
        max2d = d2sd
        max2 = d2s
      END IF
      alphaplusnewd(i) = ((alphaplusnewd(i)*d2slim+alphaplusnew(i)*&
&       d2slimd)*max2-alphaplusnew(i)*d2slim*max2d)/max2**2
      alphaplusnew(i) = alphaplusnew(i)*d2slim/max2
      pwr1 = 10.0**(-10.0)
      IF (d2s .LT. pwr1) THEN
        max3d = 0.D0
        max3 = 10.0**(-10.0)
      ELSE
        max3d = d2sd
        max3 = d2s
      END IF
      alphaminusnewd(i) = ((alphaminusnewd(i)*d2slim+alphaminusnew(i)*&
&       d2slimd)*max3-alphaminusnew(i)*d2slim*max3d)/max3**2
      alphaminusnew(i) = alphaminusnew(i)*d2slim/max3
    ELSE
      IF (alphaplusnew(i) .GE. 0.) THEN
        abs2 = alphaplusnew(i)
      ELSE
        abs2 = -alphaplusnew(i)
      END IF
      IF (alphaminusnew(i) .GE. 0.) THEN
        abs14 = alphaminusnew(i)
      ELSE
        abs14 = -alphaminusnew(i)
      END IF
      IF (abs2 .GT. 2.0*abs14) THEN
        s = SIGN(ones, alphaminusnew(i))
        deliext = -(alphaplusnew(i)**2.0/(4.0*(alphaplusnew(i)+&
&         alphaminusnew(i))))
        delsd = zd(bi) - zd(i)
        dels = z(bi) - z(i)
        IF (s*deliext .GE. s*dels) THEN
          pwr1 = 10.0**(-10.0)
          IF (s*dels - alphaminusnew(i) .GT. pwr1) THEN
            IF (dels**2.0 - dels*alphaminusnew(i) .GE. 0.) THEN
              abs3d = 2.0*dels*delsd - delsd*alphaminusnew(i) - dels*&
&               alphaminusnewd(i)
              abs3 = dels**2.0 - dels*alphaminusnew(i)
            ELSE
              abs3d = -(2.0*dels*delsd-delsd*alphaminusnew(i)-dels*&
&               alphaminusnewd(i))
              abs3 = -(dels**2.0-dels*alphaminusnew(i))
            END IF
            IF (abs3 .EQ. 0.0) THEN
              result1d = 0.D0
            ELSE
              result1d = abs3d/(2.0*SQRT(abs3))
            END IF
            result1 = SQRT(abs3)
            alphaplusnewd(i) = -(2.0*delsd) - 2.0*s*result1d
            alphaplusnew(i) = -(2.0*dels) - 2.0*s*result1
          ELSE
            alphaplusnewd(i) = -(2.0*alphaminusnewd(i))
            alphaplusnew(i) = -(2.0*alphaminusnew(i))
          END IF
        END IF
      END IF
      IF (alphaminusnew(i) .GE. 0.) THEN
        abs4 = alphaminusnew(i)
      ELSE
        abs4 = -alphaminusnew(i)
      END IF
      IF (alphaplusnew(i) .GE. 0.) THEN
        abs15 = alphaplusnew(i)
      ELSE
        abs15 = -alphaplusnew(i)
      END IF
      IF (abs4 .GT. 2.0*abs15) THEN
        s = SIGN(ones, alphaplusnew(i))
        deliext = -(alphaminusnew(i)**2.0/(4.0*(alphaplusnew(i)+&
&         alphaminusnew(i))))
        delsd = zd(ai) - zd(i)
        dels = z(ai) - z(i)
        IF (s*deliext .GE. s*dels) THEN
          pwr1 = 10.0**(-10.0)
          IF (s*dels - alphaplusnew(i) .GT. pwr1) THEN
            IF (dels**2.0 - dels*alphaplusnew(i) .GE. 0.) THEN
              abs5d = 2.0*dels*delsd - delsd*alphaplusnew(i) - dels*&
&               alphaplusnewd(i)
              abs5 = dels**2.0 - dels*alphaplusnew(i)
            ELSE
              abs5d = -(2.0*dels*delsd-delsd*alphaplusnew(i)-dels*&
&               alphaplusnewd(i))
              abs5 = -(dels**2.0-dels*alphaplusnew(i))
            END IF
            IF (abs5 .EQ. 0.0) THEN
              result1d = 0.D0
            ELSE
              result1d = abs5d/(2.0*SQRT(abs5))
            END IF
            result1 = SQRT(abs5)
            alphaminusnewd(i) = -(2.0*delsd) - 2.0*s*result1d
            alphaminusnew(i) = -(2.0*dels) - 2.0*s*result1
          ELSE
            alphaminusnewd(i) = -(2.0*alphaplusnewd(i))
            alphaminusnew(i) = -(2.0*alphaplusnew(i))
          END IF
        END IF
      END IF
    END IF
    zplusd(i) = zd(i) + alphaplusnewd(i)
    zplus(i) = z(i) + alphaplusnew(i)
    zminusd(i) = zd(i) + alphaminusnewd(i)
    zminus(i) = z(i) + alphaminusnew(i)
  END DO
  z6d = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    z6d(i) = 6.0*(zd(i)-0.5*(zminusd(i)+zplusd(i)))
    z6(i) = 6.0*(z(i)-0.5*(zminus(i)+zplus(i)))
  END DO
  fluxd = 0.D0
! calculate flux based on velocity direction
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (u .GE. 0.0) THEN
      fluxd(i) = zplusd(bi) - 0.5*(dt*ud*(-zminus(bi)+zplus(bi)-z6(bi)*(&
&       1.0-2.0/3.0*(u*dt/dx)))/dx+u*dt*(zplusd(bi)-zminusd(bi)-z6d(bi)*&
&       (1.0-2.0/3.0*(u*dt/dx))+z6(bi)*2.0*dt*ud/(3.0*dx))/dx)
      flux(i) = zplus(bi) - 0.5*(u*dt/dx)*(-zminus(bi)+zplus(bi)-z6(bi)*&
&       (1.0-2.0/3.0*(u*dt/dx)))
    ELSE
      IF (u .GE. 0.) THEN
        abs6d = ud
        abs6 = u
      ELSE
        abs6d = -ud
        abs6 = -u
      END IF
      IF (u .GE. 0.) THEN
        abs16d = ud
        abs16 = u
      ELSE
        abs16d = -ud
        abs16 = -u
      END IF
      fluxd(i) = zminusd(i) + 0.5*(dt*abs6d*(zplus(i)-zminus(i)+z6(i)*(&
&       1.0-2.0/3.0*(abs16*dt/dx)))/dx+abs6*dt*(zplusd(i)-zminusd(i)+z6d&
&       (i)*(1.0-2.0/3.0*(abs16*dt/dx))-z6(i)*2.0*dt*abs16d/(3.0*dx))/dx&
&       )
      flux(i) = zminus(i) + 0.5*(abs6*dt/dx)*(zplus(i)-zminus(i)+z6(i)*(&
&       1.0-2.0/3.0*(abs16*dt/dx)))
    END IF
  END DO
END SUBROUTINE PPMCS1D_TLM
! =========================================================

! =========================================================
SUBROUTINE ppmfromfv(flux,z,u,dt,dx,N,limchoice)

! subroutine for Lin scheme in 1D

IMPLICIT NONE
INTEGER, INTENT(in):: N, limchoice
DOUBLE PRECISION, INTENT(in) :: dx, dt, z(N), u
DOUBLE PRECISION, INTENT(out) :: flux(N)
DOUBLE PRECISION :: al(N), ar(N), a6(N), zp(N), deltaz(N), deltazmax(N), deltazmin(N), &
	minpart, dm(N), da1, da2, a6da, dl, dr, ones
INTEGER :: i, ai, bi, ci


	ones=1.0

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif
    
    		deltaz(i)=0.25*(z(ai)-z(bi))
    
    		deltazmax(i) = max(z(bi),z(i),z(ai))-z(i)
    		deltazmin(i) = z(i)-min(z(bi),z(i),z(ai))

	ENDDO

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif
    
    		minpart=min(abs(deltaz(i)),deltazmin(i),deltazmax(i))
    
    		if (deltaz(i) .ge. 0.0) then
       			dm(i)=abs(minpart)
    		else
        		dm(i)=-abs(minpart)
    		endif
    
	ENDDO


	DO i=1,N
		ci=i+2
		ai=i+1
		bi=i-1
		if (bi==0) then
			bi=N
		endif
		if (ai==N+1) then
			ai=1
		endif
		if (ci==N+1) then
			ci=1
		endif
		if (ci==N+2) then
			ci=2
		endif
	
		al(i) = 0.5*(z(bi)+z(i)) + (dm(bi) - dm(i))/3.0

	ENDDO

	DO i=1,N
		ai=i+1
		bi=i-1
		if (bi==0) then
			bi=N
		endif
		if (ai==N+1) then
			ai=1
		endif

        	ar(i) = al(ai)
		a6(i) = 6.0*(z(i) - 0.5*(al(i)+ar(i)))

	ENDDO

	IF (limchoice .eq. 0) then

		DO i=1, N
     		if (abs(dm(i)) .lt. 0.00000000001) then
         		ar(i) = z(i)
         		al(i) = z(i)
         		a6(i) = 0.0
     		else
         		da1  = ar(i) - al(i)
         		da2  = da1**2
         		a6da = a6(i)*da1
         		if (a6da .lt. -da2) then
            			a6(i) = 3.0*(al(i)-z(i))
            			ar(i) = al(i) - a6(i)
         		elseif(a6da .gt. da2) then
            			a6(i) = 3.0*(ar(i)-z(i))
            			al(i) = ar(i) - a6(i)
         		endif
     		endif
  		ENDDO


	ELSEIF (limchoice .eq. 1) then

		DO i=1, N
           	da1 = dm(i) + dm(i)
            	dl = sign(min(abs(da1),abs(al(i)-z(i))), da1)
            	dr = sign(min(abs(da1),abs(ar(i)-z(i))), da1)
         	ar(i) = z(i) + dr
         	al(i) = z(i) - dl
         	a6(i) = 3.0*(dl-dr)
      		ENDDO

	ENDIF

	! calculate flux based on velocity direction

	DO i=1,N
    		ai=i+1
		bi=i-1

		if (ai==N+1) then
			ai=1
		endif

		if (bi==0) then
			bi=N
		endif
    
    	if (u .ge. 0.0) then
        	flux(i) = ar(bi) - 0.5*(u*dt/dx)*( -al(bi)+ar(bi) - a6(bi)*(1.0-(2.0/3.0)*(u*dt/dx)) )
    	else
        	flux(i) = al(i) + 0.5*(abs(u)*dt/dx)*( ar(i)-al(i) + a6(i)*(1.0-(2.0/3.0)*(abs(u)*dt/dx)) )
    	endif
    
	ENDDO

END SUBROUTINE ppmfromfv
! =========================================================

! =========================================================
SUBROUTINE PPMFROMFV_TLM(flux, fluxd, z, zd, u, ud, dt, dx, n, limchoice)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, limchoice
  DOUBLE PRECISION, INTENT(IN) :: dx, dt, z(n), u
  DOUBLE PRECISION, INTENT(IN) :: zd(n), ud
  DOUBLE PRECISION, INTENT(OUT) :: flux(n)
  DOUBLE PRECISION, INTENT(OUT) :: fluxd(n)
  DOUBLE PRECISION :: al(n), ar(n), a6(n), zp(n), deltaz(n), deltazmax(n&
& ), deltazmin(n), minpart, dm(n), da1, da2, a6da, dl, dr, ones
  DOUBLE PRECISION :: ald(n), ard(n), a6d(n), deltazd(n), deltazmaxd(n)&
& , deltazmind(n), minpartd, dmd(n), da1d, dld, drd
  INTEGER :: i, ai, bi, ci
  INTRINSIC MAX
  INTRINSIC MIN
  INTRINSIC ABS
  INTRINSIC SIGN
  DOUBLE PRECISION :: min3
  DOUBLE PRECISION :: min2
  DOUBLE PRECISION :: min1
  DOUBLE PRECISION :: min1d
  DOUBLE PRECISION :: x3
  DOUBLE PRECISION :: x2
  DOUBLE PRECISION :: x2d
  DOUBLE PRECISION :: x1
  DOUBLE PRECISION :: abs0d
  DOUBLE PRECISION :: max1d
  DOUBLE PRECISION :: abs3d
  DOUBLE PRECISION :: x1d
  DOUBLE PRECISION :: min3d
  DOUBLE PRECISION :: y2d
  DOUBLE PRECISION :: abs3
  DOUBLE PRECISION :: abs2
  DOUBLE PRECISION :: abs2d
  DOUBLE PRECISION :: abs1
  DOUBLE PRECISION :: abs0
  DOUBLE PRECISION :: max1
  DOUBLE PRECISION :: min2d
  DOUBLE PRECISION :: y2
  DOUBLE PRECISION :: x3d
  DOUBLE PRECISION :: y1
  DOUBLE PRECISION :: y1d
  ones = 1.0
  deltazmaxd = 0.D0
  deltazmind = 0.D0
  deltazd = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    deltazd(i) = 0.25*(zd(ai)-zd(bi))
    deltaz(i) = 0.25*(z(ai)-z(bi))
    IF (z(bi) .LT. z(i)) THEN
      IF (z(i) .LT. z(ai)) THEN
        max1d = zd(ai)
        max1 = z(ai)
      ELSE
        max1d = zd(i)
        max1 = z(i)
      END IF
    ELSE IF (z(bi) .LT. z(ai)) THEN
      max1d = zd(ai)
      max1 = z(ai)
    ELSE
      max1d = zd(bi)
      max1 = z(bi)
    END IF
    deltazmaxd(i) = max1d - zd(i)
    deltazmax(i) = max1 - z(i)
    IF (z(bi) .GT. z(i)) THEN
      IF (z(i) .GT. z(ai)) THEN
        min1d = zd(ai)
        min1 = z(ai)
      ELSE
        min1d = zd(i)
        min1 = z(i)
      END IF
    ELSE IF (z(bi) .GT. z(ai)) THEN
      min1d = zd(ai)
      min1 = z(ai)
    ELSE
      min1d = zd(bi)
      min1 = z(bi)
    END IF
    deltazmind(i) = zd(i) - min1d
    deltazmin(i) = z(i) - min1
  END DO
  dmd = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (deltaz(i) .GE. 0.) THEN
      x1d = deltazd(i)
      x1 = deltaz(i)
    ELSE
      x1d = -deltazd(i)
      x1 = -deltaz(i)
    END IF
    IF (x1 .GT. deltazmin(i)) THEN
      IF (deltazmin(i) .GT. deltazmax(i)) THEN
        minpartd = deltazmaxd(i)
        minpart = deltazmax(i)
      ELSE
        minpartd = deltazmind(i)
        minpart = deltazmin(i)
      END IF
    ELSE IF (x1 .GT. deltazmax(i)) THEN
      minpartd = deltazmaxd(i)
      minpart = deltazmax(i)
    ELSE
      minpartd = x1d
      minpart = x1
    END IF
    IF (deltaz(i) .GE. 0.0) THEN
      IF (minpart .GE. 0.) THEN
        dmd(i) = minpartd
        dm(i) = minpart
      ELSE
        dmd(i) = -minpartd
        dm(i) = -minpart
      END IF
    ELSE
      IF (minpart .GE. 0.) THEN
        abs0d = minpartd
        abs0 = minpart
      ELSE
        abs0d = -minpartd
        abs0 = -minpart
      END IF
      dmd(i) = -abs0d
      dm(i) = -abs0
    END IF
  END DO
  ald = 0.D0
  DO i=1,n
    ci = i + 2
    ai = i + 1
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    IF (ai .EQ. n + 1) ai = 1
    IF (ci .EQ. n + 1) ci = 1
    IF (ci .EQ. n + 2) ci = 2
    ald(i) = 0.5*(zd(bi)+zd(i)) + (dmd(bi)-dmd(i))/3.0
    al(i) = 0.5*(z(bi)+z(i)) + (dm(bi)-dm(i))/3.0
  END DO
  ard = 0.D0
  a6d = 0.D0
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (bi .EQ. 0) bi = n
    IF (ai .EQ. n + 1) ai = 1
    ard(i) = ald(ai)
    ar(i) = al(ai)
    a6d(i) = 6.0*(zd(i)-0.5*(ald(i)+ard(i)))
    a6(i) = 6.0*(z(i)-0.5*(al(i)+ar(i)))
  END DO
  IF (limchoice .EQ. 0) THEN
    DO i=1,n
      IF (dm(i) .GE. 0.) THEN
        abs1 = dm(i)
      ELSE
        abs1 = -dm(i)
      END IF
      IF (abs1 .LT. 0.00000000001) THEN
        ard(i) = zd(i)
        ar(i) = z(i)
        ald(i) = zd(i)
        al(i) = z(i)
        a6d(i) = 0.D0
        a6(i) = 0.0
      ELSE
        da1 = ar(i) - al(i)
        da2 = da1**2
        a6da = a6(i)*da1
        IF (a6da .LT. -da2) THEN
          a6d(i) = 3.0*(ald(i)-zd(i))
          a6(i) = 3.0*(al(i)-z(i))
          ard(i) = ald(i) - a6d(i)
          ar(i) = al(i) - a6(i)
        ELSE IF (a6da .GT. da2) THEN
          a6d(i) = 3.0*(ard(i)-zd(i))
          a6(i) = 3.0*(ar(i)-z(i))
          ald(i) = ard(i) - a6d(i)
          al(i) = ar(i) - a6(i)
        END IF
      END IF
    END DO
    fluxd = 0.D0
  ELSE IF (limchoice .EQ. 1) THEN
    DO i=1,n
      da1d = 2*dmd(i)
      da1 = dm(i) + dm(i)
      IF (da1 .GE. 0.) THEN
        x2d = da1d
        x2 = da1
      ELSE
        x2d = -da1d
        x2 = -da1
      END IF
      IF (al(i) - z(i) .GE. 0.) THEN
        y1d = ald(i) - zd(i)
        y1 = al(i) - z(i)
      ELSE
        y1d = -(ald(i)-zd(i))
        y1 = -(al(i)-z(i))
      END IF
      IF (x2 .GT. y1) THEN
        min2d = y1d
        min2 = y1
      ELSE
        min2d = x2d
        min2 = x2
      END IF
      dld = min2d*SIGN(1.d0, min2*da1)
      dl = SIGN(min2, da1)
      IF (da1 .GE. 0.) THEN
        x3d = da1d
        x3 = da1
      ELSE
        x3d = -da1d
        x3 = -da1
      END IF
      IF (ar(i) - z(i) .GE. 0.) THEN
        y2d = ard(i) - zd(i)
        y2 = ar(i) - z(i)
      ELSE
        y2d = -(ard(i)-zd(i))
        y2 = -(ar(i)-z(i))
      END IF
      IF (x3 .GT. y2) THEN
        min3d = y2d
        min3 = y2
      ELSE
        min3d = x3d
        min3 = x3
      END IF
      drd = min3d*SIGN(1.d0, min3*da1)
      dr = SIGN(min3, da1)
      ard(i) = zd(i) + drd
      ar(i) = z(i) + dr
      ald(i) = zd(i) - dld
      al(i) = z(i) - dl
      a6d(i) = 3.0*(dld-drd)
      a6(i) = 3.0*(dl-dr)
    END DO
    fluxd = 0.D0
  ELSE
    fluxd = 0.D0
  END IF
! calculate flux based on velocity direction
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (u .GE. 0.0) THEN
      fluxd(i) = ard(bi) - 0.5*(dt*ud*(-al(bi)+ar(bi)-a6(bi)*(1.0-2.0/&
&       3.0*(u*dt/dx)))/dx+u*dt*(ard(bi)-ald(bi)-a6d(bi)*(1.0-2.0/3.0*(u&
&       *dt/dx))+a6(bi)*2.0*dt*ud/(3.0*dx))/dx)
      flux(i) = ar(bi) - 0.5*(u*dt/dx)*(-al(bi)+ar(bi)-a6(bi)*(1.0-2.0/&
&       3.0*(u*dt/dx)))
    ELSE
      IF (u .GE. 0.) THEN
        abs2d = ud
        abs2 = u
      ELSE
        abs2d = -ud
        abs2 = -u
      END IF
      IF (u .GE. 0.) THEN
        abs3d = ud
        abs3 = u
      ELSE
        abs3d = -ud
        abs3 = -u
      END IF
      fluxd(i) = ald(i) + 0.5*(dt*abs2d*(ar(i)-al(i)+a6(i)*(1.0-2.0/3.0*&
&       (abs3*dt/dx)))/dx+abs2*dt*(ard(i)-ald(i)+a6d(i)*(1.0-2.0/3.0*(&
&       abs3*dt/dx))-a6(i)*2.0*dt*abs3d/(3.0*dx))/dx)
      flux(i) = al(i) + 0.5*(abs2*dt/dx)*(ar(i)-al(i)+a6(i)*(1.0-2.0/3.0&
&       *(abs3*dt/dx)))
    END IF
  END DO
END SUBROUTINE PPMFROMFV_TLM
! =========================================================

! =========================================================
SUBROUTINE simplesemilagrangianschemebs(qnew,q,C,N)

! Do a simple cubic interpolation for SL scheme with BS limiter

IMPLICIT NONE
DOUBLE PRECISION, INTENT(in) :: q(N), C
INTEGER, INTENT(in) :: N
DOUBLE PRECISION, INTENT(out) :: qnew(N)
INTEGER :: i, ii, ai, bi, ci, di
INTEGER :: adder
DOUBLE PRECISION :: alpha, qmax, qmin

		adder = floor(C)
		alpha = C-adder

		! can only go up to a max Courant number of 4

		DO i=1, N

			! Indices depend on size of Courant number

                	ii=i-adder
			ai=i+1-adder
			bi=i-1-adder
			ci=i+2-adder
			di=i-2-adder

			! make sure within range 1-N

			if (ii .eq. N+1) then
				ii=1
			endif
			if (ii .eq. N+2) then
				ii=2
			endif
			if (ii .eq. N+3) then
				ii=3
			endif
			if (ii .eq. N+4) then
				ii=4
			endif
			if (ii .eq. 0) then
				ii=N
			endif
			if (ii .eq. -1) then
				ii=N-1
			endif
			if (ii .eq. -2) then
				ii=N-2
			endif
			if (ii .eq. -3) then
				ii=N-3
			endif
			if (ii .eq. -4) then
				ii=N-4
			endif

			if (ai .eq. N+1) then
				ai=1
			endif
			if (ai .eq. N+2) then
				ai=2
			endif
			if (ai .eq. N+3) then
				ai=3
			endif
			if (ai .eq. N+4) then
				ai=4
			endif
			if (ai .eq. 0) then
				ai=N
			endif
			if (ai .eq. -1) then
				ai=N-1
			endif
			if (ai .eq. -2) then
				ai=N-2
			endif
			if (ai .eq. -3) then
				ai=N-3
			endif
			if (ai .eq. -4) then
				ai=N-4
			endif

			if (bi .eq. N+1) then
				bi=1
			endif
			if (bi .eq. N+2) then
				bi=2
			endif
			if (bi .eq. N+3) then
				bi=3
			endif
			if (bi .eq. N+4) then
				bi=4
			endif
			if (bi .eq. 0) then
				bi=N
			endif
			if (bi .eq. -1) then
				bi=N-1
			endif
			if (bi .eq. -2) then
				bi=N-2
			endif
			if (bi .eq. -3) then
				bi=N-3
			endif
			if (bi .eq. -4) then
				bi=N-4
			endif

			if (ci .eq. N+1) then
				ci=1
			endif
			if (ci .eq. N+2) then
				ci=2
			endif
			if (ci .eq. N+3) then
				ci=3
			endif
			if (ci .eq. N+4) then
				ci=4
			endif
			if (ci .eq. 0) then
				ci=N
			endif
			if (ci .eq. -1) then
				ci=N-1
			endif
			if (ci .eq. -2) then
				ci=N-2
			endif
			if (ci .eq. -3) then
				ci=N-3
			endif
			if (ci .eq. -4) then
				ci=N-4
			endif

			if (di .eq. N+1) then
				di=1
			endif
			if (di .eq. N+2) then
				di=2
			endif
			if (di .eq. N+3) then
				di=3
			endif
			if (di .eq. N+4) then
				di=4
			endif
			if (di .eq. 0) then
				di=N
			endif
			if (di .eq. -1) then
				di=N-1
			endif
			if (di .eq. -2) then
				di=N-2
			endif
			if (di .eq. -3) then
				di=N-3
			endif
			if (di .eq. -4) then
				di=N-4
			endif		

			qnew(i) = -((alpha*(1.0-alpha**2.0))/6.0)*q(di) &
			+ ((alpha*(1.0+alpha)*(2.0-alpha))/2.0)*q(bi) &
			+ (((1.0-alpha**2)*(2.0-alpha))/2.0)*q(ii) &
			- ((alpha*(1.0-alpha)*(2.0-alpha))/6.0)*q(ai)

			! for positive u!!!

			qmax=max(q(ii),q(bi))
			qmin=min(q(ii),q(bi))

			if (qnew(i) .ge. qmax) then
				qnew(i)=qmax
			endif

			if (qnew(i) .le. qmin) then
				qnew(i)=qmin
			endif
                    
		ENDDO

END SUBROUTINE simplesemilagrangianschemebs
! ================================================================

! ================================================================
SUBROUTINE SIMPLESEMILAGRANGIANSCHEMEBS_TLM(qnew, qnewd, q, qd, c, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: q(n), c
  DOUBLE PRECISION, INTENT(IN) :: qd(n)
  DOUBLE PRECISION, INTENT(OUT) :: qnew(n)
  DOUBLE PRECISION, INTENT(OUT) :: qnewd(n)
  INTEGER :: i, ii, ai, bi, ci, di
  INTEGER :: adder
  DOUBLE PRECISION :: alpha, qmax, qmin
  DOUBLE PRECISION :: qmaxd, qmind
  INTRINSIC FLOOR
  INTRINSIC MAX
  INTRINSIC MIN
  adder = FLOOR(c)
  alpha = c - adder
  qnewd = 0.D0
! can only go up to a max Courant number of 4
  DO i=1,n
! Indices depend on size of Courant number
    ii = i - adder
    ai = i + 1 - adder
    bi = i - 1 - adder
    ci = i + 2 - adder
    di = i - 2 - adder
! make sure within range 1-N
    IF (ii .EQ. n + 1) ii = 1
    IF (ii .EQ. n + 2) ii = 2
    IF (ii .EQ. n + 3) ii = 3
    IF (ii .EQ. n + 4) ii = 4
    IF (ii .EQ. 0) ii = n
    IF (ii .EQ. -1) ii = n - 1
    IF (ii .EQ. -2) ii = n - 2
    IF (ii .EQ. -3) ii = n - 3
    IF (ii .EQ. -4) ii = n - 4
    IF (ai .EQ. n + 1) ai = 1
    IF (ai .EQ. n + 2) ai = 2
    IF (ai .EQ. n + 3) ai = 3
    IF (ai .EQ. n + 4) ai = 4
    IF (ai .EQ. 0) ai = n
    IF (ai .EQ. -1) ai = n - 1
    IF (ai .EQ. -2) ai = n - 2
    IF (ai .EQ. -3) ai = n - 3
    IF (ai .EQ. -4) ai = n - 4
    IF (bi .EQ. n + 1) bi = 1
    IF (bi .EQ. n + 2) bi = 2
    IF (bi .EQ. n + 3) bi = 3
    IF (bi .EQ. n + 4) bi = 4
    IF (bi .EQ. 0) bi = n
    IF (bi .EQ. -1) bi = n - 1
    IF (bi .EQ. -2) bi = n - 2
    IF (bi .EQ. -3) bi = n - 3
    IF (bi .EQ. -4) bi = n - 4
    IF (ci .EQ. n + 1) ci = 1
    IF (ci .EQ. n + 2) ci = 2
    IF (ci .EQ. n + 3) ci = 3
    IF (ci .EQ. n + 4) ci = 4
    IF (ci .EQ. 0) ci = n
    IF (ci .EQ. -1) ci = n - 1
    IF (ci .EQ. -2) ci = n - 2
    IF (ci .EQ. -3) ci = n - 3
    IF (ci .EQ. -4) ci = n - 4
    IF (di .EQ. n + 1) di = 1
    IF (di .EQ. n + 2) di = 2
    IF (di .EQ. n + 3) di = 3
    IF (di .EQ. n + 4) di = 4
    IF (di .EQ. 0) di = n
    IF (di .EQ. -1) di = n - 1
    IF (di .EQ. -2) di = n - 2
    IF (di .EQ. -3) di = n - 3
    IF (di .EQ. -4) di = n - 4
    qnewd(i) = alpha*(1.0+alpha)*(2.0-alpha)*qd(bi)/2.0 - alpha*(1.0-&
&     alpha**2.0)*qd(di)/6.0 + (1.0-alpha**2)*(2.0-alpha)*qd(ii)/2.0 - &
&     alpha*(1.0-alpha)*(2.0-alpha)*qd(ai)/6.0
    qnew(i) = -(alpha*(1.0-alpha**2.0)/6.0*q(di)) + alpha*(1.0+alpha)*(&
&     2.0-alpha)/2.0*q(bi) + (1.0-alpha**2)*(2.0-alpha)/2.0*q(ii) - &
&     alpha*(1.0-alpha)*(2.0-alpha)/6.0*q(ai)
    IF (q(ii) .LT. q(bi)) THEN
      qmaxd = qd(bi)
      qmax = q(bi)
    ELSE
      qmaxd = qd(ii)
      qmax = q(ii)
    END IF
    IF (q(ii) .GT. q(bi)) THEN
      qmind = qd(bi)
      qmin = q(bi)
    ELSE
      qmind = qd(ii)
      qmin = q(ii)
    END IF
    IF (qnew(i) .GE. qmax) THEN
      qnewd(i) = qmaxd
      qnew(i) = qmax
    END IF
    IF (qnew(i) .LE. qmin) THEN
      qnewd(i) = qmind
      qnew(i) = qmin
    END IF
  END DO
END SUBROUTINE SIMPLESEMILAGRANGIANSCHEMEBS_TLM
! ================================================================

! ================================================================
SUBROUTINE simplesemilagrangianscheme(qnew,q,C,N)

! Do a simple cubic interpolation for SL scheme

IMPLICIT NONE
DOUBLE PRECISION, INTENT(in) :: q(N), C
INTEGER, INTENT(in) :: N
DOUBLE PRECISION, INTENT(out) :: qnew(N)

INTEGER :: i, ii, ai, bi, ci, di
INTEGER :: adder
DOUBLE PRECISION :: alpha

		adder = floor(C)
		alpha = C-adder

		! can only go up to a max Courant number of 4

		DO i=1, N

			! Indices depend on size of Courant number

                	ii=i-adder
			ai=i+1-adder
			bi=i-1-adder
			ci=i+2-adder
			di=i-2-adder

			! make sure within range 1-N

			if (ii .eq. N+1) then
				ii=1
			endif
			if (ii .eq. N+2) then
				ii=2
			endif
			if (ii .eq. N+3) then
				ii=3
			endif
			if (ii .eq. N+4) then
				ii=4
			endif
			if (ii .eq. 0) then
				ii=N
			endif
			if (ii .eq. -1) then
				ii=N-1
			endif
			if (ii .eq. -2) then
				ii=N-2
			endif
			if (ii .eq. -3) then
				ii=N-3
			endif
			if (ii .eq. -4) then
				ii=N-4
			endif

			if (ai .eq. N+1) then
				ai=1
			endif
			if (ai .eq. N+2) then
				ai=2
			endif
			if (ai .eq. N+3) then
				ai=3
			endif
			if (ai .eq. N+4) then
				ai=4
			endif
			if (ai .eq. 0) then
				ai=N
			endif
			if (ai .eq. -1) then
				ai=N-1
			endif
			if (ai .eq. -2) then
				ai=N-2
			endif
			if (ai .eq. -3) then
				ai=N-3
			endif
			if (ai .eq. -4) then
				ai=N-4
			endif

			if (bi .eq. N+1) then
				bi=1
			endif
			if (bi .eq. N+2) then
				bi=2
			endif
			if (bi .eq. N+3) then
				bi=3
			endif
			if (bi .eq. N+4) then
				bi=4
			endif
			if (bi .eq. 0) then
				bi=N
			endif
			if (bi .eq. -1) then
				bi=N-1
			endif
			if (bi .eq. -2) then
				bi=N-2
			endif
			if (bi .eq. -3) then
				bi=N-3
			endif
			if (bi .eq. -4) then
				bi=N-4
			endif

			if (ci .eq. N+1) then
				ci=1
			endif
			if (ci .eq. N+2) then
				ci=2
			endif
			if (ci .eq. N+3) then
				ci=3
			endif
			if (ci .eq. N+4) then
				ci=4
			endif
			if (ci .eq. 0) then
				ci=N
			endif
			if (ci .eq. -1) then
				ci=N-1
			endif
			if (ci .eq. -2) then
				ci=N-2
			endif
			if (ci .eq. -3) then
				ci=N-3
			endif
			if (ci .eq. -4) then
				ci=N-4
			endif

			if (di .eq. N+1) then
				di=1
			endif
			if (di .eq. N+2) then
				di=2
			endif
			if (di .eq. N+3) then
				di=3
			endif
			if (di .eq. N+4) then
				di=4
			endif
			if (di .eq. 0) then
				di=N
			endif
			if (di .eq. -1) then
				di=N-1
			endif
			if (di .eq. -2) then
				di=N-2
			endif
			if (di .eq. -3) then
				di=N-3
			endif
			if (di .eq. -4) then
				di=N-4
			endif		

			qnew(i) = -((alpha*(1.0-alpha**2.0))/6.0)*q(di) &
			+ ((alpha*(1.0+alpha)*(2.0-alpha))/2.0)*q(bi) &
			+ (((1.0-alpha**2)*(2.0-alpha))/2.0)*q(ii) &
			- ((alpha*(1.0-alpha)*(2.0-alpha))/6.0)*q(ai)
                    
		ENDDO

END SUBROUTINE simplesemilagrangianscheme
! ================================================================

! ================================================================
SUBROUTINE SIMPLESEMILAGRANGIANSCHEME_TLM(qnew, qnewd, q, qd, c, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: q(n), c
  DOUBLE PRECISION, INTENT(IN) :: qd(n)
  DOUBLE PRECISION, INTENT(OUT) :: qnew(n)
  DOUBLE PRECISION, INTENT(OUT) :: qnewd(n)
  INTEGER :: i, ii, ai, bi, ci, di
  INTEGER :: adder
  DOUBLE PRECISION :: alpha
  INTRINSIC FLOOR
  adder = FLOOR(c)
  alpha = c - adder
  qnewd = 0.D0
! can only go up to a max Courant number of 4
  DO i=1,n
! Indices depend on size of Courant number
    ii = i - adder
    ai = i + 1 - adder
    bi = i - 1 - adder
    ci = i + 2 - adder
    di = i - 2 - adder
! make sure within range 1-N
    IF (ii .EQ. n + 1) ii = 1
    IF (ii .EQ. n + 2) ii = 2
    IF (ii .EQ. n + 3) ii = 3
    IF (ii .EQ. n + 4) ii = 4
    IF (ii .EQ. 0) ii = n
    IF (ii .EQ. -1) ii = n - 1
    IF (ii .EQ. -2) ii = n - 2
    IF (ii .EQ. -3) ii = n - 3
    IF (ii .EQ. -4) ii = n - 4
    IF (ai .EQ. n + 1) ai = 1
    IF (ai .EQ. n + 2) ai = 2
    IF (ai .EQ. n + 3) ai = 3
    IF (ai .EQ. n + 4) ai = 4
    IF (ai .EQ. 0) ai = n
    IF (ai .EQ. -1) ai = n - 1
    IF (ai .EQ. -2) ai = n - 2
    IF (ai .EQ. -3) ai = n - 3
    IF (ai .EQ. -4) ai = n - 4
    IF (bi .EQ. n + 1) bi = 1
    IF (bi .EQ. n + 2) bi = 2
    IF (bi .EQ. n + 3) bi = 3
    IF (bi .EQ. n + 4) bi = 4
    IF (bi .EQ. 0) bi = n
    IF (bi .EQ. -1) bi = n - 1
    IF (bi .EQ. -2) bi = n - 2
    IF (bi .EQ. -3) bi = n - 3
    IF (bi .EQ. -4) bi = n - 4
    IF (ci .EQ. n + 1) ci = 1
    IF (ci .EQ. n + 2) ci = 2
    IF (ci .EQ. n + 3) ci = 3
    IF (ci .EQ. n + 4) ci = 4
    IF (ci .EQ. 0) ci = n
    IF (ci .EQ. -1) ci = n - 1
    IF (ci .EQ. -2) ci = n - 2
    IF (ci .EQ. -3) ci = n - 3
    IF (ci .EQ. -4) ci = n - 4
    IF (di .EQ. n + 1) di = 1
    IF (di .EQ. n + 2) di = 2
    IF (di .EQ. n + 3) di = 3
    IF (di .EQ. n + 4) di = 4
    IF (di .EQ. 0) di = n
    IF (di .EQ. -1) di = n - 1
    IF (di .EQ. -2) di = n - 2
    IF (di .EQ. -3) di = n - 3
    IF (di .EQ. -4) di = n - 4
    qnewd(i) = alpha*(1.0+alpha)*(2.0-alpha)*qd(bi)/2.0 - alpha*(1.0-&
&     alpha**2.0)*qd(di)/6.0 + (1.0-alpha**2)*(2.0-alpha)*qd(ii)/2.0 - &
&     alpha*(1.0-alpha)*(2.0-alpha)*qd(ai)/6.0
    qnew(i) = -(alpha*(1.0-alpha**2.0)/6.0*q(di)) + alpha*(1.0+alpha)*(&
&     2.0-alpha)/2.0*q(bi) + (1.0-alpha**2)*(2.0-alpha)/2.0*q(ii) - &
&     alpha*(1.0-alpha)*(2.0-alpha)/6.0*q(ai)
  END DO
END SUBROUTINE SIMPLESEMILAGRANGIANSCHEME_TLM
! ================================================================

! ================================================================
SUBROUTINE sliceschemebs(qnew,q,C,x,dx,N)

! Do the SLICE scheme with bs limiter

IMPLICIT NONE
DOUBLE PRECISION, INTENT(in) :: q(N), C, x(N), dx
INTEGER, INTENT(in) :: N
DOUBLE PRECISION, INTENT(out) :: qnew(N)
DOUBLE PRECISION :: domain, qg(N), xd(N), qmax, qmin
INTEGER :: i, ai, bi

	domain=1.0

	DO i=1,N

		xd(i) = (x(i)) - C*dx

		if (xd(i) .le. 0.0) then
			xd(i)=xd(i)+domain
		endif

	ENDDO

	call slice1d(x,xd,domain,q,qnew,qg,N)
	

	DO i=1,N

		ai=i+1
		bi=i-1
	
		if (ai .eq. N+1) then
			ai=1
		endif
		if (bi .eq. 0) then
			bi=N
		endif

		! For Courant less than 1

		qmax=max(q(i),q(bi))
		qmin=min(q(i),q(bi))

		if (qnew(i) .ge. qmax) then
			qnew(i)=qmax
		endif

		if (qnew(i) .le. qmin) then
			qnew(i)=qmin
		endif
                    
	ENDDO


END SUBROUTINE sliceschemebs
! ================================================================

SUBROUTINE SLICESCHEMEBS_TLM(qnew, qnewd, q, qd, c, x, dx, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: q(n), c, x(n), dx
  DOUBLE PRECISION, INTENT(IN) :: qd(n)
  DOUBLE PRECISION, INTENT(OUT) :: qnew(n)
  DOUBLE PRECISION, INTENT(OUT) :: qnewd(n)
  DOUBLE PRECISION :: domain, qg(n), xd(n), qmax, qmin
  DOUBLE PRECISION :: qgd(n), xdd(n), qmaxd, qmind
  INTEGER :: i, ai, bi
  INTRINSIC MAX
  INTRINSIC MIN
  domain = 1.0
  DO i=1,n
    xdd(i) = 0.D0
    xd(i) = x(i) - c*dx
    IF (xd(i) .LE. 0.0) THEN
      xdd(i) = 0.D0
      xd(i) = xd(i) + domain
    END IF
  END DO
  CALL SLICE1D_D(x, xd, domain, q, qd, qnew, qnewd, qg, qgd, n)
  DO i=1,n
    ai = i + 1
    bi = i - 1
    IF (ai .EQ. n + 1) ai = 1
    IF (bi .EQ. 0) bi = n
    IF (q(i) .LT. q(bi)) THEN
      qmaxd = qd(bi)
      qmax = q(bi)
    ELSE
      qmaxd = qd(i)
      qmax = q(i)
    END IF
    IF (q(i) .GT. q(bi)) THEN
      qmind = qd(bi)
      qmin = q(bi)
    ELSE
      qmind = qd(i)
      qmin = q(i)
    END IF
    IF (qnew(i) .GE. qmax) THEN
      qnewd(i) = qmaxd
      qnew(i) = qmax
    END IF
    IF (qnew(i) .LE. qmin) THEN
      qnewd(i) = qmind
      qnew(i) = qmin
    END IF
  END DO
END SUBROUTINE SLICESCHEMEBS_TLM

! ================================================================
SUBROUTINE slicescheme(qnew,q,C,x,dx,N)

! Do the SLICE scheme

IMPLICIT NONE
DOUBLE PRECISION, INTENT(in) :: q(N), C, x(N), dx
INTEGER, INTENT(in) :: N
DOUBLE PRECISION, INTENT(out) :: qnew(N)
DOUBLE PRECISION :: domain, qg(N), xd(N)
INTEGER :: i

	domain=1.0

	DO i=1,N

		xd(i) = (x(i)) - C*dx

		if (xd(i) .le. 0.0) then
			xd(i)=xd(i)+domain
		endif

	ENDDO

	call slice1d(x,xd,domain,q,qnew,qg,N)


END SUBROUTINE slicescheme
! ================================================================

! ========================================================
SUBROUTINE trisolve(x,a,b,c,r,n)

! To solve the constant coefficient, periodic domain,
! tridiagonal linear system
! Ax = r
! where a is the value below the diagonal of A,
! b is the value on the diagonal of A,
! c is the value above the diagonal of A,
! and r is the known vector right hand side.

IMPLICIT NONE

INTEGER,INTENT(IN) :: n
REAL*8, INTENT(IN) :: a(n), b(n), c(n), r(n)
REAL*8, INTENT(OUT) :: x(n)
INTEGER :: j
REAL*8 :: q(n), s(n), rmx, p


rmx=r(n)

! Forward elimination sweep
q(1) = -c(1)/b(1)
x(1) = r(1)/b(1)
s(1) = -a(1)/b(1)
DO j = 2, n
  p = 1.0/(b(j)+a(j)*q(j-1))
  q(j) = -c(j)*p
  x(j) = (r(j)-a(j)*x(j-1))*p
  s(j) = -a(j)*s(j-1)*p
ENDDO

! Backward pass
q(n) = 0.0
s(n) = 1.0
DO j = n-1, 1, -1
  s(j) = s(j)+q(j)*s(j+1)
  q(j) = x(j)+q(j)*q(j+1)
ENDDO

! Final pass
x(n) = (rmx-c(n)*q(1)-a(n)*q(n-1))/(c(n)*s(1)+a(n)*s(n-1)+b(n))
DO j = 1, n-1
  x(j) = x(n)*s(j)+q(j)
ENDDO


END SUBROUTINE trisolve
! =========================================================

! =========================================================
SUBROUTINE slice1d(xg,xd,domain,q,qnew,qg,n)

! To use the SLICE-1D algorithm to conservatively
! advect a quantity q on a one-dimensional periodic
! domain

IMPLICIT NONE

INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: xg(n), xd(n), domain, q(n)
REAL*8, INTENT(OUT) :: qnew(n), qg(n)
INTEGER :: i, idx(n), im, ip, k, j, length, jj
REAL*8 :: a(n), b(n), c(n), r(n), dx(n), rdx(n), xi(n), &
        part(n), a0(n), a1(n), a2(n), xx, sum


! Grid intervals and reciprocals
DO i = 1, n
  ip = MODULO(i,n)+1
  dx(i) = MODULO(xg(ip) - xg(i), domain)
  rdx(i) = 1.0/dx(i)
ENDDO

! Find indices to departure cells
! and fractions of cells
DO i = 1, n
  IF (xd(i) .ge. xg(n) .or. xd(i) .le. xg(1)) THEN
    k = n
  ELSE
    k = CEILING((xd(i) - xg(1))/dx(1))
    k = MAX(MIN(k,n),1)
    ! Safety check for irregular grid
    DO WHILE(xd(i) .lt. xg(k))
      k = k - 1
    ENDDO
    DO WHILE(xd(i) .gt. xg(k+1))
      k = k + 1
    ENDDO
  ENDIF
  idx(i) = k
  xi(i) = MODULO(xd(i) - xg(k),domain)*rdx(k)
ENDDO

! Set up coefficients for tridiagonal problem
! to determine parabolic spline fit
DO i = 1, n
  im = i-1
  IF (i == 1) im = n
  a(i) = rdx(im)
  b(i) = 2.0*(rdx(im) + rdx(i))
  c(i) = rdx(i)
  r(i) = 3.0*(q(im)*rdx(im) + q(i)*rdx(i))
ENDDO

! Solve tridiagonal problem
! to obtain cell edge values qg
CALL trisolve(qg,a,b,c,r,n)

! Hence find coefficients of parabolas
DO i = 1, n
  ip = i + 1
  IF (i == n) ip = 1
  a0(i) = qg(i)
  a1(i) = -2*qg(i) - qg(ip) + 3*q(i)    ! ZWS coeff / 2
  a2(i) = qg(i) + qg(ip) - 2*q(i)       ! ZWS coeff / 3
ENDDO

! Compute partial integrals for each departure point
! and grid point value of q at departure point
DO i = 1, n
  k = idx(i)
  xx = xi(i)
  part(i) = ((((a2(k)*xx + a1(k))*xx) + a0(k))*xx)*dx(k)
  qg(i) = (3.0*a2(k)*xx + 2.0*a1(k))*xx + a0(k)
ENDDO

! Finally compute integrals between departure points
! and update values of q
DO i = 1, n
  ip = i + 1
  IF (i == n) ip = 1
  sum = part(ip) - part(i)
  length = MODULO(idx(ip) - idx(i), n)
  DO j = 1,length
    jj = MODULO(idx(i) + j - 2,n) + 1
    sum = sum + q(jj)*dx(jj)
  ENDDO
  qnew(i) = sum*rdx(i)
ENDDO



END SUBROUTINE slice1d

! ==========================================================


SUBROUTINE SLICESCHEME_TLM(qnew, qnewd, q, qd, c, x, dx, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: q(n), c, x(n), dx
  DOUBLE PRECISION, INTENT(IN) :: qd(n)
  DOUBLE PRECISION, INTENT(OUT) :: qnew(n)
  DOUBLE PRECISION, INTENT(OUT) :: qnewd(n)
  DOUBLE PRECISION :: domain, qg(n), xd(n)
  DOUBLE PRECISION :: qgd(n), xdd(n)
  INTEGER :: i
  domain = 1.0
  DO i=1,n
    xdd(i) = 0.D0
    xd(i) = x(i) - c*dx
    IF (xd(i) .LE. 0.0) THEN
      xdd(i) = 0.D0
      xd(i) = xd(i) + domain
    END IF
  END DO
  CALL SLICE1D_D(x, xd, domain, q, qd, qnew, qnewd, qg, qgd, n)
END SUBROUTINE SLICESCHEME_TLM

!  Differentiation of slice1d in forward (tangent) mode:
!   variations   of useful results: qnew
!   with respect to varying inputs: q
! =========================================================
! =========================================================
SUBROUTINE SLICE1D_D(xg, xd, domain, q, qd, qnew, qnewd, qg, qgd, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL*8, INTENT(IN) :: xg(n), xd(n), domain, q(n)
  REAL*8, INTENT(IN) :: qd(n)
  REAL*8, INTENT(OUT) :: qnew(n), qg(n)
  REAL*8, INTENT(OUT) :: qnewd(n), qgd(n)
  INTEGER :: i, idx(n), im, ip, k, j, length, jj
  REAL*8 :: a(n), b(n), c(n), r(n), dx(n), rdx(n), xi(n), part(n), a0(n)&
& , a1(n), a2(n), xx, sum
  REAL*8 :: ad(n), bd(n), cd(n), rd(n), dxd(n), rdxd(n), xid(n), partd(n&
& ), a0d(n), a1d(n), a2d(n), sumd
  INTRINSIC MIN
  INTRINSIC MAX
  INTEGER :: x1
! Grid intervals and reciprocals
  DO i=1,n
    ip = MODULO(i, n) + 1
    dxd(i) = 0.0_8
    dx(i) = MODULO(xg(ip) - xg(i), domain)
    rdxd(i) = 0.0_8
    rdx(i) = 1.0/dx(i)
  END DO
! Find indices to departure cells
! and fractions of cells
  DO i=1,n
    IF (xd(i) .GE. xg(n) .OR. xd(i) .LE. xg(1)) THEN
      k = n
    ELSE
      k = CEILING((xd(i)-xg(1))/dx(1))
      IF (k .GT. n) THEN
        x1 = n
      ELSE
        x1 = k
      END IF
      IF (x1 .LT. 1) THEN
        k = 1
      ELSE
        k = x1
      END IF
! Safety check for irregular grid
      DO WHILE (xd(i) .LT. xg(k))
        k = k - 1
      END DO
      DO WHILE (xd(i) .GT. xg(k+1))
        k = k + 1
      END DO
    END IF
    idx(i) = k
    xid(i) = 0.0_8
    xi(i) = MODULO(xd(i)-xg(k), domain)*rdx(k)
  END DO
  rd = 0.0_8
! Set up coefficients for tridiagonal problem
! to determine parabolic spline fit
  DO i=1,n
    im = i - 1
    IF (i .EQ. 1) im = n
    ad(i) = 0.0_8
    a(i) = rdx(im)
    bd(i) = 0.0_8
    b(i) = 2.0*(rdx(im)+rdx(i))
    cd(i) = 0.0_8
    c(i) = rdx(i)
    rd(i) = 3.0*(rdx(im)*qd(im)+rdx(i)*qd(i))
    r(i) = 3.0*(q(im)*rdx(im)+q(i)*rdx(i))
  END DO
! Solve tridiagonal problem
! to obtain cell edge values qg
  CALL TRISOLVE_D(qg, qgd, a, b, c, r, rd, n)
!  CALL trisolve(qg,a,b,c,r,n)
!  CALL trisolve(qgd,a,b,c,rd,n)

  a0d = 0.0_8
  a1d = 0.0_8
  a2d = 0.0_8
! Hence find coefficients of parabolas
  DO i=1,n
    ip = i + 1
    IF (i .EQ. n) ip = 1
    a0d(i) = qgd(i)
    a0(i) = qg(i)
! ZWS coeff / 2
    a1d(i) = 3*qd(i) - qgd(ip) - 2*qgd(i)
    a1(i) = -(2*qg(i)) - qg(ip) + 3*q(i)
! ZWS coeff / 3
    a2d(i) = qgd(i) + qgd(ip) - 2*qd(i)
    a2(i) = qg(i) + qg(ip) - 2*q(i)
  END DO
  partd = 0.0_8
! Compute partial integrals for each departure point
! and grid point value of q at departure point
  DO i=1,n
    k = idx(i)
    xx = xi(i)
    partd(i) = xx*dx(k)*(xx*(xx*a2d(k)+a1d(k))+a0d(k))
    part(i) = ((a2(k)*xx+a1(k))*xx+a0(k))*xx*dx(k)
    qg(i) = (3.0*a2(k)*xx+2.0*a1(k))*xx + a0(k)
  END DO
  qnewd = 0.0_8
! Finally compute integrals between departure points
! and update values of q
  DO i=1,n
    ip = i + 1
    IF (i .EQ. n) ip = 1
    sumd = partd(ip) - partd(i)
    sum = part(ip) - part(i)
    length = MODULO(idx(ip) - idx(i), n)
    DO j=1,length
      jj = MODULO(idx(i) + j - 2, n) + 1
      sumd = sumd + dx(jj)*qd(jj)
      sum = sum + q(jj)*dx(jj)
    END DO
    qnewd(i) = rdx(i)*sumd
    qnew(i) = sum*rdx(i)
  END DO
END SUBROUTINE SLICE1D_D

!  Differentiation of trisolve in forward (tangent) mode:
!   variations   of useful results: x
!   with respect to varying inputs: r
! ================================================================
! ========================================================
SUBROUTINE TRISOLVE_D(x, xd, a, b, c, r, rd, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL*8, INTENT(IN) :: a(n), b(n), c(n), r(n)
  REAL*8, INTENT(IN) :: rd(n)
  REAL*8, INTENT(OUT) :: x(n)
  REAL*8, INTENT(OUT) :: xd(n)
  INTEGER :: j
  REAL*8 :: q(n), s(n), rmx, p
  REAL*8 :: qd(n), sd(n), rmxd
  rmxd = rd(n)
  rmx = r(n)
! Forward elimination sweep
  qd(1) = 0.0_8
  q(1) = -(c(1)/b(1))
  xd = 0.0_8
  xd(1) = rd(1)/b(1)
  x(1) = r(1)/b(1)
  sd(1) = 0.0_8
  s(1) = -(a(1)/b(1))
  DO j=2,n
    p = 1.0/(b(j)+a(j)*q(j-1))
    qd(j) = 0.0_8
    q(j) = -(c(j)*p)
    xd(j) = p*(rd(j)-a(j)*xd(j-1))
    x(j) = (r(j)-a(j)*x(j-1))*p
    sd(j) = 0.0_8
    s(j) = -(a(j)*s(j-1)*p)
  END DO
! Backward pass
  qd(n) = 0.0_8
  q(n) = 0.0
  sd(n) = 0.0_8
  s(n) = 1.0
  qd = 0.0_8
  sd = 0.0_8
  DO j=n-1,1,-1
    sd(j) = sd(j) + qd(j)*s(j+1) + q(j)*sd(j+1)
    s(j) = s(j) + q(j)*s(j+1)
    qd(j) = xd(j) + qd(j)*q(j+1) + q(j)*qd(j+1)
    q(j) = x(j) + q(j)*q(j+1)
  END DO
! Final pass
  xd(n) = ((rmxd-c(n)*qd(1)-a(n)*qd(n-1))*(c(n)*s(1)+a(n)*s(n-1)+b(n))-(&
&   rmx-c(n)*q(1)-a(n)*q(n-1))*(c(n)*sd(1)+a(n)*sd(n-1)))/(c(n)*s(1)+a(n&
&   )*s(n-1)+b(n))**2
  x(n) = (rmx-c(n)*q(1)-a(n)*q(n-1))/(c(n)*s(1)+a(n)*s(n-1)+b(n))
  DO j=1,n-1
    xd(j) = xd(n)*s(j) + x(n)*sd(j) + qd(j)
    x(j) = x(n)*s(j) + q(j)
  END DO
END SUBROUTINE TRISOLVE_D

