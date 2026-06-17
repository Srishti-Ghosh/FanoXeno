      subroutine WKB(KP,Rmu,Xek,L,ax,bx,PhSh)
      ! input
      integer  KP,L
      real*8   Rmu,Xek,ax,bx
      real*8   tol/1.d-14/
      ! output
      real*8   PhSh
      ! inner variables
      real*8   RAD, RA(9001), Vra(9001)
      real*8   xx, x1, x2, drr, D1, D2, Xek2
      real*8   B, B2, B3, A , XL12, XL122, R, R1, R2, Y, xx1
	  Real*8   Pi/3.1415926535897932384626433832795d0/
      real*8   potent1
      external potent1

	  XL12=dfloat(L)+0.5d0
	  XL122=XL12**2
	  Xek2=Xek**2

      R1=ax
      R2=bx
10    R=(R1+R2)/2.0D0
      Y=potent1(R,KP)+XL122/R/R
      IF(Xek2-Y.LE.0.0D0)THEN
        R1=R
      ELSE
        R2=R
      END IF
      IF((R2-R1).GT.1.0D-10)GO TO 10
      xx=R

      drr=(bx-xx)/9000.d0
      do i=1,9001
        RAD   = xx+drr*dfloat(i-1)
        RA(i) = RAD
        VrA(i)= potent1(RAD,KP)/2.d0/Rmu
      enddo
      B=0.d0
      do i=1,9000
	    x1=RA(i)
		x2=RA(i+1)
		D1=Xek2-XL122/x1/x1-2.d0*Rmu*VrA(i)
		if(D1.lt.0.d0) D1=0.d0
		D2=Xek2-XL122/x2/x2-2.d0*Rmu*VrA(i+1)
		if(D2.lt.0.d0) D2=0.d0
		B=B+0.5d0*(dsqrt(D1)+dsqrt(D2))
	  enddo
	  B=B*drr

	  B3=Xek2*bx*bx-XL122
	  if(B3.gt.0.d0)then
	    B2=dsqrt(B3)
	    A=B2+XL12*(datan(XL12/B2)-0.5d0*Pi)
	    !write(*,*)"l=",L,B,B2,A
	    PhSh=B-A
	  else
	    PhSh=0.d0
	  endif

	  return
	  end
