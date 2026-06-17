      subroutine jac1 (KP,L,WK,neq, t, y, pd, nrowpd)
!   the subroutine calculates the Jacobian of the phase shift
!   equation at arbitrary orbital angular quantum number L
!
      integer           neq, nrowpd, KP
      double precision  t, y, pd
      dimension         y(neq), pd(nrowpd,neq)
      double precision  WK,X,Fun1,Fun2,RJ,RY,RJP,RYP,VR,DC,DS
      integer           L, IFAIL
      double precision  POTENT1
      external          POTENT1
!
      X=WK*t           ! calculate K*R
      VR=POTENT1(t,KP)
      call RICBES(X,L,RJ,RY,RJP,RYP,IFAIL)
      DC=dcos(Y(1))
      DS=dsin(Y(1))
      Fun1= DC*RJ-DS*RY
      Fun2= DC*RY+DS*RJ
      pd(1,1) = 2.d0*VR/WK*Fun1*Fun2

      return
      end subroutine

      subroutine f1(KP,L,WK,neq,R,Y,YP)
!   The subroutine calculates the derivative of phase shift at arbitrary
!   orbital angular quantum number L
!
      integer   neq, IFAIL
      double precision    WK,R,Y(neq),YP(neq),VR
      double precision    X, XIN, X1, A, Fun, B, C, RJ,RY,RJP,RYP
      integer   L, I, KP
      double precision    POTENT1
      external  POTENT1
!
      X=WK*R                ! calculate K*R
      XIN=1.D0/X            ! inverse X
      X1=X+Y(1)             ! calculate a sum of (K*R + phase shift)
      A=dsin(X1)
      Fun=A
      IF(L.GT.0)THEN
        B=A*XIN-dcos(X1)
        Fun=B
        IF (L.GT.1) THEN
          IF (L.lt.X)THEN       ! forward recurrence method
            DO I=1,L-1
             C=dfloat(2*I+1)*B*XIN-A
             A=B
             B=C
            END DO
            Fun=C
          ELSE
            call RICBES(X,L,RJ,RY,RJP,RYP, IFAIL )
            Fun= dcos(Y(1))*RJ-dsin(Y(1))*RY
          END IF       ! if L<X
        END IF         ! if L>1
      END IF           ! if L>0

      VR=POTENT1(R,KP)
      YP(1)=-VR/WK*Fun**2

      RETURN
      END subroutine

      SUBROUTINE InitPoint(Y0,Lmax,Xmax,RW)
!     The subroutine determines the starting value of RW
!     RW is R*k, i.e. the distance times the linear momentum in atomic units,
!     hbar=1 in atomic units, so linear momentum = k*hbar=k.
!     The subroutine also examines LMAX.
!     LMAX may be changed so that RMAX, as set in main program, does not lie too far inside the
!     Classically forbidden region for the maximum value of L.

      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(10000) :: RW
      INTEGER :: LMAX,L,IFAIL
      DOUBLE PRECISION :: Y0,RMIN,RL,R,R1,R2,RLmax,RMAX,RJ,RJP,RY,RYP
      DOUBLE PRECISION :: Y,XMAX

      IF(Lmax.EQ.0) RETURN
      RW(:)=0.d0
      Open(Unit=77,FILE ='RW.TXT',STATUS='UNKNOWN')
      L=1
      RW(L)=(3.d0*Y0)**(1.d0/3.d0)
!     write(77,*) L,RW(L),Y0
      do L=2,Lmax
!       write(6,*)"lmax1=",l, xmax
        Rmin=RW(L-1)
        RL=dfloat(L)
        RLmax=dfloat(2*L+1)
        Rmax=(Y0*RLmax)**(1.d0/RLmax)*RL**(2.d0*RL/RLmax)
        R1=Rmin
        R2=Rmax
10      continue
        R=(R1+R2)/2.0D0
        call RICBES(R,L,RJ,RY,RJP,RYP, IFAIL )
        Y= dabs(RJ/RY)
        IF(Y-Y0.LE.0.0D0)THEN
          R1=R
        ELSE
          R2=R
        END IF
        IF((R2-R1).GT.1.0D-10)GO TO 10
        RW(L)=R
!       write(77,*) L,RW(L),Y,RJ,RY
        if(RW(L).gt.Xmax)goto 2
      end do
!       Current Lmax value is OK
      return
!     close(77)
 2    continue
      Lmax=L-1

!       This is new value of Lmax, such that RMAX=k*R is deep in
!       classically forbidden region defined by Xmax-- but not too much so.

      RETURN
      END subroutine

      SUBROUTINE EGRID(Emin,Emax,E)
!   Input minimum and maximum kinetic energies in eV
!   Output E(1)...E(600) in atomic units (a.u.)
      IMPLICIT double precision   (A-H,O-Z)
      DIMENSION         E(600)
      data      auteV/27.211648d0/

    !Emin and Emax are given in eV, E is given in a.u.
    !enddo
      E=0.d0
      DE=(Emax-Emin)/600.d0
      do i=1,600
        E(i)=Emin+DE*dfloat(i-1)
        E(i)=E(i)/auteV
      enddo

      return
      end subroutine

      SUBROUTINE EGRID2(Emin,Emax,nestep,E)
!   Input minimum and maximum kinetic energies in eV
!   Output E(1)...E(nestep+1) in atomic units (a.u.)
      IMPLICIT double precision   (A-H,O-Z)
      DIMENSION         E(10000)
      data      auteV/27.211648d0/

    !Emin and Emax are given in eV, E is given in a.u.
    !enddo
      E=0.d0
      DE=(Emax-Emin)/nestep
      do i=1,nestep+1
        E(i)=Emin+DE*dfloat(i-1)
        E(i)=E(i)/auteV
      enddo

      return
      end subroutine

      DOUBLE PRECISION FUNCTION Ecr(RMU,Rmax, KP)
!     The function Ecr finds critical values of the  orbital angular quantum number L
!     and of the relative kinetic energy E. Only the critical value of the energy
!     is passed back to the main program. Lcr is the value of the orbital angular
!     quantum number above which the effective potential no longer has a minimum.
!     The effective potential is the real potential plus the centrifugal barrier.
!     Ecr is the critical value of the relative kinetic energy such that there are
!     no shape resonances at higher energies. Shape resonances occure when the
!     kinetic energy coincides with a quasi-bound state trapped inside effective
!     potential.
!     KP gives a kind of potential: 1 - attractive, any other - repulsive

      IMPLICIT double precision  (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      data      NR/3000/,LLMax/10000/,autev/27.21d0/
      data      CM1teV/1.2398134d-4/
      external  DPOTENT1, POTENT1

!     Examine derivative of potential, does potential have a minimum?
      DR=Rmax/dfloat(NR)
      R0=DR
      L=0
      DV0=DPOTENT1(R0,L,KP)
      IF(DV0.gt.0.0d0)THEN
        write(6,*)' R0= ',R0,' L= ',L
        write(6,*)' DV0 is positive in Ecr'
        call flush(6)
        STOP
      END IF
      do i=1,NR
        R0=dfloat(i)*DR
        if(DPOTENT1(R0,L,KP).gt.0.d0)goto 30
      enddo
!      The potential is entirely repulsive -- no minimum.
      write(6,*)' The potential is entirely repulsive'
      Ecr=0.0d0
      Return

30    CONTINUE
!     The potential does have a minimum!
!     find L at which a well disappears
      do L=1,LLmax
        do i=1,NR
          R0=dfloat(i-1)*DR
          if(DPOTENT1(R0,L,KP).gt.0.d0)goto 10
        enddo
        goto 11
  10    continue
      enddo
  11  continue
      Lcr=L-1

!     write(6,*)' Lcr= ',Lcr

!     find critical energy such that there are no more shape resonances at higher energies.

      Veff=-1.d0
      do i=1,NR
        R0=dfloat(i-1)*DR
        Veff0=DPOTENT1(R0,Lcr,KP)
        if(Veff0.lt.0.d0.and.Veff.gt.0.d0)goto 15
        Veff=Veff0
      enddo
  15  continue
      Ecr=POTENT1(R0,KP)/2.d0/RMU+dfloat(Lcr*(Lcr+1))/(2.d0*RMU*R0*R0)
  ! in a.u.
!     write(6,*) "Ecr(1/cm)=",Ecr*autev/CM1teV,Lcr,R0

      return
      end function

      subroutine angledcs(ak2,Lmin,Lmax,ym,  &
    &                     scatang,angdcs,angmin,angmax,nangstep)
! Differential cross section following Blatt Biedenharn's formula.
! J. M. Blatt and L. C. Biedenharn, Rev. mod. Phys. 24, 258 (1952).
!
      implicit real *8 (a-h,o-z)
      dimension ym(0:100), scatang(0:360), angdcs(0:360), BL(0:200)
      complex *16 wk0
      data pi/3.1415926535897932384626433832795d0/
!
! Prepare list of scattering angles in units of radian.
      do 10 i = 0, nangstep
      scatang(i) = angmin + (angmax-angmin)/nangstep*i  ! in deg.
      scatang(i) = scatang(i)*(pi/180)                ! in rad.
   10 continue   
!
! Prepare B(ll) facter in BB formula.
      do 50 ll = 0, 2*Lmax 
        wk1 = 0.0d0
          do 40 l1 = Lmin, Lmax
            lwk = abs( l1 - ll )
            l2min = max( Lmin, lwk ) 
            l2max = min ( l1 + ll, Lmax)
!     write(6,*) "Note3", lmin,lmax,l2min,l2max,ll,l1      
            do 30 l2 = l2min, l2max
              wk1 = wk1 + (2*l1+1)*(2*l2+1)*cf(l1,l2,ll,0,0,0)**2 &
     &              *sin(ym(l1))*sin(ym(l2))*cos(ym(l1)-ym(l2))
   30       continue
   40     continue
        BL(ll) = wk1
   50 continue
!
! Start scattering angle loop.
      do 80 i = 0, nangstep
        theta = scatang(i)
!        
! Start BB formula loop        
        wk2 = 0.0d0
          do 70 ll = 0, 2*Lmax
            wk2 = wk2 + BL(ll) * Pnx(ll,cos(theta))
   70     continue
        wk2 = ( 1 / ak2 ) * wk2
        angdcs(i) = wk2
! End of BB formula loop.
!        
! Start |f(theta)|^2 formla loop for comparison.        
!       wk0 = 0
!         do 75 l = lmin,lmax
!           wk0 = wk0 + (2.0d0*l+1.0d0) &
!  &             *(1.0d0-exp(2*(0.0d0,1.0d0)*ym(l)))*Pnx(l,cos(theta))
!  75     continue
!       wk3 = ( 1 / ak2 / 4 ) * abs(wk0)**2
! End of |f(theta)|^2 formula loop.
!
! Check equality of wk2 and wk3
!        write(6,*) wk2, wk3, wk2/wk3,wk3/wk2 , (wk2-wk3)/(wk2+wk3)
! Check has been performed for several combinations of l, E, angle range.
! The equality of both the formulae has been confirmed.
! The quantity |(wk2-wk3)/(wk2+wk3)| is typically 1.0x10^(-15) or less.
!
   80 continue     
      return
      end subroutine
      
      subroutine wparam (Lmin,Lmax,lr,ym,  &
    &                     scatang,angwpm,angbgd,angmin,angmax,nangstep)
! Angle dependent resonance profile parameter W(theta).
! F. Koike, J. Phys. B2, 741 (1977).                              .
!
      implicit real *8 (a-h,o-z)
      dimension ym(0:100), scatang(0:360), angbgd(0:360), ZL(0:200)
      complex *16 wk0,wk1,wk2
      complex *16, dimension(0:360) :: angwpm 
      complex *16, dimension(0:200) :: BL
      data pi/3.1415926535897932384626433832795d0/
!
! Prepare list of scattering angles in units of radian.
      do 10 i = 0, nangstep
      scatang(i) = angmin + (angmax-angmin)/nangstep*i  ! in deg.
      scatang(i) = scatang(i)*(pi/180)                ! in rad.
   10 continue   
!
! Prepare angle independent part of w-parameter.
!
      wk0 = exp(2*(0.0d0,1.0d0)*ym(lr))
!     wk1 = 0
      do 50 ll = 0, 2*Lmax 
       wk1 = 0
       do 40 l = Lmin, Lmax
!       write(6,*) l, ym(l)
         if ( l .eq. lr ) then
           wk1 = wk1 + (2*l+1)*(2*lr+1)*cf(l,lr,ll,0,0,0)**2 
                          else
           wk1 = wk1 + (2*l+1)*(2*lr+1)*cf(l,lr,ll,0,0,0)**2*(1-exp(-2*(0.0d0,1.0d0)*ym(l)))
         end if
   40  continue
      BL(ll) = wk0 * wk1                  
   50 continue
!
! Prepare angle independent part of background factor.
!
      do 55 ll = 0, 2*Lmax
        wk6 = (2*lr+1)**2*cf(lr,lr,ll,0,0,0)**2 
        do 54 l = Lmin, Lmax
          if ( l .ne. lr ) then
            wk6 = wk6 + 2*(2*l+1)*(2*lr+1)*cf(l,lr,ll,0,0,0)**2*sin(ym(l))**2
            do 53 lp = Lmin, Lmax
              if ( lp .ne. lr ) then
                wk6 = wk6 + 2*(2*l+1)*(2*lp+1)*cf(l,lp,ll,0,0,0)**2  &
   &                      * sin(ym(l))*sin(ym(lp))*cos(ym(l)-ym(lp))
              end if
  53        continue
          end if
  54    continue
        ZL(ll) = wk6
  55  continue
!  
! Start scattering angle loop.
      do 80 i = 0, nangstep
        theta = scatang(i)
! For w-parameter        
! Start L sum: ll = 0 to 2*Lmax
        wk2 = 0.0d0
        do 70 ll = 0, 2*Lmax
          wk2 = wk2 + BL(ll) * Pnx(ll,cos(theta))
   70   continue
        angwpm(i) = wk2
!       write(6,*) wk2, angwpm(i)
! End of L sum. 
! For background factor
! Start L sum: ll = 0 to 2*Lmax
        wk7 = 0.0d0
        do 75 ll = 0, 2*Lmax
          wk7 = wk7 + ZL(ll) * Pnx(ll,cos(theta))
   75   continue
        angbgd(i) = wk7   
   80 continue     
      return
      end subroutine
