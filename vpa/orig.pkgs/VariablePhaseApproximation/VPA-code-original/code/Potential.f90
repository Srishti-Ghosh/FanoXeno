
      double precision FUNCTION POTENT1(X,KP)
!     function returns the value of the potential with a well at X
!     times two and times the reduced mass.
!
      IMPLICIT double precision(A-H,O-Z)
      integer KP        ! kind of potential: 1 - attractive, any other - repulsive
      DATA amass_Ne/36443.98898270d0/   ! atomic Ne mass in electron masses
      DATA RE1/5.8456778d0/    ! position of the potential minimum in Bohr's
      DATA EPS/0.000134502d0/  ! well depth of Ne-Ne potential in Hartree

      R=RE1/X
      !if(X.le.2.5d0)R=RE1/2.5d0
      if(X.le.2.50d0)R=RE1/2.5d0
      RM=amass_Ne/2.d0          ! Reduced mass of Ne-Ne
      R6=R**6
      if(KP.eq.1)then
         POTENT=EPS*(R6-2.d0)*R6   ! L-J potential
      else
         POTENT=EPS*R6            ! repulsive 1/(R**6) potential
      endif

      POTENT1=2.d0*RM*POTENT

      RETURN
      END

      double precision FUNCTION DPOTENT1(X,L,KP)
!     DPOTENT1 returns the derivative of the effective
!      potential with a well.
!
      IMPLICIT double precision(A-H,O-Z)
      integer   L, KP   ! KP is a kind of potential: 1 - attractive, any other - repulsive

      DATA amass_Ne/36443.98898270d0/   ! atomic Ne mass in electron masses
      DATA RE1/5.8456778d0/    ! position of the potential minimum in Bohr's
      DATA EPS/0.000134502d0/  ! well depth of Ne-Ne potential in Hartree

      !if(X.le.2.5d0)then
      if(X.le.2.50d0)then
        DPOTENT1=0.d0
        return
      endif

      RM=amass_Ne/2.d0          ! Reduced mass of Ne-Ne
      R=RE1/X
      R6=R**6
      if(KP.eq.1)then
         DPOTENT=-12.d0*EPS*(R6-1.d0)*R6/X   ! derivative of L-J potential
      else
         DPOTENT=-6.0d0*EPS*R6/X   ! derivative of R^6 repulsive potential.
      endif


      DPOTENT1=DPOTENT-2.d0*dfloat(L*(L+1))/(2.0d0*RM*X**3)


      RETURN
      END



