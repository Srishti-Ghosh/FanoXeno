
      double precision FUNCTION POTENT1(X,KP)
!     function returns the value of the potential with a well at X
!     times two and times the reduced mass.
!
! Modofied by Srishti, 17 June 2026 to adjust the calculation of H + Xe collision system.
!
      IMPLICIT double precision(A-H,O-Z)
      integer KP        ! kind of potential: 1 - attractive, any other - repulsive
!      DATA amass_Ne/36443.98898270d0/   ! atomic Ne mass in electron masses
!      DATA RE1/5.8456778d0/    ! position of the potential minimum in Bohr's
!      DATA EPS/0.000134502d0/  ! well depth of Ne-Ne potential in Hartree
! Replaced by the H + Xe data
      DATA amass_H /1837.1525886135853d0/   ! atomic H mass in electron masses
      DATA amass_Xe/239332.49801983824d0/   ! atomic Xe mass in electron masses
      DATA RE1/7.219d0/    ! position of the potential minimum in Bohr (a0)
      DATA EPS/0.00026018520100363736d0/  ! well depth (7.08 meV) in Hartree
      data rs/8.957d0/ ! Switching radius (8.957 a0 = 4.74 Angstrom)
      data c6/42.298469824178895d0/ ! C6 in Hartree * a0^6
      data c8/918.3655611696183d0/ ! C8 in Hartree * a0^8
      data c10/27554.641767306115d0/ ! C10 in Hartree * a0^10

      R=RE1/X
      !if(X.le.2.5d0)R=RE1/2.5d0
      if(X.le.2.50d0)R=RE1/2.5d0
!      RM=amass_Ne/2.d0          ! Reduced mass of Ne-Ne
! Replaced by the reduced mass of H + Xe.
      RM =amass_H*amass_Xe/(amass_H + amass_Xe)
      R6=R**6
! Added the long range part potential
      if ( x .le. rs ) then
        if(KP.eq.1)then
        POTENT=EPS*(R6-2.d0)*R6   ! L-J potential
        else
        POTENT=EPS*R6            ! repulsive 1/(R**6) potential
        endif
                          else
        if(KP.eq.1)then
        POTENT= - c6/x**6 - c8/x**8 - c10/x**10 ! Long range part
        else
        POTENT=EPS*R6            ! repulsive 1/(R**6) potential
        endif
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

!      DATA amass_Ne/36443.98898270d0/   ! atomic Ne mass in electron masses
!      DATA RE1/5.8456778d0/    ! position of the potential minimum in Bohr's
!      DATA EPS/0.000134502d0/  ! well depth of Ne-Ne potential in Hartree
! Replaced by the H + Xe data
      DATA amass_H /1837.1525886135853d0/   ! atomic H mass in electron masses
      DATA amass_Xe/239332.49801983824d0/   ! atomic Xe mass in electron masses
      DATA RE1/7.219d0/    ! position of the potential minimum in Bohr (a0)
      DATA EPS/0.00026018520100363736d0/  ! well depth (7.08 meV) in Hartree
      data rs/8.957d0/ ! Switching radius (8.957 a0 = 4.74 Angstrom)
      data c6/42.298469824178895d0/ ! C6 in Hartree * a0^6
      data c8/918.3655611696183d0/ ! C8 in Hartree * a0^8
      data c10/27554.641767306115d0/ ! C10 in Hartree * a0^10

      !if(X.le.2.5d0)then
      if(X.le.2.50d0)then
        DPOTENT1=0.d0
        return
      endif

!      RM=amass_Ne/2.d0          ! Reduced mass of Ne-Ne
! Replaced by the reduced mass of H + Xe.
      RM =amass_H*amass_Xe/(amass_H + amass_Xe)
      R=RE1/X
      R6=R**6
      if ( x .le. rs ) then
        if(KP.eq.1)then
        DPOTENT=-12.d0*EPS*(R6-1.d0)*R6/X   ! derivative of L-J potential
        else
        DPOTENT=-6.0d0*EPS*R6/X   ! derivative of R^6 repulsive potential.
        endif
                        else
        if(KP.eq.1)then
        DPOTENT= 6.0d0*c6/x**7 + 8.0d0*c8/x**9 + 10.0d0*c10/x**11
                          ! derivative of the long range potentiall
        else
        DPOTENT=-6.0d0*EPS*R6/X   ! derivative of R^6 repulsive potential.
        endif
      endif

      DPOTENT1=DPOTENT-2.d0*dfloat(L*(L+1))/(2.0d0*RM*X**3)


      RETURN
      END
