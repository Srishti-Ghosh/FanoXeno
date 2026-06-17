      PROGRAM VPA

!-----------------------------------------------------------------------
! Demonstration program for the VPA package.
! This is the version of 7 May 2020.
! Written by A.P. Palov and G.G. Balint-Kurti
! ref: A.P. Palov and G.G. Balint-Kurti, Computer Physics Communications (2020)
!-----------------------------------------------------------------------------
! All output is on unit lout = 6.
!-----------------------------------------------------------------------
      IMPLICIT None
      CHARACTER*5 :: CONTR

      REAL*8 POTENT1, Ecr, Ecr1
      external f1, jac1

      INTEGER   i,iopt,istate,itask,itol,jt,liw,lrw,nerr
      INTEGER   LMIN,Lmax
      INTEGER   L,LMAX0
      INTEGER, DIMENSION(1)  :: NEQ
      INTEGER, DIMENSION(21) :: iwork
      data      neq/1/, jt/1/, nerr/0/, lrw/36/, liw/21/, iopt/0/

      REAL*8, DIMENSION(1) ::   y,RTOL,ATOL
      REAL*8, DIMENSION(36) ::  rwork



      REAL*8, DIMENSION(600) ::     E
      REAL*8, DIMENSION(10000) ::   RW1

      REAL*8 RMU,R1(2000),V1(2000), V
      REAL*8 Y0,EMAX,WVMAX,RMAX,XMAX,EMIN,YEND,XEK,FACTOR,CROSS
      REAL*8 YOLD,RW,T,TOUT,PCROSS
!      REAL*8 amass_Ne,autAng,Pi,autAng2,CM1teV,auteV
! Replaced by the H + Kr collision system.
      REAL*8 amass_H,amass_Kr,autAng,Pi,autAng2,CM1teV,auteV

      integer  KP

      ! WKB
      real*8    Pwkb

!      data      amass_Ne/36443.98898270d0/,autAng/0.529177249d0/
! Replaced by the H + Kr collision system
      data   autAng/0.529177249d0/
      DATA amass_H /1837.36d0/   ! atomic H  mass in electron masses
      DATA amass_Kr/239332.49801983824d0/   ! atomic Xe mass in electron masses (replaces Kr)
      data      Pi/3.1415926535897932384626433832795d0/
      data      CM1teV/1.2398134d-4/
      data      auteV/27.211648d0/
      data      KP/1/    ! selects type of potential , 1 if attractive, any other integer if repulsive

      autAng2=autAng*autAng      ! conversion to Angstrom sq
!      RMU=amass_Ne/2.d0          ! reduced Ne-Ne mass
! Rplaced by the reduced mass of H + K.
      RMU=amass_H*amass_Kr/(amass_H + amass_Kr)

!     open files
!      Open(Unit=5,FILE ='./VPA_input.TXT',STATUS='UNKNOWN')
      Open(Unit=8,FILE ='./Pot_Ne2.TXT',STATUS='UNKNOWN')
!      Open(Unit=6,FILE ='./output.TXT',STATUS='UNKNOWN')
      Open(Unit=7,FILE ='./LPH.TXT',STATUS='UNKNOWN')
!      Open(Unit=9,FILE ='./debug.TXT',STATUS='UNKNOWN')

!     Set some parameters
      Rmax=100.d0

!     Read in control word to determine type of scan requiored.
      read(5,*)CONTR
      write(6,*)' Control word is ',CONTR
      Call flush(6)
      
!      STOP

!     set tolerances
      itol = 1
      rtol(1)=0.d0
!     print out potential -- this step can be skipped
!     Subroutine POTENT1(R) returns the potential in atomic units (a.u.)
!     The internuclear separation, R, is in atomic units or Bohr.

!      write(7,*)'BEFORE POT loop LMAX= ',LMAX
      do i=1,2000
       R1(i)=0.1d0+0.01d0*dfloat(i-1)   ! in Bohr
       V1(i)=POTENT1(R1(i),KP)/2.d0/RMU        ! in Hartree
       V=V1(i)*autev     ! in eV
       write(8,99)R1(i), V
      enddo
      close(8)
 99   format(2f16.7)
!
!      write(7,*)'After POT loop LMAX= ',LMAX

      write(6,*)'Demonstration program for Variable Phase Approach (VPA) package'
      write(6,*)'For the calculation of the phase shift, Y, in an atom-atom scattering problem.'
      write(6,*)'References:'
      write(6,*)'F. Calogero, Variable phase approach to potential scattering '
      write(6,*)'                            (Academic Press, New York, 1967).'
      write(6,*)'A. Ronveaux, Phase Equation in Quantum Mechanics '
      write(6,*)'               (American J. Phys. 37, 135 (1969).'
      write(6,*)'G. G. Balint-Kurti and A. P. Palov, Theory of Molecular Collision, '
      write(6,*)'Royal Society of Chemistry. Theoretical and Computational Chemistry Series,2015.'


      write(6,*)' The mathematical Problem is to solve the Variable Phase Equations for the phase shift, Y:'
      write(6,*)' dY/dR + V(R)/k*(RJ(kR)*cos(Y)-RY(kR)*sin(Y))**2 = 0, '
      write(6,*)' Y(0) = 0'


      write(9,*)' Number of first order differential equations to be solved by subroutine DLSODA; neq = ',neq


      Y0=1.d-10                     ! initial phase shift
      Yend=1.d-04                   ! final phase shift

      IF(CONTR.eq.'ESCAN')THEN

!     This section performs an Energy scan at fixed value of L

       write(9,*)'Cycle over 600 values of relative kinetic energy '
       write(9,*)'Minimum Kinetic Energy/cm^-1 =  ',Emin
       write(9,*)'Maximum Kinetic Energy/cm^-1 =  ',Emin

!       write(9,*)'Kinetic Energy/cm^-1, L, Phase shift/pi, Partial Cross section'
!       write(7,*)'Kinetic Energy/cm^-1, L, Phase shift/pi, Partial Cross section'

       write(6,*)' Performing a scan over energies at fixed value of orbital angular '
       write(6,*)' momentum quantum number L'

       read(5,*) L
       write(6,*)'L ',L

       Lmax0=L

!     read in Emin and Emax in cm-1.
       read(5,*)Emin,Emax
       write(6,*)'Emin= ',Emin,' Emax= ',Emax

       Emin=Emin*CM1teV             ! in eV
       Emax=Emax*CM1teV             ! in eV

       WVmax=dsqrt(2.d0*RMU*Emax/auteV)   ! maximal wave vector in a.u.
       Xmax=WVmax*Rmax

!       Subroutine InitPoint determines the starting value of Xek*R (nondimensional),
!       for all possible values of the orbital angular momentum L.
!       It may also change LMAX0.

       call InitPoint(Y0,Lmax0,Xmax,RW1)


       IF(Lmax0.lt.L) then
        write(6,*)'input L value is too large'
        STOP
       END IF

!     RW1(L) is the starting Xek*R value for orbital angular momentum L.

       Ecr1= Ecr(RMU,Rmax,KP)           ! Ecritical in a.u.

       CALL EGRID(Emin,Emax,E)       ! E in a.u.

!     Do loop over Energy
       ENERGY: Do I=1,600,5
        Xek=dsqrt(2.d0*RMU*E(I))    ! wave vector in a.u.
        Factor=4.d0*Pi*autAng2/Xek/Xek

        if(E(i).lt.10.d0*Ecr1)then
          atol(1)=1.0d-14        ! possible shape resonance
        else
          atol(1) = 1.0d-11      ! area of no resonances or repulsive potential
        endif
        write (*,111) itol,rtol,atol
        write(9,*)' itol =',itol, ' rtol =',rtol(1),' atol =',atol(1)
        cross=0.d0
        Yold=-1.d0
        if(L.eq.0)then
          RW=Y0
        else
          RW=RW1(L)
        end if
        t=RW/Xek    ! starting R value
        y(1)=-Y0    ! starting phase shift value

        itask = 1
        istate = 1
        tout = Rmax
        call dlsoda(f1,KP,L,Xek,neq,y,t,tout,itol,rtol,atol,  &
     &      itask,istate,iopt,rwork,lrw,iwork,liw,jac1,jt)

        pcross=Factor*dfloat(2*L+1)*dsin(Y(1))**2
        cross=cross+pcross

        Write(7,778) E(I)*autev/1.24d-04, L , Y(1)/Pi, pcross   ! in 1/cm
 !       write(9,777) E(I)*autev, L, Y(1)/Pi, pcross             ! in eV
        Yold=Y(1)
       enddo ENERGY                      ! end do over energies

!Now read in control word for scan over L at fixed energy.

      ELSE IF(CONTR.eq.'LSCAN')THEN
        write(6,*)'Performing scan over quantum orbital angular momentum numbers L'
        read(5,*) Lmin,Lmax
        write(6,*)'Lmin= ',Lmin,' Lmax= ',Lmax

        read(5,*)E(1)   ! read in Energy in cm-1
        write(6,*)'Kinetic Energy/cm-1  =',E(1)
        E(1)=E(1)*CM1teV             ! convert to eV
        E(1)=E(1)/autev              ! convert to au
        Emax=E(1)

        Xek=dsqrt(2.d0*RMU*E(1))
        Factor=4.d0*Pi*autAng2/Xek/Xek

        WVmax=dsqrt(2.d0*RMU*Emax)   ! maximal wave vector
        Xmax=WVmax*Rmax

!       Subroutine InitPoint determines the starting value of Xek*R (nondimensional),
!       for all possible values of the orbital angular momentum L.
!       It may also change LMAX0.

        Lmax0=Lmax

        call InitPoint(Y0,Lmax0,Xmax,RW1)

        IF(Lmax.lt.Lmax0) then
           write(6,*)'input LMAX value is too large'
           STOP
        END IF


!     RW1(L) is the starting Xek*R value for quantum orbital angular number L.

        Ecr1= Ecr(RMU,Rmax,KP)           ! Ecritical in a.u.

        if(E(1).lt.10.d0*Ecr1)then
          atol(1)=1.0d-14        ! possible shape resonance
        else
          atol(1) = 1.0d-11      ! area of no resonances or repulsive potential
        endif
        write (6,111) itol,rtol,atol
!        write(9,*)' itol =',itol, ' rtol =',rtol(1),' atol =',atol(1)
 111    format(/' itol =',i3,//' rtol =',d10.1,//' atol =',d10.1//)

        cross=0.d0
        Yold=-1.d0

!     Cycle over orbital angular momentum L

        ORBITAL: DO L=LMIN,LMAX
          IF(L.GT.LMAX)EXIT

          if(L.eq.0)then
            RW=Y0
          else
            RW=RW1(L)
          end if
          t=RW/Xek
          y(1)=-Y0
          Pwkb=0.d0
          if(E(1).gt.10.d0*Ecr1)then
            ! WKB calculations
             call WKB(KP,Rmu,Xek,L,0.1d0,20.d0,Pwkb)
          endif
          itask = 1
          istate = 1
          tout = Rmax
          call dlsoda(f1,KP,L,Xek,neq,y,t,tout,itol,rtol,atol,  &
     &      itask,istate,iopt,rwork,lrw,iwork,liw,jac1,jt)

          if(Ecr1.gt.0.d0) then  ! potential with a well
            if(dabs(Y(1)).le.Yend.AND.Y(1).GT.0.D0)go to 22
            if(dsign(1.d0,Y(1)).le.0.d0.AND.Yold.gt.0.d0)go to 22
          else                   ! pure repulsive potential
            if(dabs(Y(1)).le.Yend)go to 22
          endif

          pcross=Factor*dfloat(2*L+1)*dsin(Y(1))**2
          cross=cross+pcross

          Write(7,778) E(1)*autev/1.24d-04, L, Y(1)/Pi, pcross, Pwkb/Pi     ! in 1/cm
!          write(9,777) E(1)*autev, L, Y(1)/Pi, pcross,Pwkb/Pi               ! in eV
          Yold=Y(1)
        enddo ORBITAL       ! end do over orbital quantum numbers L
 22     continue

      ELSE IF((CONTR.NE.'LSCAN').AND.(CONTR.NE.'ESCAN'))THEN
               write(6,*)' No valid control word input'
               STOP
      ENDIF

      Close(7)
!      Close(9)

 777  format("E(eV)=",  f11.5,2x,"L=",i5," VPA= ",e13.7," PCS= ",e12.4," WKB= ",e13.7)
 778  format("E(1/cm)=",f11.5,2x,"L=",i5," VPA= ",e13.7," PCS= ",e12.4," WKB= ",e13.7)

      stop
      end


