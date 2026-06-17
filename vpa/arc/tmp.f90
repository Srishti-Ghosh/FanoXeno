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
