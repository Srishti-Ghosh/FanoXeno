*DECK DLSODA
      SUBROUTINE DLSODA(F,KP,Lfs,Xek,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,
     1            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, JT)
      EXTERNAL F, JAC
      INTEGER KP,Lfs,NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK, LIW, JT
      DOUBLE PRECISION Xek, Y, T, TOUT, RTOL, ATOL, RWORK
      DIMENSION NEQ(*), Y(*), RTOL(*), ATOL(*), RWORK(LRW), IWORK(LIW)
C-----------------------------------------------------------------------
C This is the 12 November 2003 version of
C DLSODA: Livermore Solver for Ordinary Differential Equations, with
C         Automatic method switching for stiff and nonstiff problems.
C
C This version is in double precision.
C
C DLSODA solves the initial value problem for stiff or nonstiff
C systems of first order ODEs,
C     dy/dt = f(t,y) ,  or, in component form,
C     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
C
C This a variant version of the DLSODE package.
C It switches automatically between stiff and nonstiff methods.
C This means that the user does not have to determine whether the
C problem is stiff or not, and the solver will automatically choose the
C appropriate method.  It always starts with the nonstiff method.
C
C Authors:       Alan C. Hindmarsh
C                Center for Applied Scientific Computing, L-561
C                Lawrence Livermore National Laboratory
C                Livermore, CA 94551
C and
C                Linda R. Petzold
C                Univ. of California at Santa Barbara
C                Dept. of Computer Science
C                Santa Barbara, CA 93106
C
C References:
C 1.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
C     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
C     North-Holland, Amsterdam, 1983, pp. 55-64.
C 2.  Linda R. Petzold, Automatic Selection of Methods for Solving
C     Stiff and Nonstiff Systems of Ordinary Differential Equations,
C     Siam J. Sci. Stat. Comput. 4 (1983), pp. 136-148.
C-----------------------------------------------------------------------
C Summary of Usage.
C
C Communication between the user and the DLSODA package, for normal
C situations, is summarized here.  This summary describes only a subset
C of the full set of options available.  See the full description for
C details, including alternative treatment of the Jacobian matrix,
C optional inputs and outputs, nonstandard options, and
C instructions for special situations.  See also the example
C problem (with program and output) following this summary.
C
C A. First provide a subroutine of the form:
C               SUBROUTINE F (NEQ, T, Y, YDOT)
C               DOUBLE PRECISION T, Y(*), YDOT(*)
C which supplies the vector function f by loading YDOT(i) with f(i).
C
C B. Write a main program which calls Subroutine DLSODA once for
C each point at which answers are desired.  This should also provide
C for possible use of logical unit 6 for output of error messages
C by DLSODA.  On the first call to DLSODA, supply arguments as follows:
C F      = name of subroutine for right-hand side vector f.
C          This name must be declared External in calling program.
C NEQ    = number of first order ODEs.
C Y      = array of initial values, of length NEQ.
C T      = the initial value of the independent variable.
C TOUT   = first point where output is desired (.ne. T).
C ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
C RTOL   = relative tolerance parameter (scalar).
C ATOL   = absolute tolerance parameter (scalar or array).
C          the estimated local error in y(i) will be controlled so as
C          to be less than
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
C             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
C          Thus the local error test passes if, in each component,
C          either the absolute error is less than ATOL (or ATOL(i)),
C          or the relative error is less than RTOL.
C          Use RTOL = 0.0 for pure absolute error control, and
C          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
C          control.  Caution: actual (global) errors may exceed these
C          local tolerances, so choose them conservatively.
C ITASK  = 1 for normal computation of output values of y at t = TOUT.
C ISTATE = integer flag (input and output).  Set ISTATE = 1.
C IOPT   = 0 to indicate no optional inputs used.
C RWORK  = real work array of length at least:
C             22 + NEQ * MAX(16, NEQ + 9).
C          See also Paragraph E below.
C LRW    = declared length of RWORK (in user's dimension).
C IWORK  = integer work array of length at least  20 + NEQ.
C LIW    = declared length of IWORK (in user's dimension).
C JAC    = name of subroutine for Jacobian matrix.
C          Use a dummy name.  See also Paragraph E below.
C JT     = Jacobian type indicator.  Set JT = 2.
C          See also Paragraph E below.
C Note that the main program must declare arrays Y, RWORK, IWORK,
C and possibly ATOL.
C
C C. The output from the first call (or any call) is:
C      Y = array of computed values of y(t) vector.
C      T = corresponding value of independent variable (normally TOUT).
C ISTATE = 2  if DLSODA was successful, negative otherwise.
C          -1 means excess work done on this call (perhaps wrong JT).
C          -2 means excess accuracy requested (tolerances too small).
C          -3 means illegal input detected (see printed message).
C          -4 means repeated error test failures (check all inputs).
C          -5 means repeated convergence failures (perhaps bad Jacobian
C             supplied or wrong choice of JT or tolerances).
C          -6 means error weight became zero during problem. (Solution
C             component i vanished, and ATOL or ATOL(i) = 0.)
C          -7 means work space insufficient to finish (see messages).
C
C D. To continue the integration after a successful return, simply
C reset TOUT and call DLSODA again.  No other parameters need be reset.
C
C E. Note: If and when DLSODA regards the problem as stiff, and
C switches methods accordingly, it must make use of the NEQ by NEQ
C Jacobian matrix, J = df/dy.  For the sake of simplicity, the
C inputs to DLSODA recommended in Paragraph B above cause DLSODA to
C treat J as a full matrix, and to approximate it internally by
C difference quotients.  Alternatively, J can be treated as a band
C matrix (with great potential reduction in the size of the RWORK
C array).  Also, in either the full or banded case, the user can supply
C J in closed form, with a routine whose name is passed as the JAC
C argument.  These alternatives are described in the paragraphs on
C RWORK, JAC, and JT in the full description of the call sequence below.
C
C-----------------------------------------------------------------------
C Example Problem.
C
C The following is a simple example problem, with the coding
C needed for its solution by DLSODA.  The problem is from chemical
C kinetics, and consists of the following three rate equations:
C     dy1/dt = -.04*y1 + 1.e4*y2*y3
C     dy2/dt = .04*y1 - 1.e4*y2*y3 - 3.e7*y2**2
C     dy3/dt = 3.e7*y2**2
C on the interval from t = 0.0 to t = 4.e10, with initial conditions
C y1 = 1.0, y2 = y3 = 0.  The problem is stiff.
C
C The following coding solves this problem with DLSODA,
C printing results at t = .4, 4., ..., 4.e10.  It uses
C ITOL = 2 and ATOL much smaller for y2 than y1 or y3 because
C y2 has much smaller values.
C At the end of the run, statistical quantities of interest are
C printed (see optional outputs in the full description below).
C
C     EXTERNAL FEX
C     DOUBLE PRECISION ATOL, RTOL, RWORK, T, TOUT, Y
C     DIMENSION Y(3), ATOL(3), RWORK(70), IWORK(23)
C     NEQ = 3
C     Y(1) = 1.
C     Y(2) = 0.
C     Y(3) = 0.
C     T = 0.
C     TOUT = .4
C     ITOL = 2
C     RTOL = 1.D-4
C     ATOL(1) = 1.D-6
C     ATOL(2) = 1.D-10
C     ATOL(3) = 1.D-6
C     ITASK = 1
C     ISTATE = 1
C     IOPT = 0
C     LRW = 70
C     LIW = 23
C     JT = 2
C     DO 40 IOUT = 1,12
C       CALL DLSODA(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,
C    1     IOPT,RWORK,LRW,IWORK,LIW,JDUM,JT)
C       WRITE(6,20)T,Y(1),Y(2),Y(3)
C 20    FORMAT(' At t =',D12.4,'   Y =',3D14.6)
C       IF (ISTATE .LT. 0) GO TO 80
C 40    TOUT = TOUT*10.
C     WRITE(6,60)IWORK(11),IWORK(12),IWORK(13),IWORK(19),RWORK(15)
C 60  FORMAT(/' No. steps =',I4,'  No. f-s =',I4,'  No. J-s =',I4/
C    1   ' Method last used =',I2,'   Last switch was at t =',D12.4)
C     STOP
C 80  WRITE(6,90)ISTATE
C 90  FORMAT(///' Error halt.. ISTATE =',I3)
C     STOP
C     END
C
C     SUBROUTINE FEX (NEQ, T, Y, YDOT)
C     DOUBLE PRECISION T, Y, YDOT
C     DIMENSION Y(3), YDOT(3)
C     YDOT(1) = -.04*Y(1) + 1.D4*Y(2)*Y(3)
C     YDOT(3) = 3.D7*Y(2)*Y(2)
C     YDOT(2) = -YDOT(1) - YDOT(3)
C     RETURN
C     END
C
C The output of this program (on a CDC-7600 in single precision)
C is as follows:
C
C   At t =  4.0000e-01   y =  9.851712e-01  3.386380e-05  1.479493e-02
C   At t =  4.0000e+00   Y =  9.055333e-01  2.240655e-05  9.444430e-02
C   At t =  4.0000e+01   Y =  7.158403e-01  9.186334e-06  2.841505e-01
C   At t =  4.0000e+02   Y =  4.505250e-01  3.222964e-06  5.494717e-01
C   At t =  4.0000e+03   Y =  1.831975e-01  8.941774e-07  8.168016e-01
C   At t =  4.0000e+04   Y =  3.898730e-02  1.621940e-07  9.610125e-01
C   At t =  4.0000e+05   Y =  4.936363e-03  1.984221e-08  9.950636e-01
C   At t =  4.0000e+06   Y =  5.161831e-04  2.065786e-09  9.994838e-01
C   At t =  4.0000e+07   Y =  5.179817e-05  2.072032e-10  9.999482e-01
C   At t =  4.0000e+08   Y =  5.283401e-06  2.113371e-11  9.999947e-01
C   At t =  4.0000e+09   Y =  4.659031e-07  1.863613e-12  9.999995e-01
C   At t =  4.0000e+10   Y =  1.404280e-08  5.617126e-14  1.000000e+00
C
C   No. steps = 361  No. f-s = 693  No. J-s =  64
C   Method last used = 2   Last switch was at t =  6.0092e-03
C-----------------------------------------------------------------------
C Full description of user interface to DLSODA.
C
C The user interface to DLSODA consists of the following parts.
C
C 1.   The call sequence to Subroutine DLSODA, which is a driver
C      routine for the solver.  This includes descriptions of both
C      the call sequence arguments and of user-supplied routines.
C      following these descriptions is a description of
C      optional inputs available through the call sequence, and then
C      a description of optional outputs (in the work arrays).
C
C 2.   Descriptions of other routines in the DLSODA package that may be
C      (optionally) called by the user.  These provide the ability to
C      alter error message handling, save and restore the internal
C      Common, and obtain specified derivatives of the solution y(t).
C
C 3.   Descriptions of Common blocks to be declared in overlay
C      or similar environments, or to be saved when doing an interrupt
C      of the problem and continued solution later.
C
C 4.   Description of a subroutine in the DLSODA package,
C      which the user may replace with his/her own version, if desired.
C      this relates to the measurement of errors.
C
C-----------------------------------------------------------------------
C Part 1.  Call Sequence.
C
C The call sequence parameters used for input only are
C     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, JT,
C and those used for both input and output are
C     Y, T, ISTATE.
C The work arrays RWORK and IWORK are also used for conditional and
C optional inputs and optional outputs.  (The term output here refers
C to the return from Subroutine DLSODA to the user's calling program.)
C
C The legality of input parameters will be thoroughly checked on the
C initial call for the problem, but not checked thereafter unless a
C change in input parameters is flagged by ISTATE = 3 on input.
C
C The descriptions of the call arguments are as follows.
C
C F      = the name of the user-supplied subroutine defining the
C          ODE system.  The system must be put in the first-order
C          form dy/dt = f(t,y), where f is a vector-valued function
C          of the scalar t and the vector y.  Subroutine F is to
C          compute the function f.  It is to have the form
C               SUBROUTINE F (NEQ, T, Y, YDOT)
C               DOUBLE PRECISION T, Y(*), YDOT(*)
C          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
C          is output.  Y and YDOT are arrays of length NEQ.
C          Subroutine F should not alter Y(1),...,Y(NEQ).
C          F must be declared External in the calling program.
C
C          Subroutine F may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in F) and/or Y has length exceeding NEQ(1).
C          See the descriptions of NEQ and Y below.
C
C          If quantities computed in the F routine are needed
C          externally to DLSODA, an extra call to F should be made
C          for this purpose, for consistent and accurate results.
C          If only the derivative dy/dt is needed, use DINTDY instead.
C
C NEQ    = the size of the ODE system (number of first order
C          ordinary differential equations).  Used only for input.
C          NEQ may be decreased, but not increased, during the problem.
C          If NEQ is decreased (with ISTATE = 3 on input), the
C          remaining components of Y should be left undisturbed, if
C          these are to be accessed in F and/or JAC.
C
C          Normally, NEQ is a scalar, and it is generally referred to
C          as a scalar in this user interface description.  However,
C          NEQ may be an array, with NEQ(1) set to the system size.
C          (The DLSODA package accesses only NEQ(1).)  In either case,
C          this parameter is passed as the NEQ argument in all calls
C          to F and JAC.  Hence, if it is an array, locations
C          NEQ(2),... may be used to store other integer data and pass
C          it to F and/or JAC.  Subroutines F and/or JAC must include
C          NEQ in a Dimension statement in that case.
C
C Y      = a real array for the vector of dependent variables, of
C          length NEQ or more.  Used for both input and output on the
C          first call (ISTATE = 1), and only for output on other calls.
C          On the first call, Y must contain the vector of initial
C          values.  On output, Y contains the computed solution vector,
C          evaluated at T.  If desired, the Y array may be used
C          for other purposes between calls to the solver.
C
C          This array is passed as the Y argument in all calls to
C          F and JAC.  Hence its length may exceed NEQ, and locations
C          Y(NEQ+1),... may be used to store other real data and
C          pass it to F and/or JAC.  (The DLSODA package accesses only
C          Y(1),...,Y(NEQ).)
C
C T      = the independent variable.  On input, T is used only on the
C          first call, as the initial point of the integration.
C          on output, after each call, T is the value at which a
C          computed solution Y is evaluated (usually the same as TOUT).
C          on an error return, T is the farthest point reached.
C
C TOUT   = the next value of t at which a computed solution is desired.
C          Used only for input.
C
C          When starting the problem (ISTATE = 1), TOUT may be equal
C          to T for one call, then should .ne. T for the next call.
C          For the initial t, an input value of TOUT .ne. T is used
C          in order to determine the direction of the integration
C          (i.e. the algebraic sign of the step sizes) and the rough
C          scale of the problem.  Integration in either direction
C          (forward or backward in t) is permitted.
C
C          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
C          the first call (i.e. the first call with TOUT .ne. T).
C          Otherwise, TOUT is required on every call.
C
C          If ITASK = 1, 3, or 4, the values of TOUT need not be
C          monotone, but a value of TOUT which backs up is limited
C          to the current internal T interval, whose endpoints are
C          TCUR - HU and TCUR (see optional outputs, below, for
C          TCUR and HU).
C
C ITOL   = an indicator for the type of error control.  See
C          description below under ATOL.  Used only for input.
C
C RTOL   = a relative error tolerance parameter, either a scalar or
C          an array of length NEQ.  See description below under ATOL.
C          Input only.
C
C ATOL   = an absolute error tolerance parameter, either a scalar or
C          an array of length NEQ.  Input only.
C
C             The input parameters ITOL, RTOL, and ATOL determine
C          the error control performed by the solver.  The solver will
C          control the vector E = (E(i)) of estimated local errors
C          in y, according to an inequality of the form
C                      max-norm of ( E(i)/EWT(i) )   .le.   1,
C          where EWT = (EWT(i)) is a vector of positive error weights.
C          The values of RTOL and ATOL should all be non-negative.
C          The following table gives the types (scalar/array) of
C          RTOL and ATOL, and the corresponding form of EWT(i).
C
C             ITOL    RTOL       ATOL          EWT(i)
C              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
C              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
C              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
C              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
C
C          When either of these parameters is a scalar, it need not
C          be dimensioned in the user's calling program.
C
C          If none of the above choices (with ITOL, RTOL, and ATOL
C          fixed throughout the problem) is suitable, more general
C          error controls can be obtained by substituting a
C          user-supplied routine for the setting of EWT.
C          See Part 4 below.
C
C          If global errors are to be estimated by making a repeated
C          run on the same problem with smaller tolerances, then all
C          components of RTOL and ATOL (i.e. of EWT) should be scaled
C          down uniformly.
C
C ITASK  = an index specifying the task to be performed.
C          Input only.  ITASK has the following values and meanings.
C          1  means normal computation of output values of y(t) at
C             t = TOUT (by overshooting and interpolating).
C          2  means take one step only and return.
C          3  means stop at the first internal mesh point at or
C             beyond t = TOUT and return.
C          4  means normal computation of output values of y(t) at
C             t = TOUT but without overshooting t = TCRIT.
C             TCRIT must be input as RWORK(1).  TCRIT may be equal to
C             or beyond TOUT, but not behind it in the direction of
C             integration.  This option is useful if the problem
C             has a singularity at or beyond t = TCRIT.
C          5  means take one step, without passing TCRIT, and return.
C             TCRIT must be input as RWORK(1).
C
C          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
C          (within roundoff), it will return T = TCRIT (exactly) to
C          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
C          in which case answers at t = TOUT are returned first).
C
C ISTATE = an index used for input and output to specify the
C          the state of the calculation.
C
C          On input, the values of ISTATE are as follows.
C          1  means this is the first call for the problem
C             (initializations will be done).  See note below.
C          2  means this is not the first call, and the calculation
C             is to continue normally, with no change in any input
C             parameters except possibly TOUT and ITASK.
C             (If ITOL, RTOL, and/or ATOL are changed between calls
C             with ISTATE = 2, the new values will be used but not
C             tested for legality.)
C          3  means this is not the first call, and the
C             calculation is to continue normally, but with
C             a change in input parameters other than
C             TOUT and ITASK.  Changes are allowed in
C             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, JT, ML, MU,
C             and any optional inputs except H0, MXORDN, and MXORDS.
C             (See IWORK description for ML and MU.)
C          Note:  A preliminary call with TOUT = T is not counted
C          as a first call here, as no initialization or checking of
C          input is done.  (Such a call is sometimes useful for the
C          purpose of outputting the initial conditions.)
C          Thus the first call for which TOUT .ne. T requires
C          ISTATE = 1 on input.
C
C          On output, ISTATE has the following values and meanings.
C           1  means nothing was done; TOUT = T and ISTATE = 1 on input.
C           2  means the integration was performed successfully.
C          -1  means an excessive amount of work (more than MXSTEP
C              steps) was done on this call, before completing the
C              requested task, but the integration was otherwise
C              successful as far as T.  (MXSTEP is an optional input
C              and is normally 500.)  To continue, the user may
C              simply reset ISTATE to a value .gt. 1 and call again
C              (the excess work step counter will be reset to 0).
C              In addition, the user may increase MXSTEP to avoid
C              this error return (see below on optional inputs).
C          -2  means too much accuracy was requested for the precision
C              of the machine being used.  This was detected before
C              completing the requested task, but the integration
C              was successful as far as T.  To continue, the tolerance
C              parameters must be reset, and ISTATE must be set
C              to 3.  The optional output TOLSF may be used for this
C              purpose.  (Note: If this condition is detected before
C              taking any steps, then an illegal input return
C              (ISTATE = -3) occurs instead.)
C          -3  means illegal input was detected, before taking any
C              integration steps.  See written message for details.
C              Note:  If the solver detects an infinite loop of calls
C              to the solver with illegal input, it will cause
C              the run to stop.
C          -4  means there were repeated error test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              The problem may have a singularity, or the input
C              may be inappropriate.
C          -5  means there were repeated convergence test failures on
C              one attempted step, before completing the requested
C              task, but the integration was successful as far as T.
C              This may be caused by an inaccurate Jacobian matrix,
C              if one is being used.
C          -6  means EWT(i) became zero for some i during the
C              integration.  Pure relative error control (ATOL(i)=0.0)
C              was requested on a variable which has now vanished.
C              The integration was successful as far as T.
C          -7  means the length of RWORK and/or IWORK was too small to
C              proceed, but the integration was successful as far as T.
C              This happens when DLSODA chooses to switch methods
C              but LRW and/or LIW is too small for the new method.
C
C          Note:  Since the normal output value of ISTATE is 2,
C          it does not need to be reset for normal continuation.
C          Also, since a negative input value of ISTATE will be
C          regarded as illegal, a negative output value requires the
C          user to change it, and possibly other inputs, before
C          calling the solver again.
C
C IOPT   = an integer flag to specify whether or not any optional
C          inputs are being used on this call.  Input only.
C          The optional inputs are listed separately below.
C          IOPT = 0 means no optional inputs are being used.
C                   default values will be used in all cases.
C          IOPT = 1 means one or more optional inputs are being used.
C
C RWORK  = a real array (double precision) for work space, and (in the
C          first 20 words) for conditional and optional inputs and
C          optional outputs.
C          As DLSODA switches automatically between stiff and nonstiff
C          methods, the required length of RWORK can change during the
C          problem.  Thus the RWORK array passed to DLSODA can either
C          have a static (fixed) length large enough for both methods,
C          or have a dynamic (changing) length altered by the calling
C          program in response to output from DLSODA.
C
C                       --- Fixed Length Case ---
C          If the RWORK length is to be fixed, it should be at least
C               MAX (LRN, LRS),
C          where LRN and LRS are the RWORK lengths required when the
C          current method is nonstiff or stiff, respectively.
C
C          The separate RWORK length requirements LRN and LRS are
C          as follows:
C          IF NEQ is constant and the maximum method orders have
C          their default values, then
C             LRN = 20 + 16*NEQ,
C             LRS = 22 + 9*NEQ + NEQ**2           if JT = 1 or 2,
C             LRS = 22 + 10*NEQ + (2*ML+MU)*NEQ   if JT = 4 or 5.
C          Under any other conditions, LRN and LRS are given by:
C             LRN = 20 + NYH*(MXORDN+1) + 3*NEQ,
C             LRS = 20 + NYH*(MXORDS+1) + 3*NEQ + LMAT,
C          where
C             NYH    = the initial value of NEQ,
C             MXORDN = 12, unless a smaller value is given as an
C                      optional input,
C             MXORDS = 5, unless a smaller value is given as an
C                      optional input,
C             LMAT   = length of matrix work space:
C             LMAT   = NEQ**2 + 2              if JT = 1 or 2,
C             LMAT   = (2*ML + MU + 1)*NEQ + 2 if JT = 4 or 5.
C
C                       --- Dynamic Length Case ---
C          If the length of RWORK is to be dynamic, then it should
C          be at least LRN or LRS, as defined above, depending on the
C          current method.  Initially, it must be at least LRN (since
C          DLSODA starts with the nonstiff method).  On any return
C          from DLSODA, the optional output MCUR indicates the current
C          method.  If MCUR differs from the value it had on the
C          previous return, or if there has only been one call to
C          DLSODA and MCUR is now 2, then DLSODA has switched
C          methods during the last call, and the length of RWORK
C          should be reset (to LRN if MCUR = 1, or to LRS if
C          MCUR = 2).  (An increase in the RWORK length is required
C          if DLSODA returned ISTATE = -7, but not otherwise.)
C          After resetting the length, call DLSODA with ISTATE = 3
C          to signal that change.
C
C LRW    = the length of the array RWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C IWORK  = an integer array for work space.
C          As DLSODA switches automatically between stiff and nonstiff
C          methods, the required length of IWORK can change during
C          problem, between
C             LIS = 20 + NEQ   and   LIN = 20,
C          respectively.  Thus the IWORK array passed to DLSODA can
C          either have a fixed length of at least 20 + NEQ, or have a
C          dynamic length of at least LIN or LIS, depending on the
C          current method.  The comments on dynamic length under
C          RWORK above apply here.  Initially, this length need
C          only be at least LIN = 20.
C
C          The first few words of IWORK are used for conditional and
C          optional inputs and optional outputs.
C
C          The following 2 words in IWORK are conditional inputs:
C            IWORK(1) = ML     these are the lower and upper
C            IWORK(2) = MU     half-bandwidths, respectively, of the
C                       banded Jacobian, excluding the main diagonal.
C                       The band is defined by the matrix locations
C                       (i,j) with i-ML .le. j .le. i+MU.  ML and MU
C                       must satisfy  0 .le.  ML,MU  .le. NEQ-1.
C                       These are required if JT is 4 or 5, and
C                       ignored otherwise.  ML and MU may in fact be
C                       the band parameters for a matrix to which
C                       df/dy is only approximately equal.
C
C LIW    = the length of the array IWORK, as declared by the user.
C          (This will be checked by the solver.)
C
C Note: The base addresses of the work arrays must not be
C altered between calls to DLSODA for the same problem.
C The contents of the work arrays must not be altered
C between calls, except possibly for the conditional and
C optional inputs, and except for the last 3*NEQ words of RWORK.
C The latter space is used for internal scratch space, and so is
C available for use by the user outside DLSODA between calls, if
C desired (but not for use by F or JAC).
C
C JAC    = the name of the user-supplied routine to compute the
C          Jacobian matrix, df/dy, if JT = 1 or 4.  The JAC routine
C          is optional, but if the problem is expected to be stiff much
C          of the time, you are encouraged to supply JAC, for the sake
C          of efficiency.  (Alternatively, set JT = 2 or 5 to have
C          DLSODA compute df/dy internally by difference quotients.)
C          If and when DLSODA uses df/dy, it treats this NEQ by NEQ
C          matrix either as full (JT = 1 or 2), or as banded (JT =
C          4 or 5) with half-bandwidths ML and MU (discussed under
C          IWORK above).  In either case, if JT = 1 or 4, the JAC
C          routine must compute df/dy as a function of the scalar t
C          and the vector y.  It is to have the form
C               SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
C               DOUBLE PRECISION T, Y(*), PD(NROWPD,*)
C          where NEQ, T, Y, ML, MU, and NROWPD are input and the array
C          PD is to be loaded with partial derivatives (elements of
C          the Jacobian matrix) on output.  PD must be given a first
C          dimension of NROWPD.  T and Y have the same meaning as in
C          Subroutine F.
C               In the full matrix case (JT = 1), ML and MU are
C          ignored, and the Jacobian is to be loaded into PD in
C          columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
C               In the band matrix case (JT = 4), the elements
C          within the band are to be loaded into PD in columnwise
C          manner, with diagonal lines of df/dy loaded into the rows
C          of PD.  Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
C          ML and MU are the half-bandwidth parameters (see IWORK).
C          The locations in PD in the two triangular areas which
C          correspond to nonexistent matrix elements can be ignored
C          or loaded arbitrarily, as they are overwritten by DLSODA.
C               JAC need not provide df/dy exactly.  A crude
C          approximation (possibly with a smaller bandwidth) will do.
C               In either case, PD is preset to zero by the solver,
C          so that only the nonzero elements need be loaded by JAC.
C          Each call to JAC is preceded by a call to F with the same
C          arguments NEQ, T, and Y.  Thus to gain some efficiency,
C          intermediate quantities shared by both calculations may be
C          saved in a user Common block by F and not recomputed by JAC,
C          if desired.  Also, JAC may alter the Y array, if desired.
C          JAC must be declared External in the calling program.
C               Subroutine JAC may access user-defined quantities in
C          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
C          (dimensioned in JAC) and/or Y has length exceeding NEQ(1).
C          See the descriptions of NEQ and Y above.
C
C JT     = Jacobian type indicator.  Used only for input.
C          JT specifies how the Jacobian matrix df/dy will be
C          treated, if and when DLSODA requires this matrix.
C          JT has the following values and meanings:
C           1 means a user-supplied full (NEQ by NEQ) Jacobian.
C           2 means an internally generated (difference quotient) full
C             Jacobian (using NEQ extra calls to F per df/dy value).
C           4 means a user-supplied banded Jacobian.
C           5 means an internally generated banded Jacobian (using
C             ML+MU+1 extra calls to F per df/dy evaluation).
C          If JT = 1 or 4, the user must supply a Subroutine JAC
C          (the name is arbitrary) as described above under JAC.
C          If JT = 2 or 5, a dummy argument can be used.
C-----------------------------------------------------------------------
C Optional Inputs.
C
C The following is a list of the optional inputs provided for in the
C call sequence.  (See also Part 2.)  For each such input variable,
C this table lists its name as used in this documentation, its
C location in the call sequence, its meaning, and the default value.
C The use of any of these inputs requires IOPT = 1, and in that
C case all of these inputs are examined.  A value of zero for any
C of these optional inputs will cause the default value to be used.
C Thus to use a subset of the optional inputs, simply preload
C locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
C then set those of interest to nonzero values.
C
C Name    Location      Meaning and Default Value
C
C H0      RWORK(5)  the step size to be attempted on the first step.
C                   The default value is determined by the solver.
C
C HMAX    RWORK(6)  the maximum absolute step size allowed.
C                   The default value is infinite.
C
C HMIN    RWORK(7)  the minimum absolute step size allowed.
C                   The default value is 0.  (This lower bound is not
C                   enforced on the final step before reaching TCRIT
C                   when ITASK = 4 or 5.)
C
C IXPR    IWORK(5)  flag to generate extra printing at method switches.
C                   IXPR = 0 means no extra printing (the default).
C                   IXPR = 1 means print data on each switch.
C                   T, H, and NST will be printed on the same logical
C                   unit as used for error messages.
C
C MXSTEP  IWORK(6)  maximum number of (internally defined) steps
C                   allowed during one call to the solver.
C                   The default value is 500.
C
C MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
C                   warning that T + H = T on a step (H = step size).
C                   This must be positive to result in a non-default
C                   value.  The default value is 10.
C
C MXORDN  IWORK(8)  the maximum order to be allowed for the nonstiff
C                   (Adams) method.  the default value is 12.
C                   if MXORDN exceeds the default value, it will
C                   be reduced to the default value.
C                   MXORDN is held constant during the problem.
C
C MXORDS  IWORK(9)  the maximum order to be allowed for the stiff
C                   (BDF) method.  The default value is 5.
C                   If MXORDS exceeds the default value, it will
C                   be reduced to the default value.
C                   MXORDS is held constant during the problem.
C-----------------------------------------------------------------------
C Optional Outputs.
C
C As optional additional output from DLSODA, the variables listed
C below are quantities related to the performance of DLSODA
C which are available to the user.  These are communicated by way of
C the work arrays, but also have internal mnemonic names as shown.
C except where stated otherwise, all of these outputs are defined
C on any successful return from DLSODA, and on any return with
C ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return
C (ISTATE = -3), they will be unchanged from their existing values
C (if any), except possibly for TOLSF, LENRW, and LENIW.
C On any error return, outputs relevant to the error will be defined,
C as noted below.
C
C Name    Location      Meaning
C
C HU      RWORK(11) the step size in t last used (successfully).
C
C HCUR    RWORK(12) the step size to be attempted on the next step.
C
C TCUR    RWORK(13) the current value of the independent variable
C                   which the solver has actually reached, i.e. the
C                   current internal mesh point in t.  On output, TCUR
C                   will always be at least as far as the argument
C                   T, but may be farther (if interpolation was done).
C
C TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
C                   computed when a request for too much accuracy was
C                   detected (ISTATE = -3 if detected at the start of
C                   the problem, ISTATE = -2 otherwise).  If ITOL is
C                   left unaltered but RTOL and ATOL are uniformly
C                   scaled up by a factor of TOLSF for the next call,
C                   then the solver is deemed likely to succeed.
C                   (The user may also ignore TOLSF and alter the
C                   tolerance parameters in any other way appropriate.)
C
C TSW     RWORK(15) the value of t at the time of the last method
C                   switch, if any.
C
C NST     IWORK(11) the number of steps taken for the problem so far.
C
C NFE     IWORK(12) the number of f evaluations for the problem so far.
C
C NJE     IWORK(13) the number of Jacobian evaluations (and of matrix
C                   LU decompositions) for the problem so far.
C
C NQU     IWORK(14) the method order last used (successfully).
C
C NQCUR   IWORK(15) the order to be attempted on the next step.
C
C IMXER   IWORK(16) the index of the component of largest magnitude in
C                   the weighted local error vector ( E(i)/EWT(i) ),
C                   on an error return with ISTATE = -4 or -5.
C
C LENRW   IWORK(17) the length of RWORK actually required, assuming
C                   that the length of RWORK is to be fixed for the
C                   rest of the problem, and that switching may occur.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C LENIW   IWORK(18) the length of IWORK actually required, assuming
C                   that the length of IWORK is to be fixed for the
C                   rest of the problem, and that switching may occur.
C                   This is defined on normal returns and on an illegal
C                   input return for insufficient storage.
C
C MUSED   IWORK(19) the method indicator for the last successful step:
C                   1 means Adams (nonstiff), 2 means BDF (stiff).
C
C MCUR    IWORK(20) the current method indicator:
C                   1 means Adams (nonstiff), 2 means BDF (stiff).
C                   This is the method to be attempted
C                   on the next step.  Thus it differs from MUSED
C                   only if a method switch has just been made.
C
C The following two arrays are segments of the RWORK array which
C may also be of interest to the user as optional outputs.
C For each array, the table below gives its internal name,
C its base address in RWORK, and its description.
C
C Name    Base Address      Description
C
C YH      21             the Nordsieck history array, of size NYH by
C                        (NQCUR + 1), where NYH is the initial value
C                        of NEQ.  For j = 0,1,...,NQCUR, column j+1
C                        of YH contains HCUR**j/factorial(j) times
C                        the j-th derivative of the interpolating
C                        polynomial currently representing the solution,
C                        evaluated at T = TCUR.
C
C ACOR     LACOR         array of size NEQ used for the accumulated
C         (from Common   corrections on each step, scaled on output
C           as noted)    to represent the estimated local error in y
C                        on the last step.  This is the vector E in
C                        the description of the error control.  It is
C                        defined only on a successful return from
C                        DLSODA.  The base address LACOR is obtained by
C                        including in the user's program the
C                        following 2 lines:
C                           COMMON /DLS001/ RLS(218), ILS(37)
C                           LACOR = ILS(22)
C
C-----------------------------------------------------------------------
C Part 2.  Other Routines Callable.
C
C The following are optional calls which the user may make to
C gain additional capabilities in conjunction with DLSODA.
C (The routines XSETUN and XSETF are designed to conform to the
C SLATEC error handling package.)
C
C     Form of Call                  Function
C   CALL XSETUN(LUN)          set the logical unit number, LUN, for
C                             output of messages from DLSODA, if
C                             the default is not desired.
C                             The default value of LUN is 6.
C
C   CALL XSETF(MFLAG)         set a flag to control the printing of
C                             messages by DLSODA.
C                             MFLAG = 0 means do not print. (Danger:
C                             This risks losing valuable information.)
C                             MFLAG = 1 means print (the default).
C
C                             Either of the above calls may be made at
C                             any time and will take effect immediately.
C
C   CALL DSRCMA(RSAV,ISAV,JOB) saves and restores the contents of
C                             the internal Common blocks used by
C                             DLSODA (see Part 3 below).
C                             RSAV must be a real array of length 240
C                             or more, and ISAV must be an integer
C                             array of length 46 or more.
C                             JOB=1 means save Common into RSAV/ISAV.
C                             JOB=2 means restore Common from RSAV/ISAV.
C                                DSRCMA is useful if one is
C                             interrupting a run and restarting
C                             later, or alternating between two or
C                             more problems solved with DLSODA.
C
C   CALL DINTDY(,,,,,)        provide derivatives of y, of various
C        (see below)          orders, at a specified point t, if
C                             desired.  It may be called only after
C                             a successful return from DLSODA.
C
C The detailed instructions for using DINTDY are as follows.
C The form of the call is:
C
C   CALL DINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
C
C The input parameters are:
C
C T         = value of independent variable where answers are desired
C             (normally the same as the T last returned by DLSODA).
C             For valid results, T must lie between TCUR - HU and TCUR.
C             (See optional outputs for TCUR and HU.)
C K         = integer order of the derivative desired.  K must satisfy
C             0 .le. K .le. NQCUR, where NQCUR is the current order
C             (see optional outputs).  The capability corresponding
C             to K = 0, i.e. computing y(T), is already provided
C             by DLSODA directly.  Since NQCUR .ge. 1, the first
C             derivative dy/dt is always available with DINTDY.
C RWORK(21) = the base address of the history array YH.
C NYH       = column length of YH, equal to the initial value of NEQ.
C
C The output parameters are:
C
C DKY       = a real array of length NEQ containing the computed value
C             of the K-th derivative of y(t).
C IFLAG     = integer flag, returned as 0 if K and T were legal,
C             -1 if K was illegal, and -2 if T was illegal.
C             On an error return, a message is also written.
C-----------------------------------------------------------------------
C Part 3.  Common Blocks.
C
C If DLSODA is to be used in an overlay situation, the user
C must declare, in the primary overlay, the variables in:
C   (1) the call sequence to DLSODA, and
C   (2) the two internal Common blocks
C         /DLS001/  of length  255  (218 double precision words
C                      followed by 37 integer words),
C         /DLSA01/  of length  31    (22 double precision words
C                      followed by  9 integer words).
C
C If DLSODA is used on a system in which the contents of internal
C Common blocks are not preserved between calls, the user should
C declare the above Common blocks in the calling program to insure
C that their contents are preserved.
C
C If the solution of a given problem by DLSODA is to be interrupted
C and then later continued, such as when restarting an interrupted run
C or alternating between two or more problems, the user should save,
C following the return from the last DLSODA call prior to the
C interruption, the contents of the call sequence variables and the
C internal Common blocks, and later restore these values before the
C next DLSODA call for that problem.  To save and restore the Common
C blocks, use Subroutine DSRCMA (see Part 2 above).
C
C-----------------------------------------------------------------------
C Part 4.  Optionally Replaceable Solver Routines.
C
C Below is a description of a routine in the DLSODA package which
C relates to the measurement of errors, and can be
C replaced by a user-supplied version, if desired.  However, since such
C a replacement may have a major impact on performance, it should be
C done only when absolutely necessary, and only with great caution.
C (Note: The means by which the package version of a routine is
C superseded by the user's version may be system-dependent.)
C
C (a) DEWSET.
C The following subroutine is called just before each internal
C integration step, and sets the array of error weights, EWT, as
C described under ITOL/RTOL/ATOL above:
C     Subroutine DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
C where NEQ, ITOL, RTOL, and ATOL are as in the DLSODA call sequence,
C YCUR contains the current dependent variable vector, and
C EWT is the array of weights set by DEWSET.
C
C If the user supplies this subroutine, it must return in EWT(i)
C (i = 1,...,NEQ) a positive quantity suitable for comparing errors
C in y(i) to.  The EWT array returned by DEWSET is passed to the
C DMNORM routine, and also used by DLSODA in the computation
C of the optional output IMXER, and the increments for difference
C quotient Jacobians.
C
C In the user-supplied version of DEWSET, it may be desirable to use
C the current values of derivatives of y.  Derivatives up to order NQ
C are available from the history array YH, described above under
C optional outputs.  In DEWSET, YH is identical to the YCUR array,
C extended to NQ + 1 columns with a column length of NYH and scale
C factors of H**j/factorial(j).  On the first call for the problem,
C given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
C NYH is the initial value of NEQ.  The quantities NQ, H, and NST
C can be obtained by including in DEWSET the statements:
C     DOUBLE PRECISION RLS
C     COMMON /DLS001/ RLS(218),ILS(37)
C     NQ = ILS(33)
C     NST = ILS(34)
C     H = RLS(212)
C Thus, for example, the current value of dy/dt can be obtained as
C YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
C unnecessary when NST = 0).
C-----------------------------------------------------------------------
C
C***REVISION HISTORY  (YYYYMMDD)
C 19811102  DATE WRITTEN
C 19820126  Fixed bug in tests of work space lengths;
C           minor corrections in main prologue and comments.
C 19870330  Major update: corrected comments throughout;
C           removed TRET from Common; rewrote EWSET with 4 loops;
C           fixed t test in INTDY; added Cray directives in STODA;
C           in STODA, fixed DELP init. and logic around PJAC call;
C           combined routines to save/restore Common;
C           passed LEVEL = 0 in error message calls (except run abort).
C 19970225  Fixed lines setting JSTART = -2 in Subroutine LSODA.
C 20010425  Major update: convert source lines to upper case;
C           added *DECK lines; changed from 1 to * in dummy dimensions;
C           changed names R1MACH/D1MACH to RUMACH/DUMACH;
C           renamed routines for uniqueness across single/double prec.;
C           converted intrinsic names to generic form;
C           removed ILLIN and NTREP (data loaded) from Common;
C           removed all 'own' variables from Common;
C           changed error messages to quoted strings;
C           replaced XERRWV/XERRWD with 1993 revised version;
C           converted prologues, comments, error messages to mixed case;
C           numerous corrections to prologues and internal comments.
C 20010507  Converted single precision source to double precision.
C 20010613  Revised excess accuracy test (to match rest of ODEPACK).
C 20010808  Fixed bug in DPRJA (matrix in DBNORM call).
C 20020502  Corrected declarations in descriptions of user routines.
C 20031105  Restored 'own' variables to Common blocks, to enable
C           interrupt/restart feature.
C 20031112  Added SAVE statements for data-loaded constants.
C
C-----------------------------------------------------------------------
C Other routines in the DLSODA package.
C
C In addition to Subroutine DLSODA, the DLSODA package includes the
C following subroutines and function routines:
C  DINTDY   computes an interpolated value of the y vector at t = TOUT.
C  DSTODA   is the core integrator, which does one step of the
C           integration and the associated error control.
C  DCFODE   sets all method coefficients and test constants.
C  DPRJA    computes and preprocesses the Jacobian matrix J = df/dy
C           and the Newton iteration matrix P = I - h*l0*J.
C  DSOLSY   manages solution of linear system in chord iteration.
C  DEWSET   sets the error weight vector EWT before each step.
C  DMNORM   computes the weighted max-norm of a vector.
C  DFNORM   computes the norm of a full matrix consistent with the
C           weighted max-norm on vectors.
C  DBNORM   computes the norm of a band matrix consistent with the
C           weighted max-norm on vectors.
C  DSRCMA   is a user-callable routine to save and restore
C           the contents of the internal Common blocks.
C  DGEFA and DGESL   are routines from LINPACK for solving full
C           systems of linear algebraic equations.
C  DGBFA and DGBSL   are routines from LINPACK for solving banded
C           linear systems.
C  DUMACH   computes the unit roundoff in a machine-independent manner.
C  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH  handle the printing of all
C           error messages and warnings.  XERRWD is machine-dependent.
C Note:  DMNORM, DFNORM, DBNORM, DUMACH, IXSAV, and IUMACH are
C function routines.  All the others are subroutines.
C
C-----------------------------------------------------------------------
      EXTERNAL DPRJA, DSOLSY
      DOUBLE PRECISION DUMACH, DMNORM
      INTEGER INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER INSUFR, INSUFI, IXPR, IOWNS2, JTYP, MUSED, MXORDN, MXORDS
      INTEGER I, I1, I2, IFLAG, IMXER, KGO, LF0,
     1   LENIW, LENRW, LENWM, ML, MORD, MU, MXHNL0, MXSTP0
      INTEGER LEN1, LEN1C, LEN1N, LEN1S, LEN2, LENIWC, LENRWC
      DOUBLE PRECISION ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION TSW, ROWNS2, PDNORM
      DOUBLE PRECISION ATOLI, AYI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI,
     1   TCRIT, TDIST, TNEXT, TOL, TOLSF, TP, SIZE, SUM, W0
      DIMENSION MORD(2)
      LOGICAL IHIT
      CHARACTER*60 MSG
      SAVE MORD, MXSTP0, MXHNL0
C-----------------------------------------------------------------------
C The following two internal Common blocks contain
C (a) variables which are local to any subroutine but whose values must
C     be preserved between calls to the routine ("own" variables), and
C (b) variables which are communicated between subroutines.
C The block DLS001 is declared in subroutines DLSODA, DINTDY, DSTODA,
C DPRJA, and DSOLSY.
C The block DLSA01 is declared in subroutines DLSODA, DSTODA, and DPRJA.
C Groups of variables are replaced by dummy arrays in the Common
C declarations in routines where those variables are not used.
C-----------------------------------------------------------------------
      COMMON /DLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   INIT, MXSTEP, MXHNIL, NHNIL, NSLAST, NYH, IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
C
      COMMON /DLSA01/ TSW, ROWNS2(20), PDNORM,
     1   INSUFR, INSUFI, IXPR, IOWNS2(2), JTYP, MUSED, MXORDN, MXORDS
C
      DATA MORD(1),MORD(2)/12,5/, MXSTP0/50000/, MXHNL0/10/
C-----------------------------------------------------------------------
C Block A.
C This code block is executed on every call.
C It tests ISTATE and ITASK for legality and branches appropriately.
C If ISTATE .gt. 1 but the flag INIT shows that initialization has
C not yet been done, an error return occurs.
C If ISTATE = 1 and TOUT = T, return immediately.
C-----------------------------------------------------------------------
      IF (ISTATE .LT. 1 .OR. ISTATE .GT. 3) GO TO 601
      IF (ITASK .LT. 1 .OR. ITASK .GT. 5) GO TO 602
      IF (ISTATE .EQ. 1) GO TO 10
      IF (INIT .EQ. 0) GO TO 603
      IF (ISTATE .EQ. 2) GO TO 200
      GO TO 20
 10   INIT = 0
      IF (TOUT .EQ. T) RETURN
C-----------------------------------------------------------------------
C Block B.
C The next code block is executed for the initial call (ISTATE = 1),
C or for a continuation call with parameter changes (ISTATE = 3).
C It contains checking of all inputs and various initializations.
C
C First check legality of the non-optional inputs NEQ, ITOL, IOPT,
C JT, ML, and MU.
C-----------------------------------------------------------------------
 20   IF (NEQ(1) .LE. 0) GO TO 604
      IF (ISTATE .EQ. 1) GO TO 25
      IF (NEQ(1) .GT. N) GO TO 605
 25   N = NEQ(1)
      IF (ITOL .LT. 1 .OR. ITOL .GT. 4) GO TO 606
      IF (IOPT .LT. 0 .OR. IOPT .GT. 1) GO TO 607
      IF (JT .EQ. 3 .OR. JT .LT. 1 .OR. JT .GT. 5) GO TO 608
      JTYP = JT
      IF (JT .LE. 2) GO TO 30
      ML = IWORK(1)
      MU = IWORK(2)
      IF (ML .LT. 0 .OR. ML .GE. N) GO TO 609
      IF (MU .LT. 0 .OR. MU .GE. N) GO TO 610
 30   CONTINUE
C Next process and check the optional inputs. --------------------------
      IF (IOPT .EQ. 1) GO TO 40
      IXPR = 0
      MXSTEP = MXSTP0
      MXHNIL = MXHNL0
      HMXI = 0.0D0
      HMIN = 0.0D0
      IF (ISTATE .NE. 1) GO TO 60
      H0 = 0.0D0
      MXORDN = MORD(1)
      MXORDS = MORD(2)
      GO TO 60
 40   IXPR = IWORK(5)
      IF (IXPR .LT. 0 .OR. IXPR .GT. 1) GO TO 611
      MXSTEP = IWORK(6)
      IF (MXSTEP .LT. 0) GO TO 612
      IF (MXSTEP .EQ. 0) MXSTEP = MXSTP0
      MXHNIL = IWORK(7)
      IF (MXHNIL .LT. 0) GO TO 613
      IF (MXHNIL .EQ. 0) MXHNIL = MXHNL0
      IF (ISTATE .NE. 1) GO TO 50
      H0 = RWORK(5)
      MXORDN = IWORK(8)
      IF (MXORDN .LT. 0) GO TO 628
      IF (MXORDN .EQ. 0) MXORDN = 100
      MXORDN = MIN(MXORDN,MORD(1))
      MXORDS = IWORK(9)
      IF (MXORDS .LT. 0) GO TO 629
      IF (MXORDS .EQ. 0) MXORDS = 100
      MXORDS = MIN(MXORDS,MORD(2))
      IF ((TOUT - T)*H0 .LT. 0.0D0) GO TO 614
 50   HMAX = RWORK(6)
      IF (HMAX .LT. 0.0D0) GO TO 615
      HMXI = 0.0D0
      IF (HMAX .GT. 0.0D0) HMXI = 1.0D0/HMAX
      HMIN = RWORK(7)
      IF (HMIN .LT. 0.0D0) GO TO 616
C-----------------------------------------------------------------------
C Set work array pointers and check lengths LRW and LIW.
C If ISTATE = 1, METH is initialized to 1 here to facilitate the
C checking of work space lengths.
C Pointers to segments of RWORK and IWORK are named by prefixing L to
C the name of the segment.  E.g., the segment YH starts at RWORK(LYH).
C Segments of RWORK (in order) are denoted  YH, WM, EWT, SAVF, ACOR.
C If the lengths provided are insufficient for the current method,
C an error return occurs.  This is treated as illegal input on the
C first call, but as a problem interruption with ISTATE = -7 on a
C continuation call.  If the lengths are sufficient for the current
C method but not for both methods, a warning message is sent.
C-----------------------------------------------------------------------
 60   IF (ISTATE .EQ. 1) METH = 1
      IF (ISTATE .EQ. 1) NYH = N
      LYH = 21
      LEN1N = 20 + (MXORDN + 1)*NYH
      LEN1S = 20 + (MXORDS + 1)*NYH
      LWM = LEN1S + 1
      IF (JT .LE. 2) LENWM = N*N + 2
      IF (JT .GE. 4) LENWM = (2*ML + MU + 1)*N + 2
      LEN1S = LEN1S + LENWM
      LEN1C = LEN1N
      IF (METH .EQ. 2) LEN1C = LEN1S
      LEN1 = MAX(LEN1N,LEN1S)
      LEN2 = 3*N
      LENRW = LEN1 + LEN2
      LENRWC = LEN1C + LEN2
      IWORK(17) = LENRW
      LIWM = 1
      LENIW = 20 + N
      LENIWC = 20
      IF (METH .EQ. 2) LENIWC = LENIW
      IWORK(18) = LENIW
      IF (ISTATE .EQ. 1 .AND. LRW .LT. LENRWC) GO TO 617
      IF (ISTATE .EQ. 1 .AND. LIW .LT. LENIWC) GO TO 618
      IF (ISTATE .EQ. 3 .AND. LRW .LT. LENRWC) GO TO 550
      IF (ISTATE .EQ. 3 .AND. LIW .LT. LENIWC) GO TO 555
      LEWT = LEN1 + 1
      INSUFR = 0
      IF (LRW .GE. LENRW) GO TO 65
      INSUFR = 2
      LEWT = LEN1C + 1
      MSG='DLSODA-  Warning.. RWORK length is sufficient for now, but  '
      CALL XERRWD (MSG, 60, 103, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      may not be later.  Integration will proceed anyway.   '
      CALL XERRWD (MSG, 60, 103, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      Length needed is LENRW = I1, while LRW = I2.'
      CALL XERRWD (MSG, 50, 103, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
 65   LSAVF = LEWT + N
      LACOR = LSAVF + N
      INSUFI = 0
      IF (LIW .GE. LENIW) GO TO 70
      INSUFI = 2
      MSG='DLSODA-  Warning.. IWORK length is sufficient for now, but  '
      CALL XERRWD (MSG, 60, 104, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      may not be later.  Integration will proceed anyway.   '
      CALL XERRWD (MSG, 60, 104, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      Length needed is LENIW = I1, while LIW = I2.'
      CALL XERRWD (MSG, 50, 104, 0, 2, LENIW, LIW, 0, 0.0D0, 0.0D0)
 70   CONTINUE
C Check RTOL and ATOL for legality. ------------------------------------
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 75 I = 1,N
        IF (ITOL .GE. 3) RTOLI = RTOL(I)
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        IF (RTOLI .LT. 0.0D0) GO TO 619
        IF (ATOLI .LT. 0.0D0) GO TO 620
 75     CONTINUE
      IF (ISTATE .EQ. 1) GO TO 100
C If ISTATE = 3, set flag to signal parameter changes to DSTODA. -------
      JSTART = -1
      IF (N .EQ. NYH) GO TO 200
C NEQ was reduced.  Zero part of YH to avoid undefined references. -----
      I1 = LYH + L*NYH
      I2 = LYH + (MAXORD + 1)*NYH - 1
      IF (I1 .GT. I2) GO TO 200
      DO 95 I = I1,I2
 95     RWORK(I) = 0.0D0
      GO TO 200
C-----------------------------------------------------------------------
C Block C.
C The next block is for the initial call only (ISTATE = 1).
C It contains all remaining initializations, the initial call to F,
C and the calculation of the initial step size.
C The error weights in EWT are inverted after being loaded.
C-----------------------------------------------------------------------
 100  UROUND = DUMACH()
      TN = T
      TSW = T
      MAXORD = MXORDN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 110
      TCRIT = RWORK(1)
      IF ((TCRIT - TOUT)*(TOUT - T) .LT. 0.0D0) GO TO 625
      IF (H0 .NE. 0.0D0 .AND. (T + H0 - TCRIT)*H0 .GT. 0.0D0)
     1   H0 = TCRIT - T
 110  JSTART = 0
      NHNIL = 0
      NST = 0
      NJE = 0
      NSLAST = 0
      HU = 0.0D0
      NQU = 0
      MUSED = 0
      MITER = 0
      CCMAX = 0.3D0
      MAXCOR = 3
      MSBP = 20
      MXNCF = 10
C Initial call to F.  (LF0 points to YH(*,2).) -------------------------
      LF0 = LYH + NYH
      CALL F (KP,Lfs, Xek, NEQ, T, Y, RWORK(LF0))
      NFE = 1
C Load the initial value vector in YH. ---------------------------------
      DO 115 I = 1,N
 115    RWORK(I+LYH-1) = Y(I)
C Load and invert the EWT array.  (H is temporarily set to 1.0.) -------
      NQ = 1
      H = 1.0D0
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 120 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 621
 120    RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
C-----------------------------------------------------------------------
C The coding below computes the step size, H0, to be attempted on the
C first step, unless the user has supplied a value for this.
C First check that TOUT - T differs significantly from zero.
C A scalar tolerance quantity TOL is computed, as MAX(RTOL(i))
C if this is positive, or MAX(ATOL(i)/ABS(Y(i))) otherwise, adjusted
C so as to be between 100*UROUND and 1.0E-3.
C Then the computed value H0 is given by:
C
C   H0**(-2)  =  1./(TOL * w0**2)  +  TOL * (norm(F))**2
C
C where   w0     = MAX ( ABS(T), ABS(TOUT) ),
C         F      = the initial value of the vector f(t,y), and
C         norm() = the weighted vector norm used throughout, given by
C                  the DMNORM function routine, and weighted by the
C                  tolerances initially loaded into the EWT array.
C The sign of H0 is inferred from the initial values of TOUT and T.
C ABS(H0) is made .le. ABS(TOUT-T) in any case.
C-----------------------------------------------------------------------
      IF (H0 .NE. 0.0D0) GO TO 180
      TDIST = ABS(TOUT - T)
      W0 = MAX(ABS(T),ABS(TOUT))
      IF (TDIST .LT. 2.0D0*UROUND*W0) GO TO 622
      TOL = RTOL(1)
      IF (ITOL .LE. 2) GO TO 140
      DO 130 I = 1,N
 130    TOL = MAX(TOL,RTOL(I))
 140  IF (TOL .GT. 0.0D0) GO TO 160
      ATOLI = ATOL(1)
      DO 150 I = 1,N
        IF (ITOL .EQ. 2 .OR. ITOL .EQ. 4) ATOLI = ATOL(I)
        AYI = ABS(Y(I))
        IF (AYI .NE. 0.0D0) TOL = MAX(TOL,ATOLI/AYI)
 150    CONTINUE
 160  TOL = MAX(TOL,100.0D0*UROUND)
      TOL = MIN(TOL,0.001D0)
      SUM = DMNORM (N, RWORK(LF0), RWORK(LEWT))
      SUM = 1.0D0/(TOL*W0*W0) + TOL*SUM**2
      H0 = 1.0D0/SQRT(SUM)
      H0 = MIN(H0,TDIST)
      H0 = SIGN(H0,TOUT-T)
C Adjust H0 if necessary to meet HMAX bound. ---------------------------
 180  RH = ABS(H0)*HMXI
      IF (RH .GT. 1.0D0) H0 = H0/RH
C Load H with H0 and scale YH(*,2) by H0. ------------------------------
      H = H0
      DO 190 I = 1,N
 190    RWORK(I+LF0-1) = H0*RWORK(I+LF0-1)
      GO TO 270
C-----------------------------------------------------------------------
C Block D.
C The next code block is for continuation calls only (ISTATE = 2 or 3)
C and is to check stop conditions before taking a step.
C-----------------------------------------------------------------------
 200  NSLAST = NST
      GO TO (210, 250, 220, 230, 240), ITASK
 210  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 220  TP = TN - HU*(1.0D0 + 100.0D0*UROUND)
      IF ((TP - TOUT)*H .GT. 0.0D0) GO TO 623
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      T = TN
      GO TO 400
 230  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624
      IF ((TCRIT - TOUT)*H .LT. 0.0D0) GO TO 625
      IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 245
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      IF (IFLAG .NE. 0) GO TO 627
      T = TOUT
      GO TO 420
 240  TCRIT = RWORK(1)
      IF ((TN - TCRIT)*H .GT. 0.0D0) GO TO 624
 245  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
      IF (IHIT) T = TCRIT
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)
      IF (ISTATE .EQ. 2 .AND. JSTART .GE. 0) JSTART = -2
C-----------------------------------------------------------------------
C Block E.
C The next block is normally executed for all calls and contains
C the call to the one-step core integrator DSTODA.
C
C This is a looping point for the integration steps.
C
C First check for too many steps being taken, update EWT (if not at
C start of problem), check for too much accuracy being requested, and
C check for H below the roundoff level in T.
C-----------------------------------------------------------------------
 250  CONTINUE
      IF (METH .EQ. MUSED) GO TO 255
      IF (INSUFR .EQ. 1) GO TO 550
      IF (INSUFI .EQ. 1) GO TO 555
 255  IF ((NST-NSLAST) .GE. MXSTEP) GO TO 500
      CALL DEWSET (N, ITOL, RTOL, ATOL, RWORK(LYH), RWORK(LEWT))
      DO 260 I = 1,N
        IF (RWORK(I+LEWT-1) .LE. 0.0D0) GO TO 510
 260    RWORK(I+LEWT-1) = 1.0D0/RWORK(I+LEWT-1)
 270  TOLSF = UROUND*DMNORM (N, RWORK(LYH), RWORK(LEWT))
      IF (TOLSF .LE. 1.0D0) GO TO 280
      TOLSF = TOLSF*2.0D0
      IF (NST .EQ. 0) GO TO 626
      GO TO 520
 280  IF ((TN + H) .NE. TN) GO TO 290
      NHNIL = NHNIL + 1
      IF (NHNIL .GT. MXHNIL) GO TO 290
      MSG = 'DLSODA-  Warning..Internal T (=R1) and H (=R2) are'
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      such that in the machine, T + H = T on the next step  '
      CALL XERRWD (MSG, 60, 101, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '     (H = step size). Solver will continue anyway.'
      CALL XERRWD (MSG, 50, 101, 0, 0, 0, 0, 2, TN, H)
      IF (NHNIL .LT. MXHNIL) GO TO 290
      MSG = 'DLSODA-  Above warning has been issued I1 times.  '
      CALL XERRWD (MSG, 50, 102, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '     It will not be issued again for this problem.'
      CALL XERRWD (MSG, 50, 102, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
 290  CONTINUE
C-----------------------------------------------------------------------
C   CALL DSTODA(NEQ,Y,YH,NYH,YH,EWT,SAVF,ACOR,WM,IWM,F,JAC,DPRJA,DSOLSY)
C-----------------------------------------------------------------------
      CALL DSTODA (KP,Lfs, Xek,NEQ, Y, RWORK(LYH), NYH, RWORK(LYH),
     1   RWORK(LEWT),RWORK(LSAVF), RWORK(LACOR), RWORK(LWM),
     2   IWORK(LIWM), F, JAC, DPRJA, DSOLSY)
      KGO = 1 - KFLAG
      GO TO (300, 530, 540), KGO
C-----------------------------------------------------------------------
C Block F.
C The following block handles the case of a successful return from the
C core integrator (KFLAG = 0).
C If a method switch was just made, record TSW, reset MAXORD,
C set JSTART to -1 to signal DSTODA to complete the switch,
C and do extra printing of data if IXPR = 1.
C Then, in any case, check for stop conditions.
C-----------------------------------------------------------------------
 300  INIT = 1
      IF (METH .EQ. MUSED) GO TO 310
      TSW = TN
      MAXORD = MXORDN
      IF (METH .EQ. 2) MAXORD = MXORDS
      IF (METH .EQ. 2) RWORK(LWM) = SQRT(UROUND)
      INSUFR = MIN(INSUFR,1)
      INSUFI = MIN(INSUFI,1)
      JSTART = -1
      IF (IXPR .EQ. 0) GO TO 310
      IF (METH .EQ. 2) THEN
      MSG='DLSODA- A switch to the BDF (stiff) method has occurred     '
      CALL XERRWD (MSG, 60, 105, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      ENDIF
      IF (METH .EQ. 1) THEN
      MSG='DLSODA- A switch to the Adams (nonstiff) method has occurred'
      CALL XERRWD (MSG, 60, 106, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      ENDIF
      MSG='     at T = R1,  tentative step size H = R2,  step NST = I1 '
      CALL XERRWD (MSG, 60, 107, 0, 1, NST, 0, 2, TN, H)
 310  GO TO (320, 400, 330, 340, 350), ITASK
C ITASK = 1.  If TOUT has been reached, interpolate. -------------------
 320  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 250
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
C ITASK = 3.  Jump to exit if TOUT was reached. ------------------------
 330  IF ((TN - TOUT)*H .GE. 0.0D0) GO TO 400
      GO TO 250
C ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary.
 340  IF ((TN - TOUT)*H .LT. 0.0D0) GO TO 345
      CALL DINTDY (TOUT, 0, RWORK(LYH), NYH, Y, IFLAG)
      T = TOUT
      GO TO 420
 345  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
      IF (IHIT) GO TO 400
      TNEXT = TN + H*(1.0D0 + 4.0D0*UROUND)
      IF ((TNEXT - TCRIT)*H .LE. 0.0D0) GO TO 250
      H = (TCRIT - TN)*(1.0D0 - 4.0D0*UROUND)
      IF (JSTART .GE. 0) JSTART = -2
      GO TO 250
C ITASK = 5.  See if TCRIT was reached and jump to exit. ---------------
 350  HMX = ABS(TN) + ABS(H)
      IHIT = ABS(TN - TCRIT) .LE. 100.0D0*UROUND*HMX
C-----------------------------------------------------------------------
C Block G.
C The following block handles all successful returns from DLSODA.
C If ITASK .ne. 1, Y is loaded from YH and T is set accordingly.
C ISTATE is set to 2, and the optional outputs are loaded into the
C work arrays before returning.
C-----------------------------------------------------------------------
 400  DO 410 I = 1,N
 410    Y(I) = RWORK(I+LYH-1)
      T = TN
      IF (ITASK .NE. 4 .AND. ITASK .NE. 5) GO TO 420
      IF (IHIT) T = TCRIT
 420  ISTATE = 2
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      RWORK(15) = TSW
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = MUSED
      IWORK(20) = METH
      RETURN
C-----------------------------------------------------------------------
C Block H.
C The following block handles all unsuccessful returns other than
C those for illegal input.  First the error message routine is called.
C If there was an error test or convergence test failure, IMXER is set.
C Then Y is loaded from YH and T is set to TN.
C The optional outputs are loaded into the work arrays before returning.
C-----------------------------------------------------------------------
C The maximum number of steps was taken before reaching TOUT. ----------
 500  MSG = 'DLSODA-  At current T (=R1), MXSTEP (=I1) steps   '
      CALL XERRWD (MSG, 50, 201, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      taken on this call before reaching TOUT     '
      CALL XERRWD (MSG, 50, 201, 0, 1, MXSTEP, 0, 1, TN, 0.0D0)
      ISTATE = -1
      GO TO 580
C EWT(i) .le. 0.0 for some i (not at start of problem). ----------------
 510  EWTI = RWORK(LEWT+I-1)
      MSG = 'DLSODA-  At T (=R1), EWT(I1) has become R2 .le. 0.'
      CALL XERRWD (MSG, 50, 202, 0, 1, I, 0, 2, TN, EWTI)
      ISTATE = -6
      GO TO 580
C Too much accuracy requested for machine precision. -------------------
 520  MSG = 'DLSODA-  At T (=R1), too much accuracy requested  '
      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      for precision of machine..  See TOLSF (=R2) '
      CALL XERRWD (MSG, 50, 203, 0, 0, 0, 0, 2, TN, TOLSF)
      RWORK(14) = TOLSF
      ISTATE = -2
      GO TO 580
C KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. -----
 530  MSG = 'DLSODA-  At T(=R1) and step size H(=R2), the error'
      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      test failed repeatedly or with ABS(H) = HMIN'
      CALL XERRWD (MSG, 50, 204, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -4
      GO TO 560
C KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ----
 540  MSG = 'DLSODA-  At T (=R1) and step size H (=R2), the    '
      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      corrector convergence failed repeatedly     '
      CALL XERRWD (MSG, 50, 205, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG = '      or with ABS(H) = HMIN   '
      CALL XERRWD (MSG, 30, 205, 0, 0, 0, 0, 2, TN, H)
      ISTATE = -5
      GO TO 560
C RWORK length too small to proceed. -----------------------------------
 550  MSG = 'DLSODA-  At current T(=R1), RWORK length too small'
      CALL XERRWD (MSG, 50, 206, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      to proceed.  The integration was otherwise successful.'
      CALL XERRWD (MSG, 60, 206, 0, 0, 0, 0, 1, TN, 0.0D0)
      ISTATE = -7
      GO TO 580
C IWORK length too small to proceed. -----------------------------------
 555  MSG = 'DLSODA-  At current T(=R1), IWORK length too small'
      CALL XERRWD (MSG, 50, 207, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      to proceed.  The integration was otherwise successful.'
      CALL XERRWD (MSG, 60, 207, 0, 0, 0, 0, 1, TN, 0.0D0)
      ISTATE = -7
      GO TO 580
C Compute IMXER if relevant. -------------------------------------------
 560  BIG = 0.0D0
      IMXER = 1
      DO 570 I = 1,N
        SIZE = ABS(RWORK(I+LACOR-1)*RWORK(I+LEWT-1))
        IF (BIG .GE. SIZE) GO TO 570
        BIG = SIZE
        IMXER = I
 570    CONTINUE
      IWORK(16) = IMXER
C Set Y vector, T, and optional outputs. -------------------------------
 580  DO 590 I = 1,N
 590    Y(I) = RWORK(I+LYH-1)
      T = TN
      RWORK(11) = HU
      RWORK(12) = H
      RWORK(13) = TN
      RWORK(15) = TSW
      IWORK(11) = NST
      IWORK(12) = NFE
      IWORK(13) = NJE
      IWORK(14) = NQU
      IWORK(15) = NQ
      IWORK(19) = MUSED
      IWORK(20) = METH
      RETURN
C-----------------------------------------------------------------------
C Block I.
C The following block handles all error returns due to illegal input
C (ISTATE = -3), as detected before calling the core integrator.
C First the error message routine is called.  If the illegal input
C is a negative ISTATE, the run is aborted (apparent infinite loop).
C-----------------------------------------------------------------------
 601  MSG = 'DLSODA-  ISTATE (=I1) illegal.'
      CALL XERRWD (MSG, 30, 1, 0, 1, ISTATE, 0, 0, 0.0D0, 0.0D0)
      IF (ISTATE .LT. 0) GO TO 800
      GO TO 700
 602  MSG = 'DLSODA-  ITASK (=I1) illegal. '
      CALL XERRWD (MSG, 30, 2, 0, 1, ITASK, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 603  MSG = 'DLSODA-  ISTATE .gt. 1 but DLSODA not initialized.'
      CALL XERRWD (MSG, 50, 3, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 604  MSG = 'DLSODA-  NEQ (=I1) .lt. 1     '
      CALL XERRWD (MSG, 30, 4, 0, 1, NEQ(1), 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 605  MSG = 'DLSODA-  ISTATE = 3 and NEQ increased (I1 to I2). '
      CALL XERRWD (MSG, 50, 5, 0, 2, N, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 606  MSG = 'DLSODA-  ITOL (=I1) illegal.  '
      CALL XERRWD (MSG, 30, 6, 0, 1, ITOL, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 607  MSG = 'DLSODA-  IOPT (=I1) illegal.  '
      CALL XERRWD (MSG, 30, 7, 0, 1, IOPT, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 608  MSG = 'DLSODA-  JT (=I1) illegal.    '
      CALL XERRWD (MSG, 30, 8, 0, 1, JT, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 609  MSG = 'DLSODA-  ML (=I1) illegal: .lt.0 or .ge.NEQ (=I2) '
      CALL XERRWD (MSG, 50, 9, 0, 2, ML, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 610  MSG = 'DLSODA-  MU (=I1) illegal: .lt.0 or .ge.NEQ (=I2) '
      CALL XERRWD (MSG, 50, 10, 0, 2, MU, NEQ(1), 0, 0.0D0, 0.0D0)
      GO TO 700
 611  MSG = 'DLSODA-  IXPR (=I1) illegal.  '
      CALL XERRWD (MSG, 30, 11, 0, 1, IXPR, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 612  MSG = 'DLSODA-  MXSTEP (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 12, 0, 1, MXSTEP, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 613  MSG = 'DLSODA-  MXHNIL (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 13, 0, 1, MXHNIL, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 614  MSG = 'DLSODA-  TOUT (=R1) behind T (=R2)      '
      CALL XERRWD (MSG, 40, 14, 0, 0, 0, 0, 2, TOUT, T)
      MSG = '      Integration direction is given by H0 (=R1)  '
      CALL XERRWD (MSG, 50, 14, 0, 0, 0, 0, 1, H0, 0.0D0)
      GO TO 700
 615  MSG = 'DLSODA-  HMAX (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 15, 0, 0, 0, 0, 1, HMAX, 0.0D0)
      GO TO 700
 616  MSG = 'DLSODA-  HMIN (=R1) .lt. 0.0  '
      CALL XERRWD (MSG, 30, 16, 0, 0, 0, 0, 1, HMIN, 0.0D0)
      GO TO 700
 617  MSG='DLSODA-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)'
      CALL XERRWD (MSG, 60, 17, 0, 2, LENRW, LRW, 0, 0.0D0, 0.0D0)
      GO TO 700
 618  MSG='DLSODA-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)'
      CALL XERRWD (MSG, 60, 18, 0, 2, LENIW, LIW, 0, 0.0D0, 0.0D0)
      GO TO 700
 619  MSG = 'DLSODA-  RTOL(I1) is R1 .lt. 0.0        '
      CALL XERRWD (MSG, 40, 19, 0, 1, I, 0, 1, RTOLI, 0.0D0)
      GO TO 700
 620  MSG = 'DLSODA-  ATOL(I1) is R1 .lt. 0.0        '
      CALL XERRWD (MSG, 40, 20, 0, 1, I, 0, 1, ATOLI, 0.0D0)
      GO TO 700
 621  EWTI = RWORK(LEWT+I-1)
      MSG = 'DLSODA-  EWT(I1) is R1 .le. 0.0         '
      CALL XERRWD (MSG, 40, 21, 0, 1, I, 0, 1, EWTI, 0.0D0)
      GO TO 700
 622  MSG='DLSODA-  TOUT(=R1) too close to T(=R2) to start integration.'
      CALL XERRWD (MSG, 60, 22, 0, 0, 0, 0, 2, TOUT, T)
      GO TO 700
 623  MSG='DLSODA-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  '
      CALL XERRWD (MSG, 60, 23, 0, 1, ITASK, 0, 2, TOUT, TP)
      GO TO 700
 624  MSG='DLSODA-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   '
      CALL XERRWD (MSG, 60, 24, 0, 0, 0, 0, 2, TCRIT, TN)
      GO TO 700
 625  MSG='DLSODA-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   '
      CALL XERRWD (MSG, 60, 25, 0, 0, 0, 0, 2, TCRIT, TOUT)
      GO TO 700
 626  MSG = 'DLSODA-  At start of problem, too much accuracy   '
      CALL XERRWD (MSG, 50, 26, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
      MSG='      requested for precision of machine..  See TOLSF (=R1) '
      CALL XERRWD (MSG, 60, 26, 0, 0, 0, 0, 1, TOLSF, 0.0D0)
      RWORK(14) = TOLSF
      GO TO 700
 627  MSG = 'DLSODA-  Trouble in DINTDY.  ITASK = I1, TOUT = R1'
      CALL XERRWD (MSG, 50, 27, 0, 1, ITASK, 0, 1, TOUT, 0.0D0)
      GO TO 700
 628  MSG = 'DLSODA-  MXORDN (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 28, 0, 1, MXORDN, 0, 0, 0.0D0, 0.0D0)
      GO TO 700
 629  MSG = 'DLSODA-  MXORDS (=I1) .lt. 0  '
      CALL XERRWD (MSG, 30, 29, 0, 1, MXORDS, 0, 0, 0.0D0, 0.0D0)
C
 700  ISTATE = -3
      RETURN
C
 800  MSG = 'DLSODA-  Run aborted.. apparent infinite loop.    '
      CALL XERRWD (MSG, 50, 303, 2, 0, 0, 0, 0, 0.0D0, 0.0D0)
      RETURN
C----------------------- End of Subroutine DLSODA ----------------------
      END


*DECK DPRJA
      SUBROUTINE DPRJA(KP,Lfs,Xek,NEQ,Y,YH,NYH,EWT,FTEM,SAVF,WM,IWM,
     1   F,JAC)
      EXTERNAL F, JAC
      INTEGER KP,Lfs, NEQ, NYH, IWM
      DOUBLE PRECISION Xek, Y, YH, EWT, FTEM, SAVF, WM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), EWT(*), FTEM(*), SAVF(*),
     1   WM(*), IWM(*)
      INTEGER IOWND, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER IOWND2, IOWNS2, JTYP, MUSED, MXORDN, MXORDS
      DOUBLE PRECISION ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION ROWND2, ROWNS2, PDNORM
      COMMON /DLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   IOWND(6), IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      COMMON /DLSA01/ ROWND2, ROWNS2(20), PDNORM,
     1   IOWND2(3), IOWNS2(2), JTYP, MUSED, MXORDN, MXORDS
      INTEGER I, I1, I2, IER, II, J, J1, JJ, LENP,
     1   MBA, MBAND, MEB1, MEBAND, ML, ML3, MU, NP1
      DOUBLE PRECISION CON, FAC, HL0, R, R0, SRUR, YI, YJ, YJJ,
     1   DMNORM, DFNORM, DBNORM
C-----------------------------------------------------------------------
C DPRJA is called by DSTODA to compute and process the matrix
C P = I - H*EL(1)*J , where J is an approximation to the Jacobian.
C Here J is computed by the user-supplied routine JAC if
C MITER = 1 or 4 or by finite differencing if MITER = 2 or 5.
C J, scaled by -H*EL(1), is stored in WM.  Then the norm of J (the
C matrix norm consistent with the weighted max-norm on vectors given
C by DMNORM) is computed, and J is overwritten by P.  P is then
C subjected to LU decomposition in preparation for later solution
C of linear systems with P as coefficient matrix.  This is done
C by DGEFA if MITER = 1 or 2, and by DGBFA if MITER = 4 or 5.
C
C In addition to variables described previously, communication
C with DPRJA uses the following:
C Y     = array containing predicted values on entry.
C FTEM  = work array of length N (ACOR in DSTODA).
C SAVF  = array containing f evaluated at predicted y.
C WM    = real work space for matrices.  On output it contains the
C         LU decomposition of P.
C         Storage of matrix elements starts at WM(3).
C         WM also contains the following matrix-related data:
C         WM(1) = SQRT(UROUND), used in numerical Jacobian increments.
C IWM   = integer work space containing pivot information, starting at
C         IWM(21).   IWM also contains the band parameters
C         ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C EL0   = EL(1) (input).
C PDNORM= norm of Jacobian matrix. (Output).
C IERPJ = output error flag,  = 0 if no trouble, .gt. 0 if
C         P matrix found to be singular.
C JCUR  = output flag = 1 to indicate that the Jacobian matrix
C         (or approximation) is now current.
C This routine also uses the Common variables EL0, H, TN, UROUND,
C MITER, N, NFE, and NJE.
C-----------------------------------------------------------------------
      NJE = NJE + 1
      IERPJ = 0
      JCUR = 1
      HL0 = H*EL0
      GO TO (100, 200, 300, 400, 500), MITER
C If MITER = 1, call JAC and multiply by scalar. -----------------------
 100  LENP = N*N
      DO 110 I = 1,LENP
 110    WM(I+2) = 0.0D0
      CALL JAC (KP,Lfs, Xek, NEQ, TN, Y, WM(3), N)
      CON = -HL0
      DO 120 I = 1,LENP
 120    WM(I+2) = WM(I+2)*CON
      GO TO 240
C If MITER = 2, make N calls to F to approximate J. --------------------
 200  FAC = DMNORM (N, SAVF, EWT)
      R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0
      SRUR = WM(1)
      J1 = 2
      DO 230 J = 1,N
        YJ = Y(J)
        R = MAX(SRUR*ABS(YJ),R0/EWT(J))
        Y(J) = Y(J) + R
        FAC = -HL0/R
        CALL F (KP,Lfs, Xek, NEQ, TN, Y, FTEM)
        DO 220 I = 1,N
 220      WM(I+J1) = (FTEM(I) - SAVF(I))*FAC
        Y(J) = YJ
        J1 = J1 + N
 230    CONTINUE
      NFE = NFE + N
 240  CONTINUE
C Compute norm of Jacobian. --------------------------------------------
      PDNORM = DFNORM (N, WM(3), EWT)/ABS(HL0)
C Add identity matrix. -------------------------------------------------
      J = 3
      NP1 = N + 1
      DO 250 I = 1,N
        WM(J) = WM(J) + 1.0D0
 250    J = J + NP1
C Do LU decomposition on P. --------------------------------------------
      CALL DGEFA (WM(3), N, N, IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C Dummy block only, since MITER is never 3 in this routine. ------------
 300  RETURN
C If MITER = 4, call JAC and multiply by scalar. -----------------------
 400  ML = IWM(1)
      MU = IWM(2)
      ML3 = ML + 3
      MBAND = ML + MU + 1
      MEBAND = MBAND + ML
      LENP = MEBAND*N
      DO 410 I = 1,LENP
 410    WM(I+2) = 0.0D0
      CALL JAC (KP,Lfs, Xek, NEQ, TN, Y, WM(ML3), MEBAND)
      CON = -HL0
      DO 420 I = 1,LENP
 420    WM(I+2) = WM(I+2)*CON
      GO TO 570
C If MITER = 5, make MBAND calls to F to approximate J. ----------------
 500  ML = IWM(1)
      MU = IWM(2)
      MBAND = ML + MU + 1
      MBA = MIN(MBAND,N)
      MEBAND = MBAND + ML
      MEB1 = MEBAND - 1
      SRUR = WM(1)
      FAC = DMNORM (N, SAVF, EWT)
      R0 = 1000.0D0*ABS(H)*UROUND*N*FAC
      IF (R0 .EQ. 0.0D0) R0 = 1.0D0
      DO 560 J = 1,MBA
        DO 530 I = J,N,MBAND
          YI = Y(I)
          R = MAX(SRUR*ABS(YI),R0/EWT(I))
 530      Y(I) = Y(I) + R
        CALL F (KP,Lfs, Xek, NEQ, TN, Y, FTEM)
        DO 550 JJ = J,N,MBAND
          Y(JJ) = YH(JJ,1)
          YJJ = Y(JJ)
          R = MAX(SRUR*ABS(YJJ),R0/EWT(JJ))
          FAC = -HL0/R
          I1 = MAX(JJ-MU,1)
          I2 = MIN(JJ+ML,N)
          II = JJ*MEB1 - ML + 2
          DO 540 I = I1,I2
 540        WM(II+I) = (FTEM(I) - SAVF(I))*FAC
 550      CONTINUE
 560    CONTINUE
      NFE = NFE + MBA
 570  CONTINUE
C Compute norm of Jacobian. --------------------------------------------
      PDNORM = DBNORM (N, WM(ML+3), MEBAND, ML, MU, EWT)/ABS(HL0)
C Add identity matrix. -------------------------------------------------
      II = MBAND + 2
      DO 580 I = 1,N
        WM(II) = WM(II) + 1.0D0
 580    II = II + MEBAND
C Do LU decomposition of P. --------------------------------------------
      CALL DGBFA (WM(3), MEBAND, N, ML, MU, IWM(21), IER)
      IF (IER .NE. 0) IERPJ = 1
      RETURN
C----------------------- End of Subroutine DPRJA -----------------------
      END


*DECK DSOLSY
      SUBROUTINE DSOLSY (WM, IWM, X, TEM)
C***BEGIN PROLOGUE  DSOLSY
C***SUBSIDIARY
C***PURPOSE  ODEPACK linear system solver.
C***TYPE      DOUBLE PRECISION (SSOLSY-S, DSOLSY-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This routine manages the solution of the linear system arising from
C  a chord iteration.  It is called if MITER .ne. 0.
C  If MITER is 1 or 2, it calls DGESL to accomplish this.
C  If MITER = 3 it updates the coefficient h*EL0 in the diagonal
C  matrix, and then computes the solution.
C  If MITER is 4 or 5, it calls DGBSL.
C  Communication with DSOLSY uses the following variables:
C  WM    = real work space containing the inverse diagonal matrix if
C          MITER = 3 and the LU decomposition of the matrix otherwise.
C          Storage of matrix elements starts at WM(3).
C          WM also contains the following matrix-related data:
C          WM(1) = SQRT(UROUND) (not used here),
C          WM(2) = HL0, the previous value of h*EL0, used if MITER = 3.
C  IWM   = integer work space containing pivot information, starting at
C          IWM(21), if MITER is 1, 2, 4, or 5.  IWM also contains band
C          parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
C  X     = the right-hand side vector on input, and the solution vector
C          on output, of length N.
C  TEM   = vector of work space of length N, not used in this version.
C  IERSL = output flag (in COMMON).  IERSL = 0 if no trouble occurred.
C          IERSL = 1 if a singular matrix arose with MITER = 3.
C  This routine also uses the COMMON variables EL0, H, MITER, and N.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  DGBSL, DGESL
C***COMMON BLOCKS    DLS001
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C   010418  Reduced size of Common block /DLS001/. (ACH)
C   031105  Restored 'own' variables to Common block /DLS001/, to
C           enable interrupt/restart feature. (ACH)
C***END PROLOGUE  DSOLSY
C**End
      INTEGER IWM
      DOUBLE PRECISION WM, X, TEM
      DIMENSION WM(*), IWM(*), X(*), TEM(*)
      INTEGER IOWND, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      DOUBLE PRECISION ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      COMMON /DLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   IOWND(6), IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, MEBAND, ML, MU
      DOUBLE PRECISION DI, HL0, PHL0, R
C
C***FIRST EXECUTABLE STATEMENT  DSOLSY
      IERSL = 0
      GO TO (100, 100, 300, 400, 400), MITER
 100  CALL DGESL (WM(3), N, N, IWM(21), X, 0)
      RETURN
C
 300  PHL0 = WM(2)
      HL0 = H*EL0
      WM(2) = HL0
      IF (HL0 .EQ. PHL0) GO TO 330
      R = HL0/PHL0
      DO 320 I = 1,N
        DI = 1.0D0 - R*(1.0D0 - 1.0D0/WM(I+2))
        IF (ABS(DI) .EQ. 0.0D0) GO TO 390
 320    WM(I+2) = 1.0D0/DI
 330  DO 340 I = 1,N
 340    X(I) = WM(I+2)*X(I)
      RETURN
 390  IERSL = 1
      RETURN
C
 400  ML = IWM(1)
      MU = IWM(2)
      MEBAND = 2*ML + MU + 1
      CALL DGBSL (WM(3), MEBAND, N, ML, MU, IWM(21), X, 0)
      RETURN
      END

*DECK DSTODA
      SUBROUTINE DSTODA (KP, Lfs, Xek,NEQ,Y,YH,NYH,YH1,EWT,SAVF,ACOR,
     1   WM, IWM, F, JAC, PJAC, SLVS)
      EXTERNAL F, JAC, PJAC, SLVS
      INTEGER NEQ, NYH, IWM, Lfs, KP
      DOUBLE PRECISION Xek, Y, YH, YH1, EWT, SAVF, ACOR, WM
      DIMENSION NEQ(*), Y(*), YH(NYH,*), YH1(*), EWT(*), SAVF(*),
     1   ACOR(*), WM(*), IWM(*)
      INTEGER IOWND, IALTH, IPUP, LMAX, MEO, NQNYH, NSLP,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER IOWND2, ICOUNT, IRFLAG, JTYP, MUSED, MXORDN, MXORDS
      DOUBLE PRECISION CONIT, CRATE, EL, ELCO, HOLD, RMAX, TESCO,
     2   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      DOUBLE PRECISION ROWND2, CM1, CM2, PDEST, PDLAST, RATIO,
     1   PDNORM
      COMMON /DLS001/ CONIT, CRATE, EL(13), ELCO(13,12),
     1   HOLD, RMAX, TESCO(3,12),
     2   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     3   IOWND(6), IALTH, IPUP, LMAX, MEO, NQNYH, NSLP,
     4   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     5   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     6   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      COMMON /DLSA01/ ROWND2, CM1(12), CM2(5), PDEST, PDLAST, RATIO,
     1   PDNORM,
     2   IOWND2(3), ICOUNT, IRFLAG, JTYP, MUSED, MXORDN, MXORDS
      INTEGER I, I1, IREDO, IRET, J, JB, M, NCF, NEWQ
      INTEGER LM1, LM1P1, LM2, LM2P1, NQM1, NQM2
      DOUBLE PRECISION DCON, DDN, DEL, DELP, DSM, DUP, EXDN, EXSM, EXUP,
     1   R, RH, RHDN, RHSM, RHUP, TOLD, DMNORM
      DOUBLE PRECISION ALPHA, DM1,DM2, EXM1,EXM2,
     1   PDH, PNORM, RATE, RH1, RH1IT, RH2, RM, SM1(12)
      SAVE SM1
      DATA SM1/0.5D0, 0.575D0, 0.55D0, 0.45D0, 0.35D0, 0.25D0,
     1   0.20D0, 0.15D0, 0.10D0, 0.075D0, 0.050D0, 0.025D0/
C-----------------------------------------------------------------------
C DSTODA performs one step of the integration of an initial value
C problem for a system of ordinary differential equations.
C Note: DSTODA is independent of the value of the iteration method
C indicator MITER, when this is .ne. 0, and hence is independent
C of the type of chord method used, or the Jacobian structure.
C Communication with DSTODA is done with the following variables:
C
C Y      = an array of length .ge. N used as the Y argument in
C          all calls to F and JAC.
C NEQ    = integer array containing problem size in NEQ(1), and
C          passed as the NEQ argument in all calls to F and JAC.
C YH     = an NYH by LMAX array containing the dependent variables
C          and their approximate scaled derivatives, where
C          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate
C          j-th derivative of y(i), scaled by H**j/factorial(j)
C          (j = 0,1,...,NQ).  On entry for the first step, the first
C          two columns of YH must be set from the initial values.
C NYH    = a constant integer .ge. N, the first dimension of YH.
C YH1    = a one-dimensional array occupying the same space as YH.
C EWT    = an array of length N containing multiplicative weights
C          for local error measurements.  Local errors in y(i) are
C          compared to 1.0/EWT(i) in various error tests.
C SAVF   = an array of working storage, of length N.
C ACOR   = a work array of length N, used for the accumulated
C          corrections.  On a successful return, ACOR(i) contains
C          the estimated one-step local error in y(i).
C WM,IWM = real and integer work arrays associated with matrix
C          operations in chord iteration (MITER .ne. 0).
C PJAC   = name of routine to evaluate and preprocess Jacobian matrix
C          and P = I - H*EL0*Jac, if a chord method is being used.
C          It also returns an estimate of norm(Jac) in PDNORM.
C SLVS   = name of routine to solve linear system in chord iteration.
C CCMAX  = maximum relative change in H*EL0 before PJAC is called.
C H      = the step size to be attempted on the next step.
C          H is altered by the error control algorithm during the
C          problem.  H can be either positive or negative, but its
C          sign must remain constant throughout the problem.
C HMIN   = the minimum absolute value of the step size H to be used.
C HMXI   = inverse of the maximum absolute value of H to be used.
C          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
C          HMIN and HMXI may be changed at any time, but will not
C          take effect until the next change of H is considered.
C TN     = the independent variable. TN is updated on each step taken.
C JSTART = an integer used for input only, with the following
C          values and meanings:
C               0  perform the first step.
C           .gt.0  take a new step continuing from the last.
C              -1  take the next step with a new value of H,
C                    N, METH, MITER, and/or matrix parameters.
C              -2  take the next step with a new value of H,
C                    but with other inputs unchanged.
C          On return, JSTART is set to 1 to facilitate continuation.
C KFLAG  = a completion code with the following meanings:
C               0  the step was succesful.
C              -1  the requested error could not be achieved.
C              -2  corrector convergence could not be achieved.
C              -3  fatal error in PJAC or SLVS.
C          A return with KFLAG = -1 or -2 means either
C          ABS(H) = HMIN or 10 consecutive failures occurred.
C          On a return with KFLAG negative, the values of TN and
C          the YH array are as of the beginning of the last
C          step, and H is the last step size attempted.
C MAXORD = the maximum order of integration method to be allowed.
C MAXCOR = the maximum number of corrector iterations allowed.
C MSBP   = maximum number of steps between PJAC calls (MITER .gt. 0).
C MXNCF  = maximum number of convergence failures allowed.
C METH   = current method.
C          METH = 1 means Adams method (nonstiff)
C          METH = 2 means BDF method (stiff)
C          METH may be reset by DSTODA.
C MITER  = corrector iteration method.
C          MITER = 0 means functional iteration.
C          MITER = JT .gt. 0 means a chord iteration corresponding
C          to Jacobian type JT.  (The DLSODA/DLSODAR argument JT is
C          communicated here as JTYP, but is not used in DSTODA
C          except to load MITER following a method switch.)
C          MITER may be reset by DSTODA.
C N      = the number of first-order differential equations.
C-----------------------------------------------------------------------
      KFLAG = 0
      TOLD = TN
      NCF = 0
      IERPJ = 0
      IERSL = 0
      JCUR = 0
      ICF = 0
      DELP = 0.0D0
      IF (JSTART .GT. 0) GO TO 200
      IF (JSTART .EQ. -1) GO TO 100
      IF (JSTART .EQ. -2) GO TO 160
C-----------------------------------------------------------------------
C On the first call, the order is set to 1, and other variables are
C initialized.  RMAX is the maximum ratio by which H can be increased
C in a single step.  It is initially 1.E4 to compensate for the small
C initial H, but then is normally equal to 10.  If a failure
C occurs (in corrector convergence or error test), RMAX is set at 2
C for the next increase.
C DCFODE is called to get the needed coefficients for both methods.
C-----------------------------------------------------------------------
      LMAX = MAXORD + 1
      NQ = 1
      L = 2
      IALTH = 2
      RMAX = 10000.0D0
      RC = 0.0D0
      EL0 = 1.0D0
      CRATE = 0.7D0
      HOLD = H
      NSLP = 0
      IPUP = MITER
      IRET = 3
C Initialize switching parameters.  METH = 1 is assumed initially. -----
      ICOUNT = 20
      IRFLAG = 0
      PDEST = 0.0D0
      PDLAST = 0.0D0
      RATIO = 5.0D0
      CALL DCFODE (2, ELCO, TESCO)
      DO 10 I = 1,5
 10     CM2(I) = TESCO(2,I)*ELCO(I+1,I)
      CALL DCFODE (1, ELCO, TESCO)
      DO 20 I = 1,12
 20     CM1(I) = TESCO(2,I)*ELCO(I+1,I)
      GO TO 150
C-----------------------------------------------------------------------
C The following block handles preliminaries needed when JSTART = -1.
C IPUP is set to MITER to force a matrix update.
C If an order increase is about to be considered (IALTH = 1),
C IALTH is reset to 2 to postpone consideration one more step.
C If the caller has changed METH, DCFODE is called to reset
C the coefficients of the method.
C If H is to be changed, YH must be rescaled.
C If H or METH is being changed, IALTH is reset to L = NQ + 1
C to prevent further changes in H for that many steps.
C-----------------------------------------------------------------------
 100  IPUP = MITER
      LMAX = MAXORD + 1
      IF (IALTH .EQ. 1) IALTH = 2
      IF (METH .EQ. MUSED) GO TO 160
      CALL DCFODE (METH, ELCO, TESCO)
      IALTH = L
      IRET = 1
C-----------------------------------------------------------------------
C The el vector and related constants are reset
C whenever the order NQ is changed, or at the start of the problem.
C-----------------------------------------------------------------------
 150  DO 155 I = 1,L
 155    EL(I) = ELCO(I,NQ)
      NQNYH = NQ*NYH
      RC = RC*EL(1)/EL0
      EL0 = EL(1)
      CONIT = 0.5D0/(NQ+2)
      GO TO (160, 170, 200), IRET
C-----------------------------------------------------------------------
C If H is being changed, the H ratio RH is checked against
C RMAX, HMIN, and HMXI, and the YH array rescaled.  IALTH is set to
C L = NQ + 1 to prevent a change of H for that many steps, unless
C forced by a convergence or error test failure.
C-----------------------------------------------------------------------
 160  IF (H .EQ. HOLD) GO TO 200
      RH = H/HOLD
      H = HOLD
      IREDO = 3
      GO TO 175
 170  RH = MAX(RH,HMIN/ABS(H))
 175  RH = MIN(RH,RMAX)
      RH = RH/MAX(1.0D0,ABS(H)*HMXI*RH)
C-----------------------------------------------------------------------
C If METH = 1, also restrict the new step size by the stability region.
C If this reduces H, set IRFLAG to 1 so that if there are roundoff
C problems later, we can assume that is the cause of the trouble.
C-----------------------------------------------------------------------
      IF (METH .EQ. 2) GO TO 178
      IRFLAG = 0
      PDH = MAX(ABS(H)*PDLAST,0.000001D0)
      IF (RH*PDH*1.00001D0 .LT. SM1(NQ)) GO TO 178
      RH = SM1(NQ)/PDH
      IRFLAG = 1
 178  CONTINUE
      R = 1.0D0
      DO 180 J = 2,L
        R = R*RH
        DO 180 I = 1,N
 180      YH(I,J) = YH(I,J)*R
      H = H*RH
      RC = RC*RH
      IALTH = L
      IF (IREDO .EQ. 0) GO TO 690
C-----------------------------------------------------------------------
C This section computes the predicted values by effectively
C multiplying the YH array by the Pascal triangle matrix.
C RC is the ratio of new to old values of the coefficient  H*EL(1).
C When RC differs from 1 by more than CCMAX, IPUP is set to MITER
C to force PJAC to be called, if a Jacobian is involved.
C In any case, PJAC is called at least every MSBP steps.
C-----------------------------------------------------------------------
 200  IF (ABS(RC-1.0D0) .GT. CCMAX) IPUP = MITER
      IF (NST .GE. NSLP+MSBP) IPUP = MITER
      TN = TN + H
      I1 = NQNYH + 1
      DO 215 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 210 I = I1,NQNYH
 210      YH1(I) = YH1(I) + YH1(I+NYH)
 215    CONTINUE
      PNORM = DMNORM (N, YH1, EWT)
C-----------------------------------------------------------------------
C Up to MAXCOR corrector iterations are taken.  A convergence test is
C made on the RMS-norm of each correction, weighted by the error
C weight vector EWT.  The sum of the corrections is accumulated in the
C vector ACOR(i).  The YH array is not altered in the corrector loop.
C-----------------------------------------------------------------------
 220  M = 0
      RATE = 0.0D0
      DEL = 0.0D0
      DO 230 I = 1,N
 230    Y(I) = YH(I,1)
      CALL F (KP, Lfs, Xek, NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      IF (IPUP .LE. 0) GO TO 250
C-----------------------------------------------------------------------
C If indicated, the matrix P = I - H*EL(1)*J is reevaluated and
C preprocessed before starting the corrector iteration.  IPUP is set
C to 0 as an indicator that this has been done.
C-----------------------------------------------------------------------
      CALL PJAC (KP, Lfs, Xek,NEQ, Y, YH, NYH, EWT, ACOR, SAVF, WM, IWM,
     1  F, JAC)
      IPUP = 0
      RC = 1.0D0
      NSLP = NST
      CRATE = 0.7D0
      IF (IERPJ .NE. 0) GO TO 430
 250  DO 260 I = 1,N
 260    ACOR(I) = 0.0D0
 270  IF (MITER .NE. 0) GO TO 350
C-----------------------------------------------------------------------
C In the case of functional iteration, update Y directly from
C the result of the last function evaluation.
C-----------------------------------------------------------------------
      DO 290 I = 1,N
        SAVF(I) = H*SAVF(I) - YH(I,2)
 290    Y(I) = SAVF(I) - ACOR(I)
      DEL = DMNORM (N, Y, EWT)
      DO 300 I = 1,N
        Y(I) = YH(I,1) + EL(1)*SAVF(I)
 300    ACOR(I) = SAVF(I)
      GO TO 400
C-----------------------------------------------------------------------
C In the case of the chord method, compute the corrector error,
C and solve the linear system with that as right-hand side and
C P as coefficient matrix.
C-----------------------------------------------------------------------
 350  DO 360 I = 1,N
 360    Y(I) = H*SAVF(I) - (YH(I,2) + ACOR(I))
      CALL SLVS (WM, IWM, Y, SAVF)
      IF (IERSL .LT. 0) GO TO 430
      IF (IERSL .GT. 0) GO TO 410
      DEL = DMNORM (N, Y, EWT)
      DO 380 I = 1,N
        ACOR(I) = ACOR(I) + Y(I)
 380    Y(I) = YH(I,1) + EL(1)*ACOR(I)
C-----------------------------------------------------------------------
C Test for convergence.  If M .gt. 0, an estimate of the convergence
C rate constant is stored in CRATE, and this is used in the test.
C
C We first check for a change of iterates that is the size of
C roundoff error.  If this occurs, the iteration has converged, and a
C new rate estimate is not formed.
C In all other cases, force at least two iterations to estimate a
C local Lipschitz constant estimate for Adams methods.
C On convergence, form PDEST = local maximum Lipschitz constant
C estimate.  PDLAST is the most recent nonzero estimate.
C-----------------------------------------------------------------------
 400  CONTINUE
      IF (DEL .LE. 100.0D0*PNORM*UROUND) GO TO 450
      IF (M .EQ. 0 .AND. METH .EQ. 1) GO TO 405
      IF (M .EQ. 0) GO TO 402
      RM = 1024.0D0
      IF (DEL .LE. 1024.0D0*DELP) RM = DEL/DELP
      RATE = MAX(RATE,RM)
      CRATE = MAX(0.2D0*CRATE,RM)
 402  DCON = DEL*MIN(1.0D0,1.5D0*CRATE)/(TESCO(2,NQ)*CONIT)
      IF (DCON .GT. 1.0D0) GO TO 405
      PDEST = MAX(PDEST,RATE/ABS(H*EL(1)))
      IF (PDEST .NE. 0.0D0) PDLAST = PDEST
      GO TO 450
 405  CONTINUE
      M = M + 1
      IF (M .EQ. MAXCOR) GO TO 410
      IF (M .GE. 2 .AND. DEL .GT. 2.0D0*DELP) GO TO 410
      DELP = DEL
      CALL F (KP,Lfs, Xek, NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      GO TO 270
C-----------------------------------------------------------------------
C The corrector iteration failed to converge.
C If MITER .ne. 0 and the Jacobian is out of date, PJAC is called for
C the next try.  Otherwise the YH array is retracted to its values
C before prediction, and H is reduced, if possible.  If H cannot be
C reduced or MXNCF failures have occurred, exit with KFLAG = -2.
C-----------------------------------------------------------------------
 410  IF (MITER .EQ. 0 .OR. JCUR .EQ. 1) GO TO 430
      ICF = 1
      IPUP = MITER
      GO TO 220
 430  ICF = 2
      NCF = NCF + 1
      RMAX = 2.0D0
      TN = TOLD
      I1 = NQNYH + 1
      DO 445 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 440 I = I1,NQNYH
 440      YH1(I) = YH1(I) - YH1(I+NYH)
 445    CONTINUE
      IF (IERPJ .LT. 0 .OR. IERSL .LT. 0) GO TO 680
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 670
      IF (NCF .EQ. MXNCF) GO TO 670
      RH = 0.25D0
      IPUP = MITER
      IREDO = 1
      GO TO 170
C-----------------------------------------------------------------------
C The corrector has converged.  JCUR is set to 0
C to signal that the Jacobian involved may need updating later.
C The local error test is made and control passes to statement 500
C if it fails.
C-----------------------------------------------------------------------
 450  JCUR = 0
      IF (M .EQ. 0) DSM = DEL/TESCO(2,NQ)
      IF (M .GT. 0) DSM = DMNORM (N, ACOR, EWT)/TESCO(2,NQ)
      IF (DSM .GT. 1.0D0) GO TO 500
C-----------------------------------------------------------------------
C After a successful step, update the YH array.
C Decrease ICOUNT by 1, and if it is -1, consider switching methods.
C If a method switch is made, reset various parameters,
C rescale the YH array, and exit.  If there is no switch,
C consider changing H if IALTH = 1.  Otherwise decrease IALTH by 1.
C If IALTH is then 1 and NQ .lt. MAXORD, then ACOR is saved for
C use in a possible order increase on the next step.
C If a change in H is considered, an increase or decrease in order
C by one is considered also.  A change in H is made only if it is by a
C factor of at least 1.1.  If not, IALTH is set to 3 to prevent
C testing for that many steps.
C-----------------------------------------------------------------------
      KFLAG = 0
      IREDO = 0
      NST = NST + 1
      HU = H
      NQU = NQ
      MUSED = METH
      DO 460 J = 1,L
        DO 460 I = 1,N
 460      YH(I,J) = YH(I,J) + EL(J)*ACOR(I)
      ICOUNT = ICOUNT - 1
      IF (ICOUNT .GE. 0) GO TO 488
      IF (METH .EQ. 2) GO TO 480
C-----------------------------------------------------------------------
C We are currently using an Adams method.  Consider switching to BDF.
C If the current order is greater than 5, assume the problem is
C not stiff, and skip this section.
C If the Lipschitz constant and error estimate are not polluted
C by roundoff, go to 470 and perform the usual test.
C Otherwise, switch to the BDF methods if the last step was
C restricted to insure stability (irflag = 1), and stay with Adams
C method if not.  When switching to BDF with polluted error estimates,
C in the absence of other information, double the step size.
C
C When the estimates are OK, we make the usual test by computing
C the step size we could have (ideally) used on this step,
C with the current (Adams) method, and also that for the BDF.
C If NQ .gt. MXORDS, we consider changing to order MXORDS on switching.
C Compare the two step sizes to decide whether to switch.
C The step size advantage must be at least RATIO = 5 to switch.
C-----------------------------------------------------------------------
      IF (NQ .GT. 5) GO TO 488
      IF (DSM .GT. 100.0D0*PNORM*UROUND .AND. PDEST .NE. 0.0D0)
     1   GO TO 470
      IF (IRFLAG .EQ. 0) GO TO 488
      RH2 = 2.0D0
      NQM2 = MIN(NQ,MXORDS)
      GO TO 478
 470  CONTINUE
      EXSM = 1.0D0/L
      RH1 = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      RH1IT = 2.0D0*RH1
      PDH = PDLAST*ABS(H)
      IF (PDH*RH1 .GT. 0.00001D0) RH1IT = SM1(NQ)/PDH
      RH1 = MIN(RH1,RH1IT)
      IF (NQ .LE. MXORDS) GO TO 474
         NQM2 = MXORDS
         LM2 = MXORDS + 1
         EXM2 = 1.0D0/LM2
         LM2P1 = LM2 + 1
         DM2 = DMNORM (N, YH(1,LM2P1), EWT)/CM2(MXORDS)
         RH2 = 1.0D0/(1.2D0*DM2**EXM2 + 0.0000012D0)
         GO TO 476
 474  DM2 = DSM*(CM1(NQ)/CM2(NQ))
      RH2 = 1.0D0/(1.2D0*DM2**EXSM + 0.0000012D0)
      NQM2 = NQ
 476  CONTINUE
      IF (RH2 .LT. RATIO*RH1) GO TO 488
C THE SWITCH TEST PASSED.  RESET RELEVANT QUANTITIES FOR BDF. ----------
 478  RH = RH2
      ICOUNT = 20
      METH = 2
      MITER = JTYP
      PDLAST = 0.0D0
      NQ = NQM2
      L = NQ + 1
      GO TO 170
C-----------------------------------------------------------------------
C We are currently using a BDF method.  Consider switching to Adams.
C Compute the step size we could have (ideally) used on this step,
C with the current (BDF) method, and also that for the Adams.
C If NQ .gt. MXORDN, we consider changing to order MXORDN on switching.
C Compare the two step sizes to decide whether to switch.
C The step size advantage must be at least 5/RATIO = 1 to switch.
C If the step size for Adams would be so small as to cause
C roundoff pollution, we stay with BDF.
C-----------------------------------------------------------------------
 480  CONTINUE
      EXSM = 1.0D0/L
      IF (MXORDN .GE. NQ) GO TO 484
         NQM1 = MXORDN
         LM1 = MXORDN + 1
         EXM1 = 1.0D0/LM1
         LM1P1 = LM1 + 1
         DM1 = DMNORM (N, YH(1,LM1P1), EWT)/CM1(MXORDN)
         RH1 = 1.0D0/(1.2D0*DM1**EXM1 + 0.0000012D0)
         GO TO 486
 484  DM1 = DSM*(CM2(NQ)/CM1(NQ))
      RH1 = 1.0D0/(1.2D0*DM1**EXSM + 0.0000012D0)
      NQM1 = NQ
      EXM1 = EXSM
 486  RH1IT = 2.0D0*RH1
      PDH = PDNORM*ABS(H)
      IF (PDH*RH1 .GT. 0.00001D0) RH1IT = SM1(NQM1)/PDH
      RH1 = MIN(RH1,RH1IT)
      RH2 = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      IF (RH1*RATIO .LT. 5.0D0*RH2) GO TO 488
      ALPHA = MAX(0.001D0,RH1)
      DM1 = (ALPHA**EXM1)*DM1
      IF (DM1 .LE. 1000.0D0*UROUND*PNORM) GO TO 488
C The switch test passed.  Reset relevant quantities for Adams. --------
      RH = RH1
      ICOUNT = 20
      METH = 1
      MITER = 0
      PDLAST = 0.0D0
      NQ = NQM1
      L = NQ + 1
      GO TO 170
C
C No method switch is being made.  Do the usual step/order selection. --
 488  CONTINUE
      IALTH = IALTH - 1
      IF (IALTH .EQ. 0) GO TO 520
      IF (IALTH .GT. 1) GO TO 700
      IF (L .EQ. LMAX) GO TO 700
      DO 490 I = 1,N
 490    YH(I,LMAX) = ACOR(I)
      GO TO 700
C-----------------------------------------------------------------------
C The error test failed.  KFLAG keeps track of multiple failures.
C Restore TN and the YH array to their previous values, and prepare
C to try the step again.  Compute the optimum step size for this or
C one lower order.  After 2 or more failures, H is forced to decrease
C by a factor of 0.2 or less.
C-----------------------------------------------------------------------
 500  KFLAG = KFLAG - 1
      TN = TOLD
      I1 = NQNYH + 1
      DO 515 JB = 1,NQ
        I1 = I1 - NYH
CDIR$ IVDEP
        DO 510 I = I1,NQNYH
 510      YH1(I) = YH1(I) - YH1(I+NYH)
 515    CONTINUE
      RMAX = 2.0D0
      IF (ABS(H) .LE. HMIN*1.00001D0) GO TO 660
      IF (KFLAG .LE. -3) GO TO 640
      IREDO = 2
      RHUP = 0.0D0
      GO TO 540
C-----------------------------------------------------------------------
C Regardless of the success or failure of the step, factors
C RHDN, RHSM, and RHUP are computed, by which H could be multiplied
C at order NQ - 1, order NQ, or order NQ + 1, respectively.
C In the case of failure, RHUP = 0.0 to avoid an order increase.
C The largest of these is determined and the new order chosen
C accordingly.  If the order is to be increased, we compute one
C additional scaled derivative.
C-----------------------------------------------------------------------
 520  RHUP = 0.0D0
      IF (L .EQ. LMAX) GO TO 540
      DO 530 I = 1,N
 530    SAVF(I) = ACOR(I) - YH(I,LMAX)
      DUP = DMNORM (N, SAVF, EWT)/TESCO(3,NQ)
      EXUP = 1.0D0/(L+1)
      RHUP = 1.0D0/(1.4D0*DUP**EXUP + 0.0000014D0)
 540  EXSM = 1.0D0/L
      RHSM = 1.0D0/(1.2D0*DSM**EXSM + 0.0000012D0)
      RHDN = 0.0D0
      IF (NQ .EQ. 1) GO TO 550
      DDN = DMNORM (N, YH(1,L), EWT)/TESCO(1,NQ)
      EXDN = 1.0D0/NQ
      RHDN = 1.0D0/(1.3D0*DDN**EXDN + 0.0000013D0)
C If METH = 1, limit RH according to the stability region also. --------
 550  IF (METH .EQ. 2) GO TO 560
      PDH = MAX(ABS(H)*PDLAST,0.000001D0)
      IF (L .LT. LMAX) RHUP = MIN(RHUP,SM1(L)/PDH)
      RHSM = MIN(RHSM,SM1(NQ)/PDH)
      IF (NQ .GT. 1) RHDN = MIN(RHDN,SM1(NQ-1)/PDH)
      PDEST = 0.0D0
 560  IF (RHSM .GE. RHUP) GO TO 570
      IF (RHUP .GT. RHDN) GO TO 590
      GO TO 580
 570  IF (RHSM .LT. RHDN) GO TO 580
      NEWQ = NQ
      RH = RHSM
      GO TO 620
 580  NEWQ = NQ - 1
      RH = RHDN
      IF (KFLAG .LT. 0 .AND. RH .GT. 1.0D0) RH = 1.0D0
      GO TO 620
 590  NEWQ = L
      RH = RHUP
      IF (RH .LT. 1.1D0) GO TO 610
      R = EL(L)/L
      DO 600 I = 1,N
 600    YH(I,NEWQ+1) = ACOR(I)*R
      GO TO 630
 610  IALTH = 3
      GO TO 700
C If METH = 1 and H is restricted by stability, bypass 10 percent test.
 620  IF (METH .EQ. 2) GO TO 622
      IF (RH*PDH*1.00001D0 .GE. SM1(NEWQ)) GO TO 625
 622  IF (KFLAG .EQ. 0 .AND. RH .LT. 1.1D0) GO TO 610
 625  IF (KFLAG .LE. -2) RH = MIN(RH,0.2D0)
C-----------------------------------------------------------------------
C If there is a change of order, reset NQ, L, and the coefficients.
C In any case H is reset according to RH and the YH array is rescaled.
C Then exit from 690 if the step was OK, or redo the step otherwise.
C-----------------------------------------------------------------------
      IF (NEWQ .EQ. NQ) GO TO 170
 630  NQ = NEWQ
      L = NQ + 1
      IRET = 2
      GO TO 150
C-----------------------------------------------------------------------
C Control reaches this section if 3 or more failures have occured.
C If 10 failures have occurred, exit with KFLAG = -1.
C It is assumed that the derivatives that have accumulated in the
C YH array have errors of the wrong order.  Hence the first
C derivative is recomputed, and the order is set to 1.  Then
C H is reduced by a factor of 10, and the step is retried,
C until it succeeds or H reaches HMIN.
C-----------------------------------------------------------------------
 640  IF (KFLAG .EQ. -10) GO TO 660
      RH = 0.1D0
      RH = MAX(HMIN/ABS(H),RH)
      H = H*RH
      DO 645 I = 1,N
 645    Y(I) = YH(I,1)
      CALL F (KP, Lfs, Xek, NEQ, TN, Y, SAVF)
      NFE = NFE + 1
      DO 650 I = 1,N
 650    YH(I,2) = H*SAVF(I)
      IPUP = MITER
      IALTH = 5
      IF (NQ .EQ. 1) GO TO 200
      NQ = 1
      L = 2
      IRET = 3
      GO TO 150
C-----------------------------------------------------------------------
C All returns are made through this section.  H is saved in HOLD
C to allow the caller to change H on the next step.
C-----------------------------------------------------------------------
 660  KFLAG = -1
      GO TO 720
 670  KFLAG = -2
      GO TO 720
 680  KFLAG = -3
      GO TO 720
 690  RMAX = 10.0D0
 700  R = 1.0D0/TESCO(2,NQU)
      DO 710 I = 1,N
 710    ACOR(I) = ACOR(I)*R
 720  HOLD = H
      JSTART = 1
      RETURN
C----------------------- End of Subroutine DSTODA ----------------------
      END

*DECK JGROUP
      SUBROUTINE JGROUP (N,IA,JA,MAXG,NGRP,IGP,JGP,INCL,JDONE,IER)
      INTEGER N, IA, JA, MAXG, NGRP, IGP, JGP, INCL, JDONE, IER
      DIMENSION IA(*), JA(*), IGP(*), JGP(*), INCL(*), JDONE(*)
C-----------------------------------------------------------------------
C This subroutine constructs groupings of the column indices of
C the Jacobian matrix, used in the numerical evaluation of the
C Jacobian by finite differences.
C
C Input:
C N      = the order of the matrix.
C IA,JA  = sparse structure descriptors of the matrix by rows.
C MAXG   = length of available storage in the IGP array.
C
C Output:
C NGRP   = number of groups.
C JGP    = array of length N containing the column indices by groups.
C IGP    = pointer array of length NGRP + 1 to the locations in JGP
C          of the beginning of each group.
C IER    = error indicator.  IER = 0 if no error occurred, or 1 if
C          MAXG was insufficient.
C
C INCL and JDONE are working arrays of length N.
C-----------------------------------------------------------------------
      INTEGER I, J, K, KMIN, KMAX, NCOL, NG
C
      IER = 0
      DO 10 J = 1,N
 10     JDONE(J) = 0
      NCOL = 1
      DO 60 NG = 1,MAXG
        IGP(NG) = NCOL
        DO 20 I = 1,N
 20       INCL(I) = 0
        DO 50 J = 1,N
C Reject column J if it is already in a group.--------------------------
          IF (JDONE(J) .EQ. 1) GO TO 50
          KMIN = IA(J)
          KMAX = IA(J+1) - 1
          DO 30 K = KMIN,KMAX
C Reject column J if it overlaps any column already in this group.------
            I = JA(K)
            IF (INCL(I) .EQ. 1) GO TO 50
 30         CONTINUE
C Accept column J into group NG.----------------------------------------
          JGP(NCOL) = J
          NCOL = NCOL + 1
          JDONE(J) = 1
          DO 40 K = KMIN,KMAX
            I = JA(K)
 40         INCL(I) = 1
 50       CONTINUE
C Stop if this group is empty (grouping is complete).-------------------
        IF (NCOL .EQ. IGP(NG)) GO TO 70
 60     CONTINUE
C Error return if not all columns were chosen (MAXG too small).---------
      IF (NCOL .LE. N) GO TO 80
      NG = MAXG
 70   NGRP = NG - 1
      RETURN
 80   IER = 1
      RETURN
C----------------------- End of Subroutine JGROUP ----------------------
      END


*DECK DUMACH
      DOUBLE PRECISION FUNCTION DUMACH ()
C***BEGIN PROLOGUE  DUMACH
C***PURPOSE  Compute the unit roundoff of the machine.
C***CATEGORY  R1
C***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C *Usage:
C        DOUBLE PRECISION  A, DUMACH
C        A = DUMACH()
C
C *Function Return Values:
C     A : the unit roundoff of the machine.
C
C *Description:
C     The unit roundoff is defined as the smallest positive machine
C     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH
C     in a machine-independent manner.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DUMSUM
C***REVISION HISTORY  (YYYYMMDD)
C   19930216  DATE WRITTEN
C   19930818  Added SLATEC-format prologue.  (FNF)
C   20030707  Added DUMSUM to force normal storage of COMP.  (ACH)
C***END PROLOGUE  DUMACH
C
      DOUBLE PRECISION U, COMP
C***FIRST EXECUTABLE STATEMENT  DUMACH
      U = 1.0D0
 10   U = U*0.5D0
      CALL DUMSUM(1.0D0, U, COMP)
      IF (COMP .NE. 1.0D0) GO TO 10
      DUMACH = U*2.0D0
      RETURN
C----------------------- End of Function DUMACH ------------------------
      END
      SUBROUTINE DUMSUM(A,B,C)
C     Routine to force normal storing of A + B, for DUMACH.
      DOUBLE PRECISION A, B, C
      C = A + B
      RETURN
      END
*DECK DCFODE
      SUBROUTINE DCFODE (METH, ELCO, TESCO)
C***BEGIN PROLOGUE  DCFODE
C***SUBSIDIARY
C***PURPOSE  Set ODE integrator coefficients.
C***TYPE      DOUBLE PRECISION (SCFODE-S, DCFODE-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  DCFODE is called by the integrator routine to set coefficients
C  needed there.  The coefficients for the current method, as
C  given by the value of METH, are set for all orders and saved.
C  The maximum order assumed here is 12 if METH = 1 and 5 if METH = 2.
C  (A smaller value of the maximum order is also allowed.)
C  DCFODE is called once at the beginning of the problem,
C  and is not called again unless and until METH is changed.
C
C  The ELCO array contains the basic method coefficients.
C  The coefficients el(i), 1 .le. i .le. nq+1, for the method of
C  order nq are stored in ELCO(i,nq).  They are given by a genetrating
C  polynomial, i.e.,
C      l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
C  For the implicit Adams methods, l(x) is given by
C      dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) = 0.
C  For the BDF methods, l(x) is given by
C      l(x) = (x+1)*(x+2)* ... *(x+nq)/K,
C  where         K = factorial(nq)*(1 + 1/2 + ... + 1/nq).
C
C  The TESCO array contains test constants used for the
C  local error test and the selection of step size and/or order.
C  At order nq, TESCO(k,nq) is used for the selection of step
C  size at order nq - 1 if k = 1, at order nq if k = 2, and at order
C  nq + 1 if k = 3.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C***END PROLOGUE  DCFODE
C**End
      INTEGER METH
      INTEGER I, IB, NQ, NQM1, NQP1
      DOUBLE PRECISION ELCO, TESCO
      DOUBLE PRECISION AGAMQ, FNQ, FNQM1, PC, PINT, RAGQ,
     1   RQFAC, RQ1FAC, TSIGN, XPIN
      DIMENSION ELCO(13,12), TESCO(3,12)
      DIMENSION PC(12)
C
C***FIRST EXECUTABLE STATEMENT  DCFODE
      GO TO (100, 200), METH
C
 100  ELCO(1,1) = 1.0D0
      ELCO(2,1) = 1.0D0
      TESCO(1,1) = 0.0D0
      TESCO(2,1) = 2.0D0
      TESCO(1,2) = 1.0D0
      TESCO(3,12) = 0.0D0
      PC(1) = 1.0D0
      RQFAC = 1.0D0
      DO 140 NQ = 2,12
C-----------------------------------------------------------------------
C The PC array will contain the coefficients of the polynomial
C     p(x) = (x+1)*(x+2)*...*(x+nq-1).
C Initially, p(x) = 1.
C-----------------------------------------------------------------------
        RQ1FAC = RQFAC
        RQFAC = RQFAC/NQ
        NQM1 = NQ - 1
        FNQM1 = NQM1
        NQP1 = NQ + 1
C Form coefficients of p(x)*(x+nq-1). ----------------------------------
        PC(NQ) = 0.0D0
        DO 110 IB = 1,NQM1
          I = NQP1 - IB
 110      PC(I) = PC(I-1) + FNQM1*PC(I)
        PC(1) = FNQM1*PC(1)
C Compute integral, -1 to 0, of p(x) and x*p(x). -----------------------
        PINT = PC(1)
        XPIN = PC(1)/2.0D0
        TSIGN = 1.0D0
        DO 120 I = 2,NQ
          TSIGN = -TSIGN
          PINT = PINT + TSIGN*PC(I)/I
 120      XPIN = XPIN + TSIGN*PC(I)/(I+1)
C Store coefficients in ELCO and TESCO. --------------------------------
        ELCO(1,NQ) = PINT*RQ1FAC
        ELCO(2,NQ) = 1.0D0
        DO 130 I = 2,NQ
 130      ELCO(I+1,NQ) = RQ1FAC*PC(I)/I
        AGAMQ = RQFAC*XPIN
        RAGQ = 1.0D0/AGAMQ
        TESCO(2,NQ) = RAGQ
        IF (NQ .LT. 12) TESCO(1,NQP1) = RAGQ*RQFAC/NQP1
        TESCO(3,NQM1) = RAGQ
 140    CONTINUE
      RETURN
C
 200  PC(1) = 1.0D0
      RQ1FAC = 1.0D0
      DO 230 NQ = 1,5
C-----------------------------------------------------------------------
C The PC array will contain the coefficients of the polynomial
C     p(x) = (x+1)*(x+2)*...*(x+nq).
C Initially, p(x) = 1.
C-----------------------------------------------------------------------
        FNQ = NQ
        NQP1 = NQ + 1
C Form coefficients of p(x)*(x+nq). ------------------------------------
        PC(NQP1) = 0.0D0
        DO 210 IB = 1,NQ
          I = NQ + 2 - IB
 210      PC(I) = PC(I-1) + FNQ*PC(I)
        PC(1) = FNQ*PC(1)
C Store coefficients in ELCO and TESCO. --------------------------------
        DO 220 I = 1,NQP1
 220      ELCO(I,NQ) = PC(I)/PC(2)
        ELCO(2,NQ) = 1.0D0
        TESCO(1,NQ) = RQ1FAC
        TESCO(2,NQ) = NQP1/ELCO(1,NQ)
        TESCO(3,NQ) = (NQ+2)/ELCO(1,NQ)
        RQ1FAC = RQ1FAC/FNQ
 230    CONTINUE
      RETURN
C----------------------- END OF SUBROUTINE DCFODE ----------------------
      END
*DECK DINTDY
      SUBROUTINE DINTDY (T, K, YH, NYH, DKY, IFLAG)
C***BEGIN PROLOGUE  DINTDY
C***SUBSIDIARY
C***PURPOSE  Interpolate solution derivatives.
C***TYPE      DOUBLE PRECISION (SINTDY-S, DINTDY-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  DINTDY computes interpolated values of the K-th derivative of the
C  dependent variable vector y, and stores it in DKY.  This routine
C  is called within the package with K = 0 and T = TOUT, but may
C  also be called by the user for any K up to the current order.
C  (See detailed instructions in the usage documentation.)
C
C  The computed values in DKY are gotten by interpolation using the
C  Nordsieck history array YH.  This array corresponds uniquely to a
C  vector-valued polynomial of degree NQCUR or less, and DKY is set
C  to the K-th derivative of this polynomial at T.
C  The formula for DKY is:
C               q
C   DKY(i)  =  sum  c(j,K) * (T - tn)**(j-K) * h**(-j) * YH(i,j+1)
C              j=K
C  where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, tn = TCUR, h = HCUR.
C  The quantities  nq = NQCUR, l = nq+1, N = NEQ, tn, and h are
C  communicated by COMMON.  The above sum is done in reverse order.
C  IFLAG is returned negative if either K or T is out of bounds.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  XERRWD
C***COMMON BLOCKS    DLS001
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C   010418  Reduced size of Common block /DLS001/. (ACH)
C   031105  Restored 'own' variables to Common block /DLS001/, to
C           enable interrupt/restart feature. (ACH)
C   050427  Corrected roundoff decrement in TP. (ACH)
C***END PROLOGUE  DINTDY
C**End
      INTEGER K, NYH, IFLAG
      DOUBLE PRECISION T, YH, DKY
      DIMENSION YH(NYH,*), DKY(*)
      INTEGER IOWND, IOWNS,
     1   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     2   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     3   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      DOUBLE PRECISION ROWNS,
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND
      COMMON /DLS001/ ROWNS(209),
     1   CCMAX, EL0, H, HMIN, HMXI, HU, RC, TN, UROUND,
     2   IOWND(6), IOWNS(6),
     3   ICF, IERPJ, IERSL, JCUR, JSTART, KFLAG, L,
     4   LYH, LEWT, LACOR, LSAVF, LWM, LIWM, METH, MITER,
     5   MAXORD, MAXCOR, MSBP, MXNCF, N, NQ, NST, NFE, NJE, NQU
      INTEGER I, IC, J, JB, JB2, JJ, JJ1, JP1
      DOUBLE PRECISION C, R, S, TP
      CHARACTER*80 MSG
C
C***FIRST EXECUTABLE STATEMENT  DINTDY
      IFLAG = 0
      IF (K .LT. 0 .OR. K .GT. NQ) GO TO 80
      TP = TN - HU -  100.0D0*UROUND*SIGN(ABS(TN) + ABS(HU), HU)
      IF ((T-TP)*(T-TN) .GT. 0.0D0) GO TO 90
C
      S = (T - TN)/H
      IC = 1
      IF (K .EQ. 0) GO TO 15
      JJ1 = L - K
      DO 10 JJ = JJ1,NQ
 10     IC = IC*JJ
 15   C = IC
      DO 20 I = 1,N
 20     DKY(I) = C*YH(I,L)
      IF (K .EQ. NQ) GO TO 55
      JB2 = NQ - K
      DO 50 JB = 1,JB2
        J = NQ - JB
        JP1 = J + 1
        IC = 1
        IF (K .EQ. 0) GO TO 35
        JJ1 = JP1 - K
        DO 30 JJ = JJ1,J
 30       IC = IC*JJ
 35     C = IC
        DO 40 I = 1,N
 40       DKY(I) = C*YH(I,JP1) + S*DKY(I)
 50     CONTINUE
      IF (K .EQ. 0) RETURN
 55   R = H**(-K)
      DO 60 I = 1,N
 60     DKY(I) = R*DKY(I)
      RETURN
C
 80   MSG = 'DINTDY-  K (=I1) illegal      '
      CALL XERRWD (MSG, 30, 51, 0, 1, K, 0, 0, 0.0D0, 0.0D0)
      IFLAG = -1
      RETURN
 90   MSG = 'DINTDY-  T (=R1) illegal      '
      CALL XERRWD (MSG, 30, 52, 0, 0, 0, 0, 1, T, 0.0D0)
      MSG='      T not in interval TCUR - HU (= R1) to TCUR (=R2)      '
      CALL XERRWD (MSG, 60, 52, 0, 0, 0, 0, 2, TP, TN)
      IFLAG = -2
      RETURN
C----------------------- END OF SUBROUTINE DINTDY ----------------------
      END
*DECK DGEFA
      SUBROUTINE DGEFA (A, LDA, N, IPVT, INFO)
C***BEGIN PROLOGUE  DGEFA
C***PURPOSE  Factor a matrix using Gaussian elimination.
C***CATEGORY  D2A1
C***TYPE      DOUBLE PRECISION (SGEFA-S, DGEFA-D, CGEFA-C)
C***KEYWORDS  GENERAL MATRIX, LINEAR ALGEBRA, LINPACK,
C             MATRIX FACTORIZATION
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGEFA factors a double precision matrix by Gaussian elimination.
C
C     DGEFA is usually called by DGECO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA) .
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C     On Return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGESL or DGEDI will divide by zero
C                     if called.  Use  RCOND  in DGECO for a reliable
C                     indication of singularity.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGEFA
      INTEGER LDA,N,IPVT(*),INFO
      DOUBLE PRECISION A(LDA,*)
C
      DOUBLE PRECISION T
      INTEGER IDAMAX,J,K,KP1,L,NM1
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
C***FIRST EXECUTABLE STATEMENT  DGEFA
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        FIND L = PIVOT INDEX
C
         L = IDAMAX(N-K+1,A(K,K),1) + K - 1
         IPVT(K) = L
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/A(K,K)
            CALL DSCAL(N-K,T,A(K+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY(N-K,T,A(K+1,K),1,A(K+1,J),1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
*DECK DGESL
      SUBROUTINE DGESL (A, LDA, N, IPVT, B, JOB)
C***BEGIN PROLOGUE  DGESL
C***PURPOSE  Solve the real system A*X=B or TRANS(A)*X=B using the
C            factors computed by DGECO or DGEFA.
C***CATEGORY  D2A1
C***TYPE      DOUBLE PRECISION (SGESL-S, DGESL-D, CGESL-C)
C***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGESL solves the double precision system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGECO or DGEFA.
C
C     On Entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DGECO or DGEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGECO or DGEFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B  where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGECO has set RCOND .GT. 0.0
C        or DGEFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGESL(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DDOT
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGESL
      INTEGER LDA,N,IPVT(*),JOB
      DOUBLE PRECISION A(LDA,*),B(*)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,NM1
C***FIRST EXECUTABLE STATEMENT  DGESL
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL DAXPY(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = DDOT(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + DDOT(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
*DECK DGBFA
      SUBROUTINE DGBFA (ABD, LDA, N, ML, MU, IPVT, INFO)
C***BEGIN PROLOGUE  DGBFA
C***PURPOSE  Factor a band matrix using Gaussian elimination.
C***CATEGORY  D2A2
C***TYPE      DOUBLE PRECISION (SGBFA-S, DGBFA-D, CGBFA-C)
C***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGBFA factors a double precision band matrix by elimination.
C
C     DGBFA is usually called by DGBCO, but it can be called
C     directly with a saving in time if  RCOND  is not needed.
C
C     On Entry
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                contains the matrix in band storage.  The columns
C                of the matrix are stored in the columns of  ABD  and
C                the diagonals of the matrix are stored in rows
C                ML+1 through 2*ML+MU+1 of  ABD .
C                See the comments below for details.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C                LDA must be .GE. 2*ML + MU + 1 .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C                0 .LE. ML .LT.  N .
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C                0 .LE. MU .LT.  N .
C                More efficient if  ML .LE. MU .
C     On Return
C
C        ABD     an upper triangular matrix in band storage and
C                the multipliers which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  if  U(K,K) .EQ. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DGBSL will divide by zero if
C                     called.  Use  RCOND  in DGBCO for a reliable
C                     indication of singularity.
C
C     Band Storage
C
C           If  A  is a band matrix, the following program segment
C           will set up the input.
C
C                   ML = (band width below the diagonal)
C                   MU = (band width above the diagonal)
C                   M = ML + MU + 1
C                   DO 20 J = 1, N
C                      I1 = MAX(1, J-MU)
C                      I2 = MIN(N, J+ML)
C                      DO 10 I = I1, I2
C                         K = I - J + M
C                         ABD(K,J) = A(I,J)
C                10    CONTINUE
C                20 CONTINUE
C
C           This uses rows  ML+1  through  2*ML+MU+1  of  ABD .
C           In addition, the first  ML  rows in  ABD  are used for
C           elements generated during the triangularization.
C           The total number of rows needed in  ABD  is  2*ML+MU+1 .
C           The  ML+MU by ML+MU  upper left triangle and the
C           ML by ML  lower right triangle are not referenced.
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DSCAL, IDAMAX
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGBFA
      INTEGER LDA,N,ML,MU,IPVT(*),INFO
      DOUBLE PRECISION ABD(LDA,*)
C
      DOUBLE PRECISION T
      INTEGER I,IDAMAX,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1
C
C***FIRST EXECUTABLE STATEMENT  DGBFA
      M = ML + MU + 1
      INFO = 0
C
C     ZERO INITIAL FILL-IN COLUMNS
C
      J0 = MU + 2
      J1 = MIN(N,M) - 1
      IF (J1 .LT. J0) GO TO 30
      DO 20 JZ = J0, J1
         I0 = M + 1 - JZ
         DO 10 I = I0, ML
            ABD(I,JZ) = 0.0D0
   10    CONTINUE
   20 CONTINUE
   30 CONTINUE
      JZ = J1
      JU = 0
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 130
      DO 120 K = 1, NM1
         KP1 = K + 1
C
C        ZERO NEXT FILL-IN COLUMN
C
         JZ = JZ + 1
         IF (JZ .GT. N) GO TO 50
         IF (ML .LT. 1) GO TO 50
            DO 40 I = 1, ML
               ABD(I,JZ) = 0.0D0
   40       CONTINUE
   50    CONTINUE
C
C        FIND L = PIVOT INDEX
C
         LM = MIN(ML,N-K)
         L = IDAMAX(LM+1,ABD(M,K),1) + M - 1
         IPVT(K) = L + K - M
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF (ABD(L,K) .EQ. 0.0D0) GO TO 100
C
C           INTERCHANGE IF NECESSARY
C
            IF (L .EQ. M) GO TO 60
               T = ABD(L,K)
               ABD(L,K) = ABD(M,K)
               ABD(M,K) = T
   60       CONTINUE
C
C           COMPUTE MULTIPLIERS
C
            T = -1.0D0/ABD(M,K)
            CALL DSCAL(LM,T,ABD(M+1,K),1)
C
C           ROW ELIMINATION WITH COLUMN INDEXING
C
            JU = MIN(MAX(JU,MU+IPVT(K)),N)
            MM = M
            IF (JU .LT. KP1) GO TO 90
            DO 80 J = KP1, JU
               L = L - 1
               MM = MM - 1
               T = ABD(L,J)
               IF (L .EQ. MM) GO TO 70
                  ABD(L,J) = ABD(MM,J)
                  ABD(MM,J) = T
   70          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
   80       CONTINUE
   90       CONTINUE
         GO TO 110
  100    CONTINUE
            INFO = K
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
      IPVT(N) = N
      IF (ABD(M,N) .EQ. 0.0D0) INFO = N
      RETURN
      END
*DECK DGBSL
      SUBROUTINE DGBSL (ABD, LDA, N, ML, MU, IPVT, B, JOB)
C***BEGIN PROLOGUE  DGBSL
C***PURPOSE  Solve the real band system A*X=B or TRANS(A)*X=B using
C            the factors computed by DGBCO or DGBFA.
C***CATEGORY  D2A2
C***TYPE      DOUBLE PRECISION (SGBSL-S, DGBSL-D, CGBSL-C)
C***KEYWORDS  BANDED, LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE
C***AUTHOR  Moler, C. B., (U. of New Mexico)
C***DESCRIPTION
C
C     DGBSL solves the double precision band system
C     A * X = B  or  TRANS(A) * X = B
C     using the factors computed by DGBCO or DGBFA.
C
C     On Entry
C
C        ABD     DOUBLE PRECISION(LDA, N)
C                the output from DGBCO or DGBFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  ABD .
C
C        N       INTEGER
C                the order of the original matrix.
C
C        ML      INTEGER
C                number of diagonals below the main diagonal.
C
C        MU      INTEGER
C                number of diagonals above the main diagonal.
C
C        IPVT    INTEGER(N)
C                the pivot vector from DGBCO or DGBFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C        JOB     INTEGER
C                = 0         to solve  A*X = B ,
C                = nonzero   to solve  TRANS(A)*X = B , where
C                            TRANS(A)  is the transpose.
C
C     On Return
C
C        B       the solution vector  X .
C
C     Error Condition
C
C        A division by zero will occur if the input factor contains a
C        zero on the diagonal.  Technically this indicates singularity
C        but it is often caused by improper arguments or improper
C        setting of LDA .  It will not occur if the subroutines are
C        called correctly and if DGBCO has set RCOND .GT. 0.0
C        or DGBFA has set INFO .EQ. 0 .
C
C     To compute  INVERSE(A) * C  where  C  is a matrix
C     with  P  columns
C           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF (RCOND is too small) GO TO ...
C           DO 10 J = 1, P
C              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
C        10 CONTINUE
C
C***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
C                 Stewart, LINPACK Users' Guide, SIAM, 1979.
C***ROUTINES CALLED  DAXPY, DDOT
C***REVISION HISTORY  (YYMMDD)
C   780814  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DGBSL
      INTEGER LDA,N,ML,MU,IPVT(*),JOB
      DOUBLE PRECISION ABD(LDA,*),B(*)
C
      DOUBLE PRECISION DDOT,T
      INTEGER K,KB,L,LA,LB,LM,M,NM1
C***FIRST EXECUTABLE STATEMENT  DGBSL
      M = MU + ML + 1
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE L*Y = B
C
         IF (ML .EQ. 0) GO TO 30
         IF (NM1 .LT. 1) GO TO 30
            DO 20 K = 1, NM1
               LM = MIN(ML,N-K)
               L = IPVT(K)
               T = B(L)
               IF (L .EQ. K) GO TO 10
                  B(L) = B(K)
                  B(K) = T
   10          CONTINUE
               CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
   20       CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/ABD(M,K)
            LM = MIN(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = -B(K)
            CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            LM = MIN(K,M) - 1
            LA = M - LM
            LB = K - LM
            T = DDOT(LM,ABD(LA,K),1,B(LB),1)
            B(K) = (B(K) - T)/ABD(M,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (ML .EQ. 0) GO TO 90
         IF (NM1 .LT. 1) GO TO 90
            DO 80 KB = 1, NM1
               K = N - KB
               LM = MIN(ML,N-K)
               B(K) = B(K) + DDOT(LM,ABD(M+1,K),1,B(K+1),1)
               L = IPVT(K)
               IF (L .EQ. K) GO TO 70
                  T = B(L)
                  B(L) = B(K)
                  B(K) = T
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
*DECK DAXPY
      SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DAXPY
C***PURPOSE  Compute a constant times a vector plus a vector.
C***CATEGORY  D1A7
C***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scalar multiplier
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  double precision result (unchanged if N .LE. 0)
C
C     Overwrite double precision DY with double precision DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DAXPY
      DOUBLE PRECISION DX(*), DY(*), DA
C***FIRST EXECUTABLE STATEMENT  DAXPY
      IF (N.LE.0 .OR. DA.EQ.0.0D0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 4.
C
   20 M = MOD(N,4)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N .LT. 4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I+1) = DY(I+1) + DA*DX(I+1)
        DY(I+2) = DY(I+2) + DA*DX(I+2)
        DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DA*DX(I) + DY(I)
   70 CONTINUE
      RETURN
      END
*DECK DCOPY
      SUBROUTINE DCOPY (N, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DCOPY
C***PURPOSE  Copy a vector.
C***CATEGORY  D1A5
C***TYPE      DOUBLE PRECISION (SCOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
C***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  copy of vector DX (unchanged if N .LE. 0)
C
C     Copy double precision DX to double precision DY.
C     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DCOPY
      DOUBLE PRECISION DX(*), DY(*)
C***FIRST EXECUTABLE STATEMENT  DCOPY
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 7.
C
   20 M = MOD(N,7)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF (N .LT. 7) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I+1) = DX(I+1)
        DY(I+2) = DX(I+2)
        DY(I+3) = DX(I+3)
        DY(I+4) = DX(I+4)
        DY(I+5) = DX(I+5)
        DY(I+6) = DX(I+6)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DY(I) = DX(I)
   70 CONTINUE
      RETURN
      END
*DECK DDOT
      DOUBLE PRECISION FUNCTION DDOT (N, DX, INCX, DY, INCY)
C***BEGIN PROLOGUE  DDOT
C***PURPOSE  Compute the inner product of two vectors.
C***CATEGORY  D1A4
C***TYPE      DOUBLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
C***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C     DDOT  double precision dot product (zero if N .LE. 0)
C
C     Returns the dot product of double precision DX and DY.
C     DDOT = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DDOT
      DOUBLE PRECISION DX(*), DY(*)
C***FIRST EXECUTABLE STATEMENT  DDOT
      DDOT = 0.0D0
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. INCY) IF (INCX-1) 5,20,60
C
C     Code for unequal or nonpositive increments.
C
    5 IX = 1
      IY = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C     Code for both increments equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
      DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) +
     1              DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
   50 CONTINUE
      RETURN
C
C     Code for equal, positive, non-unit increments.
C
   60 NS = N*INCX
      DO 70 I = 1,NS,INCX
        DDOT = DDOT + DX(I)*DY(I)
   70 CONTINUE
      RETURN
      END


*DECK DSCAL
      SUBROUTINE DSCAL (N, DA, DX, INCX)
C***BEGIN PROLOGUE  DSCAL
C***PURPOSE  Multiply a vector by a constant.
C***CATEGORY  D1A6
C***TYPE      DOUBLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scale factor
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C       DX  double precision result (unchanged if N.LE.0)
C
C     Replace double precision DX by double precision DA*DX.
C     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DSCAL
      DOUBLE PRECISION DA, DX(*)
      INTEGER I, INCX, IX, M, MP1, N
C***FIRST EXECUTABLE STATEMENT  DSCAL
      IF (N .LE. 0) RETURN
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increment not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DO 10 I = 1,N
        DX(IX) = DA*DX(IX)
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increment equal to 1.
C
C     Clean-up loop so remaining vector length is a multiple of 5.
C
   20 M = MOD(N,5)
      IF (M .EQ. 0) GOTO 40
      DO 30 I = 1,M
        DX(I) = DA*DX(I)
   30 CONTINUE
      IF (N .LT. 5) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DX(I) = DA*DX(I)
        DX(I+1) = DA*DX(I+1)
        DX(I+2) = DA*DX(I+2)
        DX(I+3) = DA*DX(I+3)
        DX(I+4) = DA*DX(I+4)
   50 CONTINUE
      RETURN
      END
*DECK IDAMAX
      INTEGER FUNCTION IDAMAX (N, DX, INCX)
C***BEGIN PROLOGUE  IDAMAX
C***PURPOSE  Find the smallest index of that component of a vector
C            having the maximum magnitude.
C***CATEGORY  D1A2
C***TYPE      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C   IDAMAX  smallest index (zero if N .LE. 0)
C
C     Find smallest index of maximum magnitude of double precision DX.
C     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  IDAMAX
      DOUBLE PRECISION DX(*), DMAX, XMAG
      INTEGER I, INCX, IX, N
C***FIRST EXECUTABLE STATEMENT  IDAMAX
      IDAMAX = 0
      IF (N .LE. 0) RETURN
      IDAMAX = 1
      IF (N .EQ. 1) RETURN
C
      IF (INCX .EQ. 1) GOTO 20
C
C     Code for increments not equal to 1.
C
      IX = 1
      IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
      DMAX = ABS(DX(IX))
      IX = IX + INCX
      DO 10 I = 2,N
        XMAG = ABS(DX(IX))
        IF (XMAG .GT. DMAX) THEN
          IDAMAX = I
          DMAX = XMAG
        ENDIF
        IX = IX + INCX
   10 CONTINUE
      RETURN
C
C     Code for increments equal to 1.
C
   20 DMAX = ABS(DX(1))
      DO 30 I = 2,N
        XMAG = ABS(DX(I))
        IF (XMAG .GT. DMAX) THEN
          IDAMAX = I
          DMAX = XMAG
        ENDIF
   30 CONTINUE
      RETURN
      END

*DECK XSETF
      SUBROUTINE XSETF (MFLAG)
C***BEGIN PROLOGUE  XSETF
C***PURPOSE  Reset the error print control flag.
C***CATEGORY  R3A
C***TYPE      ALL (XSETF-A)
C***KEYWORDS  ERROR CONTROL
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C   XSETF sets the error print control flag to MFLAG:
C      MFLAG=1 means print all messages (the default).
C      MFLAG=0 means no printing.
C
C***SEE ALSO  XERRWD, XERRWV
C***REFERENCES  (NONE)
C***ROUTINES CALLED  IXSAV
C***REVISION HISTORY  (YYMMDD)
C   921118  DATE WRITTEN
C   930329  Added SLATEC format prologue. (FNF)
C   930407  Corrected SEE ALSO section. (FNF)
C   930922  Made user-callable, and other cosmetic changes. (FNF)
C***END PROLOGUE  XSETF
C
C Subroutines called by XSETF.. None
C Function routine called by XSETF.. IXSAV
C-----------------------------------------------------------------------
C**End
      INTEGER MFLAG, JUNK, IXSAV
C
C***FIRST EXECUTABLE STATEMENT  XSETF
      IF (MFLAG .EQ. 0 .OR. MFLAG .EQ. 1) JUNK = IXSAV (2,MFLAG,.TRUE.)
      RETURN
C----------------------- End of Subroutine XSETF -----------------------
      END
*DECK XSETUN
      SUBROUTINE XSETUN (LUN)
C***BEGIN PROLOGUE  XSETUN
C***PURPOSE  Reset the logical unit number for error messages.
C***CATEGORY  R3B
C***TYPE      ALL (XSETUN-A)
C***KEYWORDS  ERROR CONTROL
C***DESCRIPTION
C
C   XSETUN sets the logical unit number for error messages to LUN.
C
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***SEE ALSO  XERRWD, XERRWV
C***REFERENCES  (NONE)
C***ROUTINES CALLED  IXSAV
C***REVISION HISTORY  (YYMMDD)
C   921118  DATE WRITTEN
C   930329  Added SLATEC format prologue. (FNF)
C   930407  Corrected SEE ALSO section. (FNF)
C   930922  Made user-callable, and other cosmetic changes. (FNF)
C***END PROLOGUE  XSETUN
C
C Subroutines called by XSETUN.. None
C Function routine called by XSETUN.. IXSAV
C-----------------------------------------------------------------------
C**End
      INTEGER LUN, JUNK, IXSAV
C
C***FIRST EXECUTABLE STATEMENT  XSETUN
      IF (LUN .GT. 0) JUNK = IXSAV (1,LUN,.TRUE.)
      RETURN
C----------------------- End of Subroutine XSETUN ----------------------
      END
*DECK IXSAV
      INTEGER FUNCTION IXSAV (IPAR, IVALUE, ISET)
C***BEGIN PROLOGUE  IXSAV
C***SUBSIDIARY
C***PURPOSE  Save and recall error message control parameters.
C***CATEGORY  R3C
C***TYPE      ALL (IXSAV-A)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  IXSAV saves and recalls one of two error message parameters:
C    LUNIT, the logical unit number to which messages are printed, and
C    MESFLG, the message print flag.
C  This is a modification of the SLATEC library routine J4SAVE.
C
C  Saved local variables..
C   LUNIT  = Logical unit number for messages.  The default is obtained
C            by a call to IUMACH (may be machine-dependent).
C   MESFLG = Print control flag..
C            1 means print all messages (the default).
C            0 means no printing.
C
C  On input..
C    IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
C    IVALUE = The value to be set for the parameter, if ISET = .TRUE.
C    ISET   = Logical flag to indicate whether to read or write.
C             If ISET = .TRUE., the parameter will be given
C             the value IVALUE.  If ISET = .FALSE., the parameter
C             will be unchanged, and IVALUE is a dummy argument.
C
C  On return..
C    IXSAV = The (old) value of the parameter.
C
C***SEE ALSO  XERRWD, XERRWV
C***ROUTINES CALLED  IUMACH
C***REVISION HISTORY  (YYMMDD)
C   921118  DATE WRITTEN
C   930329  Modified prologue to SLATEC format. (FNF)
C   930915  Added IUMACH call to get default output unit.  (ACH)
C   930922  Minor cosmetic changes. (FNF)
C   010425  Type declaration for IUMACH added. (ACH)
C***END PROLOGUE  IXSAV
C
C Subroutines called by IXSAV.. None
C Function routine called by IXSAV.. IUMACH
C-----------------------------------------------------------------------
C**End
      LOGICAL ISET
      INTEGER IPAR, IVALUE
C-----------------------------------------------------------------------
      INTEGER IUMACH, LUNIT, MESFLG
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this routine.
C-----------------------------------------------------------------------
      SAVE LUNIT, MESFLG
      DATA LUNIT/-1/, MESFLG/1/
C
C***FIRST EXECUTABLE STATEMENT  IXSAV
      IF (IPAR .EQ. 1) THEN
        IF (LUNIT .EQ. -1) LUNIT = IUMACH()
        IXSAV = LUNIT
        IF (ISET) LUNIT = IVALUE
        ENDIF
C
      IF (IPAR .EQ. 2) THEN
        IXSAV = MESFLG
        IF (ISET) MESFLG = IVALUE
        ENDIF
C
      RETURN
C----------------------- End of Function IXSAV -------------------------
      END
*DECK IUMACH
      INTEGER FUNCTION IUMACH()
C***BEGIN PROLOGUE  IUMACH
C***PURPOSE  Provide standard output unit number.
C***CATEGORY  R1
C***TYPE      INTEGER (IUMACH-I)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C *Usage:
C        INTEGER  LOUT, IUMACH
C        LOUT = IUMACH()
C
C *Function Return Values:
C     LOUT : the standard logical unit for Fortran output.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   930915  DATE WRITTEN
C   930922  Made user-callable, and other cosmetic changes. (FNF)
C***END PROLOGUE  IUMACH
C
C*Internal Notes:
C  The built-in value of 6 is standard on a wide range of Fortran
C  systems.  This may be machine-dependent.
C**End
C***FIRST EXECUTABLE STATEMENT  IUMACH
      IUMACH = 6
C
      RETURN
C----------------------- End of Function IUMACH ------------------------
      END


*DECK DSRCOM
      SUBROUTINE DSRCOM (RSAV, ISAV, JOB)
C***BEGIN PROLOGUE  DSRCOM
C***SUBSIDIARY
C***PURPOSE  Save/restore ODEPACK COMMON blocks.
C***TYPE      DOUBLE PRECISION (SSRCOM-S, DSRCOM-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This routine saves or restores (depending on JOB) the contents of
C  the COMMON block DLS001, which is used internally
C  by one or more ODEPACK solvers.
C
C  RSAV = real array of length 218 or more.
C  ISAV = integer array of length 37 or more.
C  JOB  = flag indicating to save or restore the COMMON blocks:
C         JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV)
C         JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV)
C         A call with JOB = 2 presumes a prior call with JOB = 1.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  (NONE)
C***COMMON BLOCKS    DLS001
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   921116  Deleted treatment of block /EH0001/.  (ACH)
C   930801  Reduced Common block length by 2.  (ACH)
C   930809  Renamed to allow single/double precision versions. (ACH)
C   010418  Reduced Common block length by 209+12. (ACH)
C   031105  Restored 'own' variables to Common block /DLS001/, to
C           enable interrupt/restart feature. (ACH)
C   031112  Added SAVE statement for data-loaded constants.
C***END PROLOGUE  DSRCOM
C**End
      INTEGER ISAV, JOB
      INTEGER ILS
      INTEGER I, LENILS, LENRLS
      DOUBLE PRECISION RSAV,   RLS
      DIMENSION RSAV(*), ISAV(*)
      SAVE LENRLS, LENILS
      COMMON /DLS001/ RLS(218), ILS(37)
      DATA LENRLS/218/, LENILS/37/
C
C***FIRST EXECUTABLE STATEMENT  DSRCOM
      IF (JOB .EQ. 2) GO TO 100
C
      DO 10 I = 1,LENRLS
 10     RSAV(I) = RLS(I)
      DO 20 I = 1,LENILS
 20     ISAV(I) = ILS(I)
      RETURN
C
 100  CONTINUE
      DO 110 I = 1,LENRLS
 110     RLS(I) = RSAV(I)
      DO 120 I = 1,LENILS
 120     ILS(I) = ISAV(I)
      RETURN
C----------------------- END OF SUBROUTINE DSRCOM ----------------------
      END


*DECK DVNORM
      DOUBLE PRECISION FUNCTION DVNORM (N, V, W)
C***BEGIN PROLOGUE  DVNORM
C***SUBSIDIARY
C***PURPOSE  Weighted root-mean-square vector norm.
C***TYPE      DOUBLE PRECISION (SVNORM-S, DVNORM-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This function routine computes the weighted root-mean-square norm
C  of the vector of length N contained in the array V, with weights
C  contained in the array W of length N:
C    DVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 )
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C***END PROLOGUE  DVNORM
C**End
      INTEGER N,   I
      DOUBLE PRECISION V, W,   SUM
      DIMENSION V(N), W(N)
C
C***FIRST EXECUTABLE STATEMENT  DVNORM
      SUM = 0.0D0
      DO 10 I = 1,N
 10     SUM = SUM + (V(I)*W(I))**2
      DVNORM = SQRT(SUM/N)
      RETURN
C----------------------- END OF FUNCTION DVNORM ------------------------
      END


*DECK DSRCMS
      SUBROUTINE DSRCMS (RSAV, ISAV, JOB)
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLSS01, which are used
C internally by one or more ODEPACK solvers.
C
C RSAV = real array of length 224 or more.
C ISAV = integer array of length 71 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
      INTEGER ISAV, JOB
      INTEGER ILS, ILSS
      INTEGER I, LENILS, LENISS, LENRLS, LENRSS
      DOUBLE PRECISION RSAV,   RLS, RLSS
      DIMENSION RSAV(*), ISAV(*)
      SAVE LENRLS, LENILS, LENRSS, LENISS
      COMMON /DLS001/ RLS(218), ILS(37)
      COMMON /DLSS01/ RLSS(6), ILSS(34)
      DATA LENRLS/218/, LENILS/37/, LENRSS/6/, LENISS/34/
C
      IF (JOB .EQ. 2) GO TO 100
      DO 10 I = 1,LENRLS
 10     RSAV(I) = RLS(I)
      DO 15 I = 1,LENRSS
 15     RSAV(LENRLS+I) = RLSS(I)
C
      DO 20 I = 1,LENILS
 20     ISAV(I) = ILS(I)
      DO 25 I = 1,LENISS
 25     ISAV(LENILS+I) = ILSS(I)
C
      RETURN
C
 100  CONTINUE
      DO 110 I = 1,LENRLS
 110     RLS(I) = RSAV(I)
      DO 115 I = 1,LENRSS
 115     RLSS(I) = RSAV(LENRLS+I)
C
      DO 120 I = 1,LENILS
 120     ILS(I) = ISAV(I)
      DO 125 I = 1,LENISS
 125     ILSS(I) = ISAV(LENILS+I)
C
      RETURN
C----------------------- End of Subroutine DSRCMS ----------------------
      END
*DECK ODRV
      subroutine odrv
     *     (n, ia,ja,a, p,ip, nsp,isp, path, flag)
c                                                                 5/2/83
c***********************************************************************
c  odrv -- driver for sparse matrix reordering routines
c***********************************************************************
c
c  description
c
c    odrv finds a minimum degree ordering of the rows and columns
c    of a matrix m stored in (ia,ja,a) format (see below).  for the
c    reordered matrix, the work and storage required to perform
c    gaussian elimination is (usually) significantly less.
c
c    note.. odrv and its subordinate routines have been modified to
c    compute orderings for general matrices, not necessarily having any
c    symmetry.  the miminum degree ordering is computed for the
c    structure of the symmetric matrix  m + m-transpose.
c    modifications to the original odrv module have been made in
c    the coding in subroutine mdi, and in the initial comments in
c    subroutines odrv and md.
c
c    if only the nonzero entries in the upper triangle of m are being
c    stored, then odrv symmetrically reorders (ia,ja,a), (optionally)
c    with the diagonal entries placed first in each row.  this is to
c    ensure that if m(i,j) will be in the upper triangle of m with
c    respect to the new ordering, then m(i,j) is stored in row i (and
c    thus m(j,i) is not stored),  whereas if m(i,j) will be in the
c    strict lower triangle of m, then m(j,i) is stored in row j (and
c    thus m(i,j) is not stored).
c
c
c  storage of sparse matrices
c
c    the nonzero entries of the matrix m are stored row-by-row in the
c    array a.  to identify the individual nonzero entries in each row,
c    we need to know in which column each entry lies.  these column
c    indices are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
c    ja(k) = j.  to identify the individual rows, we need to know where
c    each row starts.  these row pointers are stored in the array ia.
c    i.e., if m(i,j) is the first nonzero entry (stored) in the i-th row
c    and  a(k) = m(i,j),  then  ia(i) = k.  moreover, ia(n+1) points to
c    the first location following the last element in the last row.
c    thus, the number of entries in the i-th row is  ia(i+1) - ia(i),
c    the nonzero entries in the i-th row are stored consecutively in
c
c            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
c
c    and the corresponding column indices are stored consecutively in
c
c            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
c
c    when the coefficient matrix is symmetric, only the nonzero entries
c    in the upper triangle need be stored.  for example, the matrix
c
c             ( 1  0  2  3  0 )
c             ( 0  4  0  0  0 )
c         m = ( 2  0  5  6  0 )
c             ( 3  0  6  7  8 )
c             ( 0  0  0  8  9 )
c
c    could be stored as
c
c            - 1  2  3  4  5  6  7  8  9 10 11 12 13
c         ---+--------------------------------------
c         ia - 1  4  5  8 12 14
c         ja - 1  3  4  2  1  3  4  1  3  4  5  4  5
c          a - 1  2  3  4  2  5  6  3  6  7  8  8  9
c
c    or (symmetrically) as
c
c            - 1  2  3  4  5  6  7  8  9
c         ---+--------------------------
c         ia - 1  4  5  7  9 10
c         ja - 1  3  4  2  3  4  4  5  5
c          a - 1  2  3  4  5  6  7  8  9          .
c
c
c  parameters
c
c    n    - order of the matrix
c
c    ia   - integer one-dimensional array containing pointers to delimit
c           rows in ja and a.  dimension = n+1
c
c    ja   - integer one-dimensional array containing the column indices
c           corresponding to the elements of a.  dimension = number of
c           nonzero entries in (the upper triangle of) m
c
c    a    - real one-dimensional array containing the nonzero entries in
c           (the upper triangle of) m, stored by rows.  dimension =
c           number of nonzero entries in (the upper triangle of) m
c
c    p    - integer one-dimensional array used to return the permutation
c           of the rows and columns of m corresponding to the minimum
c           degree ordering.  dimension = n
c
c    ip   - integer one-dimensional array used to return the inverse of
c           the permutation returned in p.  dimension = n
c
c    nsp  - declared dimension of the one-dimensional array isp.  nsp
c           must be at least  3n+4k,  where k is the number of nonzeroes
c           in the strict upper triangle of m
c
c    isp  - integer one-dimensional array used for working storage.
c           dimension = nsp
c
c    path - integer path specification.  values and their meanings are -
c             1  find minimum degree ordering only
c             2  find minimum degree ordering and reorder symmetrically
c                  stored matrix (used when only the nonzero entries in
c                  the upper triangle of m are being stored)
c             3  reorder symmetrically stored matrix as specified by
c                  input permutation (used when an ordering has already
c                  been determined and only the nonzero entries in the
c                  upper triangle of m are being stored)
c             4  same as 2 but put diagonal entries at start of each row
c             5  same as 3 but put diagonal entries at start of each row
c
c    flag - integer error flag.  values and their meanings are -
c               0    no errors detected
c              9n+k  insufficient storage in md
c             10n+1  insufficient storage in odrv
c             11n+1  illegal path specification
c
c
c  conversion from real to double precision
c
c    change the real declarations in odrv and sro to double precision
c    declarations.
c
c-----------------------------------------------------------------------
c
      integer  ia(*), ja(*),  p(*), ip(*),  isp(*),  path,  flag,
     *   v, l, head,  tmp, q
c...  real  a(*)
      double precision  a(*)
      logical  dflag
c
c----initialize error flag and validate path specification
      flag = 0
      if (path.lt.1 .or. 5.lt.path)  go to 111
c
c----allocate storage and find minimum degree ordering
      if ((path-1) * (path-2) * (path-4) .ne. 0)  go to 1
        max = (nsp-n)/2
        v    = 1
        l    = v     +  max
        head = l     +  max
        next = head  +  n
        if (max.lt.n)  go to 110
c
        call  md
     *     (n, ia,ja, max,isp(v),isp(l), isp(head),p,ip, isp(v), flag)
        if (flag.ne.0)  go to 100
c
c----allocate storage and symmetrically reorder matrix
   1  if ((path-2) * (path-3) * (path-4) * (path-5) .ne. 0)  go to 2
        tmp = (nsp+1) -      n
        q   = tmp     - (ia(n+1)-1)
        if (q.lt.1)  go to 110
c
        dflag = path.eq.4 .or. path.eq.5
        call sro
     *     (n,  ip,  ia, ja, a,  isp(tmp),  isp(q),  dflag)
c
   2  return
c
c ** error -- error detected in md
 100  return
c ** error -- insufficient storage
 110  flag = 10*n + 1
      return
c ** error -- illegal path specified
 111  flag = 11*n + 1
      return
      end
      subroutine md
     *     (n, ia,ja, max, v,l, head,last,next, mark, flag)
c***********************************************************************
c  md -- minimum degree algorithm (based on element model)
c***********************************************************************
c
c  description
c
c    md finds a minimum degree ordering of the rows and columns of a
c    general sparse matrix m stored in (ia,ja,a) format.
c    when the structure of m is nonsymmetric, the ordering is that
c    obtained for the symmetric matrix  m + m-transpose.
c
c
c  additional parameters
c
c    max  - declared dimension of the one-dimensional arrays v and l.
c           max must be at least  n+2k,  where k is the number of
c           nonzeroes in the strict upper triangle of m + m-transpose
c
c    v    - integer one-dimensional work array.  dimension = max
c
c    l    - integer one-dimensional work array.  dimension = max
c
c    head - integer one-dimensional work array.  dimension = n
c
c    last - integer one-dimensional array used to return the permutation
c           of the rows and columns of m corresponding to the minimum
c           degree ordering.  dimension = n
c
c    next - integer one-dimensional array used to return the inverse of
c           the permutation returned in last.  dimension = n
c
c    mark - integer one-dimensional work array (may be the same as v).
c           dimension = n
c
c    flag - integer error flag.  values and their meanings are -
c             0     no errors detected
c             9n+k  insufficient storage in md
c
c
c  definitions of internal parameters
c
c    ---------+---------------------------------------------------------
c    v(s)     - value field of list entry
c    ---------+---------------------------------------------------------
c    l(s)     - link field of list entry  (0 =) end of list)
c    ---------+---------------------------------------------------------
c    l(vi)    - pointer to element list of uneliminated vertex vi
c    ---------+---------------------------------------------------------
c    l(ej)    - pointer to boundary list of active element ej
c    ---------+---------------------------------------------------------
c    head(d)  - vj =) vj head of d-list d
c             -  0 =) no vertex in d-list d
c
c
c             -                  vi uneliminated vertex
c             -          vi in ek           -       vi not in ek
c    ---------+-----------------------------+---------------------------
c    next(vi) - undefined but nonnegative   - vj =) vj next in d-list
c             -                             -  0 =) vi tail of d-list
c    ---------+-----------------------------+---------------------------
c    last(vi) - (not set until mdp)         - -d =) vi head of d-list d
c             --vk =) compute degree        - vj =) vj last in d-list
c             - ej =) vi prototype of ej    -  0 =) vi not in any d-list
c             -  0 =) do not compute degree -
c    ---------+-----------------------------+---------------------------
c    mark(vi) - mark(vk)                    - nonneg. tag .lt. mark(vk)
c
c
c             -                   vi eliminated vertex
c             -      ei active element      -           otherwise
c    ---------+-----------------------------+---------------------------
c    next(vi) - -j =) vi was j-th vertex    - -j =) vi was j-th vertex
c             -       to be eliminated      -       to be eliminated
c    ---------+-----------------------------+---------------------------
c    last(vi) -  m =) size of ei = m        - undefined
c    ---------+-----------------------------+---------------------------
c    mark(vi) - -m =) overlap count of ei   - undefined
c             -       with ek = m           -
c             - otherwise nonnegative tag   -
c             -       .lt. mark(vk)         -
c
c-----------------------------------------------------------------------
c
      integer  ia(*), ja(*),  v(*), l(*),  head(*), last(*), next(*),
     *   mark(*),  flag,  tag, dmin, vk,ek, tail
      equivalence  (vk,ek)
c
c----initialization
      tag = 0
      call  mdi
     *   (n, ia,ja, max,v,l, head,last,next, mark,tag, flag)
      if (flag.ne.0)  return
c
      k = 0
      dmin = 1
c
c----while  k .lt. n  do
   1  if (k.ge.n)  go to 4
c
c------search for vertex of minimum degree
   2    if (head(dmin).gt.0)  go to 3
          dmin = dmin + 1
          go to 2
c
c------remove vertex vk of minimum degree from degree list
   3    vk = head(dmin)
        head(dmin) = next(vk)
        if (head(dmin).gt.0)  last(head(dmin)) = -dmin
c
c------number vertex vk, adjust tag, and tag vk
        k = k+1
        next(vk) = -k
        last(ek) = dmin - 1
        tag = tag + last(ek)
        mark(vk) = tag
c
c------form element ek from uneliminated neighbors of vk
        call  mdm
     *     (vk,tail, v,l, last,next, mark)
c
c------purge inactive elements and do mass elimination
        call  mdp
     *     (k,ek,tail, v,l, head,last,next, mark)
c
c------update degrees of uneliminated vertices in ek
        call  mdu
     *     (ek,dmin, v,l, head,last,next, mark)
c
        go to 1
c
c----generate inverse permutation from permutation
   4  do 5 k=1,n
        next(k) = -next(k)
   5    last(next(k)) = k
c
      return
      end
      subroutine mdi
     *     (n, ia,ja, max,v,l, head,last,next, mark,tag, flag)
c***********************************************************************
c  mdi -- initialization
c***********************************************************************
      integer  ia(*), ja(*),  v(*), l(*),  head(*), last(*), next(*),
     *   mark(*), tag,  flag,  sfs, vi,dvi, vj
c
c----initialize degrees, element lists, and degree lists
      do 1 vi=1,n
        mark(vi) = 1
        l(vi) = 0
   1    head(vi) = 0
      sfs = n+1
c
c----create nonzero structure
c----for each nonzero entry a(vi,vj)
      do 6 vi=1,n
        jmin = ia(vi)
        jmax = ia(vi+1) - 1
        if (jmin.gt.jmax)  go to 6
        do 5 j=jmin,jmax
          vj = ja(j)
          if (vj-vi) 2, 5, 4
c
c------if a(vi,vj) is in strict lower triangle
c------check for previous occurrence of a(vj,vi)
   2      lvk = vi
          kmax = mark(vi) - 1
          if (kmax .eq. 0) go to 4
          do 3 k=1,kmax
            lvk = l(lvk)
            if (v(lvk).eq.vj) go to 5
   3        continue
c----for unentered entries a(vi,vj)
   4        if (sfs.ge.max)  go to 101
c
c------enter vj in element list for vi
            mark(vi) = mark(vi) + 1
            v(sfs) = vj
            l(sfs) = l(vi)
            l(vi) = sfs
            sfs = sfs+1
c
c------enter vi in element list for vj
            mark(vj) = mark(vj) + 1
            v(sfs) = vi
            l(sfs) = l(vj)
            l(vj) = sfs
            sfs = sfs+1
   5      continue
   6    continue
c
c----create degree lists and initialize mark vector
      do 7 vi=1,n
        dvi = mark(vi)
        next(vi) = head(dvi)
        head(dvi) = vi
        last(vi) = -dvi
        nextvi = next(vi)
        if (nextvi.gt.0)  last(nextvi) = vi
   7    mark(vi) = tag
c
      return
c
c ** error-  insufficient storage
 101  flag = 9*n + vi
      return
      end
      subroutine mdm
     *     (vk,tail, v,l, last,next, mark)
c***********************************************************************
c  mdm -- form element from uneliminated neighbors of vk
c***********************************************************************
      integer  vk, tail,  v(*), l(*),   last(*), next(*),   mark(*),
     *   tag, s,ls,vs,es, b,lb,vb, blp,blpmax
      equivalence  (vs, es)
c
c----initialize tag and list of uneliminated neighbors
      tag = mark(vk)
      tail = vk
c
c----for each vertex/element vs/es in element list of vk
      ls = l(vk)
   1  s = ls
      if (s.eq.0)  go to 5
        ls = l(s)
        vs = v(s)
        if (next(vs).lt.0)  go to 2
c
c------if vs is uneliminated vertex, then tag and append to list of
c------uneliminated neighbors
          mark(vs) = tag
          l(tail) = s
          tail = s
          go to 4
c
c------if es is active element, then ...
c--------for each vertex vb in boundary list of element es
   2      lb = l(es)
          blpmax = last(es)
          do 3 blp=1,blpmax
            b = lb
            lb = l(b)
            vb = v(b)
c
c----------if vb is untagged vertex, then tag and append to list of
c----------uneliminated neighbors
            if (mark(vb).ge.tag)  go to 3
              mark(vb) = tag
              l(tail) = b
              tail = b
   3        continue
c
c--------mark es inactive
          mark(es) = tag
c
   4    go to 1
c
c----terminate list of uneliminated neighbors
   5  l(tail) = 0
c
      return
      end
      subroutine mdp
     *     (k,ek,tail, v,l, head,last,next, mark)
c***********************************************************************
c  mdp -- purge inactive elements and do mass elimination
c***********************************************************************
      integer  ek, tail,  v(*), l(*),  head(*), last(*), next(*),
     *   mark(*),  tag, free, li,vi,lvi,evi, s,ls,es, ilp,ilpmax
c
c----initialize tag
      tag = mark(ek)
c
c----for each vertex vi in ek
      li = ek
      ilpmax = last(ek)
      if (ilpmax.le.0)  go to 12
      do 11 ilp=1,ilpmax
        i = li
        li = l(i)
        vi = v(li)
c
c------remove vi from degree list
        if (last(vi).eq.0)  go to 3
          if (last(vi).gt.0)  go to 1
            head(-last(vi)) = next(vi)
            go to 2
   1        next(last(vi)) = next(vi)
   2      if (next(vi).gt.0)  last(next(vi)) = last(vi)
c
c------remove inactive items from element list of vi
   3    ls = vi
   4    s = ls
        ls = l(s)
        if (ls.eq.0)  go to 6
          es = v(ls)
          if (mark(es).lt.tag)  go to 5
            free = ls
            l(s) = l(ls)
            ls = s
   5      go to 4
c
c------if vi is interior vertex, then remove from list and eliminate
   6    lvi = l(vi)
        if (lvi.ne.0)  go to 7
          l(i) = l(li)
          li = i
c
          k = k+1
          next(vi) = -k
          last(ek) = last(ek) - 1
          go to 11
c
c------else ...
c--------classify vertex vi
   7      if (l(lvi).ne.0)  go to 9
            evi = v(lvi)
            if (next(evi).ge.0)  go to 9
              if (mark(evi).lt.0)  go to 8
c
c----------if vi is prototype vertex, then mark as such, initialize
c----------overlap count for corresponding element, and move vi to end
c----------of boundary list
                last(vi) = evi
                mark(evi) = -1
                l(tail) = li
                tail = li
                l(i) = l(li)
                li = i
                go to 10
c
c----------else if vi is duplicate vertex, then mark as such and adjust
c----------overlap count for corresponding element
   8            last(vi) = 0
                mark(evi) = mark(evi) - 1
                go to 10
c
c----------else mark vi to compute degree
   9            last(vi) = -ek
c
c--------insert ek in element list of vi
  10      v(free) = ek
          l(free) = l(vi)
          l(vi) = free
  11    continue
c
c----terminate boundary list
  12  l(tail) = 0
c
      return
      end
      subroutine mdu
     *     (ek,dmin, v,l, head,last,next, mark)
c***********************************************************************
c  mdu -- update degrees of uneliminated vertices in ek
c***********************************************************************
      integer  ek, dmin,  v(*), l(*),  head(*), last(*), next(*),
     *   mark(*),  tag, vi,evi,dvi, s,vs,es, b,vb, ilp,ilpmax,
     *   blp,blpmax
      equivalence  (vs, es)
c
c----initialize tag
      tag = mark(ek) - last(ek)
c
c----for each vertex vi in ek
      i = ek
      ilpmax = last(ek)
      if (ilpmax.le.0)  go to 11
      do 10 ilp=1,ilpmax
        i = l(i)
        vi = v(i)
        if (last(vi))  1, 10, 8
c
c------if vi neither prototype nor duplicate vertex, then merge elements
c------to compute degree
   1      tag = tag + 1
          dvi = last(ek)
c
c--------for each vertex/element vs/es in element list of vi
          s = l(vi)
   2      s = l(s)
          if (s.eq.0)  go to 9
            vs = v(s)
            if (next(vs).lt.0)  go to 3
c
c----------if vs is uneliminated vertex, then tag and adjust degree
              mark(vs) = tag
              dvi = dvi + 1
              go to 5
c
c----------if es is active element, then expand
c------------check for outmatched vertex
   3          if (mark(es).lt.0)  go to 6
c
c------------for each vertex vb in es
              b = es
              blpmax = last(es)
              do 4 blp=1,blpmax
                b = l(b)
                vb = v(b)
c
c--------------if vb is untagged, then tag and adjust degree
                if (mark(vb).ge.tag)  go to 4
                  mark(vb) = tag
                  dvi = dvi + 1
   4            continue
c
   5        go to 2
c
c------else if vi is outmatched vertex, then adjust overlaps but do not
c------compute degree
   6      last(vi) = 0
          mark(es) = mark(es) - 1
   7      s = l(s)
          if (s.eq.0)  go to 10
            es = v(s)
            if (mark(es).lt.0)  mark(es) = mark(es) - 1
            go to 7
c
c------else if vi is prototype vertex, then calculate degree by
c------inclusion/exclusion and reset overlap count
   8      evi = last(vi)
          dvi = last(ek) + last(evi) + mark(evi)
          mark(evi) = 0
c
c------insert vi in appropriate degree list
   9    next(vi) = head(dvi)
        head(dvi) = vi
        last(vi) = -dvi
        if (next(vi).gt.0)  last(next(vi)) = vi
        if (dvi.lt.dmin)  dmin = dvi
c
  10    continue
c
  11  return
      end
      subroutine sro
     *     (n, ip, ia,ja,a, q, r, dflag)
c***********************************************************************
c  sro -- symmetric reordering of sparse symmetric matrix
c***********************************************************************
c
c  description
c
c    the nonzero entries of the matrix m are assumed to be stored
c    symmetrically in (ia,ja,a) format (i.e., not both m(i,j) and m(j,i)
c    are stored if i ne j).
c
c    sro does not rearrange the order of the rows, but does move
c    nonzeroes from one row to another to ensure that if m(i,j) will be
c    in the upper triangle of m with respect to the new ordering, then
c    m(i,j) is stored in row i (and thus m(j,i) is not stored),  whereas
c    if m(i,j) will be in the strict lower triangle of m, then m(j,i) is
c    stored in row j (and thus m(i,j) is not stored).
c
c
c  additional parameters
c
c    q     - integer one-dimensional work array.  dimension = n
c
c    r     - integer one-dimensional work array.  dimension = number of
c            nonzero entries in the upper triangle of m
c
c    dflag - logical variable.  if dflag = .true., then store nonzero
c            diagonal elements at the beginning of the row
c
c-----------------------------------------------------------------------
c
      integer  ip(*),  ia(*), ja(*),  q(*), r(*)
c...  real  a(*),  ak
      double precision  a(*),  ak
      logical  dflag
c
c
c--phase 1 -- find row in which to store each nonzero
c----initialize count of nonzeroes to be stored in each row
      do 1 i=1,n
  1     q(i) = 0
c
c----for each nonzero element a(j)
      do 3 i=1,n
        jmin = ia(i)
        jmax = ia(i+1) - 1
        if (jmin.gt.jmax)  go to 3
        do 2 j=jmin,jmax
c
c--------find row (=r(j)) and column (=ja(j)) in which to store a(j) ...
          k = ja(j)
          if (ip(k).lt.ip(i))  ja(j) = i
          if (ip(k).ge.ip(i))  k = i
          r(j) = k
c
c--------... and increment count of nonzeroes (=q(r(j)) in that row
  2       q(k) = q(k) + 1
  3     continue
c
c
c--phase 2 -- find new ia and permutation to apply to (ja,a)
c----determine pointers to delimit rows in permuted (ja,a)
      do 4 i=1,n
        ia(i+1) = ia(i) + q(i)
  4     q(i) = ia(i+1)
c
c----determine where each (ja(j),a(j)) is stored in permuted (ja,a)
c----for each nonzero element (in reverse order)
      ilast = 0
      jmin = ia(1)
      jmax = ia(n+1) - 1
      j = jmax
      do 6 jdummy=jmin,jmax
        i = r(j)
        if (.not.dflag .or. ja(j).ne.i .or. i.eq.ilast)  go to 5
c
c------if dflag, then put diagonal nonzero at beginning of row
          r(j) = ia(i)
          ilast = i
          go to 6
c
c------put (off-diagonal) nonzero in last unused location in row
  5       q(i) = q(i) - 1
          r(j) = q(i)
c
  6     j = j-1
c
c
c--phase 3 -- permute (ja,a) to upper triangular form (wrt new ordering)
      do 8 j=jmin,jmax
  7     if (r(j).eq.j)  go to 8
          k = r(j)
          r(j) = r(k)
          r(k) = k
          jak = ja(k)
          ja(k) = ja(j)
          ja(j) = jak
          ak = a(k)
          a(k) = a(j)
          a(j) = ak
          go to 7
  8     continue
c
      return
      end
*DECK CDRV
      subroutine cdrv
     *     (n, r,c,ic, ia,ja,a, b, z, nsp,isp,rsp,esp, path, flag)
c*** subroutine cdrv
c*** driver for subroutines for solving sparse nonsymmetric systems of
c       linear equations (compressed pointer storage)
c
c
c    parameters
c    class abbreviations are--
c       n - integer variable
c       f - real variable
c       v - supplies a value to the driver
c       r - returns a result from the driver
c       i - used internally by the driver
c       a - array
c
c class - parameter
c ------+----------
c       -
c         the nonzero entries of the coefficient matrix m are stored
c    row-by-row in the array a.  to identify the individual nonzero
c    entries in each row, we need to know in which column each entry
c    lies.  the column indices which correspond to the nonzero entries
c    of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
c    ja(k) = j.  in addition, we need to know where each row starts and
c    how long it is.  the index positions in ja and a where the rows of
c    m begin are stored in the array ia.  i.e., if m(i,j) is the first
c    nonzero entry (stored) in the i-th row and a(k) = m(i,j),  then
c    ia(i) = k.  moreover, the index in ja and a of the first location
c    following the last element in the last row is stored in ia(n+1).
c    thus, the number of entries in the i-th row is given by
c    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
c    consecutively in
c            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
c    and the corresponding column indices are stored consecutively in
c            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
c    for example, the 5 by 5 matrix
c                ( 1. 0. 2. 0. 0.)
c                ( 0. 3. 0. 0. 0.)
c            m = ( 0. 4. 5. 6. 0.)
c                ( 0. 0. 0. 7. 0.)
c                ( 0. 0. 0. 8. 9.)
c    would be stored as
c               - 1  2  3  4  5  6  7  8  9
c            ---+--------------------------
c            ia - 1  3  4  7  8 10
c            ja - 1  3  2  2  3  4  4  4  5
c             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
c
c nv    - n     - number of variables/equations.
c fva   - a     - nonzero entries of the coefficient matrix m, stored
c       -           by rows.
c       -           size = number of nonzero entries in m.
c nva   - ia    - pointers to delimit the rows in a.
c       -           size = n+1.
c nva   - ja    - column numbers corresponding to the elements of a.
c       -           size = size of a.
c fva   - b     - right-hand side b.  b and z can the same array.
c       -           size = n.
c fra   - z     - solution x.  b and z can be the same array.
c       -           size = n.
c
c         the rows and columns of the original matrix m can be
c    reordered (e.g., to reduce fillin or ensure numerical stability)
c    before calling the driver.  if no reordering is done, then set
c    r(i) = c(i) = ic(i) = i  for i=1,...,n.  the solution z is returned
c    in the original order.
c         if the columns have been reordered (i.e.,  c(i).ne.i  for some
c    i), then the driver will call a subroutine (nroc) which rearranges
c    each row of ja and a, leaving the rows in the original order, but
c    placing the elements of each row in increasing order with respect
c    to the new ordering.  if  path.ne.1,  then nroc is assumed to have
c    been called already.
c
c nva   - r     - ordering of the rows of m.
c       -           size = n.
c nva   - c     - ordering of the columns of m.
c       -           size = n.
c nva   - ic    - inverse of the ordering of the columns of m.  i.e.,
c       -           ic(c(i)) = i  for i=1,...,n.
c       -           size = n.
c
c         the solution of the system of linear equations is divided into
c    three stages --
c      nsfc -- the matrix m is processed symbolically to determine where
c               fillin will occur during the numeric factorization.
c      nnfc -- the matrix m is factored numerically into the product ldu
c               of a unit lower triangular matrix l, a diagonal matrix
c               d, and a unit upper triangular matrix u, and the system
c               mx = b  is solved.
c      nnsc -- the linear system  mx = b  is solved using the ldu
c  or           factorization from nnfc.
c      nntc -- the transposed linear system  mt x = b  is solved using
c               the ldu factorization from nnf.
c    for several systems whose coefficient matrices have the same
c    nonzero structure, nsfc need be done only once (for the first
c    system).  then nnfc is done once for each additional system.  for
c    several systems with the same coefficient matrix, nsfc and nnfc
c    need be done only once (for the first system).  then nnsc or nntc
c    is done once for each additional right-hand side.
c
c nv    - path  - path specification.  values and their meanings are --
c       -           1  perform nroc, nsfc, and nnfc.
c       -           2  perform nnfc only  (nsfc is assumed to have been
c       -               done in a manner compatible with the storage
c       -               allocation used in the driver).
c       -           3  perform nnsc only  (nsfc and nnfc are assumed to
c       -               have been done in a manner compatible with the
c       -               storage allocation used in the driver).
c       -           4  perform nntc only  (nsfc and nnfc are assumed to
c       -               have been done in a manner compatible with the
c       -               storage allocation used in the driver).
c       -           5  perform nroc and nsfc.
c
c         various errors are detected by the driver and the individual
c    subroutines.
c
c nr    - flag  - error flag.  values and their meanings are --
c       -             0     no errors detected
c       -             n+k   null row in a  --  row = k
c       -            2n+k   duplicate entry in a  --  row = k
c       -            3n+k   insufficient storage in nsfc  --  row = k
c       -            4n+1   insufficient storage in nnfc
c       -            5n+k   null pivot  --  row = k
c       -            6n+k   insufficient storage in nsfc  --  row = k
c       -            7n+1   insufficient storage in nnfc
c       -            8n+k   zero pivot  --  row = k
c       -           10n+1   insufficient storage in cdrv
c       -           11n+1   illegal path specification
c
c         working storage is needed for the factored form of the matrix
c    m plus various temporary vectors.  the arrays isp and rsp should be
c    equivalenced.  integer storage is allocated from the beginning of
c    isp and real storage from the end of rsp.
c
c nv    - nsp   - declared dimension of rsp.  nsp generally must
c       -           be larger than  8n+2 + 2k  (where  k = (number of
c       -           nonzero entries in m)).
c nvira - isp   - integer working storage divided up into various arrays
c       -           needed by the subroutines.  isp and rsp should be
c       -           equivalenced.
c       -           size = lratio*nsp.
c fvira - rsp   - real working storage divided up into various arrays
c       -           needed by the subroutines.  isp and rsp should be
c       -           equivalenced.
c       -           size = nsp.
c nr    - esp   - if sufficient storage was available to perform the
c       -           symbolic factorization (nsfc), then esp is set to
c       -           the amount of excess storage provided (negative if
c       -           insufficient storage was available to perform the
c       -           numeric factorization (nnfc)).
c
c
c  conversion to double precision
c
c    to convert these routines for double precision arrays..
c    (1) use the double precision declarations in place of the real
c    declarations in each subprogram, as given in comment cards.
c    (2) change the data-loaded value of the integer  lratio
c    in subroutine cdrv, as indicated below.
c    (3) change e0 to d0 in the constants in statement number 10
c    in subroutine nnfc and the line following that.
c
      integer  r(*), c(*), ic(*),  ia(*), ja(*),  isp(*), esp,  path,
     *   flag,  d, u, q, row, tmp, ar,  umax
c     real  a(*), b(*), z(*), rsp(*)
      double precision  a(*), b(*), z(*), rsp(*)
c
c  set lratio equal to the ratio between the length of floating point
c  and integer array data.  e. g., lratio = 1 for (real, integer),
c  lratio = 2 for (double precision, integer)
c
      data lratio/2/
c
      if (path.lt.1 .or. 5.lt.path)  go to 111
c******initialize and divide up temporary storage  *******************
      il   = 1
      ijl  = il  + (n+1)
      iu   = ijl +   n
      iju  = iu  + (n+1)
      irl  = iju +   n
      jrl  = irl +   n
      jl   = jrl +   n
c
c  ******  reorder a if necessary, call nsfc if flag is set  ***********
      if ((path-1) * (path-5) .ne. 0)  go to 5
        max = (lratio*nsp + 1 - jl) - (n+1) - 5*n
        jlmax = max/2
        q     = jl   + jlmax
        ira   = q    + (n+1)
        jra   = ira  +   n
        irac  = jra  +   n
        iru   = irac +   n
        jru   = iru  +   n
        jutmp = jru  +   n
        jumax = lratio*nsp  + 1 - jutmp
        esp = max/lratio
        if (jlmax.le.0 .or. jumax.le.0)  go to 110
c
        do 1 i=1,n
          if (c(i).ne.i)  go to 2
   1      continue
        go to 3
   2    ar = nsp + 1 - n
        call  nroc
     *     (n, ic, ia,ja,a, isp(il), rsp(ar), isp(iu), flag)
        if (flag.ne.0)  go to 100
c
   3    call  nsfc
     *     (n, r, ic, ia,ja,
     *      jlmax, isp(il), isp(jl), isp(ijl),
     *      jumax, isp(iu), isp(jutmp), isp(iju),
     *      isp(q), isp(ira), isp(jra), isp(irac),
     *      isp(irl), isp(jrl), isp(iru), isp(jru),  flag)
        if(flag .ne. 0)  go to 100
c  ******  move ju next to jl  *****************************************
        jlmax = isp(ijl+n-1)
        ju    = jl + jlmax
        jumax = isp(iju+n-1)
        if (jumax.le.0)  go to 5
        do 4 j=1,jumax
   4      isp(ju+j-1) = isp(jutmp+j-1)
c
c  ******  call remaining subroutines  *********************************
   5  jlmax = isp(ijl+n-1)
      ju    = jl  + jlmax
      jumax = isp(iju+n-1)
      l     = (ju + jumax - 2 + lratio)  /  lratio    +    1
      lmax  = isp(il+n) - 1
      d     = l   + lmax
      u     = d   + n
      row   = nsp + 1 - n
      tmp   = row - n
      umax  = tmp - u
      esp   = umax - (isp(iu+n) - 1)
c
      if ((path-1) * (path-2) .ne. 0)  go to 6
        if (umax.lt.0)  go to 110
        call nnfc
     *     (n,  r, c, ic,  ia, ja, a, z, b,
     *      lmax, isp(il), isp(jl), isp(ijl), rsp(l),  rsp(d),
     *      umax, isp(iu), isp(ju), isp(iju), rsp(u),
     *      rsp(row), rsp(tmp),  isp(irl), isp(jrl),  flag)
        if(flag .ne. 0)  go to 100
c
   6  if ((path-3) .ne. 0)  go to 7
        call nnsc
     *     (n,  r, c,  isp(il), isp(jl), isp(ijl), rsp(l),
     *      rsp(d),    isp(iu), isp(ju), isp(iju), rsp(u),
     *      z, b,  rsp(tmp))
c
   7  if ((path-4) .ne. 0)  go to 8
        call nntc
     *     (n,  r, c,  isp(il), isp(jl), isp(ijl), rsp(l),
     *      rsp(d),    isp(iu), isp(ju), isp(iju), rsp(u),
     *      z, b,  rsp(tmp))
   8  return
c
c ** error.. error detected in nroc, nsfc, nnfc, or nnsc
 100  return
c ** error.. insufficient storage
 110  flag = 10*n + 1
      return
c ** error.. illegal path specification
 111  flag = 11*n + 1
      return
      end
      subroutine nroc (n, ic, ia, ja, a, jar, ar, p, flag)
c
c       ----------------------------------------------------------------
c
c               yale sparse matrix package - nonsymmetric codes
c                    solving the system of equations mx = b
c
c    i.   calling sequences
c         the coefficient matrix can be processed by an ordering routine
c    (e.g., to reduce fillin or ensure numerical stability) before using
c    the remaining subroutines.  if no reordering is done, then set
c    r(i) = c(i) = ic(i) = i  for i=1,...,n.  if an ordering subroutine
c    is used, then nroc should be used to reorder the coefficient matrix
c    the calling sequence is --
c        (       (matrix ordering))
c        (nroc   (matrix reordering))
c         nsfc   (symbolic factorization to determine where fillin will
c                  occur during numeric factorization)
c         nnfc   (numeric factorization into product ldu of unit lower
c                  triangular matrix l, diagonal matrix d, and unit
c                  upper triangular matrix u, and solution of linear
c                  system)
c         nnsc   (solution of linear system for additional right-hand
c                  side using ldu factorization from nnfc)
c    (if only one system of equations is to be solved, then the
c    subroutine trk should be used.)
c
c    ii.  storage of sparse matrices
c         the nonzero entries of the coefficient matrix m are stored
c    row-by-row in the array a.  to identify the individual nonzero
c    entries in each row, we need to know in which column each entry
c    lies.  the column indices which correspond to the nonzero entries
c    of m are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
c    ja(k) = j.  in addition, we need to know where each row starts and
c    how long it is.  the index positions in ja and a where the rows of
c    m begin are stored in the array ia.  i.e., if m(i,j) is the first
c    (leftmost) entry in the i-th row and  a(k) = m(i,j),  then
c    ia(i) = k.  moreover, the index in ja and a of the first location
c    following the last element in the last row is stored in ia(n+1).
c    thus, the number of entries in the i-th row is given by
c    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
c    consecutively in
c            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
c    and the corresponding column indices are stored consecutively in
c            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
c    for example, the 5 by 5 matrix
c                ( 1. 0. 2. 0. 0.)
c                ( 0. 3. 0. 0. 0.)
c            m = ( 0. 4. 5. 6. 0.)
c                ( 0. 0. 0. 7. 0.)
c                ( 0. 0. 0. 8. 9.)
c    would be stored as
c               - 1  2  3  4  5  6  7  8  9
c            ---+--------------------------
c            ia - 1  3  4  7  8 10
c            ja - 1  3  2  2  3  4  4  4  5
c             a - 1. 2. 3. 4. 5. 6. 7. 8. 9.         .
c
c         the strict upper (lower) triangular portion of the matrix
c    u (l) is stored in a similar fashion using the arrays  iu, ju, u
c    (il, jl, l)  except that an additional array iju (ijl) is used to
c    compress storage of ju (jl) by allowing some sequences of column
c    (row) indices to used for more than one row (column)  (n.b., l is
c    stored by columns).  iju(k) (ijl(k)) points to the starting
c    location in ju (jl) of entries for the kth row (column).
c    compression in ju (jl) occurs in two ways.  first, if a row
c    (column) i was merged into the current row (column) k, and the
c    number of elements merged in from (the tail portion of) row
c    (column) i is the same as the final length of row (column) k, then
c    the kth row (column) and the tail of row (column) i are identical
c    and iju(k) (ijl(k)) points to the start of the tail.  second, if
c    some tail portion of the (k-1)st row (column) is identical to the
c    head of the kth row (column), then iju(k) (ijl(k)) points to the
c    start of that tail portion.  for example, the nonzero structure of
c    the strict upper triangular part of the matrix
c            d 0 x x x
c            0 d 0 x x
c            0 0 d x 0
c            0 0 0 d x
c            0 0 0 0 d
c    would be represented as
c                - 1 2 3 4 5 6
c            ----+------------
c             iu - 1 4 6 7 8 8
c             ju - 3 4 5 4
c            iju - 1 2 4 3           .
c    the diagonal entries of l and u are assumed to be equal to one and
c    are not stored.  the array d contains the reciprocals of the
c    diagonal entries of the matrix d.
c
c    iii. additional storage savings
c         in nsfc, r and ic can be the same array in the calling
c    sequence if no reordering of the coefficient matrix has been done.
c         in nnfc, r, c, and ic can all be the same array if no
c    reordering has been done.  if only the rows have been reordered,
c    then c and ic can be the same array.  if the row and column
c    orderings are the same, then r and c can be the same array.  z and
c    row can be the same array.
c         in nnsc or nntc, r and c can be the same array if no
c    reordering has been done or if the row and column orderings are the
c    same.  z and b can be the same array.  however, then b will be
c    destroyed.
c
c    iv.  parameters
c         following is a list of parameters to the programs.  names are
c    uniform among the various subroutines.  class abbreviations are --
c       n - integer variable
c       f - real variable
c       v - supplies a value to a subroutine
c       r - returns a result from a subroutine
c       i - used internally by a subroutine
c       a - array
c
c class - parameter
c ------+----------
c fva   - a     - nonzero entries of the coefficient matrix m, stored
c       -           by rows.
c       -           size = number of nonzero entries in m.
c fva   - b     - right-hand side b.
c       -           size = n.
c nva   - c     - ordering of the columns of m.
c       -           size = n.
c fvra  - d     - reciprocals of the diagonal entries of the matrix d.
c       -           size = n.
c nr    - flag  - error flag.  values and their meanings are --
c       -            0     no errors detected
c       -            n+k   null row in a  --  row = k
c       -           2n+k   duplicate entry in a  --  row = k
c       -           3n+k   insufficient storage for jl  --  row = k
c       -           4n+1   insufficient storage for l
c       -           5n+k   null pivot  --  row = k
c       -           6n+k   insufficient storage for ju  --  row = k
c       -           7n+1   insufficient storage for u
c       -           8n+k   zero pivot  --  row = k
c nva   - ia    - pointers to delimit the rows of a.
c       -           size = n+1.
c nvra  - ijl   - pointers to the first element in each column in jl,
c       -           used to compress storage in jl.
c       -           size = n.
c nvra  - iju   - pointers to the first element in each row in ju, used
c       -           to compress storage in ju.
c       -           size = n.
c nvra  - il    - pointers to delimit the columns of l.
c       -           size = n+1.
c nvra  - iu    - pointers to delimit the rows of u.
c       -           size = n+1.
c nva   - ja    - column numbers corresponding to the elements of a.
c       -           size = size of a.
c nvra  - jl    - row numbers corresponding to the elements of l.
c       -           size = jlmax.
c nv    - jlmax - declared dimension of jl.  jlmax must be larger than
c       -           the number of nonzeros in the strict lower triangle
c       -           of m plus fillin minus compression.
c nvra  - ju    - column numbers corresponding to the elements of u.
c       -           size = jumax.
c nv    - jumax - declared dimension of ju.  jumax must be larger than
c       -           the number of nonzeros in the strict upper triangle
c       -           of m plus fillin minus compression.
c fvra  - l     - nonzero entries in the strict lower triangular portion
c       -           of the matrix l, stored by columns.
c       -           size = lmax.
c nv    - lmax  - declared dimension of l.  lmax must be larger than
c       -           the number of nonzeros in the strict lower triangle
c       -           of m plus fillin  (il(n+1)-1 after nsfc).
c nv    - n     - number of variables/equations.
c nva   - r     - ordering of the rows of m.
c       -           size = n.
c fvra  - u     - nonzero entries in the strict upper triangular portion
c       -           of the matrix u, stored by rows.
c       -           size = umax.
c nv    - umax  - declared dimension of u.  umax must be larger than
c       -           the number of nonzeros in the strict upper triangle
c       -           of m plus fillin  (iu(n+1)-1 after nsfc).
c fra   - z     - solution x.
c       -           size = n.
c
c       ----------------------------------------------------------------
c
c*** subroutine nroc
c*** reorders rows of a, leaving row order unchanged
c
c
c       input parameters.. n, ic, ia, ja, a
c       output parameters.. ja, a, flag
c
c       parameters used internally..
c nia   - p     - at the kth step, p is a linked list of the reordered
c       -           column indices of the kth row of a.  p(n+1) points
c       -           to the first entry in the list.
c       -           size = n+1.
c nia   - jar   - at the kth step,jar contains the elements of the
c       -           reordered column indices of a.
c       -           size = n.
c fia   - ar    - at the kth step, ar contains the elements of the
c       -           reordered row of a.
c       -           size = n.
c
      integer  ic(*), ia(*), ja(*), jar(*), p(*), flag
c     real  a(*), ar(*)
      double precision  a(*), ar(*)
c
c  ******  for each nonempty row  *******************************
      do 5 k=1,n
        jmin = ia(k)
        jmax = ia(k+1) - 1
        if(jmin .gt. jmax) go to 5
        p(n+1) = n + 1
c  ******  insert each element in the list  *********************
        do 3 j=jmin,jmax
          newj = ic(ja(j))
          i = n + 1
   1      if(p(i) .ge. newj) go to 2
            i = p(i)
            go to 1
   2      if(p(i) .eq. newj) go to 102
          p(newj) = p(i)
          p(i) = newj
          jar(newj) = ja(j)
          ar(newj) = a(j)
   3      continue
c  ******  replace old row in ja and a  *************************
        i = n + 1
        do 4 j=jmin,jmax
          i = p(i)
          ja(j) = jar(i)
   4      a(j) = ar(i)
   5    continue
      flag = 0
      return
c
c ** error.. duplicate entry in a
 102  flag = n + k
      return
      end
      subroutine nsfc
     *      (n, r, ic, ia,ja, jlmax,il,jl,ijl, jumax,iu,ju,iju,
     *       q, ira,jra, irac, irl,jrl, iru,jru, flag)
c*** subroutine nsfc
c*** symbolic ldu-factorization of nonsymmetric sparse matrix
c      (compressed pointer storage)
c
c
c       input variables.. n, r, ic, ia, ja, jlmax, jumax.
c       output variables.. il, jl, ijl, iu, ju, iju, flag.
c
c       parameters used internally..
c nia   - q     - suppose  m*  is the result of reordering  m.  if
c       -           processing of the ith row of  m*  (hence the ith
c       -           row of  u) is being done,  q(j)  is initially
c       -           nonzero if  m*(i,j) is nonzero (j.ge.i).  since
c       -           values need not be stored, each entry points to the
c       -           next nonzero and  q(n+1)  points to the first.  n+1
c       -           indicates the end of the list.  for example, if n=9
c       -           and the 5th row of  m*  is
c       -              0 x x 0 x 0 0 x 0
c       -           then  q  will initially be
c       -              a a a a 8 a a 10 5           (a - arbitrary).
c       -           as the algorithm proceeds, other elements of  q
c       -           are inserted in the list because of fillin.
c       -           q  is used in an analogous manner to compute the
c       -           ith column of  l.
c       -           size = n+1.
c nia   - ira,  - vectors used to find the columns of  m.  at the kth
c nia   - jra,      step of the factorization,  irac(k)  points to the
c nia   - irac      head of a linked list in  jra  of row indices i
c       -           such that i .ge. k and  m(i,k)  is nonzero.  zero
c       -           indicates the end of the list.  ira(i)  (i.ge.k)
c       -           points to the smallest j such that j .ge. k and
c       -           m(i,j)  is nonzero.
c       -           size of each = n.
c nia   - irl,  - vectors used to find the rows of  l.  at the kth step
c nia   - jrl       of the factorization,  jrl(k)  points to the head
c       -           of a linked list in  jrl  of column indices j
c       -           such j .lt. k and  l(k,j)  is nonzero.  zero
c       -           indicates the end of the list.  irl(j)  (j.lt.k)
c       -           points to the smallest i such that i .ge. k and
c       -           l(i,j)  is nonzero.
c       -           size of each = n.
c nia   - iru,  - vectors used in a manner analogous to  irl and jrl
c nia   - jru       to find the columns of  u.
c       -           size of each = n.
c
c  internal variables..
c    jlptr - points to the last position used in  jl.
c    juptr - points to the last position used in  ju.
c    jmin,jmax - are the indices in  a or u  of the first and last
c                elements to be examined in a given row.
c                for example,  jmin=ia(k), jmax=ia(k+1)-1.
c
      integer cend, qm, rend, rk, vj
      integer ia(*), ja(*), ira(*), jra(*), il(*), jl(*), ijl(*)
      integer iu(*), ju(*), iju(*), irl(*), jrl(*), iru(*), jru(*)
      integer r(*), ic(*), q(*), irac(*), flag
c
c  ******  initialize pointers  ****************************************
      np1 = n + 1
      jlmin = 1
      jlptr = 0
      il(1) = 1
      jumin = 1
      juptr = 0
      iu(1) = 1
      do 1 k=1,n
        irac(k) = 0
        jra(k) = 0
        jrl(k) = 0
   1    jru(k) = 0
c  ******  initialize column pointers for a  ***************************
      do 2 k=1,n
        rk = r(k)
        iak = ia(rk)
        if (iak .ge. ia(rk+1))  go to 101
        jaiak = ic(ja(iak))
        if (jaiak .gt. k)  go to 105
        jra(k) = irac(jaiak)
        irac(jaiak) = k
   2    ira(k) = iak
c
c  ******  for each column of l and row of u  **************************
      do 41 k=1,n
c
c  ******  initialize q for computing kth column of l  *****************
        q(np1) = np1
        luk = -1
c  ******  by filling in kth column of a  ******************************
        vj = irac(k)
        if (vj .eq. 0)  go to 5
   3      qm = np1
   4      m = qm
          qm =  q(m)
          if (qm .lt. vj)  go to 4
          if (qm .eq. vj)  go to 102
            luk = luk + 1
            q(m) = vj
            q(vj) = qm
            vj = jra(vj)
            if (vj .ne. 0)  go to 3
c  ******  link through jru  *******************************************
   5    lastid = 0
        lasti = 0
        ijl(k) = jlptr
        i = k
   6      i = jru(i)
          if (i .eq. 0)  go to 10
          qm = np1
          jmin = irl(i)
          jmax = ijl(i) + il(i+1) - il(i) - 1
          long = jmax - jmin
          if (long .lt. 0)  go to 6
          jtmp = jl(jmin)
          if (jtmp .ne. k)  long = long + 1
          if (jtmp .eq. k)  r(i) = -r(i)
          if (lastid .ge. long)  go to 7
            lasti = i
            lastid = long
c  ******  and merge the corresponding columns into the kth column  ****
   7      do 9 j=jmin,jmax
            vj = jl(j)
   8        m = qm
            qm = q(m)
            if (qm .lt. vj)  go to 8
            if (qm .eq. vj)  go to 9
              luk = luk + 1
              q(m) = vj
              q(vj) = qm
              qm = vj
   9        continue
            go to 6
c  ******  lasti is the longest column merged into the kth  ************
c  ******  see if it equals the entire kth column  *********************
  10    qm = q(np1)
        if (qm .ne. k)  go to 105
        if (luk .eq. 0)  go to 17
        if (lastid .ne. luk)  go to 11
c  ******  if so, jl can be compressed  ********************************
        irll = irl(lasti)
        ijl(k) = irll + 1
        if (jl(irll) .ne. k)  ijl(k) = ijl(k) - 1
        go to 17
c  ******  if not, see if kth column can overlap the previous one  *****
  11    if (jlmin .gt. jlptr)  go to 15
        qm = q(qm)
        do 12 j=jlmin,jlptr
          if (jl(j) - qm)  12, 13, 15
  12      continue
        go to 15
  13    ijl(k) = j
        do 14 i=j,jlptr
          if (jl(i) .ne. qm)  go to 15
          qm = q(qm)
          if (qm .gt. n)  go to 17
  14      continue
        jlptr = j - 1
c  ******  move column indices from q to jl, update vectors  ***********
  15    jlmin = jlptr + 1
        ijl(k) = jlmin
        if (luk .eq. 0)  go to 17
        jlptr = jlptr + luk
        if (jlptr .gt. jlmax)  go to 103
          qm = q(np1)
          do 16 j=jlmin,jlptr
            qm = q(qm)
  16        jl(j) = qm
  17    irl(k) = ijl(k)
        il(k+1) = il(k) + luk
c
c  ******  initialize q for computing kth row of u  ********************
        q(np1) = np1
        luk = -1
c  ******  by filling in kth row of reordered a  ***********************
        rk = r(k)
        jmin = ira(k)
        jmax = ia(rk+1) - 1
        if (jmin .gt. jmax)  go to 20
        do 19 j=jmin,jmax
          vj = ic(ja(j))
          qm = np1
  18      m = qm
          qm = q(m)
          if (qm .lt. vj)  go to 18
          if (qm .eq. vj)  go to 102
            luk = luk + 1
            q(m) = vj
            q(vj) = qm
  19      continue
c  ******  link through jrl,  ******************************************
  20    lastid = 0
        lasti = 0
        iju(k) = juptr
        i = k
        i1 = jrl(k)
  21      i = i1
          if (i .eq. 0)  go to 26
          i1 = jrl(i)
          qm = np1
          jmin = iru(i)
          jmax = iju(i) + iu(i+1) - iu(i) - 1
          long = jmax - jmin
          if (long .lt. 0)  go to 21
          jtmp = ju(jmin)
          if (jtmp .eq. k)  go to 22
c  ******  update irl and jrl, *****************************************
            long = long + 1
            cend = ijl(i) + il(i+1) - il(i)
            irl(i) = irl(i) + 1
            if (irl(i) .ge. cend)  go to 22
              j = jl(irl(i))
              jrl(i) = jrl(j)
              jrl(j) = i
  22      if (lastid .ge. long)  go to 23
            lasti = i
            lastid = long
c  ******  and merge the corresponding rows into the kth row  **********
  23      do 25 j=jmin,jmax
            vj = ju(j)
  24        m = qm
            qm = q(m)
            if (qm .lt. vj)  go to 24
            if (qm .eq. vj)  go to 25
              luk = luk + 1
              q(m) = vj
              q(vj) = qm
              qm = vj
  25        continue
          go to 21
c  ******  update jrl(k) and irl(k)  ***********************************
  26    if (il(k+1) .le. il(k))  go to 27
          j = jl(irl(k))
          jrl(k) = jrl(j)
          jrl(j) = k
c  ******  lasti is the longest row merged into the kth  ***************
c  ******  see if it equals the entire kth row  ************************
  27    qm = q(np1)
        if (qm .ne. k)  go to 105
        if (luk .eq. 0)  go to 34
        if (lastid .ne. luk)  go to 28
c  ******  if so, ju can be compressed  ********************************
        irul = iru(lasti)
        iju(k) = irul + 1
        if (ju(irul) .ne. k)  iju(k) = iju(k) - 1
        go to 34
c  ******  if not, see if kth row can overlap the previous one  ********
  28    if (jumin .gt. juptr)  go to 32
        qm = q(qm)
        do 29 j=jumin,juptr
          if (ju(j) - qm)  29, 30, 32
  29      continue
        go to 32
  30    iju(k) = j
        do 31 i=j,juptr
          if (ju(i) .ne. qm)  go to 32
          qm = q(qm)
          if (qm .gt. n)  go to 34
  31      continue
        juptr = j - 1
c  ******  move row indices from q to ju, update vectors  **************
  32    jumin = juptr + 1
        iju(k) = jumin
        if (luk .eq. 0)  go to 34
        juptr = juptr + luk
        if (juptr .gt. jumax)  go to 106
          qm = q(np1)
          do 33 j=jumin,juptr
            qm = q(qm)
  33        ju(j) = qm
  34    iru(k) = iju(k)
        iu(k+1) = iu(k) + luk
c
c  ******  update iru, jru  ********************************************
        i = k
  35      i1 = jru(i)
          if (r(i) .lt. 0)  go to 36
          rend = iju(i) + iu(i+1) - iu(i)
          if (iru(i) .ge. rend)  go to 37
            j = ju(iru(i))
            jru(i) = jru(j)
            jru(j) = i
            go to 37
  36      r(i) = -r(i)
  37      i = i1
          if (i .eq. 0)  go to 38
          iru(i) = iru(i) + 1
          go to 35
c
c  ******  update ira, jra, irac  **************************************
  38    i = irac(k)
        if (i .eq. 0)  go to 41
  39      i1 = jra(i)
          ira(i) = ira(i) + 1
          if (ira(i) .ge. ia(r(i)+1))  go to 40
          irai = ira(i)
          jairai = ic(ja(irai))
          if (jairai .gt. i)  go to 40
          jra(i) = irac(jairai)
          irac(jairai) = i
  40      i = i1
          if (i .ne. 0)  go to 39
  41    continue
c
      ijl(n) = jlptr
      iju(n) = juptr
      flag = 0
      return
c
c ** error.. null row in a
 101  flag = n + rk
      return
c ** error.. duplicate entry in a
 102  flag = 2*n + rk
      return
c ** error.. insufficient storage for jl
 103  flag = 3*n + k
      return
c ** error.. null pivot
 105  flag = 5*n + k
      return
c ** error.. insufficient storage for ju
 106  flag = 6*n + k
      return
      end
      subroutine nnfc
     *     (n, r,c,ic, ia,ja,a, z, b,
     *      lmax,il,jl,ijl,l, d, umax,iu,ju,iju,u,
     *      row, tmp, irl,jrl, flag)
c*** subroutine nnfc
c*** numerical ldu-factorization of sparse nonsymmetric matrix and
c      solution of system of linear equations (compressed pointer
c      storage)
c
c
c       input variables..  n, r, c, ic, ia, ja, a, b,
c                          il, jl, ijl, lmax, iu, ju, iju, umax
c       output variables.. z, l, d, u, flag
c
c       parameters used internally..
c nia   - irl,  - vectors used to find the rows of  l.  at the kth step
c nia   - jrl       of the factorization,  jrl(k)  points to the head
c       -           of a linked list in  jrl  of column indices j
c       -           such j .lt. k and  l(k,j)  is nonzero.  zero
c       -           indicates the end of the list.  irl(j)  (j.lt.k)
c       -           points to the smallest i such that i .ge. k and
c       -           l(i,j)  is nonzero.
c       -           size of each = n.
c fia   - row   - holds intermediate values in calculation of  u and l.
c       -           size = n.
c fia   - tmp   - holds new right-hand side  b*  for solution of the
c       -           equation ux = b*.
c       -           size = n.
c
c  internal variables..
c    jmin, jmax - indices of the first and last positions in a row to
c      be examined.
c    sum - used in calculating  tmp.
c
      integer rk,umax
      integer  r(*), c(*), ic(*), ia(*), ja(*), il(*), jl(*), ijl(*)
      integer  iu(*), ju(*), iju(*), irl(*), jrl(*), flag
c     real  a(*), l(*), d(*), u(*), z(*), b(*), row(*)
c     real tmp(*), lki, sum, dk
      double precision  a(*), l(*), d(*), u(*), z(*), b(*), row(*)
      double precision  tmp(*), lki, sum, dk
c
c  ******  initialize pointers and test storage  ***********************
      if(il(n+1)-1 .gt. lmax) go to 104
      if(iu(n+1)-1 .gt. umax) go to 107
      do 1 k=1,n
        irl(k) = il(k)
        jrl(k) = 0
   1    continue
c
c  ******  for each row  ***********************************************
      do 19 k=1,n
c  ******  reverse jrl and zero row where kth row of l will fill in  ***
        row(k) = 0
        i1 = 0
        if (jrl(k) .eq. 0) go to 3
        i = jrl(k)
   2    i2 = jrl(i)
        jrl(i) = i1
        i1 = i
        row(i) = 0
        i = i2
        if (i .ne. 0) go to 2
c  ******  set row to zero where u will fill in  ***********************
   3    jmin = iju(k)
        jmax = jmin + iu(k+1) - iu(k) - 1
        if (jmin .gt. jmax) go to 5
        do 4 j=jmin,jmax
   4      row(ju(j)) = 0
c  ******  place kth row of a in row  **********************************
   5    rk = r(k)
        jmin = ia(rk)
        jmax = ia(rk+1) - 1
        do 6 j=jmin,jmax
          row(ic(ja(j))) = a(j)
   6      continue
c  ******  initialize sum, and link through jrl  ***********************
        sum = b(rk)
        i = i1
        if (i .eq. 0) go to 10
c  ******  assign the kth row of l and adjust row, sum  ****************
   7      lki = -row(i)
c  ******  if l is not required, then comment out the following line  **
          l(irl(i)) = -lki
          sum = sum + lki * tmp(i)
          jmin = iu(i)
          jmax = iu(i+1) - 1
          if (jmin .gt. jmax) go to 9
          mu = iju(i) - jmin
          do 8 j=jmin,jmax
   8        row(ju(mu+j)) = row(ju(mu+j)) + lki * u(j)
   9      i = jrl(i)
          if (i .ne. 0) go to 7
c
c  ******  assign kth row of u and diagonal d, set tmp(k)  *************
  10    if (row(k) .eq. 0.0d0) go to 108
        dk = 1.0d0 / row(k)
        d(k) = dk
        tmp(k) = sum * dk
        if (k .eq. n) go to 19
        jmin = iu(k)
        jmax = iu(k+1) - 1
        if (jmin .gt. jmax)  go to 12
        mu = iju(k) - jmin
        do 11 j=jmin,jmax
  11      u(j) = row(ju(mu+j)) * dk
  12    continue
c
c  ******  update irl and jrl, keeping jrl in decreasing order  ********
        i = i1
        if (i .eq. 0) go to 18
  14    irl(i) = irl(i) + 1
        i1 = jrl(i)
        if (irl(i) .ge. il(i+1)) go to 17
        ijlb = irl(i) - il(i) + ijl(i)
        j = jl(ijlb)
  15    if (i .gt. jrl(j)) go to 16
          j = jrl(j)
          go to 15
  16    jrl(i) = jrl(j)
        jrl(j) = i
  17    i = i1
        if (i .ne. 0) go to 14
  18    if (irl(k) .ge. il(k+1)) go to 19
        j = jl(ijl(k))
        jrl(k) = jrl(j)
        jrl(j) = k
  19    continue
c
c  ******  solve  ux = tmp  by back substitution  **********************
      k = n
      do 22 i=1,n
        sum =  tmp(k)
        jmin = iu(k)
        jmax = iu(k+1) - 1
        if (jmin .gt. jmax)  go to 21
        mu = iju(k) - jmin
        do 20 j=jmin,jmax
  20      sum = sum - u(j) * tmp(ju(mu+j))
  21    tmp(k) =  sum
        z(c(k)) =  sum
  22    k = k-1
      flag = 0
      return
c
c ** error.. insufficient storage for l
 104  flag = 4*n + 1
      return
c ** error.. insufficient storage for u
 107  flag = 7*n + 1
      return
c ** error.. zero pivot
 108  flag = 8*n + k
      return
      end
      subroutine nnsc
     *     (n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)
c*** subroutine nnsc
c*** numerical solution of sparse nonsymmetric system of linear
c      equations given ldu-factorization (compressed pointer storage)
c
c
c       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
c       output variables.. z
c
c       parameters used internally..
c fia   - tmp   - temporary vector which gets result of solving  ly = b.
c       -           size = n.
c
c  internal variables..
c    jmin, jmax - indices of the first and last positions in a row of
c      u or l  to be used.
c
      integer r(*), c(*), il(*), jl(*), ijl(*), iu(*), ju(*), iju(*)
c     real l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk, sum
      double precision  l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk,sum
c
c  ******  set tmp to reordered b  *************************************
      do 1 k=1,n
   1    tmp(k) = b(r(k))
c  ******  solve  ly = b  by forward substitution  *********************
      do 3 k=1,n
        jmin = il(k)
        jmax = il(k+1) - 1
        tmpk = -d(k) * tmp(k)
        tmp(k) = -tmpk
        if (jmin .gt. jmax) go to 3
        ml = ijl(k) - jmin
        do 2 j=jmin,jmax
   2      tmp(jl(ml+j)) = tmp(jl(ml+j)) + tmpk * l(j)
   3    continue
c  ******  solve  ux = y  by back substitution  ************************
      k = n
      do 6 i=1,n
        sum = -tmp(k)
        jmin = iu(k)
        jmax = iu(k+1) - 1
        if (jmin .gt. jmax) go to 5
        mu = iju(k) - jmin
        do 4 j=jmin,jmax
   4      sum = sum + u(j) * tmp(ju(mu+j))
   5    tmp(k) = -sum
        z(c(k)) = -sum
        k = k - 1
   6    continue
      return
      end
      subroutine nntc
     *     (n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, z, b, tmp)
c*** subroutine nntc
c*** numeric solution of the transpose of a sparse nonsymmetric system
c      of linear equations given lu-factorization (compressed pointer
c      storage)
c
c
c       input variables..  n, r, c, il, jl, ijl, l, d, iu, ju, iju, u, b
c       output variables.. z
c
c       parameters used internally..
c fia   - tmp   - temporary vector which gets result of solving ut y = b
c       -           size = n.
c
c  internal variables..
c    jmin, jmax - indices of the first and last positions in a row of
c      u or l  to be used.
c
      integer r(*), c(*), il(*), jl(*), ijl(*), iu(*), ju(*), iju(*)
c     real l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk,sum
      double precision l(*), d(*), u(*), b(*), z(*), tmp(*), tmpk,sum
c
c  ******  set tmp to reordered b  *************************************
      do 1 k=1,n
   1    tmp(k) = b(c(k))
c  ******  solve  ut y = b  by forward substitution  *******************
      do 3 k=1,n
        jmin = iu(k)
        jmax = iu(k+1) - 1
        tmpk = -tmp(k)
        if (jmin .gt. jmax) go to 3
        mu = iju(k) - jmin
        do 2 j=jmin,jmax
   2      tmp(ju(mu+j)) = tmp(ju(mu+j)) + tmpk * u(j)
   3    continue
c  ******  solve  lt x = y  by back substitution  **********************
      k = n
      do 6 i=1,n
        sum = -tmp(k)
        jmin = il(k)
        jmax = il(k+1) - 1
        if (jmin .gt. jmax) go to 5
        ml = ijl(k) - jmin
        do 4 j=jmin,jmax
   4      sum = sum + l(j) * tmp(jl(ml+j))
   5    tmp(k) = -sum * d(k)
        z(r(k)) = tmp(k)
        k = k - 1
   6    continue
      return
      end


*DECK DMNORM
      DOUBLE PRECISION FUNCTION DMNORM (N, V, W)
C-----------------------------------------------------------------------
C This function routine computes the weighted max-norm
C of the vector of length N contained in the array V, with weights
C contained in the array w of length N:
C   DMNORM = MAX(i=1,...,N) ABS(V(i))*W(i)
C-----------------------------------------------------------------------
      INTEGER N,   I
      DOUBLE PRECISION V, W,   VM
      DIMENSION V(N), W(N)
      VM = 0.0D0
      DO 10 I = 1,N
 10     VM = MAX(VM,ABS(V(I))*W(I))
      DMNORM = VM
      RETURN
C----------------------- End of Function DMNORM ------------------------
      END
*DECK DFNORM
      DOUBLE PRECISION FUNCTION DFNORM (N, A, W)
C-----------------------------------------------------------------------
C This function computes the norm of a full N by N matrix,
C stored in the array A, that is consistent with the weighted max-norm
C on vectors, with weights stored in the array W:
C   DFNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
C-----------------------------------------------------------------------
      INTEGER N,   I, J
      DOUBLE PRECISION A,   W, AN, SUM
      DIMENSION A(N,N), W(N)
      AN = 0.0D0
      DO 20 I = 1,N
        SUM = 0.0D0
        DO 10 J = 1,N
 10       SUM = SUM + ABS(A(I,J))/W(J)
        AN = MAX(AN,SUM*W(I))
 20     CONTINUE
      DFNORM = AN
      RETURN
C----------------------- End of Function DFNORM ------------------------
      END
*DECK DBNORM
      DOUBLE PRECISION FUNCTION DBNORM (N, A, NRA, ML, MU, W)
C-----------------------------------------------------------------------
C This function computes the norm of a banded N by N matrix,
C stored in the array A, that is consistent with the weighted max-norm
C on vectors, with weights stored in the array W.
C ML and MU are the lower and upper half-bandwidths of the matrix.
C NRA is the first dimension of the A array, NRA .ge. ML+MU+1.
C In terms of the matrix elements a(i,j), the norm is given by:
C   DBNORM = MAX(i=1,...,N) ( W(i) * Sum(j=1,...,N) ABS(a(i,j))/W(j) )
C-----------------------------------------------------------------------
      INTEGER N, NRA, ML, MU
      INTEGER I, I1, JLO, JHI, J
      DOUBLE PRECISION A, W
      DOUBLE PRECISION AN, SUM
      DIMENSION A(NRA,N), W(N)
      AN = 0.0D0
      DO 20 I = 1,N
        SUM = 0.0D0
        I1 = I + MU + 1
        JLO = MAX(I-ML,1)
        JHI = MIN(I+MU,N)
        DO 10 J = JLO,JHI
 10       SUM = SUM + ABS(A(I1-J,J))/W(J)
        AN = MAX(AN,SUM*W(I))
 20     CONTINUE
      DBNORM = AN
      RETURN
C----------------------- End of Function DBNORM ------------------------
      END
*DECK DSRCMA
      SUBROUTINE DSRCMA (RSAV, ISAV, JOB)
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLSA01, which are used
C internally by one or more ODEPACK solvers.
C
C RSAV = real array of length 240 or more.
C ISAV = integer array of length 46 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
      INTEGER ISAV, JOB
      INTEGER ILS, ILSA
      INTEGER I, LENRLS, LENILS, LENRLA, LENILA
      DOUBLE PRECISION RSAV
      DOUBLE PRECISION RLS, RLSA
      DIMENSION RSAV(*), ISAV(*)
      SAVE LENRLS, LENILS, LENRLA, LENILA
      COMMON /DLS001/ RLS(218), ILS(37)
      COMMON /DLSA01/ RLSA(22), ILSA(9)
      DATA LENRLS/218/, LENILS/37/, LENRLA/22/, LENILA/9/
C
      IF (JOB .EQ. 2) GO TO 100
      DO 10 I = 1,LENRLS
 10     RSAV(I) = RLS(I)
      DO 15 I = 1,LENRLA
 15     RSAV(LENRLS+I) = RLSA(I)
C
      DO 20 I = 1,LENILS
 20     ISAV(I) = ILS(I)
      DO 25 I = 1,LENILA
 25     ISAV(LENILS+I) = ILSA(I)
C
      RETURN
C
 100  CONTINUE
      DO 110 I = 1,LENRLS
 110     RLS(I) = RSAV(I)
      DO 115 I = 1,LENRLA
 115     RLSA(I) = RSAV(LENRLS+I)
C
      DO 120 I = 1,LENILS
 120     ILS(I) = ISAV(I)
      DO 125 I = 1,LENILA
 125     ILSA(I) = ISAV(LENILS+I)
C
      RETURN
C----------------------- End of Subroutine DSRCMA ----------------------
      END

*DECK DSRCAR
      SUBROUTINE DSRCAR (RSAV, ISAV, JOB)
C-----------------------------------------------------------------------
C This routine saves or restores (depending on JOB) the contents of
C the Common blocks DLS001, DLSA01, DLSR01, which are used
C internally by one or more ODEPACK solvers.
C
C RSAV = real array of length 245 or more.
C ISAV = integer array of length 55 or more.
C JOB  = flag indicating to save or restore the Common blocks:
C        JOB  = 1 if Common is to be saved (written to RSAV/ISAV)
C        JOB  = 2 if Common is to be restored (read from RSAV/ISAV)
C        A call with JOB = 2 presumes a prior call with JOB = 1.
C-----------------------------------------------------------------------
      INTEGER ISAV, JOB
      INTEGER ILS, ILSA, ILSR
      INTEGER I, IOFF, LENRLS, LENILS, LENRLA, LENILA, LENRLR, LENILR
      DOUBLE PRECISION RSAV
      DOUBLE PRECISION RLS, RLSA, RLSR
      DIMENSION RSAV(*), ISAV(*)
      SAVE LENRLS, LENILS, LENRLA, LENILA, LENRLR, LENILR
      COMMON /DLS001/ RLS(218), ILS(37)
      COMMON /DLSA01/ RLSA(22), ILSA(9)
      COMMON /DLSR01/ RLSR(5), ILSR(9)
      DATA LENRLS/218/, LENILS/37/, LENRLA/22/, LENILA/9/
      DATA LENRLR/5/, LENILR/9/
C
      IF (JOB .EQ. 2) GO TO 100
      DO 10 I = 1,LENRLS
 10     RSAV(I) = RLS(I)
       DO 15 I = 1,LENRLA
 15     RSAV(LENRLS+I) = RLSA(I)
      IOFF = LENRLS + LENRLA
      DO 20 I = 1,LENRLR
 20     RSAV(IOFF+I) = RLSR(I)
C
      DO 30 I = 1,LENILS
 30     ISAV(I) = ILS(I)
      DO 35 I = 1,LENILA
 35     ISAV(LENILS+I) = ILSA(I)
      IOFF = LENILS + LENILA
      DO 40 I = 1,LENILR
 40     ISAV(IOFF+I) = ILSR(I)
C
      RETURN
C
 100  CONTINUE
      DO 110 I = 1,LENRLS
 110     RLS(I) = RSAV(I)
      DO 115 I = 1,LENRLA
 115     RLSA(I) = RSAV(LENRLS+I)
      IOFF = LENRLS + LENRLA
      DO 120 I = 1,LENRLR
 120     RLSR(I) = RSAV(IOFF+I)
C
      DO 130 I = 1,LENILS
 130     ILS(I) = ISAV(I)
      DO 135 I = 1,LENILA
 135     ILSA(I) = ISAV(LENILS+I)
      IOFF = LENILS + LENILA
      DO 140 I = 1,LENILR
 140     ILSR(I) = ISAV(IOFF+I)
C
      RETURN
C----------------------- End of Subroutine DSRCAR ----------------------
      END

*DECK DHEFA
      SUBROUTINE DHEFA (A, LDA, N, IPVT, INFO, JOB)
      INTEGER LDA, N, IPVT(*), INFO, JOB
      DOUBLE PRECISION A(LDA,*)
C-----------------------------------------------------------------------
C     This routine is a modification of the LINPACK routine DGEFA and
C     performs an LU decomposition of an upper Hessenberg matrix A.
C     There are two options available:
C
C          (1)  performing a fresh factorization
C          (2)  updating the LU factors by adding a row and a
C               column to the matrix A.
C-----------------------------------------------------------------------
C     DHEFA factors an upper Hessenberg matrix by elimination.
C
C     On entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be factored.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        JOB     INTEGER
C                JOB = 1    means that a fresh factorization of the
C                           matrix A is desired.
C                JOB .ge. 2 means that the current factorization of A
C                           will be updated by the addition of a row
C                           and a column.
C
C     On return
C
C        A       an upper triangular matrix and the multipliers
C                which were used to obtain it.
C                The factorization can be written  A = L*U  where
C                L  is a product of permutation and unit lower
C                triangular matrices and  U  is upper triangular.
C
C        IPVT    INTEGER(N)
C                an integer vector of pivot indices.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = k  if  U(k,k) .eq. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DHESL will divide by zero if called.
C
C     Modification of LINPACK, by Peter Brown, LLNL.
C     Written 7/20/83.  This version dated 6/20/01.
C
C     BLAS called: DAXPY, IDAMAX
C-----------------------------------------------------------------------
      INTEGER IDAMAX, J, K, KM1, KP1, L, NM1
      DOUBLE PRECISION T
C
      IF (JOB .GT. 1) GO TO 80
C
C A new facorization is desired.  This is essentially the LINPACK
C code with the exception that we know there is only one nonzero
C element below the main diagonal.
C
C     Gaussian elimination with partial pivoting
C
      INFO = 0
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 K = 1, NM1
         KP1 = K + 1
C
C        Find L = pivot index
C
         L = IDAMAX (2, A(K,K), 1) + K - 1
         IPVT(K) = L
C
C        Zero pivot implies this column already triangularized
C
         IF (A(L,K) .EQ. 0.0D0) GO TO 40
C
C           Interchange if necessary
C
            IF (L .EQ. K) GO TO 10
               T = A(L,K)
               A(L,K) = A(K,K)
               A(K,K) = T
   10       CONTINUE
C
C           Compute multipliers
C
            T = -1.0D0/A(K,K)
            A(K+1,K) = A(K+1,K)*T
C
C           Row elimination with column indexing
C
            DO 30 J = KP1, N
               T = A(L,J)
               IF (L .EQ. K) GO TO 20
                  A(L,J) = A(K,J)
                  A(K,J) = T
   20          CONTINUE
               CALL DAXPY (N-K, T, A(K+1,K), 1, A(K+1,J), 1)
   30       CONTINUE
         GO TO 50
   40    CONTINUE
            INFO = K
   50    CONTINUE
   60 CONTINUE
   70 CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
C
C The old factorization of A will be updated.  A row and a column
C has been added to the matrix A.
C N-1 is now the old order of the matrix.
C
  80  CONTINUE
      NM1 = N - 1
C
C Perform row interchanges on the elements of the new column, and
C perform elimination operations on the elements using the multipliers.
C
      IF (NM1 .LE. 1) GO TO 105
      DO 100 K = 2,NM1
        KM1 = K - 1
        L = IPVT(KM1)
        T = A(L,N)
        IF (L .EQ. KM1) GO TO 90
          A(L,N) = A(KM1,N)
          A(KM1,N) = T
  90    CONTINUE
        A(K,N) = A(K,N) + A(K,KM1)*T
 100    CONTINUE
 105  CONTINUE
C
C Complete update of factorization by decomposing last 2x2 block.
C
      INFO = 0
C
C        Find L = pivot index
C
         L = IDAMAX (2, A(NM1,NM1), 1) + NM1 - 1
         IPVT(NM1) = L
C
C        Zero pivot implies this column already triangularized
C
         IF (A(L,NM1) .EQ. 0.0D0) GO TO 140
C
C           Interchange if necessary
C
            IF (L .EQ. NM1) GO TO 110
               T = A(L,NM1)
               A(L,NM1) = A(NM1,NM1)
               A(NM1,NM1) = T
  110       CONTINUE
C
C           Compute multipliers
C
            T = -1.0D0/A(NM1,NM1)
            A(N,NM1) = A(N,NM1)*T
C
C           Row elimination with column indexing
C
               T = A(L,N)
               IF (L .EQ. NM1) GO TO 120
                  A(L,N) = A(NM1,N)
                  A(NM1,N) = T
  120          CONTINUE
               A(N,N) = A(N,N) + T*A(N,NM1)
         GO TO 150
  140    CONTINUE
            INFO = NM1
  150    CONTINUE
      IPVT(N) = N
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
C----------------------- End of Subroutine DHEFA -----------------------
      END
*DECK DHESL
      SUBROUTINE DHESL (A, LDA, N, IPVT, B)
      INTEGER LDA, N, IPVT(*)
      DOUBLE PRECISION A(LDA,*), B(*)
C-----------------------------------------------------------------------
C This is essentially the LINPACK routine DGESL except for changes
C due to the fact that A is an upper Hessenberg matrix.
C-----------------------------------------------------------------------
C     DHESL solves the real system A * x = b
C     using the factors computed by DHEFA.
C
C     On entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DHEFA.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                the order of the matrix  A .
C
C        IPVT    INTEGER(N)
C                the pivot vector from DHEFA.
C
C        B       DOUBLE PRECISION(N)
C                the right hand side vector.
C
C     On return
C
C        B       the solution vector  x .
C
C     Modification of LINPACK, by Peter Brown, LLNL.
C     Written 7/20/83.  This version dated 6/20/01.
C
C     BLAS called: DAXPY
C-----------------------------------------------------------------------
      INTEGER K, KB, L, NM1
      DOUBLE PRECISION T
C
      NM1 = N - 1
C
C        Solve  A * x = b
C        First solve  L*y = b
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            B(K+1) = B(K+1) + T*A(K+1,K)
   20    CONTINUE
   30    CONTINUE
C
C        Now solve  U*x = y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY (K-1, T, A(1,K), 1, B(1), 1)
   40    CONTINUE
      RETURN
C----------------------- End of Subroutine DHESL -----------------------
      END
*DECK DHEQR
      SUBROUTINE DHEQR (A, LDA, N, Q, INFO, IJOB)
      INTEGER LDA, N, INFO, IJOB
      DOUBLE PRECISION A(LDA,*), Q(*)
C-----------------------------------------------------------------------
C     This routine performs a QR decomposition of an upper
C     Hessenberg matrix A.  There are two options available:
C
C          (1)  performing a fresh decomposition
C          (2)  updating the QR factors by adding a row and a
C               column to the matrix A.
C-----------------------------------------------------------------------
C     DHEQR decomposes an upper Hessenberg matrix by using Givens
C     rotations.
C
C     On entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the matrix to be decomposed.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                A is an (N+1) by N Hessenberg matrix.
C
C        IJOB    INTEGER
C                = 1     means that a fresh decomposition of the
C                        matrix A is desired.
C                .ge. 2  means that the current decomposition of A
C                        will be updated by the addition of a row
C                        and a column.
C     On return
C
C        A       the upper triangular matrix R.
C                The factorization can be written Q*A = R, where
C                Q is a product of Givens rotations and R is upper
C                triangular.
C
C        Q       DOUBLE PRECISION(2*N)
C                the factors c and s of each Givens rotation used
C                in decomposing A.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = k  if  A(k,k) .eq. 0.0 .  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DHELS will divide by zero
C                     if called.
C
C     Modification of LINPACK, by Peter Brown, LLNL.
C     Written 1/13/86.  This version dated 6/20/01.
C-----------------------------------------------------------------------
      INTEGER I, IQ, J, K, KM1, KP1, NM1
      DOUBLE PRECISION C, S, T, T1, T2
C
      IF (IJOB .GT. 1) GO TO 70
C
C A new facorization is desired.
C
C     QR decomposition without pivoting
C
      INFO = 0
      DO 60 K = 1, N
         KM1 = K - 1
         KP1 = K + 1
C
C           Compute kth column of R.
C           First, multiply the kth column of A by the previous
C           k-1 Givens rotations.
C
            IF (KM1 .LT. 1) GO TO 20
            DO 10 J = 1, KM1
              I = 2*(J-1) + 1
              T1 = A(J,K)
              T2 = A(J+1,K)
              C = Q(I)
              S = Q(I+1)
              A(J,K) = C*T1 - S*T2
              A(J+1,K) = S*T1 + C*T2
   10         CONTINUE
C
C           Compute Givens components c and s
C
   20       CONTINUE
            IQ = 2*KM1 + 1
            T1 = A(K,K)
            T2 = A(KP1,K)
            IF (T2 .NE. 0.0D0) GO TO 30
              C = 1.0D0
              S = 0.0D0
              GO TO 50
   30       CONTINUE
            IF (ABS(T2) .LT. ABS(T1)) GO TO 40
              T = T1/T2
              S = -1.0D0/SQRT(1.0D0+T*T)
              C = -S*T
              GO TO 50
   40       CONTINUE
              T = T2/T1
              C = 1.0D0/SQRT(1.0D0+T*T)
              S = -C*T
   50       CONTINUE
            Q(IQ) = C
            Q(IQ+1) = S
            A(K,K) = C*T1 - S*T2
            IF (A(K,K) .EQ. 0.0D0) INFO = K
   60 CONTINUE
      RETURN
C
C The old factorization of A will be updated.  A row and a column
C has been added to the matrix A.
C N by N-1 is now the old size of the matrix.
C
  70  CONTINUE
      NM1 = N - 1
C
C Multiply the new column by the N previous Givens rotations.
C
      DO 100 K = 1,NM1
        I = 2*(K-1) + 1
        T1 = A(K,N)
        T2 = A(K+1,N)
        C = Q(I)
        S = Q(I+1)
        A(K,N) = C*T1 - S*T2
        A(K+1,N) = S*T1 + C*T2
 100    CONTINUE
C
C Complete update of decomposition by forming last Givens rotation,
C and multiplying it times the column vector (A(N,N), A(N+1,N)).
C
      INFO = 0
      T1 = A(N,N)
      T2 = A(N+1,N)
      IF (T2 .NE. 0.0D0) GO TO 110
        C = 1.0D0
        S = 0.0D0
        GO TO 130
 110  CONTINUE
      IF (ABS(T2) .LT. ABS(T1)) GO TO 120
        T = T1/T2
        S = -1.0D0/SQRT(1.0D0+T*T)
        C = -S*T
        GO TO 130
 120  CONTINUE
        T = T2/T1
        C = 1.0D0/SQRT(1.0D0+T*T)
        S = -C*T
 130  CONTINUE
      IQ = 2*N - 1
      Q(IQ) = C
      Q(IQ+1) = S
      A(N,N) = C*T1 - S*T2
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
C----------------------- End of Subroutine DHEQR -----------------------
      END
*DECK DHELS
      SUBROUTINE DHELS (A, LDA, N, Q, B)
      INTEGER LDA, N
      DOUBLE PRECISION A(LDA,*), B(*), Q(*)
C-----------------------------------------------------------------------
C This is part of the LINPACK routine DGESL with changes
C due to the fact that A is an upper Hessenberg matrix.
C-----------------------------------------------------------------------
C     DHELS solves the least squares problem
C
C           min (b-A*x, b-A*x)
C
C     using the factors computed by DHEQR.
C
C     On entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                the output from DHEQR which contains the upper
C                triangular factor R in the QR decomposition of A.
C
C        LDA     INTEGER
C                the leading dimension of the array  A .
C
C        N       INTEGER
C                A is originally an (N+1) by N matrix.
C
C        Q       DOUBLE PRECISION(2*N)
C                The coefficients of the N givens rotations
C                used in the QR factorization of A.
C
C        B       DOUBLE PRECISION(N+1)
C                the right hand side vector.
C
C     On return
C
C        B       the solution vector  x .
C
C     Modification of LINPACK, by Peter Brown, LLNL.
C     Written 1/13/86.  This version dated 6/20/01.
C
C     BLAS called: DAXPY
C-----------------------------------------------------------------------
      INTEGER IQ, K, KB, KP1
      DOUBLE PRECISION C, S, T, T1, T2
C
C        Minimize (b-A*x, b-A*x)
C        First form Q*b.
C
         DO 20 K = 1, N
            KP1 = K + 1
            IQ = 2*(K-1) + 1
            C = Q(IQ)
            S = Q(IQ+1)
            T1 = B(K)
            T2 = B(KP1)
            B(K) = C*T1 - S*T2
            B(KP1) = S*T1 + C*T2
   20    CONTINUE
C
C        Now solve  R*x = Q*b.
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY (K-1, T, A(1,K), 1, B(1), 1)
   40    CONTINUE
      RETURN
C----------------------- End of Subroutine DHELS -----------------------
      END
*DECK DSLSBT
      SUBROUTINE DSLSBT (WM, IWM, X, TEM)
      INTEGER IWM
      INTEGER LBLOX, LPB, LPC, MB, NB
      DOUBLE PRECISION WM, X, TEM
      DIMENSION WM(*), IWM(*), X(*), TEM(*)
C-----------------------------------------------------------------------
C This routine acts as an interface between the core integrator
C routine and the DSOLBT routine for the solution of the linear system
C arising from chord iteration.
C Communication with DSLSBT uses the following variables:
C WM    = real work space containing the LU decomposition,
C         starting at WM(3).
C IWM   = integer work space containing pivot information, starting at
C         IWM(21).  IWM also contains block structure parameters
C         MB = IWM(1) and NB = IWM(2).
C X     = the right-hand side vector on input, and the solution vector
C         on output, of length N.
C TEM   = vector of work space of length N, not used in this version.
C-----------------------------------------------------------------------
      MB = IWM(1)
      NB = IWM(2)
      LBLOX = MB*MB*NB
      LPB = 3 + LBLOX
      LPC = LPB + LBLOX
      CALL DSOLBT (MB, NB, WM(3), WM(LPB), WM(LPC), X, IWM(21))
      RETURN
C----------------------- End of Subroutine DSLSBT ----------------------
      END
*DECK DDECBT
      SUBROUTINE DDECBT (M, N, A, B, C, IP, IER)
      INTEGER M, N, IP(M,N), IER
      DOUBLE PRECISION A(M,M,N), B(M,M,N), C(M,M,N)
C-----------------------------------------------------------------------
C Block-tridiagonal matrix decomposition routine.
C Written by A. C. Hindmarsh.
C Latest revision:  November 10, 1983 (ACH)
C Reference:  UCID-30150
C             Solution of Block-Tridiagonal Systems of Linear
C             Algebraic Equations
C             A.C. Hindmarsh
C             February 1977
C The input matrix contains three blocks of elements in each block-row,
C including blocks in the (1,3) and (N,N-2) block positions.
C DDECBT uses block Gauss elimination and Subroutines DGEFA and DGESL
C for solution of blocks.  Partial pivoting is done within
C block-rows only.
C
C Note: this version uses LINPACK routines DGEFA/DGESL instead of
C of dec/sol for solution of blocks, and it uses the BLAS routine DDOT
C for dot product calculations.
C
C Input:
C     M = order of each block.
C     N = number of blocks in each direction of the matrix.
C         N must be 4 or more.  The complete matrix has order M*N.
C     A = M by M by N array containing diagonal blocks.
C         A(i,j,k) contains the (i,j) element of the k-th block.
C     B = M by M by N array containing the super-diagonal blocks
C         (in B(*,*,k) for k = 1,...,N-1) and the block in the (N,N-2)
C         block position (in B(*,*,N)).
C     C = M by M by N array containing the subdiagonal blocks
C         (in C(*,*,k) for k = 2,3,...,N) and the block in the
C         (1,3) block position (in C(*,*,1)).
C    IP = integer array of length M*N for working storage.
C Output:
C A,B,C = M by M by N arrays containing the block-LU decomposition
C         of the input matrix.
C    IP = M by N array of pivot information.  IP(*,k) contains
C         information for the k-th digonal block.
C   IER = 0  if no trouble occurred, or
C       = -1 if the input value of M or N was illegal, or
C       = k  if a singular matrix was found in the k-th diagonal block.
C Use DSOLBT to solve the associated linear system.
C
C External routines required: DGEFA and DGESL (from LINPACK) and
C DDOT (from the BLAS, or Basic Linear Algebra package).
C-----------------------------------------------------------------------
      INTEGER NM1, NM2, KM1, I, J, K
      DOUBLE PRECISION DP, DDOT
      IF (M .LT. 1 .OR. N .LT. 4) GO TO 210
      NM1 = N - 1
      NM2 = N - 2
C Process the first block-row. -----------------------------------------
      CALL DGEFA (A, M, M, IP, IER)
      K = 1
      IF (IER .NE. 0) GO TO 200
      DO 10 J = 1,M
        CALL DGESL (A, M, M, IP, B(1,J,1), 0)
        CALL DGESL (A, M, M, IP, C(1,J,1), 0)
 10     CONTINUE
C Adjust B(*,*,2). -----------------------------------------------------
      DO 40 J = 1,M
        DO 30 I = 1,M
          DP = DDOT (M, C(I,1,2), M, C(1,J,1), 1)
          B(I,J,2) = B(I,J,2) - DP
 30       CONTINUE
 40     CONTINUE
C Main loop.  Process block-rows 2 to N-1. -----------------------------
      DO 100 K = 2,NM1
        KM1 = K - 1
        DO 70 J = 1,M
          DO 60 I = 1,M
            DP = DDOT (M, C(I,1,K), M, B(1,J,KM1), 1)
            A(I,J,K) = A(I,J,K) - DP
 60         CONTINUE
 70       CONTINUE
        CALL DGEFA (A(1,1,K), M, M, IP(1,K), IER)
        IF (IER .NE. 0) GO TO 200
        DO 80 J = 1,M
 80       CALL DGESL (A(1,1,K), M, M, IP(1,K), B(1,J,K), 0)
 100    CONTINUE
C Process last block-row and return. -----------------------------------
      DO 130 J = 1,M
        DO 120 I = 1,M
          DP = DDOT (M, B(I,1,N), M, B(1,J,NM2), 1)
          C(I,J,N) = C(I,J,N) - DP
 120      CONTINUE
 130    CONTINUE
      DO 160 J = 1,M
        DO 150 I = 1,M
          DP = DDOT (M, C(I,1,N), M, B(1,J,NM1), 1)
          A(I,J,N) = A(I,J,N) - DP
 150      CONTINUE
 160    CONTINUE
      CALL DGEFA (A(1,1,N), M, M, IP(1,N), IER)
      K = N
      IF (IER .NE. 0) GO TO 200
      RETURN
C Error returns. -------------------------------------------------------
 200  IER = K
      RETURN
 210  IER = -1
      RETURN
C----------------------- End of Subroutine DDECBT ----------------------
      END
*DECK DSOLBT
      SUBROUTINE DSOLBT (M, N, A, B, C, Y, IP)
      INTEGER M, N, IP(M,N)
      DOUBLE PRECISION A(M,M,N), B(M,M,N), C(M,M,N), Y(M,N)
C-----------------------------------------------------------------------
C Solution of block-tridiagonal linear system.
C Coefficient matrix must have been previously processed by DDECBT.
C M, N, A,B,C, and IP  must not have been changed since call to DDECBT.
C Written by A. C. Hindmarsh.
C Input:
C     M = order of each block.
C     N = number of blocks in each direction of matrix.
C A,B,C = M by M by N arrays containing block LU decomposition
C         of coefficient matrix from DDECBT.
C    IP = M by N integer array of pivot information from DDECBT.
C     Y = array of length M*N containg the right-hand side vector
C         (treated as an M by N array here).
C Output:
C     Y = solution vector, of length M*N.
C
C External routines required: DGESL (LINPACK) and DDOT (BLAS).
C-----------------------------------------------------------------------
C
      INTEGER NM1, NM2, I, K, KB, KM1, KP1
      DOUBLE PRECISION DP, DDOT
      NM1 = N - 1
      NM2 = N - 2
C Forward solution sweep. ----------------------------------------------
      CALL DGESL (A, M, M, IP, Y, 0)
      DO 30 K = 2,NM1
        KM1 = K - 1
        DO 20 I = 1,M
          DP = DDOT (M, C(I,1,K), M, Y(1,KM1), 1)
          Y(I,K) = Y(I,K) - DP
 20       CONTINUE
        CALL DGESL (A(1,1,K), M, M, IP(1,K), Y(1,K), 0)
 30     CONTINUE
      DO 50 I = 1,M
        DP = DDOT (M, C(I,1,N), M, Y(1,NM1), 1)
     1     + DDOT (M, B(I,1,N), M, Y(1,NM2), 1)
        Y(I,N) = Y(I,N) - DP
 50     CONTINUE
      CALL DGESL (A(1,1,N), M, M, IP(1,N), Y(1,N), 0)
C Backward solution sweep. ---------------------------------------------
      DO 80 KB = 1,NM1
        K = N - KB
        KP1 = K + 1
        DO 70 I = 1,M
          DP = DDOT (M, B(I,1,K), M, Y(1,KP1), 1)
          Y(I,K) = Y(I,K) - DP
 70       CONTINUE
 80     CONTINUE
      DO 100 I = 1,M
        DP = DDOT (M, C(I,1,1), M, Y(1,3), 1)
        Y(I,1) = Y(I,1) - DP
 100    CONTINUE
      RETURN
C----------------------- End of Subroutine DSOLBT ----------------------
      END


*DECK XERRWD
      SUBROUTINE XERRWD (MSG, NMES, NERR, LEVEL, NI, I1, I2, NR, R1, R2)
C***BEGIN PROLOGUE  XERRWD
C***SUBSIDIARY
C***PURPOSE  Write error message with values.
C***CATEGORY  R3C
C***TYPE      DOUBLE PRECISION (XERRWV-S, XERRWD-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  Subroutines XERRWD, XSETF, XSETUN, and the function routine IXSAV,
C  as given here, constitute a simplified version of the SLATEC error
C  handling package.
C
C  All arguments are input arguments.
C
C  MSG    = The message (character array).
C  NMES   = The length of MSG (number of characters).
C  NERR   = The error number (not used).
C  LEVEL  = The error level..
C           0 or 1 means recoverable (control returns to caller).
C           2 means fatal (run is aborted--see note below).
C  NI     = Number of integers (0, 1, or 2) to be printed with message.
C  I1,I2  = Integers to be printed, depending on NI.
C  NR     = Number of reals (0, 1, or 2) to be printed with message.
C  R1,R2  = Reals to be printed, depending on NR.
C
C  Note..  this routine is machine-dependent and specialized for use
C  in limited context, in the following ways..
C  1. The argument MSG is assumed to be of type CHARACTER, and
C     the message is printed with a format of (1X,A).
C  2. The message is assumed to take only one line.
C     Multi-line messages are generated by repeated calls.
C  3. If LEVEL = 2, control passes to the statement   STOP
C     to abort the run.  This statement may be machine-dependent.
C  4. R1 and R2 are assumed to be in double precision and are printed
C     in D21.13 format.
C
C***ROUTINES CALLED  IXSAV
C***REVISION HISTORY  (YYMMDD)
C   920831  DATE WRITTEN
C   921118  Replaced MFLGSV/LUNSAV by IXSAV. (ACH)
C   930329  Modified prologue to SLATEC format. (FNF)
C   930407  Changed MSG from CHARACTER*1 array to variable. (FNF)
C   930922  Minor cosmetic change. (FNF)
C***END PROLOGUE  XERRWD
C
C*Internal Notes:
C
C For a different default logical unit number, IXSAV (or a subsidiary
C routine that it calls) will need to be modified.
C For a different run-abort command, change the statement following
C statement 100 at the end.
C-----------------------------------------------------------------------
C Subroutines called by XERRWD.. None
C Function routine called by XERRWD.. IXSAV
C-----------------------------------------------------------------------
C**End
C
C  Declare arguments.
C
      DOUBLE PRECISION R1, R2
      INTEGER NMES, NERR, LEVEL, NI, I1, I2, NR
      CHARACTER*(*) MSG
C
C  Declare local variables.
C
      INTEGER LUNIT, IXSAV, MESFLG
C
C  Get logical unit number and message print flag.
C
C***FIRST EXECUTABLE STATEMENT  XERRWD
      LUNIT = IXSAV (1, 0, .FALSE.)
      MESFLG = IXSAV (2, 0, .FALSE.)
      IF (MESFLG .EQ. 0) GO TO 100
C
C  Write the message.
C
      WRITE (LUNIT,10)  MSG
 10   FORMAT(1X,A)
      IF (NI .EQ. 1) WRITE (LUNIT, 20) I1
 20   FORMAT(6X,'In above message,  I1 =',I10)
      IF (NI .EQ. 2) WRITE (LUNIT, 30) I1,I2
 30   FORMAT(6X,'In above message,  I1 =',I10,3X,'I2 =',I10)
      IF (NR .EQ. 1) WRITE (LUNIT, 40) R1
 40   FORMAT(6X,'In above message,  R1 =',D21.13)
      IF (NR .EQ. 2) WRITE (LUNIT, 50) R1,R2
 50   FORMAT(6X,'In above,  R1 =',D21.13,3X,'R2 =',D21.13)
C
C  Abort the run if LEVEL = 2.
C
 100  IF (LEVEL .NE. 2) RETURN
      STOP
C----------------------- End of Subroutine XERRWD ----------------------
      END

*DECK DEWSET
      SUBROUTINE DEWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)
C***BEGIN PROLOGUE  DEWSET
C***SUBSIDIARY
C***PURPOSE  Set error weight vector.
C***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D)
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C
C  This subroutine sets the error weight vector EWT according to
C      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N,
C  with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
C  depending on the value of ITOL.
C
C***SEE ALSO  DLSODE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791129  DATE WRITTEN
C   890501  Modified prologue to SLATEC/LDOC format.  (FNF)
C   890503  Minor cosmetic changes.  (FNF)
C   930809  Renamed to allow single/double precision versions. (ACH)
C***END PROLOGUE  DEWSET
C**End
      INTEGER N, ITOL
      INTEGER I
      DOUBLE PRECISION RTOL, ATOL, YCUR, EWT
      DIMENSION RTOL(*), ATOL(*), YCUR(N), EWT(N)
C
C***FIRST EXECUTABLE STATEMENT  DEWSET
      GO TO (10, 20, 30, 40), ITOL
 10   CONTINUE
      DO 15 I = 1,N
 15     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
      RETURN
 20   CONTINUE
      DO 25 I = 1,N
 25     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
      RETURN
 30   CONTINUE
      DO 35 I = 1,N
 35     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
      RETURN
 40   CONTINUE
      DO 45 I = 1,N
 45     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
      RETURN
C----------------------- END OF SUBROUTINE DEWSET ----------------------
      END

