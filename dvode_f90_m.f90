  module DVODE_F90_M
! This version is the December 2005 release.
! Last change: 01/01/08
! _____________________________________________________________________
! Working Precision
  implicit none
! Define the working precision for DVODE_F90. Change D0 to E0 in the
! next statement to convert to single precision.
  integer, parameter, private :: WP = kind(1.0D0)
! ______________________________________________________________________
! Overview

! The f77 ordinary differential equation solver VODE.f is applicable to
! nonstiff systems of odes and to stiff systems having dense or banded
! Jacobians. DVODE_F90 is a Fortran 90 extension of VODE.f. While
! retaining all of the features available in VODE.f, we have
! incorporated several new options in DVODE_F90 including:
!   1. the ability to solve stiff systems with sparse Jacobians
!   2. internal management of storage and work arrays
!   3. specification of options via optional keywords
!   4. the ability to perform root finding or "event detection"
!   5. various new diagnostic and warning messages
!   6. the ability to impose solution bounds
!   7. several specialized options for dealing with sparsity
! ______________________________________________________________________
! Version Information

! This is DVODE_F90, the double precision FORTRAN 90 extension of the
! f77 DVODE.f ordinary differential equation solver. This version uses
! MA28 for sparse Jacobians. This file and related information can be
! obtained at the following support page:
!
!     http://www.radford.edu/~thompson/vodef90web/
!
! We are indebted to Richard Cox (ORNL) for providing us with his
! implementation of MA28 in LSOD28.f (a variant of Alan Hindmarsh's
! lsodes.f). We are indebted to Alan Hindmarsh for numerous contributions.
! In particular, we borrowed liberally from the f77 solvers VODE.f,
! LSODAR.f, and LSODES.f while developing DVODE_F90. We are indebted
! to Doug Salane for providing us with his JACSP Jacobian routines.
!
! If you find a bug or encounter a problem with DVODE_F90, please
! contact one of us:
!    G.D. Byrne (gbyrne@wi.rr.com)
!    S. Thompson (thompson@radford.edu)
! A set of quick start instructions is provided below.
! ______________________________________________________________________
! Note on F90/F95 Compilers

! To date we have used DVODE_F90 successfully with all F90/F95 compilers
! to which we have access. In particular, we have used it with the Lahey
! F90 and Lahey-Fujitsu F95 compilers, the Compaq Visual F90 compiler,
! the g90 compiler, and the INTEL, Portland, Salford, and SUN compilers.
! It should be noted that compilers such as Salford's FTN95 complain
! about uninitialized arrays passed as subroutine arguments and the use of
! slices of two dimensional arrays as one dimensional vectors, and will
! not run using the strictest compiler options. It is perfectly safe to
! use the /-CHECK compiler option to avoid these FTN95 runtime checks.
! DVODE_F90 does not use any variable for numerical purposes until it
! has been assigned an appropriate value.
! ______________________________________________________________________
! Quick Start Instructions

! (1) Compile this file. Then compile, link, and execute the program
!     example1.f90. The output is written to the file example1.dat.
!     Verify that the last line of the output is the string
!     'No errors occurred.'
! (2) Repeat this process for the program example2.f90.
!
! Other test programs you may wish to run to verify your installation
! of DVODE_F90 are:
!
! (3) Run the test programs nonstiffoptions.f90 and stiffoptions.f90
!     and verify that the last line in the output files produced is
!     'No errors occurred.' They solve the problems in the Toronto
!     test suites using several different error tolerances and various
!     solution options. Note that stiffoptions.f90 takes several
!     minutes to run because it performs several thousand separate
!     integrations.
! (4) Locate the file robertson.f90 in the demo programs and look at
!     how options are set using SET_OPTS, how DVODE_F90 is called to
!     obtain the solution at desired output times, and how the
!     derivative and Jacobian routines are supplied. Note too the
!     manner in which the solution is constrained to be nonnegative.
! (5) Locate demoharmonic.f90 and look at how root finding options
!     are set and how the event residual routine is supplied to
!     DVODE_F90.
! (6) The other demo programs available from the DVODE_F90 support
!     page illustrate various other solution options available in
!     DVODE_F90. The demo programs may be obtained from
!
!        http://www.radford.edu/~thompson/vodef90web/index.html
! ______________________________________________________________________
! DVODE_F90 Full Documentation Prologue

! Section 1.  Setting Options in DVODE_F90
! Section 2.  Calling DVODE_F90
! Section 3.  Choosing Error Tolerances
! Section 4.  Choosing the Method of Integration
! Section 5.  Interpolation of the Solution and Derivatives
! Section 6.  Handling Events (Root Finding)
! Section 7.  Gathering Integration Statistics
! Section 8.  Determining Jacobian Sparsity Structure Arrays
! Section 9.  Original DVODE Documentation Prologue
! Section 10. Example Usage

! Note: Search on the string 'Section' to locate these sections. You
! may wish to refer to the support page which has the sections broken
! into smaller pieces.
! ______________________________________________________________________
! Section 1.  Setting Options in DVODE_F90
!
! You can use any of three options routines:
!
! SET_NORMAL_OPTS
! SET_INTERMEDIATE_OPTS
! SET_OPTS

! OPTIONS = SET_NORMAL_OPTS(DENSE_J, BANDED_J, SPARSE_J,               &
!   USER_SUPPLIED_JACOBIAN, LOWER_BANDWIDTH, UPPER_BANDWIDTH,          &
!   RELERR, ABSERR, ABSERR_VECTOR, NEVENTS)

! OPTIONS = SET_INTERMEDIATE_OPTS(DENSE_J, BANDED_J, SPARSE_J,         &
!   USER_SUPPLIED_JACOBIAN,LOWER_BANDWIDTH, UPPER_BANDWIDTH,           &
!   RELERR, ABSERR, ABSERR_VECTOR,TCRIT, H0, HMAX, HMIN, MAXORD,       &
!   MXSTEP, MXHNIL, NZSWAG, USER_SUPPLIED_SPARSITY, MA28_RPS,          &
!   NEVENTS, CONSTRAINED, CLOWER, CUPPER, CHANGE_ONLY_f77_OPTIONS)     &

! OPTIONS = SET_OPTS(METHOD_FLAG, DENSE_J, BANDED_J, SPARSE_J,         &
!   USER_SUPPLIED_JACOBIAN, SAVE_JACOBIAN, CONSTANT_JACOBIAN,          &
!   LOWER_BANDWIDTH, UPPER_BANDWIDTH, SUB_DIAGONALS, SUP_DIAGONALS,    &
!   RELERR, RELERR_VECTOR, ABSERR, ABSERR_VECTOR, TCRIT, H0, HMAX,     &
!   HMIN, MAXORD, MXSTEP, MXHNIL, YMAGWARN, SETH, UPIVOT, NZSWAG,      &
!   USER_SUPPLIED_SPARSITY, NEVENTS, CONSTRAINED, CLOWER, CUPPER,      &
!   MA28_ELBOW_ROOM, MC19_SCALING, MA28_MESSAGES, MA28_EPS,            &
!   MA28_RPS, CHANGE_ONLY_f77_OPTIONS, JACOBIAN_BY_JACSP)

! Please refer to the documentation prologue for each of these functions
! to see what options may be used with each. Note that input to each is
! via keyword and all variables except error tolerances are optional.
! Defaults are used for unspecified options. If an option is available
! in SET_NORMAL OPTS, it is available and has the same meaning in
! SET_INTERMEDIATE_OPTS and SET_OPTS. Similarly, if an option is available
! in SET_INTERMEDIATE_OPTS, it is available and has the same meaning in
! SET_OPTS.

! The first two functions are provided merely for convenience.
! SET_NORMAL_OPTS is available simply to relieve you of reading the
! documentation for SET_OPTS and to use default values for all but
! the most common options. SET_INTERMEDIATE_OPTS is available to allow
! you more control of the integration while still using default values
! for less commonly used options. SET_OPTS allows you to specify any
! of the options available in DVODE_F90.

! Roughly, SET_NORMAL_OPTS is intended to provide for dense, banded,
! and numerical sparse Jacobians without the need to specify other
! specialized options. SET_INTERMEDIATE_OPTIONS is intended to allow
! more general sparse Jacobian options. SET_OPTS is intended to provide
! access to all options in DVODE_F90.

! Please note that SET_INTERMEDIATE_OPTS can be invoked using the same
! arguments as SET_NORMAL_OPTS; and SET_OPTS can be invoked using the
! same arguments as either SET_NORMAL_OPTS or SET_INTERMEDIATE_OPTS.
! If you wish you can simply delete SET_NORMAL_OPTS as well as
! SET_INTERMEDIATE_OPTS and use only SET_OPTS for all problems. If you
! do so, you need only include the options that you wish to use when
! you invoke SET_OPTIONS.

! In the following description any reference to SET_OPTS applies equally
! to SET_NORMAL_OPTS and SET_INTERMEDIATE OPTS.

! Before calling DVODE_F90 for the first time, SET_OPTS must be invoked.
! Typically, SET_OPTS is called once to set the desired integration
! options and parameters. DVODE_F90 is then called in an output loop to
! obtain the solution for the desired times. A detailed description of
! the DVODE_F90 arguments is given in a section below. Detailed descriptions
! of the options available via SET_OPTS are given in the documentation prologue.
! Although each option available in the f77 version of DVODE as well as
! several additional ones are available in DVODE_F90 via SET_OPTS,
! several of the available options are not relevant for most problems
! and need not be specified. Refer to the accompanying demonstration
! programs for specific examples of each usage. Note that after any call
! to DVODE_F90, you may call GET_STATS to gather relevant integration
! statistics. After your problem has completed, you may call
! RELEASE_ARRAYS to deallocate any internal arrays allocated by
! DVODE_F90 and to determine how much storage was used by DVODE_F90.
!
! To communicate with DVODE_F90 you will need to include the following
! statement in your calling program:
!    USE DVODE_F90_M
! and include the following statement in your type declarations section:
!    TYPE(VODE_OPTS) :: OPTIONS
! Below are brief summaries of typical uses of SET_OPTS.
! Nonstiff Problems:
! OPTIONS = SET_OPTS(RELERR=RTOL,ABSERR=ATOL)
!    The above use of SET_OPTS will integrate your system of odes
!    using the nonstiff Adams methods while using a relative error
!    tolerance of RTOL and an absolute error tolerance of ATOL.
!    Your subsequent call to DVODE_F90 might look like:
!    CALL DVODE_F90(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS)
! OPTIONS = SET_OPTS(RELERR=RTOL,ABSERR=ATOL,NEVENTS=NG)
!    If you wish to do root finding, SET_OPTS can be used as above.
!    Here, NEVENTS is the desired number of root finding functions.
!    Your subsequent call to DVODE_F90 might look like:
!    CALL DVODE_F90(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,G_FCN=G)
!    Here F is the name of your derivative subroutine and G is the
!    name of your subroutine to evaluate residuals of the root
!    finding functions.
! OPTIONS = SET_OPTS(RELERR=RTOL,ABSERR_VECTOR=ATOL)
!    This use of SET_OPTS indicates that a scalar relative error
!    tolerance and a vector of absolute error tolerances will be
!    used.
! Stiff Problems, internally generated dense Jacobian:
! OPTIONS = SET_OPTS(DENSE_J=.TRUE.,RELERR=RTOL,ABSERR=ATOL)
!    This use of DENSE_J=.TRUE. indicates that DVODE_F90 will
!    use the stiffly stable BDF methods and will approximate
!    the Jacobian, considered to be a dense matrix, using
!    finite differences. Your subsequent call to DVODE_F90
!    might look like:
!    CALL DVODE_F90(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS)
! OPTIONS = SET_OPTS(DENSE_J=.TRUE.,ABSERR=ATOL,RELERR=RTOL, &
!                    USER_SUPPLIED_JACOBIAN=.TRUE.)
!    If you know the Jacobian and wish to supply subroutine JAC
!    as described in the documentation for DVODE_F90, the options
!    call could look like the above.
!    Your subsequent call to DVODE_F90 might look like:
!    CALL DVODE_F90(F1,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JAC)
!    Here, JAC is the name of the subroutine that you provide to
!    evaluate the known Jacobian.
! Stiff Problems, internally generated banded Jacobian:
! OPTIONS = SET_OPTS(BANDED_J=.TRUE.,RELERR=RTOL,ABSERR=ATOL, &
!                       LOWER_BANDWIDTH=ML,UPPER_BANDWIDTH=MU)
!    This use of BANDED_J=.TRUE. indicates that DVODE_F90 will
!    use the stiffly stable BDF methods and will approximate the
!    Jacobian, considered to be a banded matrix, using finite
!    differences. Here ML is the lower bandwidth of the Jacobian
!    and ML is the upper bandwidth of the Jacobian.
! Stiff Problems, internally generated sparse Jacobian:
! OPTIONS = SET_OPTS(SPARSE_J=.TRUE.,ABSERR=ATOL,RELERR=RTOL)
!    This use of SET_OPTS indicates that the Jacobian is a sparse
!    matrix. Its structure will be approximated internally by
!    making calls to your derivative routine. If you know the
!    structure before hand, you may provide it directly in a
!    variety of ways as described in the documentation prologue
!    for SET_OPTS. In addition, several other options related
!    to sparsity are available.
! More complicated common usage:
!    Suppose you need to solve a stiff problem with a sparse Jacobian.
!    After some time, the structure of the Jacobian changes and you
!    wish to have DVODE_F90 recalculate the structure before continuing
!    the integration. Suppose that initially you want to use an absolute
!    error tolerance of 1.0D-5 and that when the Jacobian structure is
!    changed you wish to reduce the error tolerance 1.0D-7. Your calls
!    might look like this.
!    RTOL = ...
!    ATOL = 1.0D-5
!    OPTIONS = SET_OPTS(SPARSE_J=.TRUE.,ABSERR=ATOL,RELERR=RTOL)
!    Output loop:
!       CALL DVODE_F90(FCN,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS)
!       At desired time:
!       ISTATE = 3
!       ATOL = 1.0D-7
!       OPTIONS = SET_OPTS(SPARSE_J=.TRUE.,ABSERR=ATOL,RELERR=RTOL)
!    End of output loop

! In the following we have summarized and described how some of the demonstration
! programs set options and call DVODE_F90. In each case the necessary parameters
! are defined before invoking SET_OPTS. The call to DVODE_F90 is in a loop in
! which the output time is successively updated. The actual programs are available
! from the web support page
!
!    http://www.radford. edu/~thompson/vodef90web/index.html/
!
!                              Problem Summary
!
! Problem                  NEQ      Jacobian            Features Illustrated
!
! Prologue Example 1        3        Dense               Basic
!
! Prologue Example 2        3        Dense               Root finding
!
! Robertson                 3        Dense               Solution bounds
!
! Harmonic Oscillator       4        Nonstiff            Root finding
!
! Flow Equations         5-1800      Sparse              Automatic determination
!                                                        of sparsity arrays
!
! Diurnal Kinetics      50-5000      Sparse or banded    Sparsity options
!
!                       Options Used and DVODE_F90 Call
!
! Prologue Example 1
!
!    OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE., ABSERR_VECTOR=ATOL, RELERR=RTOL,     &
!              USER_SUPPLIED_JACOBIAN=.TRUE.)
!    CALL DVODE_F90(FEX,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JEX)
!
!    The problem consists of a stiff system of NEQ=3 equations. The dense
!    Jacobian option (DENSE_J) is used. A vector ATOL(*) of error tolerances
!    is used. A scalar relative error tolerance RTOL is used. Subroutine JEX
!    is provided to evaluate the analytical Jacobian. If the last argument
!    J_FCN=JEX is omitted (as in Example 2), a numerical Jacobian will
!    be used.
!
! Prologue Example 2
!
!    OPTIONS = SET_NORMAL_OPTS(DENSE_J=.TRUE., RELERR=RTOL, ABSERR_VECTOR=ATOL,     &
!              NEVENTS=NG)
!    CALL DVODE_F90(FEX,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,G_FCN=GEX)
!
!    The system in Example 1 is used to illustrate root finding. It is
!    desired to locate the times at which two of the solution components
!    attain prescribed values. NEVENTS=2 informs the solver that two such
!    functions are used. Subroutine GEX is used to calculate the residuals
!    for these two functions. A dense numerical Jacobian is used.
!
! Robertson Problem
!
!    OPTIONS = SET_INTERMEDIATE_OPTS(DENSE_J=.TRUE., RELERR_VECTOR=RTOL,            &
!              ABSERR_VECTOR=ABSERR_TOLERANCES, CONSTRAINED=BOUNDED_COMPONENTS,     &
!              CLOWER=LOWER_BOUNDS, CUPPER=UPPER_BOUNDS)
!
!    CALL DVODE_F90(DERIVS,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JACD)
!    The system used in Examples 1 and 2 is solved over a much larger
!    time interval. The solution is constrained to be nonegative. This
!    is done by prescribing the components to be constrained (BOUNDED_COMPONENTS).
!    Artificially large values are used to impose upper bounds (UPPER_BOUNDS)
!    and lower bounds of zero are used to force a nonnegative solution.
!
! Harmonic Oscillator Problem
!
!    OPTIONS = SET_NORMAL_OPTS(RELERR=RTOL, ABSERR=ATOL, NEVENTS=NG)
!    CALL DVODE_F90(DERIVS,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,G_FCN=GEVENTS)
!
!    A nonstiff system of NEQ=4 equations is solved. The nonstiff option is
!    used because neither DENSE_ nor BANDED_J nor SPARSE_J is present. It is
!    desired to find the times at which Y(2) or Y(3) is equal to 0. Residuals
!    for the two corresponding event functions are calculated in subroutine
!    GEVENTS.
!
! Flow Equations Problem
!
!    OPTIONS = SET_OPTS(SPARSE_J=SPARSE, ABSERR=ATOL(1), RELERR=RTOL(1),            &
!              MXSTEP=100000, NZSWAG=20000)
!    CALL DVODE_F90(FCN,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS)
!
!    This is a stiff system of equations resulting a method of lines
!    discretization. The Jacobian is sparse. Scalar absolute and relative
!    error tolerances are used. The Jacobian structure and a numerical
!    Jacobian are used. The solver is limited to a maximum of MXSTEP steps.
!    NZSWAG is the amount by which allocated array sizes will be increased.
!    The accompanying test program may be used to illutrate several other
!    solution options.
!
! Diurnal Kinetics Problem
!
!    OPTIONS = SET_OPTS(SPARSE_J=SPARSE, BANDED_J=BANDED, DENSE_J=DENSE,            &
!              ABSERR_VECTOR=ATOL(1:NEQ), RELERR=RTOL(1), MXSTEP=100000,            &
!              NZSWAG=50000, HMAX=MAXH, LOWER_BANDWIDTH=ML, UPPER_BANDWIDTH=MU,     &
!              MA28_ELBOW_ROOM=10, MC19_SCALING=.TRUE., MA28_MESSAGES=.FALSE.,      &
!              MA28_EPS=1.0D-4, MA28_RPS=.TRUE.,                                    &
!              USER_SUPPLIED_SPARSITY=SUPPLY_STRUCTURE)
!   CALL USERSETS_IAJA(IA, IADIM, JA, JADIM)
!   CALL DVODE_F90(FCN, NEQ, Y, T, TOUT, ITASK, ISTATE, OPTIONS)
!
!   This problem can be used to illustrate most solution options. Here, dense,
!   banded, or sparse Jacobians are used depending on the values of the first
!   three parameters. A vector error tolerance is used and a scalar relative
!   error tolerance is used. If a banded solution is desired, it is necessary
!   to supply the bandwidths ML and MU. If a sparse solution is desired,
!   several special options are used. The most important one is MA28_RPS to
!   force the solver to update the partial pivoting sequence when singular
!   iteration matrices are encountered. The sparsity pattern is determined
!   numerically if SUPPLY_STRUCTURE is FALSE. Otherwise the user will supply
!   the pattern by calling subroutine USERSETS_IAJA.
! ______________________________________________________________________
! Section 2.  Calling DVODE_F90
!
! DVODE_F90 solves the initial value problem for stiff or nonstiff
! systems of first order ODEs,
!     dy/dt = f(t,y), or, in component form,
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
! DVODE_F90 is a package based on the EPISODE and EPISODEB packages,
! and on the ODEPACK user interface standard. It was developed from
! the f77 solver DVODE developed by Brown, Byrne, and Hindmarsh.
! DVODE_F90 also provides for the solution of sparse systems in a
! fashion similar to LSODES and LSOD28. Currently, MA28 is used
! to perform the necessary sparse linear algebra. DVODE_F90 also
! contains the provision to do root finding in a fashion similar
! to the LSODAR solver from ODEPACK.

! Communication between the user and the DVODE_F90 package, for normal
! situations, is summarized here. This summary describes only a subset
! of the full set of options available. See the full description for
! details, including optional communication, nonstandard options, and
! instructions for special situations.
!    CALL DVODE_F90(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JAC,G_FCN=GEX)
!    The arguments in the call list to DVODE_F90 have the following
!    meanings.
! F       = The name of the user-supplied subroutine defining the
!           ODE system. The system must be put in the first-order
!           form dy/dt = f(t,y), where f is a vector-valued function
!           of the scalar t and the vector y. Subroutine F is to
!           compute the function f. It is to have the form
!                SUBROUTINE F(NEQ,T,Y,YDOT)
!                DOUBLE PRECISION T,Y(NEQ),YDOT(NEQ)
!           where NEQ, T, and Y are input, and the array YDOT = f(t,y)
!           is output. Y and YDOT are arrays of length NEQ.
!           Subroutine F should not alter Y(1),...,Y(NEQ).
!           If F (and JAC) are not contained in a module available to
!           your calling program, you must declare F to be EXTERNAL
!           in the calling program.
! NEQ     = The size of the ODE system (number of first order
!           ordinary differential equations).
! Y       = A double precision array for the vector of dependent variables,
!           of length NEQ or more. Used for both input and output on the
!           first call (ISTATE = 1), and only for output on other calls.
!           On the first call, Y must contain the vector of initial
!           values. In the output, Y contains the computed solution
!           evaluated at T.
! T       = The independent variable. In the input, T is used only on
!           the first call, as the initial point of the integration.
!           In the output, after each call, T is the value at which a
!           computed solution Y is evaluated (usually the same as TOUT).
!           On an error return, T is the farthest point reached.
! TOUT    = The next value of t at which a computed solution is desired.
!           TOUT is Used only for input. When starting the problem
!           (ISTATE = 1), TOUT may be equal to T for one call, then
!           should not equal T for the next call. For the initial T,
!           an input value of TOUT unequal to T is used in order to
!           determine the direction of the integration (i.e. the
!           algebraic sign of the step sizes) and the rough scale
!           of the problem. Integration in either direction (forward
!           or backward in t) is permitted. If ITASK = 2 or 5 (one-step
!           modes), TOUT is ignored after the first call (i.e. the
!           first call with TOUT \= T). Otherwise, TOUT is required
!           on every call. If ITASK = 1, 3, or 4, the values of TOUT
!           need not be monotone, but a value of TOUT which backs up
!           is limited to the current internal t interval, whose
!           endpoints are TCUR - HU and TCUR. (Refer to the description
!           of GET_STATS for a description of TCUR and HU.)
! ITASK   = An index specifying the task to be performed.
!           Input only. ITASK has the following values and meanings.
!           1  means normal computation of output values of y(t) at
!              t = TOUT (by overshooting and interpolating).
!           2  means take one step only and return.
!           3  means stop at the first internal mesh point at or
!              beyond t = TOUT and return.
!           4  means normal computation of output values of y(t) at
!              t = TOUT but without overshooting t = TCRIT.
!              TCRIT must be specified in your SET_OPTS call. TCRIT
!              may be equal to or beyond TOUT, but not behind it in
!              the direction of integration. This option is useful
!              if the problem has a singularity at or beyond t = TCRIT.
!           5  means take one step, without passing TCRIT, and return.
!              TCRIT must be specified in your SET_OPTS call.
!           If ITASK = 4 or 5 and the solver reaches TCRIT (within
!           roundoff), it will return T = TCRIT(exactly) to indicate
!           this (unless ITASK = 4 and TOUT comes before TCRIT, in
!           which case answers at T = TOUT are returned first).
! ISTATE  = an index used for input and output to specify the
!           the state of the calculation.
!           In the input, the values of ISTATE are as follows.
!           1  means this is the first call for the problem
!              (initializations will be done). See note below.
!           2  means this is not the first call, and the calculation
!              is to continue normally, with no change in any input
!              parameters except possibly TOUT and ITASK.
!           3  means this is not the first call, and the
!              calculation is to continue normally, but with
!              a change in input parameters other than
!              TOUT and ITASK. Desired changes require SET_OPTS
!              be called prior to calling DVODE_F90 again.
!           A preliminary call with TOUT = T is not counted as a
!           first call here, as no initialization or checking of
!           input is done. (Such a call is sometimes useful to
!           include the initial conditions in the output.)
!           Thus the first call for which TOUT is unequal to T
!           requires ISTATE = 1 in the input.
!           In the output, ISTATE has the following values and meanings.
!            1  means nothing was done, as TOUT was equal to T with
!               ISTATE = 1 in the input.
!            2  means the integration was performed successfully.
!            3  means a root of one of your root finding functions
!               has been located.
!           A negative value of ISTATE indicates that DVODE_F90
!           encountered an error as described in the printed error
!           message. Since the normal output value of ISTATE is 2,
!           it does not need to be reset for normal continuation.
!           Also, since a negative input value of ISTATE will be
!           regarded as illegal, a negative output value requires
!           the user to change it, and possibly other input, before
!           calling the solver again.
! OPTIONS = The options structure produced by your call to SET_OPTS.
! JAC     = The name of the user-supplied routine (MITER = 1 or 4 or 6)
!           If you do not specify that a stiff method is to be used
!           in your call to SET_OPTS, you need not include JAC in
!           your call to DVODE_F90. If you specify a stiff method and
!           that a user supplied Jacobian will be supplied, JAC must
!           compute the Jacobian matrix, df/dy, as a function of the
!           scalar t and the vector y. It is to have the form:
!              SUBROUTINE JAC(NEQ, T, Y, ML, MU, PD, NROWPD)
!              DOUBLE PRECISION T, Y(NEQ), PD(NROWPD,NEQ)
!           where NEQ, T, Y, ML, MU, and NROWPD are input and the array
!           PD is to be loaded with partial derivatives (elements of the
!           Jacobian matrix) in the output. PD must be given a first
!           dimension of NROWPD. T and Y have the same meaning as in
!           Subroutine F.
!           In the full matrix case (MITER = 1), ML and MU are
!           ignored, and the Jacobian is to be loaded into PD in
!           columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
!           In the band matrix case (MITER = 4), the elements
!           within the band are to be loaded into PD in columnwise
!           manner, with diagonal lines of df/dy loaded into the rows
!           of PD. Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
!           ML and MU are the half-bandwidth parameters. (See IUSER).
!           The locations in PD in the two triangular areas which
!           correspond to nonexistent matrix elements can be ignored
!           or loaded arbitrarily, as they are overwritten by DVODE_F90.
!           In the sparse matrix case the elements of the matrix
!           are determined by the sparsity structure given by the
!           IA and JA pointer arrays. Refer to the documentation
!           prologue for SET_OPTS for a description of the arguments
!           for JAC since they differ from the dense and banded cases.
!           JAC need not provide df/dy exactly. A crude
!           approximation (possibly with a smaller bandwidth) will do.
!           In either case, PD is preset to zero by the solver,
!           so that only the nonzero elements need be loaded by JAC.
!           In the sparse matrix case, JAC has a different form:
!                SUBROUTINE JAC (N, T, Y, IA, JA, NZ, PD)
!           Given the number of odes N, the current time T, and the
!           current solution vector Y, JAC must do the following:
!              If NZ = 0 on input:
!              Replace NZ by the number of nonzero elements in the
!              Jacobian. The diagonal of the Jacobian must be included.
!              Do NOT define the arrays IA, JA, PD at this time.
!              Once JAC has been called with NZ = 0 and you have
!              defined the value of NZ, future calls to JAC will use
!              this value of NZ.
!              When a call is made with NZ unequal to 0, you must
!              define the sparsity structure arrays IA and JA, and
!              the sparse Jacobian PD.
!                 IA defines the number of nonzeros including the
!                 diagonal in each column of the Jacobian. Define
!                 IA(1) = 1 and for J = 1,..., N,
!                 IA(J+1) = IA(J) + number of nonzeros in column J.
!                 Diagonal elements must be included even if they are
!                 zero. You should check to ensure that IA(N+1)-1 = NZ.
!                 JA defines the rows in which the nonzeros occur.
!                 For I = 1,...,NZ, JA(I) is the row in which the Ith
!                 element of the Jacobian occurs. JA must also include
!                 the diagonal elements of the Jacobian.
!                 PD defines the numerical value of the Jacobian
!                 elements. For I = 1,...,NZ, PD(I) is the numerical
!                 value of the Ith element in the Jacobian. PD must
!                 also include the diagonal elements of the Jacobian.
! GFUN    = the name of the subroutine to evaluate the residuals for
!           event functions. If you do not specify that events are
!           present (by specifying NEVENTS > 0 in SET_OPTS), you
!           need not include GFUN in your call list for DVODE_F90.
!           If GFUN is not contained in a module available to your
!           calling program, you must declare GFUN to be EXTERNAL
!           in your calling program.
! To continue the integration after a successful return, simply
! reset TOUT and call DVODE_F90 again. No other parameters need
! be reset unless ISTATE=3 in which case, reset it to 2 before
! calling DVODE_F90 again.
! ______________________________________________________________________
! Section 3.  Choosing Error Tolerances
!
! This is the most important aspect of solving odes numerically.
! You may supply any of four keywords and values. If you wish to
! use scalar error tolerances, you should supply ABSERR and RELERR.
! For a good many problems, it is advisable to supply a vector of
! absolute error tolerances ABSERR_VECTOR = desired vector. This
! allows you to use different tolerances for each component of
! the solution. If ABSERR_VECTOR is supplied, it must be a vector
! of length NEQ where NEQ is the number of odes in your system.
! Similarly, you may supply a vector of relative error tolerances,
! RELERR_VECTOR. If no tolerances are specified, DVODE_F90 will use
! default error tolerances ABSERR=1D-6 and RELERR=1D-4; but it is
! strongly recommended that you supply values that are appropriate
! for your problem. In the event you do not supply error tolerances,
! DVODE_F90 will print a reminder that the default error tolerances
! are not appropriate for all problems.
!
! RELERR can be set to a scalar value as follows.
! Determine the number of significant digits of accuracy desired, 
! which will be a positive integer, say, N.
! Then set RELERR = 10**-(N+1). 
! The authors recommend that RELERR be no larger than 10**-4.
! The authors recommend a vector valued absolute error tolerance,
! which can be set as follows.
! For the I-th component of the solution vector, Y(I), determine
! the positive number FLOOR(I) at which ABS(Y(I)) becomes
! negligible for the problem at hand. FLOOR(I) is sometimes called
! the problem zero or the floor value for the I-th component and is
! problem dependent. For a given problem that is not scaled, these
! floor values may well vary by up to 9 orders of magnitude.
! Set ABSERR(I) = FLOOR(I) or to be conservative
! ABSERR_VECTOR(I) = 0.1*FLOOR(I). There is no variable FLOOR in
! DVODE_F90. If it is difficult to divine the components of ABSERR,
! (or FLOOR) make a reasonable guess, run the problem, then set that
! ABSERR_VECTOR so for I = 1, 2,...NEQ,
! ABSERR_VECTOR(I) = 1D-6*RELERR*MAX{ABS(Y(I,TOUT): for all TOUT}.
! The correct choices for RELERR and ABSERR can and do have
! significant impact on both the quality of the solution and run
! time. Counter intuitively, error tolerances that are too loose
! can and do increase run time significantly and the quality of
! the solution may be seriously compromised.
! Examples:
! 1. OPTIONS = SET_OPTS(DENSE_J=.TRUE., ABSERR=1D-8,RELERR=1D-8)
!    This will yield MF = 22. Both the relative error tolerance
!    and the absolute error tolerance will equal 1D-8.
! 2. OPTIONS = SET_OPTS(DENSE_J=.TRUE.,RELERR=1D-5, &
!                       ABSERR_VECTOR=(/1D-6,1D-8/))
!    For a system with NEQ=2 odes, this will yield MF = 22. A scalar
!    relative error tolerance equal to 1D-5 will be used. Component 1
!    of the solution will use an absolute error tolerance of 1D-6
!    while component 2 will use an absolute error tolerance of 1D-8.
! ______________________________________________________________________
! Section 4.  Choosing the Method of Integration
!
! If you wish to specify a numerical value for METHOD_FLAG, it can equal
! any of the legal values of MF for DVODE or LSODES. (See below.) If
! you do not wish to specify a numerical value for METHOD_FLAG, you
! may specify any combination of the five logical keywords DENSE_J,
! BANDED_J, SPARSE_J, USER_SUPPLIED_JACOBIAN, SAVE_JACOBIAN that you
! wish. Appropriate values will be used in DVODE_F90 for any variables
! that are not present. The first three flags indicate the type of
! Jacobian, dense, banded, or sparse. If USER_SUPPLIED_JACOBIAN=.TRUE.,
! the Jacobian will be determined using the subroutine JAC you supply
! in your call to DVODE_F90. Otherwise, an internal Jacobian will be
! generated using finite differences.
! Examples:
! 1. OPTIONS = SET_OPTS(METHOD_FLAG=22,...)
!    DVODE will use the value MF=22 as described in the documentation
!    prologue. In this case, the stiff BDF methods will be used with
!    a dense, internally generated Jacobian.
! 2. OPTIONS = SET_OPTS(METHOD_FLAG=21,...)
!    This is the same an Example 1 except that a user supplied dense
!    Jacobian will be used and DVODE will use MF=21.
! 3. OPTIONS = SET_OPTS(DENSE_J=.TRUE.,...)
!    This will yield MF = 22 as in Example 1, provided
!    USER_SUPPLIED_JACOBIAN and SAVE_JACOBIAN are not present, or if
!    present are set to their default
!     values of .FALSE. and .TRUE., respectively.
! 4. OPTIONS = SET_OPTS(DENSE_J=.TRUE.,&
!                       USER_SUPPLIED_JACOBIAN=.TRUE.,...)
!    This will yield MF = 21 as in Example 2, provided SAVE_JACOBIAN
!    is not present, or if present is set to its default value .FALSE.
! Notes:
! 1. If you specify more than one of DENSE_J, BANDED_J, and SPARSE_J,
!    DENSE_J takes precedence over BANDED_J which in turn takes
!    precedence over SPARSE_J.
! 2. By default, DVODE_F90 saves a copy of the Jacobian and reuses the
!    copy when necessary. For problems in which storage requirements
!    are acute, you may wish to override this default and have
!    DVODE_F90 recalculate the Jacobian rather than use a saved copy.
!    You can do this by specifying SAVE_JACOBIAN=.FALSE. It is
!    recommended that you not do this unless necessary since it can
!    have a significant impact on the efficiency of DVODE_F90. (For
!    example, when solving a linear problem only one evaluation of
!    the Jacobian is required with the default option.)
! 3. If you choose BANDED_J = .TRUE. or if you supply a value of MF
!    that corresponds to a banded Jacobian, you must also supply the
!    lower  bandwidth ML and the upper bandwidth of the Jacobian MU
!    by including the keywords
!    LOWER_BANDWIDTH = value of ML and UPPER_BANDWIDTH = value of M
!                   More on Method Selection
! The keyword options available in SET_OPTS are intended to replace
! the original method indicator flag MF. However, if you wish to
! retain the flexibility of the original solver, you may specify MF
! directly in your call to SET_OPTS. This is done by using the
! keyword METHOD_FLAG=MF in your SET_OPTS where MF has any of the
! values in the following description. Refer to the demonstration
! program demosp.f90 for an example in which this is done.
! MF     = The method flag. Used only for input. The legal values of
!          MF are:
!          10, 11, 12, 13, 14, 15, 16, 17, 20, 21, 22, 23, 24, 25, 26,
!          27, -11, -12, -14, -15, -21, -22, -24, -25, -26, -27.
!          MF is a signed two-digit integer, MF = JSV*(10*METH+MITER).
!          JSV = SIGN(MF) indicates the Jacobian-saving strategy:
!            JSV =  1 means a copy of the Jacobian is saved for reuse
!                     in the corrector iteration algorithm.
!            JSV = -1 means a copy of the Jacobian is not saved
!                     (valid only for MITER = 1, 2, 4, or 5).
!          METH indicates the basic linear multistep method:
!            METH = 1 means the implicit Adams method.
!            METH = 2 means the method based on backward
!                     differentiation formulas (BDF-s).
!          MITER indicates the corrector iteration method:
!            MITER = 0 means functional iteration (no Jacobian matrix
!                      is involved).
!            MITER = 1 means chord iteration with a user-supplied
!                      full (NEQ by NEQ) Jacobian.
!            MITER = 2 means chord iteration with an internally
!                      generated (difference quotient) full Jacobian
!                      (using NEQ extra calls to F per df/dy value).
!            MITER = 4 means chord iteration with a user-supplied
!                      banded Jacobian.
!            MITER = 5 means chord iteration with an internally
!                      generated banded Jacobian (using ML+MU+1 extra
!                      calls to F per df/dy evaluation).
!            MITER = 6 means chord iteration with a user-supplied
!                      sparse Jacobian.
!            MITER = 7 means chord iteration with an internally
!                      generated sparse Jacobian
!          If MITER = 1, 4, or 6 the user must supply a subroutine
!          JAC(the name is arbitrary) as described above under JAC.
!          For other values of MITER, JAC need not be provided.
! ______________________________________________________________________
! Section 5.  Interpolation of the Solution and Derivative
!
! Following a successful return from DVODE_F90, you may call
! subroutine DVINDY to interpolate the solution or derivative.
! SUBROUTINE DVINDY(T, K, DKY, IFLAG)
! DVINDY computes interpolated values of the K-th derivative of the
! dependent variable vector y, and stores it in DKY. This routine
! is called with K = 0 or K = 1 and T = TOUT. In either case, the
! results are returned in the array DKY of length at least NEQ which
! must be declared and dimensioned in your calling program. The
! computed values in DKY are obtained by interpolation using the
! Nordsieck history array.
! ______________________________________________________________________
! Section 6.  Handling Events (Root Finding)
!
!    DVODE_F90 contains root finding provisions. It attempts to
!    locates the roots of a set of functions
!         g(i) = g(i,t,y(1),...,y(NEQ))  (i = 1,...,ng).
!    To use root finding include NEVENTS=NG in your call to SET_OPTS
!    where NG is the number of root finding functions. You must then
!    supply subroutine GFUN in your call to DVODE_F90 using
!    G_FCN=GFUN as the last argument. GFUN must has the form
!               SUBROUTINE GFUN(NEQ, T, Y, NG, GOUT)
!    where NEQ, T, Y, and NG are input, and the array GOUT is output.
!    NEQ, T, and Y have the same meaning as in the F routine, and
!    GOUT is an array of length NG. For i = 1,...,NG, this routine is
!    to load into GOUT(i) the value at (T,Y) of the i-th constraint
!    function g(i). DVODE_F90 will find roots of the g(i) of odd
!    multiplicity (i.e. sign changes) as they occur during the
!    integration. GFUN must be declared EXTERNAL in the calling
!    program. Note that because of numerical errors in the functions
!    g(i) due to roundoff and integration error, DVODE_F90 may return
!    false roots, or return the same root at two or more nearly equal
!    values of t. This is particularly true for problems in which the
!    integration is restarted (ISTATE = 1) at a root. If such false
!    roots are suspected, you should consider smaller error tolerances
!    and/or higher precision in the evaluation of the g(i). Note
!    further that if a root of some g(i) defines the end of the
!    problem, the input to DVODE_F90 should nevertheless allow
!    integration to a point slightly past that root, so that DVODE_F90
!    can locate the root by interpolation. Each time DVODE_F90 locates
!    a root of one of your event functions it makes a return to the
!    calling program with ISTATE = 3. When such a return is made and
!    you have processed the results, simply change ISTATE = 2 and call
!    DVODE_F90 again without making other changes.
! ______________________________________________________________________
! Section 7.  Gathering Integration Statistics
!
! SUBROUTINE GET_STATS(RSTATS, ISTATS, NUMEVENTS, JROOTS)
! Caution:
! RSTATS and ISTATS must be declared and dimensioned in your
! main program. The minimum dimensions are:
! DIMENSION RSTATS(22), ISTATS(31)
! This subroutine returns the user portions of the original DVODE
! RUSER and IUSER arrays, and if root finding is being done, it
! returns the original LSODAR JROOT vector. NUMEVENTS and JROOTS
! are optional parameters. NUMEVENTS is the number of root functions
! and JROOTS is an integer array of length NUMEVENTS.
! Available Integration Statistics:
! HU      RUSER(11) The step size in t last used (successfully).
! HCUR    RUSER(12) The step size to be attempted on the next step.
! TCUR    RUSER(13) The current value of the independent variable
!                   which the solver has actually reached, i.e. the
!                   current internal mesh point in t. In the output,
!                   TCUR will always be at least as far from the
!                   initial value of t as the current argument T,
!                   but may be farther (if interpolation was done).
! TOLSF   RUSER(14) A tolerance scale factor, greater than 1.0,
!                   computed when a request for too much accuracy was
!                   detected (ISTATE = -3 if detected at the start of
!                   the problem, ISTATE = -2 otherwise). If ITOL is
!                   left unaltered but RTOL and ATOL are uniformly
!                   scaled up by a factor of TOLSF for the next call,
!                   then the solver is deemed likely to succeed.
!                   (The user may also ignore TOLSF and alter the
!                   tolerance parameters in any other way appropriate.)
! NST     IUSER(11) The number of steps taken for the problem so far.
! NFE     IUSER(12) The number of f evaluations for the problem so far.
! NJE     IUSER(13) The number of Jacobian evaluations so far.
! NQU     IUSER(14) The method order last used (successfully).
! NQCUR   IUSER(15) The order to be attempted on the next step.
! IMXER   IUSER(16) The index of the component of largest magnitude in
!                   the weighted local error vector (E(i)/EWT(i)),
!                   on an error return with ISTATE = -4 or -5.
! LENRW   IUSER(17) The length of RUSER actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
! LENIW   IUSER(18) The length of IUSER actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
! NLU     IUSER(19) The number of matrix LU decompositions so far.
! NNI     IUSER(20) The number of nonlinear (Newton) iterations so far.
! NCFN    IUSER(21) The number of convergence failures of the nonlinear
!                   solver so far.
! NETF    IUSER(22) The number of error test failures of the integrator
!                   so far.
! MA28AD_CALLS      IUSER(23) The number of calls made to MA28AD
! MA28BD_CALLS      IUSER(24) The number of calls made to MA28BD
! MA28CD_CALLS      IUSER(25) The number of calls made to MA28CD
! MC19AD_CALLS      IUSER(26) The number of calls made to MC19AD
! IRNCP             IUSER(27) The number of compressions done on array JAN
! ICNCP             IUSER(28) The number of compressions done on array ICN
! MINIRN            IUSER(29) Minimum size for JAN array
! MINICN            IUSER(30) Minimum size for ICN array
! MINNZ             IUSER(31) Number of nonzeros in sparse Jacobian
! JROOTS  JROOTS    Optional array of component indices for components
!                   having a zero at the current time
! ______________________________________________________________________
! Section 8.  Determining Jacobian Sparsity Structure Arrays
!
! If you are solving a problem with a sparse Jacobian, the arrays
! that define the sparsity structure are needed. The arrays may
! be determined in any of several ways.
! 1. If you choose the default mode by indicating SPARSE=.TRUE.,
!    the sparsity arrays will be determined internally by DVODE_F90
!    by making calls to your derivative subroutine. This mode is
!    equivalent to using the integration method flag MF = 227.
! 2. The DVODE_F90 method flag MF is defined to be
!    MF = 100*MOSS + 10*METH + MITER. If you supply MF = 227 (or 217),
!    the sparse Jacobian will be determined using finite differences;
!    and the sparsity arrays will be determined by calling your
!    derivative subroutine.
! 3. If you supply MF = 126 (or 116), you must supply the Jacobian
!    subroutine JAC to define the exact Jacobian. JAC must have the
!    following form:
!           SUBROUTINE JAC (N, T, Y, IA, JA, NZ, PD)
!    Given the number of odes N, the current time T, and the current
!    solution vector Y, JAC must do the following:
!    -  If NZ = 0 on input:
!       Replace NZ by the number of nonzero elements in the Jacobian.
!       The diagonal of the Jacobian must be included.
!       Do NOT define the arrays IA, JA, PD at this time.
!       Once JAC has been called with NZ = 0 and you have defined the
!       value of NZ, future calls to JAC will use this value of NZ.
!    -  When a call is made with NZ unequal to 0, you must define the
!       sparsity structure arrays IA and JA, and the sparse Jacobian
!       PD.
!         - IA defines the number of nonzeros including the diagonal
!           in each column of the Jacobian. Define IA(1) = 1 and for
!           J = 1,..., N,
!           IA(J+1) = IA(J) + number of nonzeros in column J.
!           Diagonal elements must be include even if they are zero.
!           You should check to ensure that IA(N+1)-1 = NZ.
!         - JA defines the rows in which the nonzeros occur. For
!           I = 1,...,NZ, JA(I) is the row in which the Ith element
!           of the Jacobian occurs. JA must also include the diagonal
!           elements of the Jacobian.
!         - PD defines the numerical value of the Jacobian elements.
!           For I = 1,...,NZ, PD(I) is the numerical value of the
!           Ith element in the Jacobian. PD must also include the
!           diagonal elements of the Jacobian.
! 4. If you wish to supply the IA and JA arrays directly, use
!    MF = 27. In this case, after calling SET_OPTS, you must call
!    SET_IAJA supplying the arrays IAUSER and JAUSER described in
!    the documentation prologue for SET_IAJA. These arrays will be
!    used when approximate Jacobians are determined using finite
!    differences.
! There are two user callable sparsity structure subroutines:
! USERSETS_IAJA may be used if you wish to supply the sparsity
! structure directly.
! SUBROUTINE USERSETS_IAJA(IAUSER,NIAUSER,JAUSER,NJAUSER)
!     Caution:
!     If it is called, USERSETS_IAJA must be called after the
!     call to SET_OPTS.
!     Usage:
!     CALL SET_IAJA(IAUSER,NIAUSER,JAUSER,NJAUSER)
!       In this case, IAUSER of length NIAUSER will be used for
!       IA; and JAUSER of length NJAUSER will be used for JA.
!     Arguments:
!     IAUSER  = user supplied IA array
!     NIAUSER = length of IAUSER array
!     JAUSER  = user supplied JA vector
!     NJAUSER = length of JAUSER array
! The second subroutine allows you to approximate the sparsity
! structure using derivative differences. It allows more flexibility
! in the determination of perturbation increments used.
! SUBROUTINE SET_IAJA(DFN,NEQ,T,Y,FMIN,NTURB,DTURB,IAUSER, &
!   NIAUSER, JAUSER, NJAUSER)
!     Caution:
!     If it is called, SET_IAJA must be called after the call to
!     SET_OPTS.
!     Usage:
!     SET_IAJA may be called in one of two ways:

!     CALL SET_IAJA(DFN,NEQ,T,Y,FMIN,NTURB,DTURB)
!       In this case IA and JA will be determined using calls
!       to your derivative routine DFN.
!     CALL SET_IAJA(DFN,NEQ,T,Y,FMIN,NTURB,DTURB,IAUSER,NIAUSER, &
!       JAUSER, NJAUSER)
!       In this case, IAUSER of length NIAUSER will be used for
!       IA; and JAUSER of length NJAUSER will be used for JA.
!       T, Y, FMIN, NTURB, and DTURB will be ignored (though
!       they must be present in the argument list).
!     Arguments:
!     DFN     = DVODE derivative subroutine
!     NEQ     = Number of odes
!     T       = independent variable t
!     Y       = solution y(t)
!     FMIN    = Jacobian threshold value. Elements of the Jacobian
!               with magnitude smaller than FMIN will be ignored.
!               FMIN will be ignored if it is less than or equal
!               to ZERO.
!     NTURB   = Perturbation flag. If NTURB=1, component I of Y
!               will be perturbed by 1.01D0.
!               If NTURB=NEQ, component I of Y will be perturbed
!               by ONE + DTURB(I).
!     DTURB   = perturbation vector of length 1 or NEQ.
!     If these four optional parameters are present, IAUSER and JAUSER
!     will be copied to IA and JA rather than making derivative calls
!     to approximate IA and JA:
!        IAUSER  = user supplied IA array
!        NIAUSER = length of IAUSER array
!        JAUSER  = user supplied JA vector
!        NJAUSER = length of JAUSER array
! ______________________________________________________________________
! Section 9.  Original DVODE.F Documentation Prologue
!
! SUBROUTINE DVODE(F, NEQ, Y, T, TOUT, ITASK, ISTATE, OPTS, JAC, GFUN)
! DVODE: Variable-coefficient Ordinary Differential Equation solver,
! with fixed-leading-coefficient implementation.
! Note:
! Numerous changes have been made in the documentation and the code
! from the original Fortran 77 DVODE solver. With regard to the new
! F90 version, if you choose options that correspond to options
! available in the original f77 solver, you should obtain the same
! results. In all testing, identical results have been obtained
! between this version and a simple F90 translation of the original
! solver.
! DVODE solves the initial value problem for stiff or nonstiff
! systems of first order ODEs,
!     dy/dt = f(t,y), or, in component form,
!     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
! DVODE is a package based on the EPISODE and EPISODEB packages, and
! on the ODEPACK user interface standard, with minor modifications.
! This version is based also on elements of LSODES and LSODAR.
! Authors:
!               Peter N. Brown and Alan C. Hindmarsh
!               Center for Applied Scientific Computing, L-561
!               Lawrence Livermore National Laboratory
!               Livermore, CA 94551
!               George D. Byrne
!               Illinois Institute of Technology
!               Chicago, IL 60616
! References:
! 1. P. N. Brown, G. D. Byrne, and A. C. Hindmarsh, "VODE: A Variable
!    Coefficient ODE Solver," SIAM J. Sci. Stat. Comput., 10 (1989),
!    pp. 1038-1051. Also, LLNL Report UCRL-98412, June 1988.
! 2. G. D. Byrne and A. C. Hindmarsh, "A Polyalgorithm for the
!    Numerical Solution of Ordinary Differential Equations,"
!    ACM Trans. Math. Software, 1 (1975), pp. 71-96.
! 3. A. C. Hindmarsh and G. D. Byrne, "EPISODE: An Effective Package
!    for the Integration of Systems of Ordinary Differential
!    Equations," LLNL Report UCID/30112, Rev. 1, April 1977.
! 4. G. D. Byrne and A. C. Hindmarsh, "EPISODEB: An Experimental
!    Package for the Integration of Systems of Ordinary Differential
!    Equations with Banded Jacobians," LLNL Report UCID/30132, April
!    1976.
! 5. A. C. Hindmarsh, "ODEPACK, a Systematized Collection of ODE
!    Solvers," in Scientific Computing, R. S. Stepleman et al., eds.,
!    North-Holland, Amsterdam, 1983, pp. 55-64.
! 6. K. R. Jackson and R. Sacks-Davis, "An Alternative Implementation
!    of Variable Step-Size Multistep Formulas for Stiff ODEs," ACM
!    Trans. Math. Software, 6 (1980), pp. 295-318.
!                     Summary of Usage
! Communication between the user and the DVODE package, for normal
! situations, is summarized here. This summary describes only a subset
! of the full set of options available. See the full description for
! details, including optional communication, nonstandard options,
! and instructions for special situations. See also the example
! problem (with program and output) following this summary.
! A. First provide a subroutine of the form:
!           SUBROUTINE F(NEQ, T, Y, YDOT)
!           REAL(KIND=WP) T, Y(NEQ), YDOT(NEQ)
! which supplies the vector function f by loading YDOT(i) with f(i).
! B. Next determine (or guess) whether or not the problem is stiff.
! Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
! whose real part is negative and large in magnitude, compared to the
! reciprocal of the t span of interest. If the problem is nonstiff,
! use a method flag MF = 10. If it is stiff, there are four standard
! choices for MF(21, 22, 24, 25), and DVODE requires the Jacobian
! matrix in some form. In these cases (MF > 0), DVODE will use a
! saved copy of the Jacobian matrix. If this is undesirable because of
! storage limitations, set MF to the corresponding negative value
! (-21, -22, -24, -25). (See full description of MF below.)
! The Jacobian matrix is regarded either as full (MF = 21 or 22),
! or banded (MF = 24 or 25). In the banded case, DVODE requires two
! half-bandwidth parameters ML and MU. These are, respectively, the
! widths of the lower and upper parts of the band, excluding the main
! diagonal. Thus the band consists of the locations (i,j) with
! i-ML <= j <= i+MU, and the full bandwidth is ML+MU+1.
! C. If the problem is stiff, you are encouraged to supply the Jacobian
! directly (MF = 21 or 24), but if this is not feasible, DVODE will
! compute it internally by difference quotients (MF = 22 or 25).
! If you are supplying the Jacobian, provide a subroutine of the form:
!           SUBROUTINE JAC(NEQ, T, Y, ML, MU, PD, NROWPD)
!           REAL(KIND=WP) T, Y(NEQ), PD(NROWPD,NEQ)
! which supplies df/dy by loading PD as follows:
!     For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),
! the partial derivative of f(i) with respect to y(j). (Ignore the
! ML and MU arguments in this case.)
!     For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with
! df(i)/dy(j), i.e. load the diagonal lines of df/dy into the rows of
! PD from the top down.
!     In either case, only nonzero elements need be loaded.
! D. Write a main program which calls subroutine DVODE once for
! each point at which answers are desired. This should also provide
! for possible use of logical unit 6 for output of error messages
! by DVODE. On the first call to DVODE, supply arguments as follows:
! F      = Name of subroutine for right-hand side vector f.
!          This name must be declared external in calling program.
! NEQ    = Number of first order ODEs.
! Y      = Array of initial values, of length NEQ.
! T      = The initial value of the independent variable.
! TOUT   = First point where output is desired (/= T).
! ITOL   = 1 or 2 according as ATOL(below) is a scalar or array.
! RTOL   = Relative tolerance parameter (scalar).
! ATOL   = Absolute tolerance parameter (scalar or array).
!          The estimated local error in Y(i) will be controlled so as
!          to be roughly less (in magnitude) than
!             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
!             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
!          Thus the local error test passes if, in each component,
!          either the absolute error is less than ATOL(or ATOL(i)),
!          or the relative error is less than RTOL.
!          Use RTOL = 0.0 for pure absolute error control, and
!          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
!          control. Caution: Actual (global) errors may exceed these
!          local tolerances, so choose them conservatively.
! ITASK  = 1 for normal computation of output values of Y at t = TOUT.
! ISTATE = Integer flag (input and output). Set ISTATE = 1.
! IOPT   = 0 to indicate no optional input used.
! JAC    = Name of subroutine for Jacobian matrix (MF = 21 or 24).
!          If used, this name must be declared external in calling
!          program. If not used, pass a dummy name.
! MF     = Method flag. Standard values are:
!          10 for nonstiff (Adams) method, no Jacobian used.
!          21 for stiff (BDF) method, user-supplied full Jacobian.
!          22 for stiff method, internally generated full Jacobian.
!          24 for stiff method, user-supplied banded Jacobian.
!          25 for stiff method, internally generated banded Jacobian.
! E. The output from the first call (or any call) is:
!      Y = Array of computed values of y(t) vector.
!      T = Corresponding value of independent variable (normally TOUT).
! ISTATE = 2  if DVODE was successful, negative otherwise.
!          -1 means excess work done on this call. (Perhaps wrong MF.)
!          -2 means excess accuracy requested. (Tolerances too small.)
!          -3 means illegal input detected. (See printed message.)
!          -4 means repeated error test failures. (Check all input.)
!          -5 means repeated convergence failures. (Perhaps bad
!             Jacobian supplied or wrong choice of MF or tolerances.)
!          -6 means error weight became zero during problem. (Solution
!             component I vanished, and ATOL or ATOL(I) = 0.)
! F. To continue the integration after a successful return, simply
! reset TOUT and call DVODE again. No other parameters need be reset.
!         Full Description of User Interface to DVODE
! The user interface to DVODE consists of the following parts.
! i.  The call sequence to subroutine DVODE, which is a driver
!      routine for the solver. This includes descriptions of both
!      the call sequence arguments and of user-supplied routines.
!      Following these descriptions is
!        * a description of optional input available through the
!          call sequence,
!        * a description of optional output (in the work arrays), and
!        * instructions for interrupting and restarting a solution.
! ii. Descriptions of other routines in the DVODE package that may be
!      (optionally) called by the user. These provide the ability to
!      alter error message handling, save and restore the internal
!      PRIVATE variables, and obtain specified derivatives of the
!      solution y(t).
! iii. Descriptions of PRIVATE variables to be declared in overlay
!      or similar environments.
! iv. Description of two routines in the DVODE package, either of
!      which the user may replace with his own version, if desired.
!      these relate to the measurement of errors.
! Part i. Call Sequence.
! The call sequence parameters used for input only are
!     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, JAC, MF,
! and those used for both input and output are
!     Y, T, ISTATE.
! The work arrays RUSER and IUSER are used for conditional and
! optional input and optional output. (The term output here refers
! to the return from subroutine DVODE to the user's calling program.)
! The legality of input parameters will be thoroughly checked on the
! initial call for the problem, but not checked thereafter unless a
! change in input parameters is flagged by ISTATE = 3 in the input.
! The descriptions of the call arguments are as follows.
! F      = The name of the user-supplied subroutine defining the
!          ODE system. The system must be put in the first-order
!          form dy/dt = f(t,y), where f is a vector-valued function
!          of the scalar t and the vector y. Subroutine F is to
!          compute the function f. It is to have the form
!               SUBROUTINE F(NEQ, T, Y, YDOT)
!               REAL(KIND=WP) T, Y(NEQ), YDOT(NEQ)
!          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
!          is output. Y and YDOT are arrays of length NEQ.
!          Subroutine F should not alter Y(1),...,Y(NEQ).
!          F must be declared EXTERNAL in the calling program.

!          If quantities computed in the F routine are needed
!          externally to DVODE, an extra call to F should be made
!          for this purpose, for consistent and accurate results.
!          If only the derivative dy/dt is needed, use DVINDY instead.
! NEQ    = The size of the ODE system (number of first order
!          ordinary differential equations). Used only for input.
!          NEQ may not be increased during the problem, but
!          can be decreased (with ISTATE = 3 in the input).
! Y      = A real array for the vector of dependent variables, of
!          length NEQ or more. Used for both input and output on the
!          first call (ISTATE = 1), and only for output on other calls.
!          On the first call, Y must contain the vector of initial
!          values. In the output, Y contains the computed solution
!          evaluated at T. If desired, the Y array may be used
!          for other purposes between calls to the solver.
!          This array is passed as the Y argument in all calls to
!          F and JAC.
! T      = The independent variable. In the input, T is used only on
!          the first call, as the initial point of the integration.
!          In the output, after each call, T is the value at which a
!          computed solution Y is evaluated (usually the same as TOUT).
!          On an error return, T is the farthest point reached.
! TOUT   = The next value of t at which a computed solution is desired.
!          Used only for input.
!          When starting the problem (ISTATE = 1), TOUT may be equal
!          to T for one call, then should /= T for the next call.
!          For the initial T, an input value of TOUT /= T is used
!          in order to determine the direction of the integration
!          (i.e. the algebraic sign of the step sizes) and the rough
!          scale of the problem. Integration in either direction
!          (forward or backward in t) is permitted.
!          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
!          the first call (i.e. the first call with TOUT /= T).
!          Otherwise, TOUT is required on every call.
!          If ITASK = 1, 3, or 4, the values of TOUT need not be
!          monotone, but a value of TOUT which backs up is limited
!          to the current internal t interval, whose endpoints are
!          TCUR - HU and TCUR. (See optional output, below, for
!          TCUR and HU.)
! ITOL   = An indicator for the type of error control. See
!          description below under ATOL. Used only for input.
! RTOL   = A relative error tolerance parameter, either a scalar or
!          an array of length NEQ. See description below under ATOL.
!          Input only.
! ATOL   = An absolute error tolerance parameter, either a scalar or
!          an array of length NEQ. Input only.
!          The input parameters ITOL, RTOL, and ATOL determine
!          the error control performed by the solver. The solver will
!          control the vector e = (e(i)) of estimated local errors
!          in Y, according to an inequality of the form
!                      rms-norm of (E(i)/EWT(i)) <= 1,
!          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
!          and the rms-norm (root-mean-square norm) here is
!          rms-norm(v) = sqrt(sum v(i)**2 / NEQ). Here EWT
!          is a vector of weights which must always be positive, and
!          the values of RTOL and ATOL should all be nonnegative.
!          The following table gives the types (scalar/array) of
!          RTOL and ATOL, and the corresponding form of EWT(i).
!           ITOL    RTOL       ATOL          EWT(i)
!            1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
!            2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
!            3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
!            4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
!          When either of these parameters is a scalar, it need not
!          be dimensioned in the user's calling program.
!          If none of the above choices (with ITOL, RTOL, and ATOL
!          fixed throughout the problem) is suitable, more general
!          error controls can be obtained by substituting
!          user-supplied routines for the setting of EWT and/or for
!          the norm calculation. See Part iv below.
!          If global errors are to be estimated by making a repeated
!          run on the same problem with smaller tolerances, then all
!          components of RTOL and ATOL(i.e. of EWT) should be scaled
!          down uniformly.
! ITASK  = An index specifying the task to be performed.
!          Input only. ITASK has the following values and meanings.
!          1  means normal computation of output values of y(t) at
!             t = TOUT(by overshooting and interpolating).
!          2  means take one step only and return.
!          3  means stop at the first internal mesh point at or
!             beyond t = TOUT and return.
!          4  means normal computation of output values of y(t) at
!             t = TOUT but without overshooting t = TCRIT.
!             TCRIT must be input as RUSER(1). TCRIT may be equal to
!             or beyond TOUT, but not behind it in the direction of
!             integration. This option is useful if the problem
!             has a singularity at or beyond t = TCRIT.
!          5  means take one step, without passing TCRIT, and return.
!             TCRIT must be input as RUSER(1).
!          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
!          (within roundoff), it will return T = TCRIT(exactly) to
!          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
!          in which case answers at T = TOUT are returned first).
! ISTATE = an index used for input and output to specify the
!          the state of the calculation.
!          In the input, the values of ISTATE are as follows.
!          1  means this is the first call for the problem
!             (initializations will be done). See note below.
!          2  means this is not the first call, and the calculation
!             is to continue normally, with no change in any input
!             parameters except possibly TOUT and ITASK.
!             (If ITOL, RTOL, and/or ATOL are changed between calls
!             with ISTATE = 2, the new values will be used but not
!             tested for legality.)
!          3  means this is not the first call, and the
!             calculation is to continue normally, but with
!             a change in input parameters other than
!             TOUT and ITASK. Changes are allowed in
!             NEQ, ITOL, RTOL, ATOL, IOPT, MF, ML, MU,
!             and any of the optional input except H0.
!             (See IUSER description for ML and MU.)
!          Caution:
!          If you make a call to DVODE_F90 with ISTATE=3, you will
!          first need to call SET_OPTS again, supplying the new
!          necessary option values.
!          Note:  A preliminary call with TOUT = T is not counted
!          as a first call here, as no initialization or checking of
!          input is done. (Such a call is sometimes useful to include
!          the initial conditions in the output.)
!          Thus the first call for which TOUT /= T requires
!          ISTATE = 1 in the input.
!          In the output, ISTATE has the following values and meanings.
!           1  means nothing was done, as TOUT was equal to T with
!              ISTATE = 1 in the input.
!           2  means the integration was performed successfully.
!          -1  means an excessive amount of work (more than MXSTEP
!              steps) was done on this call, before completing the
!              requested task, but the integration was otherwise
!              successful as far as T. (MXSTEP is an optional input
!              and is normally 5000.)  To continue, the user may
!              simply reset ISTATE to a value > 1 and call again.
!              (The excess work step counter will be reset to 0.)
!              In addition, the user may increase MXSTEP to avoid
!              this error return. (See optional input below.)
!          -2  means too much accuracy was requested for the precision
!              of the machine being used. This was detected before
!              completing the requested task, but the integration
!              was successful as far as T. To continue, the tolerance
!              parameters must be reset, and ISTATE must be set
!              to 3. The optional output TOLSF may be used for this
!              purpose. (Note: If this condition is detected before
!              taking any steps, then an illegal input return
!              (ISTATE = -3) occurs instead.)
!          -3  means illegal input was detected, before taking any
!              integration steps. See written message for details.
!              Note:  If the solver detects an infinite loop of calls
!              to the solver with illegal input, it will cause
!              the run to stop.
!          -4  means there were repeated error test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              The problem may have a singularity, or the input
!              may be inappropriate.
!          -5  means there were repeated convergence test failures on
!              one attempted step, before completing the requested
!              task, but the integration was successful as far as T.
!              This may be caused by an inaccurate Jacobian matrix,
!              if one is being used.
!          -6  means EWT(i) became zero for some i during the
!              integration. Pure relative error control (ATOL(i)=0.0)
!              was requested on a variable which has now vanished.
!              The integration was successful as far as T.
!          Note:  Since the normal output value of ISTATE is 2,
!          it does not need to be reset for normal continuation.
!          Also, since a negative input value of ISTATE will be
!          regarded as illegal, a negative output value requires the
!          user to change it, and possibly other input, before
!          calling the solver again.
! IOPT   = An integer flag to specify whether or not any optional
!          input is being used on this call. Input only.
!          The optional input is listed separately below.
!          IOPT = 0 means no optional input is being used.
!                   Default values will be used in all cases.
!          IOPT = 1 means optional input is being used.
! RUSER  = A real working array (real(wp)).
!          The length of RUSER must be at least 22
!             20 + NYH * (MAXORD + 1) where
!          NYH    = the initial value of NEQ,
!          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
!                   smaller value is given as an optional input),
!          The first 22 words of RUSER are reserved for conditional
!          and optional input and optional output.

!          The following word in RUSER is a conditional input:
!            RUSER(1) = TCRIT = critical value of t which the solver
!                       is not to overshoot. Required if ITASK is
!                       4 or 5, and ignored otherwise. (See ITASK.)
! IUSER  = An integer work array. The length of IUSER must be at least 30.
!             30        if MITER = 0 or 3 (MF = 10, 13, 20, 23), or
!             30 + NEQ  otherwise (ABS(MF) = 11,12,14,15,16,17,21,22,
!             24,25,26,27).
!          The first 30 words of IUSER are reserved for conditional and
!          optional input and optional output.

!          The following 2 words in IUSER are conditional input:
!            IUSER(1) = ML  These are the lower and upper
!            IUSER(2) = MU  half-bandwidths, respectively, of the
!                       banded Jacobian, excluding the main diagonal.
!                       The band is defined by the matrix locations
!                       (i,j) with i-ML <= j <= i+MU. ML and MU
!                       must satisfy  0 <= ML,MU  <= NEQ-1.
!                       These are required if MITER is 4 or 5, and
!                       ignored otherwise. ML and MU may in fact be
!                       the band parameters for a matrix to which
!                       df/dy is only approximately equal.
! Note:  The work arrays must not be altered between calls to DVODE
! for the same problem, except possibly for the conditional and
! optional input, and except for the last 3*NEQ words of RUSER.
! The latter space is used for internal scratch space, and so is
! available for use by the user outside DVODE between calls, if
! desired (but not for use by F or JAC).
! JAC    = The name of the user-supplied routine (MITER = 1 or 4 or 6)
!          to compute the Jacobian matrix, df/dy, as a function of
!          the scalar t and the vector y. It is to have the form
!               SUBROUTINE JAC(NEQ, T, Y, ML, MU, PD, NROWPD)
!               REAL(KIND=WP) T, Y(NEQ), PD(NROWPD,NEQ)
!          where NEQ, T, Y, ML, MU, and NROWPD are input and the array
!          PD is to be loaded with partial derivatives (elements of the
!          Jacobian matrix) in the output. PD must be given a first
!          dimension of NROWPD. T and Y have the same meaning as in
!          Subroutine F.
!               In the full matrix case (MITER = 1), ML and MU are
!          ignored, and the Jacobian is to be loaded into PD in
!          columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
!               In the band matrix case (MITER = 4), the elements
!          within the band are to be loaded into PD in columnwise
!          manner, with diagonal lines of df/dy loaded into the rows
!          of PD. Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).
!          ML and MU are the half-bandwidth parameters. (See IUSER).
!          The locations in PD in the two triangular areas which
!          correspond to nonexistent matrix elements can be ignored
!          or loaded arbitrarily, as they are overwritten by DVODE.
!               In the sparse matrix case the elements of the matrix
!          are determined by the sparsity structure given by the
!          IA and JA pointer arrays. Refer to the documentation
!          prologue for SET_OPTS for a description of the arguments
!          for JAC since they differ from the dense and banded cases.
!               JAC need not provide df/dy exactly. A crude
!          approximation (possibly with a smaller bandwidth) will do.
!               In either case, PD is preset to zero by the solver,
!          so that only the nonzero elements need be loaded by JAC.
!          Each call to JAC is preceded by a call to F with the same
!          arguments NEQ, T, and Y. Thus to gain some efficiency,
!          intermediate quantities shared by both calculations may be
!          saved in a user common block by F and not recomputed by JAC,
!          if desired. Also, JAC may alter the Y array, if desired.
!          JAC must be declared external in the calling program.
! MF     = The method flag. Used only for input. The legal values of
!          MF are 10, 11, 12, 13, 14, 15, 16, 17, 20, 21, 22, 23, 24,
!          25, 26, 27, -11, -12, -14, -15, -21, -22, -24, -25, -26,
!          -27.
!          MF is a signed two-digit integer, MF = JSV*(10*METH+MITER).
!          JSV = SIGN(MF) indicates the Jacobian-saving strategy:
!            JSV =  1 means a copy of the Jacobian is saved for reuse
!                     in the corrector iteration algorithm.
!            JSV = -1 means a copy of the Jacobian is not saved
!                     (valid only for MITER = 1, 2, 4, or 5).
!          METH indicates the basic linear multistep method:
!            METH = 1 means the implicit Adams method.
!            METH = 2 means the method based on backward
!                     differentiation formulas (BDF-s).
!          MITER indicates the corrector iteration method:
!            MITER = 0 means functional iteration (no Jacobian matrix
!                      is involved).
!            MITER = 1 means chord iteration with a user-supplied
!                      full (NEQ by NEQ) Jacobian.
!            MITER = 2 means chord iteration with an internally
!                      generated (difference quotient) full Jacobian
!                      (using NEQ extra calls to F per df/dy value).
!            MITER = 3 means chord iteration with an internally
!                      generated diagonal Jacobian approximation
!                      (using 1 extra call to F per df/dy evaluation).
!            MITER = 4 means chord iteration with a user-supplied
!                      banded Jacobian.
!            MITER = 5 means chord iteration with an internally
!                      generated banded Jacobian (using ML+MU+1 extra
!                      calls to F per df/dy evaluation).
!            MITER = 6 means chord iteration with a user-supplied
!                      sparse Jacobian.
!            MITER = 7 means chord iteration with an internally
!                      generated sparse Jacobian
!          If MITER = 1, 4, or 6 the user must supply a subroutine
!          JAC(the name is arbitrary) as described above under JAC.
!          For other values of MITER, a dummy argument can be used.
!                         Optional Input
! The following is a list of the optional input provided for in the
! call sequence. (See also Part ii.)  For each such input variable,
! this table lists its name as used in this documentation, its
! location in the call sequence, its meaning, and the default value.
! The use of any of this input requires IOPT = 1, and in that
! case all of this input is examined. A value of zero for any of
! these optional input variables will cause the default value to be
! used. Thus to use a subset of the optional input, simply preload
! locations 5 to 10 in RUSER and IUSER to 0.0 and 0, respectively,
! and then set those of interest to nonzero values.
! NAME    LOCATION      MEANING AND DEFAULT VALUE
! H0      RUSER(5)  The step size to be attempted on the first step.
!                   The default value is determined by the solver.
! HMAX    RUSER(6)  The maximum absolute step size allowed.
!                   The default value is infinite.
! HMIN    RUSER(7)  The minimum absolute step size allowed.
!                   The default value is 0. (This lower bound is not
!                   enforced on the final step before reaching TCRIT
!                   when ITASK = 4 or 5.)
! MAXORD  IUSER(5)  The maximum order to be allowed. The default
!                   value is 12 if METH = 1, and 5 if METH = 2.
!                   If MAXORD exceeds the default value, it will
!                   be reduced to the default value.
!                   If MAXORD is changed during the problem, it may
!                   cause the current order to be reduced.
! MXSTEP  IUSER(6)  Maximum number of (internally defined) steps
!                   allowed during one call to the solver.
!                   The default value is 5000.
! MXHNIL  IUSER(7)  Maximum number of messages printed (per problem)
!                   warning that T + H = T on a step (H = step size).
!                   This must be positive to result in a non-default
!                   value. The default value is 10.
!                          Optional Output
! As optional additional output from DVODE, the variables listed
! below are quantities related to the performance of DVODE
! which are available to the user. These are communicated by way of
! the work arrays, but also have internal mnemonic names as shown.
! Except where stated otherwise, all of this output is defined
! on any successful return from DVODE, and on any return with
! ISTATE = -1, -2, -4, -5, or -6. On an illegal input return
! (ISTATE = -3), they will be unchanged from their existing values
! (if any), except possibly for TOLSF, LENRW, and LENIW.
! On any error return, output relevant to the error will be defined,
! as noted below.
! NAME    LOCATION      MEANING
! HU      RUSER(11) The step size in t last used (successfully).
! HCUR    RUSER(12) The step size to be attempted on the next step.
! TCUR    RUSER(13) The current value of the independent variable
!                   which the solver has actually reached, i.e. the
!                   current internal mesh point in t. In the output,
!                   TCUR will always be at least as far from the
!                   initial value of t as the current argument T,
!                   but may be farther (if interpolation was done).
! TOLSF   RUSER(14) A tolerance scale factor, greater than 1.0,
!                   computed when a request for too much accuracy was
!                   detected (ISTATE = -3 if detected at the start of
!                   the problem, ISTATE = -2 otherwise). If ITOL is
!                   left unaltered but RTOL and ATOL are uniformly
!                   scaled up by a factor of TOLSF for the next call,
!                   then the solver is deemed likely to succeed.
!                   (The user may also ignore TOLSF and alter the
!                   tolerance parameters in any other way appropriate.)
! NST     IUSER(11) The number of steps taken for the problem so far.
! NFE     IUSER(12) The number of f evaluations for the problem so far.
! NJE     IUSER(13) The number of Jacobian evaluations so far.
! NQU     IUSER(14) The method order last used (successfully).
! NQCUR   IUSER(15) The order to be attempted on the next step.
! IMXER   IUSER(16) The index of the component of largest magnitude in
!                   the weighted local error vector (e(i)/EWT(i)),
!                   on an error return with ISTATE = -4 or -5.
! LENRW   IUSER(17) The length of RUSER actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
! LENIW   IUSER(18) The length of IUSER actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
! NLU     IUSER(19) The number of matrix LU decompositions so far.
! NNI     IUSER(20) The number of nonlinear (Newton) iterations so far.
! NCFN    IUSER(21) The number of convergence failures of the nonlinear
!                   solver so far.
! NETF    IUSER(22) The number of error test failures of the integrator
!                   so far.
! The following two arrays are segments of the RUSER array which
! may also be of interest to the user as optional output.
! For each array, the table below gives its internal name,
! its base address in RUSER, and its description.
!                  Interrupting and Restarting
! If the integration of a given problem by DVODE is to be interrupted
! and then later continued, such as when restarting an interrupted run
! or alternating between two or more ODE problems, the user should save,
! following the return from the last DVODE call prior to the
! interruption, the contents of the call sequence variables and
! internal PRIVATE variables, and later restore these values before the
! next DVODE call for that problem. To save and restore the PRIVATE
! variables, use subroutine DVSRCO, as described below in part ii.
! In addition, if non-default values for either LUN or MFLAG are
! desired, an extra call to XSETUN and/or XSETF should be made just
! before continuing the integration. See Part ii below for details.
! Part ii. Other Routines Callable.
! The following are optional calls which the user may make to
! gain additional capabilities in conjunction with DVODE.
! (The routines XSETUN and XSETF are designed to conform to the
! SLATEC error handling package.)
!     FORM OF CALL                  FUNCTION
!  CALL XSETUN(LUN)           Set the logical unit number, LUN, for
!                             output of messages from DVODE, if
!                             the default is not desired.
!                             The default value of LUN is 6.
!  CALL XSETF(MFLAG)          Set a flag to control the printing of
!                             messages by DVODE.
!                             MFLAG = 0 means do not print. (Danger:
!                             This risks losing valuable information.)
!                             Either of the above calls may be made at
!                             any time and will take effect immediately.
!  CALL DVINDY(...)           Provide derivatives of y, of various
!                             orders, at a specified point T, if
!                             desired. It may be called only after
!                             a successful return from DVODE.
! The detailed instructions for using DVINDY are as follows.
! The form of the call is:
!      CALL DVINDY(T,K,DKY,IFLAG)
! The input parameters are:
! T         = Value of independent variable where answers are desired
!             (normally the same as the T last returned by DVODE).
!             For valid results, T must lie between TCUR - HU and TCUR.
!             (See optional output for TCUR and HU.)
! K         = Integer order of the derivative desired. K must satisfy
!             0 <= K <= NQCUR, where NQCUR is the current order
!             (see optional output). The capability corresponding
!             to K = 0, i.e. computing y(T), is already provided
!             by DVODE directly. Since NQCUR >= 1, the first
!             derivative dy/dt is always available with DVINDY.
! The output parameters are:
! DKY       = A real array of length NEQ containing the computed value
!             of the K-th derivative of y(t).
! IFLAG     = Integer flag, returned as 0 if K and T were legal,
!             -1 if K was illegal, and -2 if T was illegal.
!             On an error return, a message is also written.
! Part iii. Optionally Replaceable Solver Routines.
! Below are descriptions of two routines in the DVODE package which
! relate to the measurement of errors. Either routine can be
! replaced by a user-supplied version, if desired. However, since such
! a replacement may have a major impact on performance, it should be
! done only when absolutely necessary, and only with great caution.
! (Note: The means by which the package version of a routine is
! superseded by the user's version may be system-dependent.)
! (a) DEWSET.
! The following subroutine is called just before each internal
! integration step, and sets the array of error weights, EWT, as
! described under ITOL/RTOL/ATOL above:
!     SUBROUTINE DEWSET(NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
! where NEQ, ITOL, RTOL, and ATOL are as in the DVODE call sequence,
! YCUR contains the current dependent variable vector, and
! EWT is the array of weights set by DEWSET.
! If the user supplies this subroutine, it must return in EWT(i)
! (i = 1,...,NEQ) a positive quantity suitable for comparison with
! errors in Y(i). The EWT array returned by DEWSET is passed to the
! DVNORM function (See below.), and also used by DVODE in the
! computation of the optional output IMXER, the diagonal Jacobian
! approximation, and the increments for difference quotient Jacobians.
! In the user-supplied version of DEWSET, it may be desirable to use
! the current values of derivatives of y. Derivatives up to order NQ
! are available from the history array YH, described above under
! Optional Output. In DEWSET, YH is identical to the YCUR array,
! extended to NQ + 1 columns with a column length of NYH and scale
! factors of h**j/factorial(j). On the first call for the problem,
! given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
! NYH is the initial value of NEQ. Thus, for example, the current
! value of dy/dt can be obtained as YCUR(NYH+i)/H  (i=1,...,NEQ)
! (and the division by H is unnecessary when NST = 0).
! (b) DVNORM.
! The following is a function which computes the weighted
! root-mean-square norm of a vector v:
!     D = DVNORM(N, V, W)
! where:
!   N = the length of the vector,
!   V = real array of length N containing the vector,
!   W = real array of length N containing weights,
!   D = sqrt((1/N) * sum(V(i)*W(i))**2).
! DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
! EWT is as set by subroutine DEWSET.
! If the user supplies this routine, it should return a nonnegative
! value of DVNORM suitable for use in the error control in DVODE.
! None of the arguments should be altered by DVNORM.
! For example, a user-supplied DVNORM function might:
!   -substitute a max-norm of (V(i)*W(i)) for the rms-norm, or
!   -ignore some components of V in the norm, with the effect of
!    suppressing the error control on those components of Y.
!_______________________________________________________________________
! Other Routines in the DVODE Package

! In addition to subroutine DVODE, the DVODE package includes the
! following subroutines and function routines (not user callable):
!  DVHIN       computes an approximate step size for the initial step.
!  DVINDY_CORE computes an interpolated value of the y vector at t=TOUT.
!  DVINDY      computes an interpolated value of the y vector at t=TOUT.
!              (user callable)
!  DVSTEP      is the core integrator, which does one step of the
!              integration and the associated error control.
!  DVSET       sets all method coefficients and test constants.
!  DVNLSD,     solves the underlying nonlinear system -- the corrector.
!  DVNLSS28
!  DVJAC,      computes and preprocesses the Jacobian matrix J = df/dy
!  DVJACS28    and the Newton iteration matrix P = I - (h/l1)*J.
!  DVSOL,      manages solution of linear system in chord iteration.
!  DVSOLS28
!  DVJUST      adjusts the history array on a change of order.
!  DEWSET      sets the error weight vector EWT before each step.
!  DVNORM      computes the weighted r.m.s. norm of a vector.
!  DACOPY      is a routine to copy a two-dimensional array to another.
!  DGEFA_F90 and DGESL_F90 are routines from LINPACK for solving full
!              systems of linear algebraic equations.
!  DGBFA_F90 and DGBSL_F90 are routines from LINPACK for solving banded
!              linear systems.
!  DAXPY_F90, DSCAL_F90, and DCOPY_F90 are basic linear algebra modules
!              (BLAS).
!  DVCHECK     does preliminary checking for roots, and serves as an
!              interface between subroutine DVODE_F90 and subroutine
!              DVROOTS.
!  DVROOTS     finds the leftmost root of a set of functions.
! ______________________________________________________________________
! Section 10.  Example Usage
!
! MODULE example1
! The following is a simple example problem, with the coding
! needed for its solution by DVODE_F90. The problem is from
! chemical kinetics, and consists of the following three rate
! equations:
!     dy1/dt = -.04d0*y1 + 1.d4*y2*y3
!     dy2/dt = .04d0*y1 - 1.d4*y2*y3 - 3.d7*y2**2
!     dy3/dt = 3.d7*y2**2
! on the interval from t = 0.0d0 to t = 4.d10, with initial
! conditions y1 = 1.0d0, y2 = y3 = 0.0d0. The problem is stiff.
! The following coding solves this problem with DVODE_F90,
! using a user supplied Jacobian and printing results at
! t = .4, 4.,...,4.d10. It uses ITOL = 2 and ATOL much smaller
! for y2 than y1 or y3 because y2 has much smaller values. At
! the end of the run, statistical quantities of interest are
! printed. (See optional output in the full DVODE description
! below.) Output is written to the file example1.dat.
! CONTAINS
!     SUBROUTINE FEX(NEQ, T, Y, YDOT)
!     IMPLICIT NONE
!     INTEGER NEQ
!     DOUBLE PRECISION T, Y, YDOT
!     DIMENSION Y(NEQ), YDOT(NEQ)
!     YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
!     YDOT(3) = 3.E7*Y(2)*Y(2)
!     YDOT(2) = -YDOT(1) - YDOT(3)
!     RETURN
!     END SUBROUTINE FEX
!     SUBROUTINE JEX(NEQ, T, Y, ML, MU, PD, NRPD)
!     IMPLICIT NONE
!     INTEGER NEQ,ML,MU,NRPD
!     DOUBLE PRECISION PD, T, Y
!     DIMENSION Y(NEQ), PD(NRPD,NEQ)
!     PD(1,1) = -.04D0
!     PD(1,2) = 1.D4*Y(3)
!     PD(1,3) = 1.D4*Y(2)
!     PD(2,1) = .04D0
!     PD(2,3) = -PD(1,3)
!     PD(3,2) = 6.E7*Y(2)
!     PD(2,2) = -PD(1,2) - PD(3,2)
!     RETURN
!     END SUBROUTINE JEX
! END MODULE example1
!******************************************************************

!     PROGRAM runexample1
!     USE DVODE_F90_M
!     USE example1
!     IMPLICIT NONE
!     DOUBLE PRECISION ATOL, RTOL, T, TOUT, Y, RSTATS
!     INTEGER NEQ, ITASK, ISTATE, ISTATS, IOUT, IERROR, I
!     DIMENSION Y(3), ATOL(3), RSTATS(22), ISTATS(31)
!     TYPE(VODE_OPTS) :: OPTIONS
!     OPEN(UNIT=6, FILE = 'example1.dat')
!     IERROR = 0
!     NEQ = 3
!     Y(1) = 1.0D0
!     Y(2) = 0.0D0
!     Y(3) = 0.0D0
!     T = 0.0D0
!     TOUT = 0.4D0
!     RTOL = 1.D-4
!     ATOL(1) = 1.D-8
!     ATOL(2) = 1.D-14
!     ATOL(3) = 1.D-6
!     ITASK = 1
!     ISTATE = 1
!     OPTIONS = SET_OPTS(DENSE_J=.TRUE.,ABSERR_VECTOR=ATOL, &
!       RELERR=RTOL, USER_SUPPLIED_JACOBIAN=.TRUE.)
!     DO IOUT = 1,12
!       CALL DVODE_F90(FEX,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,J_FCN=JEX)
!       CALL GET_STATS(RSTATS,ISTATS)
!       WRITE(6,63)T,Y(1),Y(2),Y(3)
!       DO I = 1, NEQ
!          IF (Y(I) < 0.0D0) IERROR = 1
!       END DO
!       IF (ISTATE < 0) THEN
!          WRITE(6,64)ISTATE
!          STOP
!       END IF
!       TOUT = TOUT*10.0D0
!     END DO
!     WRITE(6,60) ISTATS(11),ISTATS(12),ISTATS(13),ISTATS(19), &
!                 ISTATS(20),ISTATS(21),ISTATS(22)
!     IF (IERROR == 1) THEN
!        WRITE(6,61)
!     ELSE
!        WRITE(6,62)
!     END IF
! 60  FORMAT(/'  No. steps =',I4,'   No. f-s =',I4,        &
!             '  No. J-s =',I4,'   No. LU-s =',I4/         &
!             '  No. nonlinear iterations =',I4/           &
!             '  No. nonlinear convergence failures =',I4/ &
!             '  No. error test failures =',I4/)
! 61  FORMAT(/' An error occurred.')
! 62  FORMAT(/' No errors occurred.')
! 63  FORMAT(' At t =',D12.4,'   y =',3D14.6)
! 64  FORMAT(///' Error halt: ISTATE =',I3)
!     STOP
!     END PROGRAM runexample1
!
! MODULE example2
! The following is a modification of the previous example
! program to illustrate root finding. The problem is from
! chemical kinetics, and consists of the following three
! rate equations:
!     dy1/dt = -.04d0*y1 + 1.d4*y2*y3
!     dy2/dt = .04d0*y1 - 1.d4*y2*y3 - 3.d7*y2**2
!     dy3/dt = 3.d7*y2**2
! on the interval from t = 0.0d0 to t = 4.d10, with initial
! conditions y1 = 1.0d0, y2 = y3 = 0.0d0. The problem is stiff.
! In addition, we want to find the values of t, y1, y2,
! and y3 at which:
!   (1) y1 reaches the value 1.d-4, and
!   (2) y3 reaches the value 1.d-2.
! The following coding solves this problem with DVODE_F90
! using an internally generated dense Jacobian and
! printing results at t = .4, 4., ..., 4.d10, and at the
! computed roots. It uses ITOL = 2 and ATOL much smaller
! for y2 than y1 or y3 because y2 has much smaller values.
! At the end of the run, statistical quantities of interest
! are printed (see optional outputs in the full description
! below). Output is written to the file example2.dat.
! CONTAINS
!     SUBROUTINE FEX (NEQ, T, Y, YDOT)
!     IMPLICIT NONE
!     INTEGER NEQ
!     DOUBLE PRECISION T, Y, YDOT
!     DIMENSION Y(3), YDOT(3)
!     YDOT(1) = -0.04D0*Y(1) + 1.0D4*Y(2)*Y(3)
!     YDOT(3) = 3.0D7*Y(2)*Y(2)
!     YDOT(2) = -YDOT(1) - YDOT(3)
!     RETURN
!     END SUBROUTINE FEX
!     SUBROUTINE GEX (NEQ, T, Y, NG, GOUT)
!     IMPLICIT NONE
!     INTEGER NEQ, NG
!     DOUBLE PRECISION T, Y, GOUT
!     DIMENSION Y(3), GOUT(2)
!     GOUT(1) = Y(1) - 1.0D-4
!     GOUT(2) = Y(3) - 1.0D-2
!     RETURN
!     END SUBROUTINE GEX
! END MODULE example2
!******************************************************************
!     PROGRAM runexample2
!     USE DVODE_F90_M
!     USE example2
!     IMPLICIT NONE
!     INTEGER ITASK, ISTATE, NG, NEQ, IOUT, JROOT, ISTATS, &
!     IERROR, I
!     DOUBLE PRECISION ATOL, RTOL, RSTATS, T, TOUT, Y
!     DIMENSION Y(3), ATOL(3), RSTATS(22), ISTATS(31), JROOT(2)
!     TYPE(VODE_OPTS) :: OPTIONS
!     OPEN (UNIT=6, FILE='example2.dat')
!     IERROR = 0
!     NEQ = 3
!     Y(1) = 1.0D0
!     Y(2) = 0.0D0
!     Y(3) = 0.0D0
!     T = 0.0D0
!     TOUT = 0.4D0
!     RTOL = 1.0D-4
!     ATOL(1) = 1.0D-8
!     ATOL(2) = 1.0D-12
!     ATOL(3) = 1.0D-8
!     ITASK = 1
!     ISTATE = 1
!     NG = 2
!     OPTIONS = SET_OPTS(DENSE_J=.TRUE.,RELERR=RTOL, &
!       ABSERR_VECTOR=ATOL,NEVENTS=NG)
!     DO 40 IOUT = 1,12
! 10    CONTINUE
!       CALL DVODE_F90(FEX,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTIONS,G_FCN=GEX)
!       CALL GET_STATS(RSTATS, ISTATS, NG, JROOT)
!       WRITE(6,20) T, Y(1), Y(2), Y(3)
!       DO I = 1, NEQ
!          IF (Y(I) < 0.0D0) IERROR = 1
!       END DO
! 20    FORMAT(' At t =',D12.4,'   Y =',3D14.6)
!       IF (ISTATE < 0) GOTO 60
!       IF (ISTATE == 2) GOTO 40
!       WRITE(6,30) JROOT(1),JROOT(2)
! 30    FORMAT(5X,' The above line is a root, JROOT =',2I5)
!       ISTATE = 2
!       GOTO 10
! 40  TOUT = TOUT*10.0D0
!     WRITE(6,50) ISTATS(11), ISTATS(12), ISTATS(13), ISTATS(10)
!     IF (IERROR == 1) THEN
!        WRITE(6,61)
!     ELSE
!        WRITE(6,62)
!     END IF
! 50  FORMAT(/' No. steps =',I4,'  No. f-s =',I4,'  No. J-s =',I4, &
!     '  No. g-s =',I4/)
!     STOP
! 60  WRITE(6,70) ISTATE
! 61  FORMAT(/' An error occurred.')
! 62  FORMAT(/' No errors occurred.')
! 70  FORMAT(///' Error halt.. ISTATE =',I3)
!     STOP
!     END PROGRAM runexample2
!_______________________________________________________________________
! BEGINNING OF DVODE_F90_M PRIVATE SECTION.
! Note: This global information is used throughout DVODE_F90.
!_______________________________________________________________________
! JACSPDB arrays and parameters.
  logical, private :: USE_JACSP, LIKE_ORIGINAL_VODE
  integer, private :: INFODS, LIWADS, MAXGRPDS, MINGRPDS, NRFJACDS,    &
    NCFJACDS, LWKDS, LIWKDS
  integer, allocatable, private :: INDROWDS(:), INDCOLDS(:),           &
    NGRPDS(:), IPNTRDS(:), JPNTRDS(:), IWADS(:), IWKDS(:), IOPTDS(:)
  real (WP), allocatable, private :: YSCALEDS(:), WKDS(:), FACDS(:)
  real (WP), private :: U125, U325
!_______________________________________________________________________
  logical, parameter, private :: USE_MA48_FOR_SPARSE=.false.
!_______________________________________________________________________
! *****MA48 build change point. Replace the above statement.
! LOGICAL, PARAMETER, PRIVATE :: USE_MA48_FOR_SPARSE=.TRUE.
!_______________________________________________________________________
! *****MA48 build change point. Insert these statements.
! MA48 type declarations:
! TYPE(ZD01_TYPE) MATRIX
! TYPE(MA48_CONTROL) CONTROL
! TYPE(MA48_FACTORS) FACTORS
! TYPE(MA48_AINFO) AINFO
! TYPE(MA48_FINFO) FINFO
! TYPE(MA48_SINFO) SINFO
!_______________________________________________________________________
! .. Parameters ..
!     IPCUTH_MAX - maximum number of times the solver will halve the
!                  stepsize to prevent an infeasible prediction if
!                  solution bounds are used
!     KFC        - maximum number of consecutive convergence failures
!                  before crashing the order
!     KFH        - maximum number of consecutive error test failures
!                  before giving up (changed from 7 to 15)
!     MAXCOR     - maximum number of corrections
!     MSBP       - maximum number of steps before forming new P matrix
!     MXNCF      - maximum number of consecutive convergence failures
!                  before giving up
!     MXHNLO     - maximum number of T+H=T messages
!     MXSTP0     - maximum number of integration steps
!     L*****     - lengths and pointers for some internal arrays
!     INTEGER, PARAMETER, PRIVATE :: KFC = -3, KFH = -7, LENIV1 = 33,    &
      integer, parameter, private :: IPCUTH_MAX = 100, KFC = -3,         &
        KFH = -15, LENIV1 = 33,                                          &
        LENIV2 = 8, LENRV1 = 48, LENRV2 = 1, LIWUSER = 30, LRWUSER = 22, &
        MAXCOR = 3, MAX_ARRAY_SIZE = 900000000, MSBP = 20, MXHNL0 = 10,  &
        MXNCF = 10, MXSTP0 = 5000
!_______________________________________________________________________
! *****LAPACK build change point. Use .TRUE. for LAPACK.
!     LOGICAL, PARAMETER, PRIVATE :: USE_LAPACK = .TRUE.
!_______________________________________________________________________
      real (WP), parameter, private :: ADDON = 1.0E-6_WP
      real (WP), parameter, private :: BIAS1 = 6.0_WP
      real (WP), parameter, private :: BIAS2 = 6.0_WP
      real (WP), parameter, private :: BIAS3 = 10.0_WP
      real (WP), parameter, private :: CCMAX = 0.3_WP
      real (WP), parameter, private :: CORTES = 0.1_WP
      real (WP), parameter, private :: CRDOWN = 0.3_WP
      real (WP), parameter, private :: ETACF = 0.25_WP
      real (WP), parameter, private :: ETAMIN = 0.1_WP
      real (WP), parameter, private :: ETAMX1 = 1.0E4_WP
      real (WP), parameter, private :: ETAMX2 = 10.0_WP
      real (WP), parameter, private :: ETAMX3 = 10.0_WP
      real (WP), parameter, private :: ETAMXF = 0.2_WP
      real (WP), parameter, private :: FIVE = 5.0_WP
      real (WP), parameter, private :: FOUR = 4.0_WP
      real (WP), parameter, private :: HALF = 0.5_WP
      real (WP), parameter, private :: HUN = 100.0_WP
      real (WP), parameter, private :: HUNDRETH = 0.01_WP
      real (WP), parameter, private :: ONE = 1.0_WP
      real (WP), parameter, private :: ONEPSM = 1.00001_WP
      real (WP), parameter, private :: PT1 = 0.1_WP
      real (WP), parameter, private :: PT2 = 0.2_WP
      real (WP), parameter, private :: RDIV = 2.0_WP
      real (WP), parameter, private :: SIX = 6.0_WP
      real (WP), parameter, private :: TEN = 10.0_WP
      real (WP), parameter, private :: TENTH = 0.1_WP
      real (WP), parameter, private :: THOU = 1000.0_WP
      real (WP), parameter, private :: THRESH = 1.5_WP
      real (WP), parameter, private :: TWO = 2.0_WP
      real (WP), parameter, private :: ZERO = 0.0_WP

! Beginning of DVODE_F90 interface.
! ..
! .. Generic Interface Blocks ..
      interface DVODE_F90

!         VODE_F90 is the interface subroutine that is actually invoked
!         when the user calls DVODE_F90. It in turn calls subroutine
!         DVODE which is the driver that directs all the work.
          module procedure VODE_F90

!         GET_STATS can be called to gather integration statistics.
          module procedure GET_STATS

!         DVINDY can be called to interpolate the solution and derivative.
          module procedure DVINDY

!         RELEASE_ARRAYS can be called to release/deallocate the work arrays.
          module procedure RELEASE_ARRAYS

!         SET_IAJA can be called to set sparse matrix descriptor arrays.
          module procedure SET_IAJA

!         USERSETS_IAJA can be called to set sparse matrix descriptor arrays.
          module procedure USERSETS_IAJA

!         CHECK_STAT can be called to stop if a storage allocation or
!         deallocation error occurs.
          module procedure CHECK_STAT

!         JACSP can be called to calculate a Jacobian using Doug Salane's
!         algoritm 
          module procedure JACSP

!         DVDSM can be called to calculate sparse pointer arrays needed
!         by JACSP
          module procedure DVDSM

      end interface
! ..
! .. Derived Type Declarations ..
      type, public :: VODE_OPTS
        real (WP), dimension (:), pointer :: ATOL, RTOL
        integer :: MF, METH, MITER, MOSS, ITOL, IOPT, NG
        logical :: DENSE, BANDED, SPARSE
      end type VODE_OPTS
! ..
! .. Local Scalars ..
!_______________________________________________________________________
! *****MA48 build change point. Insert these statements.
!     For communication with subroutine ma48_control_array:
!     REAL (WP), PUBLIC :: COPY_OF_U_PIVOT
!_______________________________________________________________________
      real (WP), private :: ACNRM, ALPHA, BIG, BIG1, CCMXJ, CGCE, CONP, CRATE,  &
        DRC, DRES, DXMAX, EPS, ERRMAX, ETA, ETAMAX, FRACINT, FRACSUB, H, HMIN,  &
        HMXI, HNEW, HSCAL, HU, MEPS, MRESID, MRMIN, PRL1, RC, RESID, RL1,       &
        RMIN, SETH, T0ST, THEMAX, TLAST, TN, TOL, TOL1, TOUTC, UMAX, UROUND,    &
        U_PIVOT, X2, WM1, WM2
      integer, private :: ADDTOJA, ADDTONNZ, CONSECUTIVE_CFAILS,                &
        CONSECUTIVE_EFAILS, ELBOW_ROOM, IADIM, IANPIV, IAVPIV,                  &
        ICF, ICNCP, IFAIL, IMAX, IMIN, INEWJ, INIT, IPUP, IRANK, IRFND, IRNCP,  &
        ISTART, ISTATC, ITASKC, JADIM, JCUR, JMIN, JSTART, JSV, KFLAG, KOUNTL,  &
        KUTH, L, LARGE, LAST, LENIGP, LICN_ALL, LIRN_ALL, LIW, LIWM, LMAX,      &
        LOCJS, LP, LRW, LWM, LWMDIM, LWMTEMP, LYH, LYHTEMP, MANPIV, MAPIV, MAXG,&
        MAXIT, MAXORD, MB28, MB48, METH, MICN, MICNCP, MINICN, MINIRN, MIRANK,  &
        MIRN, MIRNCP, MITER, MLP, MOSS, MP, MSBG, MSBJ, MXHNIL, MXSTEP, N, NZB, &
        NCFN, NDROP, NDROP1, NDX, NETF, NEWH, NEWQ, NFE, NGC, NGE, NGP, NHNIL,  &
        NJE, NLP, NLU, NNI, NNZ, NOITER, NQ, NQNYH, NQU, NQWAIT, NSLG, NSLJ,    &
        NSLP, NSRCH, NSRCH1, NST, NSUBS, NSUPS, NUM, NUMNZ, NYH, NZ_ALL,        &
        NZ_SWAG, PREVIOUS_MAXORD, WPD, WPS, MA28AD_CALLS, MA28BD_CALLS,         &
        MA28CD_CALLS, MC19AD_CALLS, MAX_MINIRN, MAX_MINICN, MAX_NNZ, BNGRP
!       MA48AD_CALLS, MA48BD_CALLS, MA48CD_CALLS
! *****MA48 build change point. Insert the above line.
      logical, private :: ABORT, ABORT1, ABORT2, ABORT3, ABORTA, ABORTB,        &
        ALLOW_DEFAULT_TOLS, BUILD_IAJA, BOUNDS, CHANGED_ACOR, GROW, IAJA_CALLED,&
        J_HAS_BEEN_COMPUTED, J_IS_CONSTANT, LBIG, LBIG1, LBLOCK, MA48_WAS_USED, &
        OK_TO_CALL_MA28, SUBS, SUPS, OPTS_CALLED, REDO_PIVOT_SEQUENCE,          &
        SCALE_MATRIX, SPARSE, USE_FAST_FACTOR, YMAXWARN  
! ..
! .. Local Arrays ..
      real (WP), allocatable, private :: ACOR(:), CSCALEX(:), EWT(:),           &
        FPTEMP(:), FTEMP(:), FTEMP1(:), G0(:), G1(:), GX(:), JMAT(:),           &
        LB(:), PMAT(:), RSCALEX(:), RWORK(:), SAVF(:), UB(:), WM(:),            &
        WMTEMP(:), WSCALEX(:,:), YHNQP2(:), YHTEMP(:), YMAX(:), YNNEG(:),       &
        YTEMP(:), DTEMP(:)
      real (WP), private :: EL(13), RUSER(22), TAU(13), TQ(5)
      integer, allocatable, private :: BIGP(:), BJGP(:), IA(:), IAB(:), IAN(:), &
        ICN(:), IDX(:), IGP(:), IKEEP28(:,:), IW28(:,:), IWORK(:), JA(:),       &
        JAB(:), JAN(:), JATEMP(:), JGP(:), JROOT(:), JVECT(:), SUBDS(:), SUPDS(:)
      integer, private :: IDISP(2), IUSER(30), LNPIV(10), LPIV(10)
      integer, private :: MORD(2) = (/ 12, 5 /)
! ..
! .. Public Subroutines and Functions ..
  public ::                                                        &
  DAXPY_F90, DCOPY_F90, DDOT_F90, DGBFA_F90, DGBSL_F90, DGEFA_F90, &
  DGESL_F90, DSCAL_F90, IDAMAX_F90
! ..
! .. Private Subroutines and Functions ..
  private ::                                                       &
  CHECK_DIAG    , DACOPY        , DEWSET        , DGROUP        ,  &
  DGROUPDS      , DVCHECK       , DVHIN         , DVINDY_BNDS   ,  &
  DVINDY_CORE   , DVJAC         , DVJACS28      , DVJUST        ,  &
  DVNLSD        , DVNLSS28      , DVNORM        , DVNRDN        ,  &
  DVNRDP        , DVNRDS        , DVODE         , DVPREPS       ,  &
  DVROOTS       , DVSET         , DVSOL         , DVSOLS28      ,  &
  DVSRCO        , DVSTEP        , GDUMMY        , IUMACH        ,  &
  IXSAV         , JACSPDB       , JDUMMY        , MA28AD        ,  &
  MA28BD        , MA28CD        , MA28DD        , MA28ID        ,  &
  MA30AD        , MA30BD        , MA30CD        , MA30DD        ,  &
  MC13E         , MC19AD        , MC20AD        , MC20BD        ,  &
  MC21A         , MC21B         , MC22AD        , MC23AD        ,  &
  MC24AD        , SET_ICN       , XERRDV        , XSETF         ,  &
  XSETUN        , DEGR          , IDO           , NUMSRT        ,  &
  SEQ           , SETR          , SLO           , SRTDAT        ,  &
  FDJS
! DVJACS48      , DVNLSS48      , DVPREPS48     , DVSOLS48
!_______________________________________________________________________
! *****MA48 build change point. Insert the above line.
!_______________________________________________________________________
! ..
! .. Intrinsic Functions ..
      intrinsic KIND
! ..
! .. Data Statements ..
      data OPTS_CALLED/ .false./
      data MP/6/, NLP/6/, MLP/6/, NSRCH/32768/, ISTART/0/, MAXIT/16/, &
        LBIG/ .false./, LBLOCK/ .true./, GROW/ .true./,               &
        TOL/0.0_WP/, CGCE/0.5_WP/, BIG/0.0_WP/, ABORT1/ .true./,      &
        ABORT2/ .true./, ABORT3/ .false./, ABORT/ .false./, MIRN/0/,  &
        MICN/0/, MIRNCP/0/, MICNCP/0/, MIRANK/0/, NDROP1/0/,          &
        MRMIN/0.0D0/, MRESID/0/, OK_TO_CALL_MA28/.false./
! ..
! END OF DVODE_F90 PRIVATE SECTION.
!_______________________________________________________________________

    contains

      subroutine VODE_F90(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTS,J_FCN,G_FCN)
! ..
! This is an interface for DVODE to allow JAC and GFUN to be
! OPTIONAL arguments.
! ..
     implicit none
! ..
! .. Structure Arguments ..
        type (VODE_OPTS) :: OPTS
! ..
! .. Scalar Arguments ..
        real (WP), intent (INOUT) :: T, TOUT
        integer, intent (INOUT) :: ISTATE
        integer, intent (IN) :: ITASK, NEQ
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: Y(*)
! ..
! .. Subroutine Arguments ..
        optional :: G_FCN, J_FCN
        external J_FCN
! ..
! .. Subroutine Interfaces ..
     interface
       subroutine F(NEQ,T,Y,YDOT)
         integer, parameter :: WP = kind(1.0D0)
         integer NEQ
         real(WP) T
         real(WP), dimension(NEQ) :: Y, YDOT
         intent(IN)  :: NEQ, T, Y
         intent(OUT) :: YDOT
       end subroutine F
     end interface

     interface
       subroutine G_FCN(NEQ,T,Y,NG,GROOT)
         integer, parameter :: WP = kind(1.0D0)
         integer NEQ, NG
         real(WP) T
         real(WP), dimension(NEQ) :: Y
         real(WP), dimension(NG) :: GROOT(NG)
         intent(IN)  :: NEQ, T, Y, NG
         intent(OUT) :: GROOT
       end subroutine G_FCN
     end interface

!    Note:
!    The best we can do here is to declare J_FCN to be
!    external. The interface for a sparse problem differs
!    from that for a banded or dense problem. The following
!    would suffuce for a banded or dense problem.
!    INTERFACE
!      SUBROUTINE J_FCN(NEQ,T,Y,ML,MU,PD,NROWPD)
!        INTEGER, PARAMETER :: WP = KIND(1.0D0)
!        INTEGER NEQ, ML, MU, NROWPD
!        REAL(WP) T
!        REAL(WP), DIMENSION(NEQ) :: Y
!        REAL(WP), DIMENSION(NEQ) :: PD(NROWPD,NEQ)
!        INTENT(IN)  :: NEQ, T, Y, ML, MU, NROWPD
!        INTENT(INOUT) :: PD
!      END SUBROUTINE J_FCN
!    END INTERFACE
!    The following would suffice for a sparse problem.
!    INTERFACE
!      SUBROUTINE J_FCN(NEQ,T,Y,IA,JA,NZ,P)
!        INTEGER, PARAMETER :: WP = KIND(1.0D0)
!        INTEGER NEQ, NZ
!        REAL(WP) T
!        REAL(WP), DIMENSION Y(*), P(*)
!        INTEGER, DIMENSION IA(*), JA(*)
!        INTENT(IN) :: NEQ, T, Y
!        INTENT(INOUT) IA, JA, NZ, P
!      END SUBROUTINE J_FCN
!    END INTERFACE
! ..
! .. Local Scalars ..
        integer :: HOWCALL, METH, MFA, MITER, MOSS, NG
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, PRESENT
! ..
! .. FIRST EXECUTABLE STATEMENT VODE_F90
! ..
!       Check that SET_OPTS has been called.
        if (.not.OPTS_CALLED) then
          MSG = 'You have not called SET_OPTS before'
          call XERRDV(MSG,10,1,0,0,0,0,ZERO,ZERO)
          MSG = 'calling DVODE_F90 the first time.'
          call XERRDV(MSG,10,2,0,0,0,0,ZERO,ZERO)
        end if

!       Check that JAC is present if it is needed.
        if (present(J_FCN)) then
        else
!         Note:
!         MOSS is irrelevant. OPTS%MF is two digits after the
!         call to SET_OPTS.
          MFA = abs(OPTS%MF)
          MOSS = MFA/100
          METH = (MFA-100*MOSS)/10
          MITER = MFA - 100*MOSS - 10*METH
          if (MITER==1 .or. MITER==4 .or. MITER==6) then
            MSG = 'You have specified a value of the integration'
            call XERRDV(MSG,20,1,0,0,0,0,ZERO,ZERO)
            MSG = 'method flag MF which requires that you supply'
            call XERRDV(MSG,20,1,0,0,0,0,ZERO,ZERO)
            MSG = 'a Jacobian subroutine JAC; but FAC is not'
            call XERRDV(MSG,20,1,0,0,0,0,ZERO,ZERO)
            MSG = 'present in the argument list.'
            call XERRDV(MSG,20,2,0,0,0,0,ZERO,ZERO)
          end if
        end if

!       Check that GFUN is present if it is needed.

        if (present(G_FCN)) then
        else
          NG = OPTS%NG
          if (NG>0) then
            MSG = 'You have indicated that events are present but'
            call XERRDV(MSG,30,1,0,0,0,0,ZERO,ZERO)
            MSG = 'you have not supplied a GFUN subroutine.'
            call XERRDV(MSG,30,2,0,0,0,0,ZERO,ZERO)
          end if
        end if

!     Determine how DVODE will be called.

!     HOWCALL = 1: JDUMMY, GDUMMY
!               2: JAC, GFUN
!               3: JAC, GDUMMY
!               4: JDUMMY, GFUN
        HOWCALL = 1
        if (present(J_FCN)) then
          if (present(G_FCN)) then
            HOWCALL = 2
          else
            HOWCALL = 3
          end if
        else
          if (present(G_FCN)) then
            HOWCALL = 4
          else
            HOWCALL = 1
          end if
        end if

!       Call DVODE to do the integration.

        if (HOWCALL==1) then
          call DVODE(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTS,JDUMMY,GDUMMY)
        else if (HOWCALL==2) then
          call DVODE(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTS,J_FCN,G_FCN)
        else if (HOWCALL==3) then
          call DVODE(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTS,J_FCN,GDUMMY)
        else if (HOWCALL==4) then
          call DVODE(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTS,JDUMMY,G_FCN)
        end if
        return

      end subroutine VODE_F90
!_______________________________________________________________________

      subroutine JDUMMY(NEQ,T,Y,ML,MU,PD,NROWPD)
! ..
! This is a dummy Jacobian subroutine for VODE_F90 (never called).
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP) :: T
        integer :: ML, MU, NEQ, NROWPD, I
        logical DUMMY
! ..
! .. Array Arguments ..
        real (WP) :: PD(NROWPD,*), Y(*)
! ..
      intent(IN) T, Y, ML, MU, NROWPD
      intent(INOUT) PD
! ..
! .. FIRST EXECUTABLE STATEMENT JDUMMY
! ..
!       Get rid of some needless compiler warning messages.
        DUMMY = .false.
        if (DUMMY) then
          I = NEQ
          I = ML
          I = MU
          I = NROWPD
          PD(1,1) = T
          PD(1,1) = Y(1)
          PD(1,1) = dble(real(I))
        end if
        return

      end subroutine JDUMMY
!_______________________________________________________________________

      subroutine GDUMMY(NEQ,T,Y,NG,GOUT)
! ..
! This is a dummy event subroutine for VODE_F90 (never called).
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP) :: T
        integer :: NEQ, NG, I
        logical DUMMY
! ..
! .. Array Arguments ..
        real (WP) :: GOUT(*), Y(*)
! ..
      intent(IN) NEQ, T, Y, NG
      intent(OUT) GOUT
! ..
! .. FIRST EXECUTABLE STATEMENT JDUMMY
! ..
!       Get rid of some needless compiler warning messages.
        DUMMY = .false.
        if (DUMMY) then
          I = NEQ
          I = NG
          GOUT(1) = T
          GOUT(1) = Y(1)
          GOUT(1) = dble(real(I))
        end if
        return

      end subroutine GDUMMY
!_______________________________________________________________________

  subroutine SET_OPTS_2(HMAX,HMIN,MXSTEP)
! ..
! Allow the maximum step size, the minimum step size, and the maximum
! number of steps to be changed without restarting the integration.
! ..
!                     Quick Summary of Options
! HMAX                   - Maximum step size in DVODE 
! HMIN                   - Minimum step size in DVODE 
! MXSTEP                 - Maximum number of integration steps in DVODE
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP), optional, intent (IN) :: HMAX, HMIN
        integer, optional, intent (IN) :: MXSTEP
! ..
! .. Local Scalars ..
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
        intrinsic ALLOCATED, PRESENT
! ..
! .. FIRST EXECUTABLE STATEMENT SET_OPTS_2
! ..
!       Check that SET_OPTS has been called:
        if (.not.OPTS_CALLED) then
          MSG = 'You have not called SET_OPTS before'
          call XERRDV(MSG,40,1,0,0,0,0,ZERO,ZERO)
          MSG = 'calling subroutine SET_OPTS_2.'
          call XERRDV(MSG,40,2,0,0,0,0,ZERO,ZERO)
        end if
        if (present(HMAX)) then
          RUSER(6) = HMAX
          MSG = 'HMAX changed in SET_OPTS_2.'
          call XERRDV(MSG,50,1,0,0,0,1,HMAX,ZERO)
        end if
        if (present(HMIN)) then
          RUSER(7) = HMIN
          MSG = 'HMIN changed in SET_OPTS_2.'
          call XERRDV(MSG,60,1,0,0,0,1,HMIN,ZERO)
        end if
        if (present(MXSTEP)) then
          IUSER(6) = MXSTEP
          MSG = 'MXSTEP changed in SET_OPTS_2.'
          call XERRDV(MSG,70,1,1,MXSTEP,0,0,ZERO,ZERO)
        end if

  end subroutine SET_OPTS_2
!_______________________________________________________________________

  function SET_NORMAL_OPTS(DENSE_J, BANDED_J, SPARSE_J,                &
    USER_SUPPLIED_JACOBIAN, LOWER_BANDWIDTH, UPPER_BANDWIDTH,          &
    RELERR, ABSERR, ABSERR_VECTOR, NEVENTS) result(OPTS)

! FUNCTION SET_NORMAL_OPTS:
!    Jacobian type:
!       DENSE_J, BANDED_J, SPARSE_J
!    Analytic Jacobian:
!       USER_SUPPLIED_JACOBIAN
!    If banded Jacobian:
!       LOWER_BANDWIDTH,UPPER_BANDWIDTH
!    Error tolerances:
!       RELERR, ABSERR, ABSERR_VECTOR
!    Rootfinding:
!       NEVENTS
! RESULT(OPTS)

! Note:
! Invoking SET_NORMAL_OPTS causes an integration restart. A common
! situation is one in which all you wish to change is one of the
! vode.f77 optional parameters HMAX, HMIN, or MXSTEP. Once the
! integration is started and this is all you wish to do, you can
! change any of these parameters without restarting the integration
! simply by calling subroutine SET_OPTS_2:
!      CALL SET_OPTS_2(HMAX,HMIN,MXSTEP)
! Each of the three arguments is optional and only the ones actually
! supplied will be used. Changes will take effect in the same manner
! as in the VODE.f77 solver.
!
! NORMAL_OPTIONS sets user parameters for DVODE via keywords.
! Values that are defined herein will be used internally by
! DVODE. All option keywords are OPTIONAL and order is not
! important. These options should be adequate for most problems.
! If you wish to use more specialized options, you must use
! SET_INTERMEDIATE_OPTS or SET_OPTS rather than NORMAL_OPTS.
! If you wish to use SET_INTERMEDIATE_OPTS or SET_OPTS, you
! may use any of the SET_NORMAL_OPTS keywords or any of the
! keywords available for these two functions. Of course, you
! may opt to simply use SET_OPTS for all problems.

! Note that DVODE_F90 requires that one of SET_NORMAL_OPTS or
! SET_INTERMEDIATE_OPTS or SET_OPTS is called before the first
! time DVODE_F90 is called.
!
! Important Note:
! If feasible, you should use the dense or banded option; but
! SET_NORMAL_OPTS allows you to use a sparse internal Jacobian
! (i.e., one that is determined using finite differences) and
! structure pointere array that are determined internally
! using finite differences. If any of the following are true
!    (1) DVODE_F90 doesn't perform satisfactorily for
!        your problem,
!    (2) you are solving a very large problem,
!    (3) you wish to supply the sparse pointer arrays
!        directly,
!    (4) you wish to supply an analytical sparse Jacobian,
! or
!    (5) you wish to use one of the specialized sparse
!        Jacobian options,
! you are encouraged to use SET_OPTS which contains several
! provisions for solving sparse problems more efficiently.

! Option Types:
! DENSE_J                - logical
! BANDED_J               - logical
! SPARSE_J               - logical
! USER_SUPPLIED_JACOBIAN - logical
! LOWER_BANDWIDTH        - integer
! UPPER_BANDWIDTH        - integer
! RELERR                 - real(wp) scalar
! ABSERR                 - real(wp) scalar
! ABSERR_VECTOR          - real(wp) vector
! NEVENTS                - integer
! Options:
! ABSERR                 = Absolute error tolerance
! ABSERR_VECTOR          = Vector of absolute error tolerances
! RELERR                 = Scalar relative error tolerance
! NEVENTS                = Number of event functions (requires
!                          user-supplied GFUN)
! DENSE_J                = Use dense linear algebra if .TRUE.
! BANDED_J               = Use banded linear algebra if .TRUE.
!   LOWER_BANDWIDTH      = Lower bandwidth of the Jacobian
!                          (required if BANDED_J = .TRUE.)
!   UPPER_BANDWIDTH      = Upper bandwidth of the Jacobian
!                          (required if BANDED_J = .TRUE.)
! SPARSE_J               = Use sparse linear algebra if .TRUE.
! USER_SUPPLIED_JACOBIAN = Exact Jacobian option
!                          (requires user-supplied JAC;
!                          ignored for SPARSE_J=.TRUE.)
!
! Note: DENSE_J takes precedence over BANDED_J which in turn
! takes precedence over SPARSE_J if more than one is supplied.
! If neither of the three flags is present, the nonstiff Adams
! option will be used. Similiarly, ABSERR_VECTOR takes
! precedence over ABSERR.
!
! Note on Jacobian Storage Formats:
!
! If you supply an analytic Jacobian PD, load the
! Jacobian elements DF(I)/DY(J), the partial
! derivative of F(I) with respect to Y(J), using
! the following formats. Here, Y is the solution,
! F is the derivative, and PD is the Jacobian.
!
! For a full Jacobian, load PD(I,J) with DF(I)/DY(J).
! Your code might look like this:
!    DO J = 1, NEQ
!       DO I = 1, NEQ
!          PD(I,J) = ... DF(I)/DY(J)
!       END DO
!    END DO       
!
! For a banded Jacobian, load PD(I-J+MU+1,J) with
! DF(I)/DY(J) where ML is the lower bandwidth
! and MU is the upper bandwidth of the Jacobian.
! Your code might look like this:
!    DO J = 1, NEQ
!       I1 = MAX(1,J-ML)
!       I2 = MIN(N,J+MU)
!       DO I = I1, I2
!          K = I-J+MU+1
!          PD(K,J) = ... DF(I)/DY(J)
!       END DO
!    END DO
! ..
     implicit none
! ..
! .. Function Return Value ..
        type (VODE_OPTS) :: OPTS
! ..
! .. Scalar Arguments ..
        real (WP), optional, intent (IN) :: ABSERR,RELERR
        integer, optional, intent (IN) :: LOWER_BANDWIDTH,   &
          NEVENTS,UPPER_BANDWIDTH
        logical, optional, intent (IN) :: BANDED_J, DENSE_J, &
          SPARSE_J, USER_SUPPLIED_JACOBIAN
! ..
! .. Array Arguments ..
        real (WP), optional, intent (IN) :: ABSERR_VECTOR(:)
! ..
! .. Local Scalars ..
        integer :: IER,IOPT,METH,MF,MFA,MFSIGN,MITER,ML,MOSS, &
          MU,NAE,NG,NRE
        logical :: BANDED,DENSE,SPARSE
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
        intrinsic ALLOCATED, IABS, MAX, MINVAL, PRESENT, SIGN, SIZE
! ..
! .. FIRST EXECUTABLE STATEMENT SET_NORMAL_OPTS
! ..
        RUSER(1:LRWUSER) = ZERO
        IUSER(1:LIWUSER) = 0
        
!       Allow default error tolerances?
        ALLOW_DEFAULT_TOLS = .false.

!       Maximum number of consecutive error test failures?
        CONSECUTIVE_EFAILS = KFH

!       Maximum number of consecutive corrector iteration failures?
        CONSECUTIVE_CFAILS = MXNCF

!       Use JACSP to approximate Jacobian?
        USE_JACSP = .false.

!       Set the flag to indicate that SET_NORMAL_OPTS has been called.
        OPTS_CALLED = .true.

!       Set the MA48 storage cleanup flag.
        MA48_WAS_USED = .false.

!       Set the fast factor option for MA48,
        USE_FAST_FACTOR = .true.

!       Set the constant Jacobian flags.
        J_IS_CONSTANT = .false.
        J_HAS_BEEN_COMPUTED = .false.

!       Determine the working precision and define the value for UMAX
!       expected by MA28. Note that it is different for single and
!       double precision.
        WPD = kind(1.0D0)
        WPS = kind(1.0E0)
        if (WPD/=WP .and. WPS/=WP) then
          MSG = 'Illegal working precision in SET_NORMAL_OPTS.'
          call XERRDV(MSG,80,2,0,0,0,0,ZERO,ZERO)
        end if
        if (WPD==WP) then
!       Working precision is double.
          UMAX = 0.999999999_WP
        else
!       Working precision is single.
          UMAX = 0.9999_WP
        end if

        MA28AD_CALLS = 0
        MA28BD_CALLS = 0
        MA28CD_CALLS = 0
        MC19AD_CALLS = 0
!_______________________________________________________________________
! *****MA48 build change point. Insert these statements.
!       MA48AD_CALLS = 0
!       MA48BD_CALLS = 0
!       MA48CD_CALLS = 0
!_______________________________________________________________________
        IRNCP = 0
        ICNCP = 0
        MINIRN = 0
        MINICN = 0
        MAX_MINIRN = 0
        MAX_MINICN = 0
        MAX_NNZ = 0

!       Set the flag to warn the user if |(y(t)| < abserr.
        YMAXWARN = .false.

!       Load defaults for the optional input arrays for DVODE.
        IUSER(1:8) = 0
        RUSER(1:8) = ZERO

!       Set the method flag.
        MF = 10
        if (present(SPARSE_J)) then
          if (SPARSE_J) then
             MF = 227
             if (present(USER_SUPPLIED_JACOBIAN)) then
               MSG = 'You have indicated you wish to supply an'
               call XERRDV(MSG,90,1,0,0,0,0,ZERO,ZERO)
               MSG = 'exact sparse Jacobian in function'
               call XERRDV(MSG,90,1,0,0,0,0,ZERO,ZERO)
               MSG = 'SET_NORMAL_OPTS. In order to do this,'
               call XERRDV(MSG,90,1,0,0,0,0,ZERO,ZERO)
               MSG = 'you must use SET_OPTS. Execution will'
               call XERRDV(MSG,90,1,0,0,0,0,ZERO,ZERO)
               MSG = 'continue.'
               call XERRDV(MSG,90,1,0,0,0,0,ZERO,ZERO)
             end if
          end if
        end if

        if (present(BANDED_J)) then
          if (BANDED_J) then
            if (present(USER_SUPPLIED_JACOBIAN)) then
              if (USER_SUPPLIED_JACOBIAN) then
                MF = 24
              else
                MF = 25
              end if
            else
              MF = 25
            end if
          end if
        end if

        if (present(DENSE_J)) then
          if (DENSE_J) then
            if (present(USER_SUPPLIED_JACOBIAN)) then
              if (USER_SUPPLIED_JACOBIAN) then
                MF = 21
              else
                MF = 22
              end if
            else
              MF = 22
            end if
          end if
        end if

!       Check for errors in MF.
        MFA = IABS(MF)
        MOSS = MFA/100
        METH = (MFA-100*MOSS)/10
        MITER = MFA - 100*MOSS - 10*METH
        if (METH<1 .or. METH>2) then
          MSG = 'Illegal value of METH in SET_NORMAL_OPTS.'
          call XERRDV(MSG,100,2,0,0,0,0,ZERO,ZERO)
        end if
        if (MITER<0 .or. MITER>7) then
          MSG = 'Illegal value of MITER in SET_NORMAL_OPTS.'
          call XERRDV(MSG,110,2,0,0,0,0,ZERO,ZERO)
        end if
        if (MOSS<0 .or. MOSS>2) then
          MSG = 'Illegal value of MOSS in SET_NORMAL_OPTS.'
          call XERRDV(MSG,120,2,0,0,0,0,ZERO,ZERO)
        end if

!       Reset MF, now that MOSS is known.
        MFSIGN = sign(1,MF)
        MF = MF - 100*MOSS*MFSIGN

        if (MITER==0) then
          DENSE = .false.
          BANDED = .false.
          SPARSE = .false.
        else if (MITER==1 .or. MITER==2) then
          DENSE = .true.
          BANDED = .false.
          SPARSE = .false.
        else if (MITER==3) then
          DENSE = .false.
          BANDED = .false.
          SPARSE = .false.
        else if (MITER==4 .or. MITER==5) then
          DENSE = .false.
          BANDED = .true.
          SPARSE = .false.
        else if (MITER==6 .or. MITER==7) then
          DENSE = .false.
          BANDED = .false.
          SPARSE = .true.
        end if

!       Define the banded Jacobian band widths.
        if (BANDED) then
          if (present(LOWER_BANDWIDTH)) then
            ML = LOWER_BANDWIDTH
            IUSER(1) = ML
          else
            MSG = 'In SET_NORMAL_OPTS you have indicated a'
            call XERRDV(MSG,130,1,0,0,0,0,ZERO,ZERO)
            MSG = 'banded Jacobian but you have not supplied'
            call XERRDV(MSG,130,1,0,0,0,0,ZERO,ZERO)
            MSG = 'the lower bandwidth.'
            call XERRDV(MSG,130,2,0,0,0,0,ZERO,ZERO)
          end if
          if (present(UPPER_BANDWIDTH)) then
            MU = UPPER_BANDWIDTH
            IUSER(2) = MU
          else
            MSG = 'In SET_NORMAL_OPTS you have indicated a'
            call XERRDV(MSG,140,1,0,0,0,0,ZERO,ZERO)
            MSG = 'banded Jacobian but you have not supplied'
            call XERRDV(MSG,140,1,0,0,0,0,ZERO,ZERO)
            MSG = 'the upper bandwidth.'
            call XERRDV(MSG,140,2,0,0,0,0,ZERO,ZERO)
          end if
        end if

!       Define the sparse Jacobian options.
        if (SPARSE) then
!         Set the MA28 message flag.
          LP = 0
!         Set the MA28 pivot sequence frequency flag.
          REDO_PIVOT_SEQUENCE = .false.
!         Set the MA28 singularity threshold.
          EPS = 1.0E-4_WP
!         Use scaling of the iteration matrix.
          SCALE_MATRIX = .true.
!         Define the elbow room factor for the MA28 sparse work arrays.
          ELBOW_ROOM = 2
!         NZSWAG is a swag for the number of nonzeros in the Jacobian.
          NZ_SWAG = 0
!         Use partial pivoting.
          U_PIVOT = ONE
!         Indicate that SET_IAJA has not yet been called successfully.
          IAJA_CALLED = .false.
!         Check for illegal method flags.
          if (MOSS==2 .and. MITER/=7) then
            MSG = 'In SET_NORMAL_OPTS MOSS=2 but MITER is not 7.'
            call XERRDV(MSG,150,2,0,0,0,0,ZERO,ZERO)
          end if
          if (MOSS==1 .and. MITER/=6) then
            MSG = 'In SET_NORMAL_OPTS MOSS=1 but MITER is not 6.'
            call XERRDV(MSG,160,2,0,0,0,0,ZERO,ZERO)
          end if
!         IF (MOSS==0 .AND. MITER/=7) THEN
!           MSG = 'In SET_NORMAL_OPTS MOSS=0 but MITER is not 7.'
!           CALL XERRDV(MSG,170,2,0,0,0,0,ZERO,ZERO)
!         END IF
        end if

!       Define the number of event functions.
        if (present(NEVENTS)) then
          if (NEVENTS>0) then
            NG = NEVENTS
          else
            NG = 0
          end if
        else
          NG = 0
        end if

!       No solution bounds will be imposed.
        BOUNDS = .false.

!       Load the user options into the solution structure.
        OPTS%MF = MF
        OPTS%METH = METH
        OPTS%MITER = MITER
        OPTS%MOSS = MOSS
        OPTS%DENSE = DENSE
        OPTS%BANDED = BANDED
        OPTS%SPARSE = SPARSE
        OPTS%NG = NG

        IOPT = 1
        OPTS%IOPT = IOPT

!       Process the error tolerances.

!       Relative error tolerance.
        NRE = 1
        allocate (OPTS%RTOL(NRE),STAT=IER)
        call CHECK_STAT(IER,10)
        if (present(RELERR)) then
          if (RELERR<ZERO) then
            MSG = 'RELERR must be nonnegative.'
            call XERRDV(MSG,180,2,0,0,0,0,ZERO,ZERO)
          end if
          OPTS%RTOL = RELERR
        else
          if (ALLOW_DEFAULT_TOLS) then
             OPTS%RTOL = 1.0E-4_WP
             MSG = 'By not specifying RELERR, you have elected to use a default'
             call XERRDV(MSG,190,1,0,0,0,0,ZERO,ZERO)
             MSG = 'relative error tolerance equal to 1.0D-4. Please be aware a'
             call XERRDV(MSG,190,1,0,0,0,0,ZERO,ZERO)
             MSG = 'tolerance this large is not appropriate for all problems.'
             call XERRDV(MSG,190,1,0,0,0,0,ZERO,ZERO)
             MSG = 'Execution will continue'
             call XERRDV(MSG,190,1,0,0,0,0,ZERO,ZERO)
          else
             MSG = 'You must specify a nonzero relative error tolerance.'
             call XERRDV(MSG,200,2,0,0,0,0,ZERO,ZERO)
          end if
        end if

!       Absolute error tolerance(s).
        if (present(ABSERR_VECTOR)) then
          if (minval(ABSERR_VECTOR)<ZERO) then
            MSG = 'All components of ABSERR_VECTOR must'
            call XERRDV(MSG,210,1,0,0,0,0,ZERO,ZERO)
            MSG = 'be nonnegative.'
            call XERRDV(MSG,210,2,0,0,0,0,ZERO,ZERO)
          end if
          NAE = size(ABSERR_VECTOR)
        else
          NAE = 1
        end if
        allocate (OPTS%ATOL(NAE),STAT=IER)
        call CHECK_STAT(IER,20)
        if (present(ABSERR_VECTOR)) then
          OPTS%ATOL = ABSERR_VECTOR
        else if (present(ABSERR)) then
          if (ABSERR<ZERO) then
            MSG = 'ABSERR must be nonnegative.'
            call XERRDV(MSG,220,2,0,0,0,0,ZERO,ZERO)
          end if
          OPTS%ATOL = ABSERR
        else
          if (ALLOW_DEFAULT_TOLS) then
             OPTS%ATOL = 1D-6
             MSG = 'By not specifying ABSERR, you have elected to use a default'
             call XERRDV(MSG,230,1,0,0,0,0,ZERO,ZERO)
             MSG = 'absolute error tolerance equal to 1.0D-6. Please be aware a'
             call XERRDV(MSG,230,1,0,0,0,0,ZERO,ZERO)
             MSG = 'tolerance this large is not appropriate for all problems.'
             call XERRDV(MSG,230,1,0,0,0,0,ZERO,ZERO)
             MSG = 'Execution will continue'
             call XERRDV(MSG,230,1,0,0,0,0,ZERO,ZERO)
          else
             MSG = 'You must specify a vector of absolute error tolerances or'
             call XERRDV(MSG,240,1,0,0,0,0,ZERO,ZERO)
             MSG = 'a scalar error tolerance. It is recommended that you use'
             call XERRDV(MSG,240,1,0,0,0,0,ZERO,ZERO)
             MSG = 'a vector of absolute error tolerances.'
             call XERRDV(MSG,240,2,0,0,0,0,ZERO,ZERO)
          end if
        end if

!       ITOL error tolerance flag.
!          ITOL   RTOL     ATOL            EWT(i)
!            1   scalar   scalar  RTOL*ABS(Y(i)) + ATOL
!            2   scalar   array   RTOL*ABS(Y(i)) + ATOL(i)
        if (present(ABSERR_VECTOR)) then
           OPTS%ITOL = 2
        else
           OPTS%ITOL = 1
        end if
        return

  end function SET_NORMAL_OPTS
!_______________________________________________________________________

  function SET_INTERMEDIATE_OPTS(DENSE_J, BANDED_J, SPARSE_J,          &
    USER_SUPPLIED_JACOBIAN,                                            &
    LOWER_BANDWIDTH, UPPER_BANDWIDTH,                                  &
    RELERR, ABSERR, ABSERR_VECTOR,                                     &
    TCRIT, H0, HMAX, HMIN, MAXORD, MXSTEP, MXHNIL,                     &
    NZSWAG, USER_SUPPLIED_SPARSITY, MA28_RPS,                          &
    NEVENTS, CONSTRAINED, CLOWER, CUPPER, CHANGE_ONLY_f77_OPTIONS)     &
  result(OPTS)

! FUNCTION SET_INTERMEDIATE_OPTS:
!    Jacobian type:
!       DENSE_J, BANDED_J, SPARSE_J
!    Analytic Jacobian:
!       USER_SUPPLIED_JACOBIAN
!    If banded Jacobian:
!       LOWER_BANDWIDTH, UPPER_BANDWIDTH
!    Error tolerances:
!       RELERR, ABSERR, ABSERR_VECTOR
!    VODE.f77 optional parameters:
!       TCRIT, H0, HMAX, HMIN, MAXORD, MXSTEP, MXHNIL
!    Sparse flags:
!       NZSWAG, USER_SUPPLIED_SPARSITY, MA28_RPS
!    Rootfinding:
!       NEVENTS
!    Impose bounds on solution:
!       CONSTRAINED, CLOWER, CUPPER
!    Change one or more of HMAX, HMIN, TCRIT, MXSTEP, MXHNIL, MAXORD
!       CHANGE_ONLY_f77_OPTIONS
! RESULT(OPTS)

! SET_OPTIONS sets user parameters for DVODE via keywords. Values
! that are defined herein will be used internally by DVODE.
! All option keywords are OPTIONAL and order is not important.

! Note that DVODE_F90 requires that one of SET_NORMAL_OPTS or
! SET_INTERMEDIATE_OPTS or SET_OPTS is called before the first
! time DVODE_F90 is called.

!                     Quick Summary of Options

! DENSE_J                - Jacobian is sparse alternative to MF
! BANDED_J               - Jacobian is banded alternative to MF 
! SPARSE_J               - Jacobian is sparse alternative to MF   
! USER_SUPPLIED_JACOBIAN - User will supply Jacobian subroutine
!                          (user supplied subroutine JAC required)
! LOWER_BANDWIDTH        - Lower bandwidth ML in DVODE
! UPPER_BANDWIDTH        - Upper bandwidth MU in DVODE
! RELERR                 - Scalar relative error tolerance in DVODE
! ABSERR                 - Scalar absolute error tolerance in DVODE
! ABSERR_VECTOR          - Vector absolute error tolerance in DVODE
! TCRIT                  - Critical time TCRIT in DVODE
! H0                     - Starting step size in DVODE
! HMAX                   - Maximum step size in DVODE 
! HMIN                   - Minimum step size in DVODE 
! MAXORD                 - Maximum integration order in DVODE
! MXSTEP                 - Maximum number of integration steps
!                          in DVODE
! MXHNIL                 - Maximum number of T+H=T messages in DVODE
! NZSWAG                 - guess for the number of nonzeros in sparse
!                          Jacobian
! USER_SUPPLIED_SPARSITY - user will supply sparsity structure
!                          arrays by calling USERSETS_IAJA 
! MA28_RPS               - Redo MA28AD pivot sequence if a singularity
!                          is encountered
! NEVENTS                - number of user defined root finding functions
!                          (user supplied subroutine G required)
! CONSTRAINED,           - array of solution component indices that
! CLOWER,                  are to be constrained by the lower and
! CUPPER                   upper bound arrays CLOWER and CUPPER so that
!                          CLOWER(I) <= Y(CONSTRAINED(I)) <= CUPPER(I)
!                     Options Types
! DENSE_J                - logical
! BANDED_J               - logical
! SPARSE_J               - logical
! USER_SUPPLIED_JACOBIAN - logical
! LOWER_BANDWIDTH        - integer
! UPPER_BANDWIDTH        - integer
! RELERR                 - real(wp) scalar
! ABSERR                 - real(wp) scalar
! ABSERR_VECTOR          - real(wp) vector
! TCRIT                  - real(wp) scalar
! H0                     - real(wp) scalar
! HMAX                   - real(wp) scalar
! HMIN                   - real(wp) scalar
! MAXORD                 - integer
! MXSTEP                 - integer
! MXHNIL                 - integer
! NZSWAG                 - integer
! USER_SUPPLIED_SPARSITY - logical
! MA28_RPS               - logical
! NEVENTS                - integer
! CONSTRAINED            - integer array
! CLOWER                 - real(wp) array
! CUPPER                 - real(wp) array
! CHANGE_ONLY_f77_OPTIONS- logical
! Argument list parameters:
! ABSERR                   = scalar absolute error tolerance
! ABSERR_VECTOR            = vector of absolute error tolerances
! RELERR                   = scalar relative error tolerance
! NEVENTS                  = Number of event functions (requires
!                            user-supplied GFUN)
! DENSE_J                  = use dense linear algebra if .TRUE.
! BANDED_J                 = use banded linear algebra if .TRUE.
!   LOWER_BANDWIDTH        = lower bandwidth of the Jacobian
!                            (required if BANDED_J = .TRUE.)
!   UPPER_BANDWIDTH        = upper bandwidth of the Jacobian
!                            (required if BANDED_J = .TRUE.)
! SPARSE_J                 = use sparse linear algebra if .TRUE.
!
!   NZSWAG                 = If you wish you may supply a guess,
!                            NZSWAG, at the number of nonzeros
!                            in the Jacobian matrix. In some cases
!                            this will speed up the determination
!                            of the necessary storage.
!   USER_SUPPLIED_SPARSITY = .TRUE. if you wish to supply the sparsity
!                            structure arrays directly by calling
!   MA28_RPS               = .TRUE. to force MA28AD to calculate a new
!                            pivot sequence for use in MA28BD if a
!                            singular iteration matrix is encountered
!                            is smaller than EPS, it will flag the
!                            Jacobian as singular and force MA28AD
!                            to calculate a new pivot sequence.
! USER_SUPPLIED_JACOBIAN   = exact Jacobian option
!                            (requires user-supplied JAC)
! TCRIT                    = critical time
! H0                       = initial step size to try
! HMAX                     = maximum allowable step size
! HMIN                     = minimum allowable step size
! MAXORD                   = maximum allowable integration order
! MXSTEP                   = maximum number of integration steps
!                            between calls to DVODE
! MXHNIL                   = maximum number of printed messages
!                            if the T+H=T condition occurs
! CONSTRAINED              = array containing the indices of
!                            solution components which are to be
!                            be constrained by CLOWER and CUPPER.
!                            The size of CONSTRAINED must be
!                            positive and not exceed NEQ, the
!                            number of ODEs. Each component of
!                            CONSTRAINED must be between 1 and NEQ.
!   CLOWER, CUPPER         = lower and upper bound arrays for the.
!                            solution. Each must be the same size
!                            as the CONSTRAINED vector.
!
!   User-callable Routines:
!   The following routines may be called by the user:
!   SET_INTERMEDIATE_OPTS      : Used to set options for DVODE_F90
!   GET_STATS     : Used to gather summary integration
!                   statistics following a successful
!                   return from DVODE_F90
!   DVINDY        : Used if the user wishes to interpolate
!                   solution or derivative following a
!                   successful return from DVODE_F90
!   USERSETS_IAJA : Used if the user wishes to supply the
!                   sparsity structure directly
!   Detailed Description of SET_INTERMEDIATE_OPTS
!   The following are defined in SET_INTERMEDIATE_OPTS:
!   OPTS%MF     = Integration method flag (MF)
!   OPTS%METH   = Integration family flag (METH)
!   OPTS%MITER  = Iteration type flag (MITER)
!   OPTS%MOSS   = Sparsity array type flag (MOSS)
!   OPTS%ITOL   = Error tolerance flag (ITOL) (ITOL)
!   OPTS%ATOL   = Absolute error tolerance(s) (ATOL)
!   OPTS%RTOL   = Relative error tolerance(s) (RTOL)
!   OPTS%DENSE  = Use dense linear algebra
!   OPTS%BANDED = Use banded linear algebra 
!   OPTS%SPARSE = Use sparse linear algebra 
!   OPTS%IOPT   = DVODE optional parameter input flag (IOPT)
!   OPTS%NG     = Number of event functions
!   RUSER(1)    = TCRIT (don't step past)
!   IUSER(1)    = Jacobian lower bandwidth (ML)
!   IUSER(2)    = Jacobian upper bandwidth (MU)
!   RUSER(5)    = Initial step size to try (H0)
!   RUSER(6)    = Maximum allowable step size (HMAX)
!   RUSER(7)    = Minimum allowable step size (HMIN)
!   IUSER(5)    = Maximum allowable integration order (MAXORD)
!   IUSER(6)    = Maximum number of integration steps (MXSTEP)
!   IUSER(7)    = Maximum number of T+H=T messages (MXHNIL)
!   NZ_SWAG     = Guess for the number of nonzeros in the
!                 sparse Jacobian
!   NG          = Number of user event functions
!   BOUNDS,     = Nonnegativity information
!   NDX, IDX
!   YMAXWARN    = Warning flag for |y(t)| < abserr
!   MA28_ELBOW_ROOM = Integer multiple by which to increase
!                     the elbow room in the MA28 sparse work
!                     arrays
!   MC19_SCALING    = logical flag to control MC19 sparse
!                     scaling of the Jacobian
!   MA28_MESSAGES   = logical flag to control printing of MA28
!                     messages.
!   MA28_EPS        = real(wp) MA28 singularity threshold
!
!                         All Options
!
! METHOD_FLAG
!
! Default:   Not used
! Change to: No need to specify if use one of next three parameters
!            but can be changed to any value of the MF method flag
!            in dvode.f77
!
! DENSE_J
!
! Default:   .FALSE.
! Change to: .TRUE. for dense linear algebra
!
! BANDED_J
!
! Default:   .FALSE.
! Change to: .TRUE. for banded linear algebra
!
! SPARSE_J
!
! Default:   .FALSE.
! Change to: .TRUE. for sparse linear algebra

! USER_SUPPLIED_JACOBIAN
!
! Default:   .FALSE.
! Change to: .TRUE. if want to supply subroutine JAC to DVODE_F90
!
! LOWER_BANDWIDTH
!
! Default:   Not used
! Change to: Lower bandwidth if BANDED_J = .TRUE.
!
! UPPER_BANDWIDTH
!
! Default:   Not used
! Change to: Upper bandwidth if BANDED_J = .TRUE.
!
! RELERR
!
! Default:   None
! Change to: Specify a nonzero relative error tolerance
!
! ABSERR
!
! Default:   None
! Change to: Specify a scalar absolute error tolerance or
!            a vector of absolute error tolerances
!
! ABSERR_VECTOR
!
! Default:   Not used
! Change to: Vector of absolute error tolerances
!
! TCRIT
!
! Default:   Not used
! Change to: Value of T which DVODE_F90 is not to step past
!
!     H0
!
! Default:   Determined by DVODE_F90
! Change to: Guess for initial step size
!
! HMAX
!
! Default:   Infinity
! Change to: Maximum allowable step size
!
! HMIN
!
! Default:   0.0D0
! Change to: Minimum allowable step size
!
! MAXORD
!
! Default:   Determined by DVODE_F90
! Change to: Maximum integration order
!
! MXSTEP
!
! Default:   5000
!     Change to: Maximum number of steps between output points
!
! MXHNIL
!
! Default:   10
! Change to: Maximum number of times T+H=T message will be printed
!
! NZSWAG
!
! Default:   Determined by DVODE_F90
! Change to: Amount by which to increment sizes of sparse arrays
!
! USER_SUPPLIED_SPARSITY
!
! Default:   .FALSE.
! Change to: .TRUE. if wish to call subroutine USERSETS_IAJA to
!            define the sparse structure array pointers directly
!
! NEVENTS
!
! Default:   0
! Change to: Number of event functions if wish to supply subroutine
!            GFUN to DVODE_F90
!
! CONSTRAINED
!
! Default:   Bounds not imposed on solution
! Change to: Array of indices for components on which to impose
!            solution bounds
!
! CLOWER
!
! Default:   Not used
! Change to: Vector containing lower bounds corresponding to
!            CONSTRAINED
!
! CUPPER
!
! Default:   Not used
! Change to: Vector containing upper bounds corresponding to
!            CONSTRAINED
!
! MA28_RPS
!
! Default:   .FALSE.
! Change to: Redo the MA28AD sparse pivoting sequence any time
!            MA28BD considers the iteration matrix to be
!            numerically singularity
!
!
! Note on Jacobian Storage Formats:
!
! If you supply an analytic Jacobian PD, load the
! Jacobian elements DF(I)/DY(J), the partial
! derivative of F(I) with respect to Y(J), using
! the following formats. Here, Y is the solution,
! F is the derivative, and PD is the Jacobian.
!
! For a full Jacobian, load PD(I,J) with DF(I)/DY(J).
! Your code might look like this:
!    DO J = 1, NEQ
!       DO I = 1, NEQ
!          PD(I,J) = ... DF(I)/DY(J)
!       END DO
!    END DO       
!
! For a banded Jacobian, load PD(I-J+MU+1,J) with
! DF(I)/DY(J) where ML is the lower bandwidth
! and MU is the upper bandwidth of the Jacobian.
! Your code might look like this:
!    DO J = 1, NEQ
!       I1 = MAX(1,J-ML)
!       I2 = MIN(N,J+MU)
!       DO I = I1, I2
!          K = I-J+MU+1
!          PD(K,J) = ... DF(I)/DY(J)
!       END DO
!    END DO
!
! For a sparse Jacobian, IA(J+1)-IA(J) is the number
! of nonzeros in column change. JA(I) indicates the
! rows in which the nonzeros occur. For column J,
! the nonzeros occurs in rows I=JA(K) for K=I1,...,I2
! where I1=IA(J) and I2= IA(J+1)-1. Load DF(I)/DY(J)
! in PD(I). Your code might look like this:
!    DO J = 1, NEQ
!       I1 = IA(J)
!       I2 = IA(J+1) - 1
!       DO K = I1, I2
!          I = JA(K)
!          PD(I) = ... DF(I)/DY(J)
!       END DO
!    END DO
!
!                    More on Sparsity Options
!
! Two facts of life regarding the use of direct sparse solvers are
! (1) significant improvements are possible, and (2) the use of
! direct sparse solvers often is more demanding of the user.
! Although SET_NORMAL_OPTS provides modest provisions for solving
! problems with sparse Jacobians, using SET_INTERMEDIATE_OPTS rather than
! SET_NORMAL_OPTS provides several advanced options. These options
! are described below. The recommended manner in which to use
! these options is also provided. Note that each of these optional
! parameters have default values and need not be specified in your
! call to SET_INTERMEDIATE_OPTS if you wish to use the default values. Note
! also that the order in which the options are specified in a call
! to SET_INTERMEDIATE_OPTS is not important.
!
!  (1) First determine if it is feasible to use either the BANDED
!      or DENSE Jacobian option.
!         Recommendation:
!         Use the DENSE or BANDED option if possible. They do
!         not require use of most of the options described here.
!    
!  (2) The Jacobian is approximated internally using differences
!      if you do not provide an analytic Jacobian. The option
!      USER_SUPPLIED_JACOBIAN=.TRUE. may be used if you wish to
!      provide an analytic Jacobian.
!         Recommendation:
!         Use an internally generated Jacobian for most problems
!         but consider providing an analytic Jacobian if it not
!         too much trouble.
!
!  (3) If you do not provide the sparse structure arrays, they
!      are approximated internally by making NEQ calls to your
!      derivative subroutine and using differences to approximate
!      the structure of the Jacobian. The option
!      USER_SUPPLIED_SPARSITY and a call to USERSETS_IAJA can
!      be used to supply the arrays directly. You can also use
!      subroutine SET_IAJA to approximate the structure using
!      different perturbation factors than those used in DVODE_F90.
!         Recommendations:
!         Although allowing DVODE_F90 to approximate the Jacobian
!         structure suffices for most problems, if you know the
!         sparsity pattern provide it directly. This eliminates
!         the possibility that an important element that happens
!         to be 0 when the sparsity pattern is approximated is
!         later nonzero; and it avoids the NEQ extra calls to
!         your derivative subroutine to approximate the structure.
!         Note that a nemesis for any sparse ode solver is a
!         problem in which the sparsity pattern actually changes
!         during the integration. An example of such a problem
!         is provided in the demohum21.f90 demonstration program.
!         If you know the sparsity pattern changes or if you
!         suspect it does because DVODE_F90 is generating
!         nonconvergence error messages, consider having DVODE_F90
!         re-approximate the structure by calling SET_INTERMEDIATE_OPTS and
!         forcing an integration restart when control is returned
!         to your calling program. Please do not do this at every
!         output point since it will be extremely time consuming
!         and inefficient if it is not needed.
!
!  (4) The optional parameter NZSWAG may/should be used to speed
!      up the determination of acceptable array sizes for the
!      internal sparse working arrays and the sparse Jacobian.
!      Initially, DVODE_F90 allocates arrays of length
!      max(10*NEQ,NSZSWAG) and increases this amount by
!      max(10*NEQ,ELBOW_ROOM*NSZSWAG) as necessary. NZSWAG
!      has a default value of 1000 and ELBOW_ROOM has a default
!      value of 2.
!         Recommendation:
!         Provide a larger value for NZSWAG particularly if NEQ
!         is large and you suspect considerable fill-in due to
!         partial pivoting.
!
!  (5) MA28BD uses the pivot sequence initially determined by
!      MA28AD. When a singularity is diagnosed, the internal
!      procedure used by the DENSE and BANDED options is used
!      to reduce the step size and re-calulate the iteration
!      matrix. This works fine for most problems; but particularly
!      for badly scaled problems, the solution may be unsuccessful
!      or it may drag along using small step sizes. This is due
!      to the fact MA28BD will continue to use the out-of-date pivot
!      sequence until singularity is again diagnosed. An option not
!      available in previous sparse ode solvers may be used to
!      instruct DVODE_F90 to force MA28AD to calculate a new pivot
!      sequence when MA28BD encounters a singularity. Although
!      MA28AD is considerably slower than MA28BD, the additional
!      calls to MA28AD ensure that the pivot sequence is more
!      up to date. This can increase both the accuracy and the
!      efficiency very dramatically. The demodirn.f90 demonstration
!      program provides an illustration of the dramatic improvement
!      that is possible. The optional parameter MA28_RPS may be set
!      .TRUE. to force these pivot sequence updates.
!         Recommendation:
!         Use the default value MA28_RPS=.FALSE; but if DVODE_F90
!         encounters nonconvergence, use MA28_RPS=.TRUE. to force
!         pivot sequence updates.
! ..
     implicit none
! ..
! .. Function Return Value ..
        type (VODE_OPTS) :: OPTS
! ..
! .. Scalar Arguments ..
        real (WP), optional, intent (IN) :: ABSERR, H0, HMAX, HMIN,         &
        RELERR, TCRIT
        integer, optional, intent (IN) :: LOWER_BANDWIDTH, MAXORD,          &
          MXHNIL, MXSTEP, NEVENTS, NZSWAG, UPPER_BANDWIDTH
        logical, optional, intent (IN) :: BANDED_J, CHANGE_ONLY_f77_OPTIONS,&
          DENSE_J, MA28_RPS, SPARSE_J, USER_SUPPLIED_JACOBIAN,              &
          USER_SUPPLIED_SPARSITY
! ..
! .. Array Arguments ..
        real (WP), optional, intent (IN) :: ABSERR_VECTOR(:)
        real (WP), optional :: CLOWER(:), CUPPER(:)
        integer, optional, intent (IN) :: CONSTRAINED(:)
! ..
! .. Local Scalars ..
        integer ::  IER, IOPT, METH, MF, MFA, MFSIGN, MITER, ML, MOSS, MU,  &
          NAE, NG, NRE
        logical :: BANDED, DENSE, SPARSE
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
        intrinsic ALLOCATED, IABS, MAX, MINVAL, PRESENT, SIGN, SIZE
! ..
! .. FIRST EXECUTABLE STATEMENT SET_INTERMEDIATE_OPTS
! ..
    RUSER(1:LRWUSER) = ZERO
    IUSER(1:LIWUSER) = 0
    
!   Allow default error tolerances?
    ALLOW_DEFAULT_TOLS = .false.

!   Maximum number of consecutive error test failures?
    CONSECUTIVE_EFAILS = KFH

!   Maximum number of consecutive corrector iteration failures?
    CONSECUTIVE_CFAILS = MXNCF

!   Use JACSP to approximate Jacobian?
    USE_JACSP = .false.

!   If only f77 options are to be changed, do it and return.
    if (present(CHANGE_ONLY_f77_OPTIONS)) then
       if (CHANGE_ONLY_f77_OPTIONS) then
          if (.not.OPTS_CALLED) then
             MSG = 'You have not previously called SET_OPTS before attempting'
             call XERRDV(MSG,250,1,0,0,0,0,ZERO,ZERO)
             MSG = 'to change one or more of the vode.f77 optional parameters.'
             call XERRDV(MSG,250,2,0,0,0,0,ZERO,ZERO)
          end if
          if (present(HMAX)) then
             IOPT = 1
             RUSER(6) = HMAX
             MSG = 'HMAX changed in SET_INTERMEDIATE_OPTS.'
             call XERRDV(MSG,260,1,0,0,0,1,HMAX,ZERO)
          end if
          if (present(HMIN)) then
             IOPT = 1
             RUSER(7) = HMIN
             MSG = 'HMIN changed in SET_INTERMEDIATE_OPTS.'
             call XERRDV(MSG,270,1,0,0,0,1,HMIN,ZERO)
          end if
          if (present(TCRIT)) then
             IOPT = 1
             RUSER(1) = TCRIT
             MSG = 'TCRIT changed in SET_INTERMEDIATE_OPTS.'
             call XERRDV(MSG,280,1,0,0,0,1,TCRIT,ZERO)
          end if
          if (present(MXSTEP)) then
             IOPT = 1
             IUSER(6) = MXSTEP
             MSG = 'MXSTEP changed in SET_INTERMEDIATE_OPTS.'
             call XERRDV(MSG,290,1,1,MXSTEP,0,0,ZERO,ZERO)
          end if
          if (present(MAXORD)) then
             IOPT = 1
             IUSER(5) = MAXORD
             MSG = 'MAXORD changed in SET_INTERMEDIATE_OPTS.'
             call XERRDV(MSG,300,1,1,MAXORD,0,0,ZERO,ZERO)
          end if
          if (present(MXHNIL)) then
             IOPT = 1
             IUSER(7) = MXHNIL
             MSG = 'MXHNIL changed in SET_INTERMEDIATE_OPTS.'
             call XERRDV(MSG,310,1,1,MXHNIL,0,0,ZERO,ZERO)
          end if
       end if
       return
    end if

!       Set the flag to indicate that SET_INTERMEDIATE_OPTS has been called.
        OPTS_CALLED = .true.

!       Set the MA48 storage cleanup flag.
        MA48_WAS_USED = .false.

!       Set the fast factor option for MA48,
        USE_FAST_FACTOR = .true.

!       Determine the working precision and define the value for UMAX
!       expected by MA28. Note that it is different for single and
!       double precision.
        WPD = kind(1.0D0)
        WPS = kind(1.0E0)
        if (WPD/=WP .and. WPS/=WP) then
          MSG = 'Illegal working precision in SET_INTERMEDIATE_OPTS.'
          call XERRDV(MSG,320,2,0,0,0,0,ZERO,ZERO)
        end if
        if (WPD==WP) then
!       Working precision is double.
          UMAX = 0.999999999_WP
        else
!       Working precision is single.
          UMAX = 0.9999_WP
        end if

!       Set the MA28 message flag.
        LP = 0

!       Set the MA28 singularity threshold.
        EPS = 1.0E-4_WP

!       Set the MA28 pivot sequence frequency flag.
        if (present(MA28_RPS)) then
           if (MA28_RPS) then
              REDO_PIVOT_SEQUENCE = MA28_RPS
           else
              REDO_PIVOT_SEQUENCE = .false.
           end if
        else
           REDO_PIVOT_SEQUENCE = .false.
        end if

        MA28AD_CALLS = 0
        MA28BD_CALLS = 0
        MA28CD_CALLS = 0
        MC19AD_CALLS = 0
!_______________________________________________________________________
! *****MA48 build change point. Insert these statements.
!       MA48AD_CALLS = 0
!       MA48BD_CALLS = 0
!       MA48CD_CALLS = 0
!_______________________________________________________________________
        IRNCP = 0
        ICNCP = 0
        MINIRN = 0
        MINICN = 0
        MAX_MINIRN = 0
        MAX_MINICN = 0
        MAX_NNZ = 0

!       Set the flag to warn the user if |(y(t)| < abserr.
        YMAXWARN = .false.

!       Load defaults for the optional input arrays for DVODE.
        IUSER(1:8) = 0
        RUSER(1:8) = ZERO

        if (present(SPARSE_J)) then
          if (SPARSE_J) then
            if (present(USER_SUPPLIED_JACOBIAN)) then
              if (USER_SUPPLIED_JACOBIAN) then
                MF = 126
              else
                MF = 227
              end if
            else
              MF = 227
            end if
            if (present(USER_SUPPLIED_SPARSITY)) then
              if (USER_SUPPLIED_SPARSITY) then
                MF = MF - 100*(MF/100)
              end if
            end if
          end if
        end if

        if (present(BANDED_J)) then
          if (BANDED_J) then
            if (present(USER_SUPPLIED_JACOBIAN)) then
              if (USER_SUPPLIED_JACOBIAN) then
                MF = 24
              else
                MF = 25
              end if
            else
              MF = 25
            end if
          end if
        end if

        if (present(DENSE_J)) then
          if (DENSE_J) then
            if (present(USER_SUPPLIED_JACOBIAN)) then
              if (USER_SUPPLIED_JACOBIAN) then
                MF = 21
              else
                MF = 22
              end if
            else
              MF = 22
            end if
          end if
        end if

!       Check for errors in MF.
        MFA = IABS(MF)
        MOSS = MFA/100
        METH = (MFA-100*MOSS)/10
        MITER = MFA - 100*MOSS - 10*METH
        if (METH<1 .or. METH>2) then
          MSG = 'Illegal value of METH in SET_INTERMEDIATE_OPTS.'
          call XERRDV(MSG,330,2,0,0,0,0,ZERO,ZERO)
        end if
        if (MITER<0 .or. MITER>7) then
          MSG = 'Illegal value of MITER in SET_INTERMEDIATE_OPTS.'
          call XERRDV(MSG,340,2,0,0,0,0,ZERO,ZERO)
        end if
        if (MOSS<0 .or. MOSS>2) then
          MSG = 'Illegal value of MOSS in SET_INTERMEDIATE_OPTS.'
          call XERRDV(MSG,350,2,0,0,0,0,ZERO,ZERO)
        end if

!       Reset MF, now that MOSS is known.
        MFSIGN = sign(1,MF)
        MF = MF - 100*MOSS*MFSIGN

        if (MITER==0) then
          DENSE = .false.
          BANDED = .false.
          SPARSE = .false.
        else if (MITER==1 .or. MITER==2) then
          DENSE = .true.
          BANDED = .false.
          SPARSE = .false.
        else if (MITER==3) then
          DENSE = .false.
          BANDED = .false.
          SPARSE = .false.
        else if (MITER==4 .or. MITER==5) then
          DENSE = .false.
          BANDED = .true.
          SPARSE = .false.
        else if (MITER==6 .or. MITER==7) then
          DENSE = .false.
          BANDED = .false.
          SPARSE = .true.
        end if

!       Define the banded Jacobian band widths.
        if (BANDED) then
          if (present(LOWER_BANDWIDTH)) then
            ML = LOWER_BANDWIDTH
            IUSER(1) = ML
          else
            MSG = 'In SET_INTERMEDIATE_OPTS you have indicated a'
            call XERRDV(MSG,360,1,0,0,0,0,ZERO,ZERO)
            MSG = 'banded Jacobian but you have not'
            call XERRDV(MSG,360,1,0,0,0,0,ZERO,ZERO)
            MSG = 'supplied the lower bandwidth.'
            call XERRDV(MSG,360,2,0,0,0,0,ZERO,ZERO)
          end if
          if (present(UPPER_BANDWIDTH)) then
            MU = UPPER_BANDWIDTH
            IUSER(2) = MU
          else
            MSG = 'In SET_INTERMEDIATE_OPTS you have indicated a'
            call XERRDV(MSG,370,1,0,0,0,0,ZERO,ZERO)
            MSG = 'banded Jacobian but you have not'
            call XERRDV(MSG,370,1,0,0,0,0,ZERO,ZERO)
            MSG = 'supplied the upper bandwidth.'
            call XERRDV(MSG,370,2,0,0,0,0,ZERO,ZERO)
          end if
!         Define the nonzero diagonals.
          BNGRP = 0
          SUBS = .false.
          SUPS = .false.
          NSUBS  = 0
          NSUPS = 0
        end if

!       Define the sparse Jacobian options.
        SCALE_MATRIX = .false.
        ELBOW_ROOM = 2
        if (SPARSE) then
!         NZSWAG for the number of nonzeros in the Jacobian.
          if (present(NZSWAG)) then
            NZ_SWAG = max(NZSWAG,0)
          else
            NZ_SWAG = 0
          end if
!         Indicate that SET_IAJA has not yet been called successfully.
          IAJA_CALLED = .false.
!         Check for illegal method flags.
          if (MOSS==2 .and. MITER/=7) then
            MSG = 'In SET_INTERMEDIATE_OPTS MOSS=2 but MITER is not 7.'
            call XERRDV(MSG,380,2,0,0,0,0,ZERO,ZERO)
          end if
          if (MOSS==1 .and. MITER/=6) then
            MSG = 'In SET_INTERMEDIATE_OPTS MOSS=1 but MITER is not 6.'
            call XERRDV(MSG,390,2,0,0,0,0,ZERO,ZERO)
          end if
!         IF (MOSS==0 .AND. MITER/=7) THEN
!           MSG = 'In SET_INTERMEDIATE_OPTS MOSS=0 but MITER is not 7.'
!           CALL XERRDV(MSG,400,2,0,0,0,0,ZERO,ZERO)
!         END IF
!         Allow MC19 scaling for the Jacobian.
          SCALE_MATRIX = .false.
        end if

!       Define the number of event functions.
        if (present(NEVENTS)) then
          if (NEVENTS>0) then
            NG = NEVENTS
          else
            NG = 0
          end if
        else
          NG = 0
        end if

!       Process the constrained solution components.
        if (present(CONSTRAINED)) then
          NDX = size(CONSTRAINED)
          if (NDX<1) then
            MSG = 'In SET_INTERMEDIATE_OPTS the size of CONSTRAINED < 1.'
            call XERRDV(MSG,410,2,0,0,0,0,ZERO,ZERO)
          end if
          if (.not.(present(CLOWER)) .or. .not.(present(CUPPER))) then
            MSG = 'In SET_INTERMEDIATE_OPTS the arrays CLOWER and CUPPER'
            call XERRDV(MSG,420,1,0,0,0,0,ZERO,ZERO)
            MSG = 'are not present.'
            call XERRDV(MSG,420,2,0,0,0,0,ZERO,ZERO)
          end if
          if (size(CLOWER)/=NDX .or. size(CUPPER)/=NDX) then
            MSG = 'In SET_INTERMEDIATE_OPTS the size of the solution bound'
            call XERRDV(MSG,430,1,0,0,0,0,ZERO,ZERO)
            MSG = 'arrays must be the same as the CONSTRAINED array.'
            call XERRDV(MSG,430,2,0,0,0,0,ZERO,ZERO)
          end if
!         Note: The contents of CONSTRAINED will be checked
!         in subroutine DVODE after NEQ is known.
          if (allocated(IDX)) then
            deallocate (IDX,LB,UB,STAT=IER)
            call CHECK_STAT(IER,30)
          end if
          allocate (IDX(NDX),LB(NDX),UB(NDX),STAT=IER)
          call CHECK_STAT(IER,40)
          IDX(1:NDX) = CONSTRAINED(1:NDX)
          LB(1:NDX) = CLOWER(1:NDX)
          UB(1:NDX) = CUPPER(1:NDX)
          BOUNDS = .true.
        else
          if (allocated(IDX)) then
            deallocate (IDX,LB,UB,STAT=IER)
            call CHECK_STAT(IER,50)
          end if
          NDX = 0
          BOUNDS = .false.
        end if

!       Is the Jacobian constant?
        J_IS_CONSTANT = .false.
        J_HAS_BEEN_COMPUTED = .false.

!       Load the user options into the solution structure.
        OPTS%MF = MF
        OPTS%METH = METH
        OPTS%MITER = MITER
        OPTS%MOSS = MOSS
        OPTS%DENSE = DENSE
        OPTS%BANDED = BANDED
        OPTS%SPARSE = SPARSE
        OPTS%NG = NG

!       Process the miscellaneous options.

!       Don't step past TCRIT variable.
        if (present(TCRIT)) then
          RUSER(1) = TCRIT
        else
          RUSER(1) = ZERO
        end if

!       DVODE optional parameters.
        IOPT = 1
        if (present(MAXORD)) then
          IUSER(5) = MAXORD
          IOPT = 1
        end if
        if (present(MXSTEP)) then
          IUSER(6) = MXSTEP
          IOPT = 1
        end if
        if (present(MXHNIL)) then
          IUSER(7) = MXHNIL
          IOPT = 1
        end if
        if (present(H0)) then
          RUSER(5) = H0
          IOPT = 1
        end if
        if (present(HMAX)) then
          RUSER(6) = HMAX
          IOPT = 1
        end if
        if (present(HMIN)) then
          RUSER(7) = HMIN
          IOPT = 1
        end if
        U_PIVOT = ONE
        OPTS%IOPT = IOPT

!       Define the error tolerances.

!       Relative error tolerances.
        NRE = 1
        allocate (OPTS%RTOL(NRE),STAT=IER)
        call CHECK_STAT(IER,60)
        if (present(RELERR)) then
          if (RELERR<ZERO) then
            MSG = 'RELERR must be nonnegative.'
            call XERRDV(MSG,440,2,0,0,0,0,ZERO,ZERO)
          end if
          OPTS%RTOL = RELERR
        else
          if (ALLOW_DEFAULT_TOLS) then
             OPTS%RTOL = 1.0E-4_WP
             MSG = 'By not specifying RELERR, you have elected to use a default'
             call XERRDV(MSG,450,1,0,0,0,0,ZERO,ZERO)
             MSG = 'relative error tolerance equal to 1.0D-4. Please be aware a'
             call XERRDV(MSG,450,1,0,0,0,0,ZERO,ZERO)
             MSG = 'tolerance this large is not appropriate for all problems.'
             call XERRDV(MSG,450,1,0,0,0,0,ZERO,ZERO)
             MSG = 'Execution will continue'
             call XERRDV(MSG,450,1,0,0,0,0,ZERO,ZERO)
          else
             MSG = 'You must specify a nonzero relative error tolerance.'
             call XERRDV(MSG,460,2,0,0,0,0,ZERO,ZERO)
          end if
        end if

!       Absolute error tolerances.
        if (present(ABSERR_VECTOR)) then
          if (minval(ABSERR_VECTOR)<ZERO) then
            MSG = 'All components of ABSERR_VECTOR must'
            call XERRDV(MSG,460,1,0,0,0,0,ZERO,ZERO)
            MSG = 'be nonnegative.'
            call XERRDV(MSG,460,2,0,0,0,0,ZERO,ZERO)
          end if
          NAE = size(ABSERR_VECTOR)
        else
          NAE = 1
        end if
        allocate (OPTS%ATOL(NAE),STAT=IER)
        call CHECK_STAT(IER,70)
        if (present(ABSERR_VECTOR)) then
          OPTS%ATOL = ABSERR_VECTOR
        else if (present(ABSERR)) then
          if (ABSERR<ZERO) then
            MSG = 'ABSERR must be nonnegative.'
            call XERRDV(MSG,470,2,0,0,0,0,ZERO,ZERO)
          end if
          OPTS%ATOL = ABSERR
        else
          if (ALLOW_DEFAULT_TOLS) then
             OPTS%ATOL = 1D-6
             MSG = 'By not specifying ABSERR, you have elected to use a default'
             call XERRDV(MSG,480,1,0,0,0,0,ZERO,ZERO)
             MSG = 'absolute error tolerance equal to 1.0D-6. Please be aware a'
             call XERRDV(MSG,480,1,0,0,0,0,ZERO,ZERO)
             MSG = 'tolerance this large is not appropriate for all problems.'
             call XERRDV(MSG,480,1,0,0,0,0,ZERO,ZERO)
             MSG = 'Execution will continue'
             call XERRDV(MSG,480,1,0,0,0,0,ZERO,ZERO)
          else
             MSG = 'You must specify a vector of absolute error tolerances or'
             call XERRDV(MSG,490,1,0,0,0,0,ZERO,ZERO)
             MSG = 'a scalar error tolerance. It is recommended that you use'
             call XERRDV(MSG,490,1,0,0,0,0,ZERO,ZERO)
             MSG = 'a vector of absolute error tolerances.'
             call XERRDV(MSG,490,2,0,0,0,0,ZERO,ZERO)
          end if
        end if

!       ITOL error tolerance flag.
!          ITOL   RTOL     ATOL            EWT(i)
!            1   scalar   scalar  RTOL*ABS(Y(i)) + ATOL
!            2   scalar   array   RTOL*ABS(Y(i)) + ATOL(i)
!            3   array    scalar  RTOL(i)*ABS(Y(i)) + ATOL
!            4   array    array   RTOL(i)*ABS(Y(i)) + ATOL(i)
        if (present(ABSERR_VECTOR)) then
          OPTS%ITOL = 2
        else
          OPTS%ITOL = 1
        end if
        return

  end function SET_INTERMEDIATE_OPTS
!_______________________________________________________________________

  function SET_OPTS(METHOD_FLAG, DENSE_J, BANDED_J, SPARSE_J,          &
    USER_SUPPLIED_JACOBIAN, SAVE_JACOBIAN, CONSTANT_JACOBIAN,          &
    LOWER_BANDWIDTH, UPPER_BANDWIDTH, SUB_DIAGONALS, SUP_DIAGONALS,    &
    RELERR, RELERR_VECTOR, ABSERR, ABSERR_VECTOR, TCRIT, H0, HMAX,     &
    HMIN, MAXORD, MXSTEP, MXHNIL, YMAGWARN, SETH, UPIVOT, NZSWAG,      &
    USER_SUPPLIED_SPARSITY, NEVENTS, CONSTRAINED, CLOWER, CUPPER,      &
    MA28_ELBOW_ROOM, MC19_SCALING, MA28_MESSAGES, MA28_EPS,            &
    MA28_RPS, CHANGE_ONLY_f77_OPTIONS,JACOBIAN_BY_JACSP)               &
  result(OPTS)

! FUNCTION SET_OPTS:
!    VODE.f77 method flag:
!       METHOD_FLAG
!    Jacobian type:
!       DENSE_J, BANDED_J, SPARSE_J
!    Analytic Jacobian:
!       USER_SUPPLIED_JACOBIAN
!    Jacobian options:
!       SAVE_JACOBIAN, CONSTANT_JACOBIAN. JACOBIAN_BY_JACSP
!    If banded Jacobian:
!       LOWER_BANDWIDTH, UPPER_BANDWIDTH
!    If specify the nonzero diagonals for banded Jacobian:
!       SUB_DIAGONALS, SUP_DIAGONALS
!    Error tolerances:
!       RELERR, RELERR_VECTOR, ABSERR, ABSERR_VECTOR
!    VODE.f77 optional parameters:
!       TCRIT, H0, HMAX, HMIN, MAXORD, MXSTEP, MXHNIL
!    Solution smaller than absolute error tolerance:
!       YMAGWARN
!    Sparse flags:
!       SETH, UPIVOT, NZSWAG, USER_SUPPLIED_SPARSITY
!       MA28_ELBOW_ROOM, MC19_SCALING, MA28_MESSAGES
!       MA28_EPS, MA28_RPS
!    Rootfinding:
!       NEVENTS
!    Impose bounds on solution:
!       CONSTRAINED, CLOWER, CUPPER
!    Change one or more of HMAX, HMIN, TCRIT, MXSTEP, MXHNIL, MAXORD
!       CHANGE_ONLY_f77_OPTIONS
! RESULT(OPTS)

! SET_OPTIONS sets user parameters for DVODE via keywords. Values that
! are defined herein will be used internally by DVODE. All option
! keywords are OPTIONAL and order is not important.

! Note that DVODE_F90 requires that one of SET_NORMAL_OPTS or
! SET_INTERMEDIATE_OPTS or SET_OPTS is called before the first
! time DVODE_F90 is called.

!                     Quick Summary of Options

! METHOD_FLAG            - any legal value of MF as in DVODE
! DENSE_J                - Jacobian is sparse alternative to MF
! BANDED_J               - Jacobian is banded alternative to MF 
! SPARSE_J               - Jacobian is sparse alternative to MF   
! USER_SUPPLIED_JACOBIAN - User will supply Jacobian subroutine
!                          (user supplied subroutine JAC required)
! SAVE_JACOBIAN          - Jacobian will be saved and reused
! CONSTANT_JACOBIAN      - Jacobian is constant
! JACOBIAN_BY_JACSP      - Use Doug Salane's approximate Jacobian
!                          algorithm
! LOWER_BANDWIDTH        - Lower bandwidth ML in DVODE
! UPPER_BANDWIDTH        - Upper bandwidth MU in DVODE
! SUB_DIAGONALS          - Nonzero sub diagonals in Jacobian
! SUP_DIAGONALS          - Nonzero super diagonals in Jacobian
! RELERR                 - Scalar relative error tolerance in DVODE
! RELERR_VECTOR          - Vector relative error tolerance in DVODE
! ABSERR                 - Scalar absolute error tolerance in DVODE
! ABSERR_VECTOR          - Vector absolute error tolerance in DVODE
! TCRIT                  - Critical time TCRIT in DVODE
! H0                     - Starting step size in DVODE
! HMAX                   - Maximum step size in DVODE 
! HMIN                   - Minimum step size in DVODE 
! MAXORD                 - Maximum integration order in DVODE
! MXSTEP                 - Maximum number of integration steps
!                          in DVODE
! MXHNIL                 - Maximum number of T+H=T messages in DVODE
! YMAGWARN               - Warn if magnitude of solution is smaller
!                          than absolute error tolerance
! SETH                   - Sparse Jacobian element threshold value
! UPIVOT                 - MA28 partial pivoting parameter
! NZSWAG                 - guess for the number of nonzeros in sparse
!                          Jacobian
! USER_SUPPLIED_SPARSITY - user will supply sparsity structure
!                          arrays by calling USERSETS_IAJA 
! MA28_ELBOW_ROOM        - Supply an integer greater than 2 if you
!                          wish to increase the elbow room in the
!                          MA28 work arrays (by MA28_ELBOW_ROOM * NZA).
! MC19_SCALING           - .TRUE. if you wish to invoke sparse scaling
!                           of the Jacobian.
! MA28_MESSAGES          - Control printing of MA28 messages
! MA28_RPS               - Redo MA28AD pivot sequence if a singularity
!                          is encountered
! NEVENTS                - number of user defined root finding functions
!                          (user supplied subroutine G required)
! CONSTRAINED,           - array of solution component indices that
! CLOWER,                  are to be constrained by the lower and
! CUPPER                   upper bound arrays CLOWER and CUPPER so that
!                          CLOWER(I) <= Y(CONSTRAINED(I)) <= CUPPER(I)
! CHANGE_ONLY_f77_OPTIONS- flag to indicate whether to only change any
!                          of the parameters MXSTEP, MXHNIL, MAXORD,
!                          HMAX, HMIN, TCRIT
!                     Options Types
! METHOD_FLAG            - integer
! DENSE_J                - logical
! BANDED_J               - logical
! SPARSE_J               - logical
! USER_SUPPLIED_JACOBIAN - logical
! SAVE_JACOBIAN          - logical
! CONSTANT_JACOBIAN      - logical
! JACOBIAN_BY_JACSP      - logical
! LOWER_BANDWIDTH        - integer
! UPPER_BANDWIDTH        - integer
! SUB_DIAGONALS          - integer array
! SUP_DIAGONALS          - integer array
! RELERR                 - real(wp) scalar
! RELERR_VECTOR          - real(wp) vector
! ABSERR                 - real(wp) scalar
! ABSERR_VECTOR          - real(wp) vector
! TCRIT                  - real(wp) scalar
! H0                     - real(wp) scalar
! HMAX                   - real(wp) scalar
! HMIN                   - real(wp) scalar
! MAXORD                 - integer
! MXSTEP                 - integer
! MXHNIL                 - integer
! YMAGWARN               - logical
! SETH                   - real(wp) scalar
! UPIVOT                 - real(wp) scalar
! NZSWAG                 - integer
! USER_SUPPLIED_SPARSITY - logical
! NEVENTS                - integer
! CONSTRAINED            - integer array
! CLOWER                 - real(wp) array
! CUPPER                 - real(wp) array
! CHANGE_ONLY_f77_OPTIONS- logical

! Argument list parameters:
! METHOD_FLAG              = integration method flag
! ABSERR                   = scalar absolute error tolerance
! ABSERR_VECTOR            = vector of absolute error tolerances
! RELERR                   = scalar relative error tolerance
! RELERR_VECTOR            = vector of relative error tolerances
! NEVENTS                  = Number of event functions (requires
!                            user-supplied GFUN)
! DENSE_J                  = use dense linear algebra if .TRUE.
! BANDED_J                 = use banded linear algebra if .TRUE.
!   LOWER_BANDWIDTH        = lower bandwidth of the Jacobian
!                            (required if BANDED_J = .TRUE.)
!   UPPER_BANDWIDTH        = upper bandwidth of the Jacobian
!                            (required if BANDED_J = .TRUE.)
!   SUB_DIAGONALS          = starting row numbers of nonzero sub
!                            diagonals counting up from the lowest
!                            sub diagonal
!   SUP_DIAGONALS          = starting column numbers of nonzero
!                            super diagonals counting up from the
!                            lowest super diagonal
! SPARSE_J                 = use sparse linear algebra if .TRUE.
!   UPIVOT                 = partial pivot control flag (default=ONE)
!                          = 0.0D0 for no partial pivoting
!                          = ONE for (full) partial pivoting
!                            Values between 0.0D0 and ONE yield
!                            MA28 (partial) partial pivoting.
!   NZSWAG                 = If you wish you may supply a guess,
!                            NZSWAG, at the number of nonzeros
!                            in the Jacobian matrix. In some cases
!                            this will speed up the determination
!                            of the necessary storage.
!   USER_SUPPLIED_SPARSITY = .TRUE. if you wish to supply the sparsity
!                            structure arrays directly by calling
!                            USERSETS_IAJA
!   MA28_ELBOW_ROOM        = Supply an integer greater than 2 if you
!                            wish to increase the elbow room in the
!                            MA28 work arrays (by MA28_ELBOW_ROOM * NZA).
!   MC19_SCALING           = .TRUE. if you wish to invoke sparse scaling
!                            of the Jacobian.
!   MA28_MESSAGES          = .TRUE. to print all MA28 messages.
!   MA28_EPS               = MA28 singularity threshold. If MA28BD
!                            determines that the ratio of the current
!                            pivot to the largest element in the row
!                            is smaller than EPS, it will flag the
!                            Jacobian as singular and force MA28AD
!                            to calculate a new pivot sequence.
!   MA28_RPS               = .TRUE. to force MA28AD to calculate a new
!                            pivot sequence for use in MA28BD if a
!                            singular iteration matrix is encountered
!                            is smaller than EPS, it will flag the
!                            Jacobian as singular and force MA28AD
!                            to calculate a new pivot sequence.
! USER_SUPPLIED_JACOBIAN   = exact Jacobian option
!                            (requires user-supplied JAC)
! SAVE_JACOBIAN            = reuse saved Jacobians
! JACOBIAN_BY_JACSP        = use Doug Salane's Jacobian algorithm
! TCRIT                    = critical time
! H0                       = initial step size to try
! HMAX                     = maximum allowable step size
! HMIN                     = minimum allowable step size
! MAXORD                   = maximum allowable integration order
! MXSTEP                   = maximum number of integration steps
!                            between calls to DVODE
! MXHNIL                   = maximum number of printed messages
!                            if the T+H=T condition occurs
! CONSTRAINED              = array containing the indices of
!                            solution components which are to be
!                            be constrained by CLOWER and CUPPER.
!                            The size of CONSTRAINED must be
!                            positive and not exceed NEQ, the
!                            number of ODEs. Each component of
!                            CONSTRAINED must be between 1 and NEQ.
! CLOWER, CUPPER           = lower and upper bound arrays for the.
!                            solution. Each must be the same size
!                            as the CONSTRAINED vector.
! YMAGWARN                 = If .TRUE., a warning will be issued
!                            before any return from DVODE for any
!                            solution component whose magnitude
!                            is smaller than the absolute error
!                            tolerance for that component.
!   User-callable Routines:
!   The following routines may be called by the user:
!   SET_OPTS      : Used to set options for DVODE_F90
!   GET_STATS     : Used to gather summary integration
!                   statistics following a successful
!                   return from DVODE_F90
!   DVINDY        : Used if the user wishes to interpolate
!                   solution or derivative following a
!                   successful return from DVODE_F90
!   USERSETS_IAJA : Used if the user wishes to supply the
!                   sparsity structure directly
!   Detailed Description of SET_OPTS
!   The following are defined in SET_OPTS:
!   OPTS%MF     = Integration method flag (MF)
!   OPTS%METH   = Integration family flag (METH)
!   OPTS%MITER  = Iteration type flag (MITER)
!   OPTS%MOSS   = Sparsity array type flag (MOSS)
!   OPTS%ITOL   = Error tolerance flag (ITOL) (ITOL)
!   OPTS%ATOL   = Absolute error tolerance(s) (ATOL)
!   OPTS%RTOL   = Relative error tolerance(s) (RTOL)
!   OPTS%DENSE  = Use dense linear algebra
!   OPTS%BANDED = Use banded linear algebra 
!   OPTS%SPARSE = Use sparse linear algebra 
!   OPTS%IOPT   = DVODE optional parameter input flag (IOPT)
!   OPTS%NG     = Number of event functions
!   RUSER(1)    = TCRIT (don't step past)
!   IUSER(1)    = Jacobian lower bandwidth (ML)
!   IUSER(2)    = Jacobian upper bandwidth (MU)
!   RUSER(5)    = Initial step size to try (H0)
!   RUSER(6)    = Maximum allowable step size (HMAX)
!   RUSER(7)    = Minimum allowable step size (HMIN)
!   IUSER(5)    = Maximum allowable integration order (MAXORD)
!   IUSER(6)    = Maximum number of integration steps (MXSTEP)
!   IUSER(7)    = Maximum number of T+H=T messages (MXHNIL)
!   NZ_SWAG     = Guess for the number of nonzeros in the
!                 sparse Jacobian
!   NG          = Number of user event functions
!   BOUNDS,     = Nonnegativity information
!   NDX, IDX
!   YMAXWARN    = Warning flag for |y(t)| < abserr
!   MA28_ELBOW_ROOM = Integer multiple by which to increase
!                     the elbow room in the MA28 sparse work
!                     arrays
!   MC19_SCALING    = logical flag to control MC19 sparse
!                     scaling of the Jacobian
!   MA28_MESSAGES   = logical flag to control printing of MA28
!                     messages.
!   MA28_EPS        = real(wp) MA28 singularity threshold
!
!                         All Options
!
! METHOD_FLAG
!
! Default:   Not used
! Change to: No need to specify if use one of next three parameters
!            but can be changed to any value of the MF method flag
!            in dvode.f77
!
! DENSE_J
!
! Default:   .FALSE.
! Change to: .TRUE. for dense linear algebra
!
! BANDED_J
!
! Default:   .FALSE.
! Change to: .TRUE. for banded linear algebra
!
! SPARSE_J
!
! Default:   .FALSE.
! Change to: .TRUE. for sparse linear algebra

! USER_SUPPLIED_JACOBIAN
!
! Default:   .FALSE.
! Change to: .TRUE. if want to supply subroutine JAC to DVODE_F90
!
! SAVE_JACOBIAN
!
! Default:   .TRUE.
! Change to: .FALSE. if do not want DVODE_F90 to reuse saved
!            Jacobians
!
! JACOBIAN_BY_JACSP
!
! Default:   .FALSE.
! Change to: .TRUE. to approximate the Jacobian using Doug Salane's
!            JACSP Jacobian subroutine
!
! CONSTANT_JACOBIAN
!
! Default:   .FALSE.
! Change to: .TRUE. if Jacobian is constant
!
! LOWER_BANDWIDTH
!
! Default:   Not used
! Change to: Lower bandwidth if BANDED_J = .TRUE.
!
! UPPER_BANDWIDTH
!
! Default:   Not used
! Change to: Upper bandwidth if BANDED_J = .TRUE.
!
! SUB_DIAGONALS
!
! Default:   Not used
! Change to: Starting row locations for sub diagonals in
!            banded Jacobian
!
! SUP_DIAGONALS
!
! Default:   Not used
! Change to: Starting columns locations for super diagonals
!            in banded Jacobian
!
! RELERR
!
! Default:   None
! Change to: Specify a nonzero relative error tolerance
!
! RELERR_VECTOR
!
! Default:   Not used
! Change to: Vector of relative error tolerances
!
! ABSERR
!
! Default:   None
! Change to: Specify a scalar absolute error tolerance or
!            a vector of absolute error tolerances
!
! ABSERR_VECTOR
!
! Default:   Not used
! Change to: Vector of absolute error tolerances
!
! TCRIT
!
! Default:   Not used
! Change to: Value of T which DVODE_F90 is not to step past
!
!     H0
!
! Default:   Determined by DVODE_F90
! Change to: Guess for initial step size
!
! HMAX
!
! Default:   Infinity
! Change to: Maximum allowable step size
!
! HMIN
!
! Default:   0.0D0
! Change to: Minimum allowable step size
!
! MAXORD
!
! Default:   Determined by DVODE_F90
! Change to: Maximum integration order
!
! MXSTEP
!
! Default:   5000
!     Change to: Maximum number of steps between output points
!
! MXHNIL
!
! Default:   10
! Change to: Maximum number of times T+H=T message will be printed
!
! YMAGWARN
!
! Default:   .FALSE.
! Change to: .TRUE. if wish to be warned when the magnitude of the
!            solution is smaller than the absolute error tolerance
!
! SETH
!
! Default:   0.0D0
! Change to: Threshold value below which Jacobian elements in the
!            sparse Jacobian are considered to be zero
!
! UPIVOT
!
! Default:   1.0D0
! Change to: partial pivoting factor between 0.0d0 and 1.0d0 to
!            control the degree of partial pivoting used in sparse
!            Gaussian Elimination
!
! NZSWAG
!
! Default:   Determined by DVODE_F90
! Change to: Amount by which to increment sizes of sparse arrays
!
! USER_SUPPLIED_SPARSITY
!
! Default:   .FALSE.
! Change to: .TRUE. if wish to call subroutine USERSETS_IAJA to
!            define the sparse structure array pointers directly
!
! NEVENTS
!
! Default:   0
! Change to: Number of event functions if wish to supply subroutine
!            GFUN to DVODE_F90
!
! CONSTRAINED
!
! Default:   Bounds not imposed on solution
! Change to: Array of indices for components on which to impose
!            solution bounds
!
! CLOWER
!
! Default:   Not used
! Change to: Vector containing lower bounds corresponding to
!            CONSTRAINED
!
! CUPPER
!
! Default:   Not used
! Change to: Vector containing upper bounds corresponding to
!            CONSTRAINED
!
! MA28_ELBOW_ROOM
!
! Default:   2
! Change to: Larger value to allow the sparse arrays more
!            elbow room
!
! MC19_SCALING
!
! Default:   .FALSE.
! Change to: .TRUE. to invoke scaling of the sparse solution
!
! MA28_MESSAGES
!
! Default:   .FALSE.
! Change to: .TRUE. to have MA28 diagnostic messages written
!            to unit 6
!
! MA28_EPS
!
! Default:   1.0D-4
! Change to: Smaller positive value to allow the ratio of the
!            magnitude of pivot elements and remaining row
!            entries to be smaller before the iteration matrix
!            is considered to be numerically singularity
!
! MA28_RPS
!
! Default:   .FALSE.
! Change to: Redo the MA28AD sparse pivoting sequence any time
!            MA28BD considers the iteration matrix to be
!            numerically singularity
!
! Note on Jacobian Storage Formats:
!
! If you supply an analytic Jacobian PD, load the
! Jacobian elements DF(I)/DY(J), the partial
! derivative of F(I) with respect to Y(J), using
! the following formats. Here, Y is the solution,
! F is the derivative, and PD is the Jacobian.
!
! For a full Jacobian, load PD(I,J) with DF(I)/DY(J).
! Your code might look like this:
!    DO J = 1, NEQ
!       DO I = 1, NEQ
!          PD(I,J) = ... DF(I)/DY(J)
!       END DO
!    END DO       
!
! For a banded Jacobian, load PD(I-J+MU+1,J) with
! DF(I)/DY(J) where ML is the lower bandwidth
! and MU is the upper bandwidth of the Jacobian.
! Your code might look like this:
!    DO J = 1, NEQ
!       I1 = MAX(1,J-ML)
!       I2 = MIN(N,J+MU)
!       DO I = I1, I2
!          K = I-J+MU+1
!          PD(K,J) = ... DF(I)/DY(J)
!       END DO
!    END DO
!
! For a sparse Jacobian, IA(J+1)-IA(J) is the number
! of nonzeros in column change. JA(I) indicates the
! rows in which the nonzeros occur. For column J,
! the nonzeros occurs in rows I=JA(K) for K=I1,...,I2
! where I1=IA(J) and I2= IA(J+1)-1. Load DF(I)/DY(J)
! in PD(I). Your code might look like this:
!    DO J = 1, NEQ
!       I1 = IA(J)
!       I2 = IA(J+1) - 1
!       DO K = I1, I2
!          I = JA(K)
!          PD(I) = ... DF(I)/DY(J)
!       END DO
!    END DO
!
!                    More on Sparsity Options
!
! Two facts of life regarding the use of direct sparse solvers are
! (1) significant improvements are possible, and (2) the use of
! direct sparse solvers often is more demanding of the user.
! Although SET_NORMAL_OPTS provides modest provisions for solving
! problems with sparse Jacobians, using SET_OPTS rather than
! SET_NORMAL_OPTS provides several advanced options. These options
! are described below. The recommended manner in which to use
! these options is also provided. Note that each of these optional
! parameters have default values and need not be specified in your
! call to SET_OPTS if you wish to use the default values. Note
! also that the order in which the options are specified in a call
! to SET_OPTS is not important.
!
!  (1) First determine if it is feasible to use either the BANDED
!      or DENSE Jacobian option.
!         Recommendation:
!         Use the DENSE or BANDED option if possible. They do
!         not require use of most of the options described here.
!    
!  (2) The Jacobian is approximated internally using differences
!      if you do not provide an analytic Jacobian. The option
!      USER_SUPPLIED_JACOBIAN=.TRUE. may be used if you wish to
!      provide an analytic Jacobian.
!         Recommendation:
!         Use an internally generated Jacobian for most problems
!         but consider providing an analytic Jacobian if it not
!         too much trouble.
!
!  (3) If you do not provide the sparse structure arrays, they
!      are approximated internally by making NEQ calls to your
!      derivative subroutine and using differences to approximate
!      the structure of the Jacobian. The option
!      USER_SUPPLIED_SPARSITY and a call to USERSETS_IAJA can
!      be used to supply the arrays directly. You can also use
!      subroutine SET_IAJA to approximate the structure using
!      different perturbation factors than those used in DVODE_F90.
!         Recommendations:
!         Although allowing DVODE_F90 to approximate the Jacobian
!         structure suffices for most problems, if you know the
!         sparsity pattern provide it directly. This eliminates
!         the possibility that an important element that happens
!         to be 0 when the sparsity pattern is approximated is
!         later nonzero; and it avoids the NEQ extra calls to
!         your derivative subroutine to approximate the structure.
!         Note that a nemesis for any sparse ode solver is a
!         problem in which the sparsity pattern actually changes
!         during the integration. An example of such a problem
!         is provided in the demohum21.f90 demonstration program.
!         If you know the sparsity pattern changes or if you
!         suspect it does because DVODE_F90 is generating
!         nonconvergence error messages, consider having DVODE_F90
!         re-approximate the structure by calling SET_OPTS and
!         forcing an integration restart when control is returned
!         to your calling program. Please do not do this at every
!         output point since it will be extremely time consuming
!         and inefficient if it is not needed.
!
!  (4) The optional parameter UPIVOT is used to control the type
!      of partial pivoting pivoting in MA28AD. UPIVOT=0.0D0
!      corresponds to no pivoting; and UPIVOT=1.0D0 corresponds
!      to partial pivoting. UPIVOT may be assigned any value
!      between 0 and 1. The pivoting strategy used in MA28AD is
!      to accept a pivot for which the magnitude of the pivot
!      element is greater than or equal to UPIVOT times the
!      magnitude of the largest remaining element in the pivot
!      row. The default value is UPIVOT=1.DO
!         Recommendation:
!         Use the default value UPIVOT=1.0D0.
!
!  (5) The optional parameter NZSWAG may/should be used to speed
!      up the determination of acceptable array sizes for the
!      internal sparse working arrays and the sparse Jacobian.
!      Initially, DVODE_F90 allocates arrays of length
!      max(10*NEQ,NSZSWAG) and increases this amount by
!      max(10*NEQ,ELBOW_ROOM*NSZSWAG) as necessary. NZSWAG
!      has a default value of 1000 and ELBOW_ROOM has a default
!      value of 2.
!         Recommendation:
!         Provide a larger value for NZSWAG particularly if NEQ
!         is large and you suspect considerable fill-in due to
!         partial pivoting.
!
!  (6) The optional parameter ELBOW_ROOM may be used to control
!      the amount of "elbow room" needed in the MA28 sparse arrays.
!      The default value is 2 but a larger value can sometimes
!      speed up the determination of sparse array sizes.
!         Recommendation:
!         MA28 error and diagnostic messages are turned off by
!         default. The optional parameter MA28_MESSAGES may be
!         used to turn them on if it is assigned a nonzero value.
!         If you see several messages that state that LIRN_ALL
!         or LICN_ALL is too small and that additional storage
!         is being allocated for another try, and you have
!         provided a large value for NZSWAG, increase the size
!         of ELBOW_ROOM.
!
!  (7) The optional parameter SETH may be assigned a nonzero value
!      that represents a threshhold value for the magitude of
!      elements in the Jacobian below which the element will be
!      treated as being zero. The default value is SETH=0.0D0.
!         Recommendation:
!         Use the default value SETH=0.0D0.
!
!  (8) MA28AD determines the pivot sequence to be used in subsequent
!      calls to MA28BD. When MA28BD is called it considers the
!      iteration matrix to be numerically singular if the magnitude
!      of the ratio of the largest remaining element in the pivot
!      row to the pivot element is less than EPS. This condition
!      can arise since MA28BD is using an out-of-date pivot sequence.
!      The default value used is EPS=1.0D-4. A smaller value
!      of EPS may be appropriate for some problems. The optional
!      parameter MA28_EPS may be used to change the default value.
!         Recommendation:
!         Use the default value EPS=1.0D-4.
!
!  (9) Badly scaled Jacobians can cause problems for sparse solvers.
!      The parameter MA28_SCALING may be set .TRUE. to instruct
!      DVODE_F90 to scale the iteration matrix using MC19AD. This
!      can positively impact the performance of MA28AD. The default
!      setting is MA28_SCALING=.FALSE.
!         Recommendation:
!         Unless the solution is not successful, use the default
!         setting for MA29_SCALING; but consider using scaling
!         if you know the problem is badly scaled (e.g., if the
!         magnitudes of the components differ greatly).
!
! (10) As mentioned above, MA28BD uses the pivot sequence initially
!      determined by MA28AD. When a singularity is diagnosed, the
!      internal procedure used by the DENSE and BANDED options is
!      used to reduce the step size and re-calulate the iteration
!      matrix. This works fine for most problems; but particularly
!      for badly scaled problems, the solution may be unsuccessful
!      or it may drag along using small step sizes. This is due
!      to the fact MA28BD will continue to use the out-of-date pivot
!      sequence until singularity is again diagnosed. An option not
!      available in previous sparse ode solvers may be used to
!      instruct DVODE_F90 to force MA28AD to calculate a new pivot
!      sequence when MA28BD encounters a singularity. Although
!      MA28AD is considerably slower than MA28BD, the additional
!      calls to MA28AD ensure that the pivot sequence is more
!      up to date. This can increase both the accuracy and the
!      efficiency very dramatically. The demodirn.f90 demonstration
!      program provides an illustration of the dramatic improvement
!      that is possible. The optional parameter MA28_RPS may be set
!      .TRUE. to force these pivot sequence updates.
!         Recommendation:
!         Use the default value MA28_RPS=.FALSE; but if DVODE_F90
!         encounters nonconvergence, use MA28_RPS=.TRUE. to force
!         pivot sequence updates. Note that use of this option
!         will usually obviate the necessity to use the other
!         options described above. Note that we decided not to
!         use MA28_RPS=.TRUE. as default simply because other
!         sparse ode solvers are averse to calling MA28AD more
!         than is absolutely necessary for a given problem.
! ..
     implicit none
! ..
! .. Function Return Value ..
        type (VODE_OPTS) :: OPTS
! ..
! .. Scalar Arguments ..
        real (WP), optional, intent (IN) :: ABSERR, H0, HMAX, HMIN,         &
        MA28_EPS, RELERR, SETH, TCRIT, UPIVOT
        integer, optional, intent (IN) :: LOWER_BANDWIDTH, MAXORD,          &
          MA28_ELBOW_ROOM, METHOD_FLAG, MXHNIL, MXSTEP, NEVENTS, NZSWAG,    &
          UPPER_BANDWIDTH
        logical, optional, intent (IN) :: BANDED_J, CHANGE_ONLY_f77_OPTIONS,&
          CONSTANT_JACOBIAN, DENSE_J, JACOBIAN_BY_JACSP, MA28_MESSAGES,     &
          MA28_RPS, MC19_SCALING, SAVE_JACOBIAN, SPARSE_J,                  &
          USER_SUPPLIED_JACOBIAN, USER_SUPPLIED_SPARSITY, YMAGWARN
! ..
! .. Array Arguments ..
        real (WP), optional, intent (IN) :: ABSERR_VECTOR(:), RELERR_VECTOR(:)
        real (WP), optional :: CLOWER(:), CUPPER(:)
        integer, optional, intent (IN) :: CONSTRAINED(:), SUB_DIAGONALS(:), &
        SUP_DIAGONALS(:)
! ..
! .. Local Scalars ..
        integer ::  IER, IOPT, METH, MF, MFA, MFSIGN, MITER, ML, MOSS, MU,  &
          NAE, NG, NRE
        logical :: BANDED, DENSE, SPARSE
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
        intrinsic ALLOCATED, IABS, MAX, MINVAL, PRESENT, SIGN, SIZE
! ..
! .. FIRST EXECUTABLE STATEMENT SET_OPTS
! ..
    RUSER(1:LRWUSER) = ZERO
    IUSER(1:LIWUSER) = 0
    
!   Allow default error tolerances?
    ALLOW_DEFAULT_TOLS = .false.

!   Maximum number of consecutive error test failures?
    CONSECUTIVE_EFAILS = KFH

!   Maximum number of consecutive corrector iteration failures?
    CONSECUTIVE_CFAILS = MXNCF

!   Use JACSP to approximate Jacobian?
    USE_JACSP = .false.
    if (present(JACOBIAN_BY_JACSP)) then
       if (JACOBIAN_BY_JACSP) USE_JACSP = .true.
    end if

!   If only f77 options are to be changed, do it and return.
    if (present(CHANGE_ONLY_f77_OPTIONS)) then
       if (CHANGE_ONLY_f77_OPTIONS) then
          if (.not.OPTS_CALLED) then
             MSG = 'You have not previously called SET_OPTS before attempting'
             call XERRDV(MSG,500,1,0,0,0,0,ZERO,ZERO)
             MSG = 'to change one or more of the vode.f77 optional parameters.'
             call XERRDV(MSG,500,2,0,0,0,0,ZERO,ZERO)
          end if
          if (present(HMAX)) then
             IOPT = 1
             RUSER(6) = HMAX
             MSG = 'HMAX changed in SET_OPTS.'
             call XERRDV(MSG,510,1,0,0,0,1,HMAX,ZERO)
          end if
          if (present(HMIN)) then
             IOPT = 1
             RUSER(7) = HMIN
             MSG = 'HMIN changed in SET_OPTS.'
             call XERRDV(MSG,520,1,0,0,0,1,HMIN,ZERO)
          end if
          if (present(TCRIT)) then
             IOPT = 1
             RUSER(1) = TCRIT
             MSG = 'TCRIT changed in SET_OPTS.'
             call XERRDV(MSG,530,1,0,0,0,1,TCRIT,ZERO)
          end if
          if (present(MXSTEP)) then
             IOPT = 1
             IUSER(6) = MXSTEP
             MSG = 'MXSTEP changed in SET_OPTS.'
             call XERRDV(MSG,530,1,1,MXSTEP,0,0,ZERO,ZERO)
          end if
          if (present(MAXORD)) then
             IOPT = 1
             IUSER(5) = MAXORD
             MSG = 'MAXORD changed in SET_OPTS.'
             call XERRDV(MSG,540,1,1,MAXORD,0,0,ZERO,ZERO)
          end if
          if (present(MXHNIL)) then
             IOPT = 1
             IUSER(7) = MXHNIL
             MSG = 'MXHNIL changed in SET_OPTS.'
             call XERRDV(MSG,550,1,1,MXHNIL,0,0,ZERO,ZERO)
          end if
       end if
       return
    end if

!   Set the flag to indicate that SET_OPTS has been called.
    OPTS_CALLED = .true.

!   Set the MA48 storage cleanup flag.
    MA48_WAS_USED = .false.

!   Set the fast factor option for MA48,
    USE_FAST_FACTOR = .true.

!   Determine the working precision and define the value for UMAX
!   expected by MA28. Note that it is different for single and
!   double precision.
    WPD = kind(1.0D0)
    WPS = kind(1.0E0)
    if (WPD/=WP .and. WPS/=WP) then
      MSG = 'Illegal working precision in SET_OPTS.'
      call XERRDV(MSG,560,2,0,0,0,0,ZERO,ZERO)
    end if
    if (WPD==WP) then
!   Working precision is double.
      UMAX = 0.999999999_WP
    else
!     Working precision is single.
      UMAX = 0.9999_WP
    end if

!   Set the MA28 message flag.
    if (present(MA28_MESSAGES)) then
       LP = 0
       if (MA28_MESSAGES) LP = 6
    else
       LP = 0
    end if

!   Set the MA28 singularity threshold.
    if (present(MA28_EPS)) then
       if (MA28_EPS > ZERO) then
          EPS = MA28_EPS
       else
          EPS = 1.0E-4_WP
       end if
     else
       EPS = 1.0E-4_WP
     end if

!   Set the MA28 pivot sequence frequency flag.
    if (present(MA28_RPS)) then
       if (MA28_RPS) then
          REDO_PIVOT_SEQUENCE = MA28_RPS
       else
          REDO_PIVOT_SEQUENCE = .false.
       end if
     else
       REDO_PIVOT_SEQUENCE = .false.
     end if

     MA28AD_CALLS = 0
     MA28BD_CALLS = 0
     MA28CD_CALLS = 0
     MC19AD_CALLS = 0
!_______________________________________________________________________
! *****MA48 build change point. Insert these statements.
!    MA48AD_CALLS = 0
!    MA48BD_CALLS = 0
!    MA48CD_CALLS = 0
!_______________________________________________________________________
     IRNCP = 0
     ICNCP = 0
     MINIRN = 0
     MINICN = 0
     MAX_MINIRN = 0
     MAX_MINICN = 0
     MAX_NNZ = 0

!   Set the flag to warn the user if |(y(t)| < abserr.
     if (present(YMAGWARN)) then
       if (YMAGWARN) then
         YMAXWARN = .true.
       else
         YMAXWARN = .false.
       end if
     else
       YMAXWARN = .false.
     end if

!    Load defaults for the optional input arrays for DVODE.
     IUSER(1:8) = 0
     RUSER(1:8) = ZERO

!    Set the method flag.
     if (.not.(present(METHOD_FLAG))) then
       MF = 10
       if (present(SPARSE_J)) then
         if (SPARSE_J) then
           if (present(USER_SUPPLIED_JACOBIAN)) then
             if (USER_SUPPLIED_JACOBIAN) then
               MF = 126
             else
               MF = 227
             end if
           else
             MF = 227
           end if
           if (present(USER_SUPPLIED_SPARSITY)) then
             if (USER_SUPPLIED_SPARSITY) then
               MF = MF - 100*(MF/100)
             end if
           end if
           if (present(SAVE_JACOBIAN)) then
             if (.not.SAVE_JACOBIAN) then
               MF = -MF
             end if
           end if
          end if
          end if

          if (present(BANDED_J)) then
            if (BANDED_J) then
              if (present(USER_SUPPLIED_JACOBIAN)) then
                if (USER_SUPPLIED_JACOBIAN) then
                  MF = 24
                else
                  MF = 25
                end if
              else
                MF = 25
              end if
              if (present(SAVE_JACOBIAN)) then
                if (.not.SAVE_JACOBIAN) then
                  MF = -MF
                end if
              end if
            end if
          end if

          if (present(DENSE_J)) then
            if (DENSE_J) then
              if (present(USER_SUPPLIED_JACOBIAN)) then
                if (USER_SUPPLIED_JACOBIAN) then
                  MF = 21
                else
                  MF = 22
                end if
              else
                MF = 22
              end if
              if (present(SAVE_JACOBIAN)) then
                if (.not.SAVE_JACOBIAN) then
                  MF = -MF
                end if
              end if
            end if
          end if
     else
       MF = METHOD_FLAG
     end if

!    Check for errors in MF.
     MFA = IABS(MF)
     MOSS = MFA/100
     METH = (MFA-100*MOSS)/10
     MITER = MFA - 100*MOSS - 10*METH
     if (METH<1 .or. METH>2) then
       MSG = 'Illegal value of METH in SET_OPTS.'
       call XERRDV(MSG,570,2,0,0,0,0,ZERO,ZERO)
     end if
     if (MITER<0 .or. MITER>7) then
       MSG = 'Illegal value of MITER in SET_OPTS.'
       call XERRDV(MSG,580,2,0,0,0,0,ZERO,ZERO)
     end if
     if (MOSS<0 .or. MOSS>2) then
       MSG = 'Illegal value of MOSS in SET_OPTS.'
       call XERRDV(MSG,580,2,0,0,0,0,ZERO,ZERO)
     end if

!    Reset MF, now that MOSS is known.
     MFSIGN = sign(1,MF)
     MF = MF - 100*MOSS*MFSIGN

     if (MITER==0) then
       DENSE = .false.
       BANDED = .false.
       SPARSE = .false.
     else if (MITER==1 .or. MITER==2) then
       DENSE = .true.
       BANDED = .false.
       SPARSE = .false.
     else if (MITER==3) then
       DENSE = .false.
       BANDED = .false.
       SPARSE = .false.
     else if (MITER==4 .or. MITER==5) then
       DENSE = .false.
       BANDED = .true.
       SPARSE = .false.
     else if (MITER==6 .or. MITER==7) then
       DENSE = .false.
       BANDED = .false.
       SPARSE = .true.
     end if

!    Define the banded Jacobian band widths.
     if (BANDED) then
       if (present(LOWER_BANDWIDTH)) then
         ML = LOWER_BANDWIDTH
         IUSER(1) = ML
       else
         MSG = 'In SET_OPTS you have indicated a'
         call XERRDV(MSG,590,1,0,0,0,0,ZERO,ZERO)
         MSG = 'banded Jacobian but you have not'
         call XERRDV(MSG,590,1,0,0,0,0,ZERO,ZERO)
         MSG = 'supplied the lower bandwidth.'
         call XERRDV(MSG,590,2,0,0,0,0,ZERO,ZERO)
       end if
       if (present(UPPER_BANDWIDTH)) then
         MU = UPPER_BANDWIDTH
         IUSER(2) = MU
       else
         MSG = 'In SET_OPTS you have indicated a'
         call XERRDV(MSG,600,1,0,0,0,0,ZERO,ZERO)
         MSG = 'banded Jacobian but you have not'
         call XERRDV(MSG,600,1,0,0,0,0,ZERO,ZERO)
         MSG = 'supplied the upper bandwidth.'
         call XERRDV(MSG,600,2,0,0,0,0,ZERO,ZERO)
       end if
!      Define the nonzero diagonals.
       BNGRP = 0
       SUBS = .false.
       SUPS = .false.
       NSUBS  = 0
       NSUPS = 0
       if (present(SUB_DIAGONALS)) then
         SUBS = .true.
         NSUBS = size(SUB_DIAGONALS)
         if (NSUBS > 0) then
            if (allocated(SUBDS)) then
               deallocate(SUBDS,STAT=IER)
               call CHECK_STAT(IER,80)
            end if
            allocate(SUBDS(NSUBS),STAT=IER)
            call CHECK_STAT(IER,90)
            SUBDS(1:NSUBS) = SUB_DIAGONALS(1:NSUBS)
         else
            if (ML > 0) then
               MSG = 'You must indicated that the lower bandwidth'
               call XERRDV(MSG,610,1,0,0,0,0,ZERO,ZERO)
               MSG = 'is positive but you have not specified the'
               call XERRDV(MSG,610,1,0,0,0,0,ZERO,ZERO)
               MSG = 'indices for the lower sub diagonals.'
               call XERRDV(MSG,610,2,0,0,0,0,ZERO,ZERO)
            end if
         end if
       end if
       if (present(SUP_DIAGONALS)) then
         SUPS = .true.
         NSUPS = size(SUP_DIAGONALS)
         if (NSUPS > 0) then
            if (allocated(SUPDS)) then
               deallocate(SUPDS,STAT=IER)
               call CHECK_STAT(IER,100)
            end if
            allocate(SUPDS(NSUPS),STAT=IER)
            call CHECK_STAT(IER,110)
            SUPDS(1:NSUPS) = SUP_DIAGONALS(1:NSUPS)
         else
            if (ML > 0) then
               MSG = 'You must indicated that the upper bandwidth'
               call XERRDV(MSG,620,1,0,0,0,0,ZERO,ZERO)
               MSG = 'is positive but you have not specified the'
               call XERRDV(MSG,620,1,0,0,0,0,ZERO,ZERO)
               MSG = 'indices for the upper sub diagonals.'
               call XERRDV(MSG,620,2,0,0,0,0,ZERO,ZERO)
            end if
         end if
       end if
     end if

!    Define the sparse Jacobian options.
     SCALE_MATRIX = .false.
     ELBOW_ROOM = 2
     if (SPARSE) then
!      NZSWAG for the number of nonzeros in the Jacobian.
       if (present(NZSWAG)) then
         NZ_SWAG = max(NZSWAG,0)
       else
         NZ_SWAG = 0
       end if
!      Indicate that SET_IAJA has not yet been called successfully.
       IAJA_CALLED = .false.
!      Check for illegal method flags.
       if (MOSS==2 .and. MITER/=7) then
         MSG = 'In SET_OPTS MOSS=2 but MITER is not 7.'
            call XERRDV(MSG,630,2,0,0,0,0,ZERO,ZERO)
         end if
         if (MOSS==1 .and. MITER/=6) then
           MSG = 'In SET_OPTS MOSS=1 but MITER is not 6.'
           call XERRDV(MSG,640,2,0,0,0,0,ZERO,ZERO)
         end if
!        IF (MOSS==0 .AND. MITER/=7) THEN
!          MSG = 'In SET_OPTS MOSS=0 but MITER is not 7.'
!          CALL XERRDV(MSG,650,2,0,0,0,0,ZERO,ZERO)
!        END IF
!        Allow the work array elbow room to be increased.
         if (present(MA28_ELBOW_ROOM)) then
            ELBOW_ROOM = max(MA28_ELBOW_ROOM, ELBOW_ROOM)
         end if
!      Allow MC19 scaling for the Jacobian.
       SCALE_MATRIX = .false.
       if (present(MC19_SCALING)) then
          if (MC19_SCALING) then
!_______________________________________________________________________
! *****MA48 build change point. Insert these statements.
!         IF (USE_MA48_FOR_SPARSE) THEN
!            MSG = 'MC29AD scaling is not available at this time.'
!            CALL XERRDV(MSG,660,1,0,0,0,0,ZERO,ZERO)
!            MSG = 'Execution will continue.'
!            CALL XERRDV(MSG,660,1,0,0,0,0,ZERO,ZERO)
!         END IF
!_______________________________________________________________________
             SCALE_MATRIX = MC19_SCALING
!            SCALE_MATRIX = .FALSE.
          end if
       end if

!_______________________________________________________________________
! *****MA48 build change point. Replace above with these statements.
!      IF (PRESENT(MC19_SCALING)) THEN
!         IF (MC19_SCALING) THEN
!            IF (USE_MA48_FOR_SPARSE) THEN
!               MSG = 'Please note that this version uses MC19AD rather'
!               CALL XERRDV(MSG,670,1,0,0,0,0,ZERO,ZERO)
!               MSG = 'than MC29AD to sacle the iteration matrix.'
!               CALL XERRDV(MSG,670,1,0,0,0,0,ZERO,ZERO)
!               MSG = 'Execution will continue.'
!               CALL XERRDV(MSG,670,1,0,0,0,0,ZERO,ZERO)
!            END IF
!            SCALE_MATRIX = MC19_SCALING
!         END IF
!       END IF
!_______________________________________________________________________

       end if

!       Define the number of event functions.
        if (present(NEVENTS)) then
          if (NEVENTS>0) then
            NG = NEVENTS
          else
            NG = 0
          end if
        else
          NG = 0
        end if

!       Process the constrained solution components.
        if (present(CONSTRAINED)) then
          NDX = size(CONSTRAINED)
          if (NDX<1) then
            MSG = 'In SET_OPTS the size of CONSTRAINED < 1.'
            call XERRDV(MSG,680,2,0,0,0,0,ZERO,ZERO)
          end if
          if (.not.(present(CLOWER)) .or. .not.(present(CUPPER))) then
            MSG = 'In SET_OPTS the arrays CLOWER and CUPPER are'
            call XERRDV(MSG,690,1,0,0,0,0,ZERO,ZERO)
            MSG = 'not present.'
            call XERRDV(MSG,690,2,0,0,0,0,ZERO,ZERO)
          end if
          if (size(CLOWER)/=NDX .or. size(CUPPER)/=NDX) then
            MSG = 'In SET_OPTS the size of the solution bound arrays'
            call XERRDV(MSG,700,1,0,0,0,0,ZERO,ZERO)
            MSG = 'must be the same as the CONSTRAINED array.'
            call XERRDV(MSG,700,2,0,0,0,0,ZERO,ZERO)
          end if
!         Note: The contents of CONSTRAINED will be checked
!         in subroutine DVODE after NEQ is known.
          if (allocated(IDX)) then
            deallocate (IDX,LB,UB,STAT=IER)
            call CHECK_STAT(IER,120)
          end if
          allocate (IDX(NDX),LB(NDX),UB(NDX),STAT=IER)
          call CHECK_STAT(IER,130)
          IDX(1:NDX) = CONSTRAINED(1:NDX)
          LB(1:NDX) = CLOWER(1:NDX)
          UB(1:NDX) = CUPPER(1:NDX)
          BOUNDS = .true.
        else
          if (allocated(IDX)) then
            deallocate (IDX,LB,UB,STAT=IER)
            call CHECK_STAT(IER,140)
          end if
          NDX = 0
          BOUNDS = .false.
        end if

!       Is the Jacobian constant?
        J_IS_CONSTANT = .false.
        J_HAS_BEEN_COMPUTED = .false.

!_______________________________________________________________________
        if (present(CONSTANT_JACOBIAN)) then
           if (CONSTANT_JACOBIAN) then
              J_IS_CONSTANT = .true.
!_______________________________________________________________________
! *****MA48 build change point. Insert these statements.
              if (USE_MA48_FOR_SPARSE .and. SPARSE) then
                 MSG = 'The constant Jacobian option is not yet available'
                 call XERRDV(MSG,710,1,0,0,0,0,ZERO,ZERO)
                 MSG = 'with the sparse MA48 solution option. Execution'
                 call XERRDV(MSG,710,1,0,0,0,0,ZERO,ZERO)
                 MSG = 'will continue.'
                 call XERRDV(MSG,710,1,0,0,0,0,ZERO,ZERO)
                 J_IS_CONSTANT = .false.
              end if
!_______________________________________________________________________
           end if
        end if
! *****MA48 build change point. Replace above with these statements.
!       IF (PRESENT(CONSTANT_JACOBIAN)) THEN
!          IF (CONSTANT_JACOBIAN) THEN
!             IF (USE_MA48_FOR_SPARSE) THEN
!                MSG = 'The constant Jacobian option is not yet available'
!                CALL XERRDV(MSG,720,1,0,0,0,0,ZERO,ZERO)
!                MSG = 'with the sparse MA48 solution option. Execution'
!                CALL XERRDV(MSG,720,1,0,0,0,0,ZERO,ZERO)
!                MSG = 'will continue.'
!                CALL XERRDV(MSG,720,1,0,0,0,0,ZERO,ZERO)
!             END IF
!          ELSE
!             J_IS_CONSTANT = .TRUE.
!          END IF
!       END IF
!_______________________________________________________________________
        if (J_IS_CONSTANT) then
           if (present(SAVE_JACOBIAN)) then
              if (.not.SAVE_JACOBIAN) then
                 MSG = 'You have specified that the Jacobian is constant.'
                 call XERRDV(MSG,730,1,0,0,0,0,ZERO,ZERO)
                 MSG = 'In this case you cannot also specify that'
                 call XERRDV(MSG,730,1,0,0,0,0,ZERO,ZERO)
                 MSG = 'SAVE_JACOBIAN=.FALSE.'
                 call XERRDV(MSG,730,2,0,0,0,0,ZERO,ZERO)
              end if
           end if
           if (present(USER_SUPPLIED_JACOBIAN)) then
              if (USER_SUPPLIED_JACOBIAN .and. SPARSE) then
                 MSG = 'You have specified that the Jacobian is constant'
                 call XERRDV(MSG,740,1,0,0,0,0,ZERO,ZERO)
                 MSG = 'and that you wish to supply an analytic Jacobian'
                 call XERRDV(MSG,740,1,0,0,0,0,ZERO,ZERO)
                 MSG = 'for a sparse problem. In this case your request'
                 call XERRDV(MSG,740,1,0,0,0,0,ZERO,ZERO)
                 MSG = 'to use a constant Jacobian will be ignored.'
                 call XERRDV(MSG,740,1,0,0,0,0,ZERO,ZERO)
                 MSG = 'Execution will continue.'
                 call XERRDV(MSG,740,1,0,0,0,0,ZERO,ZERO)
              end if
           end if
        end if

!       Load the user options into the solution structure.
        OPTS%MF = MF
        OPTS%METH = METH
        OPTS%MITER = MITER
        OPTS%MOSS = MOSS
        OPTS%DENSE = DENSE
        OPTS%BANDED = BANDED
        OPTS%SPARSE = SPARSE
        OPTS%NG = NG

!       Process the miscellaneous options.

!       Don't step past TCRIT variable.
        if (present(TCRIT)) then
          RUSER(1) = TCRIT
        else
          RUSER(1) = ZERO
        end if

!       DVODE optional parameters.
        IOPT = 1
        if (present(MAXORD)) then
          IUSER(5) = MAXORD
          IOPT = 1
        end if
        if (present(MXSTEP)) then
          IUSER(6) = MXSTEP
          IOPT = 1
        end if
        if (present(MXHNIL)) then
          IUSER(7) = MXHNIL
          IOPT = 1
        end if
        if (present(H0)) then
          RUSER(5) = H0
          IOPT = 1
        end if
        if (present(HMAX)) then
          RUSER(6) = HMAX
          IOPT = 1
        end if
        if (present(HMIN)) then
          RUSER(7) = HMIN
          IOPT = 1
        end if
        if (present(SETH)) then
          RUSER(8) = SETH
          IOPT = 1
        end if
        if (present(UPIVOT)) then
          U_PIVOT = UPIVOT
          if (U_PIVOT<ZERO) U_PIVOT = ZERO
          if (U_PIVOT>ONE) U_PIVOT = ONE
        else
          U_PIVOT = ONE
        end if
!_______________________________________________________________________
! *****MA48 build change point. Insert this statement.
!       COPY_OF_U_PIVOT = U_PIVOT
!_______________________________________________________________________
        OPTS%IOPT = IOPT

!       Define the error tolerances.

!       Relative error tolerances.
        if (present(RELERR_VECTOR)) then
          if (minval(RELERR_VECTOR)<ZERO) then
            MSG = 'All components of RELERR_VECTOR must'
            call XERRDV(MSG,750,1,0,0,0,0,ZERO,ZERO)
            MSG = 'be nonnegative.'
            call XERRDV(MSG,750,2,0,0,0,0,ZERO,ZERO)
          end if
          NRE = size(RELERR_VECTOR)
        else
          NRE = 1
        end if
        allocate (OPTS%RTOL(NRE),STAT=IER)
        call CHECK_STAT(IER,150)
        if (present(RELERR_VECTOR)) then
          OPTS%RTOL = RELERR_VECTOR
        else if (present(RELERR)) then
          if (RELERR<ZERO) then
            MSG = 'RELERR must be nonnegative.'
            call XERRDV(MSG,760,2,0,0,0,0,ZERO,ZERO)
          end if
          OPTS%RTOL = RELERR
        else
          if (ALLOW_DEFAULT_TOLS) then
             OPTS%RTOL = 1.0E-4_WP
             MSG = 'By not specifying RELERR, you have elected to use a default'
             call XERRDV(MSG,770,1,0,0,0,0,ZERO,ZERO)
             MSG = 'relative error tolerance equal to 1.0D-4. Please be aware a'
             call XERRDV(MSG,770,1,0,0,0,0,ZERO,ZERO)
             MSG = 'tolerance this large is not appropriate for all problems.'
             call XERRDV(MSG,770,1,0,0,0,0,ZERO,ZERO)
             MSG = 'Execution will continue'
             call XERRDV(MSG,770,1,0,0,0,0,ZERO,ZERO)
          else
             MSG = 'You must specify a nonzero relative error tolerance.'
             call XERRDV(MSG,780,2,0,0,0,0,ZERO,ZERO)
          end if
        end if

!       Absolute error tolerances.
        if (present(ABSERR_VECTOR)) then
          if (minval(ABSERR_VECTOR)<ZERO) then
            MSG = 'All components of ABSERR_VECTOR must'
            call XERRDV(MSG,790,1,0,0,0,0,ZERO,ZERO)
            MSG = 'be nonnegative.'
            call XERRDV(MSG,790,2,0,0,0,0,ZERO,ZERO)
          end if
          NAE = size(ABSERR_VECTOR)
        else
          NAE = 1
        end if
        allocate (OPTS%ATOL(NAE),STAT=IER)
        call CHECK_STAT(IER,160)
        if (present(ABSERR_VECTOR)) then
          OPTS%ATOL = ABSERR_VECTOR
        else if (present(ABSERR)) then
          if (ABSERR<ZERO) then
            MSG = 'ABSERR must be nonnegative.'
            call XERRDV(MSG,800,2,0,0,0,0,ZERO,ZERO)
          end if
          OPTS%ATOL = ABSERR
        else
          if (ALLOW_DEFAULT_TOLS) then
             OPTS%ATOL = 1D-6
             MSG = 'By not specifying ABSERR, you have elected to use a default'
             call XERRDV(MSG,810,1,0,0,0,0,ZERO,ZERO)
             MSG = 'absolute error tolerance equal to 1.0D-6. Please be aware a'
             call XERRDV(MSG,810,1,0,0,0,0,ZERO,ZERO)
             MSG = 'tolerance this large is not appropriate for all problems.'
             call XERRDV(MSG,810,1,0,0,0,0,ZERO,ZERO)
             MSG = 'Execution will continue'
             call XERRDV(MSG,810,1,0,0,0,0,ZERO,ZERO)
          else
             MSG = 'You must specify a vector of absolute error tolerances or'
             call XERRDV(MSG,820,1,0,0,0,0,ZERO,ZERO)
             MSG = 'a scalar error tolerance. It is recommended that you use'
             call XERRDV(MSG,820,1,0,0,0,0,ZERO,ZERO)
             MSG = 'a vector of absolute error tolerances.'
             call XERRDV(MSG,820,2,0,0,0,0,ZERO,ZERO)
          end if
        end if

!       ITOL error tolerance flag.
!          ITOL   RTOL     ATOL            EWT(i)
!            1   scalar   scalar  RTOL*ABS(Y(i)) + ATOL
!            2   scalar   array   RTOL*ABS(Y(i)) + ATOL(i)
!            3   array    scalar  RTOL(i)*ABS(Y(i)) + ATOL
!            4   array    array   RTOL(i)*ABS(Y(i)) + ATOL(i)
        if (present(ABSERR_VECTOR)) then
          if (present(RELERR_VECTOR)) then
            OPTS%ITOL = 4
          else
            OPTS%ITOL = 2
          end if
        else
          if (present(RELERR_VECTOR)) then
            OPTS%ITOL = 3
          else
            OPTS%ITOL = 1
          end if
        end if
        return

  end function SET_OPTS
!_______________________________________________________________________

      subroutine GET_STATS(RSTATS,ISTATS,NUMEVENTS,JROOTS)
! ..
! Return the user portions of the DVODE RUSER and IUSER arrays;
! and if root finding is being done, return the JROOT vector
! (not called by DVODE_F90).
! ..
! Available Integration Statistics.
! HU      RUSER(11) The step size in t last used (successfully).
! HCUR    RUSER(12) The step size to be attempted on the next step.
! TCUR    RUSER(13) The current value of the independent variable
!                   which the solver has actually reached, i.e. the
!                   current internal mesh point in t. In the output,
!                   TCUR will always be at least as far from the
!                   initial value of t as the current argument T,
!                   but may be farther (if interpolation was done).
! TOLSF   RUSER(14) A tolerance scale factor, greater than 1.0,
!                   computed when a request for too much accuracy was
!                   detected (ISTATE = -3 if detected at the start of
!                   the problem, ISTATE = -2 otherwise). If ITOL is
!                   left unaltered but RTOL and ATOL are uniformly
!                   scaled up by a factor of TOLSF for the next call,
!                   then the solver is deemed likely to succeed.
!                   (The user may also ignore TOLSF and alter the
!                   tolerance parameters in any other way appropriate.)
! NST     IUSER(11) The number of steps taken for the problem so far.
! NFE     IUSER(12) The number of f evaluations for the problem so far.
! NJE     IUSER(13) The number of Jacobian evaluations so far.
! NQU     IUSER(14) The method order last used (successfully).
! NQCUR   IUSER(15) The order to be attempted on the next step.
! IMXER   IUSER(16) The index of the component of largest magnitude in
!                   the weighted local error vector (e(i)/EWT(i)),
!                   on an error return with ISTATE = -4 or -5.
! LENRW   IUSER(17) The length of RUSER actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
! LENIW   IUSER(18) The length of IUSER actually required.
!                   This is defined on normal returns and on an illegal
!                   input return for insufficient storage.
! NLU     IUSER(19) The number of matrix LU decompositions so far.
! NNI     IUSER(20) The number of nonlinear (Newton) iterations so far.
! NCFN    IUSER(21) The number of convergence failures of the nonlinear
!                   solver so far.
! NETF    IUSER(22) The number of error test failures of the integrator
!                   so far.
! MA28AD_CALLS      IUSER(23) The number of calls made to MA28AD
! MA28BD_CALLS      IUSER(24) The number of calls made to MA28BD
! MA28CD_CALLS      IUSER(25) The number of calls made to MA28CD
! MC19AD_CALLS      IUSER(26) The number of calls made to MC19AD
! IRNCP             IUSER(27) The number of compressions done on array JAN
! ICNCP             IUSER(28) The number of compressions done on array ICN
! MINIRN            IUSER(29) Minimum size for JAN array
! MINICN            IUSER(30) Minimum size for ICN array
! JROOTS  JROOTS    Optional array of component indices for components
!                   having a zero at the current time
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, optional, intent (IN) :: NUMEVENTS
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: RSTATS(22)
        integer, intent (INOUT) :: ISTATS(31)
        integer, optional, intent (INOUT) :: JROOTS(:)
! ..
! .. Local Scalars ..
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
        intrinsic ALLOCATED, PRESENT, SIZE
! ..
!   Caution:
!   It is assumed that the size of RSTATS is 22.
!   It is assumed that the size of ISTATS is 31.
!   Refer to the documentation prologue for a description
!   the optional output contained in RUSER and IUSER.
! ..
! .. FIRST EXECUTABLE STATEMENT GET_STATS
! ..
!       Check if DVODE_F90 has been called yet.
        if (.not.OPTS_CALLED) then
          MSG = 'You have called GET_STATS before'
          call XERRDV(MSG,830,1,0,0,0,0,ZERO,ZERO)
          MSG = 'calling DVODE_F90 the first time.'
          call XERRDV(MSG,830,1,0,0,0,0,ZERO,ZERO)
          return
        end if

!       Check that the arrays are large enough to hold the statistics.
!       Some compilers don't like this use of SIZE.
!       IF (SIZE(RSTATS)<LRWUSER) THEN
!       IF (SIZE(RSTATS)<31) THEN
!         MSG = 'In GET_STATS, RSTATS array is too small.'
!         CALL XERRDV(MSG,840,1,0,0,0,0,ZERO,ZERO)
!         RETURN
!       END IF
!       IF (SIZE(ISTATS)<LIWUSER) THEN
!       IF (SIZE(ISTATS)<31) THEN
!         MSG = 'In GET_STATS, ISTATS array is too small.'
!         CALL XERRDV(MSG,850,1,0,0,0,0,ZERO,ZERO)
!         RETURN
!       END IF

!       Copy the statistics.
        RSTATS(1:LRWUSER) = RUSER(1:LRWUSER)
!       ISTATS(1:LIWUSER) = IUSER(1:LIWUSER)
        ISTATS(1:22) = IUSER(1:22)
        ISTATS(23) = MA28AD_CALLS
        ISTATS(24) = MA28BD_CALLS
        ISTATS(25) = MA28CD_CALLS
!_______________________________________________________________________
! *****MA48 build change point.Replace the three previous statements
!      with these statements.
!       IF (USE_MA48_FOR_SPARSE) THEN
!          ISTATS(23) = MA48AD_CALLS
!          ISTATS(24) = MA48BD_CALLS
!          ISTATS(25) = MA48CD_CALLS
!       ELSE
!          ISTATS(23) = MA28AD_CALLS
!          ISTATS(24) = MA28BD_CALLS
!          ISTATS(25) = MA28CD_CALLS
!       END IF
!_______________________________________________________________________
        ISTATS(26) = MC19AD_CALLS
        ISTATS(27) = IRNCP
        ISTATS(28) = ICNCP
!       ISTATS(29) = MINIRN
!       ISTATS(30) = MINICN
!       ISTATS(31) = NZ_ALL
        ISTATS(29) = MAX_MINIRN
        ISTATS(30) = MAX_MINICN
        ISTATS(31) = MAX_NNZ
!       If root finding is being used return the JROOT vector.
        if (present(NUMEVENTS)) then
          if (present(JROOTS)) then
            if (allocated(JROOT)) then
              JROOTS(1:NUMEVENTS) = JROOT(1:NUMEVENTS)
            end if
          end if
        end if
        return

      end subroutine GET_STATS
!_______________________________________________________________________

      subroutine USERSETS_IAJA(IAUSER,NIAUSER,JAUSER,NJAUSER)
! ..
! Approximate or allow the user to supply the sparse Jacobian
! structure pointer arrays IA and JA directly for DVODE_F90
! (not called by DVODE_F90).
! ..
! Used if the user wishes to supply the sparsity structure
! directly.
!     Caution:
!     If it is called, USERSETS_IAJA must be called after the
!     call to SET_OPTS.
!     Usage:
!     CALL SET_IAJA(IAUSER,NIAUSER,JAUSER,NJAUSER)
!       In this case, IAUSER of length NIAUSER will be used for
!       IA; and JAUSER of length NJAUSER will be used for JA.
!     Arguments:
!     IAUSER  = user supplied IA array
!     NIAUSER = length of IAUSER array
!     JAUSER  = user supplied JA vector
!     NJAUSER = length of JAUSER array
!     Results:
!     IA(IADIM), IADIM, JA(JADIM), JAMIN, SPARSE
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: NIAUSER, NJAUSER
! ..
! .. Array Arguments ..
        integer, intent (IN) :: IAUSER(NIAUSER), JAUSER(NJAUSER)
! ..
! .. Local Scalars ..
        integer :: JER
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
        intrinsic ALLOCATED
! ..
! .. FIRST EXECUTABLE STATEMENT USERSETS_IAJA
! ..
!       Check that SET_OPTS has been called:
        if (.not.OPTS_CALLED) then
          MSG = 'You have not called SET_OPTS before'
          call XERRDV(MSG,860,1,0,0,0,0,ZERO,ZERO)
          MSG = 'calling USERSETS_IAJA.'
          call XERRDV(MSG,860,2,0,0,0,0,ZERO,ZERO)
        end if
        SPARSE = .false.
        IAJA_CALLED = .false.
        IADIM = NIAUSER
        JADIM = NJAUSER
        if (allocated(IA)) then
          deallocate (IA,JA,STAT=JER)
          call CHECK_STAT(JER,170)
        end if
        allocate (IA(IADIM),JA(JADIM),STAT=JER)
        if (JER/=0) then
          MSG = 'Check your values of NIAUSER and NJAUSER'
          call XERRDV(MSG,870,1,0,0,0,0,ZERO,ZERO)
          MSG = 'in your call to USERSETS_IAJA.'
          call XERRDV(MSG,870,2,0,0,0,0,ZERO,ZERO)
        end if
        call CHECK_STAT(JER,180)
        IA(1:IADIM) = IAUSER(1:IADIM)
        JA(1:JADIM) = JAUSER(1:JADIM)
!       Set the flags to indicate to DVODE_F90 that IA and JA have
!       been loaded successfully.
        SPARSE = .true.
        IAJA_CALLED = .true.
        return

      end subroutine USERSETS_IAJA
!_______________________________________________________________________

      subroutine SET_IAJA(DFN,NEQ,T,Y,FMIN,NTURB,DTURB,IAUSER,NIAUSER, &
        JAUSER,NJAUSER)
! ..
!     Approximate or allow the user to supply the sparse Jacobian
!     structure pointer arrays IA and JA for DVODE_F90 (not called
!     by DVODE_F90).
! ..
!     Caution:
!     If it is called, SET_IAJA must be called after the call to
!     SET_OPTS.
!     Usage:
!     SET_IAJA may be called in one of two ways:
!     CALL SET_IAJA(DFN,NEQ,T,Y,FMIN,NTURB,DTURB)
!       In this case IA and JA will be determined using calls
!       to your derivative routine DFN.
!     CALL SET_IAJA(DFN,NEQ,T,Y,FMIN,NTURB,DTURB,IAUSER,NIAUSER, &
!       JAUSER,NJAUSER)
!       In this case, IAUSER of length NIAUSER will be used for
!       IA; and JAUSER of length NJAUSER will be used for JA.
!       T, Y, FMIN, NTURB, and DTURB will be ignored (though
!       they must be present in the argument list).
!
!     Arguments:
!     DFN     = VODE derivative subroutine
!     NEQ     = Number of odes
!     T       = independent variable t
!     Y       = solution y(t)
!     FMIN    = Jacobian threshold value. Elements of the Jacobian
!               with magnitude smaller than FMIN will be ignored.
!               FMIN will be ignored if it is less than or equal
!               to ZERO.
!     NTURB   = Perturbation flag. If NTURB=1, component I of Y
!               will be perturbed by by 1.01D0.
!               If NTURB=NEQ, component I of Y will be perturbed
!               by ONE + DTURB(I).
!     DTURB   = perturbation vector of length 1 or NEQ.
!     If these four optional parameters are present, IAUSER and JAUSER
!     will be copied to IA and JA rather than making derivative calls
!     to approximate IA and JA:
!     IAUSER  = user supplied IA array
!     NIAUSER = length of IAUSER array
!     JAUSER  = user supplied JA vector
!     NJAUSER = length of JAUSER array
!     Results:
!     IA(IADIM), IADIM, IMIN, JA(JADIM), JAMIN, JMIN, SPARSE
!     Allocated but free for further use by DVODE_F90:
!     FTEMP(NEQ), FPTEMP(NEQ)
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP), intent (IN) :: FMIN, T
        integer, intent (IN) :: NEQ, NTURB
        integer, optional, intent (IN) :: NIAUSER, NJAUSER
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: DTURB(*), Y(NEQ)
        integer, optional, intent (IN) :: IAUSER(:), JAUSER(:)
! ..
! .. Subroutine Arguments ..
        external DFN
! ..
! .. Local Scalars ..
        real (WP) :: AIJ, DTRB, EMIN, RJ, YJ, YJSAVE, YPJ
        integer :: I, J, JER, JP1, K, MTURB
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, ALLOCATED, EPSILON, MAX, MIN, PRESENT
! ..
! .. FIRST EXECUTABLE STATEMENT SET_IAJA
! ..
        SPARSE = .false.
        IAJA_CALLED = .false.
        if (NEQ<1) goto 30

!       Check if the user is supplying the IA and JA arrays.
        if (present(IAUSER)) then
!         If IAUSER is present, so must be NIAUSER, JAUSER,
!         and NJAUSER.
          if (.not.(present(NIAUSER))) goto 30
          if (.not.(present(JAUSER))) goto 30
          if (.not.(present(NJAUSER))) goto 30
          if (NIAUSER<NEQ+1 .or. NJAUSER<1) goto 30
          IADIM = min(NIAUSER,NEQ+1)
          IMIN = IADIM
          JADIM = NJAUSER
          JMIN = NJAUSER
          if (allocated(IA)) then
            deallocate (IA,JA,FTEMP,FPTEMP,STAT=JER)
            call CHECK_STAT(JER,190)
          end if
          allocate (IA(IADIM),JA(JADIM),FTEMP(NEQ),FPTEMP(NEQ),STAT=JER)
          call CHECK_STAT(JER,200)
          IA(1:IADIM) = IAUSER(1:IADIM)
          JA(1:JADIM) = JAUSER(1:JADIM)
          goto 40
        end if

!       Determine IA and JA using derivative calls.

!       Jacobian element magnitude threshold value.
        if (FMIN>ZERO) then
          EMIN = max(FMIN,ZERO)
        else
          EMIN = ZERO
        end if

!       Solution perturbation factors.
        if (NTURB>0) then
          if (NTURB<1 .or. (NTURB>1 .and. NTURB/=NEQ)) goto 30
          if (NTURB==NEQ) then
            do I = 1, NEQ
              DTURB(I) = max(DTURB(I),HUNDRETH)
            end do
          else
            MTURB = 1
            DTRB = max(DTURB(1),HUNDRETH)
          end if
        else
          MTURB = 1
          DTRB = HUNDRETH
        end if

        JADIM = min(NEQ*NEQ,max(1000,NZ_SWAG))
        ADDTOJA = max(1000,NZ_SWAG)
!       Loop point for array allocation.
10      continue
        if (allocated(IA)) then
          deallocate (IA,JA,FTEMP,FPTEMP,STAT=JER)
          call CHECK_STAT(JER,210)
        end if
        IMIN = NEQ + 1
        IADIM = NEQ + 1
        JMIN = 0
        JADIM = JADIM + ADDTOJA
        if (JADIM>MAX_ARRAY_SIZE) then
          MSG = 'Maximum array size exceeded. Stopping.'
          call XERRDV(MSG,880,2,0,0,0,0,ZERO,ZERO)
        end if
        allocate (IA(IADIM),JA(JADIM),FTEMP(NEQ),FPTEMP(NEQ),STAT=JER)
        call CHECK_STAT(JER,220)

!       f = y'(t,y).
        call DFN(NEQ,T,Y,FTEMP)
        IA(1) = 1
        K = 1

!       Calculate unit roundoff and powers of it used if JACSP
!       is used.
        UROUND = epsilon(ONE)

!       Successively perturb each of the solution components and
!       calculate the corresponding derivatives.
        do J = 1, NEQ
          if (MTURB==NEQ) DTRB = DTURB(J)
          YJ = Y(J)
          YJSAVE = YJ
          if (abs(YJ)<=ZERO) YJ = HUN*UROUND
          YPJ = YJ*(ONE+DTRB)
          RJ = abs(YPJ-YJ)
          if (abs(RJ)<=ZERO) RJ = HUN*UROUND
          if (abs(RJ)<=ZERO) goto 30
          Y(J) = YPJ
          call DFN(NEQ,T,Y,FPTEMP)
          do 20 I = 1, NEQ
!           Estimate the Jacobian element.
            AIJ = abs(FPTEMP(I)-FTEMP(I))/RJ
            if ((AIJ<=EMIN) .and. (I/=J)) goto 20
!           Need more storage for JA.
            if (K>JADIM) goto 10
            JMIN = K
            JA(K) = I
            K = K + 1
20        continue
          JP1 = J + 1
          IA(JP1) = K
          Y(J) = YJSAVE
        end do
        goto 40
30      continue
        MSG = 'An error occurred in subroutine SET_IAJA.'
        call XERRDV(MSG,890,2,0,0,0,0,ZERO,ZERO)

40      continue
!       Set the flags to indicate to DVODE_F90 that IA and JA have
!       been calculated successfully.
        SPARSE = .true.
        IAJA_CALLED = .true.

!       Trim JA.
        allocate (JATEMP(JMIN),STAT=JER)
        call CHECK_STAT(JER,230)
        JATEMP(1:JMIN) = JA(1:JMIN)
        deallocate (JA,STAT=JER)
        call CHECK_STAT(JER,240)
        allocate (JA(JMIN),STAT=JER)
        call CHECK_STAT(JER,250)
        JA(1:JMIN) = JATEMP(1:JMIN)
        deallocate (JATEMP,STAT=JER)
        call CHECK_STAT(JER,260)
        JADIM = JMIN
        return

      end subroutine SET_IAJA
!_______________________________________________________________________

      subroutine DVODE(F,NEQ,Y,T,TOUT,ITASK,ISTATE,OPTS,JAC,GFUN)
! ..
! This is the core driver (modified original DVODE.F driver).
! ..
! The documentation prologue was moved nearer the top of the file.
! ..
     implicit none
! ..
! .. Structure Arguments ..
        type (VODE_OPTS) :: OPTS
! ..
! .. Scalar Arguments ..
        real (WP), intent (INOUT) :: T, TOUT
        integer, intent (INOUT) :: ISTATE
        integer, intent (IN) :: ITASK, NEQ
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: Y(*)
! ..
! .. Subroutine Arguments ..
        external F, GFUN, JAC
! ..
! .. Local Scalars ..
        real (WP) :: ATOLI, BIG, EWTI, H0, HMAX, HMX, RH, RTOLI, SIZEST, &
          TCRIT, TNEXT, TOLSF, TP
        integer :: I, IER, IFLAG, IMXER, IOPT, IPCUTH, IRFP, IRT, ITOL, JCO, &
          JER, KGO, LENIW, LENJ, LENP, LENRW, LENWM, LF0, MBAND, MF, MFA, ML, &
          MU, NG, NITER, NSLAST
        logical :: IHIT
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, ALLOCATED, EPSILON, MAX, MIN, SIGN, SQRT
! ..
! The following internal PRIVATE variable blocks contain variables which
! are communicated between subroutines in the DVODE package, or which
! are to be saved between calls to DVODE.
! In each block, real variables precede integers.
! The variables stored in the internal PRIVATE variable blocks are as
! follows:
! ACNRM  = Weighted r.m.s. norm of accumulated correction vectors.
! CCMXJ  = Threshhold on DRC for updating the Jacobian. (See DRC.)
! CONP   = The saved value of TQ(5).
! CRATE  = Estimated corrector convergence rate constant.
! DRC    = Relative change in H*RL1 since last DVJAC call.
! EL     = Real array of integration coefficients. See DVSET.
! ETA    = Saved tentative ratio of new to old H.
! ETAMAX = Saved maximum value of ETA to be allowed.
! H      = The step size.
! HMIN   = The minimum absolute value of the step size H to be used.
! HMXI   = Inverse of the maximum absolute value of H to be used.
!          HMXI = 0.0 is allowed and corresponds to an infinite HMAX.
! HNEW   = The step size to be attempted on the next step.
! HSCAL  = Stepsize in scaling of YH array.
! PRL1   = The saved value of RL1.
! RC     = Ratio of current H*RL1 to value on last DVJAC call.
! RL1    = The reciprocal of the coefficient EL(2).
! TAU    = Real vector of past NQ step sizes, length 13.
! TQ     = A real vector of length 5 in which DVSET stores constants
!          used for the convergence test, the error test, and the
!          selection of H at a new order.
! TN     = The independent variable, updated on each step taken.
! UROUND = The machine unit roundoff. The smallest positive real number
!          such that  1.0 + UROUND /= 1.0
! ICF    = Integer flag for convergence failure in DVNLSD:
!            0 means no failures.
!            1 means convergence failure with out of date Jacobian
!                   (recoverable error).
!            2 means convergence failure with current Jacobian or
!                   singular matrix (unrecoverable error).
! INIT   = Saved integer flag indicating whether initialization of the
!          problem has been done (INIT = 1) or not.
! IPUP   = Saved flag to signal updating of Newton matrix.
! JCUR   = Output flag from DVJAC showing Jacobian status:
!            JCUR = 0 means J is not current.
!            JCUR = 1 means J is current.
! JSTART = Integer flag used as input to DVSTEP:
!            0  means perform the first step.
!            1  means take a new step continuing from the last.
!           -1  means take the next step with a new value of MAXORD,
!               HMIN, HMXI, N, METH, MITER, and/or matrix parameters.
!          On return, DVSTEP sets JSTART = 1.
! JSV    = Integer flag for Jacobian saving, = sign(MF).
! KFLAG  = A completion code from DVSTEP with the following meanings:
!               0      the step was successful.
!              -1      the requested error could not be achieved.
!              -2      corrector convergence could not be achieved.
!              -3, -4  fatal error in DVNLSD(can not occur here).
! KUTH   = Input flag to DVSTEP showing whether H was reduced by the
!          driver. KUTH = 1 if H was reduced, = 0 otherwise.
! L      = Integer variable, NQ + 1, current order plus one.
! LMAX   = MAXORD + 1 (used for dimensioning).
! LOCJS  = A pointer to the saved Jacobian, whose storage starts at
!          WM(LOCJS), if JSV = 1.
! LYH    = Saved integer pointer to segments of RWORK and IWORK.
! MAXORD = The maximum order of integration method to be allowed.
! METH/MITER = The method flags. See MF.
! MSBJ   = The maximum number of steps between J evaluations, = 50.
! MXHNIL = Saved value of optional input MXHNIL.
! MXSTEP = Saved value of optional input MXSTEP.
! N      = The number of first-order ODEs, = NEQ.
! NEWH   = Saved integer to flag change of H.
! NEWQ   = The method order to be used on the next step.
! NHNIL  = Saved counter for occurrences of T + H = T.
! NQ     = Integer variable, the current integration method order.
! NQNYH  = Saved value of NQ*NYH.
! NQWAIT = A counter controlling the frequency of order changes.
!          An order change is about to be considered if NQWAIT = 1.
! NSLJ   = The number of steps taken as of the last Jacobian update.
! NSLP   = Saved value of NST as of last Newton matrix update.
! NYH    = Saved value of the initial value of NEQ.
! HU     = The step size in t last used.
! NCFN   = Number of nonlinear convergence failures so far.
! NETF   = The number of error test failures of the integrator so far.
! NFE    = The number of f evaluations for the problem so far.
! NJE    = The number of Jacobian evaluations so far.
! NLU    = The number of matrix LU decompositions so far.
! NNI    = Number of nonlinear iterations so far.
! NQU    = The method order last used.
! NST    = The number of steps taken for the problem so far.
! Block 0.
! Retrieve the necessary flags from the OPTIONS structure and manage
! the storage allocation.
! ..
! .. FIRST EXECUTABLE STATEMENT DVODE
! ..
!       Retrieve the local flags from the options structure.
        IOPT = OPTS%IOPT
        ITOL = OPTS%ITOL
        MF = OPTS%MF
        METH = OPTS%METH
        MITER = OPTS%MITER
        MOSS = OPTS%MOSS
        NG = OPTS%NG

!       Allocate the necessary storage for RWORK and IWORK. (Assume that
!       both or neither of the arrays are allocated.)

!       If we are starting a new problem, deallocate the old RWORK and
!       IWORK arrays if they were allocated in a previous call. Also
!       manage the event residual arrays.

        if (ISTATE==1) then
          if (allocated(DTEMP)) then
            deallocate (DTEMP,YTEMP,STAT=JER)
            call CHECK_STAT(JER,270)
          end if
          if (allocated(RWORK)) then
            deallocate (RWORK,IWORK,ACOR,SAVF,EWT,WM,STAT=JER)
            call CHECK_STAT(JER,280)
          end if
          if (allocated(JROOT)) then
            deallocate (JROOT,G0,G1,GX,STAT=JER)
            call CHECK_STAT(JER,290)
          end if
          if (allocated(YMAX)) then
            deallocate (YMAX,STAT=JER)
            call CHECK_STAT(JER,300)
          end if
        end if

!       If the user has made changes and called DVODE_F90 with ISTATE=3,
!       temporarily save the necessary portion of the Nordsieck array
!       and then release RWORK and IWORK. Also, save the contents of
!       the WM array and release WM. Note: ACOR, SAVF, and EWT do not
!       need to be saved.

        if (ISTATE==3) then
          if (IOPT/=1) then
            MAXORD = MORD(METH)
          else
            MAXORD = IUSER(5)
            if (MAXORD<0) goto 520
            if (MAXORD==0) MAXORD = 100
            MAXORD = min(MAXORD,MORD(METH))
          end if
          if (MAXORD<NQ) then
!           If MAXORD is less than NQ, save column NQ+2 of YH for
!           later use.
            if (allocated(YHNQP2)) then
              deallocate (YHNQP2,STAT=JER)
              call CHECK_STAT(JER,310)
            end if
            allocate (YHNQP2(NEQ),STAT=JER)
            call CHECK_STAT(JER,320)
            call DCOPY_F90(NEQ,RWORK(LYH+NYH*(MAXORD+1)),1,YHNQP2,1)
          end if
          LYHTEMP = min((PREVIOUS_MAXORD+1)*NYH+LYH-1,(MAXORD+1)*NYH+LYH-1)
          allocate (YHTEMP(LYHTEMP),STAT=JER)
          call CHECK_STAT(JER,330)
!         Save YH.
          YHTEMP(1:LYHTEMP) = RWORK(1:LYHTEMP)
!         Save WM.
          LWMTEMP = LWMDIM
          allocate (WMTEMP(LWMTEMP),STAT=JER)
          call CHECK_STAT(JER,340)
          WMTEMP(1:LWMTEMP) = WM(1:LWMTEMP)
          if (allocated(DTEMP)) then
             deallocate (DTEMP,YTEMP,STAT=JER)
             call CHECK_STAT(JER,350)
          end if
          deallocate (RWORK,IWORK,YMAX,ACOR,SAVF,EWT,WM,STAT=JER)
          call CHECK_STAT(JER,360)
        end if

!     Allocate the RWORK and WM work arrays if they haven't already
!     been allocated in a previous call.
!     LRW = 20 + (MAXORD + 1) * NEQ
!     LWMDIM =     
!            0                        for MF = 10,
!            2 * NEQ**2               for MF = 11 or 12,
!            NEQ**2                   for MF = -11 or -12,
!            NEQ                      for MF = 13,
!            (3*ML + 2*MU + 2) * NEQ  for MF = 14 or 15,
!            (2*ML + MU + 1) * NEQ    for MF = -14 or -15,
!            0                        for MF = 16 or 17,
!            0                        for MF = -16 or -17,
!            0                        for MF = 20,
!            2 * NEQ**2               for MF = 21 or 22,
!            NEQ**2                   for MF = -21 or -22,
!            NEQ                      for MF = 23,
!            (3*ML + 2*MU + 2) * NEQ  for MF = 24 or 25,
!            (2*ML + MU + 1) * NEQ    for MF = -24 or -25,
!            0                        for MF = 26 or 27,
!            0                        for MF = -26 or -27.
        if (allocated(RWORK)) then
        else
          MFA = abs(MF)
          if (IOPT/=1) then
            MAXORD = MORD(METH)
          else
            MAXORD = IUSER(5)
            if (MAXORD<0) goto 520
            if (MAXORD==0) MAXORD = 100
            MAXORD = min(MAXORD,MORD(METH))
          end if
          LRW = LRWUSER + (MAXORD+1)*NEQ
!        IF (MITER == 0) LRW = LRW - 2
          LRW = LRW - 2
          if (MF==10) then
            LWMDIM = 0
          else if (MF==11 .or. MF==12) then
            LWMDIM = 2*NEQ**2
          else if (MF==-11 .or. MF==-12) then
            LWMDIM = NEQ**2
          else if (MF==13) then
            LWMDIM = NEQ
          else if (MF==14 .or. MF==15) then
            ML = IUSER(1)
            MU = IUSER(2)
            if (ML<0 .or. ML>=NEQ) goto 500
            if (MU<0 .or. MU>=NEQ) goto 510
            LWMDIM = (3*ML+2*MU+2)*NEQ
          else if (MF==-14 .or. MF==-15) then
            ML = IUSER(1)
            MU = IUSER(2)
            if (ML<0 .or. ML>=NEQ) goto 500
            if (MU<0 .or. MU>=NEQ) goto 510
            LWMDIM = (2*ML+MU+1)*NEQ
          else if (MF==16 .or. MF==17) then
            LWMDIM = 0
          else if (MF==-16 .or. MF==-17) then
            LWMDIM = 0
          else if (MF==20) then
            LWMDIM = 0
          else if (MF==21 .or. MF==22) then
            LWMDIM = 2*NEQ**2
          else if (MF==-21 .or. MF==-22) then
            LWMDIM = NEQ**2
          else if (MF==23) then
            LWMDIM = NEQ
          else if (MF==24 .or. MF==25) then
            ML = IUSER(1)
            MU = IUSER(2)
            if (ML<0 .or. ML>=NEQ) goto 500
            if (MU<0 .or. MU>=NEQ) goto 510
            LWMDIM = (3*ML+2*MU+2)*NEQ
          else if (MF==-24 .or. MF==-25) then
            ML = IUSER(1)
            MU = IUSER(2)
            if (ML<0 .or. ML>=NEQ) goto 500
            if (MU<0 .or. MU>=NEQ) goto 510
            LWMDIM = (2*ML+MU+1)*NEQ
          else if (MF==26 .or. MF==27) then
            LWMDIM = 0
          else if (MF==-26 .or. MF==-27) then
            LWMDIM = 0
          end if
!         LWMDIM = LWMDIM + 2
          allocate (RWORK(LRW),WM(LWMDIM),STAT=JER)
          call CHECK_STAT(JER,370)
          RWORK(1:LRW) = ZERO
          WM(1:LWMDIM) = ZERO
          if (.not.allocated(DTEMP)) then
             allocate (DTEMP(NEQ),YTEMP(NEQ),STAT=JER)
             call CHECK_STAT(JER,380)
          end if
          allocate (ACOR(NEQ),SAVF(NEQ),EWT(NEQ),STAT=JER)
          call CHECK_STAT(JER,390)
!         If necessary, reload the saved portion of the Nordsieck
!         array and the WM array, saved above.
          if (ISTATE==3) then
            RWORK(1:LYHTEMP) = YHTEMP(1:LYHTEMP)
            I = min(LWMTEMP,LWMDIM)
!           WM(1:LWMTEMP) = WMTEMP(1:LWMTEMP)
            WM(1:I) = WMTEMP(1:I)
            deallocate (YHTEMP,WMTEMP,STAT=JER)
            call CHECK_STAT(JER,400)
          end if
          RWORK(1:LRWUSER) = RUSER(1:LRWUSER)
!         IUSER: = LIWUSER if MITER = 0 or 3 (MF = 10, 13, 20, 23)
!                = LIWUSER + NEQ otherwise
!                  (ABS(MF) = 11,12,14,15,16,17,21,22,24,25,26,27).
          LIW = LIWUSER + NEQ
          if (MITER==0 .or. MITER==3) LIW = LIWUSER
          allocate (IWORK(LIW),STAT=JER)
          call CHECK_STAT(JER,410)
          IWORK(1:LIWUSER) = IUSER(1:LIWUSER)
!         Allocate the YMAX vector.
          allocate (YMAX(NEQ),STAT=JER)
          call CHECK_STAT(JER,420)
        end if
!       Allocate the event arrays if they haven't already been allocated
!       in a previous call.
        if (allocated(JROOT)) then
        else
          if (NG>0) then
            allocate (JROOT(NG),G0(NG),G1(NG),GX(NG),STAT=JER)
            call CHECK_STAT(JER,430)
          end if
        end if

! Block A.
! This code block is executed on every call. It tests ISTATE and
! ITASK for legality and branches appropriately.
! If ISTATE > 1 but the flag INIT shows that initialization has
! not yet been done, an error return occurs.
! If ISTATE = 1 and TOUT = T, return immediately.
! The user portion of RWORK and IWORK are reloaded since something
! may have changed since the previous visit.

        RWORK(1:LRWUSER) = RUSER(1:LRWUSER)
        IWORK(1:LIWUSER) = IUSER(1:LIWUSER)
        if (ISTATE<1 .or. ISTATE>3) goto 420
        if (ITASK<1 .or. ITASK>5) goto 430
        ITASKC = ITASK
        if (ISTATE==1) goto 10
        if (INIT/=1) goto 440
        if (ISTATE==2) goto 130
        goto 20
10      INIT = 0
        if (abs(TOUT-T)<=ZERO) then
          RUSER(1:LRWUSER) = RWORK(1:LRWUSER)
          IUSER(1:LIWUSER) = IWORK(1:LIWUSER)
          return
        end if

! Block B.
! The next code block is executed for the initial call (ISTATE = 1),
! or for a continuation call with parameter changes (ISTATE = 3).
! It contains checking of all input and various initializations.

!       First check legality of the non-optional input NEQ, ITOL, IOPT,
!       MF, ML, and MU.

20      if (NEQ<=0) goto 450
        if (ISTATE==1) goto 30
        if (NEQ>N) goto 460
        if (NEQ/=N) goto 465
30      N = NEQ
        if (ITOL<1 .or. ITOL>4) goto 470
        if (IOPT<0 .or. IOPT>1) goto 480
        JSV = sign(1,MF)
        INEWJ = (1-JSV)/2
        MFA = abs(MF)
        METH = MFA/10
        MITER = MFA - 10*METH
        if (METH<1 .or. METH>2) goto 490
        if (MITER<0 .or. MITER>7) goto 490
        if (MITER<=3) goto 40
        if (MITER<=5) then
          ML = IWORK(1)
          MU = IWORK(2)
          if (ML<0 .or. ML>=N) goto 500
          if (MU<0 .or. MU>=N) goto 510
        end if
40      continue
        if (NG<0) goto 700
        if (ISTATE==1) goto 50
        if (IRFND==0 .and. NG/=NGC) goto 710
50      NGC = NG
!       Next process and check the optional input.
        if (IOPT==1) goto 60
        MAXORD = MORD(METH)
        MXSTEP = MXSTP0
        MXHNIL = MXHNL0
        if (ISTATE==1) H0 = ZERO
        HMXI = ZERO
        HMIN = ZERO
        goto 80
60      MAXORD = IWORK(5)
        if (MAXORD<0) goto 520
        if (MAXORD==0) MAXORD = 100
        MAXORD = min(MAXORD,MORD(METH))
        MXSTEP = IWORK(6)
        if (MXSTEP<0) goto 530
        if (MXSTEP==0) MXSTEP = MXSTP0
        MXHNIL = IWORK(7)
        if (MXHNIL<0) goto 540
        if (MXHNIL==0) MXHNIL = MXHNL0
        if (ISTATE/=1) goto 70
        H0 = RWORK(5)
        if ((TOUT-T)*H0<ZERO) goto 550
70      HMAX = RWORK(6)
        if (HMAX<ZERO) goto 560
        HMXI = ZERO
        if (HMAX>ZERO) HMXI = ONE/HMAX
        HMIN = RWORK(7)
        if (HMIN<ZERO) goto 570
        SETH = RWORK(8)
        if (SETH<ZERO) goto 690
!       Check the nonnegativity information.
        if (BOUNDS) then
          if (NDX<1 .or. NDX>NEQ) then
            MSG = 'The size of the CONSTRAINED vector'
            call XERRDV(MSG,900,1,0,0,0,0,ZERO,ZERO)
            MSG = 'must be between 1 and NEQ.'
            call XERRDV(MSG,900,2,0,0,0,0,ZERO,ZERO)
          end if
          do I = 1, NDX
            if (IDX(I)<1 .or. IDX(I)>N) then
              MSG = 'Each component of THE CONSTRAINED'
              call XERRDV(MSG,910,1,0,0,0,0,ZERO,ZERO)
              MSG = 'vector must be between 1 and N.'
              call XERRDV(MSG,910,2,0,0,0,0,ZERO,ZERO)
            end if
          end do
        end if
!       Check the sub diagonal and super diagonal arrays.
        if (MITER==4 .or. MITER==5) then
           if (SUBS) then
              do I = 1, NSUBS
                 if (SUBDS(I) < 2 .or. SUBDS(I) > ML+1) then
                    MSG = 'Each element of SUB_DIAGONALS'
                    call XERRDV(MSG,920,1,0,0,0,0,ZERO,ZERO)
                    MSG = 'must be between 2 and ML + 1.'
                    call XERRDV(MSG,920,2,0,0,0,0,ZERO,ZERO)
                 end if
              end do
           end if
           if (SUPS) then
              do I = 1, NSUPS
                 if (SUPDS(I) < 2 .or. SUPDS(I) > MU + 1) then
                    MSG = 'Each element of SUP_DIAGONALS'
                    call XERRDV(MSG,930,1,0,0,0,0,ZERO,ZERO)
                    MSG = 'must be between 2 and MU + 1.'
                    call XERRDV(MSG,930,2,0,0,0,0,ZERO,ZERO)
                 end if
              end do
           end if
!          Compute the banded column grouping.
           if (SUBS .or. SUPS) then
             call BGROUP(N,EWT,ACOR,YMAX,ML,MU)
             BUILD_IAJA = .false.
             BUILD_IAJA = .true.
             if (BUILD_IAJA) then
                call BANDED_IAJA(N,ML,MU)
                NZB = IAB(N+1) - 1
             end if
           end if
        end if

        if ((MITER==2 .or. MITER==5) .and. USE_JACSP) then
!         Allocate the arrays needed by DVJAC/JACSPD.
          if (allocated(INDROWDS)) then
             deallocate (INDROWDS, INDCOLDS, NGRPDS, IPNTRDS, JPNTRDS, &
               IWADS, IWKDS, IOPTDS, YSCALEDS, WKDS, FACDS)
             call CHECK_STAT(IER,440)
          end if
          allocate (INDROWDS(1), INDCOLDS(1), NGRPDS(1), IPNTRDS(1),   &
            JPNTRDS(1), IWADS(1), IWKDS(50+N), IOPTDS(5), YSCALEDS(N), &
            WKDS(3*N), FACDS(N), STAT=IER)
          call CHECK_STAT(IER,450)
!         For use in DVJAC:
          IOPTDS(4) = 0
        end if

! Set work array pointers and check lengths LRW and LIW. Pointers
! to segments of RWORK and IWORK are named by prefixing L to the
! name of the segment. e.g., the segment YH starts at RWORK(LYH).
! Within WM, LOCJS is the location of the saved Jacobian (JSV > 0).

80      LYH = 21
        if (ISTATE==1) NYH = N
        LENRW = LYH + (MAXORD+1)*NYH - 1
        IWORK(17) = LENRW
        if (LENRW>LRW) goto 580
        if (LENRW/=LRW) goto 580
        LWM = 1
!       Save MAXORD in case the calling program calls with ISTATE=3.
        PREVIOUS_MAXORD = MAXORD
        JCO = max(0,JSV)
!       IF (MITER==0) LENWM = 2
        if (MITER==0) LENWM = 0
        if (MITER==1 .or. MITER==2) then
!         LENWM = 2 + (1+JCO)*N*N
!         LOCJS = N*N + 3
          LENWM = (1+JCO)*N*N
          LOCJS = N*N + 1
        end if
!       IF (MITER==3) LENWM = N + 2
        if (MITER==3) LENWM = N
        if (MITER==4 .or. MITER==5) then
          MBAND = ML + MU + 1
          LENP = (MBAND+ML)*N
          LENJ = MBAND*N
!         LENWM = 2 + LENP + JCO*LENJ
!         LOCJS = LENP + 3
          LENWM = LENP + JCO*LENJ
          LOCJS = LENP + 1
        end if
        if (MITER==6 .or. MITER==7) then
!         LENWM = 2
          LENWM = 0
        end if
        if (LENWM>LWMDIM) goto 730
        if (LENWM/=LWMDIM) goto 730
        LIWM = 1
        LENIW = 30 + N
        if (MITER==0 .or. MITER==3) LENIW = 30
        IWORK(18) = LENIW
        if (LENIW>LIW) goto 590
!       Check RTOL and ATOL for legality.
        RTOLI = OPTS%RTOL(1)
        ATOLI = OPTS%ATOL(1)
        do I = 1, N
          if (ITOL>=3) RTOLI = OPTS%RTOL(I)
          if (ITOL==2 .or. ITOL==4) ATOLI = OPTS%ATOL(I)
          if (RTOLI<ZERO) goto 600
          if (ATOLI<ZERO) goto 610
        end do
        if (ISTATE==1) goto 100
!       If ISTATE = 3, set flag to signal parameter changes to DVSTEP.
        JSTART = -1
        if (NQ<=MAXORD) goto 90
!       MAXORD was reduced below NQ. Copy YH(*,MAXORD+2) into SAVF.
!       YH(*,MAXORD+2) was copied to the YHNQP2 array at the
!       beginning of DVODE when the ISTATE=3 call was made.
        call DCOPY_F90(N,YHNQP2,1,SAVF,1)
!       Reload WM1 since LWM may have changed.
90      if (MITER>0) WM1 = sqrt(UROUND)
!       ISTATC controls the determination of the sparsity arrays
!       if the sparse solution option is used.
        ISTATC = ISTATE
        goto 130

! Block C.
! The next block is for the initial call only (ISTATE = 1).
! It contains all remaining initializations, the initial call to F,
! and the calculation of the initial step size.
! The error weights in EWT are inverted after being loaded.

100     UROUND = epsilon(ONE)
        U125 = UROUND ** 0.125_WP
        U325 = UROUND ** 0.325_WP
!       ISTATC controls the determination of the sparsity arrays if
!       the sparse solution option is used.
        ISTATC = ISTATE
        TN = T
        if (ITASK/=4 .and. ITASK/=5) goto 110
        TCRIT = RWORK(1)
        if ((TCRIT-TOUT)*(TOUT-T)<ZERO) goto 660
        if (abs(H0)>ZERO .and. (T+H0-TCRIT)*H0>ZERO) H0 = TCRIT - T
110     JSTART = 0
        if (MITER>0) WM1 = sqrt(UROUND)
        CCMXJ = PT2
        MSBJ = 50
        MSBG = 75
        NHNIL = 0
        NST = 0
        NJE = 0
        NNI = 0
        NCFN = 0
        NETF = 0
        NLU = 0
        NSLJ = 0
        NSLAST = 0
        HU = ZERO
        NQU = 0
        MB28 = 0
!_______________________________________________________________________
! *****MA48 build change point. Insert this statement.
!       MB48 = 0
!_______________________________________________________________________
        NSLG = 0
        NGE = 0
!       Initial call to F. (LF0 points to YH(*,2).)
        LF0 = LYH + NYH
        call F(N,T,Y,RWORK(LF0))
        NFE = 1
!       Load the initial value vector in YH.
        call DCOPY_F90(N,Y,1,RWORK(LYH),1)
!       Load and invert the EWT array. (H is temporarily set to 1.0.)
        NQ = 1
        H = ONE
        call DEWSET(N,ITOL,OPTS%RTOL,OPTS%ATOL,RWORK(LYH),EWT)
        do I = 1, N
          if (EWT(I)<=ZERO) goto 620
          EWT(I) = ONE/EWT(I)
        end do
        NNZ = 0
        NGP = 0
        if (OPTS%SPARSE) then
          ISTATC = ISTATE
        end if
        if (abs(H0)>ZERO) goto 120
!       Call DVHIN to set initial step size H0 to be attempted.
        call DVHIN(N,T,RWORK(LYH),RWORK(LF0),F,TOUT,EWT,ITOL,OPTS%ATOL,Y,ACOR, &
          H0,NITER,IER)
        NFE = NFE + NITER
        if (IER/=0) goto 630
!       Adjust H0 if necessary to meet HMAX bound.
120     RH = abs(H0)*HMXI
        if (RH>ONE) H0 = H0/RH
!       Load H with H0 and scale YH(*,2) by H0.
        H = H0
        call DSCAL_F90(N,H0,RWORK(LF0),1)
!       GOTO 270
!       Check for a zero of g at T.
        IRFND = 0
        TOUTC = TOUT
        if (NGC==0) goto 210
        call DVCHECK(1,GFUN,NEQ,Y,RWORK(LYH),NYH,G0,G1,GX,IRT)
        if (IRT==0) goto 210
        goto 720

! Block D.
! The next code block is for continuation calls only (ISTATE = 2 or 3)
! and is to check stop conditions before taking a step.

130     NSLAST = NST
        IRFP = IRFND
        if (NGC==0) goto 140
        if (ITASK==1 .or. ITASK==4) TOUTC = TOUT
        call DVCHECK(2,GFUN,NEQ,Y,RWORK(LYH),NYH,G0,G1,GX,IRT)
        if (IRT/=1) goto 140
        IRFND = 1
        ISTATE = 3
        T = T0ST
        goto 330
140     continue
        IRFND = 0
        if (IRFP==1 .and. (abs(TLAST-TN)>ZERO) .and. ITASK==2) goto 310
        KUTH = 0
        goto (150,200,160,170,180) ITASK
150     if ((TN-TOUT)*H<ZERO) goto 200
        call DVINDY_CORE(TOUT,0,RWORK(LYH),NYH,Y,IFLAG)
        if (IFLAG/=0) goto 680
        if (BOUNDS) then
          do I = 1, NDX
            Y(IDX(I)) = max(Y(IDX(I)),LB(I))
            Y(IDX(I)) = min(Y(IDX(I)),UB(I))
          end do
        end if
        T = TOUT
        goto 320
160     TP = TN - HU*(ONE+HUN*UROUND)
        if ((TP-TOUT)*H>ZERO) goto 640
        if ((TN-TOUT)*H<ZERO) goto 200
        goto 310
170     TCRIT = RWORK(1)
        if ((TN-TCRIT)*H>ZERO) goto 650
        if ((TCRIT-TOUT)*H<ZERO) goto 660
        if ((TN-TOUT)*H<ZERO) goto 190
        call DVINDY_CORE(TOUT,0,RWORK(LYH),NYH,Y,IFLAG)
        if (IFLAG/=0) goto 680
        if (BOUNDS) then
          do I = 1, NDX
            Y(IDX(I)) = max(Y(IDX(I)),LB(I))
            Y(IDX(I)) = min(Y(IDX(I)),UB(I))
          end do
        end if
        T = TOUT
        goto 320
180     TCRIT = RWORK(1)
        if ((TN-TCRIT)*H>ZERO) goto 650
190     HMX = abs(TN) + abs(H)
        IHIT = abs(TN-TCRIT) <= HUN*UROUND*HMX
        if (IHIT) goto 310
        TNEXT = TN + HNEW*(ONE+FOUR*UROUND)
        if ((TNEXT-TCRIT)*H<=ZERO) goto 200
        H = (TCRIT-TN)*(ONE-FOUR*UROUND)
        KUTH = 1

! Block E.
! The next block is normally executed for all calls and contains
! the call to the one-step core integrator DVSTEP. This is a
! looping point for the integration steps.
! First check for too many steps being taken, update EWT(if not
! at start of problem), check for too much accuracy being
! requested, and check for H below the roundoff level in T.

200     continue
        if ((NST-NSLAST)>=MXSTEP) goto 340
        call DEWSET(N,ITOL,OPTS%RTOL,OPTS%ATOL,RWORK(LYH),EWT)
        do I = 1, N
          if (EWT(I)<=ZERO) goto 350
          EWT(I) = ONE/EWT(I)
        end do
210     TOLSF = UROUND*DVNORM(N,RWORK(LYH),EWT)
        IPCUTH = -1
        if (TOLSF<=ONE) goto 220
        TOLSF = TOLSF*TWO
        if (NST==0) goto 670
        goto 360
220     continue
        IPCUTH = IPCUTH + 1
        if (IPCUTH>=IPCUTH_MAX) then
          MSG = 'Too many step reductions to prevent'
          call XERRDV(MSG,940,1,0,0,0,0,ZERO,ZERO)
          MSG = 'an infeasible prediction.'
          call XERRDV(MSG,940,1,0,0,0,0,ZERO,ZERO)
!         Retract the solution to TN:
          call DVNRDN(RWORK(LYH),NYH,N,NQ)
          ACOR(1:N) = ZERO
          ISTATE = -7
          goto 410
        end if
        if (abs((TN+H)-TN)>ZERO) goto 230
        NHNIL = NHNIL + 1
        if (NHNIL>MXHNIL) goto 230
        MSG = 'Warning: internal T(=R1) and H(=R2) are such that'
        call XERRDV(MSG,950,1,0,0,0,0,ZERO,ZERO)
        MSG = 'in the machine, T + H = T on the next step.'
        call XERRDV(MSG,950,1,0,0,0,0,ZERO,ZERO)
        MSG = '(H = step size). The solver will continue anyway.'
        call XERRDV(MSG,950,1,0,0,0,2,TN,H)
        if (NHNIL<MXHNIL) goto 230
        MSG = 'The above warning has been issued I1 times.'
        call XERRDV(MSG,950,1,0,0,0,0,ZERO,ZERO)
        MSG = 'It will not be issued again for this problem.'
        call XERRDV(MSG,950,1,1,MXHNIL,0,0,ZERO,ZERO)
230     continue
        if (BOUNDS) then
!         Check positive components for infeasible prediction; reduce 
!         step size if infeasible prediction will occur in DVNLSD.
          ACOR(1:N) = RWORK(LYH:LYH+N-1)
!         Predict:
          call DVNRDP(RWORK(LYH:LRW),NYH,N,NQ)
          do I = 1, NDX
            if ((ACOR(IDX(I))>LB(I) .and. RWORK(LYH+ &
                IDX(I)-1)<LB(I)) .or. (ACOR(IDX(I))<UB(I) .and. RWORK(LYH+ &
                IDX(I)-1)>UB(I))) then
!              Retract:
              call DVNRDN(RWORK(LYH),NYH,N,NQ)
              H = HALF*H
              ETA = HALF
!              Rescale:
              call DVNRDS(RWORK(LYH),NYH,N,L,ETA)
              goto 220
            end if
          end do
!         Retract:
          call DVNRDN(RWORK(LYH),NYH,N,NQ)
          ACOR(1:N) = ZERO
        end if

!       CALL DVSTEP(Y,YH,LDYH,YH1,EWT,SAVF,ACOR,WM,IWM,F,JAC, &
!         VNLS,OPTS%ATOL,ITOL)
!_______________________________________________________________________

        if (MITER/=6 .and. MITER/=7) then
          call DVSTEP(Y,RWORK(LYH),NYH,RWORK(LYH),EWT,SAVF,ACOR,WM, &
            IWORK(LIWM),F,JAC,DVNLSD,OPTS%ATOL,ITOL)
        else
          call DVSTEP(Y,RWORK(LYH),NYH,RWORK(LYH),EWT,SAVF,ACOR,WM, &
            IWORK(LIWM),F,JAC,DVNLSS28,OPTS%ATOL,ITOL)
        end if
! *****MA48 build change point. Replace above with these statements.
!     IF (MITER /= 6 .AND. MITER /= 7) THEN
!        CALL DVSTEP(Y,RWORK(LYH),NYH,RWORK(LYH),EWT,SAVF,ACOR,WM, &
!          IWORK(LIWM),F,JAC,DVNLSD,OPTS%ATOL,ITOL)
!     ELSE
!        IF (USE_MA48_FOR_SPARSE) THEN
!           CALL DVSTEP(Y,RWORK(LYH),NYH,RWORK(LYH),EWT,SAVF,ACOR,WM, &
!             IWORK(LIWM), F, JAC, DVNLSS48,OPTS%ATOL,ITOL)
!        ELSE
!           CALL DVSTEP(Y,RWORK(LYH),NYH,RWORK(LYH),EWT,SAVF,ACOR,WM, &
!           IWORK(LIWM),F,JAC,DVNLSS28,OPTS%ATOL,ITOL)
!        END IF
!     END IF
!_______________________________________________________________________

        KGO = 1 - KFLAG
!       Branch on KFLAG. Note: In this version, KFLAG can not be set to
!                            -3; KFLAG = 0, -1, -2.
        goto (240,370,380) KGO

! Block F.
! The following block handles the case of a successful return from the
! core integrator (KFLAG = 0). Test for stop conditions.

240     INIT = 1
        KUTH = 0
        goto (250,310,270,280,300) ITASK
!       ITASK = 1. If TOUT has been reached, interpolate.
250     continue
        if (NGC==0) goto 260
        call DVCHECK(3,GFUN,NEQ,Y,RWORK(LYH),NYH,G0,G1,GX,IRT)
        if (IRT/=1) goto 260
        IRFND = 1
        ISTATE = 3
        T = T0ST
        goto 330
260     continue
        if ((TN-TOUT)*H<ZERO) goto 200
        call DVINDY_CORE(TOUT,0,RWORK(LYH),NYH,Y,IFLAG)
        if (BOUNDS) then
          do I = 1, NDX
            Y(IDX(I)) = max(Y(IDX(I)),LB(I))
            Y(IDX(I)) = min(Y(IDX(I)),UB(I))
          end do
        end if
        T = TOUT
        goto 320
!       ITASK = 3. Jump to exit if TOUT was reached.
270     if ((TN-TOUT)*H>=ZERO) goto 310
        goto 200
!       ITASK = 4. See if TOUT or TCRIT was reached. Adjust H if necessary.
280     if ((TN-TOUT)*H<ZERO) goto 290
        call DVINDY_CORE(TOUT,0,RWORK(LYH),NYH,Y,IFLAG)
        if (BOUNDS) then
          do I = 1, NDX
            Y(IDX(I)) = max(Y(IDX(I)),LB(I))
            Y(IDX(I)) = min(Y(IDX(I)),UB(I))
          end do
        end if
        T = TOUT
        goto 320
290     HMX = abs(TN) + abs(H)
        IHIT = abs(TN-TCRIT) <= HUN*UROUND*HMX
        if (IHIT) goto 310
        TNEXT = TN + HNEW*(ONE+FOUR*UROUND)
        if ((TNEXT-TCRIT)*H<=ZERO) goto 200
        H = (TCRIT-TN)*(ONE-FOUR*UROUND)
        KUTH = 1
        goto 200
!       ITASK = 5. See if TCRIT was reached and jump to exit.
300     HMX = abs(TN) + abs(H)
        IHIT = abs(TN-TCRIT) <= HUN*UROUND*HMX

! Block G.
! The following block handles all successful returns from DVODE.
! If ITASK /= 1, Y is loaded from YH and T is set accordingly.
! ISTATE is set to 2, and the optional output is loaded into the
! work arrays before returning.

310     continue
        call DCOPY_F90(N,RWORK(LYH),1,Y,1)
        T = TN
        if (ITASK/=4 .and. ITASK/=5) goto 320
        if (IHIT) T = TCRIT
320     ISTATE = 2
330     continue
        RWORK(11) = HU
        RWORK(12) = HNEW
        RWORK(13) = TN
        IWORK(10) = NGE
        IWORK(11) = NST
        IWORK(12) = NFE
        IWORK(13) = NJE
        IWORK(14) = NQU
        IWORK(15) = NEWQ
        IWORK(19) = NLU
        IWORK(20) = NNI
        IWORK(21) = NCFN
        IWORK(22) = NETF
        TLAST = T
        RUSER(1:LRWUSER) = RWORK(1:LRWUSER)
        IUSER(1:LIWUSER) = IWORK(1:LIWUSER)
!       Warn the user if |y(t)| < ATOL:
        if (ISTATE==2 .or. ISTATE==3) then
          if (YMAXWARN) then
            ATOLI = OPTS%ATOL(1)
            do I = 1, N
              if (ITOL==2 .or. ITOL==4) ATOLI = OPTS%ATOL(I)
              if (abs(Y(I))<ATOLI) then
                MSG = 'Warning: Component I1 of the solution is'
                call XERRDV(MSG,960,1,0,0,0,0,ZERO,ZERO)
                MSG = 'smaller in magnitude than component I1'
                call XERRDV(MSG,960,1,0,0,0,0,ZERO,ZERO)
                MSG = 'of the absolute error tolerance vector.'
                call XERRDV(MSG,960,1,1,I,0,0,ZERO,ZERO)
              end if
            end do
          end if
        end if
        return

! Block H.
! The following block handles all unsuccessful returns other than
! those for illegal input. First the error message routine is called.
! if there was an error test or convergence test failure, IMXER is set.
! Then Y is loaded from YH, and T is set to TN. The optional output
! is loaded into the work arrays before returning.

!       The maximum number of steps was taken before reaching TOUT.
340     MSG = 'At current T(=R1), MXSTEP(=I1) steps'
        call XERRDV(MSG,970,1,0,0,0,0,ZERO,ZERO)
        MSG = 'taken on this call before reaching TOUT.'
        call XERRDV(MSG,970,1,1,MXSTEP,0,1,TN,ZERO)
        ISTATE = -1
        goto 410
!       EWT(i) <= 0.0 for some i (not at start of problem).
350     EWTI = EWT(I)
        MSG = 'At T(=R1), EWT(I1) has become R2 <= 0.'
        call XERRDV(MSG,980,1,1,I,0,2,TN,EWTI)
        ISTATE = -6
        goto 410
!       Too much accuracy requested for machine precision.
360     MSG = 'At T(=R1), too much accuracy was requested'
        call XERRDV(MSG,990,1,0,0,0,0,ZERO,ZERO)
        MSG = 'for precision of machine:   see TOLSF(=R2).'
        call XERRDV(MSG,990,1,0,0,0,2,TN,TOLSF)
        RWORK(14) = TOLSF
        ISTATE = -2
        goto 410
!       KFLAG = -1. Error test failed repeatedly or with ABS(H) = HMIN.
370     MSG = 'At T(=R1) and step size H(=R2), the error'
        call XERRDV(MSG,1000,1,0,0,0,0,ZERO,ZERO)
        MSG = 'test failed repeatedly or with ABS(H) = HMIN.'
        call XERRDV(MSG,1000,1,0,0,0,2,TN,H)
        ISTATE = -4
        goto 390
!       KFLAG = -2. Convergence failed repeatedly or with ABS(H) = HMIN.
380     MSG = 'At T(=R1) and step size H(=R2), the'
        call XERRDV(MSG,1010,1,0,0,0,0,ZERO,ZERO)
        MSG = 'corrector convergence failed repeatedly'
        call XERRDV(MSG,1010,1,0,0,0,0,ZERO,ZERO)
        MSG = 'or with ABS(H) = HMIN.'
        call XERRDV(MSG,1010,1,0,0,0,2,TN,H)
        ISTATE = -5
!       Compute IMXER if relevant.
390     BIG = ZERO
        IMXER = 1
        do 400 I = 1, N
          SIZEST = abs(ACOR(I)*EWT(I))
          if (BIG>=SIZEST) goto 400
          BIG = SIZEST
          IMXER = I
400     end do
        IWORK(16) = IMXER
!       Set Y vector, T, and optional output.
410     continue
        call DCOPY_F90(N,RWORK(LYH),1,Y,1)
        T = TN
        RWORK(11) = HU
        RWORK(12) = H
        RWORK(13) = TN
        IWORK(10) = NGE
        IWORK(11) = NST
        IWORK(12) = NFE
        IWORK(13) = NJE
        IWORK(14) = NQU
        IWORK(15) = NQ
        IWORK(19) = NLU
        IWORK(20) = NNI
        IWORK(21) = NCFN
        IWORK(22) = NETF
        TLAST = T
        RUSER(1:LRWUSER) = RWORK(1:LRWUSER)
        IUSER(1:LIWUSER) = IWORK(1:LIWUSER)
        return

! Block I.
! The following block handles all error returns due to illegal input
! (ISTATE = -3), as detected before calling the core integrator.
! First the error message routine is called. If the illegal input
! is a negative ISTATE, the run is aborted (apparent infinite loop).

420     MSG = 'ISTATE(=I1) is illegal.'
        call XERRDV(MSG,1020,1,1,ISTATE,0,0,ZERO,ZERO)
        if (ISTATE<0) goto 750
        goto 740
430     MSG = 'ITASK(=I1) is illegal.'
        call XERRDV(MSG,1030,1,1,ITASK,0,0,ZERO,ZERO)
        goto 740
440     MSG = 'ISTATE(=I1) > 1 but DVODE is not initialized.'
        call XERRDV(MSG,1040,1,1,ISTATE,0,0,ZERO,ZERO)
        goto 740
450     MSG = 'NEQ (=I1) < 1.'
        call XERRDV(MSG,1050,1,1,NEQ,0,0,ZERO,ZERO)
        goto 740
460     MSG = 'ISTATE = 3 and NEQ increased (I1 to I2).'
        call XERRDV(MSG,1060,1,2,N,NEQ,0,ZERO,ZERO)
        goto 740
465     MSG = 'This version of DVODE requires does not allow NEQ to be reduced.'
        call XERRDV(MSG,1070,2,0,0,0,0,ZERO,ZERO)
        goto 740
470     MSG = 'ITOL(=I1) is illegal.'
        call XERRDV(MSG,1080,1,1,ITOL,0,0,ZERO,ZERO)
        goto 740
480     MSG = 'IOPT(=I1) is illegal.'
        call XERRDV(MSG,1090,1,1,IOPT,0,0,ZERO,ZERO)
        goto 740
490     MSG = 'MF(=I1) is illegal.'
        call XERRDV(MSG,1100,1,1,MF,0,0,ZERO,ZERO)
        goto 740
500     MSG = 'ML(=I1) illegal: < 0 or >= NEQ (=I2)'
        call XERRDV(MSG,1110,1,2,ML,NEQ,0,ZERO,ZERO)
        goto 740
510     MSG = 'MU(=I1) illegal: < 0 or >= NEQ (=I2)'
        call XERRDV(MSG,1120,1,2,MU,NEQ,0,ZERO,ZERO)
        goto 740
520     MSG = 'MAXORD(=I1) < 0.'
        call XERRDV(MSG,1130,1,1,MAXORD,0,0,ZERO,ZERO)
        goto 740
530     MSG = 'MXSTEP(=I1) < 0.'
        call XERRDV(MSG,1140,1,1,MXSTEP,0,0,ZERO,ZERO)
        goto 740
540     MSG = 'MXHNIL(=I1) < 0.'
        call XERRDV(MSG,1150,1,1,MXHNIL,0,0,ZERO,ZERO)
        goto 740
550     MSG = 'TOUT(=R1) is behind T(=R2).'
        call XERRDV(MSG,1160,1,0,0,0,2,TOUT,T)
        MSG = 'The integration direction is given by H0 (=R1).'
        call XERRDV(MSG,1160,1,0,0,0,1,H0,ZERO)
        goto 740
560     MSG = 'HMAX(=R1) < 0.'
        call XERRDV(MSG,1170,1,0,0,0,1,HMAX,ZERO)
        goto 740
570     MSG = 'HMIN(=R1) < 0.'
        call XERRDV(MSG,1180,1,0,0,0,1,HMIN,ZERO)
        goto 740
580     continue
        MSG = 'RWORK length needed, LENRW(=I1) > LRW(=I2)'
        call XERRDV(MSG,1190,1,2,LENRW,LRW,0,ZERO,ZERO)
        goto 740
590     continue
        MSG = 'IWORK length needed, LENIW(=I1) > LIW(=I2)'
        call XERRDV(MSG,1200,1,2,LENIW,LIW,0,ZERO,ZERO)
        goto 740
600     MSG = 'RTOL(I1) is R1 < 0.'
        call XERRDV(MSG,1210,1,1,I,0,1,RTOLI,ZERO)
        goto 740
610     MSG = 'ATOL(I1) is R1 < 0.'
        call XERRDV(MSG,1220,1,1,I,0,1,ATOLI,ZERO)
        goto 740
620     EWTI = EWT(I)
        MSG = 'EWT(I1) is R1 <= 0.'
        call XERRDV(MSG,1230,1,1,I,0,1,EWTI,ZERO)
        goto 740
630     continue
        MSG = 'TOUT(=R1) too close to T(=R2) to start.'
        call XERRDV(MSG,1240,1,0,0,0,2,TOUT,T)
        goto 740
640     continue
        MSG = 'ITASK = I1 and TOUT(=R1) < TCUR - HU(=R2).'
        call XERRDV(MSG,1250,1,1,ITASK,0,2,TOUT,TP)
        goto 740
650     continue
        MSG = 'ITASK = 4 or 5 and TCRIT(=R1) < TCUR(=R2).'
        call XERRDV(MSG,1260,1,0,0,0,2,TCRIT,TN)
        goto 740
660     continue
        MSG = 'ITASK = 4 or 5 and TCRIT(=R1) < TOUT(=R2).'
        call XERRDV(MSG,1270,1,0,0,0,2,TCRIT,TOUT)
        goto 740
670     MSG = 'At the start of the problem, too much'
        call XERRDV(MSG,1280,1,0,0,0,0,ZERO,ZERO)
        MSG = 'accuracy was requested for precision'
        call XERRDV(MSG,1280,1,0,0,0,1,TOLSF,ZERO)
        MSG = 'of machine: see TOLSF(=R1).'
        call XERRDV(MSG,1280,1,0,0,0,1,TOLSF,ZERO)
        RWORK(14) = TOLSF
        goto 740
680     MSG = 'Trouble from DVINDY. ITASK = I1, TOUT = R1.'
        call XERRDV(MSG,1290,1,1,ITASK,0,1,TOUT,ZERO)
        goto 740
690     MSG = 'SETH must be nonnegative.'
        call XERRDV(MSG,1300,1,0,0,0,0,ZERO,ZERO)
        goto 740
700     MSG = 'NG(=I1) < 0.'
        call XERRDV(MSG,1310,0,1,NG,0,0,ZERO,ZERO)
        goto 740
710     MSG = 'NG changed (from I1 to I2) illegally, i.e.,'
        call XERRDV(MSG,1320,1,0,0,0,0,ZERO,ZERO)
        MSG = 'not immediately after a root was found.'
        call XERRDV(MSG,1320,1,2,NGC,NG,0,ZERO,ZERO)
        goto 740
720     MSG = 'One or more components of g has a root'
        call XERRDV(MSG,1330,1,0,0,0,0,ZERO,ZERO)
        MSG = 'too near to the initial point.'
        call XERRDV(MSG,1330,1,0,0,0,0,ZERO,ZERO)
        goto 740
730     continue
        MSG = 'WM length needed, LENWM(=I1) > LWMDIM(=I2)'
        call XERRDV(MSG,1340,1,2,LENWM,LWMDIM,0,ZERO,ZERO)

740     continue
        ISTATE = -3
        RUSER(1:LRWUSER) = RWORK(1:LRWUSER)
        IUSER(1:LIWUSER) = IWORK(1:LIWUSER)
        return

750     MSG = 'Run aborted:  apparent infinite loop.'
        call XERRDV(MSG,1350,2,0,0,0,0,ZERO,ZERO)
        RUSER(1:LRWUSER) = RWORK(1:LRWUSER)
        IUSER(1:LIWUSER) = IWORK(1:LIWUSER)
        return

      end subroutine DVODE
!_______________________________________________________________________

      subroutine DVHIN(N,T0,Y0,YDOT,F,TOUT,EWT,ITOL,ATOL,Y,TEMP,H0, &
        NITER,IER)
! ..
! Calculate the initial step size.
! ..
! This routine computes the step size, H0, to be attempted on the
! first step, when the user has not supplied a value for this.
! First we check that TOUT - T0 differs significantly from zero. Then
! an iteration is done to approximate the initial second derivative
! and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1.
! A bias factor of 1/2 is applied to the resulting h.
! The sign of H0 is inferred from the initial values of TOUT and T0.
! Communication with DVHIN is done with the following variables:
! N      = Size of ODE system, input.
! T0     = Initial value of independent variable, input.
! Y0     = Vector of initial conditions, input.
! YDOT   = Vector of initial first derivatives, input.
! F      = Name of subroutine for right-hand side f(t,y), input.
! TOUT   = First output value of independent variable
! UROUND = Machine unit roundoff
! EWT, ITOL, ATOL = Error weights and tolerance parameters
!                   as described in the driver routine, input.
! Y, TEMP = Work arrays of length N.
! H0     = Step size to be attempted, output.
! NITER  = Number of iterations (and of f evaluations) to compute H0,
!          output.
! IER    = The error flag, returned with the value
!          IER = 0  if no trouble occurred, or
!          IER = -1 if TOUT and T0 are considered too close to proceed.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP), intent (INOUT) :: H0
        real (WP), intent (IN) :: T0, TOUT
        integer :: IER
        integer, intent (IN) :: ITOL, N
        integer, intent (INOUT) :: NITER
! ..
! .. Array Arguments ..
        real (WP), intent (IN) :: ATOL(*), EWT(*), Y0(*)
        real (WP), intent (INOUT) :: TEMP(*), Y(*), YDOT(*)
! ..
! .. Subroutine Arguments ..
        external F
! ..
! .. Local Scalars ..
        real (WP) :: AFI, ATOLI, DELYI, H, HG, HLB, HNEW, HRAT, HUB, T1, &
          TDIST, TROUND, YDDNRM
        integer :: I, ITER
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, MAX, SIGN, SQRT
! ..
! .. FIRST EXECUTABLE STATEMENT DVHIN
! ..
        NITER = 0
        TDIST = abs(TOUT-T0)
        TROUND = UROUND*max(abs(T0),abs(TOUT))
        if (TDIST<TWO*TROUND) goto 40

!       Set a lower bound on H based on the roundoff level in T0
!       and TOUT.
        HLB = HUN*TROUND
!       Set an upper bound on H based on TOUT-T0 and the initial
!       Y and YDOT.
        HUB = PT1*TDIST
        ATOLI = ATOL(1)
        do I = 1, N
          if (ITOL==2 .or. ITOL==4) ATOLI = ATOL(I)
          DELYI = PT1*abs(Y0(I)) + ATOLI
          AFI = abs(YDOT(I))
          if (AFI*HUB>DELYI) HUB = DELYI/AFI
        end do
!       Set initial guess for H as geometric mean of upper and
!       lower bounds.
        ITER = 0
        HG = sqrt(HLB*HUB)
!       If the bounds have crossed, exit with the mean value.
        if (HUB<HLB) then
          H0 = HG
          goto 30
        end if

!       Looping point for iteration.
10      continue
!       Estimate the second derivative as a difference quotient in f.
        H = sign(HG,TOUT-T0)
        T1 = T0 + H
        Y(1:N) = Y0(1:N) + H*YDOT(1:N)
        call F(N,T1,Y,TEMP)
        NFE = NFE + 1
        TEMP(1:N) = (TEMP(1:N)-YDOT(1:N))/H
        YDDNRM = DVNORM(N,TEMP,EWT)
!       Get the corresponding new value of h.
        if (YDDNRM*HUB*HUB>TWO) then
          HNEW = sqrt(TWO/YDDNRM)
        else
          HNEW = sqrt(HG*HUB)
        end if
        ITER = ITER + 1

! Test the stopping conditions.
! Stop if the new and previous h values differ by a factor of < 2.
! Stop if four iterations have been done. Also, stop with previous h
! if HNEW/HG > 2 after first iteration, as this probably means that
! the second derivative value is bad because of cancellation error.

        if (ITER>=4) goto 20
        HRAT = HNEW/HG
        if ((HRAT>HALF) .and. (HRAT<TWO)) goto 20
        if ((ITER>=2) .and. (HNEW>TWO*HG)) then
          HNEW = HG
          goto 20
        end if
        HG = HNEW
        goto 10

!       Iteration done. Apply bounds, bias factor, and sign. Then exit.
20      H0 = HNEW*HALF
        if (H0<HLB) H0 = HLB
        if (H0>HUB) H0 = HUB
30      H0 = sign(H0,TOUT-T0)
        NITER = ITER
        IER = 0
        return
!       Error return for TOUT - T0 too small.
40      IER = -1
        return

      end subroutine DVHIN
!_______________________________________________________________________

      subroutine DVINDY_CORE(T,K,YH,LDYH,DKY,IFLAG)
! ..
! Interpolate the solution and derivative.
! ..
! DVINDY_CORE computes interpolated values of the K-th derivative
! of the dependent variable vector y, and stores it in DKY. This
! routine is called within the package with K = 0 and T = TOUT,
! but may also be called by the user for any K up to the current
! order. (See detailed instructions in the usage documentation.)
! The computed values in DKY are gotten by interpolation using the
! Nordsieck history array YH. This array corresponds uniquely to a
! vector-valued polynomial of degree NQCUR or less, and DKY is set
! to the K-th derivative of this polynomial at T.
! The formula for DKY is:
!              q
!  DKY(i)  =  sum  c(j,K) * (T - TN)**(j-K) * H**(-j) * YH(i,j+1)
!             j=K
! where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, TN = TCUR, H = HCUR.
! The quantities  NQ = NQCUR, L = NQ+1, N, TN, and H are
! communicated by PRIVATE variables. The above sum is done in reverse
! order.
! IFLAG is returned negative if either K or T is out of bounds.
! Discussion above and comments in driver explain all variables.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP), intent (IN) :: T
        integer, intent (INOUT) :: IFLAG
        integer, intent (IN) :: K, LDYH
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: DKY(*), YH(LDYH,*)
! ..
! .. Local Scalars ..
        real (WP) :: C, R, S, TFUZZ, TN1, TP
        integer :: IC, J, JB, JB2, JJ, JJ1, JP1
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, real, SIGN
! ..
! .. FIRST EXECUTABLE STATEMENT DVINDY_CORE
! ..
        IFLAG = 0
        if (K<0 .or. K>NQ) goto 40
!       TFUZZ = HUN * UROUND * (TN + HU)
        TFUZZ = HUN*UROUND*sign(abs(TN)+abs(HU),HU)
        TP = TN - HU - TFUZZ
        TN1 = TN + TFUZZ
        if ((T-TP)*(T-TN1)>ZERO) goto 50
        S = (T-TN)/H
        IC = 1
        if (K==0) goto 10
        JJ1 = L - K
        do JJ = JJ1, NQ
          IC = IC*JJ
        end do
10      C = real(IC)
        DKY(1:N) = C*YH(1:N,L)
        if (K==NQ) goto 30
        JB2 = NQ - K
        do JB = 1, JB2
          J = NQ - JB
          JP1 = J + 1
          IC = 1
          if (K==0) goto 20
          JJ1 = JP1 - K
          do JJ = JJ1, J
            IC = IC*JJ
          end do
20        C = real(IC)
          DKY(1:N) = C*YH(1:N,JP1) + S*DKY(1:N)
        end do
30      R = H**(-K)
        call DSCAL_F90(N,R,DKY,1)
        return
40      MSG = 'Error in DVINDY, K(=I1) is illegal.'
        call XERRDV(MSG,1360,1,1,K,0,0,ZERO,ZERO)
        IFLAG = -1
        return
50      MSG = 'Error in DVINDY, T(=R1) is illegal. T is not'
        call XERRDV(MSG,1370,1,0,0,0,1,T,ZERO)
        MSG = 'in interval TCUR - HU(= R1) to TCUR(=R2)'
        call XERRDV(MSG,1370,1,0,0,0,2,TP,TN)
        IFLAG = -2
        return

      end subroutine DVINDY_CORE
!_______________________________________________________________________

      subroutine DVINDY_BNDS(T,K,YH,LDYH,DKY,IFLAG)
! ..
! Interpolate the solution and derivative and enforce nonnegativity
! (used only if user calls DVINDY and the BOUNDS option is in
! use).
! ..
! This version of DVINDY_CORE enforces nonnegativity and is called
! only by DVINDY (which is called only by the user). It uses the
! private YNNEG array produced by a call to DVINDY_CORE from DVINDY
! to enforce nonnegativity.
! DVINDY_BNDS computes interpolated values of the K-th derivative
! of the dependent variable vector y, and stores it in DKY. This
! routine is called within the package with K = 0 and T = TOUT,
! but may also be called by the user for any K up to the current
! order. (See detailed instructions in the usage documentation.)
! The computed values in DKY are gotten by interpolation using the
! Nordsieck history array YH. This array corresponds uniquely to a
! vector-valued polynomial of degree NQCUR or less, and DKY is set
! to the K-th derivative of this polynomial at T.
! The formula for DKY is:
!              q
!  DKY(i)  =  sum  c(j,K) * (T - TN)**(j-K) * H**(-j) * YH(i,j+1)
!             j=K
! where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, TN = TCUR, H = HCUR.
! The quantities  NQ = NQCUR, L = NQ+1, N, TN, and H are
! communicated by PRIVATE variables. The above sum is done in reverse
! order.
! IFLAG is returned negative if either K or T is out of bounds.
! Discussion above and comments in driver explain all variables.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP), intent (IN) :: T
        integer, intent (INOUT) :: IFLAG
        integer, intent (IN) :: K, LDYH
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: DKY(*), YH(LDYH,*)
! ..
! .. Local Scalars ..
        real (WP) :: C, R, S, TFUZZ, TN1, TP
        integer :: I, IC, J, JB, JB2, JJ, JJ1, JP1
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, real, SIGN
! ..
! .. FIRST EXECUTABLE STATEMENT DVINDY_BNDS
! ..
        IFLAG = 0
        if (K==0) goto 50
        if (K<0 .or. K>NQ) goto 40
!       TFUZZ = HUN * UROUND * (TN + HU)
        TFUZZ = HUN*UROUND*sign(abs(TN)+abs(HU),HU)
        TP = TN - HU - TFUZZ
        TN1 = TN + TFUZZ
        if ((T-TP)*(T-TN1)>ZERO) goto 60
        S = (T-TN)/H
        IC = 1
        if (K==0) goto 10
        JJ1 = L - K
        do JJ = JJ1, NQ
          IC = IC*JJ
        end do
10      C = real(IC)
        DKY(1:N) = C*YH(1:N,L)
        if (BOUNDS) then
          do I = 1, NDX
            if (YNNEG(IDX(I))<LB(I) .or. YNNEG(IDX(I))>UB(I)) DKY(IDX(I)) &
              = ZERO
          end do
        end if
        if (K==NQ) goto 30
        JB2 = NQ - K
        do JB = 1, JB2
          J = NQ - JB
          JP1 = J + 1
          IC = 1
          if (K==0) goto 20
          JJ1 = JP1 - K
          do JJ = JJ1, J
            IC = IC*JJ
          end do
20        C = real(IC)
          DKY(1:N) = C*YH(1:N,JP1) + S*DKY(1:N)
          if (BOUNDS) then
            do I = 1, NDX
              if (YNNEG(IDX(I))<LB(I) .or. YNNEG(IDX(I))>UB(I)) DKY(IDX(I)) &
                = ZERO
            end do
          end if
        end do
30      R = H**(-K)
        call DSCAL_F90(N,R,DKY,1)
        return
40      MSG = 'Error in DVINDY, K(=I1) is illegal.'
        call XERRDV(MSG,1380,1,1,K,0,0,ZERO,ZERO)
        IFLAG = -1
        return
50      MSG = 'DVINDY_BNDS cannot be called with k = 0.'
        call XERRDV(MSG,1390,1,0,0,0,0,ZERO,ZERO)
        IFLAG = -1
        return
60      MSG = 'Error in DVINDY, T(=R1) is illegal. T is not'
        call XERRDV(MSG,1400,1,0,0,0,1,T,ZERO)
        MSG = 'not in interval TCUR - HU(= R1) to TCUR(=R2)'
        call XERRDV(MSG,1400,1,0,0,0,2,TP,TN)
        IFLAG = -2
        return

      end subroutine DVINDY_BNDS
!_______________________________________________________________________

      subroutine DVINDY(T,K,DKY,IFLAG)
! ..
! This is a dummy interface to allow the user to interpolate the
! solution and the derivative (not called by DVODE_F90).
! ..
! May be used if the user wishes to interpolate solution or
! derivative following a successful return from DVODE_F90
! DVINDY computes interpolated values of the K-th derivative of
! the dependent variable vector y, and stores it in DKY. This
! routine is called within the package with K = 0 and T = TOUT,
! but may also be called by the user for any K up to the current
! order. (See detailed instructions in the usage documentation.)
! The computed values in DKY are gotten by interpolation using the
! Nordsieck history array YH. This array corresponds uniquely to a
! vector-valued polynomial of degree NQCUR or less, and DKY is set
! to the K-th derivative of this polynomial at T.
! IFLAG is returned negative if either K or T is out of bounds.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP), intent (IN) :: T
        integer, intent (INOUT) :: IFLAG
        integer, intent (IN) :: K
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: DKY(*)
! ..
! .. Local Scalars ..
        integer :: I, IER
! ..
! .. Intrinsic Functions ..
        intrinsic ALLOCATED, SIZE
! ..
! .. FIRST EXECUTABLE STATEMENT DVINDY
! ..
        call DVINDY_CORE(T,K,RWORK(LYH),N,DKY,IFLAG)
        if (.not.BOUNDS) return

        if (K==0) then
!         Interpolate only the solution.
          call DVINDY_CORE(T,K,RWORK(LYH),N,DKY,IFLAG)
!         Enforce bounds.
          do I = 1, NDX
            if (DKY(IDX(I))<LB(I)) DKY(IDX(I)) = LB(I)
            if (DKY(IDX(I))>UB(I)) DKY(IDX(I)) = UB(I)
          end do
          return
        end if

!       k > 0 - derivatives requested.

!       Make sure space is available for the interpolated solution.
        if (allocated(YNNEG)) then
          if (size(YNNEG)<N) then
            deallocate (YNNEG,STAT=IER)
            call CHECK_STAT(IER,460)
            allocate (YNNEG(N),STAT=IER)
            call CHECK_STAT(IER,470)
          end if
        else
          allocate (YNNEG(N),STAT=IER)
          call CHECK_STAT(IER,480)
        end if
!       Interpolate the solution; do not enforce bounds.
        call DVINDY_CORE(T,0,RWORK(LYH),N,YNNEG,IFLAG)
!       Now interpolate the derivative using YNNEG to enforce
!       nonnegativity.
        call DVINDY_BNDS(T,K,RWORK(LYH),N,DKY,IFLAG)
        return

      end subroutine DVINDY
!_______________________________________________________________________

      subroutine DVSTEP(Y,YH,LDYH,YH1,EWT,SAVF,ACOR,WM,IWM,F,JAC, &
        VNLS,ATOL,ITOL)
! ..
! This is the core step integrator for nonstiff and for dense, banded,
! and sparse solutions.
! ..
! DVSTEP performs one step of the integration of an initial value
! problem for a system of ordinary differential equations. It
! calls subroutine DVNLSD for the solution of the nonlinear system
! arising in the time step. Thus it is independent of the problem
! Jacobian structure and the type of nonlinear system solution method.
! DVSTEP returns a completion flag KFLAG (in PRIVATE variables block).
! A return with KFLAG = -1 or -2 means either ABS(H) = HMIN or 10
! consecutive failures occurred. On a return with KFLAG negative,
! the values of TN and the YH array are as of the beginning of the last
! step, and H is the last step size attempted.
! Communication with DVSTEP is done with the following variables:
! Y      = An array of length N used for the dependent variable vector.
! YH     = An LDYH by LMAX array containing the dependent variables
!          and their approximate scaled derivatives, where
!          LMAX = MAXORD + 1. YH(i,j+1) contains the approximate
!          j-th derivative of y(i), scaled by H**j/factorial(j)
!          (j = 0,1,...,NQ). On entry for the first step, the first
!          two columns of YH must be set from the initial values.
! LDYH   = A constant integer >= N, the first dimension of YH.
!          N is the number of ODEs in the system.
! YH1    = A one-dimensional array occupying the same space as YH.
! EWT    = An array of length N containing multiplicative weights
!          for local error measurements. Local errors in y(i) are
!          compared to 1.0/EWT(i) in various error tests.
! SAVF   = An array of working storage, of length N.
!          also used for input of YH(*,MAXORD+2) when JSTART = -1
!          and MAXORD < the current order NQ.
! ACOR   = A work array of length N, used for the accumulated
!          corrections. On a successful return, ACOR(i) contains
!          the estimated one-step local error in y(i).
! WM,IWM = Real and integer work arrays associated with matrix
!          operations in DVNLSD.
! F      = Dummy name for the user supplied subroutine for f.
! JAC    = Dummy name for the user supplied Jacobian subroutine.
! DVNLSD = Dummy name for the nonlinear system solving subroutine,
!          whose real name is dependent on the method used.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: LDYH, ITOL
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: ACOR(*), EWT(*), SAVF(*), &
          WM(*), Y(*), YH(LDYH,*), YH1(*)
        real (WP), intent (IN) :: ATOL(*)
        integer, intent (INOUT) :: IWM(*)
! ..
! .. Subroutine Arguments ..
        external F, JAC, VNLS
! ..
! .. Local Scalars ..
        real (WP) :: CNQUOT, DDN, DSM, DUP, ETAQ, ETAQM1, ETAQP1, FLOTL, R, &
          TOLD
        integer :: I, I1, I2, IBACK, J, JB, NCF, NFLAG
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, MAX, MIN, real
! ..
! .. FIRST EXECUTABLE STATEMENT DVSTEP
! ..
        KFLAG = 0
        TOLD = TN
        NCF = 0
        JCUR = 0
        NFLAG = 0
        if (JSTART>0) goto 10
        if (JSTART==-1) goto 30

! On the first call, the order is set to 1, and other variables are
! initialized. ETAMAX is the maximum ratio by which H can be increased
! in a single step. It is normally 10, but is larger during the first
! step to compensate for the small initial H. If a failure occurs
! (in corrector convergence or error test), ETAMAX is set to 1 for
! the next increase.

        LMAX = MAXORD + 1
        NQ = 1
        L = 2
        NQNYH = NQ*LDYH
        TAU(1) = H
        PRL1 = ONE
        RC = ZERO
        ETAMAX = ETAMX1
        NQWAIT = 2
        HSCAL = H
        goto 70

! Take preliminary actions on a normal continuation step (JSTART > 0).
! If the driver changed H, then ETA must be reset and NEWH set to 1.
! If a change of order was dictated on the previous step, then it is
! done here and appropriate adjustments in the history are made.
! On an order decrease, the history array is adjusted by DVJUST.
! On an order increase, the history array is augmented by a column.
! On a change of step size H, the history array YH is rescaled.

10      continue
        if (KUTH==1) then
          ETA = min(ETA,H/HSCAL)
          NEWH = 1
        end if
20      if (NEWH==0) goto 70
        if (NEWQ==NQ) goto 60
        if (NEWQ<NQ) then
          call DVJUST(YH,LDYH,-1)
          NQ = NEWQ
          L = NQ + 1
          NQWAIT = L
          goto 60
        end if
        if (NEWQ>NQ) then
          call DVJUST(YH,LDYH,1)
          NQ = NEWQ
          L = NQ + 1
          NQWAIT = L
          goto 60
        end if

! The following block handles preliminaries needed when JSTART = -1.
! If N was reduced, zero out part of YH to avoid undefined references.
! If MAXORD was reduced to a value less than the tentative order NEWQ,
! then NQ is set to MAXORD, and a new H ratio ETA is chosen.
! Otherwise, we take the same preliminary actions as for JSTART > 0.
! In any case, NQWAIT is reset to L = NQ + 1 to prevent further
! changes in order for that many steps. The new H ratio ETA is
! limited by the input H if KUTH = 1, by HMIN if KUTH = 0, and by
! HMXI in any case. Finally, the history array YH is rescaled.

30      continue
        LMAX = MAXORD + 1
        if (N==LDYH) goto 40
        I1 = 1 + (NEWQ+1)*LDYH
        I2 = (MAXORD+1)*LDYH
        if (I1>I2) goto 40
        YH1(I1:I2) = ZERO
40      if (NEWQ<=MAXORD) goto 50
        FLOTL = real(LMAX)
        if (MAXORD<NQ-1) then
          DDN = DVNORM(N,SAVF,EWT)/TQ(1)
          ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL)+ADDON)
        end if
        if (MAXORD==NQ .and. NEWQ==NQ+1) ETA = ETAQ
        if (MAXORD==NQ-1 .and. NEWQ==NQ+1) then
          ETA = ETAQM1
          call DVJUST(YH,LDYH,-1)
        end if
        if (MAXORD==NQ-1 .and. NEWQ==NQ) then
          DDN = DVNORM(N,SAVF,EWT)/TQ(1)
          ETA = ONE/((BIAS1*DDN)**(ONE/FLOTL)+ADDON)
          call DVJUST(YH,LDYH,-1)
        end if
        ETA = min(ETA,ONE)
        NQ = MAXORD
        L = LMAX
50      if (KUTH==1) ETA = min(ETA,abs(H/HSCAL))
        if (KUTH==0) ETA = max(ETA,HMIN/abs(HSCAL))
        ETA = ETA/max(ONE,abs(HSCAL)*HMXI*ETA)
        NEWH = 1
        NQWAIT = L
        if (NEWQ<=MAXORD) goto 20
!       Rescale the history array for a change in H by a factor of ETA.
60      R = ONE
        do J = 2, L
          R = R*ETA
! Original:
!         CALL DSCAL_F90(N,R,YH(1,J),1)
          call DSCAL_F90(N,R,YH(1:N,J),1)
        end do
        H = HSCAL*ETA
        HSCAL = H
        RC = RC*ETA
        NQNYH = NQ*LDYH

! This section computes the predicted values by effectively
! multiplying the YH array by the Pascal triangle matrix.
! DVSET is called to calculate all integration coefficients.
! RC is the ratio of new to old values of the coefficient
! H/EL(2)=h/l1.

70      TN = TN + H
        I1 = NQNYH + 1
        do JB = 1, NQ
          I1 = I1 - LDYH
          do I = I1, NQNYH
            YH1(I) = YH1(I) + YH1(I+LDYH)
          end do
        end do
        call DVSET
        RL1 = ONE/EL(2)
        RC = RC*(RL1/PRL1)
        PRL1 = RL1

!       Call the nonlinear system solver.

        call VNLS(Y,YH,LDYH,SAVF,EWT,ACOR,IWM,WM,F,JAC,NFLAG, &
          ATOL,ITOL)

        if (NFLAG==0) goto 80

! The DVNLSD routine failed to achieve convergence (NFLAG /= 0).
! The YH array is retracted to its values before prediction.
! The step size H is reduced and the step is retried, if possible.
! Otherwise, an error exit is taken.

        NCF = NCF + 1
        NCFN = NCFN + 1
        ETAMAX = ONE
        TN = TOLD
        I1 = NQNYH + 1
        do JB = 1, NQ
          I1 = I1 - LDYH
          do I = I1, NQNYH
            YH1(I) = YH1(I) - YH1(I+LDYH)
          end do
        end do
        if (NFLAG<-1) goto 240
        if (abs(H)<=HMIN*ONEPSM) goto 230
!       IF (NCF==MXNCF) GOTO 230
        if (NCF==CONSECUTIVE_CFAILS) goto 230
        ETA = ETACF
        ETA = max(ETA,HMIN/abs(H))
        NFLAG = -1
        goto 60

! The corrector has converged (NFLAG = 0). The local error test is
! made and control passes to statement 500 if it fails.

80      continue
        DSM = ACNRM/TQ(2)
        if (DSM>ONE) goto 100

! After a successful step, update the YH and TAU arrays and decrement
! NQWAIT. If NQWAIT is then 1 and NQ < MAXORD, then ACOR is saved
! for use in a possible order increase on the next step.
! If ETAMAX = 1 (a failure occurred this step), keep NQWAIT >= 2.

        KFLAG = 0
        NST = NST + 1
        HU = H
        NQU = NQ
        do IBACK = 1, NQ
          I = L - IBACK
          TAU(I+1) = TAU(I)
        end do
        TAU(1) = H

        if (BOUNDS) then
! Original:
!         CALL DAXPY_F90(N,EL(1),ACOR,1,YH(1,1),1)
!         CALL DAXPY_F90(N,EL(2),ACOR,1,YH(1,2),1)
          call DAXPY_F90(N,EL(1),ACOR,1,YH(1:N,1),1)
          call DAXPY_F90(N,EL(2),ACOR,1,YH(1:N,2),1)
!         Take care of roundoff causing y(t) to be slightly unequal
!         to the constraint bound.
          do J = 1, NDX
            YH(IDX(J),1) = max(YH(IDX(J),1),LB(J))
            Y(IDX(J)) = max(Y(IDX(J)),LB(J))
            YH(IDX(J),1) = min(YH(IDX(J),1),UB(J))
            Y(IDX(J)) = min(Y(IDX(J)),UB(J))
          end do
!         Update the higher derivatives and project to zero if necessary.
          if (L>2) then
            do J = 3, L
! Original:
!             CALL DAXPY_F90(N,EL(J),ACOR,1,YH(1,J),1)
              call DAXPY_F90(N,EL(J),ACOR,1,YH(1:N,J),1)
            end do
          end if
        else
!         Proceed as usual.
          do J = 1, L
! Original:
!           CALL DAXPY_F90(N,EL(J),ACOR,1,YH(1,J),1)
            call DAXPY_F90(N,EL(J),ACOR,1,YH(1:N,J),1)
          end do
        end if

        NQWAIT = NQWAIT - 1
        if ((L==LMAX) .or. (NQWAIT/=1)) goto 90
! Original:
!       CALL DCOPY_F90(N,ACOR,1,YH(1,LMAX),1)
        call DCOPY_F90(N,ACOR,1,YH(1:N,LMAX),1)
        CONP = TQ(5)
90      if (abs(ETAMAX-ONE)>0) goto 130
        if (NQWAIT<2) NQWAIT = 2
        NEWQ = NQ
        NEWH = 0
        ETA = ONE
        HNEW = H
        goto 250

! The error test failed. KFLAG keeps track of multiple failures.
! Restore TN and the YH array to their previous values, and prepare
! to try the step again. Compute the optimum step size for the
! same order. After repeated failures, H is forced to decrease
! more rapidly.

100     KFLAG = KFLAG - 1
        NETF = NETF + 1
        NFLAG = -2
        TN = TOLD
        I1 = NQNYH + 1
        do JB = 1, NQ
          I1 = I1 - LDYH
          do I = I1, NQNYH
            YH1(I) = YH1(I) - YH1(I+LDYH)
          end do
        end do
        if (abs(H)<=HMIN*ONEPSM) goto 220
        ETAMAX = ONE
        if (KFLAG<=KFC) goto 110
!       Compute ratio of new H to current H at the current order.
        FLOTL = real(L)
        ETA = ONE/((BIAS2*DSM)**(ONE/FLOTL)+ADDON)
        ETA = max(ETA,HMIN/abs(H),ETAMIN)
        if ((KFLAG<=-2) .and. (ETA>ETAMXF)) ETA = ETAMXF
        goto 60

! Control reaches this section if 3 or more consecutive failures
! have occurred. It is assumed that the elements of the YH array
! have accumulated errors of the wrong order. The order is reduced
! by one, if possible. Then H is reduced by a factor of 0.1 and
! the step is retried. After a total of 7 consecutive failures,
! an exit is taken with KFLAG = -1.

!110     IF (KFLAG==KFH) GOTO 220
110     if (KFLAG==CONSECUTIVE_EFAILS) goto 220
        if (NQ==1) goto 120
        ETA = max(ETAMIN,HMIN/abs(H))
        call DVJUST(YH,LDYH,-1)
        L = NQ
        NQ = NQ - 1
        NQWAIT = L
        goto 60
120     ETA = max(ETAMIN,HMIN/abs(H))
        H = H*ETA
        HSCAL = H
        TAU(1) = H
        call F(N,TN,Y,SAVF)
        NFE = NFE + 1
        if (BOUNDS) then
          do I = 1, NDX
            if (abs(YH(IDX(I),1)-LB(I))<=ZERO) SAVF(IDX(I)) = &
              max(SAVF(IDX(I)),ZERO)
            if (abs(YH(IDX(I),1)-UB(I))<=ZERO) SAVF(IDX(I)) = &
              min(SAVF(IDX(I)),ZERO)
          end do
        end if
        YH(1:N,2) = H*SAVF(1:N)
        NQWAIT = 10
        goto 70

! If NQWAIT = 0, an increase or decrease in order by one is considered.
! Factors ETAQ, ETAQM1, ETAQP1 are computed by which H could be
! multiplied at order q, q-1, or q+1, respectively. The largest of
! these is determined, and the new order and step size set accordingly.
! A change of H or NQ is made only if H increases by at least a factor
! of THRESH. If an order change is considered and rejected, then NQWAIT
! is set to 2 (reconsider it after 2 steps).

!       Compute ratio of new H to current H at the current order.
130     FLOTL = real(L)
        ETAQ = ONE/((BIAS2*DSM)**(ONE/FLOTL)+ADDON)
        if (NQWAIT/=0) goto 170
        NQWAIT = 2
        ETAQM1 = ZERO
        if (NQ==1) goto 140
!       Compute ratio of new H to current H at the current order
!       less one.
! Original:
!       DDN = DVNORM(N,YH(1,L),EWT)/TQ(1)
        DDN = DVNORM(N,YH(1:N,L),EWT)/TQ(1)
        ETAQM1 = ONE/((BIAS1*DDN)**(ONE/(FLOTL-ONE))+ADDON)
140     ETAQP1 = ZERO
        if (L==LMAX) goto 150
!       Compute ratio of new H to current H at current order plus one.
        CNQUOT = (TQ(5)/CONP)*(H/TAU(2))**L
        SAVF(1:N) = ACOR(1:N) - CNQUOT*YH(1:N,LMAX)
        DUP = DVNORM(N,SAVF,EWT)/TQ(3)
        ETAQP1 = ONE/((BIAS3*DUP)**(ONE/(FLOTL+ONE))+ADDON)
150     if (ETAQ>=ETAQP1) goto 160
        if (ETAQP1>ETAQM1) goto 190
        goto 180
160     if (ETAQ<ETAQM1) goto 180
170     ETA = ETAQ
        NEWQ = NQ
        goto 200
180     ETA = ETAQM1
        NEWQ = NQ - 1
        goto 200
190     ETA = ETAQP1
        NEWQ = NQ + 1
! Original:
!       CALL DCOPY_F90(N,ACOR,1,YH(1,LMAX),1)
        call DCOPY_F90(N,ACOR,1,YH(1:N,LMAX),1)
!       Test tentative new H against THRESH, ETAMAX, and HMXI and
!       then exit.
200     if (ETA<THRESH .or. abs(ETAMAX-ONE)<=ZERO) goto 210
        ETA = min(ETA,ETAMAX)
        ETA = ETA/max(ONE,abs(H)*HMXI*ETA)
        NEWH = 1
        HNEW = H*ETA
        goto 250
210     NEWQ = NQ
        NEWH = 0
        ETA = ONE
        HNEW = H
        goto 250

! All returns are made through this section.
! On a successful return, ETAMAX is reset and ACOR is scaled.

220     KFLAG = -1
        goto 260
230     KFLAG = -2
        goto 260
240     if (NFLAG==-2) KFLAG = -3
        if (NFLAG==-3) KFLAG = -4
        goto 260
250     ETAMAX = ETAMX3
        if (NST<=10) ETAMAX = ETAMX2
        R = ONE/TQ(2)

        call DSCAL_F90(N,R,ACOR,1)
260     JSTART = 1
        return

      end subroutine DVSTEP
!_______________________________________________________________________

      subroutine DVSET
! ..
! Set the integration coefficients for DVSTEP.
! ..
! For each order NQ, the coefficients in EL are calculated by use
! of the generating polynomial lambda(x), with coefficients EL(i).
!      lambda(x) = EL(1) + EL(2)*x + ... + EL(NQ+1)*(x**NQ).
! For the backward differentiation formulas,
!                                     NQ-1
!      lambda(x) = (1 + x/xi*(NQ)) * product (1 + x/xi(i)) .
!                                     i = 1
! For the Adams formulas,
!                              NQ-1
!      (d/dx) lambda(x) = c * product (1 + x/xi(i)),
!                              i = 1
!      lambda(-1) = 0,    lambda(0) = 1,
! where c is a normalization constant.
! In both cases, xi(i) is defined by
!      H*xi(i) = t sub n - t sub (n-i)
!              = H + TAU(1) + TAU(2) + ... TAU(i-1).
! In addition to variables described previously, communication
! with DVSET uses the following:
!   TAU    = A vector of length 13 containing the past NQ values
!            of H.
!   EL     = A vector of length 13 in which vset stores the
!            coefficients for the corrector formula.
!   TQ     = A vector of length 5 in which vset stores constants
!            used for the convergence test, the error test, and the
!            selection of H at a new order.
!   METH   = The basic method indicator.
!   NQ     = The current order.
!   L      = NQ + 1, the length of the vector stored in EL, and
!            the number of columns of the YH array being used.
!   NQWAIT = A counter controlling the frequency of order changes.
!            An order change is about to be considered if NQWAIT = 1.
! ..
     implicit none
! ..
! .. Local Scalars ..
        real (WP) :: AHATN0, ALPH0, CNQM1, CSUM, ELP, EM0, FLOTI, FLOTL, &
          FLOTNQ, HSUM, RXI, RXIS, S, T1, T2, T3, T4, T5, T6, XI
        integer :: I, IBACK, J, JP1, NQM1, NQM2
! ..
! .. Local Arrays ..
        real (WP) :: EM(13)
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, real
! ..
! .. FIRST EXECUTABLE STATEMENT DVSET
! ..
        FLOTL = real(L)
        NQM1 = NQ - 1
        NQM2 = NQ - 2
        goto (10,40) METH

!       Set coefficients for Adams methods.
10      if (NQ/=1) goto 20
        EL(1) = ONE
        EL(2) = ONE
        TQ(1) = ONE
        TQ(2) = TWO
        TQ(3) = SIX*TQ(2)
        TQ(5) = ONE
        goto 60
20      HSUM = H
        EM(1) = ONE
        FLOTNQ = FLOTL - ONE
        EM(2:L) = ZERO
        do J = 1, NQM1
          if ((J/=NQM1) .or. (NQWAIT/=1)) goto 30
          S = ONE
          CSUM = ZERO
          do I = 1, NQM1
            CSUM = CSUM + S*EM(I)/real(I+1)
            S = -S
          end do
          TQ(1) = EM(NQM1)/(FLOTNQ*CSUM)
30        RXI = H/HSUM
          do IBACK = 1, J
            I = (J+2) - IBACK
            EM(I) = EM(I) + EM(I-1)*RXI
          end do
          HSUM = HSUM + TAU(J)
        end do
!       Compute integral from -1 to 0 of polynomial and of x times it.
        S = ONE
        EM0 = ZERO
        CSUM = ZERO
        do I = 1, NQ
          FLOTI = real(I)
          EM0 = EM0 + S*EM(I)/FLOTI
          CSUM = CSUM + S*EM(I)/(FLOTI+ONE)
          S = -S
        end do
!       In EL, form coefficients of normalized integrated polynomial.
        S = ONE/EM0
        EL(1) = ONE
        do I = 1, NQ
          EL(I+1) = S*EM(I)/real(I)
        end do
        XI = HSUM/H
        TQ(2) = XI*EM0/CSUM
        TQ(5) = XI/EL(L)
        if (NQWAIT/=1) goto 60
!       For higher order control constant, multiply polynomial by
!       1+x/xi(q).
        RXI = ONE/XI
        do IBACK = 1, NQ
          I = (L+1) - IBACK
          EM(I) = EM(I) + EM(I-1)*RXI
        end do
!       Compute integral of polynomial.
        S = ONE
        CSUM = ZERO
        do I = 1, L
          CSUM = CSUM + S*EM(I)/real(I+1)
          S = -S
        end do
        TQ(3) = FLOTL*EM0/CSUM
        goto 60

!       Set coefficients for BDF methods.
40      EL(3:L) = ZERO
        EL(1) = ONE
        EL(2) = ONE
        ALPH0 = -ONE
        AHATN0 = -ONE
        HSUM = H
        RXI = ONE
        RXIS = ONE
        if (NQ==1) goto 50
        do J = 1, NQM2
!       In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)).
          HSUM = HSUM + TAU(J)
          RXI = H/HSUM
          JP1 = J + 1
          ALPH0 = ALPH0 - ONE/real(JP1)
          do IBACK = 1, JP1
            I = (J+3) - IBACK
            EL(I) = EL(I) + EL(I-1)*RXI
          end do
        end do
        ALPH0 = ALPH0 - ONE/real(NQ)
        RXIS = -EL(2) - ALPH0
        HSUM = HSUM + TAU(NQM1)
        RXI = H/HSUM
        AHATN0 = -EL(2) - RXI
        do IBACK = 1, NQ
          I = (NQ+2) - IBACK
          EL(I) = EL(I) + EL(I-1)*RXIS
        end do
50      T1 = ONE - AHATN0 + ALPH0
        T2 = ONE + real(NQ)*T1
        TQ(2) = abs(ALPH0*T2/T1)
        TQ(5) = abs(T2/(EL(L)*RXI/RXIS))
        if (NQWAIT/=1) goto 60
        CNQM1 = RXIS/EL(L)
        T3 = ALPH0 + ONE/real(NQ)
        T4 = AHATN0 + RXI
        ELP = T3/(ONE-T4+T3)
        TQ(1) = abs(ELP/CNQM1)
        HSUM = HSUM + TAU(NQ)
        RXI = H/HSUM
        T5 = ALPH0 - ONE/real(NQ+1)
        T6 = AHATN0 - RXI
        ELP = T2/(ONE-T6+T5)
        TQ(3) = abs(ELP*RXI*(FLOTL+ONE)*T5)
60      TQ(4) = CORTES*TQ(2)
        return

      end subroutine DVSET
!_______________________________________________________________________

      subroutine DVJUST(YH,LDYH,IORD)
! ..
! Adjust the Nordsieck array.
! ..
! This subroutine adjusts the YH array on reduction of order, and
! also when the order is increased for the stiff option (METH = 2).
! Communication with DVJUST uses the following:
! IORD  = An integer flag used when METH = 2 to indicate an order
!         increase (IORD = +1) or an order decrease (IORD = -1).
! HSCAL = Step size H used in scaling of Nordsieck array YH.
!         (If IORD = +1, DVJUST assumes that HSCAL = TAU(1).)
! See References 1 and 2 for details.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: IORD, LDYH
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: YH(LDYH,*)
! ..
! .. Local Scalars ..
        real (WP) :: ALPH0, ALPH1, HSUM, PROD, T1, XI, XIOLD
        integer :: I, IBACK, J, JP1, LP1, NQM1, NQM2, NQP1
! ..
! .. Intrinsic Functions ..
        intrinsic real
! ..
! .. FIRST EXECUTABLE STATEMENT DVJUST
! ..
        if ((NQ==2) .and. (IORD/=1)) return
        NQM1 = NQ - 1
        NQM2 = NQ - 2
        goto (10,30) METH

!       Nonstiff option.

!       Check to see if the order is being increased or decreased.
10      continue
        if (IORD==1) goto 20
!       Order decrease.
        EL(1:LMAX) = ZERO
        EL(2) = ONE
        HSUM = ZERO
        do J = 1, NQM2
!         Construct coefficients of x*(x+xi(1))*...*(x+xi(j)).
          HSUM = HSUM + TAU(J)
          XI = HSUM/HSCAL
          JP1 = J + 1
          do IBACK = 1, JP1
            I = (J+3) - IBACK
            EL(I) = EL(I)*XI + EL(I-1)
          end do
        end do
!       Construct coefficients of integrated polynomial.
        do J = 2, NQM1
          EL(J+1) = real(NQ)*EL(J)/real(J)
        end do
!       Subtract correction terms from YH array.
        do J = 3, NQ
          do I = 1, N
            YH(I,J) = YH(I,J) - YH(I,L)*EL(J)
          end do
        end do
        return
!       Order increase.
!       Zero out next column in YH array.
20      continue
        LP1 = L + 1
        YH(1:N,LP1) = ZERO
        return

!       Stiff option.

!       Check to see if the order is being increased or decreased.
30      continue
        if (IORD==1) goto 40
!       Order decrease.
        EL(1:LMAX) = ZERO
        EL(3) = ONE
        HSUM = ZERO
        do J = 1, NQM2
!     Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)).
          HSUM = HSUM + TAU(J)
          XI = HSUM/HSCAL
          JP1 = J + 1
          do IBACK = 1, JP1
            I = (J+4) - IBACK
            EL(I) = EL(I)*XI + EL(I-1)
          end do
        end do
!       Subtract correction terms from YH array.
        do J = 3, NQ
          YH(1:N,J) = YH(1:N,J) - YH(1:N,L)*EL(J)
        end do
        return
!       Order increase.
40      EL(1:LMAX) = ZERO
        EL(3) = ONE
        ALPH0 = -ONE
        ALPH1 = ONE
        PROD = ONE
        XIOLD = ONE
        HSUM = HSCAL
        if (NQ==1) goto 50
        do J = 1, NQM1
!       Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)).
          JP1 = J + 1
          HSUM = HSUM + TAU(JP1)
          XI = HSUM/HSCAL
          PROD = PROD*XI
          ALPH0 = ALPH0 - ONE/real(JP1)
          ALPH1 = ALPH1 + ONE/XI
          do IBACK = 1, JP1
            I = (J+4) - IBACK
            EL(I) = EL(I)*XIOLD + EL(I-1)
          end do
          XIOLD = XI
        end do
50      continue
        T1 = (-ALPH0-ALPH1)/PROD
!       Load column L+1 in YH array.
        LP1 = L + 1
        YH(1:N,LP1) = T1*YH(1:N,LMAX)
!       Add correction terms to YH array.
        NQP1 = NQ + 1
        do J = 3, NQP1
! Original:
!         CALL DAXPY_F90(N,EL(J),YH(1,LP1),1,YH(1,J),1)
          call DAXPY_F90(N,EL(J),YH(1:N,LP1),1,YH(1:N,J),1)
        end do
        return

      end subroutine DVJUST
!_______________________________________________________________________

      subroutine DVNLSD(Y,YH,LDYH,SAVF,EWT,ACOR,IWM,WM,F,JAC,NFLAG,&
        ATOL,ITOL)
! ..
! This is the nonlinear system solver for dense and banded solutions.
! ..
! Subroutine DVNLSD is a nonlinear system solver which uses functional
! iteration or a chord (modified Newton) method. For the chord method
! direct linear algebraic system solvers are used. Subroutine DVNLSD
! then handles the corrector phase of this integration package.
! Communication with DVNLSD is done with the following variables. (For
! more details, please see the comments in the driver subroutine.)
! Y          = The dependent variable, a vector of length N, input.
! YH         = The Nordsieck (Taylor) array, LDYH by LMAX, input
!              and output. On input, it contains predicted values.
! LDYH       = A constant >= N, the first dimension of YH, input.
! SAVF       = A work array of length N.
! EWT        = An error weight vector of length N, input.
! ACOR       = A work array of length N, used for the accumulated
!              corrections to the predicted y vector.
! WM,IWM     = Real and integer work arrays associated with matrix
!              operations in chord iteration (MITER /= 0).
! F          = Dummy name for user supplied routine for f.
! JAC        = Dummy name for user supplied Jacobian routine.
! NFLAG      = Input/output flag, with values and meanings as follows:
!              INPUT
!                  0 first call for this time step.
!                 -1 convergence failure in previous call to DVNLSD.
!                 -2 error test failure in DVSTEP.
!              OUTPUT
!                  0 successful completion of nonlinear solver.
!                 -1 convergence failure or singular matrix.
!                 -2 unrecoverable error in matrix preprocessing
!                    (cannot occur here).
!                 -3 unrecoverable error in solution (cannot occur
!                    here).
! IPUP       = Own variable flag with values and meanings as follows:
!              0,          do not update the Newton matrix.
!              MITER /= 0  update Newton matrix, because it is the
!                          initial step, order was changed, the error
!                          test failed, or an update is indicated by
!                          the scalar RC or step counter NST.
! For more details, see comments in driver subroutine.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: ITOL, LDYH
        integer, intent (INOUT) :: NFLAG
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: ACOR(*), EWT(*), SAVF(*), &
          WM(*), Y(*), YH(LDYH,*)
        real (WP), intent (IN) :: ATOL(*)
        integer, intent (INOUT) :: IWM(*)
! ..
! .. Subroutine Arguments ..
        external F, JAC
! ..
! .. Local Scalars ..
        real (WP) :: ACNRMNEW, CSCALE, DCON, DEL, DELP
        integer :: I, IERPJ, IERSL, M
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, MAX, MIN
! ..
! .. FIRST EXECUTABLE STATEMENT DVNLSD
! ..
! On the first step, on a change of method order, or after a
! nonlinear convergence failure with NFLAG = -2, set IPUP = MITER
! to force a Jacobian update when MITER /= 0.
        if (JSTART==0) NSLP = 0
        if (NFLAG==0) ICF = 0
        if (NFLAG==-2) IPUP = MITER
        if ((JSTART==0) .or. (JSTART==-1)) IPUP = MITER
!     If this is functional iteration, set CRATE = 1 and drop to 220
        if (MITER==0) then
          CRATE = ONE
          goto 10
        end if

! RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
! When RC differs from 1 by more than CCMAX, IPUP is set to MITER to
! force DVJAC to be called, if a Jacobian is involved. In any case,
! DVJAC is called at least every MSBP steps.

        DRC = abs(RC-ONE)
        if (DRC>CCMAX .or. NST>=NSLP+MSBP) IPUP = MITER

! Up to MAXCOR corrector iterations are taken. A convergence test is
! made on the r.m.s. norm of each correction, weighted by the error
! weight vector EWT. The sum of the corrections is accumulated in the
! vector ACOR(i). The YH array is not altered in the corrector loop.

10      M = 0
        DELP = ZERO
! Original:
!       CALL DCOPY_F90(N,YH(1,1),1,Y,1)
        call DCOPY_F90(N,YH(1:N,1),1,Y(1:N),1)
        call F(N,TN,Y,SAVF)
        NFE = NFE + 1
        if (BOUNDS) then
          do I = 1, NDX
            if (abs(YH(IDX(I),1)-LB(I))<=ZERO) SAVF(IDX(I)) = &
               max(SAVF(IDX(I)),ZERO)
            if (abs(YH(IDX(I),1)-UB(I))<=ZERO) SAVF(IDX(I)) = &
               min(SAVF(IDX(I)),ZERO)
          end do
        end if
        if (IPUP<=0) goto 20

! If indicated, the matrix P = I - h*rl1*J is reevaluated and
! preprocessed before starting the corrector iteration. IPUP
! is set to 0 as an indicator that this has been done.

        call DVJAC(Y,YH,LDYH,EWT,ACOR,SAVF,WM,IWM,F,JAC,IERPJ, &
          ATOL,ITOL)
        IPUP = 0
        RC = ONE
        DRC = ZERO
        CRATE = ONE
        NSLP = NST
!       If matrix is singular, take error return to force cut in
!       step size.
        if (IERPJ/=0) goto 70
20      ACOR(1:N) = ZERO
!       This is a looping point for the corrector iteration.
30      if (MITER/=0) goto 40

! In the case of functional iteration, update Y directly from
! the result of the last function evaluation.

        SAVF(1:N) = RL1*(H*SAVF(1:N)-YH(1:N,2))
        Y(1:N) = SAVF(1:N) - ACOR(1:N)
        DEL = DVNORM(N,Y,EWT)
        Y(1:N) = YH(1:N,1) + SAVF(1:N)
        call DCOPY_F90(N,SAVF,1,ACOR,1)
        goto 50

! In the case of the chord method, compute the corrector error, and
! solve the linear system with that as right-hand side and P as
! coefficient matrix. The correction is scaled by the factor
! 2/(1+RC) to account for changes in h*rl1 since the last DVJAC call.

40      Y(1:N) = (RL1*H)*SAVF(1:N) - (RL1*YH(1:N,2)+ACOR(1:N))
        call DVSOL(WM,IWM,Y,IERSL)
        NNI = NNI + 1
        if (IERSL>0) goto 60
        if (METH==2 .and. abs(RC-ONE)>ZERO) then
          CSCALE = TWO/(ONE+RC)
          call DSCAL_F90(N,CSCALE,Y,1)
        end if
        DEL = DVNORM(N,Y,EWT)
        call DAXPY_F90(N,ONE,Y,1,ACOR,1)
        Y(1:N) = YH(1:N,1) + ACOR(1:N)

! Test for convergence. If M > 0, an estimate of the convergence
! rate constant is stored in CRATE, and this is used in the test.

50      if (M/=0) CRATE = max(CRDOWN*CRATE,DEL/DELP)
        DCON = DEL*min(ONE,CRATE)/TQ(4)
        if (DCON<=ONE) goto 80
        M = M + 1
        if (M==MAXCOR) goto 60
        if (M>=2 .and. DEL>RDIV*DELP) goto 60
        DELP = DEL
        call F(N,TN,Y,SAVF)
        NFE = NFE + 1
        if (BOUNDS) then
          do I = 1, NDX
            if (abs(YH(IDX(I),1)-LB(I))<=ZERO) SAVF(IDX(I)) = &
               max(SAVF(IDX(I)),ZERO)
            if (abs(YH(IDX(I),1)-UB(I))<=ZERO) SAVF(IDX(I)) = &
               min(SAVF(IDX(I)),ZERO)
          end do
        end if
        goto 30

60      if (MITER==0 .or. JCUR==1) goto 70
        ICF = 1
        IPUP = MITER
        goto 10

70      continue
        NFLAG = -1
        ICF = 2
        IPUP = MITER
        return

!       Return for successful step.
80      continue

!       Enforce bounds.
        if (BOUNDS) then
          CHANGED_ACOR = .false.
          if (M==0) then
            ACNRM = DEL
          else
            ACNRM = DVNORM(N,ACOR,EWT)
          end if
          if (MITER/=0) then
!           Since Y(:) = YH(:,1) + ACOR(:) ...
            do I = 1, NDX
              if (Y(IDX(I))<LB(I)) then
                CHANGED_ACOR = .true.
                ACOR(IDX(I)) = LB(I) - YH(IDX(I),1)
                SAVF(IDX(I)) = ACOR(IDX(I))
              end if
              if (Y(IDX(I))>UB(I)) then
                CHANGED_ACOR = .true.
                ACOR(IDX(I)) = UB(I) - YH(IDX(I),1)
                SAVF(IDX(I)) = ACOR(IDX(I))
              end if
            end do
          else
!           Since Y(:) = YH(:,1) + SAVF(:) and
!           since CALL DCOPY_F90(N,SAVF,1,ACOR,1) ...
            do I = 1, NDX
              if (Y(IDX(I))<LB(IDX(I))) then
                CHANGED_ACOR = .true.
                ACOR(IDX(I)) = LB(I) - YH(IDX(I),1)
              end if
              if (Y(IDX(I))>UB(IDX(I))) then
                CHANGED_ACOR = .true.
                ACOR(IDX(I)) = UB(I) - YH(IDX(I),1)
              end if
            end do
          end if
          if (CHANGED_ACOR) then
            if (M==0) then
              ACNRMNEW = DEL
            else
              ACNRMNEW = DVNORM(N,ACOR,EWT)
            end if
!           ACNRM = ACNRMNEW
            ACNRM = max(ACNRM,ACNRMNEW)
          else
          end if
          NFLAG = 0
          JCUR = 0
          ICF = 0
        else
!         No projections are required.
          NFLAG = 0
          JCUR = 0
          ICF = 0
          if (M==0) ACNRM = DEL
          if (M>0) ACNRM = DVNORM(N,ACOR,EWT)
        end if
        return

      end subroutine DVNLSD
!_______________________________________________________________________

      subroutine DVJAC(Y,YH,LDYH,EWT,FTEM,SAVF,WM,IWM,F,JAC,IERPJ, &
        ATOL,ITOL)
! ..
! Compute and process the matrix P = I - h*rl1*J, where J is an
! approximation to the Jacobian for dense and banded solutions.
! ..
! This is a version of DVJAC that allows use of the known nonzero
! diagonals if it is available.
! DVJAC is called by DVNLSD to compute and process the matrix
! P = I - h*rl1*J, where J is an approximation to the Jacobian.
! Here J is computed by the user-supplied routine JAC if
! MITER = 1 or 4, or by finite differencing if MITER = 2, 3, or 5.
! If MITER = 3, a diagonal approximation to J is used.
! If JSV = -1, J is computed from scratch in all cases.
! If JSV = 1 and MITER = 1, 2, 4, or 5, and if the saved value of J is
! considered acceptable, then P is constructed from the saved J.
! J is stored in wm and replaced by P. If MITER /= 3, P is then
! subjected to LU decomposition in preparation for later solution
! of linear systems with P as coefficient matrix. This is done
! by DGEFA_F90 if MITER = 1 or 2, and by DGBFA_F90 if MITER = 4 or 5.
! Communication with DVJAC is done with the following variables.
! (For more details, please see the comments in the driver subroutine.)
! Y          = Vector containing predicted values on entry.
! YH         = The Nordsieck array, an LDYH by LMAX array, input.
! LDYH       = A constant >= N, the first dimension of YH, input.
! EWT        = An error weight vector of length N.
! SAVF       = Array containing f evaluated at predicted y, input.
! WM         = Real work space for matrices. In the output, it contains
!              the inverse diagonal matrix if MITER = 3 and the LU
!              decomposition of P if MITER is 1, 2, 4, or 5.
!              Storage of matrix elements starts at WM(1).
!              Storage of the saved Jacobian starts at WM(LOCJS).
! IWM        = Integer work space containing pivot information,
!              starting at IWM(31), if MITER is 1, 2, 4, or 5.
!              IWM also contains band parameters ML = IWM(1) and
!              MU = IWM(2) if MITER is 4 or 5.
! F          = Dummy name for the user supplied subroutine for f.
! JAC        = Dummy name for the user supplied Jacobian subroutine.
! RL1        = 1/EL(2) (input).
! IERPJ      = Output error flag, = 0 if no trouble, 1 if the P
!              matrix is found to be singular.
! JCUR       = Output flag to indicate whether the Jacobian matrix
!              (or approximation) is now current.
!              JCUR = 0 means J is not current.
!              JCUR = 1 means J is current.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (INOUT) :: IERPJ
        integer, intent (IN) :: LDYH, ITOL
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: EWT(*), FTEM(*), SAVF(*), WM(*), &
          Y(*),YH(LDYH,*)
        real (WP), intent (IN) :: ATOL(*)
        integer, intent (INOUT) :: IWM(*)
! ..
! .. Subroutine Arguments ..
        external F, JAC
! ..
! .. Local Scalars ..
        real (WP) :: CON, DI, FAC, HRL1, R, R0, SRUR, YI, YJ, YJJ
        integer :: I, I1, I2, IER, II, J, J1, JJ, JJ1, JJ2, JOK, K,   &
          K1, K2, LENP, MBA, MBAND, MEB1, MEBAND, ML, ML1, MU, NG, NP1
!         K1, K2, LENP, MBA, MBAND, MEB1, MEBAND, ML, ML3, MU, NG, NP1
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, MAX, MIN, real
! ..
! .. FIRST EXECUTABLE STATEMENT DVJAC
! ..
        IERPJ = 0
        HRL1 = H*RL1
!       See whether J should be evaluated (JOK = -1) or not (JOK = 1).
        JOK = JSV
        if (JSV==1) then
          if (NST==0 .or. NST>NSLJ+MSBJ) JOK = -1
          if (ICF==1 .and. DRC<CCMXJ) JOK = -1
          if (ICF==2) JOK = -1
        end if
        if (J_IS_CONSTANT .and. J_HAS_BEEN_COMPUTED) JOK = 1
!       End of setting JOK.

        if (JOK==-1 .and. MITER==1) then
!         If JOK = -1 and MITER = 1, call JAC to evaluate Jacobian.
          NJE = NJE + 1
          NSLJ = NST
          JCUR = 1
          LENP = N*N
!         WM(3:LENP+2) = ZERO
          WM(1:LENP) = ZERO
!         CALL JAC(N,TN,Y,0,0,WM(3),N)
          call JAC(N,TN,Y,0,0,WM(1),N)
          if (J_IS_CONSTANT) J_HAS_BEEN_COMPUTED = .true.
!         IF (JSV==1) CALL DCOPY_F90(LENP,WM(3),1,WM(LOCJS),1)
          if (JSV==1) call DCOPY_F90(LENP,WM(1),1,WM(LOCJS),1)
        end if

!       Set flag to indicate how the YSCALE vector will be set for
!       JACSP.
        LIKE_ORIGINAL_VODE = .false.

        if (JOK==-1 .and. MITER==2) then
           if (USE_JACSP) then
!             Approximate the Jacobian using Doug Salane's JACSP.
              NSLJ = NST
              JCUR = 1
              IOPTDS(1) = 0
              IOPTDS(2) = 0
              IOPTDS(3) = 1
              IOPTDS(5) = 0
!             INFORDS(4) was initialized in DVODE.
              LWKDS  = 3 * N
              LIWKDS = 50 + N
              NRFJACDS = N
              NCFJACDS = N
              MAXGRPDS = N
!             Calculate the YSCALEDS vector for JACSPDV.
              if (LIKE_ORIGINAL_VODE) then
                 FAC = DVNORM(N,SAVF,EWT)
!                JACSPDB multiplies YSCALEDS(*) BY UROUND**0.825:
!                R0 = THOU*ABS(H)*UROUND*REAL(N)*FAC
                 R0 = THOU*abs(H)*real(N)*FAC
                 if (abs(R0)<=ZERO) R0 = ONE
                 SRUR = WM1
                 do J = 1, N
!                   JACSPDB multiplies YSCALEDS(*) BY UROUND**0.825:
!                   R = MAX(ABS(Y(J)),R0/EWT(J))
                    R = max(abs(Y(J))/U325,(R0/EWT(J))*U125)
                    YSCALEDS(J) = R
                 end do
              else
                 if (ITOL == 1 .or. ITOL == 3) then
                    do J = 1, N
                       YSCALEDS(J) = max(abs(Y(J)),ATOL(1),UROUND)
                    end do
                 else
                    do J = 1, N
                       YSCALEDS(J) = max(abs(Y(J)),ATOL(J),UROUND)
                    end do
                 end if
              end if

              call JACSPDB(F,N,TN,Y,SAVF,WM(1),NRFJACDS, &
                YSCALEDS,FACDS,IOPTDS,WKDS,LWKDS,IWKDS,LIWKDS, &
                MAXGRPDS,NGRPDS,JPNTRDS,INDROWDS)
              NFE = NFE + IWKDS(7)
              NJE = NJE + 1
           else
!             If MITER = 2, make N calls to F to approximate the Jacobian.
              NSLJ = NST
              JCUR = 1
              FAC = DVNORM(N,SAVF,EWT)
              R0 = THOU*abs(H)*UROUND*real(N)*FAC
              if (abs(R0)<=ZERO) R0 = ONE
              SRUR = WM1
!             J1 = 2
              J1 = 0
              do J = 1, N
                 YJ = Y(J)
                 R = max(SRUR*abs(YJ),R0/EWT(J))
                 Y(J) = Y(J) + R
                 FAC = ONE/R
                 call F(N,TN,Y,FTEM)
                 do I = 1, N
                   WM(I+J1) = (FTEM(I)-SAVF(I))*FAC
                 end do
                 Y(J) = YJ
                 J1 = J1 + N
              end do
              NFE = NFE + N
              NJE = NJE + 1
           end if
           LENP = N*N
!          IF (JSV==1) CALL DCOPY_F90(LENP,WM(3),1,WM(LOCJS),1)
           if (JSV==1) call DCOPY_F90(LENP,WM(1),1,WM(LOCJS),1)
           if (J_IS_CONSTANT) J_HAS_BEEN_COMPUTED = .true.
        end if

        if (JOK==1 .and. (MITER==1 .or. MITER==2)) then
          JCUR = 0
          LENP = N*N
!         CALL DCOPY_F90(LENP,WM(LOCJS),1,WM(3),1)
          call DCOPY_F90(LENP,WM(LOCJS),1,WM(1),1)
        end if

        if (MITER==1 .or. MITER==2) then
!         Multiply Jacobian by scalar, add identity, and do LU
!         decomposition.
          CON = -HRL1
!         CALL DSCAL_F90(LENP,CON,WM(3),1)
          call DSCAL_F90(LENP,CON,WM(1),1)
!         J = 3
          J = 1
          NP1 = N + 1
          do I = 1, N
            WM(J) = WM(J) + ONE
            J = J + NP1
          end do
          NLU = NLU + 1
! ______________________________________________________________________

!         CALL DGEFA_F90(WM(3),N,N,IWM(31),IER)
          call DGEFA_F90(WM(1),N,N,IWM(31),IER)
          if (IER/=0) IERPJ = 1
! *****LAPACK build change point. Replace above with these statements.
!        IF (.NOT.USE_LAPACK) THEN
!!         CALL DGEFA_f90(WM(3),N,N,IWM(31),IER)
!          CALL DGEFA_f90(WM(1),N,N,IWM(31),IER)
!          IF (IER /= 0) IERPJ = 1
!        ELSE
!!         CALL DGETRF(N,N,WM(3),N,IWM(31),IER)
!          CALL DGETRF(N,N,WM(1),N,IWM(31),IER)
!          IF (IER /= 0) IERPJ = 1
!        END IF
! ______________________________________________________________________

          return
        end if
!       End of code block for MITER = 1 or 2.

        if (MITER==3) then
!         If MITER = 3, construct a diagonal approximation to
!         J and P.
          NJE = NJE + 1
          JCUR = 1
          WM2 = HRL1
          R = RL1*PT1
          Y(1:N) = Y(1:N) + R*(H*SAVF(1:N)-YH(1:N,2))
!         CALL F(N,TN,Y,WM(3))
          call F(N,TN,Y,WM(1))
          NFE = NFE + 1
          do 10 I = 1, N
            R0 = H*SAVF(I) - YH(I,2)
!           DI = PT1*R0 - H*(WM(I+2)-SAVF(I))
            DI = PT1*R0 - H*(WM(I)-SAVF(I))
!           WM(I+2) = ONE
            WM(I) = ONE
            if (abs(R0)<UROUND/EWT(I)) goto 10
            if (abs(DI)<=ZERO) goto 20
!           WM(I+2) = PT1*R0/DI
            WM(I) = PT1*R0/DI
10        end do
          return
20        IERPJ = 1
          return
        end if
!       End of code block for MITER = 3.

!       Set constants for MITER = 4 or 5.
        ML = IWM(1)
        MU = IWM(2)
!       ML3 = ML + 3
        ML1 = ML + 1
        MBAND = ML + MU + 1
        MEBAND = MBAND + ML
        LENP = MEBAND*N

        if (JOK==-1 .and. MITER==4) then
!       If JOK = -1 and MITER = 4, call JAC to evaluate Jacobian.
          NJE = NJE + 1
          NSLJ = NST
          JCUR = 1
!         WM(3:LENP+2) = ZERO
          WM(1:LENP) = ZERO
!         CALL JAC(N,TN,Y,ML,MU,WM(ML3),MEBAND)
          call JAC(N,TN,Y,ML,MU,WM(ML1),MEBAND)
          if (J_IS_CONSTANT) J_HAS_BEEN_COMPUTED = .true.
!         IF (JSV==1) CALL DACOPY(MBAND,N,WM(ML3),MEBAND,WM(LOCJS),MBAND)
          if (JSV==1) call DACOPY(MBAND,N,WM(ML1),MEBAND,WM(LOCJS),MBAND)
        end if

        if (JOK==-1 .and. MITER==5) then
!       If MITER = 5, make ML+MU+1 calls to F to approximate the Jacobian
!       unless the user has specified which sub and super diagonals are
!       nonzero. In the latter case NSUBS+NSUPS+1 calls will be made.

!         If the subdiagonals are known, use that information to build
!         the Jacobian matrix.
          if ((SUBS .or. SUPS) .or. BNGRP >= MBAND) goto 60

!         Otherwise, use the original algorithm.
          if (USE_JACSP) then
!            Approximate the Jacobian using Doug Salane's JACSP.
             NSLJ = NST
             JCUR = 1
             WM(1:LENP) = ZERO
             IOPTDS(1) = 1
             IOPTDS(2) = MBAND
             IOPTDS(3) = 1
             IOPTDS(5) = ML
!            INFORDS(4) was initialized in DVODE.
             LWKDS  = 3 * N
             LIWKDS = 50 + N
!            NRFJACDS = MEBAND*N
!            NCFJACDS = 1
             NRFJACDS = MEBAND
             NCFJACDS = N
             MBA = min(MBAND,N)
             MAXGRPDS = MBA
!            Calculate the YSCALEDS vector for JACSPDV.
             if (LIKE_ORIGINAL_VODE) then
                FAC = DVNORM(N,SAVF,EWT)
!               JACSPDB multiplies YSCALEDS(*) BY UROUND**0.825:
!               R0 = THOU*ABS(H)*UROUND*REAL(N)*FAC
                R0 = THOU*abs(H)*real(N)*FAC
                if (abs(R0)<=ZERO) R0 = ONE
                SRUR = WM1
                do J = 1, N
!                  JACSPDB multiplies YSCALEDS(*) BY UROUND**0.825:
!                  R = MAX(ABS(Y(J)),R0/EWT(J))
                   R = max(abs(Y(J))/U325,(R0/EWT(J))*U125)
                   YSCALEDS(J) = R
                end do
             else
                if (ITOL == 1 .or. ITOL == 3) then
                   do J = 1, N
                      YSCALEDS(J) = max(abs(Y(J)),ATOL(1),UROUND)
                   end do
                else
                   do J = 1, N
                      YSCALEDS(J) = max(abs(Y(J)),ATOL(J),UROUND)
                   end do
                end if
             end if
             call JACSPDB(F,N,TN,Y,SAVF,WM(1),NRFJACDS, &
               YSCALEDS,FACDS,IOPTDS,WKDS,LWKDS,IWKDS,LIWKDS,  &
               MAXGRPDS,NGRPDS,JPNTRDS,INDROWDS)
             NFE = NFE + IWKDS(7)
             NJE = NJE + 1
          else
             NSLJ = NST
             JCUR = 1
             MBA = min(MBAND,N)
             MEB1 = MEBAND - 1
             SRUR = WM1
             FAC = DVNORM(N,SAVF,EWT)
             R0 = THOU*abs(H)*UROUND*real(N)*FAC
             if (abs(R0)<=ZERO) R0 = ONE
             do J = 1, MBA
               do I = J, N, MBAND
                 YI = Y(I)
                 R = max(SRUR*abs(YI),R0/EWT(I))
                 Y(I) = Y(I) + R
               end do
               call F(N,TN,Y,FTEM)
               do JJ = J, N, MBAND
                 Y(JJ) = YH(JJ,1)
                 YJJ = Y(JJ)
                 R = max(SRUR*abs(YJJ),R0/EWT(JJ))
                 FAC = ONE/R
                 I1 = max(JJ-MU,1)
                 I2 = min(JJ+ML,N)
!                II = JJ*MEB1 - ML + 2
                 II = JJ*MEB1 - ML
                 do I = I1, I2
                   WM(II+I) = (FTEM(I)-SAVF(I))*FAC
                 end do
               end do
             end do
             NFE = NFE + MBA
             NJE = NJE + 1
          end if
          if (J_IS_CONSTANT) J_HAS_BEEN_COMPUTED = .true.
          goto 90
   60     continue

!         User supplied diagonals information is available.
!         WM(3:LENP+2) = ZERO
          WM(1:LENP) = ZERO
          NJE = NJE + 1
          NSLJ = NST
          JCUR = 1
          MBA = min(MBAND,N)
          MEB1 = MEBAND - 1
          SRUR = WM1
          FAC = DVNORM(N,SAVF,EWT)
          R0 = THOU*abs(H)*UROUND*real(N)*FAC
          if (abs(R0)<=ZERO) R0 = ONE
!         For each group of columns...
          do NG = 1, BNGRP
!            Find the first and last columns in the group.
             JJ1 = BIGP(NG)
             JJ2 = BIGP(NG+1) - 1
!            For each column in this group...
             do JJ = JJ1, JJ2
                J = BJGP(JJ)
                R = max(SRUR*abs(Y(J)),R0/EWT(J))
                Y(J) = Y(J) + R
             end do
             call F(N,TN,Y,FTEM)
!            For each column in this group...
             do JJ = JJ1, JJ2
                 J = BJGP(JJ)
                 Y(J) = YH(J,1)
                 R = max(SRUR*abs(Y(J)),R0/EWT(J))
                 FAC = ONE/R
                 if (BUILD_IAJA) then
!                   Use the IAB, JAB sparse structure arrays to
!                   determine the first and last nonzeros in
!                   column J.
                    K1 = IAB(J)
                    K2 = IAB(J+1) - 1
                 else
!                   Extract the positions of the first and
!                   last nonzeros in this column directly.
                    call BANDED_GET_BJNZ(N,ML,MU,J,IWM(31),I)
                    do K = 1, N
                       if (IWM(30+K) /= 0) then
                          K1 = K
                          goto 70
                       end if
                    end do
   70               continue
                    do K = N, 1, -1
                       if (IWM(30+K) /= 0) then
                          K2 = K
                          goto 80
                       end if
                    end do
   80               continue
                 end if
!                Load the nonzeros for column J in the banded matrix.
                 if (BUILD_IAJA) then
                    do K = K1, K2
                       I = JAB(K)
!                      II = J * MEB1 - ML + 2
                       II = J * MEB1 - ML
                       WM(II+I) = (FTEM(I)-SAVF(I))*FAC
                    end do
                 else
                    do K = K1, K2
                       I = IWM(30+K)
                       if (I /= 0) then
!                         II = J * MEB1 - ML + 2
                          II = J * MEB1 - ML
                          WM(II+I) = (FTEM(I)-SAVF(I))*FAC
                       end if
                    end do
                 end if
             end do
          end do
          NFE = NFE + BNGRP
          if (J_IS_CONSTANT) J_HAS_BEEN_COMPUTED = .true.
   90     continue
!         IF (JSV==1) CALL DACOPY(MBAND,N,WM(ML3),MEBAND,WM(LOCJS),MBAND)
          if (JSV==1) call DACOPY(MBAND,N,WM(ML1),MEBAND,WM(LOCJS),MBAND)
        end if

        if (JOK==1) then
          JCUR = 0
!         CALL DACOPY(MBAND,N,WM(LOCJS),MBAND,WM(ML3),MEBAND)
          call DACOPY(MBAND,N,WM(LOCJS),MBAND,WM(ML1),MEBAND)
        end if

!       Multiply Jacobian by scalar, add identity, and do LU
!       decomposition.
        CON = -HRL1
!       CALL DSCAL_F90(LENP,CON,WM(3),1)
        call DSCAL_F90(LENP,CON,WM(1),1)
!       II = MBAND + 2
        II = MBAND
        do I = 1, N
          WM(II) = WM(II) + ONE
          II = II + MEBAND
        end do
        NLU = NLU + 1
! ______________________________________________________________________

!       CALL DGBFA_F90(WM(3),MEBAND,N,ML,MU,IWM(31),IER)
        call DGBFA_F90(WM(1),MEBAND,N,ML,MU,IWM(31),IER)
        if (IER/=0) IERPJ = 1
! *****LAPACK build change point. Replace above with these statements.
!       IF (.NOT.USE_LAPACK) THEN
!!        CALL DGBFA_f90(WM(3),MEBAND,N,ML,MU,IWM(31),IER)
!         CALL DGBFA_f90(WM(1),MEBAND,N,ML,MU,IWM(31),IER)
!         IF (IER /= 0) IERPJ = 1
!       ELSE
!!        CALL DGBTRF(N,N,ML,MU,WM(3),MEBAND,IWM(31),IER)
!         CALL DGBTRF(N,N,ML,MU,WM(1),MEBAND,IWM(31),IER)
!         IF (IER /= 0) IERPJ = 1
!       END IF
! ______________________________________________________________________
       return
!      End of code block for MITER = 4 or 5.

      end subroutine DVJAC
!_______________________________________________________________________

      subroutine DACOPY(NROW,NCOL,A,NROWA,B,NROWB)
! ..
! Copy one array to another.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: NCOL, NROW, NROWA, NROWB
! ..
! .. Array Arguments ..
        real (WP), intent (IN) :: A(NROWA,NCOL)
        real (WP), intent (INOUT) :: B(NROWB,NCOL)
! ..
! .. Local Scalars ..
        integer :: IC
! ..
! .. FIRST EXECUTABLE STATEMENT DACOPY
! ..
        do IC = 1, NCOL
          call DCOPY_F90(NROW,A(1,IC),1,B(1,IC),1)
        end do
        return

      end subroutine DACOPY
!_______________________________________________________________________

      subroutine DVSOL(WM,IWM,X,IERSL)
! ..
! Manage the solution of the linear system arising from a chord
! iteration for dense and banded solutions.
! ..
! This routine manages the solution of the linear system arising from
! a chord iteration. It is called if MITER /= 0.
! If MITER is 1 or 2, it calls DGESL_F90 to accomplish this.
! If MITER = 3 it updates the coefficient H*RL1 in the diagonal
! matrix, and then computes the solution.
! If MITER is 4 or 5, it calls DGBSL_F90.
! Communication with DVSOL uses the following variables:
! WM    = Real work space containing the inverse diagonal matrix if
!         MITER = 3 and the LU decomposition of the matrix otherwise.
!         Storage of matrix elements starts at WM(1).
!         WM also contains the following matrix-related data:
!         WM1 = SQRT(UROUND) (not used here),
!         WM2 = HRL1, the previous value of H*RL1, used if MITER = 3.
! IWM   = Integer work space containing pivot information, starting at
!         IWM(31), if MITER is 1, 2, 4, or 5. IWM also contains band
!         parameters ML = IWM(1) and MU = IWM(2) if MITER is 4 or 5.
! X     = The right-hand side vector on input, and the solution vector
!         on output, of length N.
! IERSL = Output flag. IERSL = 0 if no trouble occurred.
!         IERSL = 1 if a singular matrix arose with MITER = 3.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (INOUT) :: IERSL
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: WM(*), X(*)
        integer, intent (INOUT) :: IWM(*)
! ..
! .. Local Scalars ..
        real (WP) :: DI, HRL1, PHRL1, R
        integer :: I, MEBAND, ML, MU
! ..
! .. Intrinsic Functions ..
        intrinsic ABS
! ..
! ______________________________________________________________________

! *****LAPACK build change point. Insert this statement.
!       INTEGER INFO
!       CHARACTER*1 TRANS
! ______________________________________________________________________

! .. FIRST EXECUTABLE STATEMENT DVSOL
! ..
        IERSL = 0
        goto (10,10,20,50,50) MITER

10      continue
! ______________________________________________________________________

!       CALL DGESL_F90(WM(3),N,N,IWM(31),X,0)
        call DGESL_F90(WM(1),N,N,IWM(31),X,0)
! *****LAPACK build change point. Replace above with these statements.
!       IF (.NOT.USE_LAPACK) THEN
!!         CALL DGESL_f90(WM(3),N,N,IWM(31),X,0)
!          CALL DGESL_f90(WM(1),N,N,IWM(31),X,0)
!       ELSE
!          TRANS = 'N'
!!         CALL DGETRS(TRANS,N,1,WM(3),N,IWM(31),X,N,INFO)
!          CALL DGETRS(TRANS,N,1,WM(1),N,IWM(31),X,N,INFO)
!          IF (INFO /= 0) THEN
!             WRITE(6,*) 'Stopping in DVSOL with INFO = ', INFO
!             STOP
!          END IF
!       END IF
! ______________________________________________________________________

        return

20      PHRL1 = WM2
        HRL1 = H*RL1
        WM2 = HRL1
        if (abs(HRL1-PHRL1)<=ZERO) goto 30
        R = HRL1/PHRL1
        do I = 1, N
!         DI = ONE - R*(ONE-ONE/WM(I+2))
          DI = ONE - R*(ONE-ONE/WM(I))
          if (abs(DI)<=ZERO) goto 40
!         WM(I+2) = ONE/DI
          WM(I) = ONE/DI
        end do

30      do I = 1, N
!         X(I) = WM(I+2)*X(I)
          X(I) = WM(I)*X(I)
        end do
        return
40      IERSL = 1
        return

50      ML = IWM(1)
        MU = IWM(2)
        MEBAND = 2*ML + MU + 1
! ______________________________________________________________________

!       CALL DGBSL_F90(WM(3),MEBAND,N,ML,MU,IWM(31),X,0)
        call DGBSL_F90(WM(1),MEBAND,N,ML,MU,IWM(31),X,0)
! *****LAPACK build change point. Replace above with these statements.
!       IF (.NOT.USE_LAPACK) THEN
!!        CALL DGBSL_F90(WM(3),MEBAND,N,ML,MU,IWM(31),X,0)
!         CALL DGBSL_F90(WM(1),MEBAND,N,ML,MU,IWM(31),X,0)
!       ELSE
!         TRANS = 'N'
!!        CALL DGBTRS(TRANS,N,ML,MU,1,WM(3),MEBAND,IWM(31),X,N,INFO)
!         CALL DGBTRS(TRANS,N,ML,MU,1,WM(1),MEBAND,IWM(31),X,N,INFO)
!         IF (INFO /= 0) THEN
!           WRITE(6,*) 'Stopping in DVSOL with INFO = ', INFO
!           STOP
!         END IF
!       END IF
! ______________________________________________________________________

      return

      end subroutine DVSOL
!_______________________________________________________________________

      subroutine DVSRCO(RSAV,ISAV,JOB)
! ..
! Save or restore (depending on JOB) the contents of the PRIVATE
! variable blocks, which are used internally by DVODE (not called
! by DVODE_F90).
! ..
! RSAV = real array of length 49 or more.
! ISAV = integer array of length 41 or more.
! JOB  = flag indicating to save or restore the PRIVATE variable
!        blocks:
!        JOB  = 1 if PRIVATE variables is to be saved
!                 (written to RSAV/ISAV).
!        JOB  = 2 if PRIVATE variables is to be restored
!                 (read from RSAV/ISAV).
!        A call with JOB = 2 presumes a prior call with JOB = 1.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: JOB
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: RSAV(*)
        integer, intent (INOUT) :: ISAV(*)
! ..
! .. FIRST EXECUTABLE STATEMENT DVSRCO
! ..
        if (JOB/=2) then
!         Save the contents of the PRIVATE blocks.
          RSAV(1) = ACNRM
          RSAV(2) = CCMXJ
          RSAV(3) = CONP
          RSAV(4) = CRATE
          RSAV(5) = DRC
          RSAV(6:18) = EL(1:13)
          RSAV(19) = ETA
          RSAV(20) = ETAMAX
          RSAV(21) = H
          RSAV(22) = HMIN
          RSAV(23) = HMXI
          RSAV(24) = HNEW
          RSAV(25) = HSCAL
          RSAV(26) = PRL1
          RSAV(27) = RC
          RSAV(28) = RL1
          RSAV(29:41) = TAU(1:13)
          RSAV(42:46) = TQ(1:5)
          RSAV(47) = TN
          RSAV(48) = UROUND
          RSAV(LENRV1+1) = HU
          ISAV(1) = ICF
          ISAV(2) = INIT
          ISAV(3) = IPUP
          ISAV(4) = JCUR
          ISAV(5) = JSTART
          ISAV(6) = JSV
          ISAV(7) = KFLAG
          ISAV(8) = KUTH
          ISAV(9) = L
          ISAV(10) = LMAX
          ISAV(11) = LYH
          ISAV(12) = 0
          ISAV(13) = 0
          ISAV(14) = 0
          ISAV(15) = LWM
          ISAV(16) = LIWM
          ISAV(17) = LOCJS
          ISAV(18) = MAXORD
          ISAV(19) = METH
          ISAV(20) = MITER
          ISAV(21) = MSBJ
          ISAV(22) = MXHNIL
          ISAV(23) = MXSTEP
          ISAV(24) = N
          ISAV(25) = NEWH
          ISAV(26) = NEWQ
          ISAV(27) = NHNIL
          ISAV(28) = NQ
          ISAV(29) = NQNYH
          ISAV(30) = NQWAIT
          ISAV(31) = NSLJ
          ISAV(32) = NSLP
          ISAV(33) = NYH
          ISAV(LENIV1+1) = NCFN
          ISAV(LENIV1+2) = NETF
          ISAV(LENIV1+3) = NFE
          ISAV(LENIV1+4) = NJE
          ISAV(LENIV1+5) = NLU
          ISAV(LENIV1+6) = NNI
          ISAV(LENIV1+7) = NQU
          ISAV(LENIV1+8) = NST
          return
        else
!         Replace the contents of the PRIVATE blocks.
          ACNRM = RSAV(1)
          CCMXJ = RSAV(2)
          CONP = RSAV(3)
          CRATE = RSAV(4)
          DRC = RSAV(5)
          EL(1:13) = RSAV(6:18)
          ETA = RSAV(19)
          ETAMAX = RSAV(20)
          H = RSAV(21)
          HMIN = RSAV(22)
          HMXI = RSAV(23)
          HNEW = RSAV(24)
          HSCAL = RSAV(25)
          PRL1 = RSAV(26)
          RC = RSAV(27)
          RL1 = RSAV(28)
          TAU(1:13) = RSAV(29:41)
          TQ(1:5) = RSAV(42:46)
          TN = RSAV(47)
          UROUND = RSAV(48)
          HU = RSAV(LENRV1+1)
          ICF = ISAV(1)
          INIT = ISAV(2)
          IPUP = ISAV(3)
          JCUR = ISAV(4)
          JSTART = ISAV(5)
          JSV = ISAV(6)
          KFLAG = ISAV(7)
          KUTH = ISAV(8)
          L = ISAV(9)
          LMAX = ISAV(10)
          LYH = ISAV(11)
          LWM = ISAV(15)
          LIWM = ISAV(16)
          LOCJS = ISAV(17)
          MAXORD = ISAV(18)
          METH = ISAV(19)
          MITER = ISAV(20)
          MSBJ = ISAV(21)
          MXHNIL = ISAV(22)
          MXSTEP = ISAV(23)
          N = ISAV(24)
          NEWH = ISAV(25)
          NEWQ = ISAV(26)
          NHNIL = ISAV(27)
          NQ = ISAV(28)
          NQNYH = ISAV(29)
          NQWAIT = ISAV(30)
          NSLJ = ISAV(31)
          NSLP = ISAV(32)
          NYH = ISAV(33)
          NCFN = ISAV(LENIV1+1)
          NETF = ISAV(LENIV1+2)
          NFE = ISAV(LENIV1+3)
          NJE = ISAV(LENIV1+4)
          NLU = ISAV(LENIV1+5)
          NNI = ISAV(LENIV1+6)
          NQU = ISAV(LENIV1+7)
          NST = ISAV(LENIV1+8)
          return
        end if

      end subroutine DVSRCO
!_______________________________________________________________________

      subroutine DEWSET(N,ITOL,RTOL,ATOL,YCUR,EWT)
! ..
! Set the error weight vector.
! ..
! This subroutine sets the error weight vector EWT according to
! EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i), i = 1,...,N,
! with the subscript on RTOL and/or ATOL possibly replaced by 1
! above, depending on the value of ITOL.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: ITOL, N
! ..
! .. Array Arguments ..
        real (WP), intent (IN) :: ATOL(*), RTOL(*), YCUR(N)
        real (WP), intent (OUT) :: EWT(N)
! ..
! .. Intrinsic Functions ..
        intrinsic ABS
! ..
! .. FIRST EXECUTABLE STATEMENT DEWSET
! ..
        goto (10,20,30,40) ITOL
10      continue
        EWT(1:N) = RTOL(1)*abs(YCUR(1:N)) + ATOL(1)
        return
20      continue
        EWT(1:N) = RTOL(1)*abs(YCUR(1:N)) + ATOL(1:N)
        return
30      continue
        EWT(1:N) = RTOL(1:N)*abs(YCUR(1:N)) + ATOL(1)
        return
40      continue
        EWT(1:N) = RTOL(1:N)*abs(YCUR(1:N)) + ATOL(1:N)
        return

      end subroutine DEWSET
!_______________________________________________________________________

      function DVNORM(N,V,W)
! ..
! Calculate weighted root-mean-square (rms) vector norm.
! ..
! This routine computes the weighted root-mean-square norm
! of the vector of length N contained in the array V, with
! weights contained in the array W of length N.
! DVNORM = SQRT((1/N) * SUM(V(i)*W(i))**2)
! ..
     implicit none
! ..
! .. Function Return Value ..
        real (WP) :: DVNORM
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: N
! ..
! .. Array Arguments ..
        real (WP), intent (IN) :: V(N), W(N)
! ..
! .. Local Scalars ..
        real (WP) :: SUM
        integer :: I
! ..
! .. Intrinsic Functions ..
        intrinsic SQRT
! ..
! .. FIRST EXECUTABLE STATEMENT DVNORM
! ..
        SUM = ZERO
        do I = 1, N
          SUM = SUM + (V(I)*W(I))**2
        end do
        DVNORM = sqrt(SUM/N)
        return

      end function DVNORM
!_______________________________________________________________________

! The modified SLATEC error handling routines begin here.
      subroutine XERRDV(MSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
! ..
! Write error messages with values.
! ..
! This is an adaptation of subroutine XERRWD (NMES eliminated).
! Subroutines XERRDV, XSETF, XSETUN, and Functions IXSAV
! as given here, constitute a simplified version of the SLATEC
! error handling package.
!
! All arguments are input arguments.
! MSG    = The message (character array).
! NERR   = The error number (not used).
! LEVEL  = The error level
!          0 or 1 means recoverable (control returns to caller).
!          2 means fatal (run is aborted--see note below).
! NI     = Number of integers (0, 1, or 2) to be printed with message.
! I1,I2  = Integers to be printed, depending on NI.
! NR     = Number of reals (0, 1, or 2) to be printed with message.
! R1,R2  = Reals to be printed, depending on NR.
! Note: This routine is machine-dependent and specialized for use
! in limited context, in the following ways:
!  1. The argument MSG is assumed to be of type CHARACTER, and
!     the message is printed with a format of (1X,A).
!  2. The message is assumed to take only one line.
!     Multi-line messages are generated by repeated calls.
!  3. If LEVEL = 2, control passes to the statement: STOP
!     to abort the run. This statement may be machine-dependent.
!  4. R1 and R2 are assumed to be in real(wp) and are printed
!     in D21.13 format.
! Internal Notes:
! For a different default logical unit number, IXSAV(or a subsidiary
! function that it calls) will need to be modified.
! For a different run-abort command, change the statement following
! statement 100 at the end.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP) :: R1, R2
        integer :: I1, I2, LEVEL, NERR, NI, NR
        character (*) :: MSG
        logical PRINT_NERR
! ..
! .. Local Scalars ..
        integer :: LUNIT, MESFLG
! ..
! .. FIRST EXECUTABLE STATEMENT XERRDV
! ..
!       Get logical unit number and message print flag.
        LUNIT = IXSAV(1,0,.false.)
        MESFLG = IXSAV(2,0,.false.)
        if (MESFLG==0) goto 10

        PRINT_NERR = .false.
        if (PRINT_NERR) print *, MSG, 'Message number = ', NERR

!       Write the message.
        write (LUNIT,90000) MSG
90000   format (1X,A)
        if (NI==1) write (LUNIT,90001) I1
90001   format ('In the above message, I1 = ',I10)
        if (NI==2) write (LUNIT,90002) I1, I2
90002   format ('In the above message, I1 = ',I10,3X,'I2 = ',I10)
        if (NR==1) write (LUNIT,90003) R1
90003   format ('In the above message, R1 = ',D21.13)
        if (NR==2) write (LUNIT,90004) R1, R2
90004   format ('In the above message, R1 = ',D21.13,3X,'R2 = ',D21.13)

!       Abort the run if LEVEL = 2.
10      if (LEVEL/=2) return
        write (LUNIT,90005)
90005   format ('LEVEL = 2 in XERRDV. Stopping.')
        stop

      end subroutine XERRDV
!_______________________________________________________________________

      subroutine XSETF(MFLAG)
! ..
!     Reset the error print control flag.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: MFLAG
! ..
! .. Local Scalars ..
        integer :: JUNK
! ..
! .. FIRST EXECUTABLE STATEMENT XSETF
! ..
        if (MFLAG==0 .or. MFLAG==1) JUNK = IXSAV(2,MFLAG,.true.)
!       Get rid of a compiler warning message:
        if (JUNK/=JUNK) stop
        return

      end subroutine XSETF
!_______________________________________________________________________

      subroutine XSETUN(LUN)
! ..
!     Reset the logical unit number for error messages.
! ..
!     XSETUN sets the logical unit number for error messages to LUN.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer :: LUN
! ..
! .. Local Scalars ..
        integer :: JUNK
! ..
! .. FIRST EXECUTABLE STATEMENT XSETUN
! ..
        if (LUN>0) JUNK = IXSAV(1,LUN,.true.)
!       Get rid of a compiler warning message:
        if (JUNK/=JUNK) stop
        return

      end subroutine XSETUN
!_______________________________________________________________________

      function IXSAV(IPAR,IVALUE,ISET)
! ..
! Save and recall error message control parameters.
! ..
! IXSAV saves and recalls one of two error message parameters:
!  LUNIT, the logical unit number to which messages are printed,
!  and MESFLG, the message print flag.
! This is a modification of the SLATEC library routine J4SAVE.
! Saved local variables:
!  LUNIT  = Logical unit number for messages. The default is
!           obtained by a call to IUMACH(may be machine-dependent).
!  MESFLG = Print control flag:
!           1 means print all messages (the default).
!           0 means no printing.
! On input:
!  IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
!  IVALUE = The value to be set for the parameter, if ISET = .TRUE.
!  ISET   = Logical flag to indicate whether to read or write.
!           If ISET = .TRUE., the parameter will be given
!           the value IVALUE. If ISET = .FALSE., the parameter
!           will be unchanged, and IVALUE is a dummy argument.
! On return:
!   IXSAV = The (old) value of the parameter.
! ..
     implicit none
! ..
! .. Function Return Value ..
        integer :: IXSAV
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: IPAR, IVALUE
        logical, intent (IN) :: ISET
! ..
! .. Local Scalars ..
        integer, save :: LUNIT, MESFLG
! ..
! .. Data Statements ..
        data LUNIT/ -1/, MESFLG/1/
! ..
! .. FIRST EXECUTABLE STATEMENT IXSAV
! ..
        if (IPAR==1) then
          if (LUNIT==-1) LUNIT = IUMACH()
!  i       Get rid of a compiler warning message:
          if (LUNIT/=LUNIT) stop
          IXSAV = LUNIT
          if (ISET) LUNIT = IVALUE
        end if

        if (IPAR==2) then
          IXSAV = MESFLG
          if (ISET) MESFLG = IVALUE
        end if
        return

      end function IXSAV
!_______________________________________________________________________

      function IUMACH()
! ..
!     Provide the standard output unit number.
! ..
!     INTEGER LOUT, IUMACH
!     LOUT = IUMACH()
!     Function Return Values:
!     LOUT: the standard logical unit for Fortran output.
!     Internal Notes:
!     The built-in value of 6 is standard on a wide range
!     of Fortran systems. This may be machine-dependent.
! ..
! .. Function Return Value ..
        integer :: IUMACH
! ..
! .. FIRST EXECUTABLE STATEMENT IUMACH
! ..
        IUMACH = 6
        return

      end function IUMACH
! The modified error handling routines end here.
!_______________________________________________________________________

      subroutine CHECK_STAT(IER,CALLED_FROM_WHERE)
! ..
! Print an error message if a storage allocation error
! occurred.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: CALLED_FROM_WHERE, IER
! ..
! .. Local Scalars ..
        integer :: I1
        character (80) :: MSG
! ..
! .. FIRST EXECUTABLE STATEMENT CHECK_STAT
! ..
        if (IER/=0) then
          I1 = CALLED_FROM_WHERE
          MSG = 'A storage allocation error occurred.'
          call XERRDV(MSG,1410,1,0,0,0,0,ZERO,ZERO)
          MSG = 'The error occurred at location I1.'
          call XERRDV(MSG,1410,2,1,I1,0,0,ZERO,ZERO)
        end if
        return

      end subroutine CHECK_STAT
!_______________________________________________________________________

      subroutine DVPREPS(NEQ,Y,YH,LDYH,SAVF,EWT,F,JAC)
! ..
! Determine the sparsity structure and allocate the necessary arrays
! for MA28 based sparse solutions.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: LDYH, NEQ
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: EWT(*), Y(*)
        real (WP) :: SAVF(*)
        real (WP), intent (IN) :: YH(LDYH,*)
! ..
! .. Subroutine Arguments ..
        external F, JAC
! ..
! .. Local Scalars ..
        real (WP) :: DQ, DYJ, ERWT, FAC, YJ
        integer :: I, IER, J, JFOUND, K, KMAX, KMIN, KNEW, NP1, NZ
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, ALLOCATED, MAX, SIGN
! ..
! .. FIRST EXECUTABLE STATEMENT DVPREPS
! ..
        NZ_SWAG = max(max(1000,NZ_SWAG),10*N)
        NP1 = N + 1
        NZ_ALL = NZ_SWAG
!       ADDTONNZ = MAX(1000,NZ_SWAG)
        ADDTONNZ = NZ_SWAG
10      continue
        if (allocated(IAN)) then
          deallocate (IAN,JAN,IGP,JGP,FTEMP1,IKEEP28,IW28,ICN,PMAT, &
            JVECT,STAT=IER)
          call CHECK_STAT(IER,490)
          if (allocated(JMAT)) then
            deallocate (JMAT,STAT=IER)
            call CHECK_STAT(IER,500)
          end if
        end if
        NZ_ALL = NZ_ALL + ADDTONNZ
        LICN_ALL = ELBOW_ROOM * NZ_ALL
        LIRN_ALL = ELBOW_ROOM * NZ_ALL
        if (LICN_ALL>MAX_ARRAY_SIZE .or. LIRN_ALL>MAX_ARRAY_SIZE) then
          MSG = 'Maximum array size exceeded. Stopping.'
          call XERRDV(MSG,1420,2,0,0,0,0,ZERO,ZERO)
        end if
!       Note: ICN may need to be reallocated in DVJACS28.
        allocate (IAN(NP1),JAN(LIRN_ALL),IGP(NP1),JGP(N),FTEMP1(N), &
          IKEEP28(N,5),IW28(N,8),ICN(LICN_ALL),PMAT(LICN_ALL),      &
          JVECT(LIRN_ALL),STAT=IER)
        call CHECK_STAT(IER,510)
        if (JSV==1) then
          allocate (JMAT(NZ_ALL),STAT=IER)
          call CHECK_STAT(IER,520)
          JMAT(1:NZ_ALL) = ZERO
        end if

        if (MOSS==0) goto 30
        if (ISTATC==3) goto 20

!       ISTATE = 1 and MOSS /= 0.

!       Perturb Y for structure determination:
        do I = 1, N
          ERWT = ONE/EWT(I)
          FAC = ONE + ONE/(I+ONE)
          Y(I) = Y(I) + FAC*sign(ERWT,Y(I))
        end do
        goto (60,70) MOSS

20      continue
!       ISTATE = 3 and MOSS /= 0.

!       Load Y from YH(*,1):
        Y(1:N) = YH(1:N,1)
        goto (60,70) MOSS

!       MOSS = 0

!       Process user's IA,JA. Add diagonal entries if necessary:
30      continue
        if (IAJA_CALLED) then
        else
          MSG = 'You have indicated that you wish to supply the'
          call XERRDV(MSG,1430,1,0,0,0,0,ZERO,ZERO)
          MSG = 'sparsity arrays IA and JA directly but you did'
          call XERRDV(MSG,1430,1,0,0,0,0,ZERO,ZERO)
          MSG = 'not call SET_IAJA after calling SET_OPTS.'
          call XERRDV(MSG,1430,2,0,0,0,0,ZERO,ZERO)
        end if
        KNEW = 1
        KMIN = IA(1)
        IAN(1) = 1
        do J = 1, N
          JFOUND = 0
          KMAX = IA(J+1) - 1
          if (KMIN>KMAX) goto 40
          do K = KMIN, KMAX
            I = JA(K)
            if (I==J) JFOUND = 1
            if (KNEW>NZ_ALL) then
               if (LP /= 0) then
                  MSG = 'NZ_ALL (=I1) is not large enough.'
                  call XERRDV(MSG,1440,1,0,0,0,0,ZERO,ZERO)
                  MSG = 'Allocating more space for another try.'
                  call XERRDV(MSG,1440,1,1,NZ_ALL,0,0,ZERO,ZERO)
               end if
               goto 10
            end if
            JAN(KNEW) = I
            KNEW = KNEW + 1
          end do
          if (JFOUND==1) goto 50
40        if (KNEW>NZ_ALL) then
             if (LP /= 0) then
                MSG = 'NZ_ALL (=I1) is not large enough.'
                call XERRDV(MSG,1450,1,0,0,0,0,ZERO,ZERO)
                MSG = 'Allocating more space for another try.'
                call XERRDV(MSG,1450,1,1,NZ_ALL,0,0,ZERO,ZERO)
             end if
             goto 10
          end if
          JAN(KNEW) = J
          KNEW = KNEW + 1
50        IAN(J+1) = KNEW
          KMIN = KMAX + 1
        end do
        goto 90

60      continue

!       MOSS = 1.

!       Compute structure from user-supplied Jacobian routine JAC.
        NZ = 0
        call JAC(NEQ,TN,Y,IAN,JAN,NZ,PMAT)
        if (NZ<=0) then
          MSG = 'Illegal value of NZ from JAC in DVPREPS.'
          call XERRDV(MSG,1460,2,0,0,0,0,ZERO,ZERO)
        end if
        if (NZ>NZ_ALL) then
           if (LP /= 0) then
              MSG = 'NZ_ALL (=I1) is not large enough.'
              call XERRDV(MSG,1470,1,0,0,0,0,ZERO,ZERO)
              MSG = 'Allocating more space for another try.'
              call XERRDV(MSG,1470,1,1,NZ_ALL,0,0,ZERO,ZERO)
           end if
           goto 10
        end if
        call JAC(NEQ,TN,Y,IAN,JAN,NZ,PMAT)
        call SET_ICN(N,IAN,ICN)
        call CHECK_DIAG(N,IAN,JAN,ICN)
        goto 90

!       MOSS = 2.

!       Compute structure from results of N+1 calls to F.
70      K = 1
        IAN(1) = 1
        do I = 1, N
          ERWT = ONE/EWT(I)
          FAC = ONE + ONE/(I+ONE)
          Y(I) = Y(I) + FAC*sign(ERWT,Y(I))
        end do
        call F(NEQ,TN,Y,SAVF)
        NFE = NFE + 1
        do J = 1, N
          if (K>NZ_ALL) then
             if (LP /= 0) then
                MSG = 'NZ_ALL (=I1) is not large enough.'
                call XERRDV(MSG,1480,1,0,0,0,0,ZERO,ZERO)
                MSG = 'Allocating more space for another try.'
                call XERRDV(MSG,1480,1,1,NZ_ALL,0,0,ZERO,ZERO)
             end if
             goto 10
          end if
          YJ = Y(J)
          ERWT = ONE/EWT(J)
          DYJ = sign(ERWT,YJ)
          Y(J) = YJ + DYJ
          call F(NEQ,TN,Y,FTEMP1)
          NFE = NFE + 1
          Y(J) = YJ
          do 80 I = 1, N
            DQ = (FTEMP1(I)-SAVF(I))/DYJ
            if ((abs(DQ)<=SETH) .and. (I/=J)) goto 80
            JAN(K) = I
            K = K + 1
80        end do
          IAN(J+1) = K
        end do
90      continue
        if (MOSS==0 .or. ISTATC/=1) goto 100
!       If ISTATE = 1 and MOSS /= 0, restore Y from YH.
        Y(1:N) = YH(1:N,1)
100     NNZ = IAN(NP1) - 1
        LENIGP = 0
        MAXG = 0
        if (MITER==7) then
!         Compute grouping of column indices.
          MAXG = NP1
          call DGROUP(N,IAN,JAN,MAXG,NGP,IGP,JGP,IKEEP28(1,1), &
            IKEEP28(1,2),IER)
          if (IER/=0) then
            MSG = 'An error occurred in DGROUP.'
            call XERRDV(MSG,1490,2,0,0,0,0,ZERO,ZERO)
          end if
          LENIGP = NGP + 1
        end if

        if (USE_JACSP .and. MITER==7) then
!         Use Doug Salane's Jacobian routines to determine the column
!         grouping; and allocate and initialize the necessary arrays
!         for use in DVJACS48.
          if (allocated(INDROWDS)) then
             deallocate (INDROWDS, INDCOLDS, NGRPDS, IPNTRDS, JPNTRDS, IWADS, &
               IWKDS, IOPTDS, YSCALEDS, WKDS, FACDS)
             call CHECK_STAT(IER,530)
          end if
!         We could delete IWADS and use IW28 array.
          allocate (INDROWDS(NNZ), INDCOLDS(NNZ), NGRPDS(N+1), IPNTRDS(N+1), &
            JPNTRDS(N+1), IWADS(6*N), IWKDS(50+N), IOPTDS(5), YSCALEDS(N),  &
            WKDS(3*N), FACDS(N) ,STAT=IER)
          call CHECK_STAT(IER,540)
          INDROWDS(1:NNZ) = JAN(1:NNZ)
          call SET_ICN(N,IAN,INDCOLDS)
          call CHECK_DIAG(N,IAN,INDROWDS,INDCOLDS)
          LIWADS = 6 * N
          call DVDSM(N,N,NNZ,INDROWDS,INDCOLDS,NGRPDS,MAXGRPDS,MINGRPDS,      &
            INFODS,IPNTRDS,JPNTRDS,IWADS,LIWADS)
          if (INFODS /= 1) then
             MSG = 'An error occurred in subroutine DSM. INFO = I1.'
               call XERRDV(MSG,1500,2,1,INFODS,0,0,ZERO,ZERO)
          end if
!         For use in DVJACS28:
          IOPTDS(4) = 0
!         Define the IGP and JGP arrays needed by DVJACS28.
          call DGROUPDS(N,MAXGRPDS,NGRPDS,IGP,JGP)
          NGP = MAXGRPDS
          LENIGP = MAXGRPDS + 1
        end if

!       Trim the arrays to the final sizes.

        if (NZ_ALL>NNZ) then
          NZ_ALL = NNZ
          MAX_NNZ = max(MAX_NNZ,NNZ)
          LIRN_ALL = ELBOW_ROOM * NNZ
          ICN(1:NNZ) = JAN(1:NNZ)
          deallocate (JAN,STAT=IER)
          call CHECK_STAT(IER,550)
          allocate (JAN(LIRN_ALL),STAT=IER)
          call CHECK_STAT(IER,560)
          JAN(1:NNZ) = ICN(1:NNZ)
          call CHECK_STAT(IER,570)
          deallocate (ICN,PMAT,JVECT,STAT=IER)
          call CHECK_STAT(IER,580)
          if (allocated(JMAT)) then
            deallocate (JMAT,STAT=IER)
            call CHECK_STAT(IER,590)
          end if
          LICN_ALL = ELBOW_ROOM * NNZ
          allocate (ICN(LICN_ALL),PMAT(LICN_ALL),JVECT(LIRN_ALL),STAT=IER)
          call CHECK_STAT(IER,600)
          if (JSV==1) then
            allocate (JMAT(NNZ),STAT=IER)
            call CHECK_STAT(IER,610)
            JMAT(1:NNZ) = ZERO
          end if
          if (MITER==7) then
            JVECT(1:LENIGP) = IGP(1:LENIGP)
            deallocate (IGP,STAT=IER)
            call CHECK_STAT(IER,620)
            allocate (IGP(LENIGP),STAT=IER)
            call CHECK_STAT(IER,630)
            IGP(1:LENIGP) = JVECT(1:LENIGP)
          end if
        end if
        if (SCALE_MATRIX) then
           if (allocated(CSCALEX)) then
              deallocate (CSCALEX,RSCALEX,WSCALEX,STAT=IER)
              call CHECK_STAT(IER,640)
              allocate (CSCALEX(N),RSCALEX(N),WSCALEX(N,5),STAT=IER)
              call CHECK_STAT(IER,650)
           else
              allocate (CSCALEX(N),RSCALEX(N),WSCALEX(N,5),STAT=IER)
              call CHECK_STAT(IER,660)
           end if
        end if
        if (LP /= 0) then
           MSG = 'The final DVPRPEPS storage allocations are:'
           call XERRDV(MSG,1510,1,0,0,0,0,ZERO,ZERO)
           MSG = '   NZ_ALL (=I1):'
           call XERRDV(MSG,1510,1,1,NZ_ALL,0,0,ZERO,ZERO)
           MSG = '   LIRN_ALL (=I1) and LICN_ALL (=I2):'
           call XERRDV(MSG,1510,1,2,LIRN_ALL,LICN_ALL,0,ZERO,ZERO)
         end if
        return

      end subroutine DVPREPS
!_______________________________________________________________________

      subroutine DVRENEW(NEQ,Y,SAVF,EWT,F)
! ..
! In the event MA28BD encounters a zero pivot in the LU factorization
! of the iteration matrix due to an out-of-date MA28AD pivot sequence,
! re-calculate the sparsity structure using finite differences.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: NEQ
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: EWT(*), SAVF(*), Y(*)
! ..
! .. Subroutine Arguments ..
        external F
! ..
! .. Local Scalars ..
        real (WP) :: DQ, DYJ, ERWT, FAC, YJ
        integer :: ADDTONZ, I, IER, J, K, KVAL, NP1
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, ALLOCATED, MAX, SIGN
! ..
! .. FIRST EXECUTABLE STATEMENT DVRENEW
! ..
! .. Caution:
!    This routine must not be called before DVPREPS has been called.
!
!    Note:
!    On entry to DVRENEW, the allocated array sizes of arrays that
!    may change size are:
!      ICN, PMAT  = length LICN_ALL
!      JAN, JVECT = length LIRN_ALL
!      JMAT       = length NZ_ALL
!      IGP        = length LENIGP on first entry; N+1 thereafter

!       Check if a numerical Jacobian is being used and stop if
!       it is not.
        if (MITER /= 7) then
          MSG = 'DVRENEW can be used only if MITER = 7.'
          call XERRDV(MSG,1520,2,0,0,0,0,ZERO,ZERO)
        end if

!       Save Y and SAVF.
        if (.not.allocated(YTEMP)) then
           allocate (YTEMP(N),DTEMP(N),STAT=IER)
           call CHECK_STAT(IER,670)
        end if
        YTEMP(1:N) = Y(1:N)
        DTEMP(1:N) = SAVF(1:N)

!       Define the amount to be added to the array lengths
!       if necessary.
        NP1 = N + 1
        NNZ = IAN(NP1) - 1
        ADDTONZ = ELBOW_ROOM * NNZ

!       Just change the size of IGP to N+1 if have not already
!       done so.
        if (size(IGP) /= NP1) then
           deallocate (IAN,STAT=IER)
           call CHECK_STAT(IER,680)
           allocate (IAN(NP1),STAT=IER)
           call CHECK_STAT(IER,690)
        end if

!       Go to the differencing section to determine the new
!       sparsity structure.
        goto 20

10      continue

!       Reallocate the arrays if necessary.
        if (KVAL > LIRN_ALL) then
!          Note: JAN and JVECT may need to be reallocated in DVJACS28.
           LIRN_ALL = LIRN_ALL + ADDTONZ
           if (LIRN_ALL>MAX_ARRAY_SIZE) then
              MSG = 'Maximum array size exceeded. Stopping in DVRENEW.'
             call XERRDV(MSG,1530,2,0,0,0,0,ZERO,ZERO)
           end if
           deallocate (JAN,JVECT,STAT=IER)
           call CHECK_STAT(IER,700)
           allocate (JAN(LIRN_ALL),JVECT(LIRN_ALL),STAT=IER)
           call CHECK_STAT(IER,710)
        end if
        if (KVAL > LICN_ALL) then
!          Note: ICN and PMAT may need to be reallocated in DVJACS28.
           LICN_ALL = LICN_ALL + ADDTONZ
           if (LICN_ALL>MAX_ARRAY_SIZE) then
              MSG = 'Maximum array size exceeded. Stopping in DVRENEW.'
             call XERRDV(MSG,1540,2,0,0,0,0,ZERO,ZERO)
           end if
           deallocate (ICN,PMAT,STAT=IER)
           call CHECK_STAT(IER,720)
           allocate (ICN(LICN_ALL),PMAT(LICN_ALL),STAT=IER)
           call CHECK_STAT(IER,730)
        end if

   20   continue

!       Perturb Y for structure determination:
        do I = 1, N
          ERWT = ONE/EWT(I)
          FAC = ONE + ONE/(I+ONE)
          Y(I) = Y(I) + FAC*sign(ERWT,Y(I))
        end do

!       Compute structure from results of N+1 calls to F.
        K = 1
        IAN(1) = 1
        do I = 1, N
          ERWT = ONE/EWT(I)
          FAC = ONE + ONE/(I+ONE)
          Y(I) = Y(I) + FAC*sign(ERWT,Y(I))
        end do
        call F(NEQ,TN,Y,SAVF)
        NFE = NFE + 1
        do J = 1, N
          KVAL = K
          if (KVAL > LIRN_ALL) then
             if (LP /= 0) then
                MSG = 'LIRN_ALL (=I1) is not large enough.'
                call XERRDV(MSG,1550,1,0,0,0,0,ZERO,ZERO)
                MSG = 'Allocating more space for another try.'
                call XERRDV(MSG,1550,1,1,LIRN_ALL,0,0,ZERO,ZERO)
             end if
             goto 10
          end if
          if (KVAL > LICN_ALL) then
             if (LP /= 0) then
                MSG = 'LICN_ALL (=I1) is not large enough.'
                call XERRDV(MSG,1560,1,0,0,0,0,ZERO,ZERO)
                MSG = 'Allocating more space for another try.'
                call XERRDV(MSG,1560,1,1,LICN_ALL,0,0,ZERO,ZERO)
             end if
             goto 10
          end if
          YJ = Y(J)
          ERWT = ONE/EWT(J)
          DYJ = sign(ERWT,YJ)
          Y(J) = YJ + DYJ
          call F(NEQ,TN,Y,FTEMP1)
          NFE = NFE + 1
          Y(J) = YJ
          do 80 I = 1, N
            DQ = (FTEMP1(I)-SAVF(I))/DYJ
            if ((abs(DQ)<=SETH) .and. (I/=J)) goto 80
            JAN(K) = I
            K = K + 1
80        end do
          IAN(J+1) = K
        end do

        NNZ = IAN(NP1) - 1
        if (NNZ > NZ_ALL .and. JSV == 1) then
!          Increase the size of JMAT if necessary.
           NZ_ALL = NNZ
           if (LP /= 0) then
              MSG = 'NZ_ALL (=I1) is not large enough.'
              call XERRDV(MSG,1570,1,0,0,0,0,ZERO,ZERO)
              MSG = 'Allocating more space for another try.'
              call XERRDV(MSG,1570,1,1,NZ_ALL,0,0,ZERO,ZERO)
           end if
           if (NZ_ALL>MAX_ARRAY_SIZE) then
              MSG = 'Maximum array size exceeded. Stopping in DVRENEW.'
              call XERRDV(MSG,1580,2,0,0,0,0,ZERO,ZERO)
           end if
           deallocate (JMAT,STAT=IER)
           call CHECK_STAT(IER,740)
           allocate (JMAT(NZ_ALL),STAT=IER)
           call CHECK_STAT(IER,750)
           JMAT(1:NZ_ALL) = ZERO
        end if

!       Compute grouping of column indices.

        LENIGP = 0
        MAXG = 0
        MAXG = NP1
        call DGROUP(N,IAN,JAN,MAXG,NGP,IGP,JGP,IKEEP28(1,1), &
          IKEEP28(1,2),IER)
        if (IER/=0) then
          MSG = 'An error occurred in DGROUP.'
          call XERRDV(MSG,1590,2,0,0,0,0,ZERO,ZERO)
        end if
        LENIGP = NGP + 1

!       Restore Y and SAVF.
        Y(1:N) = YTEMP(1:N)
        SAVF(1:N) = DTEMP(1:N)
        return

      end subroutine DVRENEW
!_______________________________________________________________________

      subroutine DGROUP(N,IA,JA,MAXG,NGRP,IGP,JGP,INCL,JDONE,IER)
! ..
! Construct groupings of the column indices of the Jacobian matrix,
! used in the numerical evaluation of the Jacobian by finite
! differences for sparse solutions.
! ..
!     Input:
!     N      = the order of the matrix
!     IA,JA  = sparse structure descriptors of the matrix by rows
!     MAXG   = length of available storage in the IGP array
!     INCL and JDONE are working arrays of length N.
!     Output:
!     NGRP   = number of groups
!     JGP    = array of length N containing the column indices by
!              groups
!     IGP    = pointer array of length NGRP + 1 to the locations
!              in JGP of the beginning of each group
!     IER    = error indicator. IER = 0 if no error occurred, or
!              1 if MAXG was insufficient
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (INOUT) :: IER, NGRP
        integer, intent (IN) :: MAXG, N
! ..
! .. Array Arguments ..
        integer, intent (IN) :: IA(*), JA(*)
        integer, intent (INOUT) :: IGP(*), INCL(*), JDONE(*), JGP(*)
! ..
! .. Local Scalars ..
        integer :: I, J, K, KMAX, KMIN, NCOL, NG
! ..
! .. FIRST EXECUTABLE STATEMENT DGROUP
! ..
        IER = 0
        JDONE(1:N) = 0
        NCOL = 1
        do NG = 1, MAXG
          IGP(NG) = NCOL
          INCL(1:N) = 0
          do 20 J = 1, N
!           Reject column J if it is already in a group.
            if (JDONE(J)==1) goto 20
            KMIN = IA(J)
            KMAX = IA(J+1) - 1
            do 10 K = KMIN, KMAX
!           Reject column J if it overlaps any column already
!           in this group.
              I = JA(K)
              if (INCL(I)==1) goto 20
10          end do
!           Accept column J into group NG.
            JGP(NCOL) = J
            NCOL = NCOL + 1
            JDONE(J) = 1
            do K = KMIN, KMAX
              I = JA(K)
              INCL(I) = 1
            end do
20        end do
!         Stop if this group is empty (grouping is complete).
          if (NCOL==IGP(NG)) goto 30
        end do
!       Error return if not all columns were chosen (MAXG too small).
        if (NCOL<=N) goto 40
        NG = MAXG
30      NGRP = NG - 1
        return
40      IER = 1
        return

      end subroutine DGROUP
!_______________________________________________________________________

      subroutine BGROUP(N,BJA,BINCL,BDONE,ML,MU)
! ..
! Construct groupings of the column indices of the Jacobian
! matrix, used in the numerical evaluation of the Jacobian
! by finite differences for banded solutions when the nonzero
! sub and super diagonals are known. BGROUP is similar to
! DGROUP but it does not require the sparse structure arrays
! and it uses real rather than integer work arrays.
! ..
!     Input:
!
!       N       = the order of the matrix (number of odes)
!       BINCL   = real working array of length N
!       BJA     = real working array of length N
!       BDONE   = real working array of length N
!       ML      = integer lower bandwidth
!       MU      = integer upper bandwidth
!
!     Output: (PRIVATE information used in DVJAC)
!
!       BNGRP   = integer number of groups
!       BJGP    = integer array of length N containing the
!                 column indices by groups
!       BIGP    = integer pointer array of length BNGRP + 1
!                 to the locations in BJGP of the beginning
!                 of each group
!
!       Note:
!          On output:
!             For I = 1, ..., BNGRP:
!                Start of group I:
!                  J = BIGP(I)
!                Number of columns in group I:
!                  K = BIGP(I+1) - BIGP(I)
!                Columns in group I:
!                  BJGP(J-1+L), L=1, ..., K
!
!    Note:
!    The three arrays BJA, BINCL, and BDONE are REAL to avoid the
!    necessity to allocate three new INTEGER arrays. BGROUP can
!    be called only at an integration start or restart because
!    DVODE_F90 work arrays are used for these arrays.
!
!    Note:
!    The PRIVATE banded information SUBDS, NSUBDS, SUPDS, NSUPS,
!    ML, and MU must be defined before calling BGROUP.
!       ML        = lower bandwidth
!       MU        = upper bandwidth
!       NSUBS     = number of strict sub diagonals
!       SUBDS(I)  = row in which the Ith sub diagonal
!                   begins, I=1, ..., NSUBS
!       NSUPS     = number of strict super diagonals
!       SUPDS(I)  = column in which the Ith super diagonal
!                   begins, I=1, ..., NSUPS
! ..
        implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: N, ML, MU
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: BINCL(*), BJA(*), BDONE(*)
! ..
! .. Local Scalars ..
        integer :: I, IBDONE, IBINCL, IER, J, K, KBEGIN, KFINI, &
          KI, KJ, MAXG, NCOL, NG
        real (WP) :: FUDGE, ONE_PLUS_FUDGE
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
     intrinsic INT, MAX, MIN, real
! ..
! .. FIRST EXECUTABLE STATEMENT BGROUP
! ..
!       Storage for the column grouping information.
        if (allocated(BIGP)) then
           deallocate (BIGP,BJGP,STAT=IER)
           call CHECK_STAT(IER,760)
        end if
        allocate (BIGP(N+1),BJGP(N),STAT=IER)
        call CHECK_STAT(IER,770)

        FUDGE = 0.4_WP
        ONE_PLUS_FUDGE = 1.0_WP + FUDGE
        MAXG = N + 1
!       BDONE(1:N) = 0 ...
        BDONE(1:N) = FUDGE
        NCOL = 1
        do NG = 1, MAXG
           BIGP(NG) = NCOL
!          BINCL(1:N) = 0 ...
           BINCL(1:N) = FUDGE
           do 30 J = 1, N
!             Reject column J if it is already in a group.
!             IF (BDONE(J) == 1) GOTO 30 ...
              IBDONE = int(BDONE(J))
              if (IBDONE == 1) goto 30
!             Vertical extent of band = KBEGIN to KFINI.
!             KJ = number of nonzeros in column J.
!             BJA(K) = K implies nonzero at (k,j).
              KBEGIN = max(J-MU,1)
              KFINI = min(J+ML,N)
              KJ = 0
!             BJA(1:N) = 0 ...
              BJA(1:N) = FUDGE
!             Locate the row positions of the nonzeros in column J.
!             Restrict attention to the band:
              do 10 K = KBEGIN, KFINI
                 if (K < J) then
                    if (NSUPS > 0) then
                       do I = NSUPS, 1, -1
!                         KI = SUPDS(I) + J - 1
                          KI = J + 1 - SUPDS(I)
                          if (K == KI) then
                             KJ = KJ + 1
!                            BJA(K) = K ...
                             BJA(K) = real(K) + FUDGE
                          end if
                       end do
                    end if
                 elseif (K == J) then
                    KJ = KJ + 1
!                   BJA(K) = K ...
                    BJA(K) = real(K) + FUDGE
                 else
                    if (NSUBS > 0) then
                       do I = NSUBS, 1, -1
                          KI = SUBDS(I) + J -1
                          if (K == KI) then
                             KJ = KJ + 1
!                            BJA(K) = K ...
                             BJA(K) = real(K) + FUDGE
                          end if
                       end do
                    end if
              end if
 10           continue
!             At this point BJA contains the row numbers for
!             the nonzeros in column J.
              do 20 K = KBEGIN, KFINI
!                Reject column J if it overlaps any column
!                already in this group.
!                I = BJA(K)
                 I = int(BJA(K))
                 IBINCL = int(BINCL(I))
!                IF (BINCL(I) == 1 .AND. I == K) GOTO 30
                 if (IBINCL == 1 .and. I == K) goto 30
 20           end do
!             Accept column J into group NG.
              BJGP(NCOL) = J
              NCOL = NCOL + 1
!             BDONE(J) = 1 ...
              BDONE(J) = ONE_PLUS_FUDGE
              do K = 1, N
!                IF (I == K) BINCL(I) = 1 ...
                 I = int(BJA(K))
                 if (I == K) BINCL(I) = ONE_PLUS_FUDGE
              end do
 30        end do
!          Done if this group is empty (grouping is complete).
           if (NCOL == BIGP(NG)) goto 40
        end do

!       Should not get here since MAXG = N + 1.
!       Terminal error if not all columns were chosen
!       because MAXG too small.
        if (NCOL <= N) then
           MSG = 'An impossible error occurred in subroutine BGROUP.'
           call XERRDV(MSG,1600,2,0,0,0,0,ZERO,ZERO)
        end if

        NG = MAXG
 40     BNGRP = NG - 1

!       Trim BIGP to it's actual size if necessary.
        if (NG < MAXG) then
!          BJA(1:NG) = BIGP(1:NG) ...
           do I = 1, NG
              BJA(I) = real(BIGP(I)) + FUDGE
           end do
           deallocate (BIGP,STAT=IER)
           call CHECK_STAT(IER,780)
           allocate (BIGP(NG),STAT=IER)
           call CHECK_STAT(IER,790)
!          BIGP(1:NG) = BJA(1:NG) ...
           do I = 1, NG
              BIGP(I) = int(BJA(I))
           end do
        end if
        return

      end subroutine BGROUP
!_______________________________________________________________________

      subroutine BANDED_IAJA(N,ML,MU)
! ..
! Build the sparse structure descriptor arrays for a banded
! matrix if the nonzero diagonals are known.
! ..
!     Input:
!
!       N       = the order of the matrix (number of odes)
!       ML      = integer lower bandwidth
!       MU      = integer upper bandwidth
!
!     Output: (PRIVATE information used in DVJAC)
!
!       IAB     = IA descriptor array
!       JAB     = JA descriptor array
!
!    Note:
!    The PRIVATE banded information SUBDS, NSUBS, SUPDS, NSUPS,
!    ML, and MU must be defined before calling BANDED_IAJA.
!       ML        = lower bandwidth
!       MU        = upper bandwidth
!       NSUBS     = number of strict sub diagonals
!       SUBDS(I)  = row in which the Ith sub diagonal
!                   begins, I=1, ..., NSUBS
!       NSUPS     = number of strict super diagonals
!       SUPDS(I)  = column in which the Ith super diagonal
!                   begins, I=1, ..., NSUPS
! ..
        implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: N, ML, MU
! ..
! .. Local Scalars ..
        integer :: I, IER, J, K, KBEGIN, KFINI, KI, KJ, &
          NP1, NZBB
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
     intrinsic ALLOCATED, MAX, MIN
! ..
! .. FIRST EXECUTABLE STATEMENT BANDED_IAJA
! ..
!    Check for errors.

     if (.not.BUILD_IAJA) then
        MSG = 'BANDED_IAJA cannot be called with BUILD_IAJA = .FALSE.'
        call XERRDV(MSG,1610,2,0,0,0,0,ZERO,ZERO)
     end if

     if (NSUBS < 0) then
        MSG = 'NSUBS < 0 in BANDED_IAJA.'
        call XERRDV(MSG,1620,2,0,0,0,0,ZERO,ZERO)
     else
        if (NSUBS > 0) then
           if (NSUBS /= size(SUBDS)) then
              MSG = 'The size of the SUBDS array must'
              call XERRDV(MSG,1630,1,0,0,0,0,ZERO,ZERO)
              MSG = 'equal NSUBS in BANDED_IAJA.'
              call XERRDV(MSG,1630,2,0,0,0,0,ZERO,ZERO)
           end if
        end if
     end if

     if (NSUPS < 0) then
        MSG = 'NSUPS < 0 in BANDED_IAJA.'
        call XERRDV(MSG,1640,2,0,0,0,0,ZERO,ZERO)
     else
        if (NSUPS > 0) then
           if (NSUPS /= size(SUPDS)) then
              MSG = 'The size of the SUPDS array must'
              call XERRDV(MSG,1650,1,0,0,0,0,ZERO,ZERO)
              MSG = 'equal NSUPS in BANDED_IAJA.'
              call XERRDV(MSG,1650,2,0,0,0,0,ZERO,ZERO)
           end if
        end if
     end if

!    Allocate the necessary storage for the descriptor arrays.

!    Define the total number of elements in all diagonals.
     NP1 = N + 1
     NZB = (NSUBS + NSUPS + 1) * NP1 - 1
     if (NSUBS /= 0) then
        do I = 1, NSUBS
           NZB = NZB - SUBDS(I)
        end do
     end if
     if (NSUPS /= 0) then
        do I = 1, NSUPS
           NZB = NZB - SUPDS(I)
        end do
     end if
     if (allocated(IAB)) then
        deallocate(IAB,JAB,STAT=IER)
        call CHECK_STAT(IER,800)
     end if
     allocate(IAB(NP1),JAB(NZB),STAT=IER)
     call CHECK_STAT(IER,810)
     IAB(1) = 1
     NZBB = 0

!    For each column in the matrix...
     do J = 1, N
!       Vertical extent of band = KBEGIN to KFINI.
!       KJ = number of nonzeros in column J.
        KBEGIN = max(J-MU,1)
        KFINI = min(J+ML,N)
        KJ = 0
!       Locate the row positions of the nonzeros in column J.
!       (Restrict attention to the band.)
        IAB(J+1) = IAB(J)
!       For each row in the intersection of the band with
!       this column ...
        do K = KBEGIN, KFINI
!          Does column J intersect a super diagonal at (K,J)?
           if (K < J) then
              do I = NSUPS, 1, -1
                 KI = J + 1 - SUPDS(I)
                 if (K == KI) then
                    KJ = KJ + 1
                    IAB(J+1) = IAB(J+1) + 1
                    NZBB = NZBB + 1
                    JAB(NZBB) = K
                    goto 10
                 end if
              end do
           elseif (K == J) then
!             We are on the main diagonal.
              KJ = KJ + 1
              IAB(J+1) = IAB(J+1) + 1
              NZBB = NZBB + 1
              JAB(NZBB) = K
              goto 10
           else
!             Does column J intersect a sub diagonal at (K,J)?
              do I = NSUBS, 1, -1
                 KI = SUBDS(I) + J - 1
                 if (K == KI) then
                    KJ = KJ + 1
                    IAB(J+1) = IAB(J+1) + 1
                    NZBB = NZBB + 1
                    JAB(NZBB) = K
                    goto 10
                 end if
              end do
           end if
10         continue
        end do
     end do

     if (NZBB /= NZB) then
        MSG = 'NZBB (I1) is not equal to NZB (I2)'
        call XERRDV(MSG,1660,1,0,0,0,0,ZERO,ZERO)
        MSG = 'in BANDED_IAJA.'
        call XERRDV(MSG,1660,2,2,NZBB,NZB,0,ZERO,ZERO)
     end if
     return

   end subroutine BANDED_IAJA
!_______________________________________________________________________

      subroutine BANDED_GET_BJNZ(N,ML,MU,JCOL,JNZ,NZJ)
! ..
! Locate the nonzeros in a given column of a sparse banded matrix
! with known diagonals. This is a version of BANDED_IAJA modified
! to do only one column.
! ..
!     Input:
!
!       N       = the order of the matrix (number of odes)
!       ML      = integer lower bandwidth
!       MU      = integer upper bandwidth
!       JCOL    = column number between 1 and N
!       JZ      = integer array of length N
!
!     Output:
!
!       JNZ    = integer array of length N. If
!                 JNZ(K) is not 0, there is a
!                 nonzero at position (K,JCOL)
!       NZJ     = number of nozeros in column JCOL
!
!    Caution:
!    No parameter checking is done since this subroutine
!    will be called many times. Note that a number of
!    PRIVATE parameters must be set before calling this
!    subroutine.
! ..
        implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: N, ML, MU, JCOL
        integer, intent (OUT) :: NZJ
! ..
! .. Array Arguments ..
        integer, intent (OUT) :: JNZ(*)
! ..
! .. Local Scalars ..
        integer :: I, J, K, KBEGIN, KFINI, KI, KJ
! ..
! .. Intrinsic Functions ..
     intrinsic MAX, MIN
! ..
! .. FIRST EXECUTABLE STATEMENT BANDED_GET_BJNZ
! ..
        JNZ(1:N) = 0
        J = JCOL

!       Locate the row positions of the nonzeros in column J.

!       Vertical extent of band = KBEGIN to KFINI.
        KBEGIN = max(J-MU,1)
        KFINI = min(J+ML,N)
!       KJ = number of nonzeros in column J.
        KJ = 0
!       For each row in the intersection of the band with
!       this column ...
        do K = KBEGIN, KFINI
!          Does column J intersect a super diagonal at (K,J)?
           if (K < J) then
              do I = NSUPS, 1, -1
                 KI = J + 1 - SUPDS(I)
                 if (K == KI) then
                    KJ = KJ + 1
                    JNZ(KJ) = K
                    goto 10
                 end if
              end do
           elseif (K == J) then
!             We are on the main diagonal.
              KJ = KJ + 1
              JNZ(KJ) = K
              goto 10
           else
!             Does column J intersect a sub diagonal at (K,J)?
              do I = NSUBS, 1, -1
                 KI = SUBDS(I) + J - 1
                 if (K == KI) then
                    KJ = KJ + 1
                    JNZ(KJ) = K
                    goto 10
                 end if
              end do
           end if
10         continue
        end do
        NZJ = KJ
     return

   end subroutine BANDED_GET_BJNZ
!_______________________________________________________________________

! Beginning of Jacobian related routines that use MA28

      subroutine DVNLSS28(Y,YH,LDYH,SAVF,EWT,ACOR,IWM,WM,F,JAC, &
        NFLAG,ATOL,ITOL)
! ..
! This is the nonlinear system solver for MA28 based sparse solutions.
! ..
! Subroutine DVNLSS28 is a nonlinear system solver, which uses functional
! iteration or a chord (modified Newton) method. For the chord method
! direct linear algebraic system solvers are used. Subroutine DVNLSS28
! then handles the corrector phase of this integration package.
! Communication with DVNLSS28 is done with the following variables. (For
! more details, please see the comments in the driver subroutine.)
! Y          = The dependent variable, a vector of length N, input.
! YH         = The Nordsieck (Taylor) array, LDYH by LMAX, input
!              and output. On input, it contains predicted values.
! LDYH       = A constant >= N, the first dimension of YH, input.
! SAVF       = A work array of length N.
! EWT        = An error weight vector of length N, input.
! ACOR       = A work array of length N, used for the accumulated
!              corrections to the predicted y vector.
! WM,IWM     = Real and integer work arrays associated with matrix
!              operations in chord iteration (MITER /= 0).
! F          = Dummy name for user supplied routine for f.
! JAC        = Dummy name for user supplied Jacobian routine.
! NFLAG      = Input/output flag, with values and meanings as follows:
!              INPUT
!                  0 first call for this time step.
!                 -1 convergence failure in previous call to DVNLSS28.
!                 -2 error test failure in DVSTEP.
!              OUTPUT
!                  0 successful completion of nonlinear solver.
!                 -1 convergence failure or singular matrix.
!                 -2 unrecoverable error in matrix preprocessing
!                    (cannot occur here).
!                 -3 unrecoverable error in solution (cannot occur
!                    here).
! IPUP       = Own variable flag with values and meanings as follows:
!              0,          do not update the Newton matrix.
!              MITER \= 0  update Newton matrix, because it is the
!                          initial step, order was changed, the error
!                          test failed, or an update is indicated by
!                          the scalar RC or step counter NST.
! For more details, see comments in driver subroutine.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: ITOL, LDYH
        integer, intent (INOUT) :: NFLAG
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: ACOR(*), EWT(*), SAVF(*), &
          WM(*), Y(*), YH(LDYH,*)
        real (WP), intent (IN) :: ATOL(*)
        integer IWM(*)
        logical DUMMY
! ..
! .. Subroutine Arguments ..
        external F, JAC
! ..
! .. Local Scalars ..
        real (WP) :: ACNRMNEW, CSCALE, DCON, DEL, DELP
        integer :: I, IERPJ, IERSL, M
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, MAX, MIN
! ..
! .. FIRST EXECUTABLE STATEMENT DVNLSS28
! ..
! Get rid of a couple of needless compiler warning messages.
        DUMMY = .false.
        if (DUMMY) then
          WM(1) = ZERO
          IWM(1) = 0
        end if  
! On the first step, on a change of method order, or after a
! nonlinear convergence failure with NFLAG = -2, set IPUP = MITER
! to force a Jacobian update when MITER /= 0.
        if (JSTART==0) NSLP = 0
        if (NFLAG==0) ICF = 0
        if (NFLAG==-2) IPUP = MITER
        if ((JSTART==0) .or. (JSTART==-1)) IPUP = MITER
!       If this is functional iteration, set CRATE = 1 and drop
!       to 220.
        if (MITER==0) then
          CRATE = ONE
          goto 10
        end if

! RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1.
! When RC differs from 1 by more than CCMAX, IPUP is set to MITER
! to force DVJACS28 to be called, if a Jacobian is involved. In any
! case, DVJACS28 is called at least every MSBP steps.

        DRC = abs(RC-ONE)
        if (DRC>CCMAX .or. NST>=NSLP+MSBP) IPUP = MITER

! Up to MAXCOR corrector iterations are taken. A convergence test is
! made on the r.m.s. norm of each correction, weighted by the error
! weight vector EWT. The sum of the corrections is accumulated in the
! vector ACOR(i). The YH array is not altered in the corrector loop.

10      M = 0
        DELP = ZERO
! Original:
!       CALL DCOPY_F90(N,YH(1,1),1,Y,1)
        call DCOPY_F90(N,YH(1:N,1),1,Y(1:N),1)
        call F(N,TN,Y,SAVF)
        NFE = NFE + 1
        if (BOUNDS) then
          do I = 1, NDX
            if (abs(YH(IDX(I),1)-LB(I))<=ZERO) SAVF(IDX(I)) = &
              max(SAVF(IDX(I)),ZERO)
            if (abs(YH(IDX(I),1)-UB(I))<=ZERO) SAVF(IDX(I)) = &
              min(SAVF(IDX(I)),ZERO)
          end do
        end if
        if (IPUP<=0) goto 20

! If indicated, the matrix P = I - h*rl1*J is reevaluated and
! preprocessed before starting the corrector iteration. IPUP
! is set to 0 as an indicator that this has been done.

        call DVJACS28(Y,YH,LDYH,EWT,ACOR,SAVF,F,JAC,IERPJ,N, &
          ATOL,ITOL)
        IPUP = 0
        RC = ONE
        DRC = ZERO
        CRATE = ONE
        NSLP = NST
!     If matrix is singular, take error return to force cut in
!     step size.
        if (IERPJ/=0) goto 70
20      ACOR(1:N) = ZERO
!       This is a looping point for the corrector iteration.
30      if (MITER/=0) goto 40

! In the case of functional iteration, update Y directly from
! the result of the last function evaluation.

        SAVF(1:N) = RL1*(H*SAVF(1:N)-YH(1:N,2))
        Y(1:N) = SAVF(1:N) - ACOR(1:N)
        DEL = DVNORM(N,Y,EWT)
        Y(1:N) = YH(1:N,1) + SAVF(1:N)
        call DCOPY_F90(N,SAVF,1,ACOR,1)
        goto 50

! In the case of the chord method, compute the corrector error,
! and solve the linear system with that as right-hand side and
! P as coefficient matrix. The correction is scaled by the factor
! 2/(1+RC) to account for changes in h*rl1 since the last
! DVJACS28 call.

40      Y(1:N) = (RL1*H)*SAVF(1:N) - (RL1*YH(1:N,2)+ACOR(1:N))
        call DVSOLS28(Y,SAVF,IERSL)
        NNI = NNI + 1
        if (IERSL>0) goto 60
        if (METH==2 .and. abs(RC-ONE)>ZERO) then
          CSCALE = TWO/(ONE+RC)
          call DSCAL_F90(N,CSCALE,Y,1)
        end if
        DEL = DVNORM(N,Y,EWT)
        call DAXPY_F90(N,ONE,Y,1,ACOR,1)
        Y(1:N) = YH(1:N,1) + ACOR(1:N)

! Test for convergence. If M > 0, an estimate of the convergence
! rate constant is stored in CRATE, and this is used in the test.

50      if (M/=0) CRATE = max(CRDOWN*CRATE,DEL/DELP)
        DCON = DEL*min(ONE,CRATE)/TQ(4)
        if (DCON<=ONE) goto 80
        M = M + 1
        if (M==MAXCOR) goto 60
        if (M>=2 .and. DEL>RDIV*DELP) goto 60
        DELP = DEL
        call F(N,TN,Y,SAVF)
        NFE = NFE + 1
        if (BOUNDS) then
          do I = 1, NDX
            if (abs(YH(IDX(I),1)-LB(I))<=ZERO) SAVF(IDX(I)) = &
              max(SAVF(I),ZERO)
            if (abs(YH(IDX(I),1)-UB(I))<=ZERO) SAVF(IDX(I)) = &
              min(SAVF(I),ZERO)
          end do
        end if
        goto 30

60      if (MITER==0 .or. JCUR==1) goto 70
        ICF = 1
        IPUP = MITER
        goto 10

70      continue
        NFLAG = -1
        ICF = 2
        IPUP = MITER
        return

!       Return for successful step.
80      NFLAG = 0

!       Enforce bounds.
        if (BOUNDS) then
          CHANGED_ACOR = .false.
          if (M==0) then
            ACNRM = DEL
          else
            ACNRM = DVNORM(N,ACOR,EWT)
          end if
          if (MITER/=0) then
!           Since Y(:) = YH(:,1) + ACOR(:):
            do I = 1, NDX
              if (Y(IDX(I))<LB(I)) then
                CHANGED_ACOR = .true.
                ACOR(IDX(I)) = LB(I) - YH(IDX(I),1)
                SAVF(IDX(I)) = ACOR(IDX(I))
              end if
              if (Y(IDX(I))>UB(I)) then
                CHANGED_ACOR = .true.
                ACOR(IDX(I)) = UB(I) - YH(IDX(I),1)
                SAVF(IDX(I)) = ACOR(IDX(I))
              end if
            end do
          else
!           Since Y(:) = YH(:,1) + SAVF(:) and
!           since CALL DCOPY_F90(N, SAVF, 1, ACOR, 1) ...
            do I = 1, NDX
              if (Y(IDX(I))<LB(IDX(I))) then
                CHANGED_ACOR = .true.
                ACOR(IDX(I)) = LB(I) - YH(IDX(I),1)
              end if
              if (Y(IDX(I))>UB(IDX(I))) then
                CHANGED_ACOR = .true.
                ACOR(IDX(I)) = UB(I) - YH(IDX(I),1)
              end if
            end do
          end if
          if (CHANGED_ACOR) then
            if (M==0) then
              ACNRMNEW = DEL
            else
              ACNRMNEW = DVNORM(N,ACOR,EWT)
            end if
            ACNRM = max(ACNRM,ACNRMNEW)
          else
          end if
          NFLAG = 0
          JCUR = 0
          ICF = 0
        else
!         No projections are required.
          NFLAG = 0
          JCUR = 0
          ICF = 0
          if (M==0) ACNRM = DEL
          if (M>0) ACNRM = DVNORM(N,ACOR,EWT)
        end if
        return

      end subroutine DVNLSS28
!_______________________________________________________________________

      subroutine DVSOLS28(X,TEM,IERSL)
! ..
! Manage the solution of the MA28 based sparse linear system arising
! from a chord iteration.
! ..
! This routine solves the sparse linear system arising from a chord
! iteration. If MITER is 6 or 7, it calls MA28CD to accomplish this.
! Communication with DVSOLS28 uses the following variables:
! X = The right-hand side vector on input, and the solution vector
!     on output, of length N.
! TEM (=SAVF(*))
! IERSL = Output flag. IERSL = 0 if no trouble occurred.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (INOUT) :: IERSL
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: TEM(*), X(*)
! ..
! .. Local Scalars ..
        integer :: I
! ..
! .. FIRST EXECUTABLE STATEMENT DVSOLS28
! ..
        if (SCALE_MATRIX) then
           do I = 1, N
              X(I) = X(I) * RSCALEX(I)
           end do
        end if
        IERSL = 0
        call MA28CD(N,PMAT,LICN_ALL,ICN,IKEEP28,X,TEM,1)
        MA28CD_CALLS = MA28CD_CALLS + 1
        if (SCALE_MATRIX) then
           do I = 1, N
              X(I) = X(I) * CSCALEX(I)
           end do
        end if
        return

      end subroutine DVSOLS28
!_______________________________________________________________________

      subroutine DVJACS28(Y,YH,LDYH,EWT,FTEMP1,SAVF,F,JAC,IERPJ,N, &
        ATOL,ITOL)
! ..
! Compute and process P = I - H*RL1*J, where J is an approximation to
! the MA28 based sparse Jacobian.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (INOUT) :: IERPJ
        integer, intent (IN) :: LDYH, N, ITOL
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: EWT(*), FTEMP1(*), SAVF(*), Y(*), YH(LDYH,*)
        real (WP), intent (IN) :: ATOL(*)
! ..
! .. Subroutine Arguments ..
        external F, JAC
! ..
! .. Local Scalars ..
        real (WP) :: CON, FAC, HRL1, R, R0, SRUR
        integer :: I, IER, J, JER, JJ, JJ1, JJ2, K, K1, K2, MA28,       &
          MA28SAVE, MB28SAVE, NG, NZ
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, EXP, MAX, real
! ..
! .. FIRST EXECUTABLE STATEMENT DVJACS28
! ..
        IERPJ = 0

!       Structure determination

!       Calculate the sparsity structure if this is the first call to
!       DVJACS28 with ISTATE = 1 or if it is a continuation call with
!       ISTATE = 3.
        if (ISTATC==1 .or. ISTATC==3) then
          call DVPREPS(N,Y,YH,LDYH,DTEMP,EWT,F,JAC)
          ISTATC = 0
        end if
        JCUR = 0
        HRL1 = H*RL1
        CON = -HRL1

!       If MA28 = 4 the saved copy of the Jacobian will be used
!       to restore PMAT. If MA28 = 1,2,3 the Jacobian will be
!       recomputed. If MA28 = 1,2 the JVECT and ICN pointer arrays
!       will be defined and MA28AD will be called to decompose
!       PMAT. If MA28 = 3 the JVECT pointer array will be defined
!       and MA28BD will be called to decompose PMAT using the ICN
!       pointer array returned in the last call to MA28AD.

        MA28 = 4
        if (INEWJ==1 .or. MB28==0) MA28 = 3
        if (NST>=NSLJ+MSBJ) MA28 = 3
!       IF (ICF==1 .OR. ICF==2) MA28 = 3
        if (ICF==1 .and. DRC<CCMXJ) MA28 = 3
        if (ICF==2) MA28 = 3
        if (NST>=NSLG+MSBG) MA28 = 2
        if (JSTART==0 .or. JSTART==-1) MA28 = 1
        JSTART = 1
10      if (MA28<=2) NSLG = NST
        if (MA28<=3) NSLJ = NST

!       Analytical Sparse Jacobian

!       If MITER = 6, call JAC to evaluate J analytically, multiply
!       J by CON = -H*EL(1), and add the identity matrix to form P.
        if (MITER==6) then
          if (MA28==4) then
!           Reuse the saved Jacobian.
            NZ = IAN(N+1) - 1
            PMAT(1:NZ) = CON*JMAT(1:NZ)
            do K = 1, NZ
              if (JAN(K)==JVECT(K)) PMAT(K) = PMAT(K) + ONE
            end do
            goto 90
          end if
          JCUR = 1
          NJE = NJE + 1
          if (MA28==1 .or. MA28==2) then
            NZ = IAN(N+1) - 1
            call JAC(N,TN,Y,IAN,JAN,NZ,PMAT)
            NZ = IAN(N+1) - 1
            if (NZ>NZ_ALL) then
              MSG = 'DVODE_F90-- NZ > NZ_ALL in DVJACS28.'
              call XERRDV(MSG,1670,2,0,0,0,0,ZERO,ZERO)
            end if
!           Define column pointers for MA28AD.
            call SET_ICN(N,IAN,ICN)
            call CHECK_DIAG(N,IAN,JAN,ICN)
            PMAT(1:NZ) = CON*PMAT(1:NZ)
            do K = 1, NZ
              if (JAN(K)==ICN(K)) PMAT(K) = PMAT(K) + ONE
            end do
            goto 80
          else
!           MA28 = 3...
            NZ = IAN(N+1) - 1
            call JAC(N,TN,Y,IAN,JAN,NZ,PMAT)
            NZ = IAN(N+1) - 1
            if (NZ>NZ_ALL) then
              MSG = 'DVODE_F90-- NZ > NZ_ALL in DVJACS28.'
              call XERRDV(MSG,1680,2,0,0,0,0,ZERO,ZERO)
            end if
!           Define column pointers for MA28AD.
            call SET_ICN(N,IAN,JVECT)
            call CHECK_DIAG(N,IAN,JAN,JVECT)
            if (INEWJ/=1) JMAT(1:NZ) = PMAT(1:NZ)
            PMAT(1:NZ) = CON*PMAT(1:NZ)
            do K = 1, NZ
              if (JAN(K)==JVECT(K)) PMAT(K) = PMAT(K) + ONE
            end do
            goto 90
          end if
        end if

!       Finite Difference Sparse Jacobian

!       If MITER = 7, evaluate J numerically, multiply J by
!       CON, and add the identity matrix to form P.
        if (MITER==7) then
          if (MA28==4) then
!           Reuse the saved constant Jacobian.
            NZ = IAN(N+1) - 1
            PMAT(1:NZ) = CON*JMAT(1:NZ)
            do J = 1, N
              K1 = IAN(J)
              K2 = IAN(J+1) - 1
              do K = K1, K2
                I = JAN(K)
                if (I==J) PMAT(K) = PMAT(K) + ONE
              end do
            end do
            goto 90
          else
            NZ = IAN(N+1) - 1
            JCUR = 1
            if (.not.(J_IS_CONSTANT.and.J_HAS_BEEN_COMPUTED)) then
               if (USE_JACSP) then
!                 Approximate the Jacobian using Doug Salane's JACSP.
!                 The JPNTRDS and INDROWDS pointer arrays were defined
!                 in DVPREPS (and altered in DSM).
                  IOPTDS(1) = 2
                  IOPTDS(2) = 0
                  IOPTDS(3) = 1
                  IOPTDS(5) = 0
!                 INFORDS(4) was initialized in DVPREPS (and altered in
!                 the first call to JACSP).
                  LWKDS  = 3 * N
                  LIWKDS = 50 + N
                  NRFJACDS = NZ
                  NCFJACDS = 1

!                 Set flag to indicate how the YSCALE vector will be
!                 set for JACSP.
                  LIKE_ORIGINAL_VODE = .false.
!                 Calculate the YSCALEDS vector for JACSPDV.
                  if (LIKE_ORIGINAL_VODE) then
                     FAC = DVNORM(N,SAVF,EWT)
!                    JACSPDB multiplies YSCALEDS(*) BY UROUND**0.825:
!                    R0 = THOU*ABS(H)*UROUND*REAL(N)*FAC
                     R0 = THOU*abs(H)*real(N)*FAC
                     if (abs(R0)<=ZERO) R0 = ONE
!                    SRUR = WM1
                     do J = 1, N
!                       JACSPDB multiplies YSCALEDS(*) BY UROUND**0.825:
!                       R = MAX(ABS(Y(J)),R0/EWT(J))
                        R = max(abs(Y(J))/U325,(R0/EWT(J))*U125)
                        YSCALEDS(J) = R
                     end do
                  else
                     if (ITOL == 1 .or. ITOL == 3) then
                        do J = 1, N
                           YSCALEDS(J) = max(abs(Y(J)),ATOL(1),UROUND)
                        end do
                     else
                        do J = 1, N
                          YSCALEDS(J) = max(abs(Y(J)),ATOL(J),UROUND)
                        end do
                     end if
                  end if

                  call JACSPDB(F,N,TN,Y,SAVF,PMAT(1),NRFJACDS, &
                    YSCALEDS,FACDS,IOPTDS,WKDS,LWKDS,IWKDS,LIWKDS, &
                    MAXGRPDS,NGRPDS,JPNTRDS,INDROWDS)
                  NFE = NFE + IWKDS(7)
                  NJE = NJE + 1

                  do NG = 1, MAXGRPDS
                    JJ1 = IGP(NG)
                    JJ2 = IGP(NG+1) - 1
                    do JJ = JJ1, JJ2
                      J = JGP(JJ)
                      K1 = IAN(J)
                      K2 = IAN(J+1) - 1
                      do K = K1, K2
                        I = JAN(K)
                        goto (17,17,18) MA28
!                       Define the row pointers for MA28AD.
17                      JVECT(K) = I
!                       Define the column pointers for MA28AD.
                        ICN(K) = J
                        goto 19
!                       Define the column pointers for MA28AD.
18                      JVECT(K) = J
19                      continue
                        if (INEWJ==0) JMAT(K) = PMAT(K)
                        PMAT(K) = CON*PMAT(K)
                        if (I==J) PMAT(K) = PMAT(K) + ONE
                      end do
                    end do
                  end do
                  NFE = NFE + MAXGRPDS
               else
                  FAC = DVNORM(N,SAVF,EWT)
                  R0 = THOU*abs(H)*UROUND*real(N)*FAC
                  if (abs(R0)<=ZERO) R0 = ONE
                  SRUR = WM1
                  do NG = 1, NGP
                    JJ1 = IGP(NG)
                    JJ2 = IGP(NG+1) - 1
                    do JJ = JJ1, JJ2
                      J = JGP(JJ)
                      R = max(SRUR*abs(Y(J)),R0/EWT(J))
                      Y(J) = Y(J) + R
                    end do
                    call F(N,TN,Y,FTEMP1)
                    NFE = NFE + 1
                    do JJ = JJ1, JJ2
                      J = JGP(JJ)
                      Y(J) = YH(J,1)
                      R = max(SRUR*abs(Y(J)),R0/EWT(J))
                      FAC = ONE / R
                      K1 = IAN(J)
                      K2 = IAN(J+1) - 1
                      do K = K1, K2
                        I = JAN(K)
                        goto (20,20,30) MA28
!                       Define row pointers for MA28AD.
20                      JVECT(K) = I
!                       Define column pointers for MA28AD.
                        ICN(K) = J
                        goto 40
!                       Define column pointers for MA28AD.
30                      JVECT(K) = J
40                      PMAT(K) = (FTEMP1(I)-SAVF(I)) * FAC
                        if (INEWJ==0) JMAT(K) = PMAT(K)
                        PMAT(K) = CON*PMAT(K)
                        if (I==J) PMAT(K) = PMAT(K) + ONE
                      end do
                    end do
                  end do
                  NFE = NFE + NGP
                  NJE = NJE + 1
               end if
               if (J_IS_CONSTANT) J_HAS_BEEN_COMPUTED = .true.
            else
!              Do not recompute the constant Jacobian.
!              Reuse the saved Jacobian.
               NZ = IAN(N+1) - 1
               PMAT(1:NZ) = CON*JMAT(1:NZ)
               do J = 1, N
                 K1 = IAN(J)
                 K2 = IAN(J+1) - 1
                 do K = K1, K2
                   I = JAN(K)
                   if (I==J) PMAT(K) = PMAT(K) + ONE
                   goto (50,50,60) MA28
!                  Define row pointers for MA28AD.
50                 JVECT(K) = I
!                  Define column pointers for MA28AD.
                   ICN(K) = J
                   goto 70
!                  Define column pointers for MA28AD.
60                 JVECT(K) = J
70                 continue
                 end do
               end do
            end if
            goto (80,80,90) MA28
          end if
        end if

!       MA28AD does an LU factorization based on a pivotal strategy
!       designed to compromise between maintaining sparsity and
!       controlling loss of accuracy due to roundoff error. Unless
!       magnitudes of Jacobian elements change so as to invalidate
!       choice of pivots, MA28AD need only be called at beginning
!       of the integration.

80      continue
        if (SCALE_MATRIX) then
!          MA19AD computes scaling factors for the iteration matrix.
           call MC19AD(N,NZ,PMAT,JAN,ICN,RSCALEX,CSCALEX,WSCALEX)
           MC19AD_CALLS = MC19AD_CALLS + 1
           do I =1, N
              RSCALEX(I) = exp(RSCALEX(I))
              CSCALEX(I) = exp(CSCALEX(I))
           end do
           do K = 1, NZ
              I = JAN(K)
              J = ICN(K)
              PMAT(K) = PMAT(K) * RSCALEX(I) * CSCALEX(J)
           end do
        end if
        OK_TO_CALL_MA28 = .true.
        call MA28AD(N,NZ,PMAT,LICN_ALL,JAN,LIRN_ALL,ICN,U_PIVOT, &
          IKEEP28,IW28,FTEMP1,IER)
        OK_TO_CALL_MA28 = .false.
        MA28AD_CALLS = MA28AD_CALLS + 1
        MA28SAVE = MA28
        MB28SAVE = MB28
        MB28 = 0
        NLU = NLU + 1
        MAX_MINIRN = max(MAX_MINIRN,MINIRN)
        MAX_MINICN = max(MAX_MINICN,MINICN)
!       IER = -1: Numerically singular Jacobian
!       IER = -2: Structurally singular Jacobian
        if (IER==-1 .or. IER==-2) IERPJ = 1
        if (IER==-3) then
!         LIRN_ALL is not large enough.
          if (LP /= 0) then
             MSG = 'LIRN_ALL (=I1) is not large enough.'
             call XERRDV(MSG,1690,1,0,0,0,0,ZERO,ZERO)
             MSG = 'Allocating more space for another try.'
             call XERRDV(MSG,1690,1,1,LIRN_ALL,0,0,ZERO,ZERO)
          end if
!         Allocate more space for JAN and JVECT and try again.
          LIRN_ALL = LIRN_ALL + max(max(1000,ELBOW_ROOM*NZ_SWAG),10*N)
          LIRN_ALL = max(LIRN_ALL,(11*MINIRN)/10)
          if (LIRN_ALL>MAX_ARRAY_SIZE) then
            MSG = 'Maximum array size exceeded. Stopping.'
            call XERRDV(MSG,1700,2,0,0,0,0,ZERO,ZERO)
          end if
          deallocate (JAN,STAT=JER)
          call CHECK_STAT(JER,820)
          allocate (JAN(LIRN_ALL),STAT=JER)
          call CHECK_STAT(JER,830)
          if (MITER==7) then
            JAN(1:NZ) = JVECT(1:NZ)
          end if
          deallocate (JVECT,STAT=JER)
          call CHECK_STAT(JER,840)
          allocate (JVECT(LIRN_ALL),STAT=JER)
          call CHECK_STAT(JER,850)
          if (MITER==7) then
            JVECT(1:NZ) = JAN(1:NZ)
          end if
          MA28 = MA28SAVE
          MB28 = MB28SAVE
          NLU = NLU - 1
!         Since PMAT has changed, it must be restored:
          if (J_IS_CONSTANT) J_HAS_BEEN_COMPUTED = .false.
          goto 10
        end if
        if (IER==-4 .or. IER==-5 .or. IER==-6) then
!         LICN_ALL is not large enough.
          if (LP /= 0) then
             MSG = 'LICN_ALL (=I1) is not large enough.'
             call XERRDV(MSG,1710,1,0,0,0,0,ZERO,ZERO)
             MSG = 'Allocating more space for another try.'
             call XERRDV(MSG,1710,1,1,LICN_ALL,0,0,ZERO,ZERO)
          end if
!         Allocate more space for JAN and JVECT and try again.
          LICN_ALL = LICN_ALL + max(max(1000,ELBOW_ROOM*NZ_SWAG),10*N)
          LICN_ALL = max(LICN_ALL,(11*MINICN)/10)
          if (LICN_ALL>MAX_ARRAY_SIZE) then
            MSG = 'Maximum array size exceeded. Stopping.'
            call XERRDV(MSG,1720,2,0,0,0,0,ZERO,ZERO)
          end if
          deallocate (PMAT,ICN,STAT=JER)
          call CHECK_STAT(JER,860)
          allocate (PMAT(LICN_ALL),ICN(LICN_ALL),STAT=JER)
          PMAT(1:LICN_ALL) = ZERO
          call CHECK_STAT(JER,870)
          if (MITER==7) then
            JAN(1:NZ) = JVECT(1:NZ)
          end if
          MA28 = MA28SAVE
          MB28 = MB28SAVE
          NLU = NLU - 1
          if (J_IS_CONSTANT) J_HAS_BEEN_COMPUTED = .false.
          goto 10
        end if
        if (MITER/=7) return
        JAN(1:NZ) = JVECT(1:NZ)
        return

!       MA28BD uses the pivot sequence generated by an earlier call
!       to MA28AD to factor a new matrix of the same structure.

90      continue

        if (SCALE_MATRIX) then
           do K =1, NZ
              I = JAN(K)
              J = JVECT(K)
              PMAT(K) = PMAT(K) * RSCALEX(I) * CSCALEX(J)
           end do
        end if
        call MA28BD(N,NZ,PMAT,LICN_ALL,JAN,JVECT,ICN,IKEEP28,IW28, &
          FTEMP1,IER)
        MA28BD_CALLS = MA28BD_CALLS + 1
        MB28 = 1
        NLU = NLU + 1
!       IER = -2 : The matrix is numerically singular. The MA28AD
!                  pivot sequence leads to a zero pivot, that is,
!                  to one for which the ratio of it to the smallest
!                  element in the row is less than EPS.
!       IER = -13: The matrix is structurally singular.

        if (REDO_PIVOT_SEQUENCE) then
!          Force MA28AD to calculate a new pivot sequence.
           if (IER/=-13 .and. IER/=-2 .and. IER<=0) return
           if (IER==-2) MA28 = 1
           if (IER==-13) MA28 = 1
           if (J_IS_CONSTANT) J_HAS_BEEN_COMPUTED = .false.
        else
           if (IER==-2) IERPJ = 1
           if (IER/=-13 .and. IER<=0) return
           if (IER==-13) MA28 = 1
           if (J_IS_CONSTANT) J_HAS_BEEN_COMPUTED = .false.
        end if
        if (IER>0) MA28 = 2
        if (IER==-13 .and. MITER==7 .and. MOSS==2) then
!          Recompute the sparsity structure.
           call DVRENEW(N,Y,SAVF,EWT,F)
           if (J_IS_CONSTANT) J_HAS_BEEN_COMPUTED = .false.
        end if
        goto 10

      end subroutine DVJACS28
! End of Jacobian related routines that use MA28
!_______________________________________________________________________

      subroutine SET_ICN(N,IA,ICN)
! ..
! Define the column locations of nonzero elements for a sparse
! matrix.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: N
! ..
! .. Array Arguments ..
        integer, intent (IN) :: IA(*)
        integer, intent (INOUT) :: ICN(*)
! ..
! .. Local Scalars ..
        integer :: J, KMAX, KMIN
! ..
! .. FIRST EXECUTABLE STATEMENT SET_ICN
! ..
        KMIN = 1
        do J = 1, N
          KMAX = IA(J+1) - 1
          ICN(KMIN:KMAX) = J
          KMIN = KMAX + 1
        end do

      end subroutine SET_ICN
!_______________________________________________________________________

      subroutine CHECK_DIAG(N,IA,JA,ICN)
! ..
! Check that the diagonal is included in the sparse matrix
! description arrays.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: N
! ..
! .. Array Arguments ..
        integer, intent (IN) :: IA(*), ICN(*), JA(*)
! ..
! .. Local Scalars ..
        integer :: J, K, KMAX, KMIN
        character (80) :: MSG
! ..
! .. FIRST EXECUTABLE STATEMENT CHECK_DIAG
! ..
        KMIN = 1
        do J = 1, N
          KMAX = IA(J+1) - 1
          do K = KMIN, KMAX
            if (JA(K)==ICN(K)) goto 10
          end do
          MSG = 'In CHECK_DIAG, the diagonal is not present.'
          call XERRDV(MSG,1730,2,0,0,0,0,ZERO,ZERO)
10        continue
          KMIN = KMAX + 1
        end do

      end subroutine CHECK_DIAG
!_______________________________________________________________________

      subroutine DVCHECK(JOB,G,NEQ,Y,YH,NYH,G0,G1,GX,IRT)
! ..
! Check for the presence of a root in the vicinity of the current T,
! in a manner depending on the input flag JOB, and call DVROOTS to
! locate the root as precisely as possible.
! ..
! This subroutine is essentially DRCHEK from LSODAR.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (INOUT) :: IRT
        integer, intent (IN) :: JOB, NEQ, NYH
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: G0(*), G1(*), GX(*), Y(*), YH(NYH,*)
! ..
! .. Subroutine Arguments ..
        external G
! ..
! .. Local Scalars ..
        real (WP) :: HMING, T1, TEMP1, TEMP2, X
        integer :: I, IFLAG, JFLAG
        logical :: ZROOT
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, MAX, MIN, SIGN
! ..
! In addition to variables described previously, DVCHECK
! uses the following for communication:
! JOB    = integer flag indicating type of call:
!          JOB = 1 means the problem is being initialized, and DVCHECK
!                  is to look for a root at or very near the initial T.
!          JOB = 2 means a continuation call to the solver was just
!                  made, and DVCHECK is to check for a root in the
!                  relevant part of the step last taken.
!          JOB = 3 means a successful step was just taken, and DVCHECK
!                  is to look for a root in the interval of the step.
! G0     = array of length NG, containing the value of g at T = T0ST.
!          G0 is input for JOB >= 2, and output in all cases.
! G1,GX  = arrays of length NG for work space.
! IRT    = completion flag:
!          IRT = 0  means no root was found.
!          IRT = -1 means JOB = 1 and a root was found too near to T.
!          IRT = 1  means a legitimate root was found (JOB = 2 or 3).
!                   On return, T0ST is the root location, and Y is the
!                   corresponding solution vector.
! T0ST   = value of T at one endpoint of interval of interest. Only
!          roots beyond T0ST in the direction of integration are sought.
!          T0ST is input if JOB >= 2, and output in all cases.
!          T0ST is updated by DVCHECK, whether a root is found or not.
! TLAST  = last value of T returned by the solver (input only).
! TOUTC  = copy of TOUT(input only).
! IRFND  = input flag showing whether the last step taken had a root.
!          IRFND = 1 if it did, = 0 if not.
! ITASKC = copy of ITASK (input only).
! NGC    = copy of NG (input only).
! ..
! .. FIRST EXECUTABLE STATEMENT DVCHECK
! ..
        IRT = 0
        JROOT(1:NGC) = 0
        HMING = (abs(TN)+abs(H))*UROUND*HUN

        goto (10,30,80) JOB

!       Evaluate g at initial T, and check for zero values.
10      continue
        T0ST = TN
        call G(NEQ,T0ST,Y,NGC,G0)
        NGE = 1
        ZROOT = .false.
        do I = 1, NGC
          if (abs(G0(I))<=ZERO) ZROOT = .true.
        end do
        if (.not.ZROOT) goto 20
!       g has a zero at T. Look at g at T + (small increment).
!       TEMP1 = SIGN(HMING, H)
!       T0ST = T0ST + TEMP1
!       TEMP2 = TEMP1 / H
        TEMP2 = max(HMING/abs(H),TENTH)
        TEMP1 = TEMP2*H
        T0ST = T0ST + TEMP1
        Y(1:N) = Y(1:N) + TEMP2*YH(1:N,2)
        call G(NEQ,T0ST,Y,NGC,G0)
        NGE = NGE + 1
        ZROOT = .false.
        do I = 1, NGC
          if (abs(G0(I))<=ZERO) ZROOT = .true.
        end do
        if (.not.ZROOT) goto 20
!       g has a zero at T and also close to T. Take error return.
        IRT = -1
        return

20      continue
        return

30      continue
        if (IRFND==0) goto 70
!       If a root was found on the previous step, evaluate G0 = g(T0ST).
        call DVINDY_CORE(T0ST,0,YH,NYH,Y,IFLAG)
        if (BOUNDS) then
          do I = 1, NDX
            Y(IDX(I)) = max(Y(IDX(I)),LB(I))
            Y(IDX(I)) = min(Y(IDX(I)),UB(I))
          end do
        end if
        call G(NEQ,T0ST,Y,NGC,G0)
        NGE = NGE + 1
        ZROOT = .false.
        do I = 1, NGC
          if (abs(G0(I))<=ZERO) ZROOT = .true.
        end do
        if (.not.ZROOT) goto 70
!       g has a zero at T0ST. Look at g at T + (small increment).
        TEMP1 = sign(HMING,H)
        T0ST = T0ST + TEMP1
        if ((T0ST-TN)*H<ZERO) goto 40
        TEMP2 = TEMP1/H
        Y(1:N) = Y(1:N) + TEMP2*YH(1:N,2)
        goto 50
40      call DVINDY_CORE(T0ST,0,YH,NYH,Y,IFLAG)
        if (BOUNDS) then
          do I = 1, NDX
            Y(IDX(I)) = max(Y(IDX(I)),LB(I))
            Y(IDX(I)) = min(Y(IDX(I)),UB(I))
          end do
        end if
50      call G(NEQ,T0ST,Y,NGC,G0)
        NGE = NGE + 1
        ZROOT = .false.
        do 60 I = 1, NGC
          if (abs(G0(I))>ZERO) goto 60
          JROOT(I) = 1
          ZROOT = .true.
60      end do
        if (.not.ZROOT) goto 70
!       g has a zero at T0ST and also close to T0ST. Return root.
        IRT = 1
        return
!       G0 has no zero components. Proceed to check relevant interval.
70      if (abs(TN-TLAST)<=ZERO) goto 130

80      continue
!       Set T1 to TN or TOUTC, whichever comes first, and get g at T1.
        if (ITASKC==2 .or. ITASKC==3 .or. ITASKC==5) goto 90
        if ((TOUTC-TN)*H>=ZERO) goto 90
        T1 = TOUTC
        if ((T1-T0ST)*H<=ZERO) goto 130
        call DVINDY_CORE(T1,0,YH,NYH,Y,IFLAG)
        if (BOUNDS) then
          do I = 1, NDX
            Y(IDX(I)) = max(Y(IDX(I)),LB(I))
            Y(IDX(I)) = min(Y(IDX(I)),UB(I))
          end do
        end if
        goto 100
90      T1 = TN
        do I = 1, N
          Y(I) = YH(I,1)
        end do
100     call G(NEQ,T1,Y,NGC,G1)
        NGE = NGE + 1
!       Call DVROOTS to search for root in interval from T0ST to T1.
        JFLAG = 0
110     continue
        call DVROOTS(NGC,HMING,JFLAG,T0ST,T1,G0,G1,GX,X)
        if (JFLAG>1) goto 120
        call DVINDY_CORE(X,0,YH,NYH,Y,IFLAG)
        if (BOUNDS) then
          do I = 1, NDX
            Y(IDX(I)) = max(Y(IDX(I)),LB(I))
            Y(IDX(I)) = min(Y(IDX(I)),UB(I))
          end do
        end if
        call G(NEQ,X,Y,NGC,GX)
        NGE = NGE + 1
        goto 110
120     T0ST = X
        call DCOPY_F90(NGC,GX,1,G0,1)
        if (JFLAG==4) goto 130
!       Found a root. Interpolate to X and return.
        call DVINDY_CORE(X,0,YH,NYH,Y,IFLAG)
        if (BOUNDS) then
          do I = 1, NDX
            Y(IDX(I)) = max(Y(IDX(I)),LB(I))
            Y(IDX(I)) = min(Y(IDX(I)),UB(I))
          end do
        end if
        IRT = 1
        return
130     continue
        return

      end subroutine DVCHECK
!_______________________________________________________________________

      subroutine DVROOTS(NG,HMIN,JFLAG,X0,X1,G0,G1,GX,X)
! ..
! Perform root finding for DVODE_F90.
! ..
! This is essentially subroutine DROOTS from LSODAR.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP), intent (IN) :: HMIN
        real (WP), intent (INOUT) :: X, X0, X1
        integer, intent (INOUT) :: JFLAG
        integer, intent (IN) :: NG
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: G0(NG), G1(NG), GX(NG)
! ..
! .. Local Scalars ..
        real (WP) :: T2, TMAX
        integer :: I, IMXOLD, NXLAST
        logical :: SGNCHG, XROOT, ZROOT
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, SIGN
! ..
! This subroutine finds the leftmost root of a set of arbitrary
! functions gi(x) (i = 1,...,NG) in an interval (X0,X1). Only roots
! of odd multiplicity (i.e. changes of sign of the gi) are found.
! Here the sign of X1 - X0 is arbitrary, but is constant for a given
! problem, and 'leftmost' means nearest to X0.The values of the
! vector-valued function g(x) = (gi, i=1...NG) are communicated
! through the call sequence of DVROOTS. The method used is the
! Illinois algorithm.

! Reference:
! Kathie L. Hiebert and Lawrence F. Shampine, Implicitly Defined
! Output Points for Solutions of ODEs, Sandia Report SAND/80-0180,
! February 1980.

! Description of parameters.

! NG     = number of functions gi, or the number of components of
!          the vector valued function g(x). Input only.

! HMIN   = resolution parameter in X. Input only. When a root is
!          found, it is located only to within an error of HMIN in X.
!          Typically, HMIN should be set to something on the order of
!               100 * UROUND * MAX(ABS(X0),ABS(X1)),
!          where UROUND is the unit roundoff of the machine.

! JFLAG  = integer flag for input and output communication.

!          On input, set JFLAG = 0 on the first call for the problem,
!          and leave it unchanged until the problem is completed.
!          (The problem is completed when JFLAG >= 2 on return.)

!          On output, JFLAG has the following values and meanings:
!          JFLAG = 1 means DVROOTS needs a value of g(x). Set GX = g(X)
!                    and call DVROOTS again.
!          JFLAG = 2 means a root has been found. The root is
!                    at X, and GX contains g(X). (Actually, X is the
!                    rightmost approximation to the root on an interval
!                    (X0,X1) of size HMIN or less.)
!          JFLAG = 3 means X = X1 is a root, with one or more of the gi
!                    being zero at X1 and no sign changes in (X0,X1).
!                    GX contains g(X) on output.
!          JFLAG = 4 means no roots (of odd multiplicity) were
!                    found in (X0,X1) (no sign changes).

! X0,X1  = endpoints of the interval where roots are sought.
!          X1 and X0 are input when JFLAG = 0 (first call), and
!          must be left unchanged between calls until the problem is
!          completed. X0 and X1 must be distinct, but X1 - X0 may be
!          of either sign. However, the notion of 'left' and 'right'
!          will be used to mean nearer to X0 or X1, respectively.
!          When JFLAG >= 2 on return, X0 and X1 are output, and
!          are the endpoints of the relevant interval.

! G0,G1  = arrays of length NG containing the vectors g(X0) and g(X1),
!          respectively. When JFLAG = 0, G0 and G1 are input and
!          none of the G0(i) should be zero.
!          When JFLAG >= 2 on return, G0 and G1 are output.

! GX     = array of length NG containing g(X). GX is input
!          when JFLAG = 1, and output when JFLAG >= 2.

! X      = independent variable value. Output only.
!          When JFLAG = 1 on output, X is the point at which g(x)
!          is to be evaluated and loaded into GX.
!          When JFLAG = 2 or 3, X is the root.
!          When JFLAG = 4, X is the right endpoint of the interval, X1.

! JROOT  = integer array of length NG. Output only.
!          When JFLAG = 2 or 3, JROOT indicates which components
!          of g(x) have a root at X. JROOT(i) is 1 if the i-th
!          component has a root, and JROOT(i) = 0 otherwise.
! ..
! .. FIRST EXECUTABLE STATEMENT DVROOTS
! ..
        if (JFLAG==1) goto 90
!       JFLAG /= 1. Check for change in sign of g or zero at X1.
        IMAX = 0
        TMAX = ZERO
        ZROOT = .false.
        do 20 I = 1, NG
          if (abs(G1(I))>ZERO) goto 10
          ZROOT = .true.
          goto 20
!         At this point, G0(i) has been checked and cannot be zero.
10        if (abs(sign(ONE,G0(I))-sign(ONE,G1(I)))<=ZERO) goto 20
          T2 = abs(G1(I)/(G1(I)-G0(I)))
          if (T2<=TMAX) goto 20
          TMAX = T2
          IMAX = I
20      end do
        if (IMAX>0) goto 30
        SGNCHG = .false.
        goto 40
30      SGNCHG = .true.
40      if (.not.SGNCHG) goto 200
!       There is a sign change. Find the first root in the interval.
        XROOT = .false.
        NXLAST = 0
        LAST = 1

!       Repeat until the first root in the interval is found. Loop point.
50      continue
        if (XROOT) goto 170
        if (NXLAST==LAST) goto 60
        ALPHA = ONE
        goto 80
60      if (LAST==0) goto 70
        ALPHA = HALF*ALPHA
        goto 80
70      ALPHA = TWO*ALPHA
80      X2 = X1 - (X1-X0)*G1(IMAX)/(G1(IMAX)-ALPHA*G0(IMAX))
!       IF ((ABS(X2 - X0) < HMIN) .AND. (ABS(X1 - X0) > TEN * &
!       HMIN)) X2 = X0 + PT1 * (X1 - X0)

!       If X2 is too close to X0 or X1, adjust it inward,
!       by a fractional distance that is between 0.1 and 0.5.
        if (abs(X2-X0)<HALF*HMIN) then
          FRACINT = abs(X1-X0)/HMIN
          FRACSUB = TENTH
          if (FRACINT<=FIVE) FRACSUB = HALF/FRACINT
          X2 = X0 + FRACSUB*(X1-X0)
        end if
        if (abs(X1-X2)<HALF*HMIN) then
          FRACINT = abs(X1-X0)/HMIN
          FRACSUB = TENTH
          if (FRACINT<=FIVE) FRACSUB = HALF/FRACINT
          X2 = X1 - FRACSUB*(X1-X0)
        end if

        JFLAG = 1
        X = X2
!       Return to the calling routine to get a value of GX = g(X).
        return
!       Check to see in which interval g changes sign.
90      IMXOLD = IMAX
        IMAX = 0
        TMAX = ZERO
        ZROOT = .false.
        do 110 I = 1, NG
          if (abs(GX(I))>ZERO) goto 100
          ZROOT = .true.
          goto 110
!         Neither G0(i) nor GX(i) can be zero at this point.
100       if (abs(sign(ONE,G0(I))-sign(ONE,GX(I)))<=ZERO) goto 110
          T2 = abs(GX(I)/(GX(I)-G0(I)))
          if (T2<=TMAX) goto 110
          TMAX = T2
          IMAX = I
110     end do
        if (IMAX>0) goto 120
        SGNCHG = .false.
        IMAX = IMXOLD
        goto 130
120     SGNCHG = .true.
130     NXLAST = LAST
        if (.not.SGNCHG) goto 140
!       Sign change between X0 and X2, so replace X1 with X2.
        X1 = X2
        call DCOPY_F90(NG,GX,1,G1,1)
        LAST = 1
        XROOT = .false.
        goto 160
140     if (.not.ZROOT) goto 150
!       Zero value at X2 and no sign change in (X0,X2), so X2 is a root.
        X1 = X2
        call DCOPY_F90(NG,GX,1,G1,1)
        XROOT = .true.
        goto 160
!       No sign change between X0 and X2. Replace X0 with X2.
150     continue
        call DCOPY_F90(NG,GX,1,G0,1)
        X0 = X2
        LAST = 0
        XROOT = .false.
160     if (abs(X1-X0)<=HMIN) XROOT = .true.
        goto 50

!       Return with X1 as the root. Set JROOT. Set X = X1 and GX = G1.
170     JFLAG = 2
        X = X1
        call DCOPY_F90(NG,G1,1,GX,1)
        do 190 I = 1, NG
          JROOT(I) = 0
          if (abs(G1(I))>ZERO) goto 180
          JROOT(I) = 1
          goto 190
180       if (abs(sign(ONE,G0(I))-sign(ONE,G1(I)))>ZERO) JROOT(I) = 1
190     end do
        return

!       No sign change in the interval. Check for zero at right endpoint.
200     if (.not.ZROOT) goto 210

!       Zero value at X1 and no sign change in (X0,X1). Return JFLAG = 3.
        X = X1
        call DCOPY_F90(NG,G1,1,GX,1)
        do I = 1, NG
          JROOT(I) = 0
          if (abs(G1(I))<=ZERO) JROOT(I) = 1
        end do
        JFLAG = 3
        return

!       No sign changes in this interval. Set X = X1, return JFLAG = 4.
210     call DCOPY_F90(NG,G1,1,GX,1)
        X = X1
        JFLAG = 4
        return

      end subroutine DVROOTS
!_______________________________________________________________________

      subroutine DVNRDP(Y,IYDIM,NEQN,NQ)
! ..
! Retract the Nordsieck array (undo prediction).
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer :: IYDIM, NEQN, NQ
! ..
! .. Array Arguments ..
        real (WP) :: Y(*)
! ..
! .. Local Scalars ..
        integer :: I, J, J1, J2
! ..
! .. FIRST EXECUTABLE STATEMENT DVNRDP
! ..
        do J1 = 1, NQ
          do J2 = J1, NQ
            J = (NQ+J1) - J2
            do I = 1, NEQN
! Original:
!             Y(I,J) = Y(I,J) + Y(I,J+1)
              Y(I+(J-1)*IYDIM) = Y(I+(J-1)*IYDIM) + Y(I+J*IYDIM)
            end do
          end do
        end do
        return

      end subroutine DVNRDP
!_______________________________________________________________________

      subroutine DVNRDN(Y,IYDIM,NEQN,NQ)
! ..
! Apply the Nordsieck array (predict).
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer :: IYDIM, NEQN, NQ
! ..
! .. Array Arguments ..
        real (WP) :: Y(*)
! ..
! .. Local Scalars ..
        integer :: I, J, J1, J2
! ..
! .. FIRST EXECUTABLE STATEMENT DVNRDN
! ..
        do J1 = 1, NQ
          do J2 = J1, NQ
            J = (NQ+J1) - J2
            do I = 1, NEQN
! Original:
!             Y(I,J) = Y(I,J) - Y(I,J+1)
              Y(I+(J-1)*IYDIM) = Y(I+(J-1)*IYDIM) - Y(I+J*IYDIM)
            end do
          end do
        end do
        return

      end subroutine DVNRDN
!_______________________________________________________________________

      subroutine DVNRDS(Y,IYDIM,NEQN,L,RH)
! ..
! Scale the Nordsieck array.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP) :: RH
        integer :: IYDIM, L, NEQN
! ..
! .. Array Arguments ..
        real (WP) :: Y(*)
! ..
! .. Local Scalars ..
        real (WP) :: R1
        integer :: I, J
! ..
! .. FIRST EXECUTABLE STATEMENT DVNRDS
! ..
        R1 = ONE
        do J = 2, L
          R1 = R1*RH
          do I = 1, NEQN
! Original:
!           Y(I,J) = Y(I,J)*R1
            Y(I+(J-1)*IYDIM) = Y(I+(J-1)*IYDIM)*R1
          end do
        end do
        return

      end subroutine DVNRDS
!_______________________________________________________________________

      subroutine RELEASE_ARRAYS
! ..
! Deallocate any allocated arrays and determine how much storage was
! used (not called by DVODE_F90).
! ..
     implicit none
! ..
! .. Local Scalars ..
        integer :: IER, II, IR
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
        intrinsic ALLOCATED, SIZE
! ..
!_______________________________________________________________________
! *****MA48 build change point. Insert these statements.
!     INTEGER INFO
!     INTEGER ISIZE
!     COMMON /MA48SIZE/ ISIZE
!_______________________________________________________________________
! ..
! .. FIRST EXECUTABLE STATEMENT RELEASE_ARRAYS
! ..
        IR = 0
        if (allocated(ACOR)) then
          IR = IR + size(ACOR)
          deallocate (ACOR,STAT=IER)
          call CHECK_STAT(IER,880)
        end if
        if (allocated(CSCALEX)) then
          IR = IR + size(CSCALEX) + size(RSCALEX) + size(WSCALEX)
          deallocate (CSCALEX,RSCALEX,WSCALEX,STAT=IER)
          call CHECK_STAT(IER,890)
        end if
        if (allocated(DTEMP)) then
          IR = IR + size(DTEMP)
          deallocate (DTEMP,STAT=IER)
          call CHECK_STAT(IER,900)
        end if
        if (allocated(EWT)) then
          IR = IR + size(EWT)
          deallocate (EWT,STAT=IER)
          call CHECK_STAT(IER,910)
        end if
        if (allocated(FACDS)) then
          IR = IR + size(FACDS)
          deallocate (FACDS,STAT=IER)
          call CHECK_STAT(IER,920)
        end if
        if (allocated(FPTEMP)) then
          IR = IR + size(FPTEMP)
          deallocate (FPTEMP,STAT=IER)
          call CHECK_STAT(IER,930)
        end if
        if (allocated(FTEMP)) then
          IR = IR + size(FTEMP)
          deallocate (FTEMP,STAT=IER)
          call CHECK_STAT(IER,940)
        end if
        if (allocated(FTEMP1)) then
          IR = IR + size(FTEMP1)
          deallocate (FTEMP1,STAT=IER)
          call CHECK_STAT(IER,950)
        end if
        if (allocated(G0)) then
          IR = IR + size(G0)
          deallocate (G0,STAT=IER)
          call CHECK_STAT(IER,960)
        end if
        if (allocated(G1)) then
          IR = IR + size(G1)
          deallocate (G1,STAT=IER)
          call CHECK_STAT(IER,970)
        end if
        if (allocated(GX)) then
          IR = IR + size(GX)
          deallocate (GX,STAT=IER)
          call CHECK_STAT(IER,980)
        end if
        if (allocated(JMAT)) then
          IR = IR + size(JMAT)
          deallocate (JMAT,STAT=IER)
          call CHECK_STAT(IER,990)
        end if
        if (allocated(LB)) then
          IR = IR + size(LB)
          deallocate (LB,STAT=IER)
          call CHECK_STAT(IER,1000)
        end if
        if (allocated(UB)) then
          IR = IR + size(UB)
          deallocate (UB,STAT=IER)
          call CHECK_STAT(IER,1000)
        end if
        if (allocated(PMAT)) then
          IR = IR + size(PMAT)
          deallocate (PMAT,STAT=IER)
          call CHECK_STAT(IER,1010)
        end if
        if (allocated(RWORK)) then
          IR = IR + size(RWORK)
          deallocate (RWORK,STAT=IER)
          call CHECK_STAT(IER,1020)
        end if
        if (allocated(SAVF)) then
          IR = IR + size(SAVF)
          deallocate (SAVF,STAT=IER)
          call CHECK_STAT(IER,1030)
        end if
        if (allocated(YMAX)) then
          IR = IR + size(YMAX)
          deallocate (YMAX,STAT=IER)
          call CHECK_STAT(IER,1040)
        end if
        if (allocated(WM)) then
          IR = IR + size(WM)
          deallocate (WM,STAT=IER)
          call CHECK_STAT(IER,1050)
        end if
        if (allocated(YHNQP2)) then
          IR = IR + size(YHNQP2)
          deallocate (YHNQP2,STAT=IER)
          call CHECK_STAT(IER,1060)
        end if
        if (allocated(YHTEMP)) then
          IR = IR + size(YHTEMP)
          deallocate (YHTEMP,STAT=IER)
          call CHECK_STAT(IER,1070)
        end if
        if (allocated(YMAX)) then
          IR = IR + size(YMAX)
          deallocate (YMAX,STAT=IER)
          call CHECK_STAT(IER,1080)
        end if
        if (allocated(YNNEG)) then
          IR = IR + size(YNNEG)
          deallocate (YNNEG,STAT=IER)
          call CHECK_STAT(IER,1090)
        end if
        if (allocated(YSCALEDS)) then
          IR = IR + size(YSCALEDS)
          deallocate (YSCALEDS,STAT=IER)
          call CHECK_STAT(IER,1100)
        end if
        if (allocated(YTEMP)) then
          IR = IR + size(YTEMP)
          deallocate (YTEMP,STAT=IER)
          call CHECK_STAT(IER,1110)
        end if
        if (allocated(WKDS)) then
          IR = IR + size(WKDS)
          deallocate (WKDS,STAT=IER)
          call CHECK_STAT(IER,1120)
        end if
        II = 0
        if (allocated(BIGP)) then
          II = II + size(BIGP)
          deallocate (BIGP,STAT=IER)
          call CHECK_STAT(IER,1130)
        end if
        if (allocated(BJGP)) then
          II = II + size(BJGP)
          deallocate (BJGP,STAT=IER)
          call CHECK_STAT(IER,1140)
        end if
        if (allocated(IA)) then
          II = II + size(IA)
          deallocate (IA,STAT=IER)
          call CHECK_STAT(IER,1150)
        end if
        if (allocated(IAB)) then
          II = II + size(IAB)
          deallocate (IAB,STAT=IER)
          call CHECK_STAT(IER,1160)
        end if
        if (allocated(IAN)) then
          II = II + size(IAN)
          deallocate (IAN,STAT=IER)
          call CHECK_STAT(IER,1170)
        end if
        if (allocated(ICN)) then
          II = II + size(ICN)
          deallocate (ICN,STAT=IER)
          call CHECK_STAT(IER,1180)
        end if
        if (allocated(IDX)) then
          II = II + size(IDX)
          deallocate (IDX,STAT=IER)
          call CHECK_STAT(IER,1190)
        end if
        if (allocated(IGP)) then
          II = II + size(IGP)
          deallocate (IGP,STAT=IER)
          call CHECK_STAT(IER,1200)
        end if
        if (allocated(IKEEP28)) then
          II = II + size(IKEEP28,1)*size(IKEEP28,2)
          deallocate (IKEEP28,STAT=IER)
          call CHECK_STAT(IER,1210)
        end if
        if (allocated(INDCOLDS)) then
          II = II + size(INDCOLDS)
          deallocate (INDCOLDS,STAT=IER)
          call CHECK_STAT(IER,1220)
        end if
        if (allocated(INDROWDS)) then
          II = II + size(INDROWDS)
          deallocate (INDROWDS,STAT=IER)
          call CHECK_STAT(IER,1230)
        end if
        if (allocated(IOPTDS)) then
          II = II + size(IOPTDS)
          deallocate (IOPTDS,STAT=IER)
          call CHECK_STAT(IER,1240)
        end if
        if (allocated(IPNTRDS)) then
          II = II + size(IPNTRDS)
          deallocate (IPNTRDS,STAT=IER)
          call CHECK_STAT(IER,1250)
        end if
        if (allocated(IWADS)) then
          II = II + size(IWADS)
          deallocate (IWADS,STAT=IER)
          call CHECK_STAT(IER,1260)
        end if
        if (allocated(IWKDS)) then
          II = II + size(IWKDS)
          deallocate (IWKDS,STAT=IER)
          call CHECK_STAT(IER,1270)
        end if
        if (allocated(IWORK)) then
          II = II + size(IWORK)
          deallocate (IWORK,STAT=IER)
          call CHECK_STAT(IER,1280)
        end if
        if (allocated(IW28)) then
          II = II + size(IW28,1)*size(IW28,2)
          deallocate (IW28,STAT=IER)
          call CHECK_STAT(IER,1290)
        end if
        if (allocated(JA)) then
          II = II + size(JA)
          deallocate (JA,STAT=IER)
          call CHECK_STAT(IER,1300)
        end if
        if (allocated(JAB)) then
          II = II + size(JAB)
          deallocate (JAB,STAT=IER)
          call CHECK_STAT(IER,1310)
        end if
        if (allocated(JAN)) then
          II = II + size(JAN)
          deallocate (JAN,STAT=IER)
          call CHECK_STAT(IER,1320)
        end if
        if (allocated(JATEMP)) then
          II = II + size(JATEMP)
          deallocate (JATEMP,STAT=IER)
          call CHECK_STAT(IER,1330)
        end if
        if (allocated(JGP)) then
          II = II + size(JGP)
          deallocate (JGP,STAT=IER)
          call CHECK_STAT(IER,1340)
        end if
        if (allocated(JPNTRDS)) then
          II = II + size(JPNTRDS)
          deallocate (JPNTRDS,STAT=IER)
          call CHECK_STAT(IER,1350)
        end if
        if (allocated(JROOT)) then
          II = II + size(JROOT)
          deallocate (JROOT,STAT=IER)
          call CHECK_STAT(IER,1360)
        end if
        if (allocated(JVECT)) then
          II = II + size(JVECT)
          deallocate (JVECT,STAT=IER)
          call CHECK_STAT(IER,1370)
        end if
        if (allocated(NGRPDS)) then
          II = II + size(NGRPDS)
          deallocate (NGRPDS,STAT=IER)
          call CHECK_STAT(IER,1380)
        end if
        if (allocated(SUBDS)) then
          II = II + size(SUBDS)
          deallocate (SUBDS,STAT=IER)
          call CHECK_STAT(IER,1390)
        end if
        if (allocated(SUPDS)) then
          II = II + size(SUPDS)
          deallocate (SUPDS,STAT=IER)
          call CHECK_STAT(IER,1400)
        end if
!_______________________________________________________________________
! *****MA48 build change point. Insert these statements.
!       IF (MA48_WAS_USED) THEN
!          CALL MA48_FINALIZE(FACTORS,CONTROL,INFO)
!          IF (INFO /= 0) THEN
!             MSG = 'The call to MA48_FINALIZE FAILED.'
!             CALL XERRDV(MSG,1740,1,1,II,0,0,ZERO,ZERO)
!          END IF
!          MSG = 'Size of MA48 deallocated arrays (I1) = '
!          CALL XERRDV(MSG,1750,1,1,ISIZE,0,0,ZERO,ZERO)
!       END IF
!_______________________________________________________________________

!       Print the amount of storage used.
        MSG = 'I1 = Total length of REAL arrays used.'
        call XERRDV(MSG,1760,1,1,IR,0,0,ZERO,ZERO)
        MSG = 'I1 = Total length of INTEGER arrays used.'
        call XERRDV(MSG,1760,1,1,II,0,0,ZERO,ZERO)

!       In case DVODE_F90 is subsequently called:
        OPTS_CALLED = .false.
        return

      end subroutine RELEASE_ARRAYS
! End of DVODE_F90 subroutines
!_______________________________________________________________________

! Beginning of LINPACK and BLAS subroutines
      subroutine DGEFA_F90(A,LDA,N,IPVT,INFO)
! ..
! Factor a matrix using Gaussian elimination.
! ..
!     DGEFA_F90 factors a real(wp) matrix by Gaussian elimination.
!     DGEFA_F90 is usually called by DGECO, but it can be called
!     directly with a saving in time if RCOND is not needed.
!     (Time for DGECO) = (1 + 9/N)*(Time for DGEFA_F90).
!     On Entry
!        A       REAL(KIND=WP)(LDA, N)
!                the matrix to be factored.
!        LDA     INTEGER
!                the leading dimension of the array A.
!        N       INTEGER
!                the order of the matrix A.
!     On Return
!        A       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                The factorization can be written A = L*U where
!                L is a product of permutation and unit lower
!                triangular matrices and U is upper triangular.
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!        INFO    INTEGER
!                = 0  normal value.
!                = K  if U(K,K) == 0.0. This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DGESL_F90 or DGEDI will divide
!                     by zero if called. Use RCOND in DGECO for a
!                     reliable indication of singularity.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (INOUT) :: INFO
        integer, intent (IN) :: LDA, N
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: A(LDA,*)
        integer, intent (INOUT) :: IPVT(*)
! ..
! .. Local Scalars ..
        real (WP) :: T
        integer :: J, K, KP1, L, NM1
! ..
! .. Intrinsic Functions ..
        intrinsic ABS
! ..
! .. FIRST EXECUTABLE STATEMENT DGEFA_F90
! ..
        INFO = 0
        NM1 = N - 1
        if (NM1<1) goto 50
        do K = 1, NM1
          KP1 = K + 1

!         Find L = pivot index.

! Original:
!         L = IDAMAX_F90(N-K+1,A(K,K),1) + K - 1
          L = IDAMAX_F90(N-K+1,A(K:N,K),1) + K - 1
          IPVT(K) = L

!         Zero pivot implies this column already triangularized.

!         IF (A(L, K) == ZERO) GOTO 40
          if (abs(A(L,K))<=ZERO) goto 30

!         Interchange if necessary.

          if (L==K) goto 10
          T = A(L,K)
          A(L,K) = A(K,K)
          A(K,K) = T
10        continue

!         Compute multipliers.

          T = -ONE/A(K,K)
! Original:          
!         CALL DSCAL_F90(N-K,T,A(K+1,K),1)         
          call DSCAL_F90(N-K,T,A(K+1:N,K),1)

!         Row elimination with column indexing.

          do J = KP1, N
            T = A(L,J)
            if (L==K) goto 20
            A(L,J) = A(K,J)
            A(K,J) = T
20          continue
! Original:            
!           CALL DAXPY_F90(N-K,T,A(K+1,K),1,A(K+1,J),1)
            call DAXPY_F90(N-K,T,A(K+1:N,K),1,A(K+1:N,J),1)
          end do
          goto 40
30        continue
          INFO = K
40        continue
        end do
50      continue
        IPVT(N) = N
!       IF (A(N, N) == ZERO) INFO = N
        if (abs(A(N,N))<=ZERO) INFO = N
        return

      end subroutine DGEFA_F90
!_______________________________________________________________________

      subroutine DGESL_F90(A,LDA,N,IPVT,B,JOB)
! ..
! Solve the real system A*X=B or TRANS(A)*X=B using the factors
! computed by DGECO or DGEFA_F90.
! ..
!     DGESL_F90 solves the real(wp) system
!     A * X = B or TRANS(A) * X = B
!     using the factors computed by DGECO or DGEFA_F90.
!     On Entry
!        A       REAL(KIND=WP)(LDA, N)
!                the output from DGECO or DGEFA_F90.
!        LDA     INTEGER
!                the leading dimension of the array A.
!        N       INTEGER
!                the order of the matrix A.
!        IPVT    INTEGER(N)
!                the pivot vector from DGECO or DGEFA_F90.
!        B       REAL(KIND=WP)(N)
!                the right hand side vector.
!        JOB     INTEGER
!                = 0         to solve A*X = B,
!                = nonzero   to solve TRANS(A)*X = B where
!                            TRANS(A)is the transpose.
!     On Return
!        B       the solution vector X.
!     Error Condition
!        A division by zero will occur if the input factor contains a
!        zero on the diagonal. Technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of LDA. It will not occur if the subroutines are
!        called correctly and if DGECO has set RCOND > 0.0 or
!        DGEFA_F90 has set INFO == 0.
!     To compute INVERSE(A) * C where C is a matrix
!     with P columns
!           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
!           IF (RCOND is too small) GOTO
!           DO J = 1, P
!              CALL DGESL_F90(A,LDA,N,IPVT,C(1,J),0)
!           END DO
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: JOB, LDA, N
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: A(LDA,*), B(*)
        integer, intent (INOUT) :: IPVT(*)
! ..
! .. Local Scalars ..
        real (WP) :: T
        integer :: K, KB, L, NM1
! ..
! .. FIRST EXECUTABLE STATEMENT DGESL_F90
! ..
        NM1 = N - 1
        if (JOB/=0) goto 30

!       JOB = 0, solve A*X = B.
!       First solve L*Y = B.

        if (NM1<1) goto 20
        do K = 1, NM1
          L = IPVT(K)
          T = B(L)
          if (L==K) goto 10
          B(L) = B(K)
          B(K) = T
10        continue
! Original:
!         CALL DAXPY_F90(N-K,T,A(K+1,K),1,B(K+1),1)
          call DAXPY_F90(N-K,T,A(K+1:N,K),1,B(K+1:N),1)
        end do
20      continue

!       Now solve U*X = Y.

        do KB = 1, N
          K = N + 1 - KB
          B(K) = B(K)/A(K,K)
          T = -B(K)
! Original:
!         CALL DAXPY_F90(K-1,T,A(1,K),1,B(1),1)
          call DAXPY_F90(K-1,T,A(1:K-1,K),1,B(1:K-1),1)
        end do
        goto 60
30      continue

!       JOB /= 0, solve TRANS(A)*X = B.
!       First solve TRANS(U)*Y = B.

        do K = 1, N
          T = DDOT_F90(K-1,A(1,K),1,B(1),1)
          B(K) = (B(K)-T)/A(K,K)
        end do

!       Now solve TRANS(L)*X = Y.

        if (NM1<1) goto 50
        do KB = 1, NM1
          K = N - KB
          B(K) = B(K) + DDOT_F90(N-K,A(K+1,K),1,B(K+1),1)
          L = IPVT(K)
          if (L==K) goto 40
          T = B(L)
          B(L) = B(K)
          B(K) = T
40        continue
        end do
50      continue
60      continue
        return

      end subroutine DGESL_F90
!_______________________________________________________________________

      subroutine DGBFA_F90(ABD,LDA,N,ML,MU,IPVT,INFO)
! ..
! Factor a banded matrix using Gaussian elimination.
! ..
!     DGBFA_F90 factors a real(wp) band matrix by elimination.
!     DGBFA_F90 is usually called by DGBCO, but it can be called
!     directly with a saving in time if RCOND is not needed.
!     On Entry
!        ABD     REAL(KIND=WP)(LDA, N)
!                contains the matrix in band storage. The columns
!                of the matrix are stored in the columns of ABD and
!                the diagonals of the matrix are stored in rows
!                ML+1 through 2*ML+MU+1 of ABD.
!                See the comments below for details.
!        LDA     INTEGER
!                the leading dimension of the array ABD.
!                LDA must be >= 2*ML + MU + 1.
!        N       INTEGER
!                the order of the original matrix.
!        ML      INTEGER
!                number of diagonals below the main diagonal.
!                0 <= ML < N.
!        MU      INTEGER
!                number of diagonals above the main diagonal.
!                0 <= MU < N.
!                More efficient if ML <= MU.
!     On Return
!        ABD     an upper triangular matrix in band storage and
!                the multipliers which were used to obtain it.
!                The factorization can be written A = L*U where
!                L is a product of permutation and unit lower
!                triangular matrices and U is upper triangular.
!        IPVT    INTEGER(N)
!                an integer vector of pivot indices.
!        INFO    INTEGER
!                = 0  normal value.
!                = K  if U(K,K) == 0.0. This is not an error
!                     condition for this subroutine, but it does
!                     indicate that DGBSL_F90 will divide by zero
!                     if called. Use RCOND in DGBCO for a reliable
!                     indication of singularity.
!     Band Storage
!           If A is a band matrix, the following program segment
!           will set up the input.
!                   ML = (band width below the diagonal)
!                   MU = (band width above the diagonal)
!                   M = ML + MU + 1
!                   DO J = 1, N
!                      I1 = MAX(1, J-MU)
!                      I2 = MIN(N, J+ML)
!                      DO 10 I = I1, I2
!                         K = I - J + M
!                         ABD(K,J) = A(I,J)
!                      END DO
!                   END DO
!           This uses rows ML+1 through 2*ML+MU+1 of ABD.
!           In addition, the first ML rows in ABD are used for
!           elements generated during the triangularization.
!           The total number of rows needed in ABD is 2*ML+MU+1.
!           The ML+MU by ML+MU upper left triangle and the
!           ML by ML lower right triangle are not referenced.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (INOUT) :: INFO
        integer, intent (IN) :: LDA, ML, MU, N
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: ABD(LDA,*)
        integer, intent (INOUT) :: IPVT(*)
! ..
! .. Local Scalars ..
        real (WP) :: T
        integer :: I0, J, J0, J1, JU, JZ, K, KP1, L, LM, M, MM, NM1
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, MAX, MIN
! ..
! .. FIRST EXECUTABLE STATEMENT DGBFA_F90
! ..
        M = ML + MU + 1
        INFO = 0

!       Zero initial fill-in columns.

        J0 = MU + 2
        J1 = min(N,M) - 1
        if (J1<J0) goto 10
        do JZ = J0, J1
          I0 = M + 1 - JZ
          ABD(I0:ML,JZ) = ZERO
        end do
10      continue
        JZ = J1
        JU = 0

!       Gaussian elimination with partial pivoting.

        NM1 = N - 1
        if (NM1<1) goto 80
        do K = 1, NM1
          KP1 = K + 1

!         Zero next fill-in column.

          JZ = JZ + 1
          if (JZ>N) goto 20
          if (ML<1) goto 20
          ABD(1:ML,JZ) = ZERO
20        continue

!         Find L = pivot index.

          LM = min(ML,N-K)
! Original:
!         L = IDAMAX_F90(LM+1,ABD(M,K),1) + M - 1
          L = IDAMAX_F90(LM+1,ABD(M:M+LM,K),1) + M - 1

          IPVT(K) = L + K - M

!         Zero pivot implies this column already triangularized.

!         IF (ABD(L, K) == ZERO) GOTO 100
          if (abs(ABD(L,K))<=ZERO) goto 60

!         Interchange if necessary.

          if (L==M) goto 30
          T = ABD(L,K)
          ABD(L,K) = ABD(M,K)
          ABD(M,K) = T
30        continue

!         Compute multipliers.

          T = -ONE/ABD(M,K)
! Original:
!         CALL DSCAL_F90(LM,T,ABD(M+1,K),1)
          call DSCAL_F90(LM,T,ABD(M+1:M+LM,K),1)

!         Row elimination with column indexing.

          JU = min(max(JU,MU+IPVT(K)),N)
          MM = M
          if (JU<KP1) goto 50
          do J = KP1, JU
            L = L - 1
            MM = MM - 1
            T = ABD(L,J)
            if (L==MM) goto 40
            ABD(L,J) = ABD(MM,J)
            ABD(MM,J) = T
40          continue
! Original:
!           CALL DAXPY_F90(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
            call DAXPY_F90(LM,T,ABD(M+1:M+LM,K),1,ABD(MM+1:M+LM,J),1)
          end do
50        continue
          goto 70
60        continue
          INFO = K
70        continue
        end do
80      continue
        IPVT(N) = N
!       IF (ABD(M, N) == ZERO) INFO = N
        if (abs(ABD(M,N))<=ZERO) INFO = N
        return

      end subroutine DGBFA_F90
!_______________________________________________________________________

      subroutine DGBSL_F90(ABD,LDA,N,ML,MU,IPVT,B,JOB)
! ..
! Solve the real band system A*X=B or TRANS(A)*X=B using the factors
! computed by DGBCO or DGBFA_F90.
! ..
!     DGBSL_F90 solves the real(wp) band system
!     A * X = B or TRANS(A) * X = B
!     using the factors computed by DGBCO or DGBFA_F90.
!     On Entry
!        ABD     REAL(KIND=WP)(LDA, N)
!                the output from DGBCO or DGBFA_F90.
!        LDA     INTEGER
!                the leading dimension of the array ABD.
!        N       INTEGER
!                the order of the original matrix.
!        ML      INTEGER
!                number of diagonals below the main diagonal.
!        MU      INTEGER
!                number of diagonals above the main diagonal.
!        IPVT    INTEGER(N)
!                the pivot vector from DGBCO or DGBFA_F90.
!        B       REAL(KIND=WP)(N)
!                the right hand side vector.
!        JOB     INTEGER
!                = 0         to solve A*X = B,
!                = nonzero   to solve TRANS(A)*X = B, where
!                            TRANS(A) is the transpose.
!     On Return
!        B       the solution vector X.
!     Error Condition
!        A division by zero will occur if the input factor contains a
!        zero on the diagonal. Technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of LDA. It will not occur if the subroutines are
!        called correctly and if DGBCO has set RCOND > 0.0
!        or DGBFA_F90 has set INFO == 0 .
!     To compute INVERSE(A) * C where C is a matrix
!     with P columns
!           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
!           IF (RCOND is too small) GOTO ...
!           DO J = 1, P
!              CALL DGBSL_F90(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
!           END DO
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: JOB, LDA, ML, MU, N
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: ABD(LDA,*), B(*)
        integer, intent (INOUT) :: IPVT(*)
! ..
! .. Local Scalars ..
        real (WP) :: T
        integer :: K, KB, L, LA, LB, LM, M, NM1
! ..
! .. Intrinsic Functions ..
        intrinsic MIN
! ..
! .. FIRST EXECUTABLE STATEMENT DGBSL_F90
! ..
        M = MU + ML + 1
        NM1 = N - 1
        if (JOB/=0) goto 30

!       JOB = 0, solve A*X = B.
!       First solve L*Y = B.

        if (ML==0) goto 20
        if (NM1<1) goto 20
        do K = 1, NM1
          LM = min(ML,N-K)
          L = IPVT(K)
          T = B(L)
          if (L==K) goto 10
          B(L) = B(K)
          B(K) = T
10        continue
! Original:
!         CALL DAXPY_F90(LM,T,ABD(M+1,K),1,B(K+1),1)
          call DAXPY_F90(LM,T,ABD(M+1:M+LM,K),1,B(K+1:K+LM),1)
        end do
20      continue

!       Now solve U*X = Y.

        do KB = 1, N
          K = N + 1 - KB
          B(K) = B(K)/ABD(M,K)
          LM = min(K,M) - 1
          LA = M - LM
          LB = K - LM
          T = -B(K)
! Original:
!         CALL DAXPY_F90(LM,T,ABD(LA,K),1,B(LB),1)
          call DAXPY_F90(LM,T,ABD(LA:LA+LM-1,K),1,B(LB:LB+LM-1),1)
        end do
        goto 60
30      continue

!       JOB /= 0, solve TRANS(A)*X = B.
!       First solve TRANS(U)*Y = B.

        do K = 1, N
          LM = min(K,M) - 1
          LA = M - LM
          LB = K - LM
! Original:
!         T = DDOT_F90(LM,ABD(LA,K),1,B(LB),1)
          T = DDOT_F90(LM,ABD(LA:LA+LM-1,K),1,B(LB:LB+LM-1),1)
          B(K) = (B(K)-T)/ABD(M,K)
        end do

!       Now solve TRANS(L)*X = Y.

        if (ML==0) goto 50
        if (NM1<1) goto 50
        do KB = 1, NM1
          K = N - KB
          LM = min(ML,N-K)
! Original:
!         B(K) = B(K) + DDOT_F90(LM,ABD(M+1,K),1,B(K+1),1)
          B(K) = B(K) + DDOT_F90(LM,ABD(M+1:M+LM,K),1,B(K+1:K+LM),1)
          L = IPVT(K)
          if (L==K) goto 40
          T = B(L)
          B(L) = B(K)
          B(K) = T
40        continue
        end do
50      continue
60      continue
        return

      end subroutine DGBSL_F90
!_______________________________________________________________________

      subroutine DAXPY_F90(N,DA,DX,INCX,DY,INCY)
! ..
! Compute a constant times a vector plus a vector.
! ..
!     Description of Parameters
!     Input:
!     N    - number of elements in input vector(s)
!     DA   - real(wp) scalar multiplier
!     DX   - real(wp) vector with N elements
!     INCX - storage spacing between elements of DX
!     DY   - real(wp) vector with N elements
!     INCY -  storage spacing between elements of DY
!     Output:
!     DY - real(wp) result (unchanged if N <= 0)
!     Overwrite real(wp) DY with real(wp) DA*DX + DY.
!     For I = 0 to N-1, replace  DY(LY+I*INCY) with
!     DA*DX(LX+I*INCX) + DY(LY+I*INCY),
!     where LX = 1 if INCX >= 0, else LX = 1+(1-N)*INCX,
!     and LY is defined in a similar way using INCY.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP), intent (IN) :: DA
        integer, intent (IN) :: INCX, INCY, N
! ..
! .. Array Arguments ..
        real (WP), intent (IN) :: DX(*)
        real (WP), intent (INOUT) :: DY(*)
! ..
! .. Local Scalars ..
        integer :: I, IX, IY, M, MP1, NS
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, MOD
! ..
! .. FIRST EXECUTABLE STATEMENT DAXPY_F90
! ..
!       IF (N <= 0 .OR. DA == ZERO) RETURN
        if (N<=0 .or. abs(DA)<=ZERO) return
!       IF (INCX==INCY) IF (INCX-1) 10, 20, 40
        if (INCX == INCY) then
          if (INCX < 1) then
             goto 10
          elseif (INCX == 1) then
             goto 20
          else
             goto 40
          end if
        end if

!       Code for unequal or nonpositive increments.

10      IX = 1
        IY = 1
        if (INCX<0) IX = (-N+1)*INCX + 1
        if (INCY<0) IY = (-N+1)*INCY + 1
        do I = 1, N
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
        end do
        return

!       Code for both increments equal to 1.

!       Clean-up loop so remaining vector length is a multiple of 4.

20      M = mod(N,4)
        if (M==0) goto 30
        DY(1:M) = DY(1:M) + DA*DX(1:M)
        if (N<4) return
30      MP1 = M + 1
        do I = MP1, N, 4
          DY(I) = DY(I) + DA*DX(I)
          DY(I+1) = DY(I+1) + DA*DX(I+1)
          DY(I+2) = DY(I+2) + DA*DX(I+2)
          DY(I+3) = DY(I+3) + DA*DX(I+3)
        end do
        return

!       Code for equal, positive, non-unit increments.

40      NS = N*INCX
        do I = 1, NS, INCX
          DY(I) = DA*DX(I) + DY(I)
        end do
        return

      end subroutine DAXPY_F90
!_______________________________________________________________________

      subroutine DCOPY_F90(N,DX,INCX,DY,INCY)
! ..
! Copy a vector to another vector.
! ..
!     Description of Parameters
!     Input:
!     N    - number of elements in input vector(s)
!     DX   - real(wp) vector with N elements
!     INCX - storage spacing between elements of DX
!     DY   - real(wp) vector with N elements
!     INCY - storage spacing between elements of DY
!     Output:
!     DY - copy of vector DX(unchanged if N <= 0)
!     Copy real(wp) DX to real(wp) DY.
!     For I = 0 to N-1, copy DX(LX+I*INCX) to DY(LY+I*INCY),
!     where LX = 1 if INCX >= 0, else LX = 1+(1-N)*INCX,
!     and LY is defined in a similar way using INCY.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: INCX, INCY, N
! ..
! .. Array Arguments ..
        real (WP), intent (IN) :: DX(*)
        real (WP), intent (INOUT) :: DY(*)
! ..
! .. Local Scalars ..
        integer :: I, IX, IY, M, MP1, NS
! ..
! .. Intrinsic Functions ..
        intrinsic MOD
! ..
! .. FIRST EXECUTABLE STATEMENT DCOPY_F90
! ..
        if (N<=0) return
!       IF (INCX==INCY) IF (INCX-1) 10, 20, 40
        if (INCX == INCY) then
          if (INCX < 1) then
             goto 10
          elseif (INCX == 1) then
             goto 20
          else
             goto 40
          end if
        end if

!       Code for unequal or nonpositive increments.

10      IX = 1
        IY = 1
        if (INCX<0) IX = (-N+1)*INCX + 1
        if (INCY<0) IY = (-N+1)*INCY + 1
        do I = 1, N
          DY(IY) = DX(IX)
          IX = IX + INCX
          IY = IY + INCY
        end do
        return

!       Code for both increments equal to 1.

!       Clean-up loop so remaining vector length is a multiple of 7.

20      M = mod(N,7)
        if (M==0) goto 30
        do I = 1, M
          DY(I) = DX(I)
        end do
        if (N<7) return
30      MP1 = M + 1
        do I = MP1, N, 7
          DY(I) = DX(I)
          DY(I+1) = DX(I+1)
          DY(I+2) = DX(I+2)
          DY(I+3) = DX(I+3)
          DY(I+4) = DX(I+4)
          DY(I+5) = DX(I+5)
          DY(I+6) = DX(I+6)
        end do
        return

!       Code for equal, positive, non-unit increments.

40      NS = N*INCX
        do I = 1, NS, INCX
          DY(I) = DX(I)
        end do
        return

      end subroutine DCOPY_F90
!_______________________________________________________________________

      function DDOT_F90(N,DX,INCX,DY,INCY)
! ..
! Compute the inner product of two vectors.
! ..
!     Description of Parameters
!     Input:
!     N    - number of elements in input vector(s)
!     DX   - real(wp) vector with N elements
!     INCX - storage spacing between elements of DX
!     DY   - real(wp) vector with N elements
!     INCY - storage spacing between elements of DY
!     Output:
!     DDOT_F90 - real(wp) dot product (zero if N <= 0)
!     Returns the dot product of real(wp) DX and DY.
!     DDOT_F90 = sum for I = 0 to N-1 of  DX(LX+I*INCX) * DY(LY+I*INCY),
!     where LX = 1 if INCX >= 0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
! ..
     implicit none
! ..
! .. Function Return Value ..
        real (WP) :: DDOT_F90
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: INCX, INCY, N
! ..
! .. Array Arguments ..
        real (WP), intent (IN) :: DX(*), DY(*)
! ..
! .. Local Scalars ..
        integer :: I, IX, IY, M, MP1, NS
! ..
! .. Intrinsic Functions ..
        intrinsic MOD
! ..
! .. FIRST EXECUTABLE STATEMENT DDOT_F90
! ..
        DDOT_F90 = ZERO
        if (N<=0) return
!       IF (INCX==INCY) IF (INCX-1) 10, 20, 40
        if (INCX == INCY) then
          if (INCX < 1) then
             goto 10
          elseif (INCX == 1) then
             goto 20
          else
             goto 40
          end if
        end if

!       Code for unequal or nonpositive increments.

10      IX = 1
        IY = 1
        if (INCX<0) IX = (-N+1)*INCX + 1
        if (INCY<0) IY = (-N+1)*INCY + 1
        do I = 1, N
          DDOT_F90 = DDOT_F90 + DX(IX)*DY(IY)
          IX = IX + INCX
          IY = IY + INCY
        end do
        return

!       Code for both increments equal to 1.

!       Clean-up loop so remaining vector length is a multiple of 5.

20      M = mod(N,5)
        if (M==0) goto 30
        do I = 1, M
          DDOT_F90 = DDOT_F90 + DX(I)*DY(I)
        end do
        if (N<5) return
30      MP1 = M + 1
        do I = MP1, N, 5
          DDOT_F90 = DDOT_F90 + DX(I)*DY(I) + DX(I+1)*DY(I+1) + &
            DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
        end do
        return

!       Code for equal, positive, non-unit increments.

40      NS = N*INCX
        do I = 1, NS, INCX
          DDOT_F90 = DDOT_F90 + DX(I)*DY(I)
        end do
        return

      end function DDOT_F90
!_______________________________________________________________________

      subroutine DSCAL_F90(N,DA,DX,INCX)
! ..
! Multiply a vector by a constant.
! ..
!     Description of Parameters
!     Input:
!     N    - number of elements in input vector(s)
!     DA   - real(wp) scale factor
!     DX   - real(wp) vector with N elements
!     INCX - storage spacing between elements of DX
!     Output:
!     DX - real(wp) result (unchanged if N <= 0)
!     Replace real(wp) DX by real(wp) DA*DX.
!     For I = 0 to N-1, replace DX(IX+I*INCX) with DA * DX(IX+I*INCX),
!     where IX = 1 if INCX >= 0, else IX = 1+(1-N)*INCX.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP), intent (IN) :: DA
        integer, intent (IN) :: INCX, N
! ..
! .. Array Arguments ..
        real (WP), intent (INOUT) :: DX(*)
! ..
! .. Local Scalars ..
        integer :: I, IX, M, MP1
! ..
! .. Intrinsic Functions ..
        intrinsic MOD
! ..
! .. FIRST EXECUTABLE STATEMENT DSCAL_F90
! ..
        if (N<=0) return
        if (INCX==1) goto 10

!       Code for increment not equal to 1.

        IX = 1
        if (INCX<0) IX = (-N+1)*INCX + 1
        do I = 1, N
          DX(IX) = DA*DX(IX)
          IX = IX + INCX
        end do
        return

!       Code for increment equal to 1.

!       Clean-up loop so remaining vector length is a multiple of 5.

10      M = mod(N,5)
        if (M==0) goto 20
        DX(1:M) = DA*DX(1:M)
        if (N<5) return
20      MP1 = M + 1
        do I = MP1, N, 5
          DX(I) = DA*DX(I)
          DX(I+1) = DA*DX(I+1)
          DX(I+2) = DA*DX(I+2)
          DX(I+3) = DA*DX(I+3)
          DX(I+4) = DA*DX(I+4)
        end do
        return

      end subroutine DSCAL_F90
!_______________________________________________________________________

      function IDAMAX_F90(N,DX,INCX)
! ..
! Find the smallest index of that component of a vector
! having the maximum magnitude.
! ..
!     Description of Parameters
!     Input:
!     N    - number of elements in input vector(s)
!     DX   - real(wp) vector with N elements
!     INCX - storage spacing between elements of DX
!     Output:
!     IDAMAX_F90 - smallest index (zero if N <= 0)
!     Find smallest index of maximum magnitude of real(wp) DX.
!     IDAMAX_F90 = first I, I = 1 to N, to maximize
!     ABS(DX(IX+(I-1)*INCX)), where IX = 1 if INCX >= 0,
!     else IX = 1+(1-N)*INCX.
! ..
     implicit none
! ..
! .. Function Return Value ..
        integer :: IDAMAX_F90
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: INCX, N
! ..
! .. Array Arguments ..
        real (WP), intent (IN) :: DX(*)
! ..
! .. Local Scalars ..
        real (WP) :: DMAX, XMAG
        integer :: I, IX
! ..
! .. Intrinsic Functions ..
        intrinsic ABS
! ..
! .. FIRST EXECUTABLE STATEMENT IDAMAX_F90
! ..
        IDAMAX_F90 = 0
        if (N<=0) return
        IDAMAX_F90 = 1
        if (N==1) return

        if (INCX==1) goto 10

!       Code for increments not equal to 1.

        IX = 1
        if (INCX<0) IX = (-N+1)*INCX + 1
        DMAX = abs(DX(IX))
        IX = IX + INCX
        do I = 2, N
          XMAG = abs(DX(IX))
          if (XMAG>DMAX) then
            IDAMAX_F90 = I
            DMAX = XMAG
          end if
          IX = IX + INCX
        end do
        return

!       Code for increments equal to 1.

10      DMAX = abs(DX(1))
        do I = 2, N
          XMAG = abs(DX(I))
          if (XMAG>DMAX) then
            IDAMAX_F90 = I
            DMAX = XMAG
          end if
        end do
        return

      end function IDAMAX_F90
! End of LINPACK and BLAS subroutines
!_______________________________________________________________________

! Beginning of MA28 subroutines
! THIS IS A FORTRAN 90 TRANSLATION OF HSL'S F77 MA28. IT IS INTENDED
! ONLY FOR USE IN CONJUNCTION WITH THE ODE SOLVER DVODE_F90 AND IS
! NOT FUNCTIONAL IN A STANDALONE MANNER. PLEASE NOTE THAT MA28 IS NOT
! PUBLIC DOMAIN SOFTWARE BUT HAS BEEN MADE AVAILABLE TO THE NUMERICAL
! ANALYSIS COMMUNITY BY HARWELL. FOR OTHER USE PLEASE CONTACT HARWELL
! AT HTTP://WWW.HSL-LIBRARY.COM OR CONTACT HSL@HYPROTECH.COM IF YOU
! WISH TO SOLVE GENERAL SPARSE LINEAR SYSTEMS. IF YOU FIND A BUG OR
! ENCOUNTER A PROBLEM WITH THE USE OF MA28 WITH DVODE.F90, PLEASE
! CONTACT ONE OF THE AUTHORS OF DVODE_F90:
!           G.D. Byrne (gbyrne@wi.rr.com)
!           S. Thompson, thompson@radford.edu
! NUMEROUS CHANGES WERE MADE IN CONNECTION WITH DVODE_F90 USAGE. THESE
! INCLUDE USING METCALF'S CONVERTER TO TRANSLATE THE ORIGINAL f77 CODE
! TO F90, MOVING INITIALIZATIONS TO THE DVODE_F90 PRIVATE SECTION,
! ELIMINATION OF HOLLERITHS IN FORMAT STATEMENTS, ELIMINATION OF
! BLOCKDATA, CHANGES IN ARITHMETICAL OPERATOR SYNTAX, AND CONVERSION
! TO UPPER CASE. THIS VERSION OF MA28 IS INTENDED ONLY FOR USE WITH
! DVODE_F90. PLEASE DO NOT MODIFY IT FOR ANY OTHER PURPOSE. IF YOU
! HAVE LICENSED ACCESS TO THE HSL LIBRARY, AN ALTERNATE VERSION OF
! DVODE_F90 BASED ON THE SUCCESSOR TO MA28, MA48, IS AVAILABLE FROM
! THE AUTHORS. PLEASE NOTE THAT THE ALTERNATE VERSION OF DVODE_F90
! IS NOT SELF CONTAINED SINCE MA48 IS NOT DISTRIBUTED WITH DVODE_F90.
!******************************************************************
!             *****MA28 COPYRIGHT NOTICE*****
! COPYRIGHT (C) 2001 COUNCIL FOR THE CENTRAL LABORATORY
!               OF THE RESEARCH COUNCILS
! ALL RIGHTS RESERVED.

! NONE OF THE COMMENTS IN THIS COPYRIGHT NOTICE BETWEEN THE LINES
! OF ASTERISKS SHALL BE REMOVED OR ALTERED IN ANY WAY.

! THIS PACKAGE IS INTENDED FOR COMPILATION WITHOUT MODIFICATION,
! SO MOST OF THE EMBEDDED COMMENTS HAVE BEEN REMOVED.

! ALL USE IS SUBJECT TO LICENCE. IF YOU NEED FURTHER CLARIFICATION,
! PLEASE SEE HTTP://WWW.HSL-LIBRARY.COM OR CONTACT HSL@HYPROTECH.COM

! PLEASE NOTE THAT:

! 1. THE PACKAGES MAY ONLY BE USED FOR THE PURPOSES SPECIFIED IN THE
!    LICENCE AGREEMENT AND MUST NOT BE COPIED BY THE LICENSEE FOR
!    USE BY ANY OTHER PERSONS. USE OF THE PACKAGES IN ANY COMMERCIAL
!    APPLICATION SHALL BE SUBJECT TO PRIOR WRITTEN AGREEMENT BETWEEN
!    HYPROTECH UK LIMITED AND THE LICENSEE ON SUITABLE TERMS AND
!    CONDITIONS, WHICH WILL INCLUDE FINANCIAL CONDITIONS.
! 2. ALL INFORMATION ON THE PACKAGE IS PROVIDED TO THE LICENSEE ON
!    THE UNDERSTANDING THAT THE DETAILS THEREOF ARE CONFIDENTIAL.
! 3. ALL PUBLICATIONS ISSUED BY THE LICENSEE THAT INCLUDE RESULTS
!    OBTAINED WITH THE HELP OF ONE OR MORE OF THE PACKAGES SHALL
!    ACKNOWLEDGE THE USE OF THE PACKAGES. THE LICENSEE WILL NOTIFY
!    HSL@HYPROTECH.COM OR HYPROTECH UK LIMITED OF ANY SUCH PUBLICATION.
! 4. THE PACKAGES MAY BE MODIFIED BY OR ON BEHALF OF THE LICENSEE
!    FOR SUCH USE IN RESEARCH APPLICATIONS BUT AT NO TIME SHALL SUCH
!    PACKAGES OR MODIFICATIONS THEREOF BECOME THE PROPERTY OF THE
!    LICENSEE. THE LICENSEE SHALL MAKE AVAILABLE FREE OF CHARGE TO THE
!    COPYRIGHT HOLDER FOR ANY PURPOSE ALL INFORMATION RELATING TO
!    ANY MODIFICATION.
! 5. NEITHER COUNCIL FOR THE CENTRAL LABORATORY OF THE RESEARCH
!    COUNCILS NOR HYPROTECH UK LIMITED SHALL BE LIABLE FOR ANY
!    DIRECT OR CONSEQUENTIAL LOSS OR DAMAGE WHATSOEVER ARISING OUT OF
!    THE USE OF PACKAGES BY THE LICENSEE.
!******************************************************************

      subroutine MA28AD(N,NZ,A,LICN,IRN,LIRN,ICN,U,IKEEP,IW,W,IFLAG)
! ..
! This subroutine performs the LU factorization of A.
! ..
! The parameters are as follows:
! N     Order of matrix. Not altered by subroutine.
! NZ    Number of non-zeros in input matrix. Not altered by subroutine.
! A     Array of length LICN. Holds non-zeros of matrix
!       on entry and non-zeros of factors on exit. Reordered by
!       MA20AD and MC23AD and altered by MA30AD.
! LICN  Length of arrays A and ICN. Not altered by subroutine.
! IRN   Array of length LIRN. Holds row indices on input.
!       Used as workspace by MA30AD to hold column orientation of
!       matrix.
! LIRN  Length of array IRN. Not altered by the subroutine.
! ICN   Array of length LICN. Holds column indices on entry
!       and column indices of decomposed matrix on exit. Reordered
!       by MA20AD and MC23AD and altered by MA30AD.
! U     Variable set by user to control bias towards numeric or
!       sparsity pivoting. U = 1.0 gives partial pivoting
!       while U = 0. does not check multipliers at all. Values of U
!       greater than one are treated as one while negative values
!       are treated as zero. Not altered by subroutine.
! IKEEP Array of length 5*N used as workspace by MA28AD.
!       (See later comments.) It is not required to be set on entry
!       and, on exit, it contains information about the decomposition.
!       It should be preserved between this call and subsequent calls
!       to MA28BD or MA30CD.
!       IKEEP(I,1),I = 1,N holds the total length of the part of row
!       I in the diagonal block.
!       Row IKEEP(I,2),I = 1,N of the input matrix is the Ith row in
!       pivot order.
!       Column IKEEP(I,3),I = 1,N of the input matrix is the Ith
!       Column in pivot order.
!       IKEEP(I,4),I = 1,N holds the length of the part of row I in
!       the L part of the LU decomposition.
!       IKEEP(I,5),I = 1,N holds the length of the part of row I in
!       the off-diagonal blocks. If there is only one diagonal block,
!       IKEEP(1,5) will be set to -1.
! IW    Array of length 8*N. If the option NSRCH <= N is used, then
!       the length of array IW can be reduced to 7*N.
! W     Array length N. Used by MC24AD both as workspace and to return
!       growth estimate in W(1). The use of this array by MA28AD is
!       thus optional depending on logical variable GROW.
! IFLAG Variable used as error flag by subroutine. A positive or
!       zero value on exit indicates success. Possible negative
!       values are -1 through -14.

! Private Variable Information.
! LP, MP Default value 6 (line printer). Unit number for error
!     messages and duplicate element warning, respectively.
! NLP, MLP INTEGER. Unit number for messages from MA30AD and
!     MC23AD. Set by MA28AD to the value of LP.
! LBLOCK Logical variable with default value .TRUE. If .TRUE.,
!     MC23AD is used to first permute the matrix to block lower
!     triangular form.
! GROW Logical variable with default value .TRUE. If .TRUE., then
!     an estimate of the increase in size of matrix elements during
!     LU decomposition is given by MC24AD.
! EPS, RMIN, RESID. Variables not referenced by MA28AD.
! IRNCP, ICNCP INTEGER. Set to number of compresses on arrays IRN
!     and ICN/A, respectively.
! MINIRN, MINICN INTEGER. Minimum length of arrays IRN and ICN/A,
!     respectively, for success on future runs.
! IRANK INTEGER. Estimated rank of matrix.
! MIRNCP, MICNCP, MIRANK, MIRN, MICN INTEGER. Variables used to
!     communicate between MA30FD and MA28FD values of
!     abovenamed variables with somewhat similar names.
! ABORT1, ABORT2 LOGICAL. Variables with default value .TRUE.
!     If .FALSE., then decomposition will be performed even
!     if the matrix is structurally or numerically singular,
!     respectively.
! ABORTA, ABORTB LOGICAL. Variables used to communicate values
!     of ABORT1 and ABORT2 to MA30AD.
! ABORT Logical variable used to communicate value of ABORT1
!     to MC23AD.
! ABORT3 Logical variable. Not referenced by MA28AD.
! IDISP Array of length 2. Used to communicate information
!     on decomposition between this call to MA28AD and subsequent
!     calls to MA28BD and MA30CD. On exit, IDISP(1) and
!     IDISP(2) indicate position in arrays A and ICN of the
!     first and last elements in the LU decomposition of the
!     diagonal blocks, respectively.
! NUMNZ Structural rank of matrix.
! NUM   Number of diagonal blocks.
! LARGE Size of largest diagonal block.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP) :: U
        integer :: IFLAG, LICN, LIRN, N, NZ
! ..
! .. Array Arguments ..
        real (WP) :: A(LICN), W(N)
        integer :: ICN(LICN), IKEEP(N,5), IRN(LIRN), IW(N,8)
! ..
! .. Local Scalars ..
        real (WP) :: UPRIV
        integer :: I, I1, IEND, II, J, J1, J2, JAY, JJ, KNUM, LENGTH, MOVE, &
          NEWJ1, NEWPOS
        character (80) :: MSG
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, MAX
! ..
! .. FIRST EXECUTABLE STATEMENT MA28AD
! ..
!       Check that this call was made from DVODE_F90 and, if not, stop.
        if (.not.OK_TO_CALL_MA28) then
          MSG = 'This version of MA28 may be used only in conjunction with DVODE_F90.'
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = 'Please refer to the following HSL copyright notice for MA28.'
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '******************************************************************     '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '             *****MA28 COPYRIGHT NOTICE*****                           '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = ' COPYRIGHT (C) 2001 COUNCIL FOR THE CENTRAL LABORATORY                 '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '               OF THE RESEARCH COUNCILS                                '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = ' ALL RIGHTS RESERVED.                                                  '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = ' NONE OF THE COMMENTS IN THIS COPYRIGHT NOTICE BETWEEN THE LINES       '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = ' OF ASTERISKS SHALL BE REMOVED OR ALTERED IN ANY WAY.                  '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = ' THIS PACKAGE IS INTENDED FOR COMPILATION WITHOUT MODIFICATION,        '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = ' SO MOST OF THE EMBEDDED COMMENTS HAVE BEEN REMOVED.                   '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = ' ALL USE IS SUBJECT TO LICENCE. IF YOU NEED FURTHER CLARIFICATION,     '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = ' PLEASE SEE HTTP://WWW.HSL-LIBRARY.COM OR CONTACT HSL@HYPROTECH.COM    '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = ' PLEASE NOTE THAT:                                                     '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = ' 1. THE PACKAGES MAY ONLY BE USED FOR THE PURPOSES SPECIFIED IN THE    '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    LICENCE AGREEMENT AND MUST NOT BE COPIED BY THE LICENSEE FOR       '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    USE BY ANY OTHER PERSONS. USE OF THE PACKAGES IN ANY COMMERCIAL    '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    APPLICATION SHALL BE SUBJECT TO PRIOR WRITTEN AGREEMENT BETWEEN    '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    HYPROTECH UK LIMITED AND THE LICENSEE ON SUITABLE TERMS AND        '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    CONDITIONS, WHICH WILL INCLUDE FINANCIAL CONDITIONS.               '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = ' 2. ALL INFORMATION ON THE PACKAGE IS PROVIDED TO THE LICENSEE ON      '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    THE UNDERSTANDING THAT THE DETAILS THEREOF ARE CONFIDENTIAL.       '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = ' 3. ALL PUBLICATIONS ISSUED BY THE LICENSEE THAT INCLUDE RESULTS       '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    OBTAINED WITH THE HELP OF ONE OR MORE OF THE PACKAGES SHALL        '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    ACKNOWLEDGE THE USE OF THE PACKAGES. THE LICENSEE WILL NOTIFY      '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    HSL@HYPROTECH.COM OR HYPROTECH UK LIMITED OF ANY SUCH PUBLICATION. '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = ' 4. THE PACKAGES MAY BE MODIFIED BY OR ON BEHALF OF THE LICENSEE       '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    FOR SUCH USE IN RESEARCH APPLICATIONS BUT AT NO TIME SHALL SUCH    '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    PACKAGES OR MODIFICATIONS THEREOF BECOME THE PROPERTY OF THE       '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    LICENSEE. THE LICENSEE SHALL MAKE AVAILABLE FREE OF CHARGE TO THE  '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    COPYRIGHT HOLDER FOR ANY PURPOSE ALL INFORMATION RELATING TO       '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    ANY MODIFICATION.                                                  '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = ' 5. NEITHER COUNCIL FOR THE CENTRAL LABORATORY OF THE RESEARCH         '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    COUNCILS NOR HYPROTECH UK LIMITED SHALL BE LIABLE FOR ANY          '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    DIRECT OR CONSEQUENTIAL LOSS OR DAMAGE WHATSOEVER ARISING OUT OF   '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '    THE USE OF PACKAGES BY THE LICENSEE.                               '
          call XERRDV(MSG,1770,1,0,0,0,0,ZERO,ZERO)
          MSG = '******************************************************************     '
          call XERRDV(MSG,1770,2,0,0,0,0,ZERO,ZERO)
        end if

!       Some initialization and transfer of information between
!       common blocks (see earlier comments).
        IFLAG = 0
        ABORTA = ABORT1
        ABORTB = ABORT2
        ABORT = ABORT1
        MLP = LP
        NLP = LP
        TOL1 = TOL
        LBIG1 = LBIG
        NSRCH1 = NSRCH
!       UPRIV private copy of U is used in case it is outside
!       range zero to one and is thus altered by MA30AD.
        UPRIV = U
!       Simple data check on input variables and array dimensions.
        if (N>0) goto 10
        IFLAG = -8
        if (LP/=0) write (LP,90000) N
        goto 170
10      if (NZ>0) goto 20
        IFLAG = -9
        if (LP/=0) write (LP,90001) NZ
        goto 170
20      if (LICN>=NZ) goto 30
        IFLAG = -10
        if (LP/=0) write (LP,90002) LICN
        goto 170
30      if (LIRN>=NZ) goto 40
        IFLAG = -11
        if (LP/=0) write (LP,90003) LIRN
        goto 170

!       Data check to see if all indices lie between 1 and N.
40      do 50 I = 1, NZ
          if (IRN(I)>0 .and. IRN(I)<=N .and. ICN(I)>0 .and. ICN(I)<=N) &
            goto 50
          if (IFLAG==0 .and. LP/=0) write (LP,90004)
          IFLAG = -12
          if (LP/=0) write (LP,90005) I, A(I), IRN(I), ICN(I)
50      end do
        if (IFLAG<0) goto 180

!       Sort matrix into row order.
        call MC20AD(N,NZ,A,ICN,IW(1,1),IRN,0)
!       Part of IKEEP is used here as a work-array. IKEEP(I,2) is
!       the last row to have a non-zero in column I. IKEEP(I,3)
!       is the off-set of column I from the start of the row.
        IKEEP(1:N,2) = 0
        IKEEP(1:N,1) = 0

!       Check for duplicate elements summing any such entries
!       and printing a warning message on unit MP.
!       MOVE is equal to the number of duplicate elements found.
        MOVE = 0
!       The loop also calculates the largest element in the matrix,
!       THEMAX.
        THEMAX = ZERO
!       J1 is position in arrays of first non-zero in row.
        J1 = IW(1,1)
        do 90 I = 1, N
          IEND = NZ + 1
          if (I/=N) IEND = IW(I+1,1)
          LENGTH = IEND - J1
          if (LENGTH==0) goto 90
          J2 = IEND - 1
          NEWJ1 = J1 - MOVE
          do 80 JJ = J1, J2
            J = ICN(JJ)
            THEMAX = max(THEMAX,abs(A(JJ)))
            if (IKEEP(J,2)==I) goto 70
!           First time column has ocurred in current row.
            IKEEP(J,2) = I
            IKEEP(J,3) = JJ - MOVE - NEWJ1
            if (MOVE==0) goto 80
!           Shift necessary because of previous duplicate element.
            NEWPOS = JJ - MOVE
            A(NEWPOS) = A(JJ)
            ICN(NEWPOS) = ICN(JJ)
            goto 80
!           Duplicate element.
70          MOVE = MOVE + 1
            LENGTH = LENGTH - 1
            JAY = IKEEP(J,3) + NEWJ1
            if (MP/=0) write (MP,90006) I, J, A(JJ)
            A(JAY) = A(JAY) + A(JJ)
            THEMAX = max(THEMAX,abs(A(JAY)))
80        end do
          IKEEP(I,1) = LENGTH
          J1 = IEND
90      end do

!       KNUM is actual number of non-zeros in matrix with any multiple
!       entries counted only once.
        KNUM = NZ - MOVE
        if (.not.LBLOCK) goto 100

!       Perform block triangularisation.
        call MC23AD(N,ICN,A,LICN,IKEEP(1,1),IDISP,IKEEP(1,2),IKEEP(1,3), &
          IKEEP(1,5),IW(1,3),IW)
        if (IDISP(1)>0) goto 130
        IFLAG = -7
        if (IDISP(1)==-1) IFLAG = -1
        if (LP/=0) write (LP,90007)
        goto 170

!       Block triangularization not requested.
!       Move structure to end of data arrays in preparation for MA30AD.
!       Also set LENOFF(1) to -1 and set permutation arrays.
100     do I = 1, KNUM
          II = KNUM - I + 1
          NEWPOS = LICN - I + 1
          ICN(NEWPOS) = ICN(II)
          A(NEWPOS) = A(II)
        end do
        IDISP(1) = 1
        IDISP(2) = LICN - KNUM + 1
        do I = 1, N
          IKEEP(I,2) = I
          IKEEP(I,3) = I
        end do
        IKEEP(1,5) = -1
130     if (LBIG) BIG1 = THEMAX
        if (NSRCH<=N) goto 140

!       Perform LU decomposition on diagonal blocks.
        call MA30AD(N,ICN,A,LICN,IKEEP(1,1),IKEEP(1,4),IDISP,IKEEP(1,2),       &
          IKEEP(1,3),IRN,LIRN,IW(1,2),IW(1,3),IW(1,4),IW(1,5),IW(1,6),IW(1,7), &
          IW(1,8),IW,UPRIV,IFLAG)
        goto 150
!       This call if used if NSRCH has been set less than or equal to N.
!       In this case, two integer work arrays of length can be saved.
140     call MA30AD(N,ICN,A,LICN,IKEEP(1,1),IKEEP(1,4),IDISP,IKEEP(1,2),       &
          IKEEP(1,3),IRN,LIRN,IW(1,2),IW(1,3),IW(1,4),IW(1,5),IW,IW,IW(1,6),IW,&
          UPRIV,IFLAG)

!       Transfer private variable information.
!150    MINIRN = MAX(MIRN,NZ)
!       MINICN = MAX(MICN,NZ)
150     MINIRN = max(MINIRN,NZ)
        MINICN = max(MINICN,NZ)
!       IRNCP = MIRNCP
!       ICNCP = MICNCP
!       IRANK = MIRANK
!       NDROP = NDROP1
        if (LBIG) BIG = BIG1
        if (IFLAG>=0) goto 160
        if (LP/=0) write (LP,90008)
        goto 170

!       Reorder off-diagonal blocks according to pivot permutation.
160     I1 = IDISP(1) - 1
        if (I1/=0) call MC22AD(N,ICN,A,I1,IKEEP(1,5),IKEEP(1,2), &
          IKEEP(1,3),IW,IRN)
        I1 = IDISP(1)
        IEND = LICN - I1 + 1

!       Optionally calculate element growth estimate.
        if (GROW) call MC24AD(N,ICN,A(I1),IEND,IKEEP,IKEEP(1,4),W)
!       Increment growth estimate by original maximum element.
        if (GROW) W(1) = W(1) + THEMAX
        if (GROW .and. N>1) W(2) = THEMAX
!       Set flag if the only error is due to duplicate elements.
        if (IFLAG>=0 .and. MOVE/=0) IFLAG = -14
        goto 180
170     if (LP/=0) write (LP,90009)
180     return
90000   format (' N is out of range = ',I10)
90001   format (' NZ is non positive = ',I10)
90002   format (' LICN is too small = ',I10)
90003   format (' LIRN is too small = ',I10)
90004   format (' Error return from MA28AD because indices found out', &
          ' of range')
90005   format (1X,I6,'The element with value ',1P,D22.14, &
          ' is out of range with indices ',I8,',',I8)
90006   format (' Duplicate element in position ',I8,',',I8,' with value ',1P, &
          D22.14)
90007   format (' Error return from MC23AD')
90008   format (' Error return from MA30AD')
90009   format (' Error return from MA28AD')

      end subroutine MA28AD
!_______________________________________________________________________

      subroutine MA28BD(N,NZ,A,LICN,IVECT,JVECT,ICN,IKEEP,IW,W,IFLAG)
! ..
! This subroutine factorizes a matrix of a similar sparsity
! pattern to that previously factorized by MA28AD.
! ..
! The parameters are as follows:
! N            Order of matrix. Not altered by subroutine.
! NZ           Number of non-zeros in input matrix. Not
!              altered by subroutine.
! A            Array of length LICN. Holds non-zeros of
!              matrix on entry and non-zeros of factors on exit.
!              Reordered by MA28DD and altered by MA30BD.
! LICN         Length of arrays A and ICN. Not altered by
!              subroutine.
! IVECT, JVECT Arrays of length NZ. Hold row and column
!              indices of non-zeros, respectively. Not altered
!              by subroutine.
! ICN          Array of length LICN. Same array as output
!              from MA28AD. Unchanged by MA28BD.
! IKEEP        Array of length 5*N. Same array as output
!              from MA28AD. Unchanged by MA28BD.
! IW           Array length 5*N. Used as workspace by
!              MA28DD and MA30BD.
! W            Array of length N. Used as workspace by
!              MA28DD, MA30BD and (optionally) MC24AD.
! IFLAG        Integer used as error flag, with positive
!              or zero value indicating success.

! Private Variable Information.
! Unless otherwise stated private variables are as in MA28AD.
! Those variables referenced by MA28BD are mentioned below.
! LP, MP     Integers used as in MA28AD as unit number for error
!            and warning messages, respectively.
! NLP        Integer variable used to give value of LP to MA30ED.
! EPS        MA30BD will output a positive value
!            for IFLAG if any modulus of the ratio of pivot element
!            to the largest element in its row (U part only) is less
!            than EPS (unless EPS is greater than 1.0 when no action
!            takes place).
! RMIN       Variable equal to the value of this minimum ratio in
!            cases where EPS is less than or equal to 1.0.
! MEPS,MRMIN Variables used by the subroutine to communicate between
!            MA28FD and MA30GD.
! IDISP      Integer array of length 2. The same as that used by
!            MA28AD. Unchanged by MA28BD.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer :: IFLAG, LICN, N, NZ
! ..
! .. Array Arguments ..
        real (WP) :: A(LICN), W(N)
        integer :: ICN(LICN), IKEEP(N,5), IVECT(NZ), IW(N,5), JVECT(NZ)
! ..
! .. Local Scalars ..
        integer :: I1, IDUP, IEND
! ..
! .. FIRST EXECUTABLE STATEMENT MA28BD
! ..
!       Check to see if elements were dropped in previous MA28AD call.
        if (NDROP==0) goto 10
        IFLAG = -15
        if (LP/=0) write (LP,90000) IFLAG, NDROP
        goto 70
10      IFLAG = 0
        MEPS = EPS
        NLP = LP
!       Simple data check on variables.
        if (N>0) goto 20
        IFLAG = -11
        if (LP/=0) write (LP,90001) N
        goto 60
20      if (NZ>0) goto 30
        IFLAG = -10
        if (LP/=0) write (LP,90002) NZ
        goto 60
30      if (LICN>=NZ) goto 40
        IFLAG = -9
        if (LP/=0) write (LP,90003) LICN
        goto 60

40      call MA28DD(N,A,LICN,IVECT,JVECT,NZ,ICN,IKEEP(1,1),IKEEP(1,4), &
          IKEEP(1,5),IKEEP(1,2),IKEEP(1,3),IW(1,3),IW,W(1),IFLAG)
!       THEMAX is largest element in matrix.
        THEMAX = W(1)
        if (LBIG) BIG1 = THEMAX
!       IDUP equals one if there were duplicate elements, zero otherwise.
        IDUP = 0
        if (IFLAG==(N+1)) IDUP = 1
        if (IFLAG<0) goto 60

!       Perform row Gauss elimination on the structure received from MA28DD.
        call MA30BD(N,ICN,A,LICN,IKEEP(1,1),IKEEP(1,4),IDISP,IKEEP(1,2), &
          IKEEP(1,3),W,IW,IFLAG)

!       Transfer private variable information.
        if (LBIG) BIG1 = BIG
!       RMIN = MRMIN
        if (IFLAG>=0) goto 50
        IFLAG = -2
        if (LP/=0) write (LP,90004)
        goto 60

!       Optionally calculate the growth parameter.
50      I1 = IDISP(1)
        IEND = LICN - I1 + 1
        if (GROW) call MC24AD(N,ICN,A(I1),IEND,IKEEP,IKEEP(1,4),W)
!       Increment estimate by largest element in input matrix.
        if (GROW) W(1) = W(1) + THEMAX
        if (GROW .and. N>1) W(2) = THEMAX
!       Set flag if the only error is due to duplicate elements.
        if (IDUP==1 .and. IFLAG>=0) IFLAG = -14
        goto 70
60      if (LP/=0) write (LP,90005)
70      return
90000   format (' Error return from MA28BD with IFLAG = ',I4/I7, &
          ' Entries dropped from structure by MA28AD')
90001   format (' N is out of range = ',I10)
90002   format (' NZ is non positive = ',I10)
90003   format (' LICN is too small = ',I10)
90004   format (' Error return from MA30BD')
90005   format (' + Error return from MA28BD')

      end subroutine MA28BD
!_______________________________________________________________________

      subroutine MA28DD(N,A,LICN,IVECT,JVECT,NZ,ICN,LENR,LENRL,LENOFF, &
        IP,IQ,IW1,IW,W1,IFLAG)
! ..
! This subroutine need never be called by the user directly. It sorts
! the user's matrix into the structure of the decomposed form and checks
! for the presence of duplicate entries or non-zeros lying outside the
! sparsity pattern of the decomposition. It also calculates the largest
! element in the input matrix.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP) :: W1
        integer :: IFLAG, LICN, N, NZ
! ..
! .. Array Arguments ..
        real (WP) :: A(LICN)
        integer :: ICN(LICN), IP(N), IQ(N), IVECT(NZ), IW(N,2), IW1(N,3), &
          JVECT(NZ), LENOFF(N), LENR(N), LENRL(N)
! ..
! .. Local Scalars ..
        real (WP) :: AA
        integer :: I, IBLOCK, IDISP2, IDUMMY, II, INEW, IOLD, J1, J2, JCOMP, &
          JDUMMY, JJ, JNEW, JOLD, MIDPT
        logical :: BLOCKL
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, IABS, MAX
! ..
! .. FIRST EXECUTABLE STATEMENT MA28DD
! ..
        BLOCKL = LENOFF(1) >= 0
!       IW1(I,3) is set to the block in which row I lies and the
!       inverse permutations to IP and IQ are set in IW1(:,1) and
!       IW1(:,2), respectively.
!       Pointers to beginning of the part of row I in diagonal and
!       off-diagonal blocks are set in IW(I,2) and IW(I,1),
!       respectively.
        IBLOCK = 1
        IW(1,1) = 1
        IW(1,2) = IDISP(1)
        do 10 I = 1, N
          IW1(I,3) = IBLOCK
          if (IP(I)<0) IBLOCK = IBLOCK + 1
          II = IABS(IP(I)+0)
!         II = IABS(IP(I))
          IW1(II,1) = I
          JJ = IQ(I)
          JJ = IABS(JJ)
          IW1(JJ,2) = I
          if (I==1) goto 10
          if (BLOCKL) IW(I,1) = IW(I-1,1) + LENOFF(I-1)
          IW(I,2) = IW(I-1,2) + LENR(I-1)
10      end do
!       Place each non-zero in turn into its correct location in
!       the A/ICN array.
        IDISP2 = IDISP(2)
        do 170 I = 1, NZ
!         Necessary to avoid reference to unassigned element of ICN.
          if (I>IDISP2) goto 20
          if (ICN(I)<0) goto 170
20        IOLD = IVECT(I)
          JOLD = JVECT(I)
          AA = A(I)
!         This is a dummy loop for following a chain of interchanges.
!         It will be executed NZ TIMES in total.
          do IDUMMY = 1, NZ
!           Perform some validity checks on IOLD and JOLD.
            if (IOLD<=N .and. IOLD>0 .and. JOLD<=N .and. JOLD>0) goto 30
            if (LP/=0) write (LP,90000) I, A(I), IOLD, JOLD
            IFLAG = -12
            goto 180
30          INEW = IW1(IOLD,1)
            JNEW = IW1(JOLD,2)
!           Are we in a valid block and is it diagonal or off-diagonal?
!           IF (IW1(INEW,3)-IW1(JNEW,3)) 40, 60, 50
!40         IFLAG = -13
            if (IW1(INEW,3)-IW1(JNEW,3) == 0) goto 60
            if (IW1(INEW,3)-IW1(JNEW,3) > 0) goto 50
            IFLAG = -13
            if (LP/=0) write (LP,90001) IOLD, JOLD
            goto 180
50          J1 = IW(INEW,1)
            J2 = J1 + LENOFF(INEW) - 1
            goto 110
!           Element is in diagonal block.
60          J1 = IW(INEW,2)
            if (INEW>JNEW) goto 70
            J2 = J1 + LENR(INEW) - 1
            J1 = J1 + LENRL(INEW)
            goto 110
70          J2 = J1 + LENRL(INEW)
!           Binary search of ordered list. Element in L part of row.
            do 100 JDUMMY = 1, N
              MIDPT = (J1+J2)/2
              JCOMP = IABS(ICN(MIDPT)+0)
!             JCOMP = IABS(ICN(MIDPT))
!             IF (JNEW-JCOMP) 80, 130, 90
!80           J2 = MIDPT
              if (JNEW-JCOMP == 0) goto 130
              if (JNEW-JCOMP > 0) goto 90
              J2 = MIDPT
              goto 100
90            J1 = MIDPT
100         end do
            IFLAG = -13
            if (LP/=0) write (LP,90002) IOLD, JOLD
            goto 180
!           Linear search. Element in L part of row or off-diagonal blocks.
110         do MIDPT = J1, J2
              if (IABS(ICN(MIDPT)+0)==JNEW) goto 130
!             IF (IABS(ICN(MIDPT))==JNEW) GOTO 130
            end do
            IFLAG = -13
            if (LP/=0) write (LP,90002) IOLD, JOLD
            goto 180
!           Equivalent element of ICN is in position MIDPT.
130         if (ICN(MIDPT)<0) goto 160
            if (MIDPT>NZ .or. MIDPT<=I) goto 150
            W1 = A(MIDPT)
            A(MIDPT) = AA
            AA = W1
            IOLD = IVECT(MIDPT)
            JOLD = JVECT(MIDPT)
            ICN(MIDPT) = -ICN(MIDPT)
          end do
150       A(MIDPT) = AA
          ICN(MIDPT) = -ICN(MIDPT)
          goto 170
160       A(MIDPT) = A(MIDPT) + AA
!        Set flag for duplicate elements.
          iflag = n + 1
170     end do
!       Reset ICN array and zero elements in LU but not in A.
!       Also calculate the maximum element of A.
180     W1 = ZERO
        do 200 I = 1, IDISP2
          if (ICN(I)<0) goto 190
          A(I) = ZERO
          goto 200
190       ICN(I) = -ICN(I)
          W1 = max(W1,abs(A(I)))
200     end do
        return
90000   format (' Element ',I6,' with value ',D22.14,' has indices ',I8, &
          ','/I8,' indices out of range')
90001   format (' Non-zero ',I7,',',I6,' in zero off-diagonal',' block')
90002   format (' Element ',I6,',',I6,' was not in LU pattern')

      end subroutine MA28DD
!_______________________________________________________________________

      subroutine MA28CD(N,A,LICN,ICN,IKEEP,RHS,W,MTYPE)
! ..
! This subroutine uses the factors from MA28AD or MA28BD to solve a
! system of equations without iterative refinement.
! ..
! The parameters are:
! N     Order of matrix. Not altered by subroutine.
! A     array of length LICN. the same array as
!       was used in the most recent call to MA28AD or MA28BD.
! LICN  Length of arrays A and ICN. Not altered by subroutine.
! ICN   Array of length LICN. Same array as output from
!       MA28AD. Unchanged by MA30CD.
! IKEEP Array of length 5*N. Same array as output from
!       MA28AD. Unchanged by MA30CD.
! RHS   Array of length N. On entry, it holds the
!       right hand side, on exit, the solution vector.
! W     Array of length N. Used as workspace by MA30CD.
! MTYPE Integer used to tell MA30CD to solve the direct
!       equation (MTYPE = 1) or its transpose (MTYPE /= 1).

! Private Variable Information.
! Unless otherwise stated private variables are as in MA28AD.
! Those variables referenced by MA30CD are mentioned below.
! RESID   Variable which returns maximum residual of
!         equations where pivot was zero.
! MRESID  Variable used by MA30CD to communicate between
!         MA28FD and MA30HD.
! IDISP   Integer array of length 2. The same as that used
!         by MA28AD. Unchanged by MA28BD.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer :: LICN, MTYPE, N
! ..
! .. Array Arguments ..
        real (WP) :: A(LICN), RHS(N), W(N)
        integer :: ICN(LICN), IKEEP(N,5)
! ..
! .. FIRST EXECUTABLE STATEMENT MA28CD
! ..
!       This call performs the solution of the set of equations.
        call MA30CD(N,ICN,A,LICN,IKEEP(1,1),IKEEP(1,4),IKEEP(1,5),IDISP, &
          IKEEP(1,2),IKEEP(1,3),RHS,W,MTYPE)
!       Transfer private variable information.
        RESID = MRESID
        return

      end subroutine MA28CD
!_______________________________________________________________________

      subroutine MA28ID(N,NZ,AORG,IRNORG,ICNORG,LICN,A,ICN,IKEEP,RHS, &
        X,R,W,MTYPE,PREC,IFLAG)
! ..
! This subroutine uses the factors from an earlier call to MA28AD or
! MA28BD to solve the system of equations with iterative refinement.
! ..
! The parameters are:
! N      Order of the matrix. Not altered by the subroutine.
! NZ     Number of entries in the original matrix. Not altered by
!        the subroutine.
!        For this entry the original matrix must have been saved in
!        AORG, IRNORG, ICNORG where entry AORG(K) is in row IRNORG(K)
!        and column ICNORG(K),K = 1,...,NZ. Information about the
!        factors of A is communicated to this subroutine via the
!        parameters LICN, A, ICN and IKEEP where:
! AORG   Array of length NZ. Not altered by MA28ID.
! IRNORG Array of length NZ. Not altered by MA28ID.
! ICNORG Array of length NZ. Not altered by MA28ID.
! LICN   Length of arrays A and ICN. Not altered by the subroutine.
! A      Array of length LICN. It must be unchanged since the
!        last call to MA28AD or MA28BD. Not altered by the
!        subroutine.
! ICN, IKEEP are the arrays (of lengths LICN and 5*N, respectively)
!        of the same names as in the previous all to MA28AD. They
!        should be unchanged since this earlier call and they are
!        not altered by MA28ID.
! The other parameters are as follows:
! RHS    Array of length N. The user must set RHS(I) to contain
!        the value of the Ith component of the right hand side.
!        Not altered by MA28ID.
! X      Array of length N. If an initial guess of the solution
!        is given (ISTART equal to 1), then the user must set X(I)
!        to contain the value of the Ith component of the estimated
!        solution. On exit, X(I) contains the Ith component of the
!        solution vector.
! R      Array of length N. It need not be set on entry. On exit,
!        R(I) contains the Ith component of an estimate of the error
!        if MAXIT is greater than 0.
! W      Array of length N. Used as workspace by MA28ID.
! MTYPE  Must be set to determine whether MA28ID will solve A*X = RHS
!        (MTYPE = 1) or AT*X = RHS (MTYPE /= 1). Not altered by MA28ID.
! PREC   Should be set by the user to the relative accuracy required.
!        The iterative refinement will terminate if the magnitude of
!        the largest component of the estimated error relative to the
!        largest component in the solution is less than PREC.
!        Not altered by MA28ID.
! IFLAG  Diagnostic flag which will be set to zero on successful
!        exit from MA28ID. Otherwise it will have a non-zero value.
!        The non-zero value IFLAG can have on exit from MA28ID are
!        -16 indicating that more than MAXIT iteartions are required.
!        -17 indicating that more convergence was too slow.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP) :: PREC
        integer :: IFLAG, LICN, MTYPE, N, NZ
! ..
! .. Array Arguments ..
        real (WP) :: A(LICN), AORG(NZ), R(N), RHS(N), W(N), X(N)
        integer :: ICN(LICN), ICNORG(NZ), IKEEP(N,5), IRNORG(NZ)
! ..
! .. Local Scalars ..
        real (WP) :: CONVER, D, DD
        integer :: I, ITERAT, NCOL, NROW
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, MAX
! ..
! .. FIRST EXECUTABLE STATEMENT MA28ID
! ..
!       Initialization of NOITER, ERRMAX, and IFLAG.
        NOITER = 0
        ERRMAX = ZERO
        IFLAG = 0

!       Jump if a starting vector has been supplied by the user.

        if (ISTART==1) goto 20

!       Make a copy of the right-hand side vector.

        X(1:N) = RHS(1:N)

!       Find the first solution.

        call MA28CD(N,A,LICN,ICN,IKEEP,X,W,MTYPE)

!       Stop the computations if MAXIT = 0.

20      if (MAXIT==0) goto 160

!       Calculate the max-norm of the first solution.

        DD = 0.0
        do I = 1, N
          DD = max(DD,abs(X(I)))
        end do
        DXMAX = DD

!       Begin the iterative process.

        do 120 ITERAT = 1, MAXIT
          D = DD

!         Calculate the residual vector.

          R(1:N) = RHS(1:N)
          if (MTYPE==1) goto 60
          do I = 1, NZ
            NROW = IRNORG(I)
            NCOL = ICNORG(I)
            R(NCOL) = R(NCOL) - AORG(I)*X(NROW)
          end do
          goto 80
!         MTYPE = 1.
60        do I = 1, NZ
            NROW = IRNORG(I)
            NCOL = ICNORG(I)
            R(NROW) = R(NROW) - AORG(I)*X(NCOL)
          end do
80        DRES = 0.0

!         Find the max-norm of the residual vector.

          do I = 1, N
            DRES = max(DRES,abs(R(I)))
          end do

!         Stop the calculations if the max-norm of
!         the residual vector is zero.

!         IF (DRES == 0.0) GOTO 150
          if (abs(DRES)<=ZERO) goto 150

!         Calculate the correction vector.

          NOITER = NOITER + 1
          call MA28CD(N,A,LICN,ICN,IKEEP,R,W,MTYPE)

!         Find the max-norm of the correction vector.

          DD = 0.0
          do I = 1, N
            DD = max(DD,abs(R(I)))
          end do

!         Check the convergence.

          if (DD>D*CGCE .and. ITERAT>=2) goto 130
          if (abs((DXMAX*TEN+DD)-(DXMAX*TEN))<=ZERO) goto 140

!         Attempt to improve the solution.

          DXMAX = 0.0
          do I = 1, N
            X(I) = X(I) + R(I)
            DXMAX = max(DXMAX,abs(X(I)))
          end do

!         Check the stopping criterion.
 
          if (DD<PREC*DXMAX) goto 140
120     end do
!       More than MAXIT iterations required.
        IFLAG = -16
        write (LP,90000) IFLAG, MAXIT
        goto 140
!       Convergence rate unacceptably slow.
130     IFLAG = -17
        CONVER = DD/D
        write (LP,90001) IFLAG, CONVER, CGCE

!       The iterative process is terminated.

140     ERRMAX = DD
150     continue
160     return
90000   format (' Error return from MA28ID with IFLAG = ',I3/' More than ', &
          I5,' iterations required')
90001   format (' Error return from MA28I with IFLAG = ', &
          I3/' Convergence rate of ',1P,E9.2,' too slow'/ &
          ' Maximum acceptable rate set to ',1P,E9.2)

      end subroutine MA28ID
!_______________________________________________________________________

! Private Variable Information.
! LP, MP are used by the subroutine as the unit numbers for its warning
!     and diagnostic messages. Default value for both is 6 (for line
!     printer output). The user can either reset them to a different
!     stream number or suppress the output by setting them to zero.
!     While LP directs the output of error diagnostics from the
!     principal subroutines and internally called subroutines, MP
!     controls only the output of a message which warns the user that
!     he has input two or more non-zeros A(I),. .,A(K) with the same
!     row and column indices. The action taken in this case is to
!     proceed using a numerical value of A(I)+...+A(K). In the absence
!     of other errors, IFLAG will equal -14 on exit.
! LBLOCK Is a logical variable which controls an option of first
!     preordering the matrix to block lower triangular form (using
!     Harwell subroutine MC23A). The preordering is performed if LBLOCK
!     is equal to its default value of .TRUE. If LBLOCK is set to
!     .FALSE., the option is not invoked and the space allocated to
!     IKEEP can be reduced to 4*N+1.
! GROW is a logical variable. If it is left at its default value of
!     .TRUE., then on return from MA28AD or MA28BD, W(1) will give
!     an estimate (an upper bound) of the increase in size of elements
!     encountered during the decomposition. If the matrix is well
!     scaled, then a high value for W(1), relative to the largest entry
!     in the input matrix, indicates that the LU decomposition may be
!     inaccurate and the user should be wary of his results and perhaps
!     increase U for subsequent runs. We would like to emphasise that
!     this value only relates to the accuracy of our LU decomposition
!     and gives no indication as to the singularity of the matrix or the
!     accuracy of the solution. This upper bound can be a significant
!     overestimate particularly if the matrix is badly scaled. If an
!     accurate value for the growth is required, LBIG(Q.V.) should be
!     set to .TRUE.
! EPS, RMIN If, on entry to MA28BD, EPS is less than one, then RMIN
!     will give the smallest ratio of the pivot to the largest element
!     in the corresponding row of the upper triangular factor thus
!     monitoring the stability of successive factorizations. If RMIN
!     becomes very large and W(1) from MA28BD is also very large, it
!     may be advisable to perform a new decomposition using MA28AD.
! RESID Variable which on exit from MA30CD gives the value of the
!     maximum residual over all the equations unsatisfied because
!     of dependency (zero pivots).
! IRNCP, ICNCP Variables which monitor the adequacy of "elbow room"
!     in IRN and A/ICN, respectively. If either is quite large (say
!     greater than N/10), it will probably pay to increase the size of
!     the corresponding array for subsequent runs. If either is very low
!     or zero, then one can perhaps save storage by reducing the size of
!     the corresponding array.
! MINIRN, MINICN Integer variables which, in the event of a successful
!     return (IFLAG >= 0 or IFLAG = -14) give the minimum size of IRN
!     and A/ICN, respectively which would enable a successful run on
!     an identical matrix. On an exit with IFLAG equal to -5, MINICN
!     gives the minimum value of ICN for success on subsequent runs on
!     an identical matrix. In the event of failure with IFLAG = -6,-4,
!     -3,-2, OR -1, then MINICN and MINIRN give the minimum value of
!     LICN and LIRN, respectively which would be required for a
!     successful decomposition up to the point at which the failure
!     occurred.
! IRANK Integer variable which gives an upper bound on the rank of
!     the matrix.
! ABORT1 is a logical variable with default value .TRUE. If ABORT1 is
!     set to .FALSE., then MA28AD will decompose structurally singular
!     matrices (including rectangular ones).
! ABORT2 is a logical variable with default value .TRUE. If ABORT2 is
!     set to .FALSE., then MA28AD will decompose numerically singular
!     matrices.
! IDISP is an integer array of length 2. On output from MA28AD, The
!     indices of the diagonal blocks of the factors lie in positions
!     IDISP(1) to IDISP(2) of A/ICN. This array must be preserved
!     between a call to MA28AD and subsequent calls to MA28BD,
!     MA30CD or MA28ID.
! TOL If set to a positive value, then any non-zero whose modulus is
!     less than TOL will be dropped from the factorization. The
!     factorization will then require less storage but will be
!     inaccurate. After a run of MA28AD with TOL positive it is not
!     possible to use MA28BD and the user is recommended to use
!     MA28ID to obtain the solution. The default value for TOL is 0.0.
! THEMAX On exit from MA28AD, THEMAX will hold the largest entry of
!     the original matrix.
! BIG If LBIG has been set to .TRUE., BIG will hold the largest entry
!     encountered during the factorization by MA28AD or MA28BD.
! DXMAX On exit from MA28ID, DXMAX will be set to the largest component
!     of the solution.
! ERRMAX If MAXIT is positive, ERRMAX will be set to the largest
!     component in the estimate of the error.
! DRES On exit from MA28ID, if MAXIT is positive, DRES will be set to
!     the largest component of the residual.
! CGCE Used by MA28ID to check the convergence rate. If the ratio of
!     successive corrections is not less than CGCE, then we terminate
!     since the convergence rate is adjudged too slow.
! NDROP If TOL has been set positive on exit from MA28AD, NDROP will
!     hold the number of entries dropped from the data structure.
! MAXIT Maximum number of iterations performed by MA28ID. Default = 16.
! NOITER Set by MA28ID to the number of iterative refinement iterations
!     actually used.
! NSRCH If NSRCH is set to a value less than N, then a different pivot
!     option will be employed by MA28AD. This may result in different
!     fill-in and execution time for MA28AD. If NSRCH is less than or
!     equal to N, the workspace array IW can be reduced in length. The
!     default value for NSRCH is 32768.
! ISTART If ISTART is set to a value other than zero, then the user
!     must supply an estimate of the solution to MA28ID. The default
!     value for istart is zero.
! LBIG If LBIG is set to .TRUE., the value of the largest element
!     encountered in the factorization by MA28AD or MA28BD is returned
!     in BIG. Setting LBIG to .TRUE. will increase the time for MA28AD
!     marginally and that for MA28BD by about 20%. The default value
!     for LBIG is .FALSE.
!_______________________________________________________________________

      subroutine MA30AD(NN,ICN,A,LICN,LENR,LENRL,IDISP,IP,IQ,IRN,LIRN, &
        LENC,IFIRST,LASTR,NEXTR,LASTC,NEXTC,IPTR,IPC,U,IFLAG)
! ..
!     If the user requires a more convenient data interface then the
!     MA28 package should be used. The MA28 subroutines call the MA30
!     subroutines after checking the user's input data and optionally
!     using MC23AD to permute the matrix to block triangular form.
!     This package of subroutines (MA30AD, MA30BD, MA30CD, and MA30DD)
!     performs operations pertinent to the solution of a general
!     sparse N by N system of linear equations (i.e., solve
!     AX = B). Structually singular matrices are permitted including
!     those with row or columns consisting entirely of zeros (i.e.,
!     including rectangular matrices). It is assumed that the
!     non-zeros of the matrix A do not differ widely in size. If
!     necessary a prior call of the scaling subroutine MA19AD may be
!     made.
!     A discussion of the design of these subroutines is given by Duff
!     and Reid (ACM TRANS MATH SOFTWARE 5 PP 18-35, 1979(CSS 48)) while
!     fuller details of the implementation are given in duff (HARWELL
!     REPORT AERE-R 8730, 1977). The additional pivoting option in
!     ma30ad and the use of drop tolerances (see private variables for
!     MA30ID) were added to the package after joint work with Reid,
!     Schaumburg, Wasniewski and Zlatev (Duff, Reid, Schaumburg,
!     Wasniewski and Zlatev, HARWELL REPORT CSS 135, 1983).

!     MA30AD performs the LU decomposition of the diagonal blocks of
!     the permutation PAQ of a sparse matrix A, where input permutations
!     P1 and Q1 are used to define the diagonal blocks. There may be
!     non-zeros in the off-diagonal blocks but they are unaffected by
!     MA30AD. P and P1 differ only within blocks as do Q and Q1. The
!     permutations P1 and Q1 may be found by calling MC23AD or the
!     matrix may be treated as a single block by using P1 = Q1 = I. The
!     matrix non-zeros should be held compactly by rows, although it
!     should be noted that the user can supply the matrix by columns
!     to get the LU decomposition of A transpose.

!     This description should also be consulted for further information
!     on most of the parameters of MA30BD and MA30CD.
!     The parameters are:
! N   is an integer variable which must be set by the user to the order
!     of the matrix. It is not altered by MA30AD.
! ICN is an integer array of length LICN. Positions IDISP(2) to
!     LICN must be set by the user to contain the column indices of
!     the non-zeros in the diagonal blocks of P1*A*Q1. Those belonging
!     to a single row must be contiguous but the ordering of column
!     indices with each row is unimportant. The non-zeros of row I
!     precede those of row I+1,I = 1,...,N-1 and no wasted space is
!     allowed between the rows. On output the column indices of the
!     LU decomposition of PAQ are held in positions IDISP(1) to
!     IDISP(2), the rows are in pivotal order, and the column indices
!     of the L Part of each row are in pivotal order and precede those
!     of U. Again there is no wasted space either within a row or
!     between the rows. ICN(1) to ICN(IDISP(1)-1), are neither
!     required nor altered. If MC23AD has been called, these will hold
!     information about the off-diagonal blocks.
! A   Array of length LICN whose entries IDISP(2) to LICN must be set
!     by the user to the values of the non-zero entries of the matrix
!     in the order indicated by ICN. On output A will hold the LU
!     factors of the matrix where again the position in the matrix
!     is determined by the corresponding values in ICN.
!     A(1) to A(IDISP(1)-1) are neither required nor altered.
! LICN is an integer variable which must be set by the user to the
!     length of arrays ICN and A. It must be big enough for A and ICN
!     to hold all the non-zeros of L and U and leave some "elbow
!     room". It is possible to calculate a minimum value for LICN by
!     a preliminary run of MA30AD. The adequacy of the elbow room
!     can be judged by the size of the private variable ICNCP.
!     It is not altered by MA30AD.
! LENR is an integer array of length N. On input, LENR(I) should
!     equal the number of non-zeros in row I,I = 1,...,N of the
!     diagonal blocks of P1*A*Q1. On output, LENR(I) will equal the
!     total number of non-zeros in row I of L and row I of U.
! LENRL is an integer array of length N. On output from MA30AD,
!     LENRL(I) will hold the number of non-zeros in row I of L.
! IDISP is an integer array of length 2. The user should set IDISP(1)
!     to be the first available position in A/ICN for the LU
!     decomposition while IDISP(2) is set to the position in A/ICN of
!     the first non-zero in the diagonal blocks of P1*A*Q1. On output,
!     IDISP(1) will be unaltered while IDISP(2) will be set to the
!     position in A/ICN of the last non-zero of the LU decomposition.
! IP  is an integer array of length N which holds a permutation of
!     the integers 1 to N. On input to MA30AD, the absolute value of
!     IP(I) must be set to the row of A which is row I of P1*A*Q1. A
!     negative value for IP(I) indicates that row I is at the end of a
!     diagonal block. On output from MA30AD, IP(I) indicates the row
!     of A which is the Ith row in PAQ. IP(I) will still be negative
!     for the last row of each block (except the last).
! IQ  is an integer array of length N which again holds a
!     permutation of the integers 1 to N. On input to MA30AD, IQ(J)
!     must be set to the column of A which is column J of P1*A*Q1. On
!     output from MA30AD, the absolute value of IQ(J) indicates the
!     column of A which is the Jth in PAQ. For rows, I say, in which
!     structural or numerical singularity is detected IQ(I) is
!     negated.
! IRN is an integer array of length LIRN used as workspace by
!     MA30AD.
! LIRN is an integer variable. It should be greater than the
!     largest number of non-zeros in a diagonal block of P1*A*Q1 but
!     need not be as large as LICN. It is the length of array IRN and
!     should be large enough to hold the active part of any block,
!     plus some "elbow room", the a posteriori adequacy of which can
!     be estimated by examining the size of private variable
!     IRNCP.
! LENC, FIRST, LASTR, NEXTR, LASTC, NEXTC are all integer arrays of
!     length N which are used as workspace by MA30AD. If NSRCH is
!     set to a value less than or equal to N, then arrays LASTC and
!     NEXTC are not referenced by MA30AD and so can be dummied in
!     the call to MA30AD.
! IPTR, IPC are integer arrays of length N which are used as
!     workspace by MA30AD.
! U   is a real(kind=wp) variable which should be set by the user
!     to a value between 0. and 1.0. If less than zero it is reset
!     to zero and if its value is 1.0 or greater it is reset to
!     0.9999 (0.999999999 in D version). It determines the balance
!     between pivoting for sparsity and for stability, values near
!     zero emphasizing sparsity and values near one emphasizing
!     stability. We recommend U = 0.1 as a posible first trial value.
!     The stability can be judged by a later call to MC24AD or by
!     setting LBIG to .TRUE.
! IFLAG is an integer variable. It will have a non-negative value
!     if MA30AD is successful. Negative values indicate error
!     conditions while positive values indicate that the matrix has
!     been successfully decomposed but is singular. For each non-zero
!     value, an appropriate message is output on unit LP. Possible
!     non-zero values for IFLAG are
! -1  THe matrix is structurally singular with rank given by IRANK
!     in private variables for MA30FD.
! +1  If, however, the user wants the LU decomposition of a
!     structurally singular matrix and sets the private variable
!     ABORT1 to .FALSE., then, in the event of singularity and a
!     successful decomposition, Iflag is returned with the value +1
!     and no message is output.
! -2  The matrix is numerically singular (it may also be structurally
!     singular) with estimated rank given by IRANK in private variables
!     for MA30FD.
! +2  THE user can choose to continue the decomposition even when a
!     zero pivot is encountered by setting private variable
!     ABORT2 TO .FALSE. If a singularity is encountered, IFLAG will
!     then return with a value of +2, and no message is output if
!     the decomposition has been completed successfully.
! -3  LIRN has not been large enough to continue with the
!     decomposition. If the stage was zero then private variable
!     MINIRN gives the length sufficient to start the decomposition on
!     this block. FOr a successful decomposition on this block the
!     user should make LIRN slightly (say about N/2) greater than this
!     value.
! -4  LICN is not large enough to continue with the decomposition.
! -5  The decomposition has been completed but some of the LU factors
!     have been discarded to create enough room in A/ICN to continue
!     the decomposition. The variable MINICN in private variables for
!     MA30FD. Then gives the size that LICN should be to enable the
!     factorization to be successful. If the user sets private
!     variable ABORT3 to .TRUE., then the subroutine will exit
!     immediately instead of destroying any factors and continuing.
! -6  Both LICN and LIRN are too small. Termination has been caused
!     by lack of space in IRN (see error IFLAG = -3), but already
!     some of the LU factors in A/ICN have been lost (see error
!     IFLAG = -5). MINICN gives the minimum amount of space
!     required in A/ICN for decomposition up to this point.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        real (WP) :: U
        integer :: IFLAG, LICN, LIRN, NN
! ..
! .. Array Arguments ..
        real (WP) :: A(LICN)
        integer :: ICN(LICN), IDISP(2), IFIRST(NN), IP(NN), IPC(NN), &
          IPTR(NN), IQ(NN), IRN(LIRN), LASTC(NN), LASTR(NN),         &
          LENC(NN), LENR(NN), LENRL(NN), NEXTC(NN), NEXTR(NN)
! ..
! .. Local Scalars ..
        real (WP) :: AANEW, AMAX, ANEW, AU, PIVR, PIVRAT, SCALE
        integer :: COLUPD, DISPC, I, I1, I2, IACTIV, IBEG, IDISPC,   &
          IDROP, IDUMMY, IEND, IFILL, IFIR, II, III, IJFIR, IJP1,    &
          IJPOS, ILAST, INDROW, IOP, IPIV, IPOS, IROWS, ISING, ISRCH,&
          ISTART, ISW, ISW1, ITOP, J, J1, J2, JBEG, JCOST, JCOUNT,   &
          JDIFF, JDUMMY, JEND, JJ, JMORE, JNEW, JNPOS, JOLD, JPIV,   &
          JPOS, JROOM, JVAL, JZER, JZERO, K, KCOST, KDROP, L, LC,    &
          LENPIV, LENPP, LL, LR, MOREI, MSRCH, N, NBLOCK, NC, NNM1,  &
          NR, NUM, NZ, NZ2, NZCOL, NZMIN, NZPC, NZROW, OLDEND,       &
          OLDPIV, PIVEND, PIVOT, PIVROW, ROWI
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, IABS, MAX, MIN
! ..
! .. FIRST EXECUTABLE STATEMENT MA30AD
! ..
        MSRCH = NSRCH
        NDROP = 0
        LNPIV(1:10) = 0
        LPIV(1:10) = 0
        MAPIV = 0
        MANPIV = 0
        IAVPIV = 0
        IANPIV = 0
        KOUNTL = 0
        MINIRN = 0
        MINICN = IDISP(1) - 1
        MOREI = 0
        IRANK = NN
        IRNCP = 0
        ICNCP = 0
        IFLAG = 0
!       Reset U if necessary.
        U = min(U,UMAX)
!       IBEG is the position of the next pivot row after elimination
!       step using it.
        U = max(U,ZERO)
        IBEG = IDISP(1)
!       IACTIV is the position of the first entry in the active part
!       of A/ICN.
        IACTIV = IDISP(2)
!       NZROW is current number of non-zeros in active and unprocessed
!       part of row file ICN.
        NZROW = LICN - IACTIV + 1
        MINICN = NZROW + MINICN

!       Count the number of diagonal blocks and set up pointers to the
!       beginnings of the rows. NUM is the number of diagonal blocks.
        NUM = 1
        IPTR(1) = IACTIV
        if (NN==1) goto 30
        NNM1 = NN - 1
        do I = 1, NNM1
          if (IP(I)<0) NUM = NUM + 1
          IPTR(I+1) = IPTR(I) + LENR(I)
        end do
!       ILAST is the last row in the previous block.
30      ILAST = 0

! ***********************************************
! ****    LU decomposition of block NBLOCK   ****
! ***********************************************

!       Each pass through this loop performs LU decomposition on one
!       of the diagonal blocks.
        do 1070 NBLOCK = 1, NUM
          ISTART = ILAST + 1
          do IROWS = ISTART, NN
            if (IP(IROWS)<0) goto 50
          end do
          IROWS = NN
50        ILAST = IROWS
!         N is the number of rows in the current block.
!         ISTART is the index of the first row in the current block.
!         ILAST is the index of the last row in the current block.
!         IACTIV is the position of the first entry in the block.
!         ITOP is the position of the last entry in the block.
          N = ILAST - ISTART + 1
          if (N/=1) goto 100

!         Code for dealing WITH 1x1 block.
          LENRL(ILAST) = 0
          ISING = ISTART
          if (LENR(ILAST)/=0) goto 60
!         Block is structurally singular.
          IRANK = IRANK - 1
          ISING = -ISING
          if (IFLAG/=2 .and. IFLAG/=-5) IFLAG = 1
          if (.not.ABORT1) goto 90
          IDISP(2) = IACTIV
          IFLAG = -1
          if (LP/=0) write (LP,90000)
!         RETURN
          goto 1190
60        SCALE = abs(A(IACTIV))
          if (abs(SCALE)<=ZERO) goto 70
          if (LBIG) BIG = max(BIG,SCALE)
          goto 80
70        ISING = -ISING
          IRANK = IRANK - 1
          IPTR(ILAST) = 0
          if (IFLAG/=-5) IFLAG = 2
          if (.not.ABORT2) goto 80
          IDISP(2) = IACTIV
          IFLAG = -2
          if (LP/=0) write (LP,90001)
          goto 1190
80        A(IBEG) = A(IACTIV)
          ICN(IBEG) = ICN(IACTIV)
          IACTIV = IACTIV + 1
          IPTR(ISTART) = 0
          IBEG = IBEG + 1
          NZROW = NZROW - 1
90        LASTR(ISTART) = ISTART
          IPC(ISTART) = -ISING
          goto 1070

!         Non-trivial block.
100       ITOP = LICN
          if (ILAST/=NN) ITOP = IPTR(ILAST+1) - 1

!         Set up column oriented storage.
          LENRL(ISTART:ILAST) = 0
          LENC(ISTART:ILAST) = 0
          if (ITOP-IACTIV<LIRN) goto 120
          MINIRN = ITOP - IACTIV + 1
          PIVOT = ISTART - 1
          goto 1170

!         Calculate column counts.
120       do II = IACTIV, ITOP
            I = ICN(II)
            LENC(I) = LENC(I) + 1
          end do
!         Set up column pointers so that IPC(J) points to position
!         after end of column J in column file.
          IPC(ILAST) = LIRN + 1
          J1 = ISTART + 1
          do JJ = J1, ILAST
            J = ILAST - JJ + J1 - 1
            IPC(J) = IPC(J+1) - LENC(J+1)
          end do
          do 160 INDROW = ISTART, ILAST
            J1 = IPTR(INDROW)
            J2 = J1 + LENR(INDROW) - 1
            if (J1>J2) goto 160
            do JJ = J1, J2
              J = ICN(JJ)
              IPOS = IPC(J) - 1
              IRN(IPOS) = INDROW
              IPC(J) = IPOS
            end do
160       end do
!         DISPC is the lowest indexed active location in the column file.
          DISPC = IPC(ISTART)
          NZCOL = LIRN - DISPC + 1
          MINIRN = max(NZCOL,MINIRN)
          NZMIN = 1

!         Initialize array IFIRST. IFIRST(I) = +/- K indicates that
!         row/col K has I non-zeros. If IFIRST(I) = 0, there is no
!         row or column with I non-zeros.
          IFIRST(1:N) = 0

!         Compute ordering of row and column counts.
!         First run through columns (from column N to column 1).
          do 190 JJ = ISTART, ILAST
            J = ILAST - JJ + ISTART
            NZ = LENC(J)
            if (NZ/=0) goto 180
            IPC(J) = 0
            goto 190
180         if (NSRCH<=NN) goto 190
            ISW = IFIRST(NZ)
            IFIRST(NZ) = -J
            LASTC(J) = 0
            NEXTC(J) = -ISW
            ISW1 = IABS(ISW)
            if (ISW/=0) LASTC(ISW1) = J
190       end do
!        Now run through rows (again from N to 1).
          do 220 II = ISTART, ILAST
               I = ILAST - II + ISTART
            NZ = LENR(I)
            if (NZ/=0) goto 200
            IPTR(I) = 0
            LASTR(I) = 0
            goto 220
200         ISW = IFIRST(NZ)
            IFIRST(NZ) = I
            if (ISW>0) goto 210
            NEXTR(I) = 0
            LASTR(I) = ISW
            goto 220
210         NEXTR(I) = ISW
            LASTR(I) = LASTR(ISW)
            LASTR(ISW) = I
220       end do

! **********************************************
! ****    Start of main elimination loop    ****
! **********************************************
          do 1050 PIVOT = ISTART, ILAST

!           First find the pivot using MARKOWITZ criterion with
!           stability control.
!           JCOST is the Markowitz cost of the best pivot so far.
!           This pivot is in row IPIV and column JPIV.
            NZ2 = NZMIN
            JCOST = N*N

!           Examine rows/columns in order of ascending count.
            do L = 1, 2
              PIVRAT = ZERO
              ISRCH = 1
              LL = L
!             A pass with L equal to 2 is only performed in the
!             case of singularity.
              do 340 NZ = NZ2, N
                if (JCOST<=(NZ-1)**2) goto 430
                IJFIR = IFIRST(NZ)
!               IF (IJFIR) 240, 230, 250
!230            IF (LL==1) NZMIN = NZ + 1
                if (IJFIR < 0) goto 240
                if (IJFIR > 0) goto 250
                if (LL==1) NZMIN = NZ + 1
                goto 340
240             LL = 2
                IJFIR = -IJFIR
                goto 300
250             LL = 2
!               Scan rows with NZ non-zeros.
                do IDUMMY = 1, N
                  if (JCOST<=(NZ-1)**2) goto 430
                  if (ISRCH>MSRCH) goto 430
                  if (IJFIR==0) goto 290
!                 Row IJFIR is now examined.
                  I = IJFIR
                  IJFIR = NEXTR(I)
!                 First calculate multiplier threshold level.
                  AMAX = ZERO
                  J1 = IPTR(I) + LENRL(I)
                  J2 = IPTR(I) + LENR(I) - 1
                  do JJ = J1, J2
                    AMAX = max(AMAX,abs(A(JJ)))
                  end do
                  AU = AMAX*U
                  ISRCH = ISRCH + 1
!                 Scan row for possible pivots.
                  do 270 JJ = J1, J2
                    if (abs(A(JJ))<=AU .and. L==1) goto 270
                    J = ICN(JJ)
                    KCOST = (NZ-1)*(LENC(J)-1)
                    if (KCOST>JCOST) goto 270
                    PIVR = ZERO
                    if (abs(AMAX)>ZERO) PIVR = abs(A(JJ))/AMAX
                    if (KCOST==JCOST .and. (PIVR<=PIVRAT .or. NSRCH>NN+1)) &
                      goto 270
!                   Best pivot so far is found.
                    JCOST = KCOST
                    IJPOS = JJ
                    IPIV = I
                    JPIV = J
                    if (MSRCH>NN+1 .and. JCOST<=(NZ-1)**2) goto 430
                    PIVRAT = PIVR
270               end do
                end do

!               Columns with NZ non-zeros now examined.
290             IJFIR = IFIRST(NZ)
                IJFIR = -LASTR(IJFIR)
300             if (JCOST<=NZ*(NZ-1)) goto 430
                if (MSRCH<=NN) goto 340
                do 330 IDUMMY = 1, N
                  if (IJFIR==0) goto 340
                  J = IJFIR
                  IJFIR = NEXTC(IJFIR)
                  I1 = IPC(J)
                  I2 = I1 + NZ - 1
!                 Scan column J.
                  do 320 II = I1, I2
                    I = IRN(II)
                    KCOST = (NZ-1)*(LENR(I)-LENRL(I)-1)
                    if (KCOST>=JCOST) goto 320
!                   Pivot has best Markowitz count so far. Now
!                   check its suitability on numeric grounds by
!                   examining the other non-zeros in its row.
                    J1 = IPTR(I) + LENRL(I)
                    J2 = IPTR(I) + LENR(I) - 1
!                   We need a stability check on singleton columns
!                   because of possible problems with
!                   underdetermined systems.
                    AMAX = ZERO
                    do JJ = J1, J2
                      AMAX = max(AMAX,abs(A(JJ)))
                      if (ICN(JJ)==J) JPOS = JJ
                    end do
                    if (abs(A(JPOS))<=AMAX*U .and. L==1) goto 320
                    JCOST = KCOST
                    IPIV = I
                    JPIV = J
                    IJPOS = JPOS
                    if (abs(AMAX)>ZERO) PIVRAT = abs(A(JPOS))/AMAX
                    if (JCOST<=NZ*(NZ-1)) goto 430
320               end do
330             end do
340           end do
!             In the event of singularity, we must make sure all
!             rows and columns are tested.
              MSRCH = N

!             Matrix is numerically or structurally singular. Which
!             it is will be diagnosed later.
              IRANK = IRANK - 1
            end do
!           Assign rest of rows and columns to ordering array.
!           Matrix is structurally singular.
            if (IFLAG/=2 .and. IFLAG/=-5) IFLAG = 1
            IRANK = IRANK - ILAST + PIVOT + 1
            if (.not.ABORT1) goto 360
            IDISP(2) = IACTIV
            IFLAG = -1
            if (LP/=0) write (LP,90000)
            goto 1190
360         K = PIVOT - 1
            do 400 I = ISTART, ILAST
              if (LASTR(I)/=0) goto 400
              K = K + 1
              LASTR(I) = K
              if (LENRL(I)==0) goto 390
              MINICN = max(MINICN,NZROW+IBEG-1+MOREI+LENRL(I))
              if (IACTIV-IBEG>=LENRL(I)) goto 370
              call MA30DD(A,ICN,IPTR(ISTART),N,IACTIV,ITOP,.true.)
!             Check now to see if MA30DD has created enough
!             available space.
              if (IACTIV-IBEG>=LENRL(I)) goto 370
!             Create more space by destroying previously created
!             LU factors.
              MOREI = MOREI + IBEG - IDISP(1)
              IBEG = IDISP(1)
              if (LP/=0) write (LP,90002)
              IFLAG = -5
              if (ABORT3) goto 1160
370           J1 = IPTR(I)
              J2 = J1 + LENRL(I) - 1
              IPTR(I) = 0
              do JJ = J1, J2
                A(IBEG) = A(JJ)
                ICN(IBEG) = ICN(JJ)
                ICN(JJ) = 0
                IBEG = IBEG + 1
              end do
              NZROW = NZROW - LENRL(I)
390           if (K==ILAST) goto 410
400         end do
410         K = PIVOT - 1
            do 420 I = ISTART, ILAST
              if (IPC(I)/=0) goto 420
              K = K + 1
              IPC(I) = K
              if (K==ILAST) goto 1060
420         end do

!           The pivot has now been found in position (IPIV,JPIV)
!           in location IJPOS in row file.
!           Update column and row ordering arrays to correspond
!           with removal of the active part of the matrix.
430         ISING = PIVOT
            if (abs(A(IJPOS))>ZERO) goto 440
!           Numerical singularity is recorded here.
            ISING = -ISING
            if (IFLAG/=-5) IFLAG = 2
            if (.not.ABORT2) goto 440
            IDISP(2) = IACTIV
            IFLAG = -2
            if (LP/=0) write (LP,90001)
            goto 1190
440         OLDPIV = IPTR(IPIV) + LENRL(IPIV)
            OLDEND = IPTR(IPIV) + LENR(IPIV) - 1
!           Changes to column ordering.
            if (NSRCH<=NN) goto 470
            COLUPD = NN + 1
            LENPP = OLDEND - OLDPIV + 1
            if (LENPP<4) LPIV(1) = LPIV(1) + 1
            if (LENPP>=4 .and. LENPP<=6) LPIV(2) = LPIV(2) + 1
            if (LENPP>=7 .and. LENPP<=10) LPIV(3) = LPIV(3) + 1
            if (LENPP>=11 .and. LENPP<=15) LPIV(4) = LPIV(4) + 1
            if (LENPP>=16 .and. LENPP<=20) LPIV(5) = LPIV(5) + 1
            if (LENPP>=21 .and. LENPP<=30) LPIV(6) = LPIV(6) + 1
            if (LENPP>=31 .and. LENPP<=50) LPIV(7) = LPIV(7) + 1
            if (LENPP>=51 .and. LENPP<=70) LPIV(8) = LPIV(8) + 1
            if (LENPP>=71 .and. LENPP<=100) LPIV(9) = LPIV(9) + 1
            if (LENPP>=101) LPIV(10) = LPIV(10) + 1
            MAPIV = max(MAPIV,LENPP)
            IAVPIV = IAVPIV + LENPP
            do 460 JJ = OLDPIV, OLDEND
              J = ICN(JJ)
              LC = LASTC(J)
              NC = NEXTC(J)
              NEXTC(J) = -COLUPD
              if (JJ/=IJPOS) COLUPD = J
              if (NC/=0) LASTC(NC) = LC
              if (LC==0) goto 450
              NEXTC(LC) = NC
              goto 460
450           NZ = LENC(J)
              ISW = IFIRST(NZ)
              if (ISW>0) LASTR(ISW) = -NC
              if (ISW<0) IFIRST(NZ) = -NC
460         end do
!           Changes to row ordering.
470         I1 = IPC(JPIV)
            I2 = I1 + LENC(JPIV) - 1
            do 490 II = I1, I2
              I = IRN(II)
              LR = LASTR(I)
              NR = NEXTR(I)
              if (NR/=0) LASTR(NR) = LR
              if (LR<=0) goto 480
              NEXTR(LR) = NR
              goto 490
480           NZ = LENR(I) - LENRL(I)
              if (NR/=0) IFIRST(NZ) = NR
              if (NR==0) IFIRST(NZ) = LR
490         end do

!           Move pivot to position LENRL+1 in pivot row and move pivot
!           row to the beginning of the available storage. The L part
!           and the pivot in the old copy of the pivot row is nullified
!           while, in the strictly upper triangular part, the column
!           indices, J say, are overwritten by the corresponding entry
!           of IQ (IQ(J)) and IQ(J) is set to the negative of the
!           displacement of the column index from the pivot entry.
            if (OLDPIV==IJPOS) goto 500
            AU = A(OLDPIV)
            A(OLDPIV) = A(IJPOS)
            A(IJPOS) = AU
            ICN(IJPOS) = ICN(OLDPIV)
            ICN(OLDPIV) = JPIV
!           Check to see if there is space immediately available in
!           A/ICN to hold new copy of pivot row.
500         MINICN = max(MINICN,NZROW+IBEG-1+MOREI+LENR(IPIV))
            if (IACTIV-IBEG>=LENR(IPIV)) goto 510
            call MA30DD(A,ICN,IPTR(ISTART),N,IACTIV,ITOP,.true.)
            OLDPIV = IPTR(IPIV) + LENRL(IPIV)
            OLDEND = IPTR(IPIV) + LENR(IPIV) - 1
!           Check now to see if MA30DD has created enough
!           available space.
            if (IACTIV-IBEG>=LENR(IPIV)) goto 510
!           Create more space by destroying previously created
!           LU factors.
            MOREI = MOREI + IBEG - IDISP(1)
            IBEG = IDISP(1)
            if (LP/=0) write (LP,90002)
            IFLAG = -5
            if (ABORT3) goto 1160
            if (IACTIV-IBEG>=LENR(IPIV)) goto 510
!           There is still not enough room in A/ICN.
            IFLAG = -4
            goto 1160
!           Copy pivot row and set up IQ array.
510         IJPOS = 0
            J1 = IPTR(IPIV)

            do JJ = J1, OLDEND
              A(IBEG) = A(JJ)
              ICN(IBEG) = ICN(JJ)
              if (IJPOS/=0) goto 520
              if (ICN(JJ)==JPIV) IJPOS = IBEG
              ICN(JJ) = 0
              goto 530
520           K = IBEG - IJPOS
              J = ICN(JJ)
              ICN(JJ) = IQ(J)
              IQ(J) = -K
530           IBEG = IBEG + 1
            end do

            IJP1 = IJPOS + 1
            PIVEND = IBEG - 1
            LENPIV = PIVEND - IJPOS
            NZROW = NZROW - LENRL(IPIV) - 1
            IPTR(IPIV) = OLDPIV + 1
            if (LENPIV==0) IPTR(IPIV) = 0

!           Remove pivot row (including pivot) from column
!           oriented file.
            do JJ = IJPOS, PIVEND
              J = ICN(JJ)
              I1 = IPC(J)
              LENC(J) = LENC(J) - 1
!             I2 is last position in new column.
              I2 = IPC(J) + LENC(J) - 1
              if (I2<I1) goto 560
              do 550 II = I1, I2
                if (IRN(II)/=IPIV) goto 550
                IRN(II) = IRN(I2+1)
                goto 560
550           end do
560           IRN(I2+1) = 0
            end do
            NZCOL = NZCOL - LENPIV - 1

!           Go down the pivot column and for each row with a
!           non-zero add the appropriate multiple of the pivot
!           row to it. We loop on the number of non-zeros in
!           the pivot column since MA30DD may change its
!           actual position.

            NZPC = LENC(JPIV)
            if (NZPC==0) goto 940
            do 880 III = 1, NZPC
              II = IPC(JPIV) + III - 1
              I = IRN(II)
!             Search row I For non-zero to be eliminated,
!             calculate multiplier, and place it in position
!             LENRL+1 in its row. IDROP is the number of
!             non-zero entries dropped from row I because
!             these fall beneath tolerance level.
              IDROP = 0
              J1 = IPTR(I) + LENRL(I)
              IEND = IPTR(I) + LENR(I) - 1
              do 580 JJ = J1, IEND
                if (ICN(JJ)/=JPIV) goto 580
!               IF pivot is zero, rest of column is and so
!               multiplier is zero.
                AU = ZERO
                if (abs(A(IJPOS))>ZERO) AU = -A(JJ)/A(IJPOS)
                if (LBIG) BIG = max(BIG,abs(AU))
                A(JJ) = A(J1)
                A(J1) = AU
                ICN(JJ) = ICN(J1)
                ICN(J1) = JPIV
                LENRL(I) = LENRL(I) + 1
                goto 590
580           end do
!             JUMP if pivot row is a singleton.
590           if (LENPIV==0) goto 880
!             Now perform necessary operations on rest of non-pivot
!             row I.
              ROWI = J1 + 1
              IOP = 0
!             Jump if all the pivot row causes fill-in.
              if (ROWI>IEND) goto 670
!             Perform operations on current non-zeros in row I.
!             Innermost loop.
              LENPP = IEND - ROWI + 1
              if (LENPP<4) LNPIV(1) = LNPIV(1) + 1
              if (LENPP>=4 .and. LENPP<=6) LNPIV(2) = LNPIV(2) + 1
              if (LENPP>=7 .and. LENPP<=10) LNPIV(3) = LNPIV(3) + 1
              if (LENPP>=11 .and. LENPP<=15) LNPIV(4) = LNPIV(4) + 1
              if (LENPP>=16 .and. LENPP<=20) LNPIV(5) = LNPIV(5) + 1
              if (LENPP>=21 .and. LENPP<=30) LNPIV(6) = LNPIV(6) + 1
              if (LENPP>=31 .and. LENPP<=50) LNPIV(7) = LNPIV(7) + 1
              if (LENPP>=51 .and. LENPP<=70) LNPIV(8) = LNPIV(8) + 1
              if (LENPP>=71 .and. LENPP<=100) LNPIV(9) = LNPIV(9) + 1
              if (LENPP>=101) LNPIV(10) = LNPIV(10) + 1
              MANPIV = max(MANPIV,LENPP)
              IANPIV = IANPIV + LENPP
              KOUNTL = KOUNTL + 1
              do 600 JJ = ROWI, IEND
                J = ICN(JJ)
                if (IQ(J)>0) goto 600
                IOP = IOP + 1
                PIVROW = IJPOS - IQ(J)
                A(JJ) = A(JJ) + AU*A(PIVROW)
                if (LBIG) BIG = max(abs(A(JJ)),BIG)
                ICN(PIVROW) = -ICN(PIVROW)
                if (abs(A(JJ))<TOL) IDROP = IDROP + 1
600           end do

!             JUMP if no non-zeros in non-pivot row have been removed
!             because these are beneath the drop-tolerance TOL.

              if (IDROP==0) goto 670

!             Run through non-pivot row compressing row so that only
!             non-zeros greater than TOL are stored. All non-zeros
!             less than TOL are also removed from the column structure.

              JNEW = ROWI
              do 650 JJ = ROWI, IEND
                if (abs(A(JJ))<TOL) goto 610
                A(JNEW) = A(JJ)
                ICN(JNEW) = ICN(JJ)
                JNEW = JNEW + 1
                goto 650

!               Remove non-zero entry from column structure.

610             J = ICN(JJ)
                I1 = IPC(J)
                I2 = I1 + LENC(J) - 1
                do II = I1, I2
                  if (IRN(II)==I) goto 630
                end do
630             IRN(II) = IRN(I2)
                IRN(I2) = 0
                LENC(J) = LENC(J) - 1
                if (NSRCH<=NN) goto 650
!               Remove column from column chain and place in
!               update chain.
                if (NEXTC(J)<0) goto 650
!               JUMP if column already in update chain.
                LC = LASTC(J)
                NC = NEXTC(J)
                NEXTC(J) = -COLUPD
                COLUPD = J
                if (NC/=0) LASTC(NC) = LC
                if (LC==0) goto 640
                NEXTC(LC) = NC
                goto 650
640             NZ = LENC(J) + 1
                ISW = IFIRST(NZ)
                if (ISW>0) LASTR(ISW) = -NC
                if (ISW<0) IFIRST(NZ) = -NC
650           end do
              ICN(JNEW:IEND) = 0
!             The value of IDROP might be different from that
!             calculated earlier because we may now have dropped
!             some non-zeros which were not modified by the pivot
!             row.
              IDROP = IEND + 1 - JNEW
              IEND = JNEW - 1
              LENR(I) = LENR(I) - IDROP
              NZROW = NZROW - IDROP
              NZCOL = NZCOL - IDROP
              NDROP = NDROP + IDROP
670           IFILL = LENPIV - IOP
!             Jump is if there is no fill-in.
              if (IFILL==0) goto 770
!             Now for the fill-in.
              MINICN = max(MINICN,MOREI+IBEG-1+NZROW+IFILL+LENR(I))
!             See if there is room for fill-in.
!             Get maximum space for row I in situ.
              do JDIFF = 1, IFILL
                JNPOS = IEND + JDIFF
                if (JNPOS>LICN) goto 690
                if (ICN(JNPOS)/=0) goto 690
              end do
!             There is room for all the fill-in after the end of the
!             row so it can be left in situ.
!             Next available space for fill-in.
              IEND = IEND + 1
              goto 770
!             JMORE spaces for fill-in are required in front of row.
690           JMORE = IFILL - JDIFF + 1
              I1 = IPTR(I)
!             We now look in front of the row to see if there is space
!             for the rest of the fill-in.
              do JDIFF = 1, JMORE
                JNPOS = I1 - JDIFF
                if (JNPOS<IACTIV) goto 710
                if (ICN(JNPOS)/=0) goto 720
              end do
710           JNPOS = I1 - JMORE
              goto 730
!             Whole row must be moved to the beginning of available
!             storage.
720           JNPOS = IACTIV - LENR(I) - IFILL
!             Jump if there is space immediately available for the
!             shifted row.
730           if (JNPOS>=IBEG) goto 750
              call MA30DD(A,ICN,IPTR(ISTART),N,IACTIV,ITOP,.true.)
              I1 = IPTR(I)
              IEND = I1 + LENR(I) - 1
              JNPOS = IACTIV - LENR(I) - IFILL
              if (JNPOS>=IBEG) goto 750
!             No space available so try to create some by throwing
!             away previous LU decomposition.
              MOREI = MOREI + IBEG - IDISP(1) - LENPIV - 1
              if (LP/=0) write (LP,90002)
              IFLAG = -5
              if (ABORT3) goto 1160
!             Keep record of current pivot row.
              IBEG = IDISP(1)
              ICN(IBEG) = JPIV
              A(IBEG) = A(IJPOS)
              IJPOS = IBEG
              do JJ = IJP1, PIVEND
                IBEG = IBEG + 1
                A(IBEG) = A(JJ)
                ICN(IBEG) = ICN(JJ)
              end do
              IJP1 = IJPOS + 1
              PIVEND = IBEG
              IBEG = IBEG + 1
              if (JNPOS>=IBEG) goto 750
!             This still does not give enough room.
              IFLAG = -4
              goto 1160
750           IACTIV = min(IACTIV,JNPOS)
!             Move non-pivot row I.
              IPTR(I) = JNPOS
              do JJ = I1, IEND
                A(JNPOS) = A(JJ)
                ICN(JNPOS) = ICN(JJ)
                JNPOS = JNPOS + 1
                ICN(JJ) = 0
              end do
!             First new available space.
              IEND = JNPOS
770           NZROW = NZROW + IFILL
!             Innermost fill-in loop which also resets ICN.
              IDROP = 0
              do 850 JJ = IJP1, PIVEND
                J = ICN(JJ)
                if (J<0) goto 840
                ANEW = AU*A(JJ)
                AANEW = abs(ANEW)
                if (AANEW>=TOL) goto 780
                IDROP = IDROP + 1
                NDROP = NDROP + 1
                NZROW = NZROW - 1
                MINICN = MINICN - 1
                IFILL = IFILL - 1
                goto 850
780             if (LBIG) BIG = max(AANEW,BIG)
                A(IEND) = ANEW
                ICN(IEND) = J
                IEND = IEND + 1

!               Put new entry in column file.
                MINIRN = max(MINIRN,NZCOL+LENC(J)+1)
                JEND = IPC(J) + LENC(J)
                JROOM = NZPC - III + 1 + LENC(J)
                if (JEND>LIRN) goto 790
                if (IRN(JEND)==0) goto 830
790             if (JROOM<DISPC) goto 800
!               Compress column file to obtain space for new
!               copy of column.
                call MA30DD(A,IRN,IPC(ISTART),N,DISPC,LIRN,.false.)
                if (JROOM<DISPC) goto 800
                JROOM = DISPC - 1
                if (JROOM>=LENC(J)+1) goto 800
!               Column file is not large enough.
                goto 1170
!               Copy column to beginning of file.
800             JBEG = IPC(J)
                JEND = IPC(J) + LENC(J) - 1
                JZERO = DISPC - 1
                DISPC = DISPC - JROOM
                IDISPC = DISPC
                do II = JBEG, JEND
                  IRN(IDISPC) = IRN(II)
                  IRN(II) = 0
                  IDISPC = IDISPC + 1
                end do
                IPC(J) = DISPC
                JEND = IDISPC
                IRN(JEND:JZERO) = 0
830             IRN(JEND) = I
                NZCOL = NZCOL + 1
                LENC(J) = LENC(J) + 1
!               End of adjustment to column file.
                goto 850
840             ICN(JJ) = -J
850           end do
              if (IDROP==0) goto 870
              do KDROP = 1, IDROP
                ICN(IEND) = 0
                IEND = IEND + 1
              end do
870           LENR(I) = LENR(I) + IFILL
!             End of scan of pivot column.
880         end do

!           Remove pivot column from column oriented storage
!           and update row ordering arrays.
            I1 = IPC(JPIV)
            I2 = IPC(JPIV) + LENC(JPIV) - 1
            NZCOL = NZCOL - LENC(JPIV)
            do 930 II = I1, I2
              I = IRN(II)
              IRN(II) = 0
              NZ = LENR(I) - LENRL(I)
              if (NZ/=0) goto 890
              LASTR(I) = 0
              goto 930
890           IFIR = IFIRST(NZ)
              IFIRST(NZ) = I
!             IF (IFIR) 900, 920, 910
!900          LASTR(I) = IFIR
              if (IFIR == 0) goto 920
              if (IFIR > 0) goto 910
              LASTR(I) = IFIR
              NEXTR(I) = 0
              goto 930
910           LASTR(I) = LASTR(IFIR)
              NEXTR(I) = IFIR
              LASTR(IFIR) = I
              goto 930
920           LASTR(I) = 0
              NEXTR(I) = 0
              NZMIN = min(NZMIN,NZ)
930         end do
!           Restore IQ and nullify U part of old pivot row.
!           Record the column permutation in LASTC(JPIV) and
!           the row permutation in LASTR(IPIV).
940         IPC(JPIV) = -ISING
            LASTR(IPIV) = PIVOT
            if (LENPIV==0) goto 1050
            NZROW = NZROW - LENPIV
            JVAL = IJP1
            JZER = IPTR(IPIV)
            IPTR(IPIV) = 0
            do JCOUNT = 1, LENPIV
              J = ICN(JVAL)
              IQ(J) = ICN(JZER)
              ICN(JZER) = 0
              JVAL = JVAL + 1
              JZER = JZER + 1
            end do
!           Adjust column ordering arrays.
            if (NSRCH>NN) goto 980
            do 970 JJ = IJP1, PIVEND
              J = ICN(JJ)
              NZ = LENC(J)
              if (NZ/=0) goto 960
              IPC(J) = 0
              goto 970
960           NZMIN = min(NZMIN,NZ)
970         end do
            goto 1050
980         JJ = COLUPD
            do 1040 JDUMMY = 1, NN
              J = JJ
              if (J==NN+1) goto 1050
              JJ = -NEXTC(J)
              NZ = LENC(J)
              if (NZ/=0) goto 990
              IPC(J) = 0
              goto 1040
990           IFIR = IFIRST(NZ)
              LASTC(J) = 0
!             IF (IFIR) 1000, 1010, 1020
!1000         IFIRST(NZ) = -J
              if (IFIR == 0) goto 1010
              if (IFIR > 0) goto 1020
              IFIRST(NZ) = -J
              IFIR = -IFIR
              LASTC(IFIR) = J
              NEXTC(J) = IFIR
              goto 1040
1010          IFIRST(NZ) = -J
              NEXTC(J) = 0
              goto 1030
1020          LC = -LASTR(IFIR)
              LASTR(IFIR) = -J
              NEXTC(J) = LC
              if (LC/=0) LASTC(LC) = J
1030          NZMIN = min(NZMIN,NZ)
1040        end do
1050      end do
! ********************************************
! ****    End of main elimination loop    ****
! ********************************************

!         Reset IACTIV to point to the beginning of the
!         next block.
1060      if (ILAST/=NN) IACTIV = IPTR(ILAST+1)
1070    end do

! ********************************************
! ****    End of deomposition of block    ****
! ********************************************

!     Record singularity (if any) in IQ array.
        if (IRANK==NN) goto 1090
        do 1080 I = 1, NN
          if (IPC(I)<0) goto 1080
          ISING = IPC(I)
          IQ(ISING) = -IQ(ISING)
          IPC(I) = -ISING
1080    end do

!       Run through LU decomposition changing column indices
!       to that of new order and permuting LENR and LENRL
!       arrays according to pivot permutations.
1090    ISTART = IDISP(1)
        IEND = IBEG - 1
        if (IEND<ISTART) goto 1110
        do JJ = ISTART, IEND
          JOLD = ICN(JJ)
          ICN(JJ) = -IPC(JOLD)
        end do
1110    do II = 1, NN
          I = LASTR(II)
          NEXTR(I) = LENR(II)
          IPTR(I) = LENRL(II)
        end do
        LENRL(1:NN) = IPTR(1:NN)
        LENR(1:NN) = NEXTR(1:NN)

!       Update permutation arrays IP and IQ.
        do II = 1, NN
          I = LASTR(II)
          J = -IPC(II)
          NEXTR(I) = IABS(IP(II)+0)
          IPTR(J) = IABS(IQ(II)+0)
!         NEXTR(I) = IABS(IP(II))
!         IPTR(J) = IABS(IQ(II))
        end do
        do I = 1, NN
          if (IP(I)<0) NEXTR(I) = -NEXTR(I)
          IP(I) = NEXTR(I)
          if (IQ(I)<0) IPTR(I) = -IPTR(I)
          IQ(I) = IPTR(I)
        end do
        IP(NN) = IABS(IP(NN)+0)
!       IP(NN) = IABS(IP(NN))
        IDISP(2) = IEND
        goto 1190

!       Error returns
1160    IDISP(2) = IACTIV
        if (LP==0) goto 1190
        write (LP,90003)
        goto 1180
1170    if (IFLAG==-5) IFLAG = -6
        if (IFLAG/=-6) IFLAG = -3
        IDISP(2) = IACTIV
        if (LP==0) goto 1190
        if (IFLAG==-3) write (LP,90004)
        if (IFLAG==-6) write (LP,90005)
1180    PIVOT = PIVOT - ISTART + 1
        write (LP,90006) PIVOT, NBLOCK, ISTART, ILAST
        if (PIVOT==0) write (LP,90007) MINIRN

1190    return
90000   format (' Error return from MA30AD because matrix is', &
          ' structurally singular')
90001   format (' Error return from MA30AD because matrix is', &
          ' numerically singular')
90002   format (' LU decomposition destroyed to create more space')
90003   format (' Error return from MA30AD because LICN is not big enough')
90004   format (' Error return from MA30AD because LIRN is not big enough')
90005   format (' Error return from MA30AD LIRN and LICN are too small')
90006   format (' At stage ',I5,' in block ',I5,' with first row ',I5, &
          ' and last row ',I5)
90007   format (' To continue set LIRN to at least ',I8)

      end subroutine MA30AD
!_______________________________________________________________________

      subroutine MA30DD(A,ICN,IPTR,N,IACTIV,ITOP,REALS)
! ..
! This subroutine performs garbage collection operations on the arrays
! A, ICN, and IRN.
! ..
! IACTIV is the first position in arrays A/ICN from which the compress
! starts. On exit, IACTIV equals the position of the first entry.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer :: IACTIV, ITOP, N
        logical :: REALS
! ..
! .. Array Arguments ..
        real (WP) :: A(ITOP)
        integer :: ICN(ITOP), IPTR(N)
! ..
! .. Local Scalars ..
        integer :: J, JPOS, K, KL, KN
! ..
! .. FIRST EXECUTABLE STATEMENT MA30DD
! ..
        if (REALS) ICNCP = ICNCP + 1
        if (.not.REALS) IRNCP = IRNCP + 1
!       Set the first non-zero entry in each row to the negative of the
!       row/col number and hold this row/col index in the row/col
!       pointer. This is so that the beginning of each row/col can
!       be recognized in the subsequent scan.
        do 10 J = 1, N
          K = IPTR(J)
          if (K<IACTIV) goto 10
          IPTR(J) = ICN(K)
          ICN(K) = -J
10      end do
        KN = ITOP + 1
        KL = ITOP - IACTIV + 1
!       Go through arrays in reverse order compressing to the back so
!       that there are no zeros held in positions IACTIV to ITOP in ICN.
!       Reset first entry of each row/col and pointer array IPTR.
        do 30 K = 1, KL
          JPOS = ITOP - K + 1
          if (ICN(JPOS)==0) goto 30
          KN = KN - 1
          if (REALS) A(KN) = A(JPOS)
          if (ICN(JPOS)>=0) goto 20
!         First non-zero of row/col has been located.
          J = -ICN(JPOS)
          ICN(JPOS) = IPTR(J)
          IPTR(J) = KN
20        ICN(KN) = ICN(JPOS)
30      end do
        IACTIV = KN
        return

      end subroutine MA30DD
!_______________________________________________________________________

      subroutine MA30BD(N,ICN,A,LICN,LENR,LENRL,IDISP,IP,IQ,W,IW,IFLAG)
! ..
! MA30BD performs the LU decomposition of the diagonal blocks of
! a new matrix PAQ of the same sparsity pattern, using information
! from a previous call to MA30AD. THe entries of the input matrix
! must already be in their final positions in the LU decomposition
! structure. This routine executes about five times faster than
! MA30AD.
! ..
! We now describe the argument list for MA30BD. Consult MA30AD for
!     further information on these parameters.
! N   is an integer variable set to the order of the matrix.
! ICN is an integer array of length LICN. It should be unchanged
!     since the last call to MA30AD. It is not altered by MA30BD.
! A   is a real(kind=wp) array of length LICN. The user must set
!     entries IDISP(1) to IDISP(2) to contain the entries in the
!     diagonal blocks of the matrix PAQ whose column numbers are
!     held in ICN, using corresponding positions. Note that some
!     zeros may need to be held explicitly. On output entries
!     IDISP(1) to IDISP(2) of array A contain the LU decomposition
!     of the diagonal blocks of PAQ. Entries A(1) to A(IDISP(1)-1)
!     are neither required nor altered by MA30BD.
! LICN is an integer variable which must be set by the user to the
!     length of arrays A and ICN. it is not altered by MA30BD.
! LENR, LENRL are integer arrays of length N. They should be
!     unchanged since the last call to MA30AD. They are not
!     altered by MA30BD.
! IDISP is an integer array of length 2. It should be unchanged
!     since the last call to MA30AD. It is not altered by MA30BD.
! IP, IQ are integer arrays of length N. They should be unchanged
!     since the last call to MA30AD. They are not altered by
!     MA30BD.
! W   is a REAL(KIND=WP) array of length N which is used as
!     workspace by MA30BD.
! IW  is an integer array of length N which is used as workspace
!     by MA30BD.
! IFLAG is an integer variable. On output from MA30BD, IFLAG has
!     the value zero if the factorization was successful, has the
!     value I if pivot I was very small and has the value -I if
!     an unexpected singularity was detected at stage I of the
!     decomposition.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer :: IFLAG, LICN, N
! ..
! .. Array Arguments ..
        real (WP) :: A(LICN), W(N)
        integer :: ICN(LICN), IDISP(2), IP(N), IQ(N), IW(N), LENR(N), LENRL(N)
! ..
! .. Local Scalars ..
        real (WP) :: AU, ROWMAX
        integer :: I, IFIN, ILEND, IPIVJ, ISING, ISTART, J, JAY, JAYJAY, JFIN, &
          JJ, PIVPOS
        logical :: STAB
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, MAX
! ..
! .. FIRST EXECUTABLE STATEMENT MA30BD
! ..
        STAB = EPS <= ONE
        RMIN = EPS
        ISING = 0
        IFLAG = 0
        W(1:N) = ZERO
!       Set up pointers to the beginning of the rows.
        IW(1) = IDISP(1)
        if (N==1) goto 30
        do I = 2, N
          IW(I) = IW(I-1) + LENR(I-1)
        end do

!       Start of main loop.
!       At step I, row I of A is transformed to row I of LU by
!       adding appropriate multiples of rows 1 TO I-1 using row
!       Gauss elimination.
30      do 170 I = 1, N
!         ISTART is beginning of row I of A and row I of L.
          ISTART = IW(I)
!         IFIN is end of row I of A and row I of U.
          IFIN = ISTART + LENR(I) - 1
!         ILEND is end of row I of L.
          ILEND = ISTART + LENRL(I) - 1
          if (ISTART>ILEND) goto 100
!         Load row I of A into vector W.
          do JJ = ISTART, IFIN
            J = ICN(JJ)
            W(J) = A(JJ)
          end do

!         Add multiples of appropriate rows of I to I-1 to row I.
          do 80 JJ = ISTART, ILEND
            J = ICN(JJ)
!           IPIVJ is position of pivot in row J.
            IPIVJ = IW(J) + LENRL(J)
!           Form multiplier AU.
            AU = -W(J)/A(IPIVJ)
            if (LBIG) BIG = max(abs(AU),BIG)
            W(J) = AU
!           AU * ROW J (U part) is added to row I.
            IPIVJ = IPIVJ + 1
            JFIN = IW(J) + LENR(J) - 1
            if (IPIVJ>JFIN) goto 80
!           Innermost loop.
            if (LBIG) goto 60
            do JAYJAY = IPIVJ, JFIN
              JAY = ICN(JAYJAY)
              W(JAY) = W(JAY) + AU*A(JAYJAY)
            end do
            goto 80
60          do JAYJAY = IPIVJ, JFIN
              JAY = ICN(JAYJAY)
              W(JAY) = W(JAY) + AU*A(JAYJAY)
              BIG = max(abs(W(JAY)),BIG)
            end do
80        end do

!         Reload W back into A (now LU).
          do JJ = ISTART, IFIN
            J = ICN(JJ)
            A(JJ) = W(J)
            W(J) = ZERO
          end do
!         We now perform the stability checks.
100       PIVPOS = ILEND + 1
          if (IQ(I)>0) goto 150
!         Matrix had singularity at this point in MA30AD.
!         Is it the first such pivot in current block?
          if (ISING==0) ISING = I
!         Does current matrix have a singularity in the same place?
          if (PIVPOS>IFIN) goto 110
          if (abs(A(PIVPOS))>ZERO) goto 180
!         It does, so set ising if it is not the end of the current
!         block. Check to see that appropriate part of LU is zero
!         or null.
110       if (ISTART>IFIN) goto 130
          do 120 JJ = ISTART, IFIN
            if (ICN(JJ)<ISING) goto 120
            if (abs(A(JJ))>ZERO) goto 180
120       end do
130       if (PIVPOS<=IFIN) A(PIVPOS) = ONE
          if (IP(I)>0 .and. I/=N) goto 170
!         End of current block. Reset zero pivots and ISING.
          do 140 J = ISING, I
            if ((LENR(J)-LENRL(J))==0) goto 140
            JJ = IW(J) + LENRL(J)
            A(JJ) = ZERO
140       end do
          ISING = 0
          goto 170
!         Matrix had non-zero pivot in MA30AD at this stage.
150       if (PIVPOS>IFIN) goto 180
          if (abs(A(PIVPOS))<=ZERO) goto 180
          if (.not.STAB) goto 170
          ROWMAX = ZERO
          do JJ = PIVPOS, IFIN
            ROWMAX = max(ROWMAX,abs(A(JJ)))
          end do
          if (abs(A(PIVPOS))/ROWMAX>=RMIN) goto 170
          IFLAG = I
          RMIN = abs(A(PIVPOS))/ROWMAX
!         End of main loop.
170     end do

        goto 190
!       Error return
180     if (LP/=0) write (LP,90000) I
        IFLAG = -I

190     return
90000   format (' Error return from MA30BD; Singularity detected in',' row', &
          I8)

      end subroutine MA30BD
!_______________________________________________________________________

      subroutine MA30CD(N,ICN,A,LICN,LENR,LENRL,LENOFF,IDISP,IP,IQ,X, &
        W,MTYPE)
! ..
! MA30CD uses the factors produced by MA30AD or MA30BD to solve
! AX = B or A TRANSPOSE X = B when the matrix P1*A*Q1(PAQ) is
! block lower triangular (including the case of only one diagonal
! block).
! ..
! We now describe the argument list for MA30CD.
! N   is an integer variable set to the order of the matrix. It is
!     not altered by the subroutine.
! ICN is an integer array of length LICN. Entries IDISP(1) to
!     IDISP(2) should be unchanged since the last call to MA30AD. If
!     the matrix has more than one diagonal block, then column indices
!     corresponding to non-zeros in sub-diagonal blocks of PAQ must
!     appear in positions 1 to IDISP(1)-1. For the same row those
!     entries must be contiguous, with those in row I preceding those
!     in row I+1 (I=1,...,N-1) and no wasted space between rows.
!     Entries may be in any order within each row. It is not altered
!     by MA30CD.
! A   is a REAL(KIND=WP) array of length LICN. Entries IDISP(1) to
!     IDISP(2) should be unchanged since the last call to MA30AD or
!     MA30BD. If the matrix has more than one diagonal block, then
!     the values of the non-zeros in sub-diagonal blocks must be in
!     positions 1 to IDISP(1)-1 in the order given by ICN. It is not
!     altered by MA30CD.
! LICN is an integer variable set to the size of arrays ICN and A.
!     It is not altered by MA30CD.
! LENR, LENRL are integer arrays of length N which should be
!     unchanged since the last call to MA30AD. They are not
!     altered by MA30CD.
! LENOFF is an integer array of length N. If the matrix PAQ (or
!     P1*A*Q1) has more than one diagonal block, then LENOFF(I),
!     I=1,...,N should be set to the number of non-zeros in row I of
!     the matrix PAQ which are in sub-diagonal blocks. If there is
!     only one diagonal block then LENOFF(1) may be set TO -1, in
!     which case the other entries of LENOFF are never accessed.
!     It is not altered by ma30cd.
! IDISP is an integer array of length 2 which should be unchanged
!     since the last call to MA30AD. It is not altered by MA30CD.
! IP, IQ are integer arrays of length N which should be unchanged
!     since the last call to MA30AD. They are not altered by
!     MA30CD.
! X   is a REAL(KIND=WP) array of length N. It must be set by the
!     user to the values of the right hand side vector B for the
!     equations being solved. ON exit from MA30CD it will be equal
!     to the solution X required.
! W   is a REAL(KIND=WP) array of length N which is used as
!     workspace by MA30CD.
! MTYPE is an integer variable which must be set by the user. If
!     MTYPE = 1, then the solution to the system AX = B is returned.
!     Any other value for mtype will return the solution to the
!     system A TRANSPOSE X = B. It is not altered by MA30CD.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer :: LICN, MTYPE, N
! ..
! .. Array Arguments ..
        real (WP) :: A(LICN), W(N), X(N)
        integer :: ICN(LICN), IDISP(2), IP(N), IQ(N), LENOFF(N), LENR(N), &
          LENRL(N)
! ..
! .. Local Scalars ..
        real (WP) :: WI, WII
        integer :: I, IB, IBACK, IBLEND, IBLOCK, IEND, IFIRST, II, III,   &
          ILAST, J, J1, J2, J3, JJ, JPIV, JPIVP1, K, LJ1, LJ2, LT, LTEND, &
          NUMBLK
        logical :: NEG, NOBLOC
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, IABS, MAX
! ..
! .. FIRST EXECUTABLE STATEMENT MA30CD
! ..
!       The final value of RESID is the maximum residual for an
!       inconsistent set of equations.
        RESID = ZERO
!       NOBLOC is .TRUE. if subroutine block has been used
!       previously and is .FALSE. otherwise. The value .FALSE.
!       means that LENOFF will not be subsequently accessed.
        NOBLOC = LENOFF(1) < 0
        if (MTYPE/=1) goto 140

!       We now solve  A * X = B.
!       NEG is used to indicate when the last row in a block
!       has been reached. It is then set to .TRUE. whereafter
!       back substitution is performed on the block.
        NEG = .false.
!       IP(N) is negated so that the last row of the last
!       block can be recognised. It is reset to its positive
!       value on exit.
        IP(N) = -IP(N)
!       Preorder VECTOR ... W(I) = X(IP(I))
        do II = 1, N
          I = IP(II)
          I = IABS(I)
          W(II) = X(I)
        end do
!       LT holds the position of the first non-zero in the current
!       row of the off-diagonal blocks.
        LT = 1
!       IFIRST holds the index of the first row in the current block.
        IFIRST = 1
!       IBLOCK holds the position of the first non-zero in the current
!       row of the LU decomposition of the diagonal blocks.
        IBLOCK = IDISP(1)
!       If I is not the last row of a block, then a pass through this
!       loop adds the inner product of row I of the off-diagonal blocks
!       and W to W and performs forward elimination using row I of the
!       LU decomposition. If I is the last row of a block then, after
!       performing these aforementioned operations, back substitution
!       is performed using the rows of the block.
        do 120 I = 1, N
          WI = W(I)
          if (NOBLOC) goto 30
          if (LENOFF(I)==0) goto 30
!         Operations using lower triangular blocks.
!         LTEND is the end of row I in the off-diagonal blocks.
          LTEND = LT + LENOFF(I) - 1
          do JJ = LT, LTEND
            J = ICN(JJ)
            WI = WI - A(JJ)*W(J)
          end do
!         LT is set the beginning of the next off-diagonal row.
          LT = LTEND + 1
!         Set NEG to .TRUE. if we are on the last row of the block.
30        if (IP(I)<0) NEG = .true.
          if (LENRL(I)==0) goto 50
!         Forward elimination phase.
!         IEND is the end of the L part of row I in the LU decomposition.
          IEND = IBLOCK + LENRL(I) - 1
          do JJ = IBLOCK, IEND
            J = ICN(JJ)
            WI = WI + A(JJ)*W(J)
          end do
!         IBLOCK is adjusted to point to the start of the next row.
50        IBLOCK = IBLOCK + LENR(I)
          W(I) = WI
          if (.not.NEG) goto 120
!         Back substitution phase.
!         J1 is position in A/ICN after end of block beginning in
!         row IFIRST and ending in row I.
          J1 = IBLOCK
!         Are there any singularities in this block? If not, continue
!         with the back substitution.
          IB = I
          if (IQ(I)>0) goto 70
          do III = IFIRST, I
            IB = I - III + IFIRST
            if (IQ(IB)>0) goto 70
            J1 = J1 - LENR(IB)
            RESID = max(RESID,abs(W(IB)))
            W(IB) = ZERO
          end do
!         Entire block is singular.
          goto 110
!         Each pass through this loop performs the back substitution
!         operations for a single row, starting at the end of the
!         block and working through it in reverse order.
70        do III = IFIRST, IB
            II = IB - III + IFIRST
!           J2 is end of row II.
            J2 = J1 - 1
!           J1 is beginning of row II.
            J1 = J1 - LENR(II)
!           JPIV is the position of the pivot in row II.
            JPIV = J1 + LENRL(II)
            JPIVP1 = JPIV + 1
!           JUMP if row II of U has no non-zeros.
            if (J2<JPIVP1) goto 90
            WII = W(II)
            do JJ = JPIVP1, J2
              J = ICN(JJ)
              WII = WII - A(JJ)*W(J)
            end do
            W(II) = WII
90          W(II) = W(II)/A(JPIV)
          end do
110       IFIRST = I + 1
          NEG = .false.
120     end do

!       Reorder solution vector, X(I) = W(IQINVERSE(I)).
        do II = 1, N
          I = IQ(II)
          I = IABS(I)
          X(I) = W(II)
        end do
        IP(N) = -IP(N)
        goto 320

!       WE now solve ATRANSPOSE * X = B.
!       Preorder vector, W(I) = X(IQ(I)).
140     do II = 1, N
          I = IQ(II)
          I = IABS(I)
          W(II) = X(I)
        end do
!       LJ1 points to the beginning the current row in the off-diagonal
!       blocks.
        LJ1 = IDISP(1)
!       IBLOCK is initialized to point to the beginning of the block
!       after the last one.
        IBLOCK = IDISP(2) + 1
!       ILAST is the last row in the current block.
        ILAST = N
!       IBLEND points to the position after the last non-zero in the
!       current block.
        IBLEND = IBLOCK
!       Each pass through this loop operates with one diagonal block and
!       the off-diagonal part of the matrix corresponding to the rows
!       of this block. The blocks are taken in reverse order and the
!       number of times the loop is entered is MIN(N, NO. BLOCKS+1).
        do NUMBLK = 1, N
          if (ILAST==0) goto 300
          IBLOCK = IBLOCK - LENR(ILAST)
!         This loop finds the index of the first row in the current
!         block. It is FIRST and IBLOCK is set to the position of
!         the beginning of this first row.
          do K = 1, N
            II = ILAST - K
            if (II==0) goto 170
            if (IP(II)<0) goto 170
            IBLOCK = IBLOCK - LENR(II)
          end do
170       IFIRST = II + 1
!         J1 points to the position of the beginning of row I (LT part)
!         or pivot.
          J1 = IBLOCK
!         Forward elimination.
!         Each pass through this loop performs the operations for one row
!         of the block. If the corresponding entry of W is zero then the
!         operations can be avoided.
          do I = IFIRST, ILAST
            if (abs(W(I))<=ZERO) goto 200
!           JUMP IF ROW I SINGULAR.
            if (IQ(I)<0) goto 220
!           J2 first points to the pivot in row I and then is made
!           to point to the first non-zero in the U TRANSPOSE part
!           of the row.
            J2 = J1 + LENRL(I)
            WI = W(I)/A(J2)
            if (LENR(I)-LENRL(I)==1) goto 190
            J2 = J2 + 1
!           J3 points to the end of row I.
            J3 = J1 + LENR(I) - 1
            do JJ = J2, J3
              J = ICN(JJ)
              W(J) = W(J) - A(JJ)*WI
            end do
190         W(I) = WI
200         J1 = J1 + LENR(I)
          end do
          goto 240
!         Deal with rest of block which is singular.
220       do II = I, ILAST
            RESID = max(RESID,abs(W(II)))
            W(II) = ZERO
          end do
!         Back substitution.
!         This loop does the back substitution on the rows of the
!         block in the reverse order doing it simultaneously on
!         the L TRANSPOSE part of the diagonal blocks and the
!         off-diagonal blocks.
240       J1 = IBLEND
          do 280 IBACK = IFIRST, ILAST
            I = ILAST - IBACK + IFIRST
!           J1 points to the beginning of row I.
            J1 = J1 - LENR(I)
            if (LENRL(I)==0) goto 260
!           J2 points to the end of the L TRANSPOSE part of row I.
            J2 = J1 + LENRL(I) - 1
            do JJ = J1, J2
              J = ICN(JJ)
              W(J) = W(J) + A(JJ)*W(I)
            end do
260         if (NOBLOC) goto 280
!           Operations using lower triangular blocks.
            if (LENOFF(I)==0) goto 280
!           LJ2 points to the end of row I of the off-diagonal blocks.
            LJ2 = LJ1 - 1
!           LJ1 points to the beginning of row I of the off-diagonal
!           blocks.
            LJ1 = LJ1 - LENOFF(I)
            do JJ = LJ1, LJ2
              J = ICN(JJ)
              W(J) = W(J) - A(JJ)*W(I)
            end do
280       end do
          IBLEND = J1
          ILAST = IFIRST - 1
        end do
!       Reorder solution vector, X(I) = W(IPINVERSE(I)).
300     do II = 1, N
          I = IP(II)
          I = IABS(I)
          X(I) = W(II)
        end do

320     return

      end subroutine MA30CD
!_______________________________________________________________________

! Private Variable Information.
! Private variables for MA30ED hold control parameters.
!     Original : COMMON /MA30ED/ LP,ABORT1,ABORT2,ABORT3
! The integer LP is the unit number to which the error messages are
!     sent. LP has a default value of 6. This default value can be
!     reset by the user, if desired. A value of 0 suppresses all
!     messages.
! The logical variables ABORT1, ABORT2, ABORT3 are used to control
!     the conditions under which the subroutine will terminate.
! If ABORT1 iS .TRUE. then the subroutine will exit immediately on
!     detecting structural singularity.
! If ABORT2 iS .TRUE. then the subroutine will exit immediately on
!     detecting numerical singularity.
! If ABORT3 iS .TRUE. then the subroutine will exit immediately when
!     the available space in A/ICN is filled up by the previously
!     decomposed, active, and undecomposed parts of the matrix.
! The default values for ABORT1, ABORT2, ABORT3 are set to .TRUE.,
!     .TRUE., and .FALSE., respectively.

! The private variables for MA30FD are used to provide the user with
! information on the decomposition.
!     Original : COMMON /MA30FD/ IRNCP,ICNCP,IRANK,MINIRN,MINICN
! IRNCP AND ICNCP are integer variables used to monitor the adequacy
!     of the allocated space in arrays IRN and A/ICN, respectively,
!     by taking account of the number of data management compresses
!     required on these arrays. If IRNCP or ICNCP is fairly large (say
!     greater than N/10), It may be advantageous to increase the size
!     of the corresponding array(s). IRNCP and ICNCP are initialized
!     to zero on entry to MA30AD and are incremented each time the
!     compressing routine MA30DD is entered.
! ICNCP is the number of compresses on A/ICN.
! IRNCP is the number of compresses on IRN.
! IRANK is an integer variable which gives an estimate (actually an
!     upper bound) of the rank of the matrix. On an exit with IFLAG
!     equal to 0, this will be equal to N.
! MINIRN is an integer variable which, after a successful call to
!     MA30AD, indicates the minimum length to which IRN can be
!     reduced while still permitting a successful decomposition
!     of the same matrix. If, however, the user were to decrease the
!     length of IRN to that size, the number of compresses (IRNCP)
!     may be very high and quite costly. If LIRN is not large enough
!     to begin the decomposition on a diagonal block, MINIRN will be
!     equal to the value required to continue the decomposition and
!     IFLAG will be set to -3 or -6. A value of LIRN slightly greater
!     than this (say about N/2) will usually provide enough space to
!     complete the decomposition on that block. In the event of any
!     other failure minirn gives the minimum size of IRN required
!     for a successful decomposition up to that point.
! MINICN is an integer variable which after a successful call to
!     MA30AD, indicates the minimum size of LICN required to enable
!     a successful decomposition. In the event of failure with
!     IFLAG = -5, MINICN will, if ABORT3 is left set to .FALSE.,
!     indicate the minimum length that would be sufficient to
!     prevent this error in a subsequent run on an identical matrix.
!     Again the user may prefer to use a value of icn slightly
!     greater than MINICN for subsequent runs to avoid too many
!     conpresses (ICNCP). In the event of failure with IFLAG equal
!     to any negative value except -4, minicn will give the minimum
!     length to which licn could be reduced to enable a successful
!     decomposition to the point at which failure occurred. Notice
!     that, on a successful entry IDISP(2) gives the amount of space
!     in A/ICN required for the decomposition while MINICN will
!     usually be slightly greater because of the need for
!     "elbow room". If the user is very unsure how large to make
!     LICN, the variable MINICN can be used to provide that
!     information. A preliminary run should be performed with ABORT3
!     left set to .FALSE. and LICN about 3/2 times as big as the
!     number of non-zeros in the original matrix. Unless the initial
!     problem is very sparse (when the run will be successful) or
!     fills in extremely badly (giving an error return with
!     IFLAG = -4), an error return with IFLAG = -5 should result and
!     MINICN will give the amount of space required for a successful
!     decomposition.

! The private variables for MA30GD are used by the MA30BD entry only.
!     Original : COMMON /MA30GD/ EPS,RMIN
! EPS is a REAL(KIND=WP) variable. It is used to test for small
!     pivots. Its default value is 1.0D-4 (1.0E-4 in E version).
!     If the user sets EPS to any value greater than 1.0, then no
!     check is made on the size of the pivots. Although the absence
!     of such a check would fail to warn the user of bad instability,
!     its absence will enable MA30BD to run slightly faster. An a
!     posteriori check on the stability of the factorization can
!     be obtained from MC24AD.
! RMIN is a REAL(KIND=WP) variable which gives the user some
!     information about the stability of the decomposition. At each
!     stage of the LU decomposition the magnitude of the pivot APIV
!     is compared with the largest off-diagonal entry currently in
!     its row (row of U), ROWMAX say. If the ratio
!                       MIN(APIV/ROWMAX)
!     where the minimum is taken over all the rows, is less than EPS
!     then RMIN is set to this minimum value and IFLAG is returned
!     with the value +I where I is the row in which this minimum
!     occurs. If the user sets EPS greater than one, then this test
!     is not performed. In this case, and when there are no small
!     pivots rmin will be set equal to EPS.

! The private variables for MA30HD are used by MA30CD only.
!     Original : COMMON /MA30HD/ RESID
! RESID is a REAL(KIND=WP) variable. In the case of singular or
!     rectangular matrices its final value will be equal to the
!     maximum residual for the unsatisfied equations; otherwise
!     its value will be set to zero.

!  MA30ID private variables control the use of drop tolerances,
!     the modified pivot option and the calculation of the largest
!     entry in the factorization process. These private variables
!     were added to the MA30 package in February, 1983.
!     Original : COMMON /MA30ID/ TOL,BIG,NDROP,NSRCH,LBIG
! TOL is a REAL(KIND=WP) variable. If it is set to a positive
!     value, then MA30AD will drop from the factors any non-zero
!     whose modulus is less than TOL. The factorization will then
!     require less storage but will be inaccurate. After a run of
!     MA30AD where entries have been dropped, MA30BD should not
!     be called. The default value for TOL is 0.0.
! BIG is a REAL(KIND=WP) variable. If LBIG has been set to .TRUE.,
!     BIG will be set to the largest entry encountered during the
!     factorization.
! NDROP is an integer variable. If TOL has been set positive, on
!     exit from MA30AD, NDROP will hold the number of entries
!     dropped from the data structure.
! NSRCH is an integer variable. If NSRCH is set to a value less
!     than or equal to N, then a different pivot option will be
!     employed by MA30AD. This may result in different fill-in
!     and execution time for MA30AD. If NSRCH is less than or
!     equal to N, the workspace arrays LASTC and NEXTC are not
!     referenced by MA30AD. The default value for NSRCH is 32768.
! LBIG is a logical variable. If LBIG is set to .TRUE., the value
!     of the largest entry encountered in the factorization by
!     MA30AD is returned in BIG. Setting LBIG to .TRUE. will
!     marginally increase the factorization time for MA30AD and
!     will increase that for MA30BD by about 20%. The default
!     value for LBIG is .FALSE.
!_______________________________________________________________________

      subroutine MC13D(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW)
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer :: LICN, N, NUM
! ..
! .. Array Arguments ..
        integer :: IB(N), ICN(LICN), ior(N), IP(N), IW(N,3), LENR(N)
! ..
! .. FIRST EXECUTABLE STATEMENT MA13D
! ..
        call MC13E(N,ICN,LICN,IP,LENR,IOR,IB,NUM,IW(1,1),IW(1,2),IW(1,3))
        return

      end subroutine MC13D
!_______________________________________________________________________

      subroutine MC13E(N,ICN,LICN,IP,LENR,ARP,IB,NUM,LOWL,NUMB,PREV)
! ..
! ARP(I) is one less than the number of unsearched edges leaving
!     node I. At the end of the algorithm it is set to a
!     permutation which puts the matrix in block lower
!     triangular form.
! IB(I) is the position in the ordering of the start of the Ith
!     block. IB(N+1-I) holds the node number of the Ith node
!     on the stack.
! LOWL(I) is the smallest stack position of any node to which a
!     path from node I has been found. It is set to N+1 when
!     node I is removed from the stack.
! NUMB(I) is the position of node I in the stack if it is on
!     the stack. It is the permuted order of node I for those
!     nodes whose final position has been found and is otherwise
!     zero.
! PREV(I) is the node at the end of the path when node I was
!     placed on the stack.
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer :: LICN, N, NUM
! ..
! .. Array Arguments ..
        integer :: ARP(N), IB(N), ICN(LICN), IP(N), LENR(N), LOWL(N),  &
          NUMB(N), PREV(N)
! ..
! .. Local Scalars ..
        integer :: DUMMY, I, I1, I2, ICNT, II, ISN, IST, IST1, IV, IW, &
          J, K, LCNT, NNM1, STP
! ..
! .. Intrinsic Functions ..
        intrinsic MIN
! ..
! .. FIRST EXECUTABLE STATEMENT MA13E
! ..
!       ICNT is the number of nodes whose positions in final ordering
!       have been found.
        ICNT = 0
!       NUM is the number of blocks that have been found.
        NUM = 0
        NNM1 = N + N - 1

!       Initialization of arrays.
        do J = 1, N
          NUMB(J) = 0
          ARP(J) = LENR(J) - 1
        end do

        do 90 ISN = 1, N
!         Look for a starting node.
          if (NUMB(ISN)/=0) goto 90
          IV = ISN
!         IST is the number of nodes on the stack. It is the
!         stack pointer.
          IST = 1
!         Put node IV at beginning of stack.
          LOWL(IV) = 1
          NUMB(IV) = 1
          IB(N) = IV

!         The body of this loop puts a new node on the stack
!         or backtracks.
          do 80 DUMMY = 1, NNM1
            I1 = ARP(IV)
!           Have all edges leaving node iv been searched?
            if (I1<0) goto 30
            I2 = IP(IV) + LENR(IV) - 1
            I1 = I2 - I1

!           Look at edges leaving node iv until one enters a
!           new node or all edges are exhausted.
            do 20 II = I1, I2
              IW = ICN(II)
!             Has node iw been on stack already?
              if (NUMB(IW)==0) goto 70
!             Update value of LOWL(IV) if necessary.
20          LOWL(IV) = min(LOWL(IV),LOWL(IW))

!           There are no more edges leaving node IV.
            ARP(IV) = -1
!           Is node IV the root of a block.
30          if (LOWL(IV)<NUMB(IV)) goto 60

!           Order nodes in a block.
            NUM = NUM + 1
            IST1 = N + 1 - IST
            LCNT = ICNT + 1
!           Peel block off the top of the stack starting at the
!           top and working down to the root of the block.
            do STP = IST1, N
              IW = IB(STP)
              LOWL(IW) = N + 1
              ICNT = ICNT + 1
              NUMB(IW) = ICNT
              if (IW==IV) goto 50
            end do
50          IST = N - STP
            IB(NUM) = LCNT
!           Are there any nodes left on the stack?
            if (IST/=0) goto 60
!           Have all the nodes been ordered?
            if (ICNT<N) goto 90
            goto 100

!           Backtrack to previous node on path.
60          IW = IV
            IV = PREV(IV)
!           Update value of LOWL(IV) if necessary.
            LOWL(IV) = min(LOWL(IV),LOWL(IW))
            goto 80

!           Put new node on the stack.
70          ARP(IV) = I2 - II - 1
            PREV(IW) = IV
            IV = IW
            IST = IST + 1
            LOWL(IV) = IST
            NUMB(IV) = IST
            K = N + 1 - IST
            IB(K) = IV
80        end do

90      end do

!       Put permutation in the required form.
100     do I = 1, N
          II = NUMB(I)
          ARP(II) = I
        end do
        return

      end subroutine MC13E
!_______________________________________________________________________

      subroutine MC20AD(NC,MAXA,A,INUM,JPTR,JNUM,JDISP)
! ..
     implicit none
! ..
! .. Scalar Arguments ..
! ..
        integer :: JDISP, MAXA, NC
! ..
! .. Array Arguments ..
        real (WP) :: A(MAXA)
        integer :: INUM(MAXA), JNUM(MAXA), JPTR(NC)
! ..
! .. Local Scalars ..
        real (WP) :: ACE, ACEP
        integer :: I, ICE, ICEP, J, JA, JB, JCE, JCEP, K, KR, LOC, MYNULL
! ..
! .. FIRST EXECUTABLE STATEMENT MA20AD
! ..
        MYNULL = -JDISP
!       CLEAR JPTR
        JPTR(1:NC) = 0
!       Count the number of elements in each column.
        do K = 1, MAXA
          J = JNUM(K) + JDISP
          JPTR(J) = JPTR(J) + 1
        end do
!       SET THE JPTR ARRAY.
        K = 1
        do J = 1, NC
          KR = K + JPTR(J)
          JPTR(J) = K
          K = KR
        end do

!       Reorder the elements into column order. The algorithm is an
!       in-place sort and is of order MAXA.
        do 50 I = 1, MAXA
!         Establish the current entry.
          JCE = JNUM(I) + JDISP
          if (JCE==0) goto 50
          ACE = A(I)
          ICE = INUM(I)
!         Clear the location vacated.
          JNUM(I) = MYNULL
!         Chain from current entry to store items.
          do 40 J = 1, MAXA
!           Current entry not in correct position. Determine the
!           correct position to store entry.
            LOC = JPTR(JCE)
            JPTR(JCE) = JPTR(JCE) + 1
!           Save contents of that location.
            ACEP = A(LOC)
            ICEP = INUM(LOC)
            JCEP = JNUM(LOC)
!           Store current entry.
            A(LOC) = ACE
            INUM(LOC) = ICE
            JNUM(LOC) = MYNULL
!           Check if next current entry needs to be processed.
            if (JCEP==MYNULL) goto 50
!           It does. Copy into current entry.
            ACE = ACEP
            ICE = ICEP
40        JCE = JCEP + JDISP

50      end do

!       Reset JPTR vector.
        JA = 1
        do J = 1, NC
          JB = JPTR(J)
          JPTR(J) = JA
          JA = JB
        end do
        return

      end subroutine MC20AD
!_______________________________________________________________________

      subroutine MC20BD(NC,MAXA,A,INUM,JPTR)
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer :: MAXA, NC
! ..
! .. Array Arguments ..
        real (WP) :: A(MAXA)
        integer :: INUM(MAXA), JPTR(NC)
! ..
! .. Local Scalars ..
        real (WP) :: ACE
        integer :: ICE, IK, J, JJ, K, KDUMMY, KLO, KMAX, KOR
! ..
! .. Intrinsic Functions ..
        intrinsic IABS
! ..
! .. FIRST EXECUTABLE STATEMENT MA20BD
! ..
        KMAX = MAXA
        do JJ = 1, NC
          J = NC + 1 - JJ
          KLO = JPTR(J) + 1
          if (KLO>KMAX) goto 40
          KOR = KMAX
          do KDUMMY = KLO, KMAX
!           Items KOR,KOR+1,...,KMAX are in order.
            ACE = A(KOR-1)
            ICE = INUM(KOR-1)
            do K = KOR, KMAX
              IK = INUM(K)
              if (IABS(ICE)<=IABS(IK)) goto 20
              INUM(K-1) = IK
              A(K-1) = A(K)
            end do
            K = KMAX + 1
20          INUM(K-1) = ICE
            A(K-1) = ACE
            KOR = KOR - 1
          end do
!         Next column.
40        KMAX = KLO - 2
        end do
        return

      end subroutine MC20BD
!_______________________________________________________________________

      subroutine MC21A(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer :: LICN, N, NUMNZ
! ..
! .. Array Arguments ..
        integer :: ICN(LICN), IP(N), IPERM(N), IW(N,4), LENR(N)
! ..
! .. FIRST EXECUTABLE STATEMENT MA21A
! ..
        call MC21B(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2), &
          IW(1,3),IW(1,4))
        return

      end subroutine MC21A
!_______________________________________________________________________

      subroutine MC21B(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
! ..
     implicit none
! ..
! .. Scalar Arguments ..
! ..
        integer :: LICN, N, NUMNZ
! ..
! .. Array Arguments ..
        integer :: ARP(N), CV(N), ICN(LICN), IP(N), IPERM(N), LENR(N), &
          OUT(N), PR(N)
! ..
! .. Local Scalars ..
        integer :: I, II, IN1, IN2, IOUTK, J, J1, JORD, K, KK
! ..
! .. FIRST EXECUTABLE STATEMENT MA21B
! ..
! PR(I) is the previous row to I in the depth first search.
!    It is used as a work array in the sorting algorithm.
!    Elements (IPERM(I),I) I=1,...,N are non-zero at the
!    end of the algorithm unless N assignments have not
!    been made, in which case (IPERM(I),I) will be zero
!    for N-NUMNZ entries.
! CV(I) is the most recent row extension at which column I
!    was visited.
! ARP(I) is one less than the number of non-zeros in row I
!    which have not been scanned when looking for a cheap
!    assignment.
! OUT(I) is one less than the number of non-zeros in row I
!    which have not been scanned during one pass through
!    the main loop.

!       Initialization of arrays.
        do I = 1, N
          ARP(I) = LENR(I) - 1
          CV(I) = 0
          IPERM(I) = 0
        end do
        NUMNZ = 0

!       Main loop.
!       Each pass round this loop either results in a new
!       assignment or gives a row with no assignment.
        do 100 JORD = 1, N
          J = JORD
          PR(J) = -1
          do 70 K = 1, JORD
!           Look for a cheap assignment.
            IN1 = ARP(J)
            if (IN1<0) goto 30
            IN2 = IP(J) + LENR(J) - 1
            IN1 = IN2 - IN1
            do II = IN1, IN2
              I = ICN(II)
              if (IPERM(I)==0) goto 80
            end do
!           No cheap assignment in row.
            ARP(J) = -1
!           Begin looking for assignment chain starting with row J.
30          OUT(J) = LENR(J) - 1
!           Inner loop. Extends chain by one or backtracks.
            do KK = 1, JORD
              IN1 = OUT(J)
              if (IN1<0) goto 50
              IN2 = IP(J) + LENR(J) - 1
              IN1 = IN2 - IN1
!             Forward scan.
              do 40 II = IN1, IN2
                I = ICN(II)
                if (CV(I)==JORD) goto 40
!               Column I has not yet been accessed during this pass.
                J1 = J
                J = IPERM(I)
                CV(I) = JORD
                PR(J) = J1
                OUT(J1) = IN2 - II - 1
                goto 70
40            end do

!             Backtracking step.
50            J = PR(J)
              if (J==-1) goto 100
            end do

70        end do

!         New assignment is made.
80        IPERM(I) = J
          ARP(J) = IN2 - II - 1
          NUMNZ = NUMNZ + 1
          do K = 1, JORD
            J = PR(J)
            if (J==-1) goto 100
            II = IP(J) + LENR(J) - OUT(J) - 2
            I = ICN(II)
            IPERM(I) = J
          end do

100     end do

!       If matrix is structurally singular, we now complete the
!       permutation IPERM.
        if (NUMNZ==N) return
        ARP(1:N) = 0
        K = 0
        do 130 I = 1, N
          if (IPERM(I)/=0) goto 120
          K = K + 1
          OUT(K) = I
          goto 130
120       J = IPERM(I)
          ARP(J) = I
130     end do
        K = 0
        do 140 I = 1, N
          if (ARP(I)/=0) goto 140
          K = K + 1
          IOUTK = OUT(K)
          IPERM(IOUTK) = I
140     end do
        return

      end subroutine MC21B
!_______________________________________________________________________

      subroutine MC22AD(N,ICN,A,NZ,LENROW,IP,IQ,IW,IW1)
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer :: N, NZ
! ..
! .. Array Arguments ..
        real (WP) :: A(NZ)
        integer :: ICN(NZ), IP(N), IQ(N), IW(N,2), IW1(NZ), LENROW(N)
! ..
! .. Local Scalars ..
        real (WP) :: AVAL
        integer :: I, ICHAIN, IOLD, IPOS, J, J2, JJ, JNUM, JVAL, LENGTH, &
          NEWPOS
! ..
! .. Intrinsic Functions ..
        intrinsic IABS
! ..
! .. FIRST EXECUTABLE STATEMENT MA22AD
! ..
        if (NZ<=0) goto 90
        if (N<=0) goto 90
!       Set start of row I in IW(I,1) and LENROW(I) in IW(I,2)
        IW(1,1) = 1
        IW(1,2) = LENROW(1)
        do I = 2, N
          IW(I,1) = IW(I-1,1) + LENROW(I-1)
          IW(I,2) = LENROW(I)
        end do
!       Permute LENROW according to IP. Set off-sets for new
!       position of row IOLD in IW(IOLD,1) and put old row
!       indices in IW1 in positions corresponding to the new
!       position of this row in A/ICN.
        JJ = 1
        do 30 I = 1, N
          IOLD = IP(I)
          IOLD = IABS(IOLD)
          LENGTH = IW(IOLD,2)
          LENROW(I) = LENGTH
          if (LENGTH==0) goto 30
          IW(IOLD,1) = IW(IOLD,1) - JJ
          J2 = JJ + LENGTH - 1
          do 20 J = JJ, J2
20        IW1(J) = IOLD
          JJ = J2 + 1
30      end do
!       Set inverse permutation to IQ in IW(:,2).
        do I = 1, N
          IOLD = IQ(I)
          IOLD = IABS(IOLD)
          IW(IOLD,2) = I
        end do
!       Permute A and ICN in place, changing to new column numbers.
!       Main loop.
!       Each pass through this loop places a closed chain of
!       column indices in their new (and final) positions.
!       This is recorded by setting the iw1 entry to zero so
!       that any which are subsequently encountered during
!       this major scan can be bypassed.
        do 80 I = 1, NZ
          IOLD = IW1(I)
          if (IOLD==0) goto 80
          IPOS = I
          JVAL = ICN(I)
!         If row IOLD is in same positions after permutation,
!         GOTO 150.
          if (IW(IOLD,1)==0) goto 70
          AVAL = A(I)
!         Chain loop.
!         Each pass through this loop places one(permuted) column
!         index in its final position, viz. IPOS.
          do ICHAIN = 1, NZ
!           NEWPOS is the original position in A/ICN of the
!           element to be placed in position IPOS. It is also
!           the position of the next element in the chain.
            NEWPOS = IPOS + IW(IOLD,1)
!           Is chain complete?
            if (NEWPOS==I) goto 60
            A(IPOS) = A(NEWPOS)
            JNUM = ICN(NEWPOS)
            ICN(IPOS) = IW(JNUM,2)
            IPOS = NEWPOS
            IOLD = IW1(IPOS)
            IW1(IPOS) = 0
!           End of chain loop.
          end do
60        A(IPOS) = AVAL
70        ICN(IPOS) = IW(JVAL,2)
!         END OF MAIN LOOP
80      end do
90      return

      end subroutine MC22AD
!_______________________________________________________________________

      subroutine MC23AD(N,ICN,A,LICN,LENR,IDISP,IP,IQ,LENOFF,IW,IW1)
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer :: LICN, N
! ..
! .. Array Arguments ..
        real (WP) :: A(LICN)
        integer :: ICN(LICN), IDISP(2), IP(N), IQ(N), IW(N,5), IW1(N,2), &
          LENOFF(N), LENR(N)
! ..
! .. Local Scalars ..
        integer :: I, I1, I2, IBEG, IBLOCK, IEND, II, ILEND, INEW, IOLD, &
          IROWB, IROWE, J, JJ, JNEW, JNPOS, JOLD, K, LENI, NZ
! ..
! .. Intrinsic Functions ..
        intrinsic MAX, MIN
! ..
! .. FIRST EXECUTABLE STATEMENT MA23AD
! ..
!     Input ... N,ICN ... A,ICN,LENR ...
!     Set up pointers IW(:,1) to the beginning of the rows
!     and set LENOFF equal to LENR.
        IW1(1,1) = 1
        LENOFF(1) = LENR(1)
        if (N==1) goto 20
        do I = 2, N
          LENOFF(I) = LENR(I)
          IW1(I,1) = IW1(I-1,1) + LENR(I-1)
        end do
!       IDISP(1) points to the first position in A/ICN after
!       the off-diagonal blocks and untreated rows.
20      IDISP(1) = IW1(N,1) + LENR(N)

!       Find row permutation ip to make diagonal zero-free.
        call MC21A(N,ICN,LICN,IW1(1,1),LENR,IP,NUMNZ,IW(1,1))

!       Possible error return for structurally singular matrices.
        if (NUMNZ/=N .and. ABORT) goto 170

!       IW1(:,2) and LENR are permutations of IW1(:,1) and
!       LENR/LENOFF suitable for entry to MC13D since matrix
!       with these row pointer and length arrays has maximum
!       number of non-zeros on the diagonal.
        do II = 1, N
          I = IP(II)
          IW1(II,2) = IW1(I,1)
          LENR(II) = LENOFF(I)
        end do

!       Find symmetric permutation IQ to block lower triangular form.
        call MC13D(N,ICN,LICN,IW1(1,2),LENR,IQ,IW(1,4),NUM,IW)

        if (NUM/=1) goto 60

!       Action taken if matrix is irreducible: The
!       whole matrix is just moved to the end of the storage.
        do I = 1, N
          LENR(I) = LENOFF(I)
          IP(I) = I
          IQ(I) = I
        end do
        LENOFF(1) = -1
!       IDISP(1) is the first position after the last element in
!       the off-diagonal blocks and untreated rows.
        NZ = IDISP(1) - 1
        IDISP(1) = 1
!       IDISP(2) Is the position in A/ICN of the first element
!       in the diagonal blocks.
        IDISP(2) = LICN - NZ + 1
        LARGE = N
        if (NZ==LICN) goto 200
        do K = 1, NZ
          J = NZ - K + 1
          JJ = LICN - K + 1
          A(JJ) = A(J)
          ICN(JJ) = ICN(J)
        end do
        goto 200

!       Data structure reordered.

!       Form composite row permutation, IP(I) = IP(IQ(I)).
60      do II = 1, N
          I = IQ(II)
          IW(II,1) = IP(I)
        end do
        do I = 1, N
          IP(I) = IW(I,1)
        end do

!       Run through blocks in reverse order separating diagonal
!       blocks which are moved to the end of the storage. Elements
!       in off-diagonal blocks are left in place unless a compress
!       is necessary.

!       IBEG indicates the lowest value of J for which ICN(J) has
!       been set to zero when element in position J was moved to
!       the diagonal block part of storage.
        IBEG = LICN + 1
!       IEND is the position of the first element of those treated
!       rows which are in diagonal blocks.
        IEND = LICN + 1
!       LARGE is the dimension of the largest block encountered
!       so far.
        LARGE = 0

!       NUM is the number of diagonal blocks.
        do K = 1, NUM
          IBLOCK = NUM - K + 1
!         I1 is first row (in permuted form) of block IBLOCK.
!         I2 is last row (in permuted form) of block IBLOCK.
          I1 = IW(IBLOCK,4)
          I2 = N
          if (K/=1) I2 = IW(IBLOCK+1,4) - 1
          LARGE = max(LARGE,I2-I1+1)
!         Go through the rows of block IBLOCK in the reverse order.
          do II = I1, I2
            INEW = I2 - II + I1
!           We now deal with row inew in permuted form (row IOLD
!           in original matrix).
            IOLD = IP(INEW)
!           If there is space to move up diagonal block portion
!           of row GOTO 110.
            if (IEND-IDISP(1)>=LENOFF(IOLD)) goto 110

!           In-line compress.
!           Moves separated off-diagonal elements and untreated
!           rows to front of storage.
            JNPOS = IBEG
            ILEND = IDISP(1) - 1
            if (ILEND<IBEG) goto 180
            do 90 J = IBEG, ILEND
              if (ICN(J)==0) goto 90
              ICN(JNPOS) = ICN(J)
              A(JNPOS) = A(J)
              JNPOS = JNPOS + 1
90          end do
            IDISP(1) = JNPOS
            if (IEND-JNPOS<LENOFF(IOLD)) goto 180
            IBEG = LICN + 1
!           Reset pointers to the beginning of the rows.
            do 100 I = 2, N
100         IW1(I,1) = IW1(I-1,1) + LENOFF(I-1)

!           Row IOLD is now split into diagonal and off-diagonal parts.
110         IROWB = IW1(IOLD,1)
            LENI = 0
            IROWE = IROWB + LENOFF(IOLD) - 1
!           Backward scan of whole of row IOLD (in original matrix).
            if (IROWE<IROWB) goto 130
            do 120 JJ = IROWB, IROWE
              J = IROWE - JJ + IROWB
              JOLD = ICN(J)
!             IW(:,2) holds the inverse permutation to IQ.
!             It was set to this in MC13D.
              JNEW = IW(JOLD,2)
!             IF (JNEW < I1) THEN ...
!             Element is in off-diagonal block and so is
!             left in situ.
              if (JNEW<I1) goto 120
!             Element is in diagonal block and is moved to
!             the end of the storage.
              IEND = IEND - 1
              A(IEND) = A(J)
              ICN(IEND) = JNEW
              IBEG = min(IBEG,J)
              ICN(J) = 0
              LENI = LENI + 1
120         end do

            LENOFF(IOLD) = LENOFF(IOLD) - LENI
130         LENR(INEW) = LENI
          end do

          IP(I2) = -IP(I2)
        end do
!       Resets IP(N) to positive value.
        IP(N) = -IP(N)
!       IDISP(2) is position of first element in diagonal blocks.
        IDISP(2) = IEND

!       This compress is used to move all off-diagonal elements
!       to the front of the storage.
        if (IBEG>LICN) goto 200
        JNPOS = IBEG
        ILEND = IDISP(1) - 1
        do 160 J = IBEG, ILEND
          if (ICN(J)==0) goto 160
          ICN(JNPOS) = ICN(J)
          A(JNPOS) = A(J)
          JNPOS = JNPOS + 1
160     end do
!       IDISP(1) is first position after last element of
!       off-diagonal blocks.
        IDISP(1) = JNPOS
        goto 200

!       Error return
170     if (LP/=0) write (LP,90000) NUMNZ
90000   format (' Matrix is structurally singular, rank = ',I6)
        IDISP(1) = -1
        goto 190
180     if (LP/=0) write (LP,90001) N
90001   format (' LICN is not big enough; increase by ',I6)
        IDISP(1) = -2
190     if (LP/=0) write (LP,90002)
90002   format (' + Error return from MC23AD')
200     return

      end subroutine MC23AD
!_______________________________________________________________________

      subroutine MC24AD(N,ICN,A,LICN,LENR,LENRL,W)
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer :: LICN, N
! ..
! .. Array Arguments ..
        real (WP) :: A(LICN), W(N)
        integer :: ICN(LICN), LENR(N), LENRL(N)
! ..
! .. Local Scalars ..
        real (WP) :: AMAXL, AMAXU, WROWL
        integer :: I, J, J0, J1, J2, JJ
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, MAX
! ..
! .. FIRST EXECUTABLE STATEMENT MA24AD
! ..
        AMAXL = ZERO
        W(1:N) = ZERO
        J0 = 1
        do 60 I = 1, N
          if (LENR(I)==0) goto 60
          J2 = J0 + LENR(I) - 1
          if (LENRL(I)==0) goto 30
!         Calculation of 1-norm of L.
          J1 = J0 + LENRL(I) - 1
          WROWL = ZERO
          do 20 JJ = J0, J1
20        WROWL = WROWL + abs(A(JJ))
!         AMAXL is the maximum norm of columns of L so far found.
          AMAXL = max(AMAXL,WROWL)
          J0 = J1 + 1
!         Calculation of norms of columns of U(MAX-NORMS).
30        J0 = J0 + 1
          if (J0>J2) goto 50
          do 40 JJ = J0, J2
            J = ICN(JJ)
40        W(J) = max(abs(A(JJ)),W(J))
50        J0 = J2 + 1
60      end do
!       AMAXU is set to maximum max-norm of columns of U.
        AMAXU = ZERO
        do I = 1, N
          AMAXU = max(AMAXU,W(I))
        end do
!       GROFAC is MAX U max-norm times MAX L 1-norm.
        W(1) = AMAXL*AMAXU
        return

      end subroutine MC24AD
!_______________________________________________________________________

      subroutine MC19AD(N,NA,A,IRN,ICN,R,C,W)
! ..
!     MC19A was altered to use same precision for R, C,
!     and W as is used for other variables in program.
! ..
!     REAL A(NA),R(N),C(N),W(N,5)
!     R(I) is used to return log(scaling factor for row I).
!     C(J) is used to return log(scaling factor for col J).
!     W(I,1), W(I,2) hold row,col non-zero counts.
!     W(J,3) holds - COL J LOG during execution.
!     W(J,4) holds 2-iteration change in W(J,3).
!     W(I,5) is used to save average element log for row I.
!     INTEGER*2 IRN(NA), ICN(NA)
!     IRN(K) gives row number of element in A(K).
!     ICN(K) gives col number of element in A(K).
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer :: N, NA
! ..
! .. Array Arguments ..
        real (WP) :: A(NA), C(N), R(N), W(N,5)
        integer :: ICN(NA), IRN(NA)
! ..
! .. Local Scalars ..
        real (WP) :: E, E1, EM, Q, Q1, QM, S, S1, SM, SMIN, U, V
        integer :: I, I1, I2, ITER, J, K, MAXIT
! ..
! .. Intrinsic Functions ..
        intrinsic ABS, LOG, real
! ..
! .. Data Statements ..
!       MAXIT is the maximal permitted number of iterations.
!       SMIN is used in a convergence test on (residual norm)**2.
        data MAXIT/100/, SMIN/0.1/
! ..
! .. FIRST EXECUTABLE STATEMENT MA19AD
! ..
!       Check scalar data.
        IFAIL = 1
        if (N<1) goto 210
        IFAIL = 2
!       IF (N > 32767)GOTO 230
        IFAIL = 2
        IFAIL = 0

!       Initialise for accumulation of sums and products.
        C(1:N) = ZERO
        R(1:N) = ZERO
        W(1:N,1:4) = ZERO
        if (NA<=0) goto 220
        do 40 K = 1, NA
          U = abs(A(K))
!         IF (U == ZERO) GOTO 30
          if (abs(U-ZERO)<=ZERO) goto 40
          U = log(U)
          I1 = IRN(K)
          I2 = ICN(K)
          if (I1>=1 .and. I1<=N .and. I2>=1 .and. I2<=N) goto 30
          if (LP>0) write (LP,90000) K, I1, I2
90000     format (' MC19 error. Element ',I5,' is in row ',I5,' and column ', &
            I5)
          IFAIL = 3
          goto 40
!         Count row/col non-zeros and compute rhs vectors.
30        W(I1,1) = W(I1,1) + 1.
          W(I2,2) = W(I2,2) + 1.
          R(I1) = R(I1) + U
          W(I2,3) = W(I2,3) + U
40      end do
        if (IFAIL==3) goto 210

!       Divide rhs by diagonal matrices.
        do I = 1, N
!        IF (W(I,1) == ZERO) W(I,1) = 1.
          if (abs(W(I,1))<=ZERO) W(I,1) = ONE
          R(I) = R(I)/W(I,1)
!         SAVE R(I) FOR USE AT END.
          W(I,5) = R(I)
!         IF (W(I,2) == ZERO) W(I,2) = 1.
          if (abs(W(I,2))<=ZERO) W(I,2) = ONE
          W(I,3) = W(I,3)/W(I,2)
        end do
!       SM = SMIN*FLOAT(NA)
        SM = SMIN*real(NA)
!       Sweep to compute initial residual vector.
        do 60 K = 1, NA
!         IF (A(K) == ZERO) GOTO 80
          if (abs(A(K))<=ZERO) goto 60
          I = IRN(K)
          J = ICN(K)
          R(I) = R(I) - W(J,3)/W(I,1)
60      end do

!       Initialise iteration.
        E = ZERO
        Q = 1.
        S = ZERO
        do I = 1, N
          S = S + W(I,1)*R(I)**2
        end do
        if (S<=SM) goto 160

!       Iteration loop.
        do ITER = 1, MAXIT
!         Sweep through matrix to update residual vector.
          do 80 K = 1, NA
!           IF (A(K) == ZERO) GOTO 130
            if (abs(A(K))<=ZERO) goto 80
            I = ICN(K)
            J = IRN(K)
            C(I) = C(I) + R(J)
80        end do
          S1 = S
          S = ZERO
          do I = 1, N
            V = -C(I)/Q
            C(I) = V/W(I,2)
            S = S + V*C(I)
          end do
          E1 = E
          E = Q*S/S1
          Q = 1. - E
          if (S<=SM) E = ZERO
!         Update residual.
          do I = 1, N
            R(I) = R(I)*E*W(I,1)
          end do
          if (S<=SM) goto 180
          EM = E*E1
!         Sweep through matrix to update residual vector.
          do 110 K = 1, NA
!           IF (A(K) == ZERO) GOTO 152
            if (abs(A(K))<=ZERO) goto 110
            I = IRN(K)
            J = ICN(K)
            R(I) = R(I) + C(J)
110       end do
          S1 = S
          S = ZERO
          do I = 1, N
            V = -R(I)/Q
            R(I) = V/W(I,1)
            S = S + V*R(I)
          end do
          E1 = E
          E = Q*S/S1
          Q1 = Q
          Q = 1. - E
!         Special fixup for last iteration.
          if (S<=SM) Q = 1.
!         Update colulm scaling powers.
          QM = Q*Q1
          do I = 1, N
            W(I,4) = (EM*W(I,4)+C(I))/QM
            W(I,3) = W(I,3) + W(I,4)
          end do
          if (S<=SM) goto 160
!         Update residual.
          do I = 1, N
            C(I) = C(I)*E*W(I,2)
          end do
        end do
160     do I = 1, N
          R(I) = R(I)*W(I,1)
        end do

!       Sweep through matrix to prepare to get row scaling powers.
180     do 190 K = 1, NA
!         IF (A(K) == ZERO) GOTO 200
          if (abs(A(K))<=ZERO) goto 190
          I = IRN(K)
          J = ICN(K)
          R(I) = R(I) + W(J,3)
190     end do

!       Final conversion to output values.
        do I = 1, N
          R(I) = R(I)/W(I,1) - W(I,5)
          C(I) = -W(I,3)
        end do
        goto 220
210     if (LP>0) write (LP,90001) IFAIL
90001   format (' Error return ',I2,' from MC19')
220     return
      end subroutine MC19AD
! End of  MA28 subroutines.
!_______________________________________________________________________
! Beginning of JACSP routines.
! Change Record:
! ST 09-01-05
!   Convert to F90 using Metcalf converter
!   Trim trailing blanks using trimem.pl
!   Convert arithmetic operators to F90
!   Replace R1MACH by EPSILON
!   Run through nag tools suite
! ST 09-17-05
!   Delete RDUM and IDUM
!   Change FCN argument list
! ST 09-18-05
!   Add subroutine DGROUPDS to define the IGP and JGP arrays for JACSP
! ST 09-20-05
!   Modify JACSP to produce JACSPDB for sparse, dense, and
!   banded Jacobians
!_______________________________________________________________________

    subroutine JACSP(FCN,N,T,Y,F,FJAC,NRFJAC,YSCALE,FAC,IOPT,WK, &
      LWK,IWK,LIWK,MAXGRP,NGRP,JPNTR,INDROW)

! BEGIN PROLOGUE JACSP
! DATE WRITTEN 850415
! CATEGORY NO. D4
! KEYWORDS NUMERICAL DIFFERENCING, SPARSE JACOBIANS
! AUTHOR SALANE, DOUGLAS E., SANDIA NATIONAL LABORATORIES
!        NUMERICAL MATHEMATICS DIVISION,1642
!        ALBUQUERQUE, NM 87185

! SUBROUTINE JACSP USES FINITE DIFFERENCES TO COMPUTE THE JACOBIAN OF
! A SPARSE SYSTEM OF N EQUATIONS AND N UNKNOWNS. JACSP IS DESIGNED FOR
! USE IN NUMERICAL METHODS FOR SOLVING NONLINEAR PROBLEMS WHERE THE
! JACOBIAN IS EVALUATED REPEATEDLY AND OFTEN AT NEIGHBORING ARGUMENTS
! (E.G., NEWTON'S METHOD OR A BDF METHOD FOR SOVING STIFF ORDINARY
! DIFFERENTIAL EQUATIONS). JACSP IS INTENDED FOR APPLICATIONS IN WHICH
! THE REQUIRED JACOBIANS ARE LARGE AND SPARSE.

! TAKING ADVANTAGE OF SPARSITY.

! SUBROUTINE JACSP TAKES ADVANTAGE OF THE SPARSITY OF A MATRIX TO REDUCE
! THE NUMBER OF FUNCTION EVALUATIONS REQUIRED TO COMPUTE THE JACOBIAN.
! TO REALIZE THIS ADVANTAGE, THE USER MUST PROVIDE JACSP WITH A COLUMN
! GROUPING. THIS MEANS THE USER MUST ASSIGN COLUMNS OF THE JACOBIAN TO
! MUTUALLY EXCLUSIVE GROUPS BASED ON THE SPARSITY OF THE COLUMNS.THE
! DIFFERENCES REQUIRED TO COMPUTE THE NONZERO ELEMENTS IN ALL COLUMNS IN
! A GROUP ARE FORMED USING ONLY ONE ADDITIONAL FUNCTION EVALUATION. FOR
! MORE DETAILS, THE USER IS REFERRED TO THE REPORT BY D.E. SALANE AND
! L.F. SHAMPINE (REF. 1).

! THE SUBROUTINE DVDSM (REF 2.) IS THE WAY MOST USERS WILL DETERMINE A
! COLUMN GROUPING FOR JACSP. THE USE OF DSM AND JACSP IS DESCRIBED
! LATER IN THE PROLOGUE.

! STORAGE.

! JACSP REQUIRES THE USER TO PROVIDE A SPARSE DATA STRUCTURE FOR THE
! JACOBIAN TO BE COMPUTED. JACSP REQUIRES THE USER TO SPECIFY THE
! INDICES OF THE NONZERO ELEMENTS IN A COLUMN PACKED SPARSE DATA
! STRUCTURE (SEE REF.1). THE SUBROUTINE DVDSM( REF 2.) CAN BE USED TO
! ESTABLISH THE SPARSE DATA STRUCTURE REQUIRED BY JACSP. THE USE OF
! JACSP AND DSM IS DESCRIBED LATER IN THE PROLOGUE.

! ON OUTPUT, SUBROUTINE JACSP WILL RETURN THE JACOBIAN IN ANY ONE OF
! THE FOLLOWING THREE FORMATS SPECIFIED BY THE USER.

! (1) FULL STORAGE FORMAT.......THE COMPUTED JACOBIAN IS STORED IN
!      AN ARRAY WHOSE ROW AND COLUMN DIMENSIONS ARE THE NUMBER OF
!      EQUATIONS.

! (2) BANDED STORAGE FORMAT.....THE COMPUTED MATRIX IS STORED IN A TWO
!      DIMENSIONAL ARRAY WHOSE ROW DIMENSION IS
!      2*(NUMBER OF LOWER DIAGONALS) + (THE NUMBER OF UPPER DIAGONALS)
!      + 1. THE COLUMN DIMENSION OF THIS ARRAY IS THE NUMBER OF
!      EQUATIONS. THIS STORAGE FORMAT IS COMPATIBLE WITH THE LINPACK
!      GENERAL BAND MATRIX EQUATION SOLVER.

! (3) SPARSE STORAGE FORMAT.....THE COMPUTED JACOBIAN IS STORED IN A
!      ONE DIMENSIONAL ARRAY WHOSE LENGTH IS THE NUMBER OF NONZERO
!      ELEMENTS IN THE JACOBIAN.

! DESCRIPTION

!  SUBROUTINE PARAMETERS

!  FCN()......A USER-PROVIDED FUNCTION (SEE SUBROUTINE DESCRIPTION).
!  N..........THE NUMBER OF EQUATIONS.
!  T..........A SCALAR VARIABLE. T IS PROVIDED SO USERS CAN PASS THE
!             VALUE OF AN INDEPENDENT VARIALBE TO THE FUNCTION
!             EVALUATION ROUTINE.
!  Y(*).......AN ARRAY OF DIMENSION N. THE POINT AT WHICH THE
!             JACOBIAN IS TO BE EVALUATED.
!  F(*).......AN ARRAY OF DIMENSION N. THE EQUATIONS EVALUATED AT
!             THE POINT Y.
!  FJAC(*,*)..AN ARRAY OF WHOSE DIMENSIONS DEPEND ON THE STORAGE
!      FORMAT SELECTED BY THE USER. THE ROW AND COLUMN DIMENSIONS ARE
!      SET BY THE USER. IF THE SPARSE MATRIX FORMAT IS SELECTED, FJAC
!      SHOULD BE TREATED AS A ONE DIMENSIONAL ARRAY BY THE CALLING
!      PROGRAM. FOR FURTHER DETAILS, SEE THE DESCRIPTION OF THE
!      PARAMETER IOPT(1) WHICH CONTROLS THE STORAGE FORMAT.

!      NOTE THAT IF THE BANDED OR FULL OPTION IS USED,THE USER SHOULD
!      ZERO OUT THOSE POSITIONS OF FJAC THAT WILL NOT BE ACCESSED BY
!      JACSP. JACSP DOES NOT ZERO OUT THE MATRIX FJAC BEFORE
!      COMPUTING THE JACOBIAN. POSITIONS OF FJAC THAT ARE NOT ASSIGNED
!      VALUES BY JACSP WILL BE THE SAME ON EXIT AS ON ENTRY TO JACSP.

!  NRFJAC.....THE NUMBER OF ROWS IN FJAC AND THE LEADING DIMENSION
!             OF FJAC (SEE IOPT(1)).
!  NCFJAC.....THE NUMBER OF COLUMNS IN FJAC AND THE COLUMN DIMENSION
!             OF FJAC (SEE IOPT(1)).
!  YSCALE(*)..AN ARRAY OF DIMENSION N. YSCALE(I) CONTAINS A
!     REPRESENTATIVE MAGNITUDE FOR Y(I). YSCALE(I) MUST BE POSITIVE.
!     YSCALE IS AN OPTIONAL FEATURE OF THE ROUTINE. IF THE USER DOES
!     NOT WISH TO PROVIDE YSCALE, IT CAN BE TREATED AS A DUMMY SCALAR
!     VARIABLE (SEE IOPT(3)).

!  FAC(*).....AN ARRAY OF DIMENSION N. FAC CONTAINS A PERCENTAGE FOR
!     USE IN COMPUTING THE INCREMENT.

!     THE NORMAL WAY TO USE THE CODE IS TO LET JACSP ADJUST FAC VALUES.
!     IN THIS CASE FAC IS INITIALIZED BY JACSP. ALSO, THE USER MUST
!     NOT ALTER FAC VALUES BETWEEN CALLS TO JACSP.

!     THE USER MAY PROVIDE FAC VALUES IF DESIRED (SEE IOPT(4) FOR
!     DETAILS). IF THE USER PROVIDES FAC VALUES, FAC(I) SHOULD BE
!     SHOULD BE SET TO A VALUE BETWEEN O AND 1. JACSP WILL NOT PERMIT
!     FAC(I) TO BE SET TO A VALUE THAT RESULTS IN TOO SMALL AN
!     INCREMENT. JACSP ENSURES THAT
!                     FACMIN <= FAC(I) <= FACMAX.
!     FOR FURTHER DETAILS ON FACMIN AND FACMAX SEE
!     THE REPORT(REF.3).

!  IOPT(*)....AN INTEGER ARRAY OF LENGTH 5 FOR USER SELECTED OPTIONS.

!          IOPT(1) CONTROLS THE STORAGE FORMAT.

!             IOPT(1) = 0 INDICATES FULL STORAGE FORMAT. SET BOTH
!             NRFJAC = N AMD NCFJAC = N.

!             IOPT(1) = 1 INDICATES BANDED STORAGE FORMAT. SET
!             NRFJAC = 2 * ML + MU + 1 WHERE
!             ML = NUMBER OF SUB-DIAGONALS BELOW THE MAIN DIAGONAL AND
!             MU = NUMBER OF SUPER-DIAGONALS ABOVE THE MAIN DIAGONAL.
!             SET NCFJAC = N.

!             IOPT(1) = 2 INDICATES SPARSE STORAGE FORMAT. SET
!             NRFJAC = TO THE NUMBER OF NONZEROS IN THE JACOBIAN AND
!             SET NCFJAC = 1.

!          IOPT(2) MUST BE SET TO THE BANDWIDTH OF THE MATRIX.
!          IOPT(2) NEED ONLY BE PROVIDED IF BANDED STORAGE
!          FORMAT IS REQUESTED(I.E., IOPT(1) = 1).

!          IOPT(3) ALLOWS THE USER TO PROVIDE TYPICAL VALUES
!          TO BE USED IN COMPUTING INCREMENTS FOR DIFFERENCING.

!             IOPT(3) = 0 INDICATES Y VALUES ARE USED IN COMPUTING
!             INCREMENTS. IF IOPT(3) = 0, NO STORAGE IS REQUIRED FOR
!             YSCALE AND IT MAY BE TREATED AS A DUMMY VARIABLE.

!             IOPT(3) = 1 INDICATES THAT YSCALE VALUES ARE TO BE USED.
!             IF IOPT(3) = 1, THE USER MUST PROVIDE AN ARRAY OF NONZERO
!             VALUES IN YSCALE.

!          IOPT(4) ALLOWS THE USER TO PROVIDE THE VALUES USED
!          IN THE FAC ARRAY TO COMPUTE INCREMENTS FOR DIFFERENCING.

!             IF IOPT(4) = 0, EACH COMPONENT OF FAC IS
!             SET TO THE SQUARE ROOT OF MACHINE UNIT ROUNDOFF ON THE
!             FIRST CALL TO JACSP. IOPT(4) IS SET TO ONE ON RETURN.

!             IF IOPT(4) = 1, EACH COMPONENT OF FAC MUST BE SET
!             BY THE CALLING ROUTINE. UNLESS THE USER WISHES TO
!             INITIALIZE FAC, THE FAC ARRAY SHOULD NOT BE ALTERED
!             BETWEEN SUBSEQUENT CALLS TO JACSP. ALSO, THE USER
!             SHOULD NOT CHANGE THE VALUE OF IOPT(4) RETURNED BY
!             JACSP.

!          IOPT(5) IS NOT USED IN JACSP.

!  WK(*)......A WORK ARRAY OF DIMENSION AT LEAST 3*N
!  LWK........THE LENGTH OF THE WORK ARRAY. LWK IS AT LEAST 3*N.
!  IWK(*)...AN INTEGER ARRAY OF LENGTH LIWK = 50 + N WHICH GIVES
!      DIAGNOSTIC INFORMATION IN POSITIONS 1 THROUGH 50. POSITIONS 51
!      THROUGH 50 + N ARE USED AS INTEGER WORKSPACE.

!      IWK(1) GIVES THE NUMBER OF TIMES THE INCREMENT FOR DIFFERENCING
!      (DEL) WAS COMPUTED AND HAD TO BE INCREASED BECAUSE (Y(JCOL)+DEL)
!      -Y(JCOL)) WAS TOO SMALL RELATIVE TO Y(JCOL) OR YSCALE(JCOL).

!      IWK(2) GIVES THE NUMBER OF COLUMNS IN WHICH THREE ATTEMPTS WERE
!      MADE TO INCREASE A PERCENTAGE FACTOR FOR DIFFERENCING (I.E., A
!      COMPONENT IN THE FAC ARRAY) BUT THE COMPUTED DEL REMAINED
!      UNACCEPTABLY SMALL RELATIVE TO Y(JCOL) OR YSCALE(JCOL). IN SUCH
!      CASES THE PERCENTAGE FACTOR IS SET TO THE SQUARE ROOT OF THE UNIT
!      ROUNDOFF OF THE MACHINE. THE FIRST 10 COLUMNS ARE GIVEN IN
!      IWK(21),...,IWK(30).

!      IWK(3) GIVES THE NUMBER OF COLUMNS IN WHICH THE COMPUTED DEL WAS
!      ZERO TO MACHINE PRECISION BECAUSE Y(JCOL) OR YSCALE(JCOL) WAS
!      ZERO. IN SUCH CASES DEL IS SET TO THE SQUARE ROOT OF THE UNIT
!      ROUNDOFF.

!      IWK(4) GIVES THE NUMBER OF COLUMNS WHICH HAD TO BE RECOMPUTED
!      BECAUSE THE LARGEST DIFFERENCE FORMED IN THE COLUMN WAS VERY
!      CLOSE TO ZERO RELATIVE TO SCALE WHERE

!                    SCALE = MAX(F(Y),F(Y+DEL))
!                                 I     I

!      AND I DENOTES THE ROW INDEX OF THE LARGEST DIFFERENCE IN THE
!      COLUMN CURRENTLY BEING PROCESSED. IWK(31),...,IWK(40) GIVES THE
!      FIRST 10 OF THESE COLUMNS.

!      IWK(5) GIVES THE NUMBER OF COLUMNS WHOSE LARGEST DIFFERENCE IS
!      CLOSE TO ZERO RELATIVE TO SCALE AFTER THE COLUMN HAS BEEN
!      RECOMPUTED. THE FIRST 10 OF THESE COLUMNS ARE GIVEN IN POSITIONS
!      IWK(41)...IWK(50).
!      IWK(6) GIVES THE NUMBER OF TIMES SCALE INFORMATION WAS NOT
!      AVAILABLE FOR USE IN THE ROUNDOFF AND TRUNCATION ERROR TESTS.
!      THIS OCCURS WHEN
!                    MIN(F(Y),F(Y+DEL)) = 0.
!                          I     I
!      WHERE I IS THE INDEX OF THE LARGEST DIFFERENCE FOR THE COLUMN
!      CURRENTLY BEING PROCESSED.
!      IWK(7) GIVES THE NUMBER OF TIMES THE FUNCTION EVALUATION ROUTINE
!      WAS CALLED.
!      IWK(8) GIVES THE NUMBER OF TIMES A COMPONENT OF THE FAC ARRAY WAS
!      REDUCED BECAUSE CHANGES IN FUNCTION VALUES WERE LARGE AND EXCESS
!      TRUNCATION ERROR WAS SUSPECTED. IWK(11),...,IWK(20) GIVES THE FIRST
!      10 OF THESE COLUMNS.
!      IWK(9) AND IWK(10) ARE NOT USED IN JACSP.
!  LIWK.....THE LENGTH OF THE ARRAY IWK. LIWK = 50 + N.
!      THE FOLLOWING PARAMETERS MAY BE PROVIDED BY THE USER OR
!      INITIALIZED BY THE SUBROUTINE DVDSM (SEE LONG DESCRIPTION
!      SECTION OF THE PROLOGUE).
!  PARAMETERS CONTAINING COLUMN GROUPING:
!  MAXGRP.....THE NUMBER OF DIFFERENT GROUPS TO WHICH THE COLUMNS
!             HAVE BEEN ASSIGNED.
!  NGRP(*)....AN ARRAY OF LENGTH N THAT CONTAINS THE COLUMN GROUPING.
!             NGRP(I) IS THE GROUP TO WHICH THE I-TH COLUMN HAS BEEN
!             ASSIGNED.
!             (SEE USE OF DSM AND JACSP FOR DETERMINING MAXGRP AND NGRP)
!  PARAMETERS CONTAINING THE SPARSE DATA STRUCTURE:
!  JPNTR(*)...AN ARRAY OF LENGTH N+1 THAT CONTAINS POINTERS TO
!             THE ROW INDICES IN INDROW (SEE INDROW).
!  INDROW(*)..AN ARRAY WHOSE LENGTH IS THE NUMBER OF NONZERO ELEMENTS
!             IN THE JACOBIAN. INDROW CONTAINS THE ROW INDICES
!             STORED IN COLUMN PACKED FORMAT. IN OTHER WORDS, THE
!             ROW INDICES OF THE NONZERO ELEMENTS OF A GIVEN
!             COLUMN OF THE JACOBIAN ARE STORED CONTIGUOUSLY IN
!             INDROW. JPNTR(JCOL) POINTS TO THE ROW INDEX OF THE FIRST
!             NONZERO ELEMENT IN THE COLUMN JCOL. JPNTR(JCOL+1) - 1
!             POINTS TO THE ROW INDEX OF THE LAST NONZERO ELEMENT IN
!             THE COLUMN JCOL.
!             (SEE USE OF DSM AND JACSP FOR DETERMINING INDROW AND JPNTR)
!  REQUIRED SUBROUTINES: SUBROUTINE FCN(N,T,Y,F)
!     T......AN INDEPENDENT SCALAR VARIABLE WHICH MAY BE USED
!           IN EVALUATING F(E.G., F(Y(T),T)).
!     Y(*)...AN ARRAY OF DIMENSION N WHICH CONTAINS THE POINT AT
!            WHICH THE EQUATIONS ARE TO BE EVALUATED.
!     F(*)...AN ARRAY OF DIMENSION N WHICH ON RETURN FROM FCN
!            CONTAINS THE EQUATIONS EVALUATED AT Y.

! LONG DESCRIPTION

! ROUNDOFF AND TRUNCATION ERRORS.
! SUBROUTINE JACSP TAKES ADVANTAGE OF THE WAY IN WHICH THE JACOBIAN
! IS EVALUATED TO ADJUST INCREMENTS FOR DIFFERENCING TO CONTROL
! ROUNDOFF AND TRUNCATION ERRORS. THE ROUTINE SELDOM REQUIRES MORE
! THAN ONE ADDITIONAL FUNCTION EVALUATION TO COMPUTE A COLUMN OF THE
! JACOBIAN. ALSO, THE ROUTINE RETURNS A VARIETY OF ERROR DIAGNOSTICS
! TO WARN USERS WHEN COMPUTED DERIVATIVES MAY NOT BE ACCURATE.
! WARNING: JACSP CAN NOT GUARANTEE THE ACCURACY OF THE COMPUTED
! DERIVATIVES. IN ORDER O SAVE ON FUNCTION EVALUATIONS, HEURISTIC
! TECNIQUES FOR INCREMENT ADJUSTMENT AND SAFEGUARDING INCREMENTS ARE
! USED. THESE USUALLY WORK WELL.
! WARNING: SOME OF THE DIAGNOSTICS RETURNED CAN ONLY BE INTERPRETED
! WITH A DETAILED KNOWLEDGE OF THE ROUTINE. NEVERTHELESS, THEY ARE
! PROVIDED TO GIVE USERS FULL ACCESS TO THE INFORMATION PRODUCED BY
! THE SUBROUTINE.
! USE OF DSM AND JACSP.
! SUBROUTINE DVDSM CAN BE USED TO DETERMINE THE COLUMN GROUPING (MAXGRP
! AND NGRP(*)) AND THE SPARSE DATA STRUCTURE VARIABLES (JPNTR(*) AND
! INDROW(*)). THE USER CAN CALL DVDSM ONCE TO INITIALIZE
! MAXGRP,NGRP,JPNTR AND INDROW. JACSP MAY THEN BE CALLED REPEATEDLY TO
! EVALUATE THE JACOBIAN. THE FOLLOWING ARE THE IMPORTANT VARIABLES IN
! THE DSM CALLING SEQUENCE.
! SUBROUTINE DVDSM(...,INDROW,INDCOL,NGRP,MAXGRP,..,JPNTR,...)
! ON INPUT, THE USER MUST PROVIDE DSM WITH THE INTEGER ARRAYS INDROW AND
! INDCOL. THE PAIR
!                  (INDROW(I),INDCOL(I))
! PROVIDES THE INDEX OF A NONZERO ELEMENT OF THE JACOBIAN. THE LENGTH OF
! INDROW AND INDCOL IS THE NUMBER OF NONZERO ELEMENTS IN THE JACOBIAN.
! NO ORDERING OF THE INDICES FOR NONZERO ELEMENTS IS REQUIRED IN THE
! ARRAYS INDROW AND INDCOL.
! ON RETURN FROM DSM, MAXGRP,NGRP,INDROW,JPNTR ARE INITIALIZED
! FOR INPUT TO JACSP. THE USER MUST NOT CHANGE MAXGRP,NGRP,JPNTR OR
! INDROW AFTER CALLING DSM.
! WE REFER THE READER TO THE IN CODE DOCUMENTATION FOR DSM OR THE REPORT
! BY MORE AND COLEMAN (REF 2.) FOR FURTHER DETAILS.
! REFERENCES
!  (1) D.E. SALANE AND L. F. SHAMPINE
!      "AN ECONOMICAL AND EFFICIENT ROUTINE FOR COMPUTING
!      SPARSE JACOBIANS", REPORT NO. SAND85-0977, SANDIA NATIONAL
!      LABORATORIES, ALBUQUERQUE,NM,87185.
!  (2) T.F. COLEMAN AND J.J. MORE
!      "SOFTWARE FOR ESTIMATING SPARSE JACOBIAN MATRICES" ACM TOMS,
!      V10.,N0.3,SEPT. 1984.
!  (3) D.E. SALANE AND L. F. SHAMPINE
!      "THREE ADAPTIVE ROUTINES FOR FORMING JACOBIANS NUMERICALLY,"
!       REPORT NO. SAND86- ****, SANDIA NATIONAL LABORATORIES,
!       ALBUQUERQUE, NM, 87185.
!  REQUIRED FORTRAN INTRINSIC FUNCTIONS: MAX,MIN,ABS,SIGN
!  MACHINE DEPENDENT CONSTANT: U (MACHINE UNIT ROUNDOFF)
!  JACSP SETS THE REQUIRED MACHINE CONSTANT USING THE F90
!  INTRINSIC EPSILON.
!  END PROLOGUE JACSP

     implicit none

! .. Parameters ..
      integer, parameter :: WP = kind(0.0D0)
! ..
! .. Scalar Arguments ..
      real (WP) :: T
      integer :: LIWK, LWK, MAXGRP, N, NRFJAC
! ..
! .. Array Arguments ..
      real (WP) :: F(N), FAC(N), FJAC(NRFJAC,*), WK(LWK), Y(N), YSCALE(*)
      integer :: INDROW(*), IOPT(5), IWK(LIWK), JPNTR(*), NGRP(N)
! ..
! .. Subroutine Arguments ..
      external FCN
! ..
! .. Local Scalars ..
      real (WP) :: ADIFF, AY, DEL, DELM, DFMJ, DIFF, DMAX, EXPFMN, FACMAX, &
        FACMIN, FJACL, FMJ, ONE, P125, P25, P75, P875, PERT, RDEL, RMNFDF, &
        RMXFDF, SDF, SF, SGN, T1, T2, U, U3QRT, U7EGT, UEGT, UMEGT, UQRT, &
        USQT, ZERO
      integer :: IDXL, IDXU, IFLAG1, IFLAG2, IRCMP, IRDEL, IROW, IROWB, &
        IROWMX, ITRY, J, JCOL, KT1, KT2, KT3, KT4, KT5, L, NID1, NID2, NID3, &
        NID4, NID5, NID6, NIFAC, NT2, NUMGRP
! ..
! .. Intrinsic Functions ..
!     INTRINSIC ABS, EPSILON, KIND, MAX, MIN, SIGN, SQRT
      intrinsic ABS, EPSILON, MAX, MIN, SIGN, SQRT
! ..
! .. Data Statements ..
      data PERT/2.0E+0_WP/, FACMAX/1.E-1_WP/, EXPFMN/.75E+0_WP/
      data ONE/1.0E0_WP/, ZERO/0.0E0_WP/
      data P125/.125E+0_WP/, P25/.25E+0_WP/, P75/.75E+0_WP/, P875/.875E+0_WP/
      data NIFAC/3/, NID1/10/, NID2/20/, NID3/30/, NID4/40/, NID5/50/
! ..
!     COMPUTE ALGORITHM AND MACHINE CONSTANTS.
! ..
! ..  FIRST EXECUTABLE STATEMENT JACSP
! ..
      U = epsilon(ONE)
      USQT = sqrt(U)
      UEGT = U**P125
      UMEGT = ONE/UEGT
      UQRT = U**P25
      U3QRT = U**P75
      U7EGT = U**P875
      FACMIN = U**EXPFMN

      if (IOPT(4) == 0) then
        IOPT(4) = 1
        do 10 J = 1, N
          FAC(J) = USQT
10      continue
      end if
      do 20 J = 1, 50
        IWK(J) = 0
20    continue
      KT1 = NID1
      KT2 = NID2
      KT3 = NID3
      KT4 = NID4
      KT5 = NID5
      NID6 = LIWK
      NT2 = 2*N
      do NUMGRP = 1, MAXGRP
!       COMPUTE AND SAVE THE INCREMENTS FOR THE COLUMNS IN GROUP NUMGRP.
        IRCMP = 0
        ITRY = 0
        do 30 J = NID5 + 1, NID6
          IWK(J) = 0
30      continue
40      continue
        do JCOL = 1, N
          if (NGRP(JCOL) == NUMGRP) then
            WK(N+JCOL) = Y(JCOL)
!           COMPUTE DEL. IF DEL IS TOO SMALL INCREASE FAC(JCOL) AND RECOMPUTE
!           DEL. NIFAC ATTEMPTS ARE MADE TO INCREASE FAC(JCOL) AND FIND AN
!           APPROPRIATE DEL. IF DEL CANT BE FOUND IN THIS MANNER, DEL IS COMPUTED
!           WITH FAC(JCOL) SET TO THE SQUARE ROOT OF THE MACHINE PRECISION (USQT).
!           IF DEL IS ZERO TO MACHINE PRECISION BECAUSE Y(JCOL) IS ZERO OR
!           YSCALE(JCOL) IS ZERO, DEL IS SET TO USQT.
            SGN = sign(ONE,F(JCOL))
            IRDEL = 0
            if (IOPT(3) == 1) then
              AY = abs(YSCALE(JCOL))
            else
              AY = abs(Y(JCOL))
            end if
            DELM = U7EGT*AY

            do 50 J = 1, NIFAC
              DEL = FAC(JCOL)*AY*SGN
!             IF (DEL == ZERO) THEN
              if (abs(DEL) <= ZERO) then

                DEL = USQT*SGN
                if (ITRY == 0) IWK(3) = IWK(3) + 1
              end if
              T1 = Y(JCOL) + DEL
              DEL = T1 - Y(JCOL)
              if (abs(DEL) < DELM) then
                if (J >= NIFAC) goto 50
                if (IRDEL == 0) then
                  IRDEL = 1
                  IWK(1) = IWK(1) + 1
                end if
                T1 = FAC(JCOL)*UMEGT
                FAC(JCOL) = min(T1,FACMAX)
              else
                goto 60
              end if
50          end do

            FAC(JCOL) = USQT
            DEL = USQT*AY*SGN
            IWK(2) = IWK(2) + 1
            if (KT2 < NID3) then
              KT2 = KT2 + 1
              IWK(KT2) = JCOL
            end if

60          continue
            WK(NT2+JCOL) = DEL
            Y(JCOL) = Y(JCOL) + DEL
          end if
        end do

        IWK(7) = IWK(7) + 1
        call FCN(N,T,Y,WK)
        do JCOL = 1, N
          if (NGRP(JCOL) == NUMGRP) Y(JCOL) = WK(N+JCOL)
        end do

!       COMPUTE THE JACOBIAN ENTRIES FOR ALL COLUMNS IN NUMGRP.
!       STORE ENTRIES ACCORDING TO SELECTED STORAGE FORMAT.
!       USE LARGEST ELEMENTS IN A COLUMN TO DETERMINE SCALING
!       INFORMATION FOR ROUNDOFF AND TRUNCATION ERROR TESTS.
        do JCOL = 1, N
          if (NGRP(JCOL) == NUMGRP) then
            IDXL = JPNTR(JCOL)
            IDXU = JPNTR(JCOL+1) - 1
            DMAX = ZERO
            RDEL = ONE/WK(NT2+JCOL)
            IROWMX = 1
            do L = IDXL, IDXU
              IROW = INDROW(L)
              DIFF = WK(IROW) - F(IROW)
              ADIFF = abs(DIFF)
              if (ADIFF >= DMAX) then
                IROWMX = IROW
                DMAX = ADIFF
                SF = F(IROW)
                SDF = WK(IROW)
              end if
              FJACL = DIFF*RDEL
              if (ITRY == 1) WK(IROW) = FJACL

              if (IOPT(1) == 0) then
                if (ITRY == 1) WK(IROW+N) = FJAC(IROW,JCOL)
                FJAC(IROW,JCOL) = FJACL
              end if
              if (IOPT(1) == 1) then
                IROWB = IROW - JCOL + IOPT(2)
                if (ITRY == 1) WK(IROW+N) = FJAC(IROWB,JCOL)
                FJAC(IROWB,JCOL) = FJACL
              end if
              if (IOPT(1) == 2) then
                if (ITRY == 1) WK(IROW+N) = FJAC(L,1)
                FJAC(L,1) = FJACL
              end if
            end do

!           IF A COLUMN IS BEING RECOMPUTED (ITRY=1),THIS SECTION OF THE
!           CODE PERFORMS AN EXTRAPOLATION TEST TO ENABLE THE CODE TO
!           COMPUTE SMALL DERIVATIVES MORE ACCURATELY. THIS TEST IS ONLY
!           PERFORMED ON THOSE COLUMNS WHOSE LARGEST DIFFERENCE IS CLOSE
!           TO ZERO RELATIVE TO SCALE.
            if (ITRY == 1) then
              IFLAG1 = 0
              IFLAG2 = 0
              do 100 J = NID5 + 1, NID6
                if (IWK(J) == JCOL) IFLAG1 = 1
100           continue
              if (IFLAG1 == 1) then
                IFLAG1 = 0
                T1 = WK(IROWMX+N)
                T2 = WK(IROWMX)*FAC(JCOL)
                if (abs(T2) < abs(T1)*PERT) IFLAG2 = 1
              end if

              if (IFLAG2 == 1) then
                IFLAG2 = 0
                T1 = FAC(JCOL)*FAC(JCOL)
                FAC(JCOL) = max(T1,FACMIN)
                do L = IDXL, IDXU
                  IROW = INDROW(L)
                  FJACL = WK(IROW+N)
                  if (IOPT(1) == 0) FJAC(IROW,JCOL) = FJACL
                  if (IOPT(1) == 1) then
                    IROWB = IROW - JCOL + IOPT(2)
                    FJAC(IROWB,JCOL) = FJACL
                  end if
                  if (IOPT(1) == 2) FJAC(L,1) = FJACL
                end do
              end if
            end if

            FMJ = abs(SF)
            DFMJ = abs(SDF)
            RMXFDF = max(FMJ,DFMJ)
            RMNFDF = min(FMJ,DFMJ)

!           IF SCALE INFORMATION IS NOT AVAILABLE, PERFORM NO ROUNDOFF
!           OR TRUNCATION ERROR TESTS. IF THE EXTRAPOLATION TEST HAS
!           CAUSED FAC(JCOL) TO BE RESET TO ITS PREVIOUS VALUE (IAJAC=1)
!           THEN NO FURTHER ROUNDOFF OR TRUNCATION ERROR TESTS ARE
!           PERFORMED.
!           IF (RMNFDF/=ZERO) THEN
            if (abs(RMNFDF) > ZERO) then
!             TEST FOR POSSIBLE ROUNDOFF ERROR (FIRST TEST)
!             AND ALSO FOR POSSIBLE SERIOUS ROUNDOFF ERROR (SECOND TEST).
              if (DMAX <= (U3QRT*RMXFDF)) then
                if (DMAX <= (U7EGT*RMXFDF)) then
                  if (ITRY == 0) then
                    T1 = sqrt(FAC(JCOL))
                    FAC(JCOL) = min(T1,FACMAX)
                    IRCMP = 1
                    if (KT5 < NID6) then
                      KT5 = KT5 + 1
                      IWK(KT5) = JCOL
                    end if
                    IWK(4) = IWK(4) + 1
                    if (KT3 < NID4) then
                      KT3 = KT3 + 1
                      IWK(KT3) = JCOL
                    end if
                  else
                    IWK(5) = IWK(5) + 1
                    if (KT4 < NID5) then
                      KT4 = KT4 + 1
                      IWK(KT4) = JCOL
                    end if
                  end if
                else
                  T1 = UMEGT*FAC(JCOL)
                  FAC(JCOL) = min(T1,FACMAX)
                end if
              end if
!             TEST FOR POSSIBLE TRUNCATION ERROR.
              if (DMAX > UQRT*RMXFDF) then
                T1 = FAC(JCOL)*UEGT
                FAC(JCOL) = max(T1,FACMIN)
                IWK(8) = IWK(8) + 1
                if (KT1 < NID2) then
                  KT1 = KT1 + 1
                  IWK(KT1) = JCOL
                end if
              end if
            else
              IWK(6) = IWK(6) + 1
            end if
          end if
        end do

!       IF SERIOUS ROUNDOFF ERROR IS SUSPECTED, RECOMPUTE ALL
!       COLUMNS IN GROUP NUMGRP.
        if (IRCMP == 1) then
          IRCMP = 0
          ITRY = 1
          goto 40
        end if
        ITRY = 0
      end do
      return

    end subroutine JACSP
!_______________________________________________________________________

    subroutine DEGR(N,INDROW,JPNTR,INDCOL,IPNTR,NDEG,IWA)

!     GIVEN THE SPARSITY PATTERN OF AN M BY N MATRIX A,
!     THIS SUBROUTINE DETERMINES THE DEGREE SEQUENCE FOR
!     THE INTERSECTION GRAPH OF THE COLUMNS OF A.
!     IN GRAPH-THEORY TERMINOLOGY, THE INTERSECTION GRAPH OF
!     THE COLUMNS OF A IS THE LOOPLESS GRAPH G WITH VERTICES
!     A(J), J = 1,2,...,N WHERE A(J) IS THE J-TH COLUMN OF A
!     AND WITH EDGE (A(I),A(J)) IF AND ONLY IF COLUMNS I AND J
!     HAVE A NON-ZERO IN THE SAME ROW POSITION.
!     NOTE THAT THE VALUE OF M IS NOT NEEDED BY DEGR AND IS
!     THEREFORE NOT PRESENT IN THE SUBROUTINE STATEMENT.
!     THE SUBROUTINE STATEMENT IS
!       SUBROUTINE DEGR(N,INDROW,JPNTR,INDCOL,IPNTR,NDEG,IWA)
!     WHERE
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF COLUMNS OF A.
!       INDROW IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW
!         INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!       JPNTR IS AN INTEGER INPUT ARRAY OF LENGTH N + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW.
!         THE ROW INDICES FOR COLUMN J ARE
!               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1.
!         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!       INDCOL IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE
!         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!       IPNTR IS AN INTEGER INPUT ARRAY OF LENGTH M + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL.
!         THE COLUMN INDICES FOR ROW I ARE
!               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1.
!         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!       NDEG IS AN INTEGER OUTPUT ARRAY OF LENGTH N WHICH
!         SPECIFIES THE DEGREE SEQUENCE. THE DEGREE OF THE
!         J-TH COLUMN OF A IS NDEG(J).
!       IWA IS AN INTEGER WORK ARRAY OF LENGTH N.
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983.
!     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE'

     implicit none

! .. Scalar Arguments ..
      integer :: N
! ..
! .. Array Arguments ..
      integer :: INDCOL(*), INDROW(*), IPNTR(*), IWA(N), JPNTR(N+1), NDEG(N)
! ..
! .. Local Scalars ..
      integer :: IC, IP, IR, JCOL, JP
! ..
! ..  FIRST EXECUTABLE STATEMENT DEGR
! ..
!     INITIALIZATION BLOCK.
      NDEG(1:N) = 0
      IWA(1:N) = 0

!     COMPUTE THE DEGREE SEQUENCE BY DETERMINING THE CONTRIBUTIONS
!     TO THE DEGREES FROM THE CURRENT (JCOL) COLUMN AND FURTHER
!     COLUMNS WHICH HAVE NOT YET BEEN CONSIDERED.
      do JCOL = 2, N
        IWA(JCOL) = N
!       DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND
!       TO NON-ZEROES IN THE MATRIX.
        do JP = JPNTR(JCOL), JPNTR(JCOL+1) - 1
          IR = INDROW(JP)
!         FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC)
!         WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX.
          do IP = IPNTR(IR), IPNTR(IR+1) - 1
            IC = INDCOL(IP)
!           ARRAY IWA MARKS COLUMNS WHICH HAVE CONTRIBUTED TO
!           THE DEGREE COUNT OF COLUMN JCOL. UPDATE THE DEGREE
!           COUNTS OF THESE COLUMNS AS WELL AS COLUMN JCOL.
            if (IWA(IC) < JCOL) then
              IWA(IC) = JCOL
              NDEG(IC) = NDEG(IC) + 1
              NDEG(JCOL) = NDEG(JCOL) + 1
            end if
          end do
        end do
      end do
      return

    end subroutine DEGR
!_______________________________________________________________________

    subroutine IDO(M,N,INDROW,JPNTR,INDCOL,IPNTR,NDEG,LIST,MAXCLQ, &
      IWA1,IWA2,IWA3,IWA4)

!     SUBROUTINE IDO
!     GIVEN THE SPARSITY PATTERN OF AN M BY N MATRIX A, THIS
!     SUBROUTINE DETERMINES AN INCIDENCE-DEGREE ORDERING OF THE
!     COLUMNS OF A.
!     THE INCIDENCE-DEGREE ORDERING IS DEFINED FOR THE LOOPLESS
!     GRAPH G WITH VERTICES A(J), J = 1,2,...,N WHERE A(J) IS THE
!     J-TH COLUMN OF A AND WITH EDGE (A(I),A(J)) IF AND ONLY IF
!     COLUMNS I AND J HAVE A NON-ZERO IN THE SAME ROW POSITION.
!     THE INCIDENCE-DEGREE ORDERING IS DETERMINED RECURSIVELY BY
!     LETTING LIST(K), K = 1,...,N BE A COLUMN WITH MAXIMAL
!     INCIDENCE TO THE SUBGRAPH SPANNED BY THE ORDERED COLUMNS.
!     AMONG ALL THE COLUMNS OF MAXIMAL INCIDENCE, IDO CHOOSES A
!     COLUMN OF MAXIMAL DEGREE.
!     THE SUBROUTINE STATEMENT IS
!       SUBROUTINE IDO(M,N,INDROW,JPNTR,INDCOL,IPNTR,NDEG,LIST,
!                      MAXCLQ,IWA1,IWA2,IWA3,IWA4)
!     WHERE
!       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF ROWS OF A.
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF COLUMNS OF A.
!       INDROW IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW
!         INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!       JPNTR IS AN INTEGER INPUT ARRAY OF LENGTH N + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW.
!         THE ROW INDICES FOR COLUMN J ARE
!               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1.
!         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!       INDCOL IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE
!         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!       IPNTR IS AN INTEGER INPUT ARRAY OF LENGTH M + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL.
!         THE COLUMN INDICES FOR ROW I ARE
!               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1.
!         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!       NDEG IS AN INTEGER INPUT ARRAY OF LENGTH N WHICH SPECIFIES
!         THE DEGREE SEQUENCE. THE DEGREE OF THE J-TH COLUMN
!         OF A IS NDEG(J).
!       LIST IS AN INTEGER OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES
!         THE INCIDENCE-DEGREE ORDERING OF THE COLUMNS OF A. THE J-TH
!         COLUMN IN THIS ORDER IS LIST(J).
!       MAXCLQ IS AN INTEGER OUTPUT VARIABLE SET TO THE SIZE
!         OF THE LARGEST CLIQUE FOUND DURING THE ORDERING.
!       IWA1,IWA2,IWA3, AND IWA4 ARE INTEGER WORK ARRAYS OF LENGTH N.
!     SUBPROGRAMS CALLED
!       MINPACK-SUPPLIED ... NUMSRT
!       INTRINSIC ... MAX
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983.
!     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE'

     implicit none

! .. Scalar Arguments ..
      integer :: M, MAXCLQ, N
! ..
! .. Array Arguments ..
      integer :: INDCOL(*), INDROW(*), IPNTR(M+1), IWA1(0:N-1), IWA2(N), &
        IWA3(N), IWA4(N), JPNTR(N+1), LIST(N), NDEG(N)
! ..
! .. Local Scalars ..
      integer :: IC, IP, IR, JCOL, JP, MAXINC, MAXLST, NCOMP, NUMINC, NUMLST, &
        NUMORD, NUMWGT
! ..
! .. External Subroutines ..
!     EXTERNAL NUMSRT
! ..
! .. Intrinsic Functions ..
      intrinsic MAX
! ..
! ..  FIRST EXECUTABLE STATEMENT IDO
! ..
!     SORT THE DEGREE SEQUENCE.
      call NUMSRT(N,N-1,NDEG,-1,IWA4,IWA2,IWA3)

!     INITIALIZATION BLOCK.

!     CREATE A DOUBLY-LINKED LIST TO ACCESS THE INCIDENCES OF THE
!     COLUMNS. THE POINTERS FOR THE LINKED LIST ARE AS FOLLOWS.

!     EACH UNORDERED COLUMN IC IS IN A LIST (THE INCIDENCE LIST)
!     OF COLUMNS WITH THE SAME INCIDENCE.

!     IWA1(NUMINC) IS THE FIRST COLUMN IN THE NUMINC LIST
!     UNLESS IWA1(NUMINC) = 0. IN THIS CASE THERE ARE
!     NO COLUMNS IN THE NUMINC LIST.

!     IWA2(IC) IS THE COLUMN BEFORE IC IN THE INCIDENCE LIST
!     UNLESS IWA2(IC) = 0. IN THIS CASE IC IS THE FIRST
!     COLUMN IN THIS INCIDENCE LIST.

!     IWA3(IC) IS THE COLUMN AFTER IC IN THE INCIDENCE LIST
!     UNLESS IWA3(IC) = 0. IN THIS CASE IC IS THE LAST
!     COLUMN IN THIS INCIDENCE LIST.

!     IF IC IS AN UN-ORDERED COLUMN, THEN LIST(IC) IS THE
!     INCIDENCE OF IC TO THE GRAPH INDUCED BY THE ORDERED
!     COLUMNS. IF JCOL IS AN ORDERED COLUMN, THEN LIST(JCOL)
!     IS THE INCIDENCE-DEGREE ORDER OF COLUMN JCOL.

      MAXINC = 0
      do JP = N, 1, -1
        IC = IWA4(JP)
        IWA1(N-JP) = 0
        IWA2(IC) = 0
        IWA3(IC) = IWA1(0)
        if (IWA1(0) > 0) IWA2(IWA1(0)) = IC
        IWA1(0) = IC
        IWA4(JP) = 0
        LIST(JP) = 0
      end do

!     DETERMINE THE MAXIMAL SEARCH LENGTH FOR THE LIST
!     OF COLUMNS OF MAXIMAL INCIDENCE.
      MAXLST = 0
      do IR = 1, M
        MAXLST = MAXLST + (IPNTR(IR+1)-IPNTR(IR))**2
      end do
      MAXLST = MAXLST/N
      MAXCLQ = 0
      NUMORD = 1

!     BEGINNING OF ITERATION LOOP.

30    continue

!     UPDATE THE SIZE OF THE LARGEST CLIQUE
!     FOUND DURING THE ORDERING.
      if (MAXINC == 0) NCOMP = 0
      NCOMP = NCOMP + 1
      if (MAXINC+1 == NCOMP) MAXCLQ = max(MAXCLQ,NCOMP)

!     CHOOSE A COLUMN JCOL OF MAXIMAL DEGREE AMONG THE
!     COLUMNS OF MAXIMAL INCIDENCE MAXINC.
40    continue
      JP = IWA1(MAXINC)
      if (JP > 0) goto 50
      MAXINC = MAXINC - 1
      goto 40
50    continue
      NUMWGT = -1
      do NUMLST = 1, MAXLST
        if (NDEG(JP) > NUMWGT) then
          NUMWGT = NDEG(JP)
          JCOL = JP
        end if
        JP = IWA3(JP)
        if (JP <= 0) goto 70
      end do
70    continue
      LIST(JCOL) = NUMORD
      NUMORD = NUMORD + 1

!     TERMINATION TEST.
      if (NUMORD > N) goto 100

!     DELETE COLUMN JCOL FROM THE MAXINC LIST.
      if (IWA2(JCOL) == 0) then
        IWA1(MAXINC) = IWA3(JCOL)
      else
        IWA3(IWA2(JCOL)) = IWA3(JCOL)
      end if
      if (IWA3(JCOL) > 0) IWA2(IWA3(JCOL)) = IWA2(JCOL)

!     FIND ALL COLUMNS ADJACENT TO COLUMN JCOL.
      IWA4(JCOL) = N

!     DETERMINE ALL POSITIONS(IR,JCOL) WHICH CORRESPOND
!     TO NON-ZEROES IN THE MATRIX.
      do JP = JPNTR(JCOL), JPNTR(JCOL+1) - 1
        IR = INDROW(JP)
!       FOR EACH ROW IR, DETERMINE ALL POSITIONS(IR,IC)
!       WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX.
        do IP = IPNTR(IR), IPNTR(IR+1) - 1
          IC = INDCOL(IP)
!         ARRAY IWA4 MARKS COLUMNS WHICH ARE ADJACENT TO
!         COLUMN JCOL.
          if (IWA4(IC) < NUMORD) then
            IWA4(IC) = NUMORD
!           UPDATE THE POINTERS TO THE CURRENT INCIDENCE LISTS.
            NUMINC = LIST(IC)
            LIST(IC) = LIST(IC) + 1
            MAXINC = max(MAXINC,LIST(IC))
!           DELETE COLUMN IC FROM THE NUMINC LIST.
            if (IWA2(IC) == 0) then
              IWA1(NUMINC) = IWA3(IC)
            else
              IWA3(IWA2(IC)) = IWA3(IC)
            end if
            if (IWA3(IC) > 0) IWA2(IWA3(IC)) = IWA2(IC)
!           ADD COLUMN IC TO THE NUMINC+1 LIST.
            IWA2(IC) = 0
            IWA3(IC) = IWA1(NUMINC+1)
            if (IWA1(NUMINC+1) > 0) IWA2(IWA1(NUMINC+1)) = IC
            IWA1(NUMINC+1) = IC
          end if
        end do
      end do
!     END OF ITERATION LOOP.
      goto 30
100   continue
!     INVERT THE ARRAY LIST.
      do JCOL = 1, N
        IWA2(LIST(JCOL)) = JCOL
      end do
      LIST(1:N) = IWA2(1:N)
      return

    end subroutine IDO
!_______________________________________________________________________

    subroutine NUMSRT(N,NMAX,NUM,MODE,INDEX,LAST,NEXT)

!     GIVEN A SEQUENCE OF INTEGERS, THIS SUBROUTINE GROUPS TOGETHER THOSE
!     INDICES WITH THE SAME SEQUENCE VALUE AND, OPTIONALLY, SORTS THE
!     SEQUENCE INTO EITHER ASCENDING OR DESCENDING ORDER. THE SEQUENCE
!     OF INTEGERS IS DEFINED BY THE ARRAY NUM, AND IT IS ASSUMED THAT THE
!     INTEGERS ARE EACH FROM THE SET 0,1,...,NMAX. ON OUTPUT THE INDICES
!     K SUCH THAT NUM(K) = L FOR ANY L = 0,1,...,NMAX CAN BE OBTAINED
!     FROM THE ARRAYS LAST AND NEXT AS FOLLOWS.
!           K = LAST(L)
!           WHILE(K /= 0) K = NEXT(K)
!     OPTIONALLY, THE SUBROUTINE PRODUCES AN ARRAY INDEX SO THAT
!     THE SEQUENCE NUM(INDEX(I)), I = 1,2,...,N IS SORTED.
!     THE SUBROUTINE STATEMENT IS
!       SUBROUTINE NUMSRT(N,NMAX,NUM,MODE,INDEX,LAST,NEXT)
!     WHERE
!       N IS A POSITIVE INTEGER INPUT VARIABLE.
!       NMAX IS A POSITIVE INTEGER INPUT VARIABLE.
!       NUM IS AN INPUT ARRAY OF LENGTH N WHICH CONTAINS THE
!         SEQUENCE OF INTEGERS TO BE GROUPED AND SORTED. IT
!         IS ASSUMED THAT THE INTEGERS ARE EACH FROM THE SET
!         0,1,...,NMAX.
!       MODE IS AN INTEGER INPUT VARIABLE. THE SEQUENCE NUM IS
!         SORTED IN ASCENDING ORDER IF MODE IS POSITIVE AND IN
!         DESCENDING ORDER IF MODE IS NEGATIVE. IF MODE IS 0,
!         NO SORTING IS DONE.
!       INDEX IS AN INTEGER OUTPUT ARRAY OF LENGTH N SET SO
!         THAT THE SEQUENCE
!               NUM(INDEX(I)), I = 1,2,...,N
!         IS SORTED ACCORDING TO THE SETTING OF MODE. IF MODE
!         IS 0, INDEX IS NOT REFERENCED.
!       LAST IS AN INTEGER OUTPUT ARRAY OF LENGTH NMAX + 1. THE
!         INDEX OF NUM FOR THE LAST OCCURRENCE OF L IS LAST(L)
!         FOR ANY L = 0,1,...,NMAX UNLESS LAST(L) = 0. IN
!         THIS CASE L DOES NOT APPEAR IN NUM.
!       NEXT IS AN INTEGER OUTPUT ARRAY OF LENGTH N. IF
!         NUM(K) = L, THEN THE INDEX OF NUM FOR THE PREVIOUS
!         OCCURRENCE OF L IS NEXT(K) FOR ANY L = 0,1,...,NMAX
!         UNLESS NEXT(K) = 0. IN THIS CASE THERE IS NO PREVIOUS
!         OCCURRENCE OF L IN NUM.
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983.
!     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE'

      implicit none

! .. Scalar Arguments ..
      integer :: MODE, N, NMAX
! ..
! .. Array Arguments ..
      integer :: index(N), LAST(0:NMAX), NEXT(N), NUM(N)
! ..
! .. Local Scalars ..
      integer :: I, J, JINC, JL, JU, K, L
! ..
! ..  FIRST EXECUTABLE STATEMENT NUMSRT
! ..
!     DETERMINE THE ARRAYS NEXT AND LAST.
      LAST(0:NMAX) = 0
      do K = 1, N
        L = NUM(K)
        NEXT(K) = LAST(L)
        LAST(L) = K
      end do
      if (MODE == 0) return

!     STORE THE POINTERS TO THE SORTED ARRAY IN INDEX.
      I = 1
      if (MODE > 0) then
        JL = 0
        JU = NMAX
        JINC = 1
      else
        JL = NMAX
        JU = 0
        JINC = -1
      end if
      do J = JL, JU, JINC
        K = LAST(J)
30      continue
        if (K == 0) goto 40
        index(I) = K
        I = I + 1
        K = NEXT(K)
        goto 30
40      continue
      end do
      return

    end subroutine NUMSRT
!_______________________________________________________________________

    subroutine SEQ(N,INDROW,JPNTR,INDCOL,IPNTR,LIST,NGRP,MAXGRP,IWA)

!     GIVEN THE SPARSITY PATTERN OF AN M BY N MATRIX A, THIS
!     SUBROUTINE DETERMINES A CONSISTENT PARTITION OF THE
!     COLUMNS OF A BY A SEQUENTIAL ALGORITHM.
!     A CONSISTENT PARTITION IS DEFINED IN TERMS OF THE LOOPLESS
!     GRAPH G WITH VERTICES A(J), J = 1,2,...,N WHERE A(J) IS THE
!     J-TH COLUMN OF A AND WITH EDGE(A(I),A(J)) IF AND ONLY IF
!     COLUMNS I AND J HAVE A NON-ZERO IN THE SAME ROW POSITION.
!     A PARTITION OF THE COLUMNS OF A INTO GROUPS IS CONSISTENT
!     IF THE COLUMNS IN ANY GROUP ARE NOT ADJACENT IN THE GRAPH G.
!     IN GRAPH-THEORY TERMINOLOGY, A CONSISTENT PARTITION OF THE
!     COLUMNS OF A CORRESPONDS TO A COLORING OF THE GRAPH G.
!     THE SUBROUTINE EXAMINES THE COLUMNS IN THE ORDER SPECIFIED
!     BY THE ARRAY LIST, AND ASSIGNS THE CURRENT COLUMN TO THE
!     GROUP WITH THE SMALLEST POSSIBLE NUMBER.
!     NOTE THAT THE VALUE OF M IS NOT NEEDED BY SEQ AND IS
!     THEREFORE NOT PRESENT IN THE SUBROUTINE STATEMENT.
!     THE SUBROUTINE STATEMENT IS
!     SUBROUTINE SEQ(N,INDROW,JPNTR,INDCOL,IPNTR,LIST,NGRP,MAXGRP,IWA)
!     WHERE
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF COLUMNS OF A.
!       INDROW IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW
!         INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!       JPNTR IS AN INTEGER INPUT ARRAY OF LENGTH N + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW.
!         THE ROW INDICES FOR COLUMN J ARE
!               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1.
!         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!       INDCOL IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE
!         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!       IPNTR IS AN INTEGER INPUT ARRAY OF LENGTH M + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL.
!         THE COLUMN INDICES FOR ROW I ARE
!               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1.
!         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!       LIST IS AN INTEGER INPUT ARRAY OF LENGTH N WHICH SPECIFIES
!         THE ORDER TO BE USED BY THE SEQUENTIAL ALGORITHM.
!         THE J-TH COLUMN IN THIS ORDER IS LIST(J).
!       NGRP IS AN INTEGER OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES
!         THE PARTITION OF THE COLUMNS OF A. COLUMN JCOL BELONGS
!         TO GROUP NGRP(JCOL).
!       MAXGRP IS AN INTEGER OUTPUT VARIABLE WHICH SPECIFIES THE
!         NUMBER OF GROUPS IN THE PARTITION OF THE COLUMNS OF A.
!       IWA IS AN INTEGER WORK ARRAY OF LENGTH N.
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983.
!     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE'

     implicit none

! .. Scalar Arguments ..
      integer :: MAXGRP, N
! ..
! .. Array Arguments ..
      integer :: INDCOL(*), INDROW(*), IPNTR(*), IWA(N), JPNTR(N+1), &
        LIST(N), NGRP(N)
! ..
! .. Local Scalars ..
      integer :: IC, IP, IR, J, JCOL, JP
! ..
! ..  FIRST EXECUTABLE STATEMENT SEQ
! ..
!     INITIALIZATION BLOCK.
      MAXGRP = 0
      NGRP(1:N) = N
      IWA(1:N) = 0

!     BEGINNING OF ITERATION LOOP.

      do J = 1, N
        JCOL = LIST(J)
!       FIND ALL COLUMNS ADJACENT TO COLUMN JCOL.
!       DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND
!       TO NON-ZEROES IN THE MATRIX.
        do JP = JPNTR(JCOL), JPNTR(JCOL+1) - 1
          IR = INDROW(JP)
!         FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC)
!         WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX.
          do IP = IPNTR(IR), IPNTR(IR+1) - 1
            IC = INDCOL(IP)
!           ARRAY IWA MARKS THE GROUP NUMBERS OF THE
!           COLUMNS WHICH ARE ADJACENT TO COLUMN JCOL.
            IWA(NGRP(IC)) = J
          end do
        end do
!       ASSIGN THE SMALLEST UN-MARKED GROUP NUMBER TO JCOL.
        do JP = 1, MAXGRP
          if (IWA(JP)/=J) goto 50
        end do
        MAXGRP = MAXGRP + 1
50      continue
        NGRP(JCOL) = JP
      end do

!     END OF ITERATION LOOP.

      return

    end subroutine SEQ
!_______________________________________________________________________

    subroutine SETR(M,N,INDROW,JPNTR,INDCOL,IPNTR,IWA)

!     GIVEN A COLUMN-ORIENTED DEFINITION OF THE SPARSITY PATTERN
!     OF AN M BY N MATRIX A, THIS SUBROUTINE DETERMINES A
!     ROW-ORIENTED DEFINITION OF THE SPARSITY PATTERN OF A.
!     ON INPUT THE COLUMN-ORIENTED DEFINITION IS SPECIFIED BY
!     THE ARRAYS INDROW AND JPNTR. ON OUTPUT THE ROW-ORIENTED
!     DEFINITION IS SPECIFIED BY THE ARRAYS INDCOL AND IPNTR.
!     THE SUBROUTINE STATEMENT IS
!       SUBROUTINE SETR(M,N,INDROW,JPNTR,INDCOL,IPNTR,IWA)
!     WHERE
!       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF ROWS OF A.
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF COLUMNS OF A.
!       INDROW IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW
!         INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!       JPNTR IS AN INTEGER INPUT ARRAY OF LENGTH N + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW.
!         THE ROW INDICES FOR COLUMN J ARE
!               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1.
!         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!       INDCOL IS AN INTEGER OUTPUT ARRAY WHICH CONTAINS THE
!         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!       IPNTR IS AN INTEGER OUTPUT ARRAY OF LENGTH M + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL.
!         THE COLUMN INDICES FOR ROW I ARE
!               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1.
!         NOTE THAT IPNTR(1) IS SET TO 1 AND THAT IPNTR(M+1)-1 IS
!         THEN THE NUMBER OF NON-ZERO ELEMENTS OF THE MATRIX A.
!       IWA IS AN INTEGER WORK ARRAY OF LENGTH M.
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983.
!     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE'

     implicit none

! .. Scalar Arguments ..
      integer :: M, N
! ..
! .. Array Arguments ..
      integer :: INDCOL(*), INDROW(*), IPNTR(M+1), IWA(M), JPNTR(N+1)
! ..
! .. Local Scalars ..
      integer :: IR, JCOL, JP
! ..
! ..  FIRST EXECUTABLE STATEMENT SETR
! ..
!     STORE IN ARRAY IWA THE COUNTS OF NON-ZEROES IN THE ROWS.
      IWA(1:M) = 0
      do JP = 1, JPNTR(N+1) - 1
        IWA(INDROW(JP)) = IWA(INDROW(JP)) + 1
      end do

!     SET POINTERS TO THE START OF THE ROWS IN INDCOL.
      IPNTR(1) = 1
      do IR = 1, M
        IPNTR(IR+1) = IPNTR(IR) + IWA(IR)
        IWA(IR) = IPNTR(IR)
      end do

!     FILL INDCOL.
      do JCOL = 1, N
        do JP = JPNTR(JCOL), JPNTR(JCOL+1) - 1
          IR = INDROW(JP)
          INDCOL(IWA(IR)) = JCOL
          IWA(IR) = IWA(IR) + 1
        end do
      end do
      return

    end subroutine SETR
!_______________________________________________________________________

    subroutine SLO(N,INDROW,JPNTR,INDCOL,IPNTR,NDEG,LIST,MAXCLQ,IWA1,IWA2, &
        IWA3,IWA4)

!     GIVEN THE SPARSITY PATTERN OF AN M BY N MATRIX A, THIS
!     SUBROUTINE DETERMINES THE SMALLEST-LAST ORDERING OF THE
!     COLUMNS OF A.
!     THE SMALLEST-LAST ORDERING IS DEFINED FOR THE LOOPLESS
!     GRAPH G WITH VERTICES A(J), J = 1,2,...,N WHERE A(J) IS THE
!     J-TH COLUMN OF A AND WITH EDGE(A(I),A(J)) IF AND ONLY IF
!     COLUMNS I AND J HAVE A NON-ZERO IN THE SAME ROW POSITION.
!     THE SMALLEST-LAST ORDERING IS DETERMINED RECURSIVELY BY
!     LETTING LIST(K), K = N,...,1 BE A COLUMN WITH LEAST DEGREE
!     IN THE SUBGRAPH SPANNED BY THE UN-ORDERED COLUMNS.
!     NOTE THAT THE VALUE OF M IS NOT NEEDED BY SLO AND IS
!     THEREFORE NOT PRESENT IN THE SUBROUTINE STATEMENT.
!     THE SUBROUTINE STATEMENT IS
!     SUBROUTINE SLO(N,INDROW,JPNTR,INDCOL,IPNTR,NDEG,LIST, &
!                    MAXCLQ,IWA1,IWA2,IWA3,IWA4)
!     WHERE
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF COLUMNS OF A.
!       INDROW IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW
!         INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!       JPNTR IS AN INTEGER INPUT ARRAY OF LENGTH N + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW.
!         THE ROW INDICES FOR COLUMN J ARE
!               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1.
!         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!       INDCOL IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE
!         COLUMN INDICES FOR THE NON-ZEROES IN THE MATRIX A.
!       IPNTR IS AN INTEGER INPUT ARRAY OF LENGTH M + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL.
!         THE COLUMN INDICES FOR ROW I ARE
!               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1.
!         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!       NDEG IS AN INTEGER INPUT ARRAY OF LENGTH N WHICH SPECIFIES
!         THE DEGREE SEQUENCE. THE DEGREE OF THE J-TH COLUMN
!         OF A IS NDEG(J).
!       LIST IS AN INTEGER OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES
!         THE SMALLEST-LAST ORDERING OF THE COLUMNS OF A. THE J-TH
!         COLUMN IN THIS ORDER IS LIST(J).
!       MAXCLQ IS AN INTEGER OUTPUT VARIABLE SET TO THE SIZE
!         OF THE LARGEST CLIQUE FOUND DURING THE ORDERING.
!       IWA1,IWA2,IWA3, AND IWA4 ARE INTEGER WORK ARRAYS OF LENGTH N.
!     SUBPROGRAMS CALLED
!       INTRINSIC ... MIN
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983.
!     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE'

     implicit none

! .. Scalar Arguments ..
      integer :: MAXCLQ, N
! ..
! .. Array Arguments ..
      integer :: INDCOL(*), INDROW(*), IPNTR(*), IWA1(0:N-1), IWA2(N), &
        IWA3(N), IWA4(N), JPNTR(N+1), LIST(N), NDEG(N)
! ..
! .. Local Scalars ..
      integer :: IC, IP, IR, JCOL, JP, MINDEG, NUMDEG, NUMORD
! ..
! .. Intrinsic Functions ..
      intrinsic MIN
! ..
! ..  FIRST EXECUTABLE STATEMENT SLO
! ..
!     INITIALIZATION BLOCK.
      MINDEG = N
      do JP = 1, N
        IWA1(JP-1) = 0
        IWA4(JP) = N
        LIST(JP) = NDEG(JP)
        MINDEG = min(MINDEG,NDEG(JP))
      end do

!     CREATE A DOUBLY-LINKED LIST TO ACCESS THE DEGREES OF THE
!     COLUMNS. THE POINTERS FOR THE LINKED LIST ARE AS FOLLOWS.

!     EACH UN-ORDERED COLUMN IC IS IN A LIST (THE DEGREE LIST)
!     OF COLUMNS WITH THE SAME DEGREE.

!     IWA1(NUMDEG) IS THE FIRST COLUMN IN THE NUMDEG LIST
!     UNLESS IWA1(NUMDEG) = 0. IN THIS CASE THERE ARE
!     NO COLUMNS IN THE NUMDEG LIST.

!     IWA2(IC) IS THE COLUMN BEFORE IC IN THE DEGREE LIST
!     UNLESS IWA2(IC) = 0. IN THIS CASE IC IS THE FIRST
!     COLUMN IN THIS DEGREE LIST.

!     IWA3(IC) IS THE COLUMN AFTER IC IN THE DEGREE LIST
!     UNLESS IWA3(IC) = 0. IN THIS CASE IC IS THE LAST
!     COLUMN IN THIS DEGREE LIST.

!     IF IC IS AN UN-ORDERED COLUMN, THEN LIST(IC) IS THE
!     DEGREE OF IC IN THE GRAPH INDUCED BY THE UN-ORDERED
!     COLUMNS. IF JCOL IS AN ORDERED COLUMN, THEN LIST(JCOL)
!     IS THE SMALLEST-LAST ORDER OF COLUMN JCOL.

      do JP = 1, N
        NUMDEG = NDEG(JP)
        IWA2(JP) = 0
        IWA3(JP) = IWA1(NUMDEG)
        if (IWA1(NUMDEG) > 0) IWA2(IWA1(NUMDEG)) = JP
        IWA1(NUMDEG) = JP
      end do
      MAXCLQ = 0
      NUMORD = N

!     BEGINNING OF ITERATION LOOP.

30    continue

!     MARK THE SIZE OF THE LARGEST CLIQUE
!     FOUND DURING THE ORDERING.
      if (MINDEG+1 == NUMORD .and. MAXCLQ == 0) MAXCLQ = NUMORD

!     CHOOSE A COLUMN JCOL OF MINIMAL DEGREE MINDEG.
40    continue
      JCOL = IWA1(MINDEG)
      if (JCOL > 0) goto 50
      MINDEG = MINDEG + 1
      goto 40
50    continue
      LIST(JCOL) = NUMORD
      NUMORD = NUMORD - 1

!     TERMINATION TEST.
      if (NUMORD == 0) goto 80

!     DELETE COLUMN JCOL FROM THE MINDEG LIST.
      IWA1(MINDEG) = IWA3(JCOL)
      if (IWA3(JCOL) > 0) IWA2(IWA3(JCOL)) = 0

!     FIND ALL COLUMNS ADJACENT TO COLUMN JCOL.
      IWA4(JCOL) = 0

!     DETERMINE ALL POSITIONS (IR,JCOL) WHICH CORRESPOND
!     TO NON-ZEROES IN THE MATRIX.
      do JP = JPNTR(JCOL), JPNTR(JCOL+1) - 1
        IR = INDROW(JP)
!       FOR EACH ROW IR, DETERMINE ALL POSITIONS (IR,IC)
!       WHICH CORRESPOND TO NON-ZEROES IN THE MATRIX.
        do IP = IPNTR(IR), IPNTR(IR+1) - 1
          IC = INDCOL(IP)
!         ARRAY IWA4 MARKS COLUMNS WHICH ARE ADJACENT TO
!         COLUMN JCOL.

          if (IWA4(IC) > NUMORD) then
            IWA4(IC) = NUMORD
!           UPDATE THE POINTERS TO THE CURRENT DEGREE LISTS.
            NUMDEG = LIST(IC)
            LIST(IC) = LIST(IC) - 1
            MINDEG = min(MINDEG,LIST(IC))
!           DELETE COLUMN IC FROM THE NUMDEG LIST.
            if (IWA2(IC) == 0) then
              IWA1(NUMDEG) = IWA3(IC)
            else
              IWA3(IWA2(IC)) = IWA3(IC)
            end if
            if (IWA3(IC) > 0) IWA2(IWA3(IC)) = IWA2(IC)
!           ADD COLUMN IC TO THE NUMDEG-1 LIST.
            IWA2(IC) = 0
            IWA3(IC) = IWA1(NUMDEG-1)
            if (IWA1(NUMDEG-1) > 0) IWA2(IWA1(NUMDEG-1)) = IC
            IWA1(NUMDEG-1) = IC
          end if
        end do
      end do
!     END OF ITERATION LOOP.

      goto 30
80    continue

!     INVERT THE ARRAY LIST.
      do JCOL = 1, N
        IWA2(LIST(JCOL)) = JCOL
      end do
      LIST(1:N) = IWA2(1:N)
      return

    end subroutine SLO
!_______________________________________________________________________

    subroutine SRTDAT(N,NNZ,INDROW,INDCOL,JPNTR,IWA)

!     GIVEN THE NON-ZERO ELEMENTS OF AN M BY N MATRIX A IN ARBITRARY
!     ORDER AS SPECIFIED BY THEIR ROW AND COLUMN INDICES, THIS SUBROUTINE
!     PERMUTES THESE ELEMENTS SO THAT THEIR COLUMN INDICES ARE IN
!     NON-DECREASING ORDER. ON INPUT IT IS ASSUMED THAT THE ELEMENTS ARE
!     SPECIFIED IN
!           INDROW(K),INDCOL(K), K = 1,...,NNZ.
!     ON OUTPUT THE ELEMENTS ARE PERMUTED SO THAT INDCOL IS IN
!     NON-DECREASING ORDER. IN ADDITION, THE ARRAY JPNTR IS SET SO THAT
!     THE ROW INDICES FOR COLUMN J ARE
!           INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1.
!     NOTE THAT THE VALUE OF M IS NOT NEEDED BY SRTDAT AND IS THEREFORE
!     NOT PRESENT IN THE SUBROUTINE STATEMENT.
!     THE SUBROUTINE STATEMENT IS
!       SUBROUTINE SRTDAT(N,NNZ,INDROW,INDCOL,JPNTR,IWA)
!     WHERE
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF COLUMNS OF A.
!       NNZ IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF NON-ZERO ELEMENTS OF A.
!       INDROW IS AN INTEGER ARRAY OF LENGTH NNZ. ON INPUT INDROW
!         MUST CONTAIN THE ROW INDICES OF THE NON-ZERO ELEMENTS OF A.
!         ON OUTPUT INDROW IS PERMUTED SO THAT THE CORRESPONDING
!         COLUMN INDICES OF INDCOL ARE IN NON-DECREASING ORDER.
!       INDCOL IS AN INTEGER ARRAY OF LENGTH NNZ. ON INPUT INDCOL
!         MUST CONTAIN THE COLUMN INDICES OF THE NON-ZERO ELEMENTS
!         OF A. ON OUTPUT INDCOL IS PERMUTED SO THAT THESE INDICES
!         ARE IN NON-DECREASING ORDER.
!       JPNTR IS AN INTEGER OUTPUT ARRAY OF LENGTH N + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN THE OUTPUT
!         INDROW. THE ROW INDICES FOR COLUMN J ARE
!               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1.
!         NOTE THAT JPNTR(1) IS SET TO 1 AND THAT JPNTR(N+1)-1
!         IS THEN NNZ.
!       IWA IS AN INTEGER WORK ARRAY OF LENGTH N.
!     SUBPROGRAMS CALLED - NONE
!     INTRINSIC - MAX
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983.
!     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE'

     implicit none

! .. Scalar Arguments ..
      integer :: N, NNZ
! ..
! .. Array Arguments ..
      integer :: INDCOL(NNZ), INDROW(NNZ), IWA(N), JPNTR(N+1)
! ..
! .. Local Scalars ..
      integer :: I, J, K, L
! ..
! .. Intrinsic Functions ..
      intrinsic MAX
! ..
! ..  FIRST EXECUTABLE STATEMENT SRTDAT
! ..
!     STORE IN ARRAY IWA THE COUNTS OF NON-ZEROES IN THE COLUMNS.
      IWA(1:N) = 0
      do K = 1, NNZ
        IWA(INDCOL(K)) = IWA(INDCOL(K)) + 1
      end do

!     SET POINTERS TO THE START OF THE COLUMNS IN INDROW.
      JPNTR(1) = 1
      do J = 1, N
        JPNTR(J+1) = JPNTR(J) + IWA(J)
        IWA(J) = JPNTR(J)
      end do
      K = 1

!     BEGIN IN-PLACE SORT.
40    continue
      J = INDCOL(K)
      if (K >= JPNTR(J)) then
!        CURRENT ELEMENT IS IN POSITION. NOW EXAMINE THE
!        NEXT ELEMENT OR THE FIRST UN-SORTED ELEMENT IN
!        THE J-TH GROUP.
        K = max(K+1,IWA(J))
      else
!       CURRENT ELEMENT IS NOT IN POSITION. PLACE ELEMENT
!       IN POSITION AND MAKE THE DISPLACED ELEMENT THE
!       CURRENT ELEMENT.
        L = IWA(J)
        IWA(J) = IWA(J) + 1
        I = INDROW(K)
        INDROW(K) = INDROW(L)
        INDCOL(K) = INDCOL(L)
        INDROW(L) = I
        INDCOL(L) = J
      end if
      if (K <= NNZ) goto 40
      return

    end subroutine SRTDAT
!_______________________________________________________________________

    subroutine FDJS(M,N,COL,IND,NPNTR,NGRP,NUMGRP,D,FJACD,FJAC)

!     GIVEN A CONSISTENT PARTITION OF THE COLUMNS OF AN M BY N
!     JACOBIAN MATRIX INTO GROUPS, THIS SUBROUTINE COMPUTES
!     APPROXIMATIONS TO THOSE COLUMNS IN A GIVEN GROUP. THE
!     APPROXIMATIONS ARE STORED INTO EITHER A COLUMN-ORIENTED
!     OR A ROW-ORIENTED PATTERN.
!     A PARTITION IS CONSISTENT IF THE COLUMNS IN ANY GROUP
!     DO NOT HAVE A NON-ZERO IN THE SAME ROW POSITION.
!     APPROXIMATIONS TO THE COLUMNS OF THE JACOBIAN MATRIX IN A
!     GIVEN GROUP CAN BE OBTAINED BY SPECIFYING A DIFFERENCE
!     PARAMETER ARRAY D WITH D(JCOL) NON-ZERO IF AND ONLY IF
!     JCOL IS A COLUMN IN THE GROUP, AND AN APPROXIMATION TO
!     JAC*D WHERE JAC DENOTES THE JACOBIAN MATRIX OF A MAPPING F.
!     D CAN BE DEFINED WITH THE FOLLOWING SEGMENT OF CODE.
!           DO 10 JCOL = 1, N
!             D(JCOL) = 0.0
!             IF (NGRP(JCOL) == NUMGRP) D(JCOL) = ETA(JCOL)
!        10 CONTINUE
!     IN THE ABOVE CODE NUMGRP IS THE GIVEN GROUP NUMBER,
!     NGRP(JCOL) IS THE GROUP NUMBER OF COLUMN JCOL, AND
!     ETA(JCOL) IS THE DIFFERENCE PARAMETER USED TO
!     APPROXIMATE COLUMN JCOL OF THE JACOBIAN MATRIX.
!     SUITABLE VALUES FOR THE ARRAY ETA MUST BE PROVIDED.
!     AS MENTIONED ABOVE, AN APPROXIMATION TO JAC*D MUST
!     ALSO BE PROVIDED. FOR EXAMPLE, THE APPROXIMATION
!           F(X+D) - F(X)
!     CORRESPONDS TO THE FORWARD DIFFERENCE FORMULA AT X.
!     THE SUBROUTINE STATEMENT IS
!       SUBROUTINE FDJS(M,N,COL,IND,NPNTR,NGRP,NUMGRP,D,FJACD,FJAC)
!     WHERE
!       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF ROWS OF THE JACOBIAN MATRIX.
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF COLUMNS OF THE JACOBIAN MATRIX.
!       COL IS A LOGICAL INPUT VARIABLE. IF COL IS SET TRUE, THEN THE
!         JACOBIAN APPROXIMATIONS ARE STORED INTO A COLUMN-ORIENTED
!         PATTERN. IF COL IS SET FALSE, THEN THE JACOBIAN
!         APPROXIMATIONS ARE STORED INTO A ROW-ORIENTED PATTERN.
!       IND IS AN INTEGER INPUT ARRAY WHICH CONTAINS THE ROW
!         INDICES FOR THE NON-ZEROES IN THE JACOBIAN MATRIX
!         IF COL IS TRUE, AND CONTAINS THE COLUMN INDICES FOR
!         THE NON-ZEROES IN THE JACOBIAN MATRIX IF COL IS FALSE.
!       NPNTR IS AN INTEGER INPUT ARRAY WHICH SPECIFIES THE
!         LOCATIONS OF THE ROW INDICES IN IND IF COL IS TRUE, AND
!         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN IND IF
!         COL IS FALSE. IF COL IS TRUE, THE INDICES FOR COLUMN J ARE
!               IND(K), K = NPNTR(J),...,NPNTR(J+1)-1.
!         IF COL IS FALSE, THE INDICES FOR ROW I ARE
!               IND(K), K = NPNTR(I),...,NPNTR(I+1)-1.
!         NOTE THAT NPNTR(N+1)-1 IF COL IS TRUE, OR NPNTR(M+1)-1
!         IF COL IS FALSE, IS THEN THE NUMBER OF NON-ZERO ELEMENTS
!         OF THE JACOBIAN MATRIX.
!       NGRP IS AN INTEGER INPUT ARRAY OF LENGTH N WHICH SPECIFIES
!         THE PARTITION OF THE COLUMNS OF THE JACOBIAN MATRIX.
!         COLUMN JCOL BELONGS TO GROUP NGRP(JCOL).
!       NUMGRP IS A POSITIVE INTEGER INPUT VARIABLE SET TO A GROUP
!         NUMBER IN THE PARTITION. THE COLUMNS OF THE JACOBIAN
!         MATRIX IN THIS GROUP ARE TO BE ESTIMATED ON THIS CALL.
!       D IS AN INPUT ARRAY OF LENGTH N WHICH CONTAINS THE
!         DIFFERENCE PARAMETER VECTOR FOR THE ESTIMATE OF
!         THE JACOBIAN MATRIX COLUMNS IN GROUP NUMGRP.
!       FJACD IS AN INPUT ARRAY OF LENGTH M WHICH CONTAINS
!         AN APPROXIMATION TO THE DIFFERENCE VECTOR JAC*D,
!         WHERE JAC DENOTES THE JACOBIAN MATRIX.
!       FJAC IS AN OUTPUT ARRAY OF LENGTH NNZ, WHERE NNZ IS THE
!         NUMBER OF ITS NON-ZERO ELEMENTS. AT EACH CALL OF FDJS,
!         FJAC IS UPDATED TO INCLUDE THE NON-ZERO ELEMENTS OF THE
!         JACOBIAN MATRIX FOR THOSE COLUMNS IN GROUP NUMGRP. FJAC
!         SHOULD NOT BE ALTERED BETWEEN SUCCESSIVE CALLS TO FDJS.
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983.
!     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE'

     implicit none

! .. Parameters ..
      integer, parameter :: WP = kind(0.0D0)
! ..
! .. Scalar Arguments ..
      integer :: M, N, NUMGRP
      logical :: COL
! ..
! .. Array Arguments ..
      real (WP) :: D(N), FJAC(*), FJACD(M)
      integer :: IND(*), NGRP(N), NPNTR(*)
! ..
! .. Local Scalars ..
      integer :: IP, IROW, JCOL, JP
! ..
! .. Intrinsic Functions ..
!     INTRINSIC KIND
! ..
! ..  FIRST EXECUTABLE STATEMENT FDJS
! ..
!     COMPUTE ESTIMATES OF JACOBIAN MATRIX COLUMNS IN GROUP
!     NUMGRP. THE ARRAY FJACD MUST CONTAIN AN APPROXIMATION
!     TO JAC*D, WHERE JAC DENOTES THE JACOBIAN MATRIX AND D
!     IS A DIFFERENCE PARAMETER VECTOR WITH D(JCOL) NON-ZERO
!     IF AND ONLY IF JCOL IS A COLUMN IN GROUP NUMGRP.
      if (COL) then
!       COLUMN ORIENTATION.
        do JCOL = 1, N
          if (NGRP(JCOL) == NUMGRP) then
            do JP = NPNTR(JCOL), NPNTR(JCOL+1) - 1
              IROW = IND(JP)
              FJAC(JP) = FJACD(IROW)/D(JCOL)
            end do
          end if
        end do
      else
!       ROW ORIENTATION.
        do IROW = 1, M
          do IP = NPNTR(IROW), NPNTR(IROW+1) - 1
            JCOL = IND(IP)
            if (NGRP(JCOL) == NUMGRP) then
              FJAC(IP) = FJACD(IROW)/D(JCOL)
              goto 40
            end if
          end do
40        continue
        end do
      end if
      return

    end subroutine FDJS
!_______________________________________________________________________

    subroutine DVDSM(M,N,NPAIRS,INDROW,INDCOL,NGRP,MAXGRP,MINGRP,INFO, &
        IPNTR, JPNTR,IWA,LIWA)

!     THE PURPOSE OF DSM IS TO DETERMINE AN OPTIMAL OR NEAR-
!     OPTIMAL CONSISTENT PARTITION OF THE COLUMNS OF A SPARSE
!     M BY N MATRIX A.
!     THE SPARSITY PATTERN OF THE MATRIX A IS SPECIFIED BY
!     THE ARRAYS INDROW AND INDCOL. ON INPUT THE INDICES
!     FOR THE NON-ZERO ELEMENTS OF A ARE
!           INDROW(K),INDCOL(K), K = 1,2,...,NPAIRS.
!     THE(INDROW,INDCOL) PAIRS MAY BE SPECIFIED IN ANY ORDER.
!     DUPLICATE INPUT PAIRS ARE PERMITTED, BUT THE SUBROUTINE
!     ELIMINATES THEM.
!     THE SUBROUTINE PARTITIONS THE COLUMNS OF A INTO GROUPS
!     SUCH THAT COLUMNS IN THE SAME GROUP DO NOT HAVE A
!     NON-ZERO IN THE SAME ROW POSITION. A PARTITION OF THE
!     COLUMNS OF A WITH THIS PROPERTY IS CONSISTENT WITH THE
!     DIRECT DETERMINATION OF A.
!     THE SUBROUTINE STATEMENT IS
!       SUBROUTINE DVDSM(M,N,NPAIRS,INDROW,INDCOL,NGRP,MAXGRP,MINGRP,
!                      INFO,IPNTR,JPNTR,IWA,LIWA)
!     WHERE
!       M IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF ROWS OF A.
!       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!         OF COLUMNS OF A.
!       NPAIRS IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE
!         NUMBER OF (INDROW,INDCOL) PAIRS USED TO DESCRIBE THE
!         SPARSITY PATTERN OF A.
!       INDROW IS AN INTEGER ARRAY OF LENGTH NPAIRS. ON INPUT INDROW
!         MUST CONTAIN THE ROW INDICES OF THE NON-ZERO ELEMENTS OF A.
!         ON OUTPUT INDROW IS PERMUTED SO THAT THE CORRESPONDING
!         COLUMN INDICES ARE IN NON-DECREASING ORDER. THE COLUMN
!         INDICES CAN BE RECOVERED FROM THE ARRAY JPNTR.
!       INDCOL IS AN INTEGER ARRAY OF LENGTH NPAIRS. ON INPUT INDCOL
!         MUST CONTAIN THE COLUMN INDICES OF THE NON-ZERO ELEMENTS OF
!         A. ON OUTPUT INDCOL IS PERMUTED SO THAT THE CORRESPONDING
!         ROW INDICES ARE IN NON-DECREASING ORDER. THE ROW INDICES
!         CAN BE RECOVERED FROM THE ARRAY IPNTR.
!       NGRP IS AN INTEGER OUTPUT ARRAY OF LENGTH N WHICH SPECIFIES
!         THE PARTITION OF THE COLUMNS OF A. COLUMN JCOL BELONGS
!         TO GROUP NGRP(JCOL).
!       MAXGRP IS AN INTEGER OUTPUT VARIABLE WHICH SPECIFIES THE
!         NUMBER OF GROUPS IN THE PARTITION OF THE COLUMNS OF A.
!       MINGRP IS AN INTEGER OUTPUT VARIABLE WHICH SPECIFIES A LOWER
!         BOUND FOR THE NUMBER OF GROUPS IN ANY CONSISTENT PARTITION
!         OF THE COLUMNS OF A.
!       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS. FOR
!         NORMAL TERMINATION INFO = 1. IF M, N, OR NPAIRS IS NOT
!         POSITIVE OR LIWA IS LESS THAN MAX(M,6*N), THEN INFO = 0.
!         IF THE K-TH ELEMENT OF INDROW IS NOT AN INTEGER BETWEEN
!         1 AND M OR THE K-TH ELEMENT OF INDCOL IS NOT AN INTEGER
!         BETWEEN 1 AND N, THEN INFO = -K.
!       IPNTR IS AN INTEGER OUTPUT ARRAY OF LENGTH M + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE COLUMN INDICES IN INDCOL.
!         THE COLUMN INDICES FOR ROW I ARE
!               INDCOL(K), K = IPNTR(I),...,IPNTR(I+1)-1.
!         NOTE THAT IPNTR(M+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!       JPNTR IS AN INTEGER OUTPUT ARRAY OF LENGTH N + 1 WHICH
!         SPECIFIES THE LOCATIONS OF THE ROW INDICES IN INDROW.
!         THE ROW INDICES FOR COLUMN J ARE
!               INDROW(K), K = JPNTR(J),...,JPNTR(J+1)-1.
!         NOTE THAT JPNTR(N+1)-1 IS THEN THE NUMBER OF NON-ZERO
!         ELEMENTS OF THE MATRIX A.
!       IWA IS AN INTEGER WORK ARRAY OF LENGTH LIWA.
!       LIWA IS A POSITIVE INTEGER INPUT VARIABLE NOT LESS THAN
!         MAX(M,6*N).
!       MINPACK-SUPPLIED ... DEGR,IDO,NUMSRT,SEQ,SETR,SLO,SRTDAT
!       INTRINSIC - MAX
!     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JULY 1983.
!     THOMAS F. COLEMAN, BURTON S. GARBOW, JORGE J. MORE'

     implicit none

! .. Scalar Arguments ..
      integer :: INFO, LIWA, M, MAXGRP, MINGRP, N, NPAIRS
! ..
! .. Array Arguments ..
      integer :: INDCOL(NPAIRS), INDROW(NPAIRS), IPNTR(M+1), IWA(LIWA), &
        JPNTR(N+1), NGRP(N)
! ..
! .. Local Scalars ..
      integer :: I, IR, J, JP, K, MAXCLQ, NNZ, NUMGRP
! ..
! .. External Subroutines ..
!     EXTERNAL DEGR, IDO, NUMSRT, SEQ, SETR, SLO, SRTDAT
! ..
! .. Intrinsic Functions ..
      intrinsic MAX
! ..
! ..  FIRST EXECUTABLE STATEMENT DSM
! ..
!     CHECK THE INPUT DATA.
      INFO = 0
      if (M<1 .or. N<1 .or. NPAIRS<1 .or. LIWA<max(M,6*N)) return
      do K = 1, NPAIRS
        INFO = -K
        if (INDROW(K)<1 .or. INDROW(K)>M .or. INDCOL(K)<1 .or. INDCOL(K)>N) &
          return
      end do
      INFO = 1

!     SORT THE DATA STRUCTURE BY COLUMNS.
      call SRTDAT(N,NPAIRS,INDROW,INDCOL,JPNTR,IWA)

!     COMPRESS THE DATA AND DETERMINE THE NUMBER OF
!     NON-ZERO ELEMENTS OF A.
      IWA(1:M) =  0
      NNZ = 1
      do J = 1, N
        K = NNZ
        do JP = JPNTR(J), JPNTR(J+1) - 1
          IR = INDROW(JP)
          if (IWA(IR)/=J) then
            INDROW(NNZ) = IR
            NNZ = NNZ + 1
            IWA(IR) = J
          end if
        end do
        JPNTR(J) = K
      end do
      JPNTR(N+1) = NNZ

!     EXTEND THE DATA STRUCTURE TO ROWS.
      call SETR(M,N,INDROW,JPNTR,INDCOL,IPNTR,IWA)

!     DETERMINE A LOWER BOUND FOR THE NUMBER OF GROUPS.
      MINGRP = 0
      do I = 1, M
        MINGRP = max(MINGRP,IPNTR(I+1)-IPNTR(I))
      end do

!     DETERMINE THE DEGREE SEQUENCE FOR THE INTERSECTION
!     GRAPH OF THE COLUMNS OF A.
      call DEGR(N,INDROW,JPNTR,INDCOL,IPNTR,IWA(5*N+1),IWA(N+1))

!     COLOR THE INTERSECTION GRAPH OF THE COLUMNS OF A
!     WITH THE SMALLEST-LAST (SL) ORDERING.
      call SLO(N,INDROW,JPNTR,INDCOL,IPNTR,IWA(5*N+1),IWA(4*N+1),MAXCLQ, &
        IWA(1),IWA(N+1),IWA(2*N+1),IWA(3*N+1))
      call SEQ(N,INDROW,JPNTR,INDCOL,IPNTR,IWA(4*N+1),NGRP,MAXGRP,IWA(N+1))
      MINGRP = max(MINGRP,MAXCLQ)

!     EXIT IF THE SMALLEST-LAST ORDERING IS OPTIMAL.

      if (MAXGRP == MINGRP) return

!     COLOR THE INTERSECTION GRAPH OF THE COLUMNS OF A
!     WITH THE INCIDENCE-DEGREE(ID) ORDERING.
      call IDO(M,N,INDROW,JPNTR,INDCOL,IPNTR,IWA(5*N+1),IWA(4*N+1),MAXCLQ, &
        IWA(1),IWA(N+1),IWA(2*N+1),IWA(3*N+1))
      call SEQ(N,INDROW,JPNTR,INDCOL,IPNTR,IWA(4*N+1),IWA(1),NUMGRP,IWA(N+1))
      MINGRP = max(MINGRP,MAXCLQ)

!     RETAIN THE BETTER OF THE TWO ORDERINGS SO FAR.
      if (NUMGRP < MAXGRP) then
        MAXGRP = NUMGRP
        NGRP(1:N) = IWA(1:N)

!       EXIT IF THE INCIDENCE-DEGREE ORDERING IS OPTIMAL.
        if (MAXGRP == MINGRP) return
      end if

!     COLOR THE INTERSECTION GRAPH OF THE COLUMNS OF A
!     WITH THE LARGEST-FIRST (LF) ORDERING.
      call NUMSRT(N,N-1,IWA(5*N+1),-1,IWA(4*N+1),IWA(2*N+1),IWA(N+1))
      call SEQ(N,INDROW,JPNTR,INDCOL,IPNTR,IWA(4*N+1),IWA(1),NUMGRP,IWA(N+1))

!     RETAIN THE BEST OF THE THREE ORDERINGS AND EXIT.
      if (NUMGRP < MAXGRP) then
        MAXGRP = NUMGRP
        NGRP(1:N) = IWA(1:N)
      end if
      return

    end subroutine DVDSM
!_______________________________________________________________________

      subroutine DGROUPDS(N,MAXGRPDS,NGRPDS,IGP,JGP)
! ..
! Construct the column grouping arrays IGP anf JGP needed by DVJACS28
! from the DSM array NGRPDS.
! ..
!     Input:
!
!     N        = the order of the matrix
!     MAXGRPDS = number of groups (from DSM)
!     NGRPDS   = DSM output array
!     Output:
!     JGP      = array of length N containing the column
!                indices by groups
!     IGP      = pointer array of length NGRP + 1 to the
!                locations in JGP of the beginning of
!                each group
! ..
     implicit none
! ..
! .. Scalar Arguments ..
        integer, intent (IN) :: N, MAXGRPDS
! ..
! .. Array Arguments ..
        integer, intent (IN) :: NGRPDS(MAXGRPDS)
        integer, intent (OUT) :: IGP(MAXGRPDS+1), JGP(N)
! ..
! .. Local Scalars ..
        integer :: IGRP, INDEX, JCOL
! ..
! .. FIRST EXECUTABLE STATEMENT DGROUPDS
! ..
        IGP(1) = 1
        INDEX =  0
        do IGRP = 1, MAXGRPDS
           IGP(IGRP+1) = IGP(IGRP)
           do JCOL = 1, N
              if (NGRPDS(JCOL) == IGRP) then
                 IGP(IGRP+1) = IGP(IGRP+1) + 1
                 INDEX = INDEX + 1
                 JGP(INDEX) = JCOL
              end if
           end do
        end do
        return

      end subroutine DGROUPDS
!_______________________________________________________________________

    subroutine JACSPDB(FCN,N,T,Y,F,FJAC,NRFJAC,YSCALE,FAC,IOPT, &
      WK,LWK,IWK,LIWK,MAXGRP,NGRP,JPNTR,INDROW)

! This is a modified version of JACSP which does not require the NGRP,
! JPNTR and INDROW sparse pointer arrays in the event a dense or a
! banded matrix is being processed.

! Refer to the documentation for JACSP for a description of the
! parameters. If the banded option is used, IOPT(5) is used to
! input the lower bandwidth ML in this version.

     implicit none

! .. Parameters ..
      integer, parameter :: WP = kind(0.0D0)
! ..
! .. Scalar Arguments ..
      real (WP) :: T
      integer :: LIWK, LWK, MAXGRP, N, NRFJAC
! ..
! .. Array Arguments ..
      real (WP) :: F(N), FAC(N), FJAC(NRFJAC,*), WK(LWK), Y(N), YSCALE(*)
      integer :: INDROW(*), IOPT(5), IWK(LIWK), JPNTR(*), NGRP(N)
! ..
! .. Subroutine Arguments ..
      external FCN
! ..
! .. Local Scalars ..
      real (WP) :: ADIFF, AY, DEL, DELM, DFMJ, DIFF, DMAX, EXPFMN, FACMAX, &
        FACMIN, FJACL, FMJ, ONE, P125, P25, P75, P875, PERT, RDEL, RMNFDF, &
        RMXFDF, SDF, SF, SGN, T1, T2, U, U3QRT, U7EGT, UEGT, UMEGT, UQRT,  &
        USQT, ZERO
      integer :: IDXL, IDXU, IFLAG1, IFLAG2, IRCMP, IRDEL, IROW, IROWB,    &
        IROWMX, ITRY, J, JCOL, JFIRST, JINC, JLAST, KT1, KT2, KT3, KT4,    &
        KT5, L, MBAND, ML, MU, NID1, NID2, NID3, NID4, NID5, NID6,         &
        NIFAC, NT2, NUMGRP
!       KT5, L, MBAND, MEB1, ML, MU, NID1, NID2, NID3, NID4, NID5, NID6,   &
      logical :: DOTHISBLOCK
      character (80) :: MSG
! ..
! .. Intrinsic Functions ..
!     INTRINSIC ABS, EPSILON, KIND, MAX, MIN, SIGN, SQRT
      intrinsic ABS, EPSILON, MAX, MIN, SIGN, SQRT
! ..
! .. Data Statements ..
      data PERT/2.0E+0_WP/, FACMAX/1.E-1_WP/, EXPFMN/.75E+0_WP/
      data ONE/1.0E0_WP/, ZERO/0.0E0_WP/
      data P125/.125E+0_WP/, P25/.25E+0_WP/, P75/.75E+0_WP/, P875/.875E+0_WP/
      data NIFAC/3/, NID1/10/, NID2/20/, NID3/30/, NID4/40/, NID5/50/
! ..
!     COMPUTE ALGORITHM AND MACHINE CONSTANTS.
! ..
! ..  FIRST EXECUTABLE STATEMENT JACSPDB
! ..
      if (IOPT(1) == 0 .and. MAXGRP /= N) then
         MSG = 'JACSPDB requires that MAXGRP=N for a dense matrix.'
         call XERRDV(MSG,1800,2,0,0,0,0,ZERO,ZERO)
      end if

      if (IOPT(1) == 1) then
         MBAND = IOPT(2)
         ML = IOPT(5)
         MU = MBAND - ML - 1
!        MEB1 = 2*ML + MU
      end if

      U = epsilon(ONE)
      USQT = sqrt(U)
      UEGT = U**P125
      UMEGT = ONE/UEGT
      UQRT = U**P25
      U3QRT = U**P75
      U7EGT = U**P875
      FACMIN = U**EXPFMN

      if (IOPT(4) == 0) then
        IOPT(4) = 1
        do 10 J = 1, N
          FAC(J) = USQT
10      continue
      end if
      do 20 J = 1, 50
        IWK(J) = 0
20    continue
      KT1 = NID1
      KT2 = NID2
      KT3 = NID3
      KT4 = NID4
      KT5 = NID5
      NID6 = LIWK
      NT2 = 2*N

      do NUMGRP = 1, MAXGRP
!       COMPUTE AND SAVE THE INCREMENTS FOR THE COLUMNS IN GROUP NUMGRP.
        IRCMP = 0
        ITRY = 0
        do 30 J = NID5 + 1, NID6
          IWK(J) = 0
30      continue
40      continue

!       Note: For banded in DVJAC:
!       mba = min(mband,n)
!          Groups: j=1,...,mba
!             Columns in group j: jj=j,...,n by mband
!                Rows in column jj: i1=max(jj-mu,1),...,i2=min(jj+ml,n)       

        if (IOPT(1) == 0 .or. IOPT(1) == 2) then
          JFIRST = 1
          JLAST = N
          JINC = 1
        else
          JFIRST = NUMGRP
          JLAST = N
          JINC = MBAND
        end if
        do JCOL = JFIRST, JLAST, JINC
          DOTHISBLOCK = .false.
          if (IOPT(1) == 2) then
            if (NGRP(JCOL) == NUMGRP) DOTHISBLOCK = .true.
          else
            if (IOPT(1) == 0) then
              if (JCOL == NUMGRP) DOTHISBLOCK = .true.
            else
              if (IOPT(1) == 1) then
                DOTHISBLOCK = .true.
              end if
            end if
          end if

          if (DOTHISBLOCK) then
            WK(N+JCOL) = Y(JCOL)
!           COMPUTE DEL. IF DEL IS TOO SMALL INCREASE FAC(JCOL) AND RECOMPUTE
!           DEL. NIFAC ATTEMPTS ARE MADE TO INCREASE FAC(JCOL) AND FIND AN
!           APPROPRIATE DEL. IF DEL CANT BE FOUND IN THIS MANNER, DEL IS COMPUTED
!           WITH FAC(JCOL) SET TO THE SQUARE ROOT OF THE MACHINE PRECISION (USQT).
!           IF DEL IS ZERO TO MACHINE PRECISION BECAUSE Y(JCOL) IS ZERO OR
!           YSCALE(JCOL) IS ZERO, DEL IS SET TO USQT.
            SGN = sign(ONE,F(JCOL))
            IRDEL = 0
            if (IOPT(3) == 1) then
              AY = abs(YSCALE(JCOL))
            else
              AY = abs(Y(JCOL))
            end if
            DELM = U7EGT*AY

            do 50 J = 1, NIFAC
              DEL = FAC(JCOL)*AY*SGN
              if (abs(DEL) <= ZERO) then
                DEL = USQT*SGN
                if (ITRY == 0) IWK(3) = IWK(3) + 1
              end if
              T1 = Y(JCOL) + DEL
              DEL = T1 - Y(JCOL)
              if (abs(DEL) < DELM) then
                if (J >= NIFAC) goto 50
                if (IRDEL == 0) then
                  IRDEL = 1
                  IWK(1) = IWK(1) + 1
                end if
                T1 = FAC(JCOL)*UMEGT
                FAC(JCOL) = min(T1,FACMAX)
              else
                goto 60
              end if
50          end do

            FAC(JCOL) = USQT
            DEL = USQT*AY*SGN
            IWK(2) = IWK(2) + 1
            if (KT2 < NID3) then
              KT2 = KT2 + 1
              IWK(KT2) = JCOL
            end if

60          continue
            WK(NT2+JCOL) = DEL
            Y(JCOL) = Y(JCOL) + DEL
          end if
        end do

        IWK(7) = IWK(7) + 1
        call FCN(N,T,Y,WK)

        do JCOL = JFIRST, JLAST, JINC
          if (IOPT(1) == 2) then
             if (NGRP(JCOL) == NUMGRP) Y(JCOL) = WK(N+JCOL)
          else
             if (IOPT(1) == 0) then
                if (JCOL == NUMGRP) Y(JCOL) = WK(N+JCOL)
             else
                if (IOPT(1) == 1) Y(JCOL) = WK(N+JCOL) 
             end if
          end if
        end do

!       COMPUTE THE JACOBIAN ENTRIES FOR ALL COLUMNS IN NUMGRP.
!       STORE ENTRIES ACCORDING TO SELECTED STORAGE FORMAT.
!       USE LARGEST ELEMENTS IN A COLUMN TO DETERMINE SCALING
!       INFORMATION FOR ROUNDOFF AND TRUNCATION ERROR TESTS.

        do JCOL = JFIRST, JLAST, JINC
          DOTHISBLOCK = .false.
          if (IOPT(1) == 2) then
             if (NGRP(JCOL) == NUMGRP) DOTHISBLOCK = .true.
          else
             if (IOPT(1) == 0) then
                if (JCOL == NUMGRP) DOTHISBLOCK = .true.
             else
                if(IOPT(1) == 1) DOTHISBLOCK = .true. 
             end if
          end if
          if (DOTHISBLOCK) then
            if (IOPT(1) == 2) then
              IDXL = JPNTR(JCOL)
              IDXU = JPNTR(JCOL+1) - 1
            else
              if (IOPT(1) == 0) then
                IDXL = 1
                IDXU = N
              else
                if (IOPT(1) == 1) then
                  IDXL = max(JCOL-MU,1)
                  IDXU = min(JCOL+ML,N)
                end if
              end if
            end if
            DMAX = ZERO
            RDEL = ONE/WK(NT2+JCOL)
            IROWMX = 1
            do L = IDXL, IDXU
              if (IOPT(1) == 2) then
                 IROW = INDROW(L)
              else
                 if (IOPT(1)==0 .or. IOPT(1)==1) IROW = L
              end if
              DIFF = WK(IROW) - F(IROW)
              ADIFF = abs(DIFF)
              if (ADIFF >= DMAX) then
                IROWMX = IROW
                DMAX = ADIFF
                SF = F(IROW)
                SDF = WK(IROW)
              end if
              FJACL = DIFF*RDEL
              if (ITRY == 1) WK(IROW) = FJACL

              if (IOPT(1) == 0) then
                if (ITRY == 1) WK(IROW+N) = FJAC(IROW,JCOL)
                FJAC(IROW,JCOL) = FJACL
              end if
              if (IOPT(1) == 1) then
                IROWB = IROW - JCOL + IOPT(2)
                if (ITRY == 1) WK(IROW+N) = FJAC(IROWB,JCOL)
                FJAC(IROWB,JCOL) = FJACL
!               IROWB = JCOL * MEB1 - ML + L
!               IF (ITRY == 1) WK(IROW+N) = FJAC(IROWB,1)
!               FJAC(IROWB,1) = FJACL
              end if
              if (IOPT(1) == 2) then
                if (ITRY == 1) WK(IROW+N) = FJAC(L,1)
                FJAC(L,1) = FJACL
              end if
            end do

!           IF A COLUMN IS BEING RECOMPUTED (ITRY=1),THIS SECTION OF THE
!           CODE PERFORMS AN EXTRAPOLATION TEST TO ENABLE THE CODE TO
!           COMPUTE SMALL DERIVATIVES MORE ACCURATELY. THIS TEST IS ONLY
!           PERFORMED ON THOSE COLUMNS WHOSE LARGEST DIFFERENCE IS CLOSE
!           TO ZERO RELATIVE TO SCALE.
            if (ITRY == 1) then
              IFLAG1 = 0
              IFLAG2 = 0
              do 100 J = NID5 + 1, NID6
                if (IWK(J) == JCOL) IFLAG1 = 1
100           continue
              if (IFLAG1 == 1) then
                IFLAG1 = 0
                T1 = WK(IROWMX+N)
                T2 = WK(IROWMX)*FAC(JCOL)
                if (abs(T2) < abs(T1)*PERT) IFLAG2 = 1
              end if

              if (IFLAG2 == 1) then
                IFLAG2 = 0
                T1 = FAC(JCOL)*FAC(JCOL)
                FAC(JCOL) = max(T1,FACMIN)
                do L = IDXL, IDXU
                  if (IOPT(1) == 2) then
                     IROW = INDROW(L)
                  else
                     if (IOPT(1)==0 .or. IOPT(1)==1) IROW = L
                  end if
                  FJACL = WK(IROW+N)
                  if (IOPT(1) == 0) FJAC(IROW,JCOL) = FJACL
                  if (IOPT(1) == 1) then
                    IROWB = IROW - JCOL + IOPT(2)
                    FJAC(IROWB,JCOL) = FJACL
!                   IROWB = JCOL * MEB1 - ML + L
!                   FJAC(IROWB,1) = FJACL
                  end if
                  if (IOPT(1) == 2) FJAC(L,1) = FJACL
              end do
              end if
            end if

            FMJ = abs(SF)
            DFMJ = abs(SDF)
            RMXFDF = max(FMJ,DFMJ)
            RMNFDF = min(FMJ,DFMJ)

!           IF SCALE INFORMATION IS NOT AVAILABLE, PERFORM NO ROUNDOFF
!           OR TRUNCATION ERROR TESTS. IF THE EXTRAPOLATION TEST HAS
!           CAUSED FAC(JCOL) TO BE RESET TO ITS PREVIOUS VALUE (IAJAC=1)
!           THEN NO FURTHER ROUNDOFF OR TRUNCATION ERROR TESTS ARE
!           PERFORMED.
!           IF (RMNFDF/=ZERO) THEN
            if (abs(RMNFDF) > ZERO) then
!             TEST FOR POSSIBLE ROUNDOFF ERROR (FIRST TEST)
!             AND ALSO FOR POSSIBLE SERIOUS ROUNDOFF ERROR (SECOND TEST).
              if (DMAX <= (U3QRT*RMXFDF)) then
                if (DMAX <= (U7EGT*RMXFDF)) then
                  if (ITRY == 0) then
                    T1 = sqrt(FAC(JCOL))
                    FAC(JCOL) = min(T1,FACMAX)
                    IRCMP = 1
                    if (KT5 < NID6) then
                      KT5 = KT5 + 1
                      IWK(KT5) = JCOL
                    end if
                    IWK(4) = IWK(4) + 1
                    if (KT3 < NID4) then
                      KT3 = KT3 + 1
                      IWK(KT3) = JCOL
                    end if
                  else
                    IWK(5) = IWK(5) + 1
                    if (KT4 < NID5) then
                      KT4 = KT4 + 1
                      IWK(KT4) = JCOL
                    end if
                  end if
                else
                  T1 = UMEGT*FAC(JCOL)
                  FAC(JCOL) = min(T1,FACMAX)
                end if
              end if
!             TEST FOR POSSIBLE TRUNCATION ERROR.
              if (DMAX > UQRT*RMXFDF) then
                T1 = FAC(JCOL)*UEGT
                FAC(JCOL) = max(T1,FACMIN)
                IWK(8) = IWK(8) + 1
                if (KT1 < NID2) then
                  KT1 = KT1 + 1
                  IWK(KT1) = JCOL
                end if
              end if
            else
              IWK(6) = IWK(6) + 1
            end if
          end if
        end do

!       IF SERIOUS ROUNDOFF ERROR IS SUSPECTED, RECOMPUTE ALL
!       COLUMNS IN GROUP NUMGRP.
        if (IRCMP == 1) then
          IRCMP = 0
          ITRY = 1
          goto 40
        end if
        ITRY = 0
      end do
      return

    end subroutine JACSPDB

! End of JACSP routines.
!_______________________________________________________________________
! *****MA48 build change point. Insert the MA48 Jacobian related
! routines DVNLSS48, DVSOLS48, DVPREPS48, and DVJACS48 here.
! filename = jacobian_for_ma48.f90. Insert the following line
! after the first line of this file:
! USE hsl_ma48_double
!_______________________________________________________________________

! *****LAPACK build change point. Insert the following line after the
! first line of this file:
! USE lapackd_f90_m
! Include the module lapackd_f90_m.f90 at the beginning of this file
! or as an external module.
!_______________________________________________________________________

    end module DVODE_F90_M
!_______________________________________________________________________
