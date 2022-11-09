use libc::{c_double, c_int};

#[link(name = "gfortran")]
extern "C" {
/// Call `DLSODE` subroutine from ODEPACK
///
///***PURPOSE  Livermore Solver for Ordinary Differential Equations.
///            DLSODE solves the initial-value problem for stiff or
///            nonstiff systems of first-order ODE's,
///               dy/dt = f(t,y),   or, in component form,
///               dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(N)),  i=1,...,N.
///***CATEGORY  I1A
///***TYPE      DOUBLE PRECISION (SLSODE-S, DLSODE-D)
///***KEYWORDS  ORDINARY DIFFERENTIAL EQUATIONS, INITIAL VALUE PROBLEM,
///             STIFF, NONSTIFF
///***AUTHOR  Hindmarsh, Alan C., (LLNL)
///             Center for Applied Scientific Computing, L-561
///             Lawrence Livermore National Laboratory
///             Livermore, CA 94551.
///***DESCRIPTION
///
///     NOTE: The "Usage" and "Arguments" sections treat only a subset of
///           available options, in condensed fashion.  The options
///           covered and the information supplied will support most
///           standard uses of DLSODE.
///
///           For more sophisticated uses, full details on all options are
///           given in the concluding section, headed "Long Description."
///           A synopsis of the DLSODE Long Description is provided at the
///           beginning of that section; general topics covered are:
///           - Elements of the call sequence; optional input and output
///           - Optional supplemental routines in the DLSODE package
///           - internal COMMON block
///
/// *Usage:
///     Communication between the user and the DLSODE package, for normal
///     situations, is summarized here.  This summary describes a subset
///     of the available options.  See "Long Description" for complete
///     details, including optional communication, nonstandard options,
///     and instructions for special situations.
///
///     A sample program is given in the "Examples" section.
///
///     Refer to the argument descriptions for the definitions of the
///     quantities that appear in the following sample declarations.
///
///     For MF = 10,
///        PARAMETER  (LRW = 20 + 16*NEQ,           LIW = 20)
///     For MF = 21 or 22,
///        PARAMETER  (LRW = 22 +  9*NEQ + NEQ**2,  LIW = 20 + NEQ)
///     For MF = 24 or 25,
///        PARAMETER  (LRW = 22 + 10*NEQ + (2*ML+MU)*NEQ,
///       *                                         LIW = 20 + NEQ)
///
///        EXTERNAL F, JAC
///        INTEGER  NEQ, ITOL, ITASK, ISTATE, IOPT, LRW, IWORK(LIW),
///       *         LIW, MF
///        DOUBLE PRECISION Y(NEQ), T, TOUT, RTOL, ATOL(ntol), RWORK(LRW)
///
///        CALL DLSODE (F, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
///       *            ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF)
///
/// *Arguments:
///     F     :EXT    Name of subroutine for right-hand-side vector f.
///                   This name must be declared EXTERNAL in calling
///                   program.  The form of F must be:
///
///                   SUBROUTINE  F (NEQ, T, Y, YDOT)
///                   INTEGER  NEQ
///                   DOUBLE PRECISION  T, Y(*), YDOT(*)
///
///                   The inputs are NEQ, T, Y.  F is to set
///
///                   YDOT(i) = f(i,T,Y(1),Y(2),...,Y(NEQ)),
///                                                     i = 1, ..., NEQ .
///
///     NEQ   :IN     Number of first-order ODE's.
///
///     Y     :INOUT  Array of values of the y(t) vector, of length NEQ.
///                   Input:  For the first call, Y should contain the
///                           values of y(t) at t = T. (Y is an input
///                           variable only if ISTATE = 1.)
///                   Output: On return, Y will contain the values at the
///                           new t-value.
///
///     T     :INOUT  Value of the independent variable.  On return it
///                   will be the current value of t (normally TOUT).
///
///     TOUT  :IN     Next point where output is desired (.NE. T).
///
///     ITOL  :IN     1 or 2 according as ATOL (below) is a scalar or
///                   an array.
///
///     RTOL  :IN     Relative tolerance parameter (scalar).
///
///     ATOL  :IN     Absolute tolerance parameter (scalar or array).
///                   If ITOL = 1, ATOL need not be dimensioned.
///                   If ITOL = 2, ATOL must be dimensioned at least NEQ.
///
///                   The estimated local error in Y(i) will be controlled
///                   so as to be roughly less (in magnitude) than
///
///                   EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
///                   EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
///
///                   Thus the local error test passes if, in each
///                   component, either the absolute error is less than
///                   ATOL (or ATOL(i)), or the relative error is less
///                   than RTOL.
///
///                   Use RTOL = 0.0 for pure absolute error control, and
///                   use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative
///                   error control.  Caution:  Actual (global) errors may
///                   exceed these local tolerances, so choose them
///                   conservatively.
///
///     ITASK :IN     Flag indicating the task DLSODE is to perform.
///                   Use ITASK = 1 for normal computation of output
///                   values of y at t = TOUT.
///
///     ISTATE:INOUT  Index used for input and output to specify the state
///                   of the calculation.
///                   Input:
///                    1   This is the first call for a problem.
///                    2   This is a subsequent call.
///                   Output:
///                    1   Nothing was done, because TOUT was equal to T.
///                    2   DLSODE was successful (otherwise, negative).
///                        Note that ISTATE need not be modified after a
///                        successful return.
///                   -1   Excess work done on this call (perhaps wrong
///                        MF).
///                   -2   Excess accuracy requested (tolerances too
///                        small).
///                   -3   Illegal input detected (see printed message).
///                   -4   Repeated error test failures (check all
///                        inputs).
///                   -5   Repeated convergence failures (perhaps bad
///                        Jacobian supplied or wrong choice of MF or
///                        tolerances).
///                   -6   Error weight became zero during problem
///                        (solution component i vanished, and ATOL or
///                        ATOL(i) = 0.).
///
///     IOPT  :IN     Flag indicating whether optional inputs are used:
///                   0   No.
///                   1   Yes.  (See "Optional inputs" under "Long
///                       Description," Part 1.)
///
///     RWORK :WORK   Real work array of length at least:
///                   20 + 16*NEQ                    for MF = 10,
///                   22 +  9*NEQ + NEQ**2           for MF = 21 or 22,
///                   22 + 10*NEQ + (2*ML + MU)*NEQ  for MF = 24 or 25.
///
///     LRW   :IN     Declared length of RWORK (in user's DIMENSION
///                   statement).
///
///     IWORK :WORK   Integer work array of length at least:
///                   20        for MF = 10,
///                   20 + NEQ  for MF = 21, 22, 24, or 25.
///
///                   If MF = 24 or 25, input in IWORK(1),IWORK(2) the
///                   lower and upper Jacobian half-bandwidths ML,MU.
///
///                   On return, IWORK contains information that may be
///                   of interest to the user:
///
///            Name   Location   Meaning
///            -----  ---------  -----------------------------------------
///            NST    IWORK(11)  Number of steps taken for the problem so
///                              far.
///            NFE    IWORK(12)  Number of f evaluations for the problem
///                              so far.
///            NJE    IWORK(13)  Number of Jacobian evaluations (and of
///                              matrix LU decompositions) for the problem
///                              so far.
///            NQU    IWORK(14)  Method order last used (successfully).
///            LENRW  IWORK(17)  Length of RWORK actually required.  This
///                              is defined on normal returns and on an
///                              illegal input return for insufficient
///                              storage.
///            LENIW  IWORK(18)  Length of IWORK actually required.  This
///                              is defined on normal returns and on an
///                              illegal input return for insufficient
///                              storage.
///
///     LIW   :IN     Declared length of IWORK (in user's DIMENSION
///                   statement).
///
///     JAC   :EXT    Name of subroutine for Jacobian matrix (MF =
///                   21 or 24).  If used, this name must be declared
///                   EXTERNAL in calling program.  If not used, pass a
///                   dummy name.  The form of JAC must be:
///
///                   SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
///                   INTEGER  NEQ, ML, MU, NROWPD
///                   DOUBLE PRECISION  T, Y(*), PD(NROWPD,*)
///
///                   See item c, under "Description" below for more
///                   information about JAC.
///
///     MF    :IN     Method flag.  Standard values are:
///                   10  Nonstiff (Adams) method, no Jacobian used.
///                   21  Stiff (BDF) method, user-supplied full Jacobian.
///                   22  Stiff method, internally generated full
///                       Jacobian.
///                   24  Stiff method, user-supplied banded Jacobian.
///                   25  Stiff method, internally generated banded
///                       Jacobian.
///
/// *Description:
///     DLSODE solves the initial value problem for stiff or nonstiff
///     systems of first-order ODE's,
///
///        dy/dt = f(t,y) ,
///
///     or, in component form,
///
///        dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ))
///                                                  (i = 1, ..., NEQ) .
///
///     DLSODE is a package based on the GEAR and GEARB packages, and on
///     the October 23, 1978, version of the tentative ODEPACK user
///     interface standard, with minor modifications.
///
///     The steps in solving such a problem are as follows.
///
///     a. First write a subroutine of the form
///
///           SUBROUTINE  F (NEQ, T, Y, YDOT)
///           INTEGER  NEQ
///           DOUBLE PRECISION  T, Y(*), YDOT(*)
///
///        which supplies the vector function f by loading YDOT(i) with
///        f(i).
///
///     b. Next determine (or guess) whether or not the problem is stiff.
///        Stiffness occurs when the Jacobian matrix df/dy has an
///        eigenvalue whose real part is negative and large in magnitude
///        compared to the reciprocal of the t span of interest.  If the
///        problem is nonstiff, use method flag MF = 10.  If it is stiff,
///        there are four standard choices for MF, and DLSODE requires the
///        Jacobian matrix in some form.  This matrix is regarded either
///        as full (MF = 21 or 22), or banded (MF = 24 or 25).  In the
///        banded case, DLSODE requires two half-bandwidth parameters ML
///        and MU. These are, respectively, the widths of the lower and
///        upper parts of the band, excluding the main diagonal.  Thus the
///        band consists of the locations (i,j) with
///
///           i - ML <= j <= i + MU ,
///
///        and the full bandwidth is ML + MU + 1 .
///
///     c. If the problem is stiff, you are encouraged to supply the
///        Jacobian directly (MF = 21 or 24), but if this is not feasible,
///        DLSODE will compute it internally by difference quotients (MF =
///        22 or 25).  If you are supplying the Jacobian, write a
///        subroutine of the form
///
///           SUBROUTINE  JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
///           INTEGER  NEQ, ML, MU, NRWOPD
///           DOUBLE PRECISION  T, Y(*), PD(NROWPD,*)
///
///        which provides df/dy by loading PD as follows:
///        - For a full Jacobian (MF = 21), load PD(i,j) with df(i)/dy(j),
///          the partial derivative of f(i) with respect to y(j).  (Ignore
///          the ML and MU arguments in this case.)
///        - For a banded Jacobian (MF = 24), load PD(i-j+MU+1,j) with
///          df(i)/dy(j); i.e., load the diagonal lines of df/dy into the
///          rows of PD from the top down.
///        - In either case, only nonzero elements need be loaded.
///
///     d. Write a main program that calls subroutine DLSODE once for each
///        point at which answers are desired.  This should also provide
///        for possible use of logical unit 6 for output of error messages
///        by DLSODE.
///
///        Before the first call to DLSODE, set ISTATE = 1, set Y and T to
///        the initial values, and set TOUT to the first output point.  To
///        continue the integration after a successful return, simply
///        reset TOUT and call DLSODE again.  No other parameters need be
///        reset.
///
/// *Examples:
///     The following is a simple example problem, with the coding needed
///     for its solution by DLSODE. The problem is from chemical kinetics,
///     and consists of the following three rate equations:
///
///        dy1/dt = -.04*y1 + 1.E4*y2*y3
///        dy2/dt = .04*y1 - 1.E4*y2*y3 - 3.E7*y2**2
///        dy3/dt = 3.E7*y2**2
///
///     on the interval from t = 0.0 to t = 4.E10, with initial conditions
///     y1 = 1.0, y2 = y3 = 0. The problem is stiff.
///
///     The following coding solves this problem with DLSODE, using 
///     MF = 21 and printing results at t = .4, 4., ..., 4.E10.  It uses 
///     ITOL = 2 and ATOL much smaller for y2 than for y1 or y3 because y2 
///     has much smaller values.  At the end of the run, statistical 
///     quantities of interest are printed.
///
///        EXTERNAL  FEX, JEX
///        INTEGER  IOPT, IOUT, ISTATE, ITASK, ITOL, IWORK(23), LIW, LRW,
///       *         MF, NEQ
///        DOUBLE PRECISION  ATOL(3), RTOL, RWORK(58), T, TOUT, Y(3)
///        NEQ = 3
///        Y(1) = 1.D0
///        Y(2) = 0.D0
///        Y(3) = 0.D0
///        T = 0.D0
///        TOUT = .4D0
///        ITOL = 2
///        RTOL = 1.D-4
///        ATOL(1) = 1.D-6
///        ATOL(2) = 1.D-10
///        ATOL(3) = 1.D-6
///        ITASK = 1
///        ISTATE = 1
///        IOPT = 0
///        LRW = 58
///        LIW = 23
///        MF = 21
///        DO 40 IOUT = 1,12
///          CALL DLSODE (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
///       *               ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
///          WRITE(6,20)  T, Y(1), Y(2), Y(3)
///    20    FORMAT(' At t =',D12.4,'   y =',3D14.6)
///          IF (ISTATE .LT. 0)  GO TO 80
///    40    TOUT = TOUT*10.D0
///        WRITE(6,60)  IWORK(11), IWORK(12), IWORK(13)
///    60  FORMAT(/' No. steps =',i4,',  No. f-s =',i4,',  No. J-s =',i4)
///        STOP
///    80  WRITE(6,90)  ISTATE
///    90  FORMAT(///' Error halt.. ISTATE =',I3)
///        STOP
///        END
///
///        SUBROUTINE  FEX (NEQ, T, Y, YDOT)
///        INTEGER  NEQ
///        DOUBLE PRECISION  T, Y(3), YDOT(3)
///        YDOT(1) = -.04D0*Y(1) + 1.D4*Y(2)*Y(3)
///        YDOT(3) = 3.D7*Y(2)*Y(2)
///        YDOT(2) = -YDOT(1) - YDOT(3)
///        RETURN
///        END
///
///        SUBROUTINE  JEX (NEQ, T, Y, ML, MU, PD, NRPD)
///        INTEGER  NEQ, ML, MU, NRPD
///        DOUBLE PRECISION  T, Y(3), PD(NRPD,3)
///        PD(1,1) = -.04D0
///        PD(1,2) = 1.D4*Y(3)
///        PD(1,3) = 1.D4*Y(2)
///        PD(2,1) = .04D0
///        PD(2,3) = -PD(1,3)
///        PD(3,2) = 6.D7*Y(2)
///        PD(2,2) = -PD(1,2) - PD(3,2)
///        RETURN
///        END
///
///     The output from this program (on a Cray-1 in single precision)
///     is as follows.
///
///     At t =  4.0000e-01   y =  9.851726e-01  3.386406e-05  1.479357e-02
///     At t =  4.0000e+00   y =  9.055142e-01  2.240418e-05  9.446344e-02
///     At t =  4.0000e+01   y =  7.158050e-01  9.184616e-06  2.841858e-01
///     At t =  4.0000e+02   y =  4.504846e-01  3.222434e-06  5.495122e-01
///     At t =  4.0000e+03   y =  1.831701e-01  8.940379e-07  8.168290e-01
///     At t =  4.0000e+04   y =  3.897016e-02  1.621193e-07  9.610297e-01
///     At t =  4.0000e+05   y =  4.935213e-03  1.983756e-08  9.950648e-01
///     At t =  4.0000e+06   y =  5.159269e-04  2.064759e-09  9.994841e-01
///     At t =  4.0000e+07   y =  5.306413e-05  2.122677e-10  9.999469e-01
///     At t =  4.0000e+08   y =  5.494530e-06  2.197825e-11  9.999945e-01
///     At t =  4.0000e+09   y =  5.129458e-07  2.051784e-12  9.999995e-01
///     At t =  4.0000e+10   y = -7.170603e-08 -2.868241e-13  1.000000e+00
///
///     No. steps = 330,  No. f-s = 405,  No. J-s = 69
///
/// *Accuracy:
///     The accuracy of the solution depends on the choice of tolerances
///     RTOL and ATOL.  Actual (global) errors may exceed these local
///     tolerances, so choose them conservatively.
///
/// *Cautions:
///     The work arrays should not be altered between calls to DLSODE for
///     the same problem, except possibly for the conditional and optional
///     inputs.
///
/// *Portability:
///     Since NEQ is dimensioned inside DLSODE, some compilers may object
///     to a call to DLSODE with NEQ a scalar variable.  In this event, 
///     use DIMENSION NEQ(1).  Similar remarks apply to RTOL and ATOL.
///
///     Note to Cray users:
///     For maximum efficiency, use the CFT77 compiler.  Appropriate
///     compiler optimization directives have been inserted for CFT77.
///
/// *Reference:
///     Alan C. Hindmarsh, "ODEPACK, A Systematized Collection of ODE
///     Solvers," in Scientific Computing, R. S. Stepleman, et al., Eds.
///     (North-Holland, Amsterdam, 1983), pp. 55-64.
///
/// *Long Description:
///     The following complete description of the user interface to
///     DLSODE consists of four parts:
///
///     1.  The call sequence to subroutine DLSODE, which is a driver
///         routine for the solver.  This includes descriptions of both
///         the call sequence arguments and user-supplied routines.
///         Following these descriptions is a description of optional
///         inputs available through the call sequence, and then a
///         description of optional outputs in the work arrays.
///
///     2.  Descriptions of other routines in the DLSODE package that may
///         be (optionally) called by the user.  These provide the ability
///         to alter error message handling, save and restore the internal
///         COMMON, and obtain specified derivatives of the solution y(t).
///
///     3.  Descriptions of COMMON block to be declared in overlay or
///         similar environments, or to be saved when doing an interrupt
///         of the problem and continued solution later.
///
///     4.  Description of two routines in the DLSODE package, either of
///         which the user may replace with his own version, if desired.
///         These relate to the measurement of errors.
///
///
///                         Part 1.  Call Sequence
///                         ----------------------
///
///     Arguments
///     ---------
///     The call sequence parameters used for input only are
///
///        F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF,
///
///     and those used for both input and output are
///
///        Y, T, ISTATE.
///
///     The work arrays RWORK and IWORK are also used for conditional and
///     optional inputs and optional outputs.  (The term output here
///     refers to the return from subroutine DLSODE to the user's calling
///     program.)
///
///     The legality of input parameters will be thoroughly checked on the
///     initial call for the problem, but not checked thereafter unless a
///     change in input parameters is flagged by ISTATE = 3 on input.
///
///     The descriptions of the call arguments are as follows.
///
///     F        The name of the user-supplied subroutine defining the ODE
///              system.  The system must be put in the first-order form
///              dy/dt = f(t,y), where f is a vector-valued function of
///              the scalar t and the vector y. Subroutine F is to compute
///              the function f. It is to have the form
///
///                 SUBROUTINE F (NEQ, T, Y, YDOT)
///                 DOUBLE PRECISION  T, Y(*), YDOT(*)
///
///              where NEQ, T, and Y are input, and the array YDOT =
///              f(T,Y) is output.  Y and YDOT are arrays of length NEQ.
///              Subroutine F should not alter Y(1),...,Y(NEQ).  F must be
///              declared EXTERNAL in the calling program.
///
///              Subroutine F may access user-defined quantities in
///              NEQ(2),... and/or in Y(NEQ(1)+1),..., if NEQ is an array
///              (dimensioned in F) and/or Y has length exceeding NEQ(1).
///              See the descriptions of NEQ and Y below.
///
///              If quantities computed in the F routine are needed
///              externally to DLSODE, an extra call to F should be made
///              for this purpose, for consistent and accurate results.
///              If only the derivative dy/dt is needed, use DINTDY
///              instead.
///
///     NEQ      The size of the ODE system (number of first-order
///              ordinary differential equations).  Used only for input.
///              NEQ may be decreased, but not increased, during the
///              problem.  If NEQ is decreased (with ISTATE = 3 on input),
///              the remaining components of Y should be left undisturbed,
///              if these are to be accessed in F and/or JAC.
///
///              Normally, NEQ is a scalar, and it is generally referred
///              to as a scalar in this user interface description.
///              However, NEQ may be an array, with NEQ(1) set to the
///              system size.  (The DLSODE package accesses only NEQ(1).)
///              In either case, this parameter is passed as the NEQ
///              argument in all calls to F and JAC.  Hence, if it is an
///              array, locations NEQ(2),... may be used to store other
///              integer data and pass it to F and/or JAC.  Subroutines
///              F and/or JAC must include NEQ in a DIMENSION statement
///              in that case.
///
///     Y        A real array for the vector of dependent variables, of
///              length NEQ or more.  Used for both input and output on
///              the first call (ISTATE = 1), and only for output on
///              other calls.  On the first call, Y must contain the
///              vector of initial values.  On output, Y contains the
///              computed solution vector, evaluated at T. If desired,
///              the Y array may be used for other purposes between
///              calls to the solver.
///
///              This array is passed as the Y argument in all calls to F
///              and JAC.  Hence its length may exceed NEQ, and locations
///              Y(NEQ+1),... may be used to store other real data and
///              pass it to F and/or JAC.  (The DLSODE package accesses
///              only Y(1),...,Y(NEQ).)
///
///     T        The independent variable.  On input, T is used only on
///              the first call, as the initial point of the integration.
///              On output, after each call, T is the value at which a
///              computed solution Y is evaluated (usually the same as
///              TOUT).  On an error return, T is the farthest point
///              reached.
///
///     TOUT     The next value of T at which a computed solution is
///              desired.  Used only for input.
///
///              When starting the problem (ISTATE = 1), TOUT may be equal
///              to T for one call, then should not equal T for the next
///              call.  For the initial T, an input value of TOUT .NE. T
///              is used in order to determine the direction of the
///              integration (i.e., the algebraic sign of the step sizes)
///              and the rough scale of the problem.  Integration in
///              either direction (forward or backward in T) is permitted.
///
///              If ITASK = 2 or 5 (one-step modes), TOUT is ignored
///              after the first call (i.e., the first call with
///              TOUT .NE. T).  Otherwise, TOUT is required on every call.
///
///              If ITASK = 1, 3, or 4, the values of TOUT need not be
///              monotone, but a value of TOUT which backs up is limited
///              to the current internal T interval, whose endpoints are
///              TCUR - HU and TCUR.  (See "Optional Outputs" below for
///              TCUR and HU.)
///
///
///     ITOL     An indicator for the type of error control.  See
///              description below under ATOL.  Used only for input.
///
///     RTOL     A relative error tolerance parameter, either a scalar or
///              an array of length NEQ.  See description below under
///              ATOL.  Input only.
///
///     ATOL     An absolute error tolerance parameter, either a scalar or
///              an array of length NEQ.  Input only.
///
///              The input parameters ITOL, RTOL, and ATOL determine the
///              error control performed by the solver.  The solver will
///              control the vector e = (e(i)) of estimated local errors
///              in Y, according to an inequality of the form
///
///                 rms-norm of ( e(i)/EWT(i) ) <= 1,
///
///              where
///
///                 EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
///
///              and the rms-norm (root-mean-square norm) here is
///
///                 rms-norm(v) = SQRT(sum v(i)**2 / NEQ).
///
///              Here EWT = (EWT(i)) is a vector of weights which must
///              always be positive, and the values of RTOL and ATOL
///              should all be nonnegative.  The following table gives the
///              types (scalar/array) of RTOL and ATOL, and the
///              corresponding form of EWT(i).
///
///              ITOL    RTOL      ATOL      EWT(i)
///              ----    ------    ------    -----------------------------
///              1       scalar    scalar    RTOL*ABS(Y(i)) + ATOL
///              2       scalar    array     RTOL*ABS(Y(i)) + ATOL(i)
///              3       array     scalar    RTOL(i)*ABS(Y(i)) + ATOL
///              4       array     array     RTOL(i)*ABS(Y(i)) + ATOL(i)
///
///              When either of these parameters is a scalar, it need not
///              be dimensioned in the user's calling program.
///
///              If none of the above choices (with ITOL, RTOL, and ATOL
///              fixed throughout the problem) is suitable, more general
///              error controls can be obtained by substituting
///              user-supplied routines for the setting of EWT and/or for
///              the norm calculation.  See Part 4 below.
///
///              If global errors are to be estimated by making a repeated
///              run on the same problem with smaller tolerances, then all
///              components of RTOL and ATOL (i.e., of EWT) should be
///              scaled down uniformly.
///
///     ITASK    An index specifying the task to be performed.  Input
///              only.  ITASK has the following values and meanings:
///              1   Normal computation of output values of y(t) at
///                  t = TOUT (by overshooting and interpolating).
///              2   Take one step only and return.
///              3   Stop at the first internal mesh point at or beyond
///                  t = TOUT and return.
///              4   Normal computation of output values of y(t) at
///                  t = TOUT but without overshooting t = TCRIT.  TCRIT
///                  must be input as RWORK(1).  TCRIT may be equal to or
///                  beyond TOUT, but not behind it in the direction of
///                  integration.  This option is useful if the problem
///                  has a singularity at or beyond t = TCRIT.
///              5   Take one step, without passing TCRIT, and return.
///                  TCRIT must be input as RWORK(1).
///
///              Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
///              (within roundoff), it will return T = TCRIT (exactly) to
///              indicate this (unless ITASK = 4 and TOUT comes before
///              TCRIT, in which case answers at T = TOUT are returned
///              first).
///
///     ISTATE   An index used for input and output to specify the state
///              of the calculation.
///
///              On input, the values of ISTATE are as follows:
///              1   This is the first call for the problem
///                  (initializations will be done).  See "Note" below.
///              2   This is not the first call, and the calculation is to
///                  continue normally, with no change in any input
///                  parameters except possibly TOUT and ITASK.  (If ITOL,
///                  RTOL, and/or ATOL are changed between calls with
///                  ISTATE = 2, the new values will be used but not
///                  tested for legality.)
///              3   This is not the first call, and the calculation is to
///                  continue normally, but with a change in input
///                  parameters other than TOUT and ITASK.  Changes are
///                  allowed in NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF,
///                  ML, MU, and any of the optional inputs except H0.
///                  (See IWORK description for ML and MU.)
///
///              Note:  A preliminary call with TOUT = T is not counted as
///              a first call here, as no initialization or checking of
///              input is done.  (Such a call is sometimes useful for the
///              purpose of outputting the initial conditions.)  Thus the
///              first call for which TOUT .NE. T requires ISTATE = 1 on
///              input.
///
///              On output, ISTATE has the following values and meanings:
///               1  Nothing was done, as TOUT was equal to T with
///                  ISTATE = 1 on input.
///               2  The integration was performed successfully.
///              -1  An excessive amount of work (more than MXSTEP steps)
///                  was done on this call, before completing the
///                  requested task, but the integration was otherwise
///                  successful as far as T. (MXSTEP is an optional input
///                  and is normally 500.)  To continue, the user may
///                  simply reset ISTATE to a value >1 and call again (the
///                  excess work step counter will be reset to 0).  In
///                  addition, the user may increase MXSTEP to avoid this
///                  error return; see "Optional Inputs" below.
///              -2  Too much accuracy was requested for the precision of
///                  the machine being used.  This was detected before
///                  completing the requested task, but the integration
///                  was successful as far as T. To continue, the
///                  tolerance parameters must be reset, and ISTATE must
///                  be set to 3. The optional output TOLSF may be used
///                  for this purpose.  (Note:  If this condition is
///                  detected before taking any steps, then an illegal
///                  input return (ISTATE = -3) occurs instead.)
///              -3  Illegal input was detected, before taking any
///                  integration steps.  See written message for details.
///                  (Note:  If the solver detects an infinite loop of
///                  calls to the solver with illegal input, it will cause
///                  the run to stop.)
///              -4  There were repeated error-test failures on one
///                  attempted step, before completing the requested task,
///                  but the integration was successful as far as T.  The
///                  problem may have a singularity, or the input may be
///                  inappropriate.
///              -5  There were repeated convergence-test failures on one
///                  attempted step, before completing the requested task,
///                  but the integration was successful as far as T. This
///                  may be caused by an inaccurate Jacobian matrix, if
///                  one is being used.
///              -6  EWT(i) became zero for some i during the integration.
///                  Pure relative error control (ATOL(i)=0.0) was
///                  requested on a variable which has now vanished.  The
///                  integration was successful as far as T.
///
///              Note:  Since the normal output value of ISTATE is 2, it
///              does not need to be reset for normal continuation.  Also,
///              since a negative input value of ISTATE will be regarded
///              as illegal, a negative output value requires the user to
///              change it, and possibly other inputs, before calling the
///              solver again.
///
///     IOPT     An integer flag to specify whether any optional inputs
///              are being used on this call.  Input only.  The optional
///              inputs are listed under a separate heading below.
///              0   No optional inputs are being used.  Default values
///                  will be used in all cases.
///              1   One or more optional inputs are being used.
///
///     RWORK    A real working array (double precision).  The length of
///              RWORK must be at least
///
///                 20 + NYH*(MAXORD + 1) + 3*NEQ + LWM
///
///              where
///                 NYH = the initial value of NEQ,
///              MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
///                       smaller value is given as an optional input),
///                 LWM = 0           if MITER = 0,
///                 LWM = NEQ**2 + 2  if MITER = 1 or 2,
///                 LWM = NEQ + 2     if MITER = 3, and
///                 LWM = (2*ML + MU + 1)*NEQ + 2
///                                   if MITER = 4 or 5.
///              (See the MF description below for METH and MITER.)
///
///              Thus if MAXORD has its default value and NEQ is constant,
///              this length is:
///              20 + 16*NEQ                    for MF = 10,
///              22 + 16*NEQ + NEQ**2           for MF = 11 or 12,
///              22 + 17*NEQ                    for MF = 13,
///              22 + 17*NEQ + (2*ML + MU)*NEQ  for MF = 14 or 15,
///              20 +  9*NEQ                    for MF = 20,
///              22 +  9*NEQ + NEQ**2           for MF = 21 or 22,
///              22 + 10*NEQ                    for MF = 23,
///              22 + 10*NEQ + (2*ML + MU)*NEQ  for MF = 24 or 25.
///
///              The first 20 words of RWORK are reserved for conditional
///              and optional inputs and optional outputs.
///
///              The following word in RWORK is a conditional input:
///              RWORK(1) = TCRIT, the critical value of t which the
///                         solver is not to overshoot.  Required if ITASK
///                         is 4 or 5, and ignored otherwise.  See ITASK.
///
///     LRW      The length of the array RWORK, as declared by the user.
///              (This will be checked by the solver.)
///
///     IWORK    An integer work array.  Its length must be at least
///              20       if MITER = 0 or 3 (MF = 10, 13, 20, 23), or
///              20 + NEQ otherwise (MF = 11, 12, 14, 15, 21, 22, 24, 25).
///              (See the MF description below for MITER.)  The first few
///              words of IWORK are used for conditional and optional
///              inputs and optional outputs.
///
///              The following two words in IWORK are conditional inputs:
///              IWORK(1) = ML   These are the lower and upper half-
///              IWORK(2) = MU   bandwidths, respectively, of the banded
///                              Jacobian, excluding the main diagonal.
///                         The band is defined by the matrix locations
///                         (i,j) with i - ML <= j <= i + MU. ML and MU
///                         must satisfy 0 <= ML,MU <= NEQ - 1. These are
///                         required if MITER is 4 or 5, and ignored
///                         otherwise.  ML and MU may in fact be the band
///                         parameters for a matrix to which df/dy is only
///                         approximately equal.
///
///     LIW      The length of the array IWORK, as declared by the user.
///              (This will be checked by the solver.)
///
///     Note:  The work arrays must not be altered between calls to DLSODE
///     for the same problem, except possibly for the conditional and
///     optional inputs, and except for the last 3*NEQ words of RWORK.
///     The latter space is used for internal scratch space, and so is
///     available for use by the user outside DLSODE between calls, if
///     desired (but not for use by F or JAC).
///
///     JAC      The name of the user-supplied routine (MITER = 1 or 4) to
///              compute the Jacobian matrix, df/dy, as a function of the
///              scalar t and the vector y.  (See the MF description below
///              for MITER.)  It is to have the form
///
///                 SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD)
///                 DOUBLE PRECISION T, Y(*), PD(NROWPD,*)
///
///              where NEQ, T, Y, ML, MU, and NROWPD are input and the
///              array PD is to be loaded with partial derivatives
///              (elements of the Jacobian matrix) on output.  PD must be
///              given a first dimension of NROWPD.  T and Y have the same
///              meaning as in subroutine F.
///
///              In the full matrix case (MITER = 1), ML and MU are
///              ignored, and the Jacobian is to be loaded into PD in
///              columnwise manner, with df(i)/dy(j) loaded into PD(i,j).
///
///              In the band matrix case (MITER = 4), the elements within
///              the band are to be loaded into PD in columnwise manner,
///              with diagonal lines of df/dy loaded into the rows of PD.
///              Thus df(i)/dy(j) is to be loaded into PD(i-j+MU+1,j).  ML
///              and MU are the half-bandwidth parameters (see IWORK).
///              The locations in PD in the two triangular areas which
///              correspond to nonexistent matrix elements can be ignored
///              or loaded arbitrarily, as they are overwritten by DLSODE.
///
///              JAC need not provide df/dy exactly. A crude approximation
///              (possibly with a smaller bandwidth) will do.
///
///              In either case, PD is preset to zero by the solver, so
///              that only the nonzero elements need be loaded by JAC.
///              Each call to JAC is preceded by a call to F with the same
///              arguments NEQ, T, and Y. Thus to gain some efficiency,
///              intermediate quantities shared by both calculations may
///              be saved in a user COMMON block by F and not recomputed
///              by JAC, if desired.  Also, JAC may alter the Y array, if
///              desired.  JAC must be declared EXTERNAL in the calling
///              program.
///
///              Subroutine JAC may access user-defined quantities in
///              NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
///              (dimensioned in JAC) and/or Y has length exceeding
///              NEQ(1).  See the descriptions of NEQ and Y above.
///
///     MF       The method flag.  Used only for input.  The legal values
///              of MF are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24,
///              and 25.  MF has decimal digits METH and MITER:
///                 MF = 10*METH + MITER .
///
///              METH indicates the basic linear multistep method:
///              1   Implicit Adams method.
///              2   Method based on backward differentiation formulas
///                  (BDF's).
///
///              MITER indicates the corrector iteration method:
///              0   Functional iteration (no Jacobian matrix is
///                  involved).
///              1   Chord iteration with a user-supplied full (NEQ by
///                  NEQ) Jacobian.
///              2   Chord iteration with an internally generated
///                  (difference quotient) full Jacobian (using NEQ
///                  extra calls to F per df/dy value).
///              3   Chord iteration with an internally generated
///                  diagonal Jacobian approximation (using one extra call
///                  to F per df/dy evaluation).
///              4   Chord iteration with a user-supplied banded Jacobian.
///              5   Chord iteration with an internally generated banded
///                  Jacobian (using ML + MU + 1 extra calls to F per
///                  df/dy evaluation).
///
///              If MITER = 1 or 4, the user must supply a subroutine JAC
///              (the name is arbitrary) as described above under JAC.
///              For other values of MITER, a dummy argument can be used.
///
///     Optional Inputs
///     ---------------
///     The following is a list of the optional inputs provided for in the
///     call sequence.  (See also Part 2.)  For each such input variable,
///     this table lists its name as used in this documentation, its
///     location in the call sequence, its meaning, and the default value.
///     The use of any of these inputs requires IOPT = 1, and in that case
///     all of these inputs are examined.  A value of zero for any of
///     these optional inputs will cause the default value to be used.
///     Thus to use a subset of the optional inputs, simply preload
///     locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively,
///     and then set those of interest to nonzero values.
///
///     Name    Location   Meaning and default value
///     ------  ---------  -----------------------------------------------
///     H0      RWORK(5)   Step size to be attempted on the first step.
///                        The default value is determined by the solver.
///     HMAX    RWORK(6)   Maximum absolute step size allowed.  The
///                        default value is infinite.
///     HMIN    RWORK(7)   Minimum absolute step size allowed.  The
///                        default value is 0.  (This lower bound is not
///                        enforced on the final step before reaching
///                        TCRIT when ITASK = 4 or 5.)
///     MAXORD  IWORK(5)   Maximum order to be allowed.  The default value
///                        is 12 if METH = 1, and 5 if METH = 2. (See the
///                        MF description above for METH.)  If MAXORD
///                        exceeds the default value, it will be reduced
///                        to the default value.  If MAXORD is changed
///                        during the problem, it may cause the current
///                        order to be reduced.
///     MXSTEP  IWORK(6)   Maximum number of (internally defined) steps
///                        allowed during one call to the solver.  The
///                        default value is 500.
///     MXHNIL  IWORK(7)   Maximum number of messages printed (per
///                        problem) warning that T + H = T on a step
///                        (H = step size).  This must be positive to
///                        result in a nondefault value.  The default
///                        value is 10.
///
///     Optional Outputs
///     ----------------
///     As optional additional output from DLSODE, the variables listed
///     below are quantities related to the performance of DLSODE which 
///     are available to the user.  These are communicated by way of the
///     work arrays, but also have internal mnemonic names as shown. 
///     Except where stated otherwise, all of these outputs are defined on
///     any successful return from DLSODE, and on any return with ISTATE =
///     -1, -2, -4, -5, or -6.  On an illegal input return (ISTATE = -3),
///     they will be unchanged from their existing values (if any), except
///     possibly for TOLSF, LENRW, and LENIW.  On any error return,
///     outputs relevant to the error will be defined, as noted below.
///
///     Name   Location   Meaning
///     -----  ---------  ------------------------------------------------
///     HU     RWORK(11)  Step size in t last used (successfully).
///     HCUR   RWORK(12)  Step size to be attempted on the next step.
///     TCUR   RWORK(13)  Current value of the independent variable which
///                       the solver has actually reached, i.e., the
///                       current internal mesh point in t. On output,
///                       TCUR will always be at least as far as the
///                       argument T, but may be farther (if interpolation
///                       was done).
///     TOLSF  RWORK(14)  Tolerance scale factor, greater than 1.0,
///                       computed when a request for too much accuracy
///                       was detected (ISTATE = -3 if detected at the
///                       start of the problem, ISTATE = -2 otherwise).
///                       If ITOL is left unaltered but RTOL and ATOL are
///                       uniformly scaled up by a factor of TOLSF for the
///                       next call, then the solver is deemed likely to
///                       succeed.  (The user may also ignore TOLSF and
///                       alter the tolerance parameters in any other way
///                       appropriate.)
///     NST    IWORK(11)  Number of steps taken for the problem so far.
///     NFE    IWORK(12)  Number of F evaluations for the problem so far.
///     NJE    IWORK(13)  Number of Jacobian evaluations (and of matrix LU
///                       decompositions) for the problem so far.
///     NQU    IWORK(14)  Method order last used (successfully).
///     NQCUR  IWORK(15)  Order to be attempted on the next step.
///     IMXER  IWORK(16)  Index of the component of largest magnitude in
///                       the weighted local error vector ( e(i)/EWT(i) ),
///                       on an error return with ISTATE = -4 or -5.
///     LENRW  IWORK(17)  Length of RWORK actually required.  This is
///                       defined on normal returns and on an illegal
///                       input return for insufficient storage.
///     LENIW  IWORK(18)  Length of IWORK actually required.  This is
///                       defined on normal returns and on an illegal
///                       input return for insufficient storage.
///
///     The following two arrays are segments of the RWORK array which may
///     also be of interest to the user as optional outputs.  For each
///     array, the table below gives its internal name, its base address
///     in RWORK, and its description.
///
///     Name  Base address  Description
///     ----  ------------  ----------------------------------------------
///     YH    21            The Nordsieck history array, of size NYH by
///                         (NQCUR + 1), where NYH is the initial value of
///                         NEQ.  For j = 0,1,...,NQCUR, column j + 1 of
///                         YH contains HCUR**j/factorial(j) times the jth
///                         derivative of the interpolating polynomial
///                         currently representing the solution, evaluated
///                         at t = TCUR.
///     ACOR  LENRW-NEQ+1   Array of size NEQ used for the accumulated
///                         corrections on each step, scaled on output to
///                         represent the estimated local error in Y on
///                         the last step.  This is the vector e in the
///                         description of the error control.  It is
///                         defined only on successful return from DLSODE.
///
///
///                    Part 2.  Other Callable Routines
///                    --------------------------------
///
///     The following are optional calls which the user may make to gain
///     additional capabilities in conjunction with DLSODE.
///
///     Form of call              Function
///     ------------------------  ----------------------------------------
///     CALL XSETUN(LUN)          Set the logical unit number, LUN, for
///                               output of messages from DLSODE, if the
///                               default is not desired.  The default
///                               value of LUN is 6. This call may be made
///                               at any time and will take effect
///                               immediately.
///     CALL XSETF(MFLAG)         Set a flag to control the printing of
///                               messages by DLSODE.  MFLAG = 0 means do
///                               not print.  (Danger:  this risks losing
///                               valuable information.)  MFLAG = 1 means
///                               print (the default).  This call may be
///                               made at any time and will take effect
///                               immediately.
///     CALL DSRCOM(RSAV,ISAV,JOB)  Saves and restores the contents of the
///                               internal COMMON blocks used by DLSODE
///                               (see Part 3 below).  RSAV must be a
///                               real array of length 218 or more, and
///                               ISAV must be an integer array of length
///                               37 or more.  JOB = 1 means save COMMON
///                               into RSAV/ISAV.  JOB = 2 means restore
///                               COMMON from same.  DSRCOM is useful if
///                               one is interrupting a run and restarting
///                               later, or alternating between two or
///                               more problems solved with DLSODE.
///     CALL DINTDY(,,,,,)        Provide derivatives of y, of various
///     (see below)               orders, at a specified point t, if
///                               desired.  It may be called only after a
///                               successful return from DLSODE.  Detailed
///                               instructions follow.
///
///     Detailed instructions for using DINTDY
///     --------------------------------------
///     The form of the CALL is:
///
///           CALL DINTDY (T, K, RWORK(21), NYH, DKY, IFLAG)
///
///     The input parameters are:
///
///     T          Value of independent variable where answers are
///                desired (normally the same as the T last returned by
///                DLSODE).  For valid results, T must lie between
///                TCUR - HU and TCUR.  (See "Optional Outputs" above
///                for TCUR and HU.)
///     K          Integer order of the derivative desired.  K must
///                satisfy 0 <= K <= NQCUR, where NQCUR is the current
///                order (see "Optional Outputs").  The capability
///                corresponding to K = 0, i.e., computing y(t), is
///                already provided by DLSODE directly.  Since
///                NQCUR >= 1, the first derivative dy/dt is always
///                available with DINTDY.
///     RWORK(21)  The base address of the history array YH.
///     NYH        Column length of YH, equal to the initial value of NEQ.
///
///     The output parameters are:
///
///     DKY        Real array of length NEQ containing the computed value
///                of the Kth derivative of y(t).
///     IFLAG      Integer flag, returned as 0 if K and T were legal,
///                -1 if K was illegal, and -2 if T was illegal.
///                On an error return, a message is also written.
///
///
///                          Part 3.  Common Blocks
///                          ----------------------
///
///     If DLSODE is to be used in an overlay situation, the user must
///     declare, in the primary overlay, the variables in:
///     (1) the call sequence to DLSODE,
///     (2) the internal COMMON block /DLS001/, of length 255 
///         (218 double precision words followed by 37 integer words).
///
///     If DLSODE is used on a system in which the contents of internal
///     COMMON blocks are not preserved between calls, the user should
///     declare the above COMMON block in his main program to insure that
///     its contents are preserved.
///
///     If the solution of a given problem by DLSODE is to be interrupted
///     and then later continued, as when restarting an interrupted run or
///     alternating between two or more problems, the user should save,
///     following the return from the last DLSODE call prior to the
///     interruption, the contents of the call sequence variables and the
///     internal COMMON block, and later restore these values before the
///     next DLSODE call for that problem.   In addition, if XSETUN and/or
///     XSETF was called for non-default handling of error messages, then
///     these calls must be repeated.  To save and restore the COMMON
///     block, use subroutine DSRCOM (see Part 2 above).
///
///
///              Part 4.  Optionally Replaceable Solver Routines
///              -----------------------------------------------
///
///     Below are descriptions of two routines in the DLSODE package which
///     relate to the measurement of errors.  Either routine can be
///     replaced by a user-supplied version, if desired.  However, since
///     such a replacement may have a major impact on performance, it
///     should be done only when absolutely necessary, and only with great
///     caution.  (Note:  The means by which the package version of a
///     routine is superseded by the user's version may be system-
///     dependent.)
///
///     DEWSET
///     ------
///     The following subroutine is called just before each internal
///     integration step, and sets the array of error weights, EWT, as
///     described under ITOL/RTOL/ATOL above:
///
///           SUBROUTINE DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
///
///     where NEQ, ITOL, RTOL, and ATOL are as in the DLSODE call
///     sequence, YCUR contains the current dependent variable vector,
///     and EWT is the array of weights set by DEWSET.
///
///     If the user supplies this subroutine, it must return in EWT(i)
///     (i = 1,...,NEQ) a positive quantity suitable for comparing errors
///     in Y(i) to.  The EWT array returned by DEWSET is passed to the
///     DVNORM routine (see below), and also used by DLSODE in the
///     computation of the optional output IMXER, the diagonal Jacobian
///     approximation, and the increments for difference quotient
///     Jacobians.
///
///     In the user-supplied version of DEWSET, it may be desirable to use
///     the current values of derivatives of y. Derivatives up to order NQ
///     are available from the history array YH, described above under
///     optional outputs.  In DEWSET, YH is identical to the YCUR array,
///     extended to NQ + 1 columns with a column length of NYH and scale
///     factors of H**j/factorial(j).  On the first call for the problem,
///     given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
///     NYH is the initial value of NEQ.  The quantities NQ, H, and NST
///     can be obtained by including in SEWSET the statements:
///           DOUBLE PRECISION RLS
///           COMMON /DLS001/ RLS(218),ILS(37)
///           NQ = ILS(33)
///           NST = ILS(34)
///           H = RLS(212)
///     Thus, for example, the current value of dy/dt can be obtained as
///     YCUR(NYH+i)/H (i=1,...,NEQ) (and the division by H is unnecessary
///     when NST = 0).
///
///     DVNORM
///     ------
///     DVNORM is a real function routine which computes the weighted
///     root-mean-square norm of a vector v:
///
///        d = DVNORM (n, v, w)
///
///     where:
///     n = the length of the vector,
///     v = real array of length n containing the vector,
///     w = real array of length n containing weights,
///     d = SQRT( (1/n) * sum(v(i)*w(i))**2 ).
///
///     DVNORM is called with n = NEQ and with w(i) = 1.0/EWT(i), where
///     EWT is as set by subroutine DEWSET.
///
///     If the user supplies this function, it should return a nonnegative
///     value of DVNORM suitable for use in the error control in DLSODE.
///     None of the arguments should be altered by DVNORM.  For example, a
///     user-supplied DVNORM routine might:
///     - Substitute a max-norm of (v(i)*w(i)) for the rms-norm, or
///     - Ignore some components of v in the norm, with the effect of
///       suppressing the error control on those components of Y.
    pub fn dlsode_(
        f: extern "C" fn(*const c_int, *const c_double, *mut c_double, *mut c_double),
        neq: &c_int,
        y: *mut c_double,
        t: &mut c_double,
        tout: &c_double,
        itol: &c_int,
        rtol: &c_double,
        atol: &c_double,
        itask: &c_int,
        istate: &mut c_int,
        iopt: &c_int,
        rwork: *mut c_double,
        lrw: &c_int,
        iwork: *mut c_int,
        liw: &c_int,
        jac: extern "C" fn(
            *const c_int,
            *const c_double,
            *const c_double,
            *const c_int,
            *const c_int,
            *mut c_double,
            *const c_int,
        ),
        mf: &c_int,
    );

/// Call `DLSODES` subroutine from ODEPACK
///
/// DLSODES solves the initial value problem for stiff or nonstiff
/// systems of first order ODEs,
///     dy/dt = f(t,y) ,  or, in component form,
///     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(NEQ)) (i = 1,...,NEQ).
/// DLSODES is a variant of the DLSODE package, and is intended for
/// problems in which the Jacobian matrix df/dy has an arbitrary
/// sparse structure (when the problem is stiff).
///
/// Authors:       Alan C. Hindmarsh
///                Center for Applied Scientific Computing, L-561
///                Lawrence Livermore National Laboratory
///                Livermore, CA 94551
/// and
///                Andrew H. Sherman
///                J. S. Nolen and Associates
///                Houston, TX 77084
///-----------------------------------------------------------------------
/// References:
/// 1.  Alan C. Hindmarsh,  ODEPACK, A Systematized Collection of ODE
///     Solvers, in Scientific Computing, R. S. Stepleman et al. (Eds.),
///     North-Holland, Amsterdam, 1983, pp. 55-64.
///
/// 2.  S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman,
///     Yale Sparse Matrix Package: I. The Symmetric Codes,
///     Int. J. Num. Meth. Eng., 18 (1982), pp. 1145-1151.
///
/// 3.  S. C. Eisenstat, M. C. Gursky, M. H. Schultz, and A. H. Sherman,
///     Yale Sparse Matrix Package: II. The Nonsymmetric Codes,
///     Research Report No. 114, Dept. of Computer Sciences, Yale
///     University, 1977.
///-----------------------------------------------------------------------
/// Summary of Usage.
///
/// Communication between the user and the DLSODES package, for normal
/// situations, is summarized here.  This summary describes only a subset
/// of the full set of options available.  See the full description for
/// details, including optional communication, nonstandard options,
/// and instructions for special situations.  See also the example
/// problem (with program and output) following this summary.
///
/// A. First provide a subroutine of the form:
///               SUBROUTINE F (NEQ, T, Y, YDOT)
///               DOUBLE PRECISION T, Y(*), YDOT(*)
/// which supplies the vector function f by loading YDOT(i) with f(i).
///
/// B. Next determine (or guess) whether or not the problem is stiff.
/// Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue
/// whose real part is negative and large in magnitude, compared to the
/// reciprocal of the t span of interest.  If the problem is nonstiff,
/// use a method flag MF = 10.  If it is stiff, there are two standard
/// choices for the method flag, MF = 121 and MF = 222.  In both cases,
/// DLSODES requires the Jacobian matrix in some form, and it treats this
/// matrix in general sparse form, with sparsity structure determined
/// internally.  (For options where the user supplies the sparsity
/// structure, see the full description of MF below.)
///
/// C. If the problem is stiff, you are encouraged to supply the Jacobian
/// directly (MF = 121), but if this is not feasible, DLSODES will
/// compute it internally by difference quotients (MF = 222).
/// If you are supplying the Jacobian, provide a subroutine of the form:
///               SUBROUTINE JAC (NEQ, T, Y, J, IAN, JAN, PDJ)
///               DOUBLE PRECISION T, Y(*), IAN(*), JAN(*), PDJ(*)
/// Here NEQ, T, Y, and J are input arguments, and the JAC routine is to
/// load the array PDJ (of length NEQ) with the J-th column of df/dy.
/// I.e., load PDJ(i) with df(i)/dy(J) for all relevant values of i.
/// The arguments IAN and JAN should be ignored for normal situations.
/// DLSODES will call the JAC routine with J = 1,2,...,NEQ.
/// Only nonzero elements need be loaded.  Usually, a crude approximation
/// to df/dy, possibly with fewer nonzero elements, will suffice.
///
/// D. Write a main program which calls Subroutine DLSODES once for
/// each point at which answers are desired.  This should also provide
/// for possible use of logical unit 6 for output of error messages by
/// DLSODES.  On the first call to DLSODES, supply arguments as follows:
/// F      = name of subroutine for right-hand side vector f.
///          This name must be declared External in calling program.
/// NEQ    = number of first order ODEs.
/// Y      = array of initial values, of length NEQ.
/// T      = the initial value of the independent variable t.
/// TOUT   = first point where output is desired (.ne. T).
/// ITOL   = 1 or 2 according as ATOL (below) is a scalar or array.
/// RTOL   = relative tolerance parameter (scalar).
/// ATOL   = absolute tolerance parameter (scalar or array).
///          The estimated local error in Y(i) will be controlled so as
///          to be roughly less (in magnitude) than
///             EWT(i) = RTOL*ABS(Y(i)) + ATOL     if ITOL = 1, or
///             EWT(i) = RTOL*ABS(Y(i)) + ATOL(i)  if ITOL = 2.
///          Thus the local error test passes if, in each component,
///          either the absolute error is less than ATOL (or ATOL(i)),
///          or the relative error is less than RTOL.
///          Use RTOL = 0.0 for pure absolute error control, and
///          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error
///          control.  Caution: actual (global) errors may exceed these
///          local tolerances, so choose them conservatively.
/// ITASK  = 1 for normal computation of output values of Y at t = TOUT.
/// ISTATE = integer flag (input and output).  Set ISTATE = 1.
/// IOPT   = 0 to indicate no optional inputs used.
/// RWORK  = real work array of length at least:
///             20 + 16*NEQ            for MF = 10,
///             20 + (2 + 1./LENRAT)*NNZ + (11 + 9./LENRAT)*NEQ
///                                    for MF = 121 or 222,
///          where:
///          NNZ    = the number of nonzero elements in the sparse
///                   Jacobian (if this is unknown, use an estimate), and
///          LENRAT = the real to integer wordlength ratio (usually 1 in
///                   single precision and 2 in double precision).
///          In any case, the required size of RWORK cannot generally
///          be predicted in advance if MF = 121 or 222, and the value
///          above is a rough estimate of a crude lower bound.  Some
///          experimentation with this size may be necessary.
///          (When known, the correct required length is an optional
///          output, available in IWORK(17).)
/// LRW    = declared length of RWORK (in user dimension).
/// IWORK  = integer work array of length at least 30.
/// LIW    = declared length of IWORK (in user dimension).
/// JAC    = name of subroutine for Jacobian matrix (MF = 121).
///          If used, this name must be declared External in calling
///          program.  If not used, pass a dummy name.
/// MF     = method flag.  Standard values are:
///          10  for nonstiff (Adams) method, no Jacobian used
///          121 for stiff (BDF) method, user-supplied sparse Jacobian
///          222 for stiff method, internally generated sparse Jacobian
/// Note that the main program must declare arrays Y, RWORK, IWORK,
/// and possibly ATOL.
///
/// E. The output from the first call (or any call) is:
///      Y = array of computed values of y(t) vector.
///      T = corresponding value of independent variable (normally TOUT).
/// ISTATE = 2  if DLSODES was successful, negative otherwise.
///          -1 means excess work done on this call (perhaps wrong MF).
///          -2 means excess accuracy requested (tolerances too small).
///          -3 means illegal input detected (see printed message).
///          -4 means repeated error test failures (check all inputs).
///          -5 means repeated convergence failures (perhaps bad Jacobian
///             supplied or wrong choice of MF or tolerances).
///          -6 means error weight became zero during problem. (Solution
///             component i vanished, and ATOL or ATOL(i) = 0.)
///          -7 means a fatal error return flag came from sparse solver
///             CDRV by way of DPRJS or DSOLSS.  Should never happen.
///          A return with ISTATE = -1, -4, or -5 may result from using
///          an inappropriate sparsity structure, one that is quite
///          different from the initial structure.  Consider calling
///          DLSODES again with ISTATE = 3 to force the structure to be
///          reevaluated.  See the full description of ISTATE below.
///
/// F. To continue the integration after a successful return, simply
/// reset TOUT and call DLSODES again.  No other parameters need be reset.
///
///-----------------------------------------------------------------------
/// Example Problem.
///
/// The following is a simple example problem, with the coding
/// needed for its solution by DLSODES.  The problem is from chemical
/// kinetics, and consists of the following 12 rate equations:
///    dy1/dt  = -rk1*y1
///    dy2/dt  = rk1*y1 + rk11*rk14*y4 + rk19*rk14*y5
///                - rk3*y2*y3 - rk15*y2*y12 - rk2*y2
///    dy3/dt  = rk2*y2 - rk5*y3 - rk3*y2*y3 - rk7*y10*y3
///                + rk11*rk14*y4 + rk12*rk14*y6
///    dy4/dt  = rk3*y2*y3 - rk11*rk14*y4 - rk4*y4
///    dy5/dt  = rk15*y2*y12 - rk19*rk14*y5 - rk16*y5
///    dy6/dt  = rk7*y10*y3 - rk12*rk14*y6 - rk8*y6
///    dy7/dt  = rk17*y10*y12 - rk20*rk14*y7 - rk18*y7
///    dy8/dt  = rk9*y10 - rk13*rk14*y8 - rk10*y8
///    dy9/dt  = rk4*y4 + rk16*y5 + rk8*y6 + rk18*y7
///    dy10/dt = rk5*y3 + rk12*rk14*y6 + rk20*rk14*y7
///                + rk13*rk14*y8 - rk7*y10*y3 - rk17*y10*y12
///                - rk6*y10 - rk9*y10
///    dy11/dt = rk10*y8
///    dy12/dt = rk6*y10 + rk19*rk14*y5 + rk20*rk14*y7
///                - rk15*y2*y12 - rk17*y10*y12
///
/// with rk1 = rk5 = 0.1,  rk4 = rk8 = rk16 = rk18 = 2.5,
///      rk10 = 5.0,  rk2 = rk6 = 10.0,  rk14 = 30.0,
///      rk3 = rk7 = rk9 = rk11 = rk12 = rk13 = rk19 = rk20 = 50.0,
///      rk15 = rk17 = 100.0.
///
/// The t interval is from 0 to 1000, and the initial conditions
/// are y1 = 1, y2 = y3 = ... = y12 = 0.  The problem is stiff.
///
/// The following coding solves this problem with DLSODES, using MF = 121
/// and printing results at t = .1, 1., 10., 100., 1000.  It uses
/// ITOL = 1 and mixed relative/absolute tolerance controls.
/// During the run and at the end, statistical quantities of interest
/// are printed (see optional outputs in the full description below).
///
///     EXTERNAL FEX, JEX
///     DOUBLE PRECISION ATOL, RTOL, RWORK, T, TOUT, Y
///     DIMENSION Y(12), RWORK(500), IWORK(30)
///     DATA LRW/500/, LIW/30/
///     NEQ = 12
///     DO 10 I = 1,NEQ
/// 10    Y(I) = 0.0D0
///     Y(1) = 1.0D0
///     T = 0.0D0
///     TOUT = 0.1D0
///     ITOL = 1
///     RTOL = 1.0D-4
///     ATOL = 1.0D-6
///     ITASK = 1
///     ISTATE = 1
///     IOPT = 0
///     MF = 121
///     DO 40 IOUT = 1,5
///       CALL DLSODES (FEX, NEQ, Y, T, TOUT, ITOL, RTOL, ATOL,
///    1     ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JEX, MF)
///       WRITE(6,30)T,IWORK(11),RWORK(11),(Y(I),I=1,NEQ)
/// 30    FORMAT(//' At t =',D11.3,4X,
///    1    ' No. steps =',I5,4X,' Last step =',D11.3/
///    2    '  Y array =  ',4D14.5/13X,4D14.5/13X,4D14.5)
///       IF (ISTATE .LT. 0) GO TO 80
///       TOUT = TOUT*10.0D0
/// 40    CONTINUE
///     LENRW = IWORK(17)
///     LENIW = IWORK(18)
///     NST = IWORK(11)
///     NFE = IWORK(12)
///     NJE = IWORK(13)
///     NLU = IWORK(21)
///     NNZ = IWORK(19)
///     NNZLU = IWORK(25) + IWORK(26) + NEQ
///     WRITE (6,70) LENRW,LENIW,NST,NFE,NJE,NLU,NNZ,NNZLU
/// 70  FORMAT(//' Required RWORK size =',I4,'   IWORK size =',I4/
///    1   ' No. steps =',I4,'   No. f-s =',I4,'   No. J-s =',I4,
///    2   '   No. LU-s =',I4/' No. of nonzeros in J =',I5,
///    3   '   No. of nonzeros in LU =',I5)
///     STOP
/// 80  WRITE(6,90)ISTATE
/// 90  FORMAT(///' Error halt.. ISTATE =',I3)
///     STOP
///     END
///
///     SUBROUTINE FEX (NEQ, T, Y, YDOT)
///     DOUBLE PRECISION T, Y, YDOT
///     DOUBLE PRECISION RK1, RK2, RK3, RK4, RK5, RK6, RK7, RK8, RK9,
///    1   RK10, RK11, RK12, RK13, RK14, RK15, RK16, RK17
///     DIMENSION Y(12), YDOT(12)
///     DATA RK1/0.1D0/, RK2/10.0D0/, RK3/50.0D0/, RK4/2.5D0/, RK5/0.1D0/,
///    1   RK6/10.0D0/, RK7/50.0D0/, RK8/2.5D0/, RK9/50.0D0/, RK10/5.0D0/,
///    2   RK11/50.0D0/, RK12/50.0D0/, RK13/50.0D0/, RK14/30.0D0/,
///    3   RK15/100.0D0/, RK16/2.5D0/, RK17/100.0D0/, RK18/2.5D0/,
///    4   RK19/50.0D0/, RK20/50.0D0/
///     YDOT(1)  = -RK1*Y(1)
///     YDOT(2)  = RK1*Y(1) + RK11*RK14*Y(4) + RK19*RK14*Y(5)
///    1           - RK3*Y(2)*Y(3) - RK15*Y(2)*Y(12) - RK2*Y(2)
///     YDOT(3)  = RK2*Y(2) - RK5*Y(3) - RK3*Y(2)*Y(3) - RK7*Y(10)*Y(3)
///    1           + RK11*RK14*Y(4) + RK12*RK14*Y(6)
///     YDOT(4)  = RK3*Y(2)*Y(3) - RK11*RK14*Y(4) - RK4*Y(4)
///     YDOT(5)  = RK15*Y(2)*Y(12) - RK19*RK14*Y(5) - RK16*Y(5)
///     YDOT(6)  = RK7*Y(10)*Y(3) - RK12*RK14*Y(6) - RK8*Y(6)
///     YDOT(7)  = RK17*Y(10)*Y(12) - RK20*RK14*Y(7) - RK18*Y(7)
///     YDOT(8)  = RK9*Y(10) - RK13*RK14*Y(8) - RK10*Y(8)
///     YDOT(9)  = RK4*Y(4) + RK16*Y(5) + RK8*Y(6) + RK18*Y(7)
///     YDOT(10) = RK5*Y(3) + RK12*RK14*Y(6) + RK20*RK14*Y(7)
///    1           + RK13*RK14*Y(8) - RK7*Y(10)*Y(3) - RK17*Y(10)*Y(12)
///    2           - RK6*Y(10) - RK9*Y(10)
///     YDOT(11) = RK10*Y(8)
///     YDOT(12) = RK6*Y(10) + RK19*RK14*Y(5) + RK20*RK14*Y(7)
///    1           - RK15*Y(2)*Y(12) - RK17*Y(10)*Y(12)
///     RETURN
///     END
///
///     SUBROUTINE JEX (NEQ, T, Y, J, IA, JA, PDJ)
///     DOUBLE PRECISION T, Y, PDJ
///     DOUBLE PRECISION RK1, RK2, RK3, RK4, RK5, RK6, RK7, RK8, RK9,
///    1   RK10, RK11, RK12, RK13, RK14, RK15, RK16, RK17
///     DIMENSION Y(12), IA(*), JA(*), PDJ(12)
///     DATA RK1/0.1D0/, RK2/10.0D0/, RK3/50.0D0/, RK4/2.5D0/, RK5/0.1D0/,
///    1   RK6/10.0D0/, RK7/50.0D0/, RK8/2.5D0/, RK9/50.0D0/, RK10/5.0D0/,
///    2   RK11/50.0D0/, RK12/50.0D0/, RK13/50.0D0/, RK14/30.0D0/,
///    3   RK15/100.0D0/, RK16/2.5D0/, RK17/100.0D0/, RK18/2.5D0/,
///    4   RK19/50.0D0/, RK20/50.0D0/
///     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), J
/// 1   PDJ(1) = -RK1
///     PDJ(2) = RK1
///     RETURN
/// 2   PDJ(2) = -RK3*Y(3) - RK15*Y(12) - RK2
///     PDJ(3) = RK2 - RK3*Y(3)
///     PDJ(4) = RK3*Y(3)
///     PDJ(5) = RK15*Y(12)
///     PDJ(12) = -RK15*Y(12)
///     RETURN
/// 3   PDJ(2) = -RK3*Y(2)
///     PDJ(3) = -RK5 - RK3*Y(2) - RK7*Y(10)
///     PDJ(4) = RK3*Y(2)
///     PDJ(6) = RK7*Y(10)
///     PDJ(10) = RK5 - RK7*Y(10)
///     RETURN
/// 4   PDJ(2) = RK11*RK14
///     PDJ(3) = RK11*RK14
///     PDJ(4) = -RK11*RK14 - RK4
///     PDJ(9) = RK4
///     RETURN
/// 5   PDJ(2) = RK19*RK14
///     PDJ(5) = -RK19*RK14 - RK16
///     PDJ(9) = RK16
///     PDJ(12) = RK19*RK14
///     RETURN
/// 6   PDJ(3) = RK12*RK14
///     PDJ(6) = -RK12*RK14 - RK8
///     PDJ(9) = RK8
///     PDJ(10) = RK12*RK14
///     RETURN
/// 7   PDJ(7) = -RK20*RK14 - RK18
///     PDJ(9) = RK18
///     PDJ(10) = RK20*RK14
///     PDJ(12) = RK20*RK14
///     RETURN
/// 8   PDJ(8) = -RK13*RK14 - RK10
///     PDJ(10) = RK13*RK14
///     PDJ(11) = RK10
/// 9   RETURN
/// 10  PDJ(3) = -RK7*Y(3)
///     PDJ(6) = RK7*Y(3)
///     PDJ(7) = RK17*Y(12)
///     PDJ(8) = RK9
///     PDJ(10) = -RK7*Y(3) - RK17*Y(12) - RK6 - RK9
///     PDJ(12) = RK6 - RK17*Y(12)
/// 11  RETURN
/// 12  PDJ(2) = -RK15*Y(2)
///     PDJ(5) = RK15*Y(2)
///     PDJ(7) = RK17*Y(10)
///     PDJ(10) = -RK17*Y(10)
///     PDJ(12) = -RK15*Y(2) - RK17*Y(10)
///     RETURN
///     END
///
/// The output of this program (on a Cray-1 in single precision)
/// is as follows:
///
///
/// At t =  1.000e-01     No. steps =   12     Last step =  1.515e-02
///  Y array =     9.90050e-01   6.28228e-03   3.65313e-03   7.51934e-07
///                1.12167e-09   1.18458e-09   1.77291e-12   3.26476e-07
///                5.46720e-08   9.99500e-06   4.48483e-08   2.76398e-06
///
///
/// At t =  1.000e+00     No. steps =   33     Last step =  7.880e-02
///  Y array =     9.04837e-01   9.13105e-03   8.20622e-02   2.49177e-05
///                1.85055e-06   1.96797e-06   1.46157e-07   2.39557e-05
///                3.26306e-05   7.21621e-04   5.06433e-05   3.05010e-03
///
///
/// At t =  1.000e+01     No. steps =   48     Last step =  1.239e+00
///  Y array =     3.67876e-01   3.68958e-03   3.65133e-01   4.48325e-05
///                6.10798e-05   4.33148e-05   5.90211e-05   1.18449e-04
///                3.15235e-03   3.56531e-03   4.15520e-03   2.48741e-01
///
///
/// At t =  1.000e+02     No. steps =   91     Last step =  3.764e+00
///  Y array =     4.44981e-05   4.42666e-07   4.47273e-04  -3.53257e-11
///                2.81577e-08  -9.67741e-11   2.77615e-07   1.45322e-07
///                1.56230e-02   4.37394e-06   1.60104e-02   9.52246e-01
///
///
/// At t =  1.000e+03     No. steps =  111     Last step =  4.156e+02
///  Y array =    -2.65492e-13   2.60539e-14  -8.59563e-12   6.29355e-14
///               -1.78066e-13   5.71471e-13  -1.47561e-12   4.58078e-15
///                1.56314e-02   1.37878e-13   1.60184e-02   9.52719e-01
///
///
/// Required RWORK size = 442   IWORK size =  30
/// No. steps = 111   No. f-s = 142   No. J-s =   2   No. LU-s =  20
/// No. of nonzeros in J =   44   No. of nonzeros in LU =   50
///
///-----------------------------------------------------------------------
/// Full Description of User Interface to DLSODES.
///
/// The user interface to DLSODES consists of the following parts.
///
/// 1.   The call sequence to Subroutine DLSODES, which is a driver
///      routine for the solver.  This includes descriptions of both
///      the call sequence arguments and of user-supplied routines.
///      Following these descriptions is a description of
///      optional inputs available through the call sequence, and then
///      a description of optional outputs (in the work arrays).
///
/// 2.   Descriptions of other routines in the DLSODES package that may be
///      (optionally) called by the user.  These provide the ability to
///      alter error message handling, save and restore the internal
///      Common, and obtain specified derivatives of the solution y(t).
///
/// 3.   Descriptions of Common blocks to be declared in overlay
///      or similar environments, or to be saved when doing an interrupt
///      of the problem and continued solution later.
///
/// 4.   Description of two routines in the DLSODES package, either of
///      which the user may replace with his/her own version, if desired.
///      These relate to the measurement of errors.
///
///-----------------------------------------------------------------------
/// Part 1.  Call Sequence.
///
/// The call sequence parameters used for input only are
///     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, JAC, MF,
/// and those used for both input and output are
///     Y, T, ISTATE.
/// The work arrays RWORK and IWORK are also used for conditional and
/// optional inputs and optional outputs.  (The term output here refers
/// to the return from Subroutine DLSODES to the user's calling program.)
///
/// The legality of input parameters will be thoroughly checked on the
/// initial call for the problem, but not checked thereafter unless a
/// change in input parameters is flagged by ISTATE = 3 on input.
///
/// The descriptions of the call arguments are as follows.
///
/// F      = the name of the user-supplied subroutine defining the
///          ODE system.  The system must be put in the first-order
///          form dy/dt = f(t,y), where f is a vector-valued function
///          of the scalar t and the vector y.  Subroutine F is to
///          compute the function f.  It is to have the form
///               SUBROUTINE F (NEQ, T, Y, YDOT)
///               DOUBLE PRECISION T, Y(*), YDOT(*)
///          where NEQ, T, and Y are input, and the array YDOT = f(t,y)
///          is output.  Y and YDOT are arrays of length NEQ.
///          Subroutine F should not alter y(1),...,y(NEQ).
///          F must be declared External in the calling program.
///
///          Subroutine F may access user-defined quantities in
///          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
///          (dimensioned in F) and/or Y has length exceeding NEQ(1).
///          See the descriptions of NEQ and Y below.
///
///          If quantities computed in the F routine are needed
///          externally to DLSODES, an extra call to F should be made
///          for this purpose, for consistent and accurate results.
///          If only the derivative dy/dt is needed, use DINTDY instead.
///
/// NEQ    = the size of the ODE system (number of first order
///          ordinary differential equations).  Used only for input.
///          NEQ may be decreased, but not increased, during the problem.
///          If NEQ is decreased (with ISTATE = 3 on input), the
///          remaining components of Y should be left undisturbed, if
///          these are to be accessed in F and/or JAC.
///
///          Normally, NEQ is a scalar, and it is generally referred to
///          as a scalar in this user interface description.  However,
///          NEQ may be an array, with NEQ(1) set to the system size.
///          (The DLSODES package accesses only NEQ(1).)  In either case,
///          this parameter is passed as the NEQ argument in all calls
///          to F and JAC.  Hence, if it is an array, locations
///          NEQ(2),... may be used to store other integer data and pass
///          it to F and/or JAC.  Subroutines F and/or JAC must include
///          NEQ in a Dimension statement in that case.
///
/// Y      = a real array for the vector of dependent variables, of
///          length NEQ or more.  Used for both input and output on the
///          first call (ISTATE = 1), and only for output on other calls.
///          on the first call, Y must contain the vector of initial
///          values.  On output, Y contains the computed solution vector,
///          evaluated at T.  If desired, the Y array may be used
///          for other purposes between calls to the solver.
///
///          This array is passed as the Y argument in all calls to
///          F and JAC.  Hence its length may exceed NEQ, and locations
///          Y(NEQ+1),... may be used to store other real data and
///          pass it to F and/or JAC.  (The DLSODES package accesses only
///          Y(1),...,Y(NEQ).)
///
/// T      = the independent variable.  On input, T is used only on the
///          first call, as the initial point of the integration.
///          on output, after each call, T is the value at which a
///          computed solution Y is evaluated (usually the same as TOUT).
///          On an error return, T is the farthest point reached.
///
/// TOUT   = the next value of t at which a computed solution is desired.
///          Used only for input.
///
///          When starting the problem (ISTATE = 1), TOUT may be equal
///          to T for one call, then should .ne. T for the next call.
///          For the initial T, an input value of TOUT .ne. T is used
///          in order to determine the direction of the integration
///          (i.e. the algebraic sign of the step sizes) and the rough
///          scale of the problem.  Integration in either direction
///          (forward or backward in t) is permitted.
///
///          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after
///          the first call (i.e. the first call with TOUT .ne. T).
///          Otherwise, TOUT is required on every call.
///
///          If ITASK = 1, 3, or 4, the values of TOUT need not be
///          monotone, but a value of TOUT which backs up is limited
///          to the current internal T interval, whose endpoints are
///          TCUR - HU and TCUR (see optional outputs, below, for
///          TCUR and HU).
///
/// ITOL   = an indicator for the type of error control.  See
///          description below under ATOL.  Used only for input.
///
/// RTOL   = a relative error tolerance parameter, either a scalar or
///          an array of length NEQ.  See description below under ATOL.
///          Input only.
///
/// ATOL   = an absolute error tolerance parameter, either a scalar or
///          an array of length NEQ.  Input only.
///
///             The input parameters ITOL, RTOL, and ATOL determine
///          the error control performed by the solver.  The solver will
///          control the vector E = (E(i)) of estimated local errors
///          in y, according to an inequality of the form
///                      RMS-norm of ( E(i)/EWT(i) )   .le.   1,
///          where       EWT(i) = RTOL(i)*ABS(Y(i)) + ATOL(i),
///          and the RMS-norm (root-mean-square norm) here is
///          RMS-norm(v) = SQRT(sum v(i)**2 / NEQ).  Here EWT = (EWT(i))
///          is a vector of weights which must always be positive, and
///          the values of RTOL and ATOL should all be non-negative.
///          The following table gives the types (scalar/array) of
///          RTOL and ATOL, and the corresponding form of EWT(i).
///
///             ITOL    RTOL       ATOL          EWT(i)
///              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
///              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
///              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
///              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)
///
///          When either of these parameters is a scalar, it need not
///          be dimensioned in the user's calling program.
///
///          If none of the above choices (with ITOL, RTOL, and ATOL
///          fixed throughout the problem) is suitable, more general
///          error controls can be obtained by substituting
///          user-supplied routines for the setting of EWT and/or for
///          the norm calculation.  See Part 4 below.
///
///          If global errors are to be estimated by making a repeated
///          run on the same problem with smaller tolerances, then all
///          components of RTOL and ATOL (i.e. of EWT) should be scaled
///          down uniformly.
///
/// ITASK  = an index specifying the task to be performed.
///          Input only.  ITASK has the following values and meanings.
///          1  means normal computation of output values of y(t) at
///             t = TOUT (by overshooting and interpolating).
///          2  means take one step only and return.
///          3  means stop at the first internal mesh point at or
///             beyond t = TOUT and return.
///          4  means normal computation of output values of y(t) at
///             t = TOUT but without overshooting t = TCRIT.
///             TCRIT must be input as RWORK(1).  TCRIT may be equal to
///             or beyond TOUT, but not behind it in the direction of
///             integration.  This option is useful if the problem
///             has a singularity at or beyond t = TCRIT.
///          5  means take one step, without passing TCRIT, and return.
///             TCRIT must be input as RWORK(1).
///
///          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT
///          (within roundoff), it will return T = TCRIT (exactly) to
///          indicate this (unless ITASK = 4 and TOUT comes before TCRIT,
///          in which case answers at t = TOUT are returned first).
///
/// ISTATE = an index used for input and output to specify the
///          the state of the calculation.
///
///          On input, the values of ISTATE are as follows.
///          1  means this is the first call for the problem
///             (initializations will be done).  See note below.
///          2  means this is not the first call, and the calculation
///             is to continue normally, with no change in any input
///             parameters except possibly TOUT and ITASK.
///             (If ITOL, RTOL, and/or ATOL are changed between calls
///             with ISTATE = 2, the new values will be used but not
///             tested for legality.)
///          3  means this is not the first call, and the
///             calculation is to continue normally, but with
///             a change in input parameters other than
///             TOUT and ITASK.  Changes are allowed in
///             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF,
///             the conditional inputs IA and JA,
///             and any of the optional inputs except H0.
///             In particular, if MITER = 1 or 2, a call with ISTATE = 3
///             will cause the sparsity structure of the problem to be
///             recomputed (or reread from IA and JA if MOSS = 0).
///          Note:  a preliminary call with TOUT = T is not counted
///          as a first call here, as no initialization or checking of
///          input is done.  (Such a call is sometimes useful for the
///          purpose of outputting the initial conditions.)
///          Thus the first call for which TOUT .ne. T requires
///          ISTATE = 1 on input.
///
///          On output, ISTATE has the following values and meanings.
///           1  means nothing was done; TOUT = T and ISTATE = 1 on input.
///           2  means the integration was performed successfully.
///          -1  means an excessive amount of work (more than MXSTEP
///              steps) was done on this call, before completing the
///              requested task, but the integration was otherwise
///              successful as far as T.  (MXSTEP is an optional input
///              and is normally 500.)  To continue, the user may
///              simply reset ISTATE to a value .gt. 1 and call again
///              (the excess work step counter will be reset to 0).
///              In addition, the user may increase MXSTEP to avoid
///              this error return (see below on optional inputs).
///          -2  means too much accuracy was requested for the precision
///              of the machine being used.  This was detected before
///              completing the requested task, but the integration
///              was successful as far as T.  To continue, the tolerance
///              parameters must be reset, and ISTATE must be set
///              to 3.  The optional output TOLSF may be used for this
///              purpose.  (Note: If this condition is detected before
///              taking any steps, then an illegal input return
///              (ISTATE = -3) occurs instead.)
///          -3  means illegal input was detected, before taking any
///              integration steps.  See written message for details.
///              Note:  If the solver detects an infinite loop of calls
///              to the solver with illegal input, it will cause
///              the run to stop.
///          -4  means there were repeated error test failures on
///              one attempted step, before completing the requested
///              task, but the integration was successful as far as T.
///              The problem may have a singularity, or the input
///              may be inappropriate.
///          -5  means there were repeated convergence test failures on
///              one attempted step, before completing the requested
///              task, but the integration was successful as far as T.
///              This may be caused by an inaccurate Jacobian matrix,
///              if one is being used.
///          -6  means EWT(i) became zero for some i during the
///              integration.  Pure relative error control (ATOL(i)=0.0)
///              was requested on a variable which has now vanished.
///              The integration was successful as far as T.
///          -7  means a fatal error return flag came from the sparse
///              solver CDRV by way of DPRJS or DSOLSS (numerical
///              factorization or backsolve).  This should never happen.
///              The integration was successful as far as T.
///
///          Note: an error return with ISTATE = -1, -4, or -5 and with
///          MITER = 1 or 2 may mean that the sparsity structure of the
///          problem has changed significantly since it was last
///          determined (or input).  In that case, one can attempt to
///          complete the integration by setting ISTATE = 3 on the next
///          call, so that a new structure determination is done.
///
///          Note:  since the normal output value of ISTATE is 2,
///          it does not need to be reset for normal continuation.
///          Also, since a negative input value of ISTATE will be
///          regarded as illegal, a negative output value requires the
///          user to change it, and possibly other inputs, before
///          calling the solver again.
///
/// IOPT   = an integer flag to specify whether or not any optional
///          inputs are being used on this call.  Input only.
///          The optional inputs are listed separately below.
///          IOPT = 0 means no optional inputs are being used.
///                   Default values will be used in all cases.
///          IOPT = 1 means one or more optional inputs are being used.
///
/// RWORK  = a work array used for a mixture of real (double precision)
///          and integer work space.
///          The length of RWORK (in real words) must be at least
///             20 + NYH*(MAXORD + 1) + 3*NEQ + LWM    where
///          NYH    = the initial value of NEQ,
///          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a
///                   smaller value is given as an optional input),
///          LWM = 0                                    if MITER = 0,
///          LWM = 2*NNZ + 2*NEQ + (NNZ+9*NEQ)/LENRAT   if MITER = 1,
///          LWM = 2*NNZ + 2*NEQ + (NNZ+10*NEQ)/LENRAT  if MITER = 2,
///          LWM = NEQ + 2                              if MITER = 3.
///          In the above formulas,
///          NNZ    = number of nonzero elements in the Jacobian matrix.
///          LENRAT = the real to integer wordlength ratio (usually 1 in
///                   single precision and 2 in double precision).
///          (See the MF description for METH and MITER.)
///          Thus if MAXORD has its default value and NEQ is constant,
///          the minimum length of RWORK is:
///             20 + 16*NEQ        for MF = 10,
///             20 + 16*NEQ + LWM  for MF = 11, 111, 211, 12, 112, 212,
///             22 + 17*NEQ        for MF = 13,
///             20 +  9*NEQ        for MF = 20,
///             20 +  9*NEQ + LWM  for MF = 21, 121, 221, 22, 122, 222,
///             22 + 10*NEQ        for MF = 23.
///          If MITER = 1 or 2, the above formula for LWM is only a
///          crude lower bound.  The required length of RWORK cannot
///          be readily predicted in general, as it depends on the
///          sparsity structure of the problem.  Some experimentation
///          may be necessary.
///
///          The first 20 words of RWORK are reserved for conditional
///          and optional inputs and optional outputs.
///
///          The following word in RWORK is a conditional input:
///            RWORK(1) = TCRIT = critical value of t which the solver
///                       is not to overshoot.  Required if ITASK is
///                       4 or 5, and ignored otherwise.  (See ITASK.)
///
/// LRW    = the length of the array RWORK, as declared by the user.
///          (This will be checked by the solver.)
///
/// IWORK  = an integer work array.  The length of IWORK must be at least
///             31 + NEQ + NNZ   if MOSS = 0 and MITER = 1 or 2, or
///             30               otherwise.
///          (NNZ is the number of nonzero elements in df/dy.)
///
///          In DLSODES, IWORK is used only for conditional and
///          optional inputs and optional outputs.
///
///          The following two blocks of words in IWORK are conditional
///          inputs, required if MOSS = 0 and MITER = 1 or 2, but not
///          otherwise (see the description of MF for MOSS).
///            IWORK(30+j) = IA(j)     (j=1,...,NEQ+1)
///            IWORK(31+NEQ+k) = JA(k) (k=1,...,NNZ)
///          The two arrays IA and JA describe the sparsity structure
///          to be assumed for the Jacobian matrix.  JA contains the row
///          indices where nonzero elements occur, reading in columnwise
///          order, and IA contains the starting locations in JA of the
///          descriptions of columns 1,...,NEQ, in that order, with
///          IA(1) = 1.  Thus, for each column index j = 1,...,NEQ, the
///          values of the row index i in column j where a nonzero
///          element may occur are given by
///            i = JA(k),  where   IA(j) .le. k .lt. IA(j+1).
///          If NNZ is the total number of nonzero locations assumed,
///          then the length of the JA array is NNZ, and IA(NEQ+1) must
///          be NNZ + 1.  Duplicate entries are not allowed.
///
/// LIW    = the length of the array IWORK, as declared by the user.
///          (This will be checked by the solver.)
///
/// Note:  The work arrays must not be altered between calls to DLSODES
/// for the same problem, except possibly for the conditional and
/// optional inputs, and except for the last 3*NEQ words of RWORK.
/// The latter space is used for internal scratch space, and so is
/// available for use by the user outside DLSODES between calls, if
/// desired (but not for use by F or JAC).
///
/// JAC    = name of user-supplied routine (MITER = 1 or MOSS = 1) to
///          compute the Jacobian matrix, df/dy, as a function of
///          the scalar t and the vector y.  It is to have the form
///               SUBROUTINE JAC (NEQ, T, Y, J, IAN, JAN, PDJ)
///               DOUBLE PRECISION T, Y(*), IAN(*), JAN(*), PDJ(*)
///          where NEQ, T, Y, J, IAN, and JAN are input, and the array
///          PDJ, of length NEQ, is to be loaded with column J
///          of the Jacobian on output.  Thus df(i)/dy(J) is to be
///          loaded into PDJ(i) for all relevant values of i.
///          Here T and Y have the same meaning as in Subroutine F,
///          and J is a column index (1 to NEQ).  IAN and JAN are
///          undefined in calls to JAC for structure determination
///          (MOSS = 1).  otherwise, IAN and JAN are structure
///          descriptors, as defined under optional outputs below, and
///          so can be used to determine the relevant row indices i, if
///          desired.
///               JAC need not provide df/dy exactly.  A crude
///          approximation (possibly with greater sparsity) will do.
///               In any case, PDJ is preset to zero by the solver,
///          so that only the nonzero elements need be loaded by JAC.
///          Calls to JAC are made with J = 1,...,NEQ, in that order, and
///          each such set of calls is preceded by a call to F with the
///          same arguments NEQ, T, and Y.  Thus to gain some efficiency,
///          intermediate quantities shared by both calculations may be
///          saved in a user Common block by F and not recomputed by JAC,
///          if desired.  JAC must not alter its input arguments.
///          JAC must be declared External in the calling program.
///               Subroutine JAC may access user-defined quantities in
///          NEQ(2),... and/or in Y(NEQ(1)+1),... if NEQ is an array
///          (dimensioned in JAC) and/or Y has length exceeding NEQ(1).
///          See the descriptions of NEQ and Y above.
///
/// MF     = the method flag.  Used only for input.
///          MF has three decimal digits-- MOSS, METH, MITER--
///             MF = 100*MOSS + 10*METH + MITER.
///          MOSS indicates the method to be used to obtain the sparsity
///          structure of the Jacobian matrix if MITER = 1 or 2:
///            MOSS = 0 means the user has supplied IA and JA
///                     (see descriptions under IWORK above).
///            MOSS = 1 means the user has supplied JAC (see below)
///                     and the structure will be obtained from NEQ
///                     initial calls to JAC.
///            MOSS = 2 means the structure will be obtained from NEQ+1
///                     initial calls to F.
///          METH indicates the basic linear multistep method:
///            METH = 1 means the implicit Adams method.
///            METH = 2 means the method based on Backward
///                     Differentiation Formulas (BDFs).
///          MITER indicates the corrector iteration method:
///            MITER = 0 means functional iteration (no Jacobian matrix
///                      is involved).
///            MITER = 1 means chord iteration with a user-supplied
///                      sparse Jacobian, given by Subroutine JAC.
///            MITER = 2 means chord iteration with an internally
///                      generated (difference quotient) sparse Jacobian
///                      (using NGP extra calls to F per df/dy value,
///                      where NGP is an optional output described below.)
///            MITER = 3 means chord iteration with an internally
///                      generated diagonal Jacobian approximation
///                      (using 1 extra call to F per df/dy evaluation).
///          If MITER = 1 or MOSS = 1, the user must supply a Subroutine
///          JAC (the name is arbitrary) as described above under JAC.
///          Otherwise, a dummy argument can be used.
///
///          The standard choices for MF are:
///            MF = 10  for a nonstiff problem,
///            MF = 21 or 22 for a stiff problem with IA/JA supplied
///                     (21 if JAC is supplied, 22 if not),
///            MF = 121 for a stiff problem with JAC supplied,
///                     but not IA/JA,
///            MF = 222 for a stiff problem with neither IA/JA nor
///                     JAC supplied.
///          The sparseness structure can be changed during the
///          problem by making a call to DLSODES with ISTATE = 3.
///-----------------------------------------------------------------------
/// Optional Inputs.
///
/// The following is a list of the optional inputs provided for in the
/// call sequence.  (See also Part 2.)  For each such input variable,
/// this table lists its name as used in this documentation, its
/// location in the call sequence, its meaning, and the default value.
/// The use of any of these inputs requires IOPT = 1, and in that
/// case all of these inputs are examined.  A value of zero for any
/// of these optional inputs will cause the default value to be used.
/// Thus to use a subset of the optional inputs, simply preload
/// locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and
/// then set those of interest to nonzero values.
///
/// Name    Location      Meaning and Default Value
///
/// H0      RWORK(5)  the step size to be attempted on the first step.
///                   The default value is determined by the solver.
///
/// HMAX    RWORK(6)  the maximum absolute step size allowed.
///                   The default value is infinite.
///
/// HMIN    RWORK(7)  the minimum absolute step size allowed.
///                   The default value is 0.  (This lower bound is not
///                   enforced on the final step before reaching TCRIT
///                   when ITASK = 4 or 5.)
///
/// SETH    RWORK(8)  the element threshhold for sparsity determination
///                   when MOSS = 1 or 2.  If the absolute value of
///                   an estimated Jacobian element is .le. SETH, it
///                   will be assumed to be absent in the structure.
///                   The default value of SETH is 0.
///
/// MAXORD  IWORK(5)  the maximum order to be allowed.  The default
///                   value is 12 if METH = 1, and 5 if METH = 2.
///                   If MAXORD exceeds the default value, it will
///                   be reduced to the default value.
///                   If MAXORD is changed during the problem, it may
///                   cause the current order to be reduced.
///
/// MXSTEP  IWORK(6)  maximum number of (internally defined) steps
///                   allowed during one call to the solver.
///                   The default value is 500.
///
/// MXHNIL  IWORK(7)  maximum number of messages printed (per problem)
///                   warning that T + H = T on a step (H = step size).
///                   This must be positive to result in a non-default
///                   value.  The default value is 10.
///-----------------------------------------------------------------------
/// Optional Outputs.
///
/// As optional additional output from DLSODES, the variables listed
/// below are quantities related to the performance of DLSODES
/// which are available to the user.  These are communicated by way of
/// the work arrays, but also have internal mnemonic names as shown.
/// Except where stated otherwise, all of these outputs are defined
/// on any successful return from DLSODES, and on any return with
/// ISTATE = -1, -2, -4, -5, or -6.  On an illegal input return
/// (ISTATE = -3), they will be unchanged from their existing values
/// (if any), except possibly for TOLSF, LENRW, and LENIW.
/// On any error return, outputs relevant to the error will be defined,
/// as noted below.
///
/// Name    Location      Meaning
///
/// HU      RWORK(11) the step size in t last used (successfully).
///
/// HCUR    RWORK(12) the step size to be attempted on the next step.
///
/// TCUR    RWORK(13) the current value of the independent variable
///                   which the solver has actually reached, i.e. the
///                   current internal mesh point in t.  On output, TCUR
///                   will always be at least as far as the argument
///                   T, but may be farther (if interpolation was done).
///
/// TOLSF   RWORK(14) a tolerance scale factor, greater than 1.0,
///                   computed when a request for too much accuracy was
///                   detected (ISTATE = -3 if detected at the start of
///                   the problem, ISTATE = -2 otherwise).  If ITOL is
///                   left unaltered but RTOL and ATOL are uniformly
///                   scaled up by a factor of TOLSF for the next call,
///                   then the solver is deemed likely to succeed.
///                   (The user may also ignore TOLSF and alter the
///                   tolerance parameters in any other way appropriate.)
///
/// NST     IWORK(11) the number of steps taken for the problem so far.
///
/// NFE     IWORK(12) the number of f evaluations for the problem so far,
///                   excluding those for structure determination
///                   (MOSS = 2).
///
/// NJE     IWORK(13) the number of Jacobian evaluations for the problem
///                   so far, excluding those for structure determination
///                   (MOSS = 1).
///
/// NQU     IWORK(14) the method order last used (successfully).
///
/// NQCUR   IWORK(15) the order to be attempted on the next step.
///
/// IMXER   IWORK(16) the index of the component of largest magnitude in
///                   the weighted local error vector ( E(i)/EWT(i) ),
///                   on an error return with ISTATE = -4 or -5.
///
/// LENRW   IWORK(17) the length of RWORK actually required.
///                   This is defined on normal returns and on an illegal
///                   input return for insufficient storage.
///
/// LENIW   IWORK(18) the length of IWORK actually required.
///                   This is defined on normal returns and on an illegal
///                   input return for insufficient storage.
///
/// NNZ     IWORK(19) the number of nonzero elements in the Jacobian
///                   matrix, including the diagonal (MITER = 1 or 2).
///                   (This may differ from that given by IA(NEQ+1)-1
///                   if MOSS = 0, because of added diagonal entries.)
///
/// NGP     IWORK(20) the number of groups of column indices, used in
///                   difference quotient Jacobian aproximations if
///                   MITER = 2.  This is also the number of extra f
///                   evaluations needed for each Jacobian evaluation.
///
/// NLU     IWORK(21) the number of sparse LU decompositions for the
///                   problem so far.
///
/// LYH     IWORK(22) the base address in RWORK of the history array YH,
///                   described below in this list.
///
/// IPIAN   IWORK(23) the base address of the structure descriptor array
///                   IAN, described below in this list.
///
/// IPJAN   IWORK(24) the base address of the structure descriptor array
///                   JAN, described below in this list.
///
/// NZL     IWORK(25) the number of nonzero elements in the strict lower
///                   triangle of the LU factorization used in the chord
///                   iteration (MITER = 1 or 2).
///
/// NZU     IWORK(26) the number of nonzero elements in the strict upper
///                   triangle of the LU factorization used in the chord
///                   iteration (MITER = 1 or 2).
///                   The total number of nonzeros in the factorization
///                   is therefore NZL + NZU + NEQ.
///
/// The following four arrays are segments of the RWORK array which
/// may also be of interest to the user as optional outputs.
/// For each array, the table below gives its internal name,
/// its base address, and its description.
/// For YH and ACOR, the base addresses are in RWORK (a real array).
/// The integer arrays IAN and JAN are to be obtained by declaring an
/// integer array IWK and identifying IWK(1) with RWORK(21), using either
/// an equivalence statement or a subroutine call.  Then the base
/// addresses IPIAN (of IAN) and IPJAN (of JAN) in IWK are to be obtained
/// as optional outputs IWORK(23) and IWORK(24), respectively.
/// Thus IAN(1) is IWK(IPIAN), etc.
///
/// Name    Base Address      Description
///
/// IAN    IPIAN (in IWK)  structure descriptor array of size NEQ + 1.
/// JAN    IPJAN (in IWK)  structure descriptor array of size NNZ.
///         (see above)    IAN and JAN together describe the sparsity
///                        structure of the Jacobian matrix, as used by
///                        DLSODES when MITER = 1 or 2.
///                        JAN contains the row indices of the nonzero
///                        locations, reading in columnwise order, and
///                        IAN contains the starting locations in JAN of
///                        the descriptions of columns 1,...,NEQ, in
///                        that order, with IAN(1) = 1.  Thus for each
///                        j = 1,...,NEQ, the row indices i of the
///                        nonzero locations in column j are
///                        i = JAN(k),  IAN(j) .le. k .lt. IAN(j+1).
///                        Note that IAN(NEQ+1) = NNZ + 1.
///                        (If MOSS = 0, IAN/JAN may differ from the
///                        input IA/JA because of a different ordering
///                        in each column, and added diagonal entries.)
///
/// YH      LYH            the Nordsieck history array, of size NYH by
///          (optional     (NQCUR + 1), where NYH is the initial value
///           output)      of NEQ.  For j = 0,1,...,NQCUR, column j+1
///                        of YH contains HCUR**j/factorial(j) times
///                        the j-th derivative of the interpolating
///                        polynomial currently representing the solution,
///                        evaluated at t = TCUR.  The base address LYH
///                        is another optional output, listed above.
///
/// ACOR     LENRW-NEQ+1   array of size NEQ used for the accumulated
///                        corrections on each step, scaled on output
///                        to represent the estimated local error in y
///                        on the last step.  This is the vector E  in
///                        the description of the error control.  It is
///                        defined only on a successful return from
///                        DLSODES.
///
///-----------------------------------------------------------------------
/// Part 2.  Other Routines Callable.
///
/// The following are optional calls which the user may make to
/// gain additional capabilities in conjunction with DLSODES.
/// (The routines XSETUN and XSETF are designed to conform to the
/// SLATEC error handling package.)
///
///     Form of Call                  Function
///   CALL XSETUN(LUN)          Set the logical unit number, LUN, for
///                             output of messages from DLSODES, if
///                             the default is not desired.
///                             The default value of LUN is 6.
///
///   CALL XSETF(MFLAG)         Set a flag to control the printing of
///                             messages by DLSODES.
///                             MFLAG = 0 means do not print. (Danger:
///                             This risks losing valuable information.)
///                             MFLAG = 1 means print (the default).
///
///                             Either of the above calls may be made at
///                             any time and will take effect immediately.
///
///   CALL DSRCMS(RSAV,ISAV,JOB) saves and restores the contents of
///                             the internal Common blocks used by
///                             DLSODES (see Part 3 below).
///                             RSAV must be a real array of length 224
///                             or more, and ISAV must be an integer
///                             array of length 71 or more.
///                             JOB=1 means save Common into RSAV/ISAV.
///                             JOB=2 means restore Common from RSAV/ISAV.
///                                DSRCMS is useful if one is
///                             interrupting a run and restarting
///                             later, or alternating between two or
///                             more problems solved with DLSODES.
///
///   CALL DINTDY(,,,,,)        Provide derivatives of y, of various
///        (see below)          orders, at a specified point t, if
///                             desired.  It may be called only after
///                             a successful return from DLSODES.
///
/// The detailed instructions for using DINTDY are as follows.
/// The form of the call is:
///
///   LYH = IWORK(22)
///   CALL DINTDY (T, K, RWORK(LYH), NYH, DKY, IFLAG)
///
/// The input parameters are:
///
/// T         = value of independent variable where answers are desired
///             (normally the same as the T last returned by DLSODES).
///             For valid results, T must lie between TCUR - HU and TCUR.
///             (See optional outputs for TCUR and HU.)
/// K         = integer order of the derivative desired.  K must satisfy
///             0 .le. K .le. NQCUR, where NQCUR is the current order
///             (See optional outputs).  The capability corresponding
///             to K = 0, i.e. computing y(T), is already provided
///             by DLSODES directly.  Since NQCUR .ge. 1, the first
///             derivative dy/dt is always available with DINTDY.
/// LYH       = the base address of the history array YH, obtained
///             as an optional output as shown above.
/// NYH       = column length of YH, equal to the initial value of NEQ.
///
/// The output parameters are:
///
/// DKY       = a real array of length NEQ containing the computed value
///             of the K-th derivative of y(t).
/// IFLAG     = integer flag, returned as 0 if K and T were legal,
///             -1 if K was illegal, and -2 if T was illegal.
///             On an error return, a message is also written.
///-----------------------------------------------------------------------
/// Part 3.  Common Blocks.
///
/// If DLSODES is to be used in an overlay situation, the user
/// must declare, in the primary overlay, the variables in:
///   (1) the call sequence to DLSODES, and
///   (2) the two internal Common blocks
///         /DLS001/  of length  255  (218 double precision words
///                      followed by 37 integer words),
///         /DLSS01/  of length  40  (6 double precision words
///                      followed by 34 integer words),
///
/// If DLSODES is used on a system in which the contents of internal
/// Common blocks are not preserved between calls, the user should
/// declare the above Common blocks in the calling program to insure
/// that their contents are preserved.
///
/// If the solution of a given problem by DLSODES is to be interrupted
/// and then later continued, such as when restarting an interrupted run
/// or alternating between two or more problems, the user should save,
/// following the return from the last DLSODES call prior to the
/// interruption, the contents of the call sequence variables and the
/// internal Common blocks, and later restore these values before the
/// next DLSODES call for that problem.  To save and restore the Common
/// blocks, use Subroutine DSRCMS (see Part 2 above).
///
///-----------------------------------------------------------------------
/// Part 4.  Optionally Replaceable Solver Routines.
///
/// Below are descriptions of two routines in the DLSODES package which
/// relate to the measurement of errors.  Either routine can be
/// replaced by a user-supplied version, if desired.  However, since such
/// a replacement may have a major impact on performance, it should be
/// done only when absolutely necessary, and only with great caution.
/// (Note: The means by which the package version of a routine is
/// superseded by the user's version may be system-dependent.)
///
/// (a) DEWSET.
/// The following subroutine is called just before each internal
/// integration step, and sets the array of error weights, EWT, as
/// described under ITOL/RTOL/ATOL above:
///     Subroutine DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT)
/// where NEQ, ITOL, RTOL, and ATOL are as in the DLSODES call sequence,
/// YCUR contains the current dependent variable vector, and
/// EWT is the array of weights set by DEWSET.
///
/// If the user supplies this subroutine, it must return in EWT(i)
/// (i = 1,...,NEQ) a positive quantity suitable for comparing errors
/// in y(i) to.  The EWT array returned by DEWSET is passed to the DVNORM
/// routine (see below), and also used by DLSODES in the computation
/// of the optional output IMXER, the diagonal Jacobian approximation,
/// and the increments for difference quotient Jacobians.
///
/// In the user-supplied version of DEWSET, it may be desirable to use
/// the current values of derivatives of y.  Derivatives up to order NQ
/// are available from the history array YH, described above under
/// optional outputs.  In DEWSET, YH is identical to the YCUR array,
/// extended to NQ + 1 columns with a column length of NYH and scale
/// factors of H**j/factorial(j).  On the first call for the problem,
/// given by NST = 0, NQ is 1 and H is temporarily set to 1.0.
/// NYH is the initial value of NEQ.  The quantities NQ, H, and NST
/// can be obtained by including in DEWSET the statements:
///     DOUBLE PRECISION RLS
///     COMMON /DLS001/ RLS(218),ILS(37)
///     NQ = ILS(33)
///     NST = ILS(34)
///     H = RLS(212)
/// Thus, for example, the current value of dy/dt can be obtained as
/// YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is
/// unnecessary when NST = 0).
///
/// (b) DVNORM.
/// The following is a real function routine which computes the weighted
/// root-mean-square norm of a vector v:
///     D = DVNORM (N, V, W)
/// where
///   N = the length of the vector,
///   V = real array of length N containing the vector,
///   W = real array of length N containing weights,
///   D = SQRT( (1/N) * sum(V(i)*W(i))**2 ).
/// DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where
/// EWT is as set by Subroutine DEWSET.
///
/// If the user supplies this function, it should return a non-negative
/// value of DVNORM suitable for use in the error control in DLSODES.
/// None of the arguments should be altered by DVNORM.
/// For example, a user-supplied DVNORM routine might:
///   -substitute a max-norm of (V(i)*W(i)) for the RMS-norm, or
///   -ignore some components of V in the norm, with the effect of
///    suppressing the error control on those components of y.
    pub fn dlsodes_(
        f: extern "C" fn(*const c_int, *const c_double, *mut c_double, *mut c_double),
        neq: &c_int,
        y: *mut c_double,
        t: &mut c_double,
        tout: &c_double,
        itol: &c_int,
        rtol: &c_double,
        atol: &c_double,
        itask: &c_int,
        istate: &mut c_int,
        iopt: &c_int,
        rwork: *mut c_double,
        lrw: &c_int,
        iwork: *mut c_int,
        liw: &c_int,
        jac: extern "C" fn(
            *const c_int,    // neq
            *const c_double, // t
            *const c_double, // y
            *const c_int,    // j
            *const c_int,    // ian
            *const c_int,    // jan
            *mut c_double,   //pdj
        ),
        mf: &c_int,
    );
}
