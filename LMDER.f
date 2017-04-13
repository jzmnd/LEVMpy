C
C     Levenberg–Marquardt algorithm from MINPACK
C
C       The purpose of LMDER is to minimize the sum of the squares of
C       M nonlinear functions in N variables by a modification of
C       the Levenberg-Marquardt algorithm.
C       The user must provide a subroutine which calculates the
C       functions and the Jacobian.
C
C       FUNC: the name of the user-supplied subroutine which calculates
C             the functions and the Jacobian.
C
C       M: a positive integer input variable set to the number
C          of functions.
C
C       N: a positive integer input variable set to the number
C          of variables. n must not exceed m.
C
C       X: an array of length n. on input x must contain
C          an initial estimate of the solution vector. on output x
C          contains the final estimate of the solution vector.
C
C       FVEC: an output array of length m which contains
C             the functions evaluated at the output x.
C
C       FJAC: an output m by n array. the upper n by n submatrix
C             of fjac contains an upper triangular matrix r with
C             diagonal elements of nonincreasing magnitude such that
C
C                t     t           t
C               p *(jac *jac)*p = r *r,
C
C             where p is a permutation matrix and jac is the final
C             calculated jacobian. column j of p is column ipvt(j)
C             (see below) of the identity matrix. the lower trapezoidal
C             part of fjac contains information generated during
C             the computation of r.
C
C       LDFJAC: a positive integer input variable not less than m
C               which specifies the leading dimension of the array fjac.
C
C       FTOL: a nonnegative input variable. termination
C             occurs when both the actual and predicted relative
C             reductions in the sum of squares are at most ftol.
C             therefore, ftol measures the relative error desired
C             in the sum of squares.
C
C       XTOL: a nonnegative input variable. termination
C             occurs when the relative error between two consecutive
C             iterates is at most xtol. therefore, xtol measures the
C             relative error desired in the approximate solution.
C
C       GTOL: a nonnegative input variable. termination
C         occurs when the cosine of the angle between fvec and
C         any column of the jacobian is at most gtol in absolute
C         value. therefore, gtol measures the orthogonality
C         desired between the function vector and the columns
C         of the jacobian.
C
C       MAXFEV: a positive integer input variable. termination
C               occurs when the number of calls to fcn with iflag = 1
C               has reached maxfev.
C
C       DIAG: an array of length n. if mode = 1 (see
C             below), diag is internally set. if mode = 2, diag
C             must contain positive entries that serve as
C             multiplicative scale factors for the variables.
C
C       MODE: an integer input variable. if mode = 1, the
C             variables will be scaled internally. if mode = 2,
C             the scaling is specified by the input diag. other
C             values of mode are equivalent to mode = 1.
C
C       FACTOR: a positive input variable used in determining the
C         initial step bound. this bound is set to the product of
C         factor and the euclidean norm of diag*x if nonzero, or else
C         to factor itself. in most cases factor should lie in the
C         interval (.1,100.).100. is a generally recommended value.
C
C       NPRINT: an integer input variable that enables controlled
C         printing of iterates if it is positive. in this case,
C         fcn is called with iflag = 0 at the beginning of the first
C         iteration and every nprint iterations thereafter and
C         immediately prior to return, with x, fvec, and fjac
C         available for printing. fvec and fjac should not be
C         altered. if nprint is not positive, no special calls
C         of fcn with iflag = 0 are made.
C
C       INFO: an integer output variable. if the user has
C         terminated execution, info is set to the (negative)
C         value of iflag. see description of fcn. otherwise,
C         info is set as follows.
C
C         info = 0  improper input parameters.
C
C         info = 1  both actual and predicted relative reductions
C                   in the sum of squares are at most ftol.
C
C         info = 2  relative error between two consecutive iterates
C                   is at most xtol.
C
C         info = 3  conditions for info = 1 and info = 2 both hold.
C
C         info = 4  the cosine of the angle between fvec and any
C                   column of the jacobian is at most gtol in
C                   absolute value.
C
C         info = 5  number of calls to fcn with iflag = 1 has
C                   reached maxfev.
C
C         info = 6  ftol is too small. no further reduction in
C                   the sum of squares is possible.
C
C         info = 7  xtol is too small. no further improvement in
C                   the approximate solution x is possible.
C
C         info = 8  gtol is too small. fvec is orthogonal to the
C                   columns of the jacobian to machine precision.
C
C       NFEV: an integer output variable set to the number of
C             calls to fcn with iflag = 1.
C
C       NJEV: an integer output variable set to the number of
C             calls to fcn with iflag = 2.
C
C       IPVT: an integer output array of length n. ipvt
C             defines a permutation matrix p such that jac*p = q*r,
C             where jac is the final calculated jacobian, q is
C             orthogonal (not stored), and r is upper triangular
C             with diagonal elements of nonincreasing magnitude.
C             column j of p is column ipvt(j) of the identity matrix.
C
C       QTF: an output array of length n which contains
C            the first n elements of the vector (q transpose)*fvec.
C
C       WA1, WA2, and WA3 are work arrays of length n.
C       WA4 is a work array of length m.
C
      SUBROUTINE LMDER(FUNC,M,N,X,FVEC,FJAC,LDFJAC,FTOL,XTOL,GTOL,
     + MAXFEV,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,
     + IPVT,QTF,WA1,WA2,WA3,WA4)
      EXTERNAL FUNC
      INTEGER M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV
      INTEGER IPVT(N)
      DOUBLE PRECISION FTOL,XTOL,GTOL,FACTOR
      DOUBLE PRECISION X(*),FVEC(*),FJAC(LDFJAC,N),DIAG(*),QTF(*),
     + WA1(*),WA2(*),WA3(*),WA4(M)
      INTEGER I,IFLAG,ITER,J,L
      DOUBLE PRECISION ACTRED,DELTA,DIRDER,EPSMCH,FNORM,FNORM1,GNORM,
     + ONE,PAR,PNORM,PRERED,P1,P5,P25,P75,P0001,RATIO,
     + SUM,TEMP,TEMP1,TEMP2,XNORM,ZERO
      DOUBLE PRECISION DPMPAR,ENORM
      COMMON /CM4/ FNORM
      DATA ONE,P1,P5,P25,P75,P0001,ZERO
     + /1.0D0,1.0D-1,5.0D-1,2.5D-1,7.5D-1,1.0D-4,0.0D0/
C
      WRITE(*,*) '      == LMDER =='
C
C     EPSMCH IS THE MACHINE PRECISION.
C
      EPSMCH = DPMPAR(1)
C
      INFO = 0
      IFLAG = 0
      NFEV = 0
      NJEV = 0
C
C     CHECK THE INPUT PARAMETERS FOR ERRORS.
C
      IF (N .LE. 0 .OR. M .LT. N .OR. LDFJAC .LT. M
     +    .OR. FTOL .LT. ZERO .OR. XTOL .LT. ZERO .OR. GTOL .LT. ZERO
     +    .OR. MAXFEV .LE. 0 .OR. FACTOR .LE. ZERO) GO TO 290
      IF (MODE .NE. 2) GO TO 20
      DO 10 J = 1, N
         IF (DIAG(J) .LE. ZERO) GO TO 290
10    CONTINUE
20    CONTINUE
C
C     EVALUATE THE FUNCTION AT THE STARTING POINT
C     AND CALCULATE ITS NORM.
C
      IFLAG = 1
      CALL FUNC(M,N,X,FVEC,FJAC,LDFJAC,IFLAG)
      NFEV = 1
      IF (IFLAG .LT. 0) GO TO 290
      FNORM = ENORM(M,FVEC)
C
C     INITIALIZE LEVENBERG-MARQUARDT PARAMETER AND ITERATION COUNTER.
C
      PAR = ZERO
      ITER = 0
C
C     BEGINNING OF THE OUTER LOOP.
C
30    CONTINUE
      ITER = ITER + 1
C
C     CALCULATE THE JACOBIAN MATRIX.
C
      IFLAG = 2
      CALL FUNC(M,N,X,FVEC,FJAC,LDFJAC,IFLAG)
      NJEV = NJEV + 1
      IF (IFLAG .LT. 0) GO TO 290
C
C     IF REQUESTED, CALL FUNC TO ENABLE PRINTING OF ITERATES.
C
      IF (NPRINT .LE. 0) GO TO 35
      IFLAG = 0
      IF(NPRINT*(ITER/NPRINT).EQ.ITER)
     +CALL FUNC(M,N,X,FVEC,FJAC,LDFJAC,IFLAG)
      IF (IFLAG .LT. 0) GO TO 290
35    CONTINUE
C
C     COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
C
      CALL QRFAC(M,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3)
C
C     ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
C     TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
C
      IF (ITER .NE. 1) GO TO 70
      IF (MODE .EQ. 2) GO TO 50
      DO 40 J = 1, N
        DIAG(J) = WA2(J)
        IF (WA2(J) .EQ. ZERO) DIAG(J) = ONE
40    CONTINUE
50    CONTINUE
C
C     ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED X
C     AND INITIALIZE THE STEP BOUND DELTA.
C
      DO 60 J = 1, N
        WA3(J) = DIAG(J)*X(J)
60    CONTINUE
      XNORM = ENORM(N,WA3)
      DELTA = FACTOR*XNORM
      IF (DELTA .EQ. ZERO) DELTA = FACTOR
70    CONTINUE
C
C     FORM (Q TRANSPOSE)*FVEC AND STORE THE FIRST N COMPONENTS IN QTF.
C
      DO 80 I = 1, M
        WA4(I) = FVEC(I)
80    CONTINUE
      DO 120 J = 1, N
        IF (FJAC(J,J) .EQ. ZERO) GO TO 110
        SUM = ZERO
        DO 90 I = J, M
          SUM = SUM + FJAC(I,J)*WA4(I)
90      CONTINUE
        TEMP = -SUM/FJAC(J,J)
        DO 100 I = J, M
          WA4(I) = WA4(I) + FJAC(I,J)*TEMP
100     CONTINUE
110   CONTINUE
      FJAC(J,J) = WA1(J)
      QTF(J) = WA4(J)
120   CONTINUE
C
C     COMPUTE THE NORM OF THE SCALED GRADIENT.
C
      GNORM = ZERO
      IF (FNORM .EQ. ZERO) GO TO 160
        DO 150 J = 1, N
          L = IPVT(J)
          IF (WA2(L) .EQ. ZERO) GO TO 140
          SUM = ZERO
          DO 130 I = 1, J
            SUM = SUM + FJAC(I,J)*(QTF(I)/FNORM)
130       CONTINUE
          GNORM = DMAX1(GNORM,DABS(SUM/WA2(L)))
140     CONTINUE
150     CONTINUE
160   CONTINUE
C
C     TEST FOR CONVERGENCE OF THE GRADIENT NORM.
C
      IF (GNORM .LE. GTOL) INFO = 4
      IF (INFO .NE. 0) GO TO 290
C
C     RESCALE IF NECESSARY.
C
      IF (MODE .EQ. 2) GO TO 180
      DO 170 J = 1, N
        DIAG(J) = DMAX1(DIAG(J),WA2(J))
170   CONTINUE
180   CONTINUE
C
C        BEGINNING OF THE INNER LOOP.
C
190   CONTINUE
C
C        DETERMINE THE LEVENBERG-MARQUARDT PARAMETER.
C
      CALL LMPAR(N,FJAC,LDFJAC,IPVT,DIAG,QTF,DELTA,PAR,WA1,WA2,WA3,WA4)
C       STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
        DO 200 J = 1, N
               WA1(J) = -WA1(J)
               WA2(J) = X(J) + WA1(J)
               WA3(J) = DIAG(J)*WA1(J)
200     CONTINUE
        PNORM = ENORM(N,WA3)
C
C       ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
C
        IF (ITER .EQ. 1) DELTA = DMIN1(DELTA,PNORM)
C
C       EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
C
        IFLAG = 1
        CALL FUNC(M,N,WA2,WA4,FJAC,LDFJAC,IFLAG)
        NFEV = NFEV + 1
        IF (IFLAG .LT. 0) GO TO 290
        FNORM1 = ENORM(M,WA4)
C
C         COMPUTE THE SCALED ACTUAL REDUCTION.
C
        ACTRED = -ONE
        IF (P1*FNORM1 .LT. FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
C
C         COMPUTE THE SCALED PREDICTED REDUCTION AND
C         THE SCALED DIRECTIONAL DERIVATIVE.
C
        DO 220 J = 1, N
          WA3(J) = ZERO
          L = IPVT(J)
          TEMP = WA1(L)
          DO 210 I = 1, J
            WA3(I) = WA3(I) + FJAC(I,J)*TEMP
210       CONTINUE
220     CONTINUE
        TEMP1 = ENORM(N,WA3)/FNORM
        TEMP2 = (DSQRT(PAR)*PNORM)/FNORM
        PRERED = TEMP1**2 + TEMP2**2/P5
        DIRDER = -(TEMP1**2 + TEMP2**2)
C
C       COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED REDUCTION.
C
        RATIO = ZERO
        IF (PRERED .NE. ZERO) RATIO = ACTRED/PRERED
C
C       UPDATE THE STEP BOUND.
C
        IF (RATIO .GT. P25) GO TO 230
        IF (ACTRED .GE. ZERO) TEMP = P5
        IF (ACTRED .LT. ZERO) TEMP = P5*DIRDER/(DIRDER + P5*ACTRED)
        IF (P1*FNORM1 .GE. FNORM .OR. TEMP .LT. P1) TEMP = P1
        DELTA = TEMP*DMIN1(DELTA,PNORM/P1)
        PAR = PAR/TEMP
        GO TO 250
230     CONTINUE
        IF (PAR .NE. ZERO .AND. RATIO .LT. P75) GO TO 240
        DELTA = PNORM/P5
        PAR = P5*PAR
240     CONTINUE
250   CONTINUE
C
C     TEST FOR SUCCESSFUL ITERATION.
C
      IF (RATIO .LT. P0001) GO TO 280
C
C     SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
C
      DO 260 J = 1, N
         X(J) = WA2(J)
         WA2(J) = DIAG(J)*X(J)
260   CONTINUE
      DO 270 I = 1, M
         FVEC(I) = WA4(I)
270   CONTINUE
      XNORM = ENORM(N,WA2)
      FNORM = FNORM1
280   CONTINUE
C
C     TESTS FOR CONVERGENCE.
C
      IF (DABS(ACTRED) .LE. FTOL .AND. PRERED .LE. FTOL
     +      .AND. P5*RATIO .LE. ONE) INFO = 1
      IF (DELTA .LE. XTOL*XNORM) INFO = 2
      IF (DABS(ACTRED) .LE. FTOL .AND. PRERED .LE. FTOL
     +      .AND. P5*RATIO .LE. ONE .AND. INFO .EQ. 2) INFO = 3
      IF (INFO .NE. 0) GO TO 290
C
C     TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
C     Levm will stop iterating when NFEV .GE. MAXFEV
C
      IF (NFEV .GE. MAXFEV) INFO = 5
      IF (DABS(ACTRED) .LE. EPSMCH .AND. PRERED .LE. EPSMCH
     +      .AND. P5*RATIO .LE. ONE) INFO = 6
      IF (DELTA .LE. EPSMCH*XNORM) INFO = 7
      IF (GNORM .LE. EPSMCH) INFO = 8
      IF (INFO .NE. 0) GO TO 290
C
C     END OF THE INNER LOOP. REPEAT IF ITERATION UNSUCCESSFUL.
C
      IF (RATIO .LT. P0001) GO TO 190
C
C     END OF THE OUTER LOOP.
C
      GO TO 30
290   CONTINUE
C
C     TERMINATION, EITHER NORMAL OR USER IMPOSED.
C
      IF (IFLAG .LT. 0) INFO = IFLAG
      RETURN
      END
