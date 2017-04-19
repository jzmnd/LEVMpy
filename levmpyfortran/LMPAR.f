C
C     LEVENBERG–MARQUARDT ALGORITHM FROM MINPACK
C
C       N: a positive integer input variable set to the order of r.
C
C       R: an n by n array. on input the full upper triangle
C          must contain the full upper triangle of the matrix r.
C          on output the full upper triangle is unaltered, and the
C          strict lower triangle contains the strict upper triangle
C          (transposed) of the upper triangular matrix s.
C
C       LDR: a positive integer input variable not less than n
C            which specifies the leading dimension of the array r.
C
C       IPVT: an integer input array of length n which defines the
C             permutation matrix p such that a*p = q*r. column j of p
C             is column ipvt(j) of the identity matrix.
C
C       DIAG: an input array of length n which must contain the
C             diagonal elements of the matrix d.
C
C       QTB: an input array of length n which must contain the first
C            n elements of the vector (q transpose)*b.
C
C       DELTA: a positive input variable which specifies an upper
C              bound on the euclidean norm of d*x.
C
C       PAR: a nonnegative variable. on input par contains an
C            initial estimate of the levenberg-marquardt parameter.
C            on output par contains the final estimate.
C
C       X: an output array of length n which contains the least
C          squares solution of the system a*x = b, sqrt(par)*d*x = 0,
C          for the output par.
C
C       SDIAG: is an output array of length n which contains the
C              diagonal elements of the upper triangular matrix s.
C
C       WA1 and WA2 are work arrays of length n.
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE LMPAR(N,R,LDR,IPVT,DIAG,QTB,DELTA,PAR,X,SIGMA,WA1,WA2)
      INTEGER N,LDR
      INTEGER IPVT(N)
      DOUBLE PRECISION DELTA,PAR
      DOUBLE PRECISION R(LDR,N),DIAG(N),QTB(N),X(N),SIGMA(N),WA1(N),
     + WA2(N)
      INTEGER I,ITER,J,JM1,JP1,K,L,NSING
      DOUBLE PRECISION DXNORM,FP,GNORM,PARC,PARL,PARU,P1,P001,SUM,
     + TEMP,ZERO
      DOUBLE PRECISION ENORM
      DATA P1,P001,ZERO /1.0D-1,1.0D-3,0.0D0/
C
C     COMPUTE AND STORE IN X THE GAUSS-NEWTON DIRECTION. IF THE
C     JACOBIAN IS RANK-DEFICIENT, OBTAIN A LEAST SQUARES SOLUTION.
C
      NSING = N
      DO 10 J = 1, N
        WA1(J) = QTB(J)
        IF (R(J,J) .EQ. ZERO .AND. NSING .EQ. N) NSING = J - 1
        IF (NSING .LT. N) WA1(J) = ZERO
10    CONTINUE
      IF (NSING .LT. 1) GO TO 50
      DO 40 K = 1, NSING
         J = NSING - K + 1
         WA1(J) = WA1(J)/R(J,J)
         TEMP = WA1(J)
         JM1 = J - 1
         IF (JM1 .LT. 1) GO TO 30
         DO 20 I = 1, JM1
            WA1(I) = WA1(I) - R(I,J)*TEMP
20       CONTINUE
30    CONTINUE
40    CONTINUE
50    CONTINUE
      DO 60 J = 1, N
         L = IPVT(J)
         X(L) = WA1(J)
60    CONTINUE
C
C     INITIALIZE THE ITERATION COUNTER.
C     EVALUATE THE FUNCTION AT THE ORIGIN, AND TEST
C     FOR ACCEPTANCE OF THE GAUSS-NEWTON DIRECTION.
C
      ITER = 0
      DO 70 J = 1, N
        WA2(J) = DIAG(J)*X(J)
70    CONTINUE
      DXNORM = ENORM(N,WA2)
      FP = DXNORM - DELTA
      IF (FP .LE. P1*DELTA) GO TO 220
C
C     IF THE JACOBIAN IS NOT RANK DEFICIENT, THE NEWTON
C     STEP PROVIDES A LOWER BOUND, PARL, FOR THE ZERO OF
C     THE FUNCTION. OTHERWISE SET THIS BOUND TO ZERO.
C
      PARL = ZERO
      IF (NSING .LT. N) GO TO 120
      DO 80 J = 1, N
        L = IPVT(J)
        WA1(J) = DIAG(L)*(WA2(L)/DXNORM)
80    CONTINUE
      DO 110 J = 1, N
        SUM = ZERO
        JM1 = J - 1
        IF (JM1 .LT. 1) GO TO 100
        DO 90 I = 1, JM1
          SUM = SUM + R(I,J)*WA1(I)
90       CONTINUE
100   CONTINUE
         WA1(J) = (WA1(J) - SUM)/R(J,J)
110   CONTINUE
      TEMP = ENORM(N,WA1)
      PARL = ((FP/DXNORM)/TEMP)/TEMP
120   CONTINUE
C
C     CALCULATE AN UPPER BOUND, PARU, FOR THE ZERO OF THE FUNCTION.
C
      DO 140 J = 1, N
        SUM = ZERO
        DO 130 I = 1, J
          SUM = SUM + R(I,J)*QTB(I)
130     CONTINUE
        L = IPVT(J)
        WA1(J) = SUM/DIAG(L)
140   CONTINUE
      GNORM = ENORM(N,WA1)
      PARU = GNORM/DELTA
C
C     IF THE INPUT PAR IS ZERO, REPLACE IT WITH AN IMPROVED ESTIMATE.
C
      IF (PAR .EQ. ZERO .AND. PARL .GT. ZERO) 
     + PAR = ((FP + DELTA)/DELTA)*PARL
      IF (PAR .EQ. ZERO .AND. PARL .EQ. ZERO) PAR = GNORM/DXNORM
C
C     BEGINNING OF AN ITERATION.
C
150   CONTINUE
      ITER = ITER + 1
C
C     IF PAR LIES OUTSIDE THE INTERVAL (PARL,PARU), SET PAR TO THE
C     MAXIMUM OF 0.001*PARU AND THE GEOMETRIC MEAN OF PARL AND PARU.
C
      IF (PAR .LE. PARL .OR. PAR .GE. PARU)
     + PAR = DMAX1(P001,DSQRT(PARL/PARU))*PARU
C
C     EVALUATE THE FUNCTION AT THE CURRENT VALUE OF PAR.
C
      TEMP = DSQRT(PAR)
      DO 160 J = 1, N
        WA1(J) = TEMP*DIAG(J)
160   CONTINUE
      CALL QRSOLV(N,R,LDR,IPVT,WA1,QTB,X,SIGMA,WA2)
      DO 170 J = 1, N
        WA2(J) = DIAG(J)*X(J)
170   CONTINUE
      DXNORM = ENORM(N,WA2)
      TEMP = FP
      FP = DXNORM - DELTA
C
C     IF THE FUNCTION IS SMALL ENOUGH, ACCEPT THE CURRENT VALUE
C     OF PAR. ALSO TEST FOR THE EXCEPTIONAL CASES WHERE PARL
C     IS ZERO OR THE NUMBER OF ITERATIONS HAS REACHED 10.
C
      IF (DABS(FP) .LE. P1*DELTA
     +       .OR. PARL .EQ. ZERO
     +      .AND. DABS(TEMP-FP) .LE. P001*DABS(FP)
     +       .OR. ITER .EQ. 10) GO TO 220
C
C        COMPUTE THE NEWTON CORRECTION.
C
      DO 180 J = 1, N
        L = IPVT(J)
        WA1(J) = DIAG(L)*(WA2(L)/DXNORM)
180   CONTINUE
      DO 210 J = 1, N
        WA1(J) = WA1(J)/SIGMA(J)
        TEMP = WA1(J)
        JP1 = J + 1
        IF (N .LT. JP1) GO TO 200
        DO 190 I = JP1, N
          WA1(I) = WA1(I) - R(I,J)*TEMP
190     CONTINUE
200   CONTINUE
210   CONTINUE
      TEMP = ENORM(N,WA1)
      PARC = ((FP/DXNORM)/TEMP)/TEMP
C
C     DEPENDING ON THE SIGN OF THE FUNCTION, UPDATE PARL OR PARU.
C
      IF (FP .GT. ZERO) PARL = PAR
      IF (FP .LT. ZERO) PARU = PAR
C
C     COMPUTE AN IMPROVED ESTIMATE FOR PAR.
C
      PAR = PAR + ((FP + DELTA)/DELTA)*PARC
C
C     END OF AN ITERATION.
C
      GO TO 150
220   CONTINUE
C
C     TERMINATION.
C
      IF (ITER .EQ. 0) PAR = ZERO
      RETURN
      END
