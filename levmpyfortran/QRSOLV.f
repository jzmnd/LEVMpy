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
C       X: an output array of length n which contains the least
C          squares solution of the system a*x = b, d*x = 0.
C
C       SDIAG: an output array of length n which contains the
C              diagonal elements of the upper triangular matrix s.
C
C       WA is a work array of length n.
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE QRSOLV(N,R,LDR,IPVT,DIAG,QTB,X,SIGMA,WA)
      INTEGER N,LDR
      INTEGER IPVT(N)
      DOUBLE PRECISION R(LDR,N),DIAG(N),QTB(N),X(N),SIGMA(N),WA(N)
      INTEGER I,J,JP1,K,KP1,L,NSING
      DOUBLE PRECISION COS,COTAN,P5,P25,QTBPJ,SIN,SUM,TAN,TEMP,ZERO
      DATA P5,P25,ZERO /5.0D-1,2.5D-1,0.0D0/
C
C     COPY R AND (Q TRANSPOSE)*B TO PRESERVE INPUT AND INITIALIZE S.
C     IN PARTICULAR, SAVE THE DIAGONAL ELEMENTS OF R IN X.
C
      DO 20 J = 1, N
        DO 10 I = J, N
          R(I,J) = R(J,I)
10      CONTINUE
        X(J) = R(J,J)
        WA(J) = QTB(J)
20    CONTINUE
C
C     ELIMINATE THE DIAGONAL MATRIX D USING A GIVENS ROTATION.
C
      DO 100 J = 1, N
C
C       PREPARE THE ROW OF D TO BE ELIMINATED, LOCATING THE
C       DIAGONAL ELEMENT USING P FROM THE QR FACTORIZATION.
C
        L = IPVT(J)
        IF (DIAG(L) .EQ. ZERO) GO TO 90
        DO 30 K = J, N
          SIGMA(K) = ZERO
30      CONTINUE
        SIGMA(J) = DIAG(L)
C
C       THE TRANSFORMATIONS TO ELIMINATE THE ROW OF D
C       MODIFY ONLY A SINGLE ELEMENT OF (Q TRANSPOSE)*B
C       BEYOND THE FIRST N, WHICH IS INITIALLY ZERO.
C
        QTBPJ = ZERO
        DO 80 K = J, N
C
C         DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
C         APPROPRIATE ELEMENT IN THE CURRENT ROW OF D.
C
          IF (SIGMA(K) .EQ. ZERO) GO TO 70
          IF (DABS(R(K,K)) .GE. DABS(SIGMA(K))) GO TO 40
          COTAN = R(K,K)/SIGMA(K)
          SIN = P5/DSQRT(P25+P25*COTAN**2)
          COS = SIN*COTAN
          GO TO 50
40      CONTINUE
          TAN = SIGMA(K)/R(K,K)
          COS = P5/DSQRT(P25+P25*TAN**2)
          SIN = COS*TAN
50    CONTINUE
C
C     COMPUTE THE MODIFIED DIAGONAL ELEMENT OF R AND
C     THE MODIFIED ELEMENT OF ((Q TRANSPOSE)*B,0).
C
      R(K,K) = COS*R(K,K) + SIN*SIGMA(K)
      TEMP = COS*WA(K) + SIN*QTBPJ
      QTBPJ = -SIN*WA(K) + COS*QTBPJ
      WA(K) = TEMP
C
C     ACCUMULATE THE TRANFORMATION IN THE ROW OF S.
C
      KP1 = K + 1
      IF (N .LT. KP1) GO TO 70
      DO 60 I = KP1, N
        TEMP = COS*R(I,K) + SIN*SIGMA(I)
        SIGMA(I) = -SIN*R(I,K) + COS*SIGMA(I)
        R(I,K) = TEMP
60    CONTINUE
70    CONTINUE
80    CONTINUE
90    CONTINUE
C
C     STORE THE DIAGONAL ELEMENT OF S AND RESTORE
C     THE CORRESPONDING DIAGONAL ELEMENT OF R.
C
      SIGMA(J) = R(J,J)
      R(J,J) = X(J)
100   CONTINUE
C
C     SOLVE THE TRIANGULAR SYSTEM FOR Z. IF THE SYSTEM IS
C     SINGULAR, THEN OBTAIN A LEAST SQUARES SOLUTION.
C
      NSING = N
      DO 110 J = 1, N
        IF (SIGMA(J) .EQ. ZERO .AND. NSING .EQ. N) NSING = J - 1
        IF (NSING .LT. N) WA(J) = ZERO
110   CONTINUE
      IF (NSING .LT. 1) GO TO 150
      DO 140 K = 1, NSING
        J = NSING - K + 1
        SUM = ZERO
        JP1 = J + 1
        IF (NSING .LT. JP1) GO TO 130
        DO 120 I = JP1, NSING
        SUM = SUM + R(I,J)*WA(I)
120   CONTINUE
130   CONTINUE
        WA(J) = (WA(J) - SUM)/SIGMA(J)
140   CONTINUE
150   CONTINUE
C
C     PERMUTE THE COMPONENTS OF Z BACK TO COMPONENTS OF X.
C
      DO 160 J = 1, N
        L = IPVT(J)
        X(L) = WA(J)
160   CONTINUE
      RETURN
      END
