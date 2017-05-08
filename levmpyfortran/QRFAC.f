C
C     LEVENBERG–MARQUARDT ALGORITHM FROM MINPACK
C
C       M: a positive integer input variable set to the number
C          of rows of a.
C
C       N: a positive integer input variable set to the number
C          of columns of a.
C
C       A: an m by n array. on input a contains the matrix for
C          which the qr factorization is to be computed. on output
C          the strict upper trapezoidal part of a contains the strict
C          upper trapezoidal part of r, and the lower trapezoidal
C          part of a contains a factored form of q (the non-trivial
C          elements of the u vectors described above).
C
C       LDA: a positive integer input variable not less than m
C            which specifies the leading dimension of the array a.
C
C       PIVOT: a logical input variable. if pivot is set true,
C              then column pivoting is enforced. if pivot is set false,
C              then no column pivoting is done.
C
C       IPVT: an integer output array of length lipvt. ipvt
C             defines the permutation matrix p such that a*p = q*r.
C             column j of p is column ipvt(j) of the identity matrix.
C             if pivot is false, ipvt is not referenced.
C
C       LIPVT: a positive integer input variable. if pivot is false,
C              then lipvt may be as small as 1. if pivot is true, then
C              lipvt must be at least n.
C
C       RDIAG: an output array of length n which contains the
C              diagonal elements of r.
C
C       ACNORM: an output array of length n which contains the
C         norms of the corresponding columns of the input matrix a.
C         if this information is not needed, then acnorm can coincide
C         with rdiag.
C
C       WA is a work array of length n.
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE QRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,SIGMA,ACNORM,WA)
      INTEGER M,N,LDA,LIPVT
      INTEGER IPVT(LIPVT)
      LOGICAL PIVOT
      DOUBLE PRECISION A(LDA,N),SIGMA(N),ACNORM(N),WA(N)
      INTEGER I,J,JP1,K,KMAX,MINMN
      DOUBLE PRECISION AJNORM,EPSMCH,ONE,P05,SUM,TEMP,ZERO
      DOUBLE PRECISION DPMPAR,ENORM
      DATA ONE,P05,ZERO /1.0D0,5.0D-2,0.0D0/
C
C     EPSMCH IS THE MACHINE PRECISION.
C
      EPSMCH = DPMPAR(1)
C
C     COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS.
C
      DO 10 J = 1, N
        ACNORM(J) = ENORM(M,A(1,J))
        SIGMA(J) = ACNORM(J)
        WA(J) = SIGMA(J)
        IF (PIVOT) IPVT(J) = J
10    CONTINUE
C
C     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS.
C
      MINMN = MIN0(M,N)
      DO 110 J = 1, MINMN
        IF (.NOT.PIVOT) GO TO 40
C
C        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION.
C
        KMAX = J
        DO 20 K = J, N
          IF (SIGMA(K) .GT. SIGMA(KMAX)) KMAX = K
20      CONTINUE
        IF (KMAX .EQ. J) GO TO 40
          DO 30 I = 1, M
            TEMP = A(I,J)
            A(I,J) = A(I,KMAX)
            A(I,KMAX) = TEMP
30        CONTINUE
          SIGMA(KMAX) = SIGMA(J)
          WA(KMAX) = WA(J)
          K = IPVT(J)
          IPVT(J) = IPVT(KMAX)
          IPVT(KMAX) = K
40      CONTINUE
C
C     COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE
C     J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR.
C
        AJNORM = ENORM(M-J+1,A(J,J))
        IF (AJNORM .EQ. ZERO) GO TO 100
        IF (A(J,J) .LT. ZERO) AJNORM = -AJNORM
        DO 50 I = J, M
          A(I,J) = A(I,J)/AJNORM
50      CONTINUE
        A(J,J) = A(J,J) + ONE
C
C     APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS
C     AND UPDATE THE NORMS.
C
        JP1 = J + 1
        IF (N .LT. JP1) GO TO 100
        DO 90 K = JP1, N
          SUM = ZERO
          DO 60 I = J, M
            SUM = SUM + A(I,J)*A(I,K)
60        CONTINUE
          TEMP = SUM/A(J,J)
          DO 70 I = J, M
            A(I,K) = A(I,K) - TEMP*A(I,J)
70        CONTINUE
          IF (.NOT.PIVOT .OR. SIGMA(K) .EQ. ZERO) GO TO 80
          TEMP = A(J,K)/SIGMA(K)
          SIGMA(K) = SIGMA(K)*DSQRT(DMAX1(ZERO,ONE-TEMP**2))
          IF (P05*(SIGMA(K)/WA(K))**2 .GT. EPSMCH) GO TO 80
          SIGMA(K) = ENORM(M-J,A(JP1,K))
          WA(K) = SIGMA(K)
80      CONTINUE
90      CONTINUE
100     CONTINUE
        SIGMA(J) = -AJNORM
110   CONTINUE
      RETURN
      END
