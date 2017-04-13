C
C     Calculates Euclidean Norm from MINPACK
C
C     N: a positive integer input variable
C     X: an input array of length n
C
      DOUBLE PRECISION FUNCTION ENORM(N,X)
      INTEGER N,I
      DOUBLE PRECISION ONE,SCALE,SNORMX,SUM,TEMP,TERM,ZERO,X(N)
      DATA ONE,ZERO /1.0D0,0.0D0/
C
      SNORMX = ZERO
      DO 10 I = 1, N
        SNORMX = DMAX1(SNORMX,DABS(X(I)))
10    CONTINUE
      ENORM = SNORMX
      IF (SNORMX.EQ.ZERO) GOTO 30
      SCALE = SNORMX
      IF (SNORMX.GT.ONE) SCALE = DSQRT(SNORMX)
      SUM = ZERO
      DO 20 I = 1, N
        TERM = ZERO
C
C       THE FOLLOWING TESTS PREVENT UNDERFLOWS IN THE
C       CALCULATION OF TERM AND TERM**2.
C
        TEMP = SCALE + DABS(X(I))
        IF (TEMP.NE.SCALE) TERM = X(I)/SNORMX
        IF (TEMP.NE.ONE) SUM = SUM + TERM**2
        TEMP = ONE + TERM
20    CONTINUE
      ENORM = SNORMX*DSQRT(SUM)
30    CONTINUE
      RETURN
      END
