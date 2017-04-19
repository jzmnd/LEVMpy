C
C     POLYNOMIAL INTERPOLATION OR EXTRAPOLATION SUBROUTINE
C     The unique polynomial interpolation for N points of degree N-1
C       XA: array of x values
C       YA: array of y values
C       N: length of arrays
C       X: value of x to interpolate to
C       Y: interpolated y value at x
C       DY: error estimate in Y
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N
      REAL*8 XA,YA,X,Y,DY
      PARAMETER (NMAX=10)
      DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
C
      NS = 1
      DIF = DABS(X-XA(1))
      DO 11 I=1,N
        DIFT = DABS(X-XA(I))
        IF (DIFT.LT.DIF) THEN
           NS = 1
          DIF = DIFT
        ENDIF
          C(I) = YA(I)
          D(I) = YA(I)
11    CONTINUE
       Y = YA(NS)
      NS = NS-1
      DO 13 M=1,N-1
        DO 12 I=1,N-M
          HO = XA(I)-X
          HP = XA(I+M)-X
           W = C(I+1)-D(I)
          DEN = HO-HP
          IF(DEN.EQ.0.D0) THEN
            WRITE(*,*) 'DEN IS ZERO'
            RETURN
C            STOP
          ENDIF
           DEN = W/DEN
          D(I) = HP*DEN
          C(I) = HO*DEN
12      CONTINUE
        IF(2*NS.LT.N-M) THEN
          DY = C(NS+1)
        ELSE
          DY = D(NS)
          NS = NS-1
        ENDIF
        Y = Y+DY
13    CONTINUE
      RETURN
      END
