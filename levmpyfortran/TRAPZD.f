C
C     TRAPEZOIDAL METHOD SUBROUTINE
C       FUNC: function to be integrated
C       A: starting point
C       B: final point
C       S: area
C       N: iteration number
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE TRAPZD(FUNC,A,B,S,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N
      REAL*8 A,B,S
      EXTERNAL FUNC
      INTEGER IT,J
      REAL*8 X
C
      IF (N.EQ.1) THEN
C     First iteration (IT = 1)
        S = 0.5D0*(B-A)*(FUNC(A)+FUNC(B))
      ELSE
C     Subsequent iterations (IT = 2, 4, 8 etc.)
         IT = 2**(N-2)
        TNM = IT
        DEL = (B-A)/TNM       ! spacing of midpoints to be added
          X = A+0.5D0*DEL     ! midpoint
        SUM = 0.D0
        DO 11 J=1,IT
          SUM = SUM+FUNC(X)
            X = X+DEL
11      CONTINUE
         S = 0.5D0*(S + (B-A)*SUM/TNM)
      ENDIF
      RETURN
      END
