C
C     Extended midpoint rule subroutine
C       FUNC: function to be integrated
C       A: starting point
C       B: final point
C       S: area
C       N: iteration number
C
      SUBROUTINE MIDPNT(FUNC,A,B,S,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N
      REAL*8 A,B,S
      EXTERNAL FUNC
      INTEGER IT,J
      REAL*8 X,XX
C
      IF (N.EQ.1) THEN
C     First iteration (IT = 1)
        XX = 0.5D0*(A+B)
        S = (B-A)*FUNC(XX)
      ELSE
C     Subsequent iterations
          IT = 3**(N-2)
         TNM = IT
         DEL = (B-A)/(3.D0*TNM)
        DDEL = DEL+DEL
           X = A+0.5D0*DEL
         SUM = 0.D0
        DO 11 J=1,IT
          SUM = SUM+FUNC(X)
            X = X+DDEL
          SUM = SUM+FUNC(X)
            X = X+DEL
11      CONTINUE
        S = (S+(B-A)*SUM/TNM)/3.D0
      ENDIF
      RETURN
      END
