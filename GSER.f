C
C     Returns the incomplete gamma function P(A,X) evaluated by
C       its series representation as gamser.
C       Also returns ln(Gamma(a)) as gln.
C
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER ITMAX
      REAL*8 GAMSER,A,X,GLN
      PARAMETER (ITMAX=100,EPS=3.E-8)
      INTEGER N
      GLN = GAMMLN(A)
      IF(X.LE.0.D0) THEN
        IF(X.LT.0.D0) THEN
          WRITE(*,*) 'X IS LESS THAN ZERO IN GSER'
          RETURN
C          STOP
        ENDIF
        GAMSER = 0.D0
        RETURN
      ENDIF
      AP = A
      SUM = 1.D0/A
      DEL = SUM
      DO 11 N=1,ITMAX
        AP = AP+1.D0
        DEL = DEL*X/AP
        SUM = SUM+DEL
        IF(DABS(DEL).LT.DABS(SUM)*EPS) GOTO 1
11    CONTINUE
1     GAMSER = SUM*DEXP(-X+A*DLOG(X)-GLN)
      RETURN
      END
