C
C     Returns the incomplete gamma function Q(A,X) evaluated by
C       its continued fraction representation as gammcf.
C       Also returns ln(Gamma(a)) as gln.
C
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER ITMAX
      REAL*8 GAMMCF,A,X,GLN
      PARAMETER (ITMAX=32700,EPS=3.E-8,FPMIN=1.D-30)
      INTEGER I
      GLN = GAMMLN(A)
      B = X+1.D0-A
      C = 1.D0/FPMIN
      D = 1.D0/B
      H = D
      DO 11 I=1,ITMAX
        AN = -I*(I-A)
        B = B+2.D0
        D = AN*D+B
        IF(ABS(D).LT.FPMIN) D = FPMIN
        C = B+AN/C
        IF(ABS(C).LT.FPMIN) C = FPMIN
        D = 1.D0/D
        DEL = D*C
        H = H*DEL
        IF(ABS(DEL-1.D0).LT.EPS) GOTO 1
11    CONTINUE
1     GAMMCF = EXP(-X+A*LOG(X)-GLN)*H
      RETURN
      END
