      SUBROUTINE HYPER(XW,LL,S4,S5,AA,BB,CC,DD,D6,IC)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER LL,IC
      REAL*8 XW,S4,S5,AA,BB,CC,DD,D6
      DIMENSION SS1(38),SS2(38)
      DATA NMAX/11/,EPS/1.D-8/,WLL/0.69D0/
C
      ID = IC
1330  D6 = D6*AA*BB*XW/(CC*DD)
      IF(LL.EQ.0) GOTO 1380
      S5 = S5 + D6
      IC = IC + 1
      IF(IC.LE.NMAX) SS1(IC) = S5
      LL = 0
      GOTO 1410
1380   D6 = -D6
      S4 = S4 + D6
      ID = ID + 1
      IF(ID.LE.NMAX) SS2(ID) = S4
      LL = 1
1410  IF(DABS(D6).LT.EPS) GOTO 1480
C
      AA = AA + 1.D0
      BB = BB + 1.D0
      CC = CC + 1.D0
      DD = DD + 1.D0
      IF(IC.LE.NMAX.OR.DABS(XW).LT.WLL) GOTO 1330
C
      IF(IC.LT.NMAX) THEN
        RETURN
      ENDIF
      CALL EPSALG(NMAX,SS1,SF)
        S5 = SF
      CALL EPSALG(NMAX,SS2,SF)
        S4 = SF
      GOTO 1480
1480  CONTINUE
      RETURN
      END
