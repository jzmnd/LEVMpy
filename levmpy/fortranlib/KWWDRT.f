C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE KWWDRT(TAU,PHI,SUM,ATEMP)
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL GAMMLN
      REAL*8 TAU,PHI,SUM,ATEMP
      DIMENSION XLH(133)
      DATA PIIN/0.318309886183791D0/,PPI/3.141592653589793D0/,
     +    TPIIN/0.398942280D0/,P24I/.04166666666666667D0/,P1152I/
     +    8.6805555555555556D-4/
      SAVE ISAVE,XLH,BTX,BETA,NEPS,ONEMP,PHDIMP,PHQU,ST1,ST2,ST0,ITW,
     + SACC
C
C   CALCULATE MOMENTS AND PHI QUANTITIES FOR FINAL OUTPUT
C
      IF(ISAVE.EQ.0.OR.PHI.NE.BETA) THEN
C
        NEPS = 133
        PHIM = 1.D0/PHI
        ONEMP = 1.D0 - PHI
        ONEMI = 1.D0/ONEMP
        PHPRI = PHIM*ONEMI
        PHP5 = DSQRT(PHPRI)
C        PHPI2 = PHPRI*PHPRI
        TWMP = 2.D0 - PHI
        ONM2P = 1.D0 - 2.D0*PHI
        PHDIMP = PHI*ONEMI
        PHQU = 0.5D0*TWMP*ONEMI
        PPOL = 2.D0 + PHI*(19.0D0 + 2.D0*PHI)
        SS1 = TWMP*ONM2P*PHPRI
        SS2 = SS1*PPOL*PHPRI
        ST1 = -P24I*SS1
        ST2 = P1152I*SS2
        ST0 = TPIIN*PHP5
C
        DO IJ = 1,NEPS
          AGNL = GAMMLN(DFLOAT(IJ + 1))
          AGNH = GAMMLN(1.D0 + DFLOAT(IJ)*PHI)
          XLH(IJ) = DEXP(AGNH - AGNL)
        ENDDO
          ISAVE = 1
          ITW = IDINT(DABS(ATEMP) + 1.D-3)
          SACC = 1.D-12
C
C   CALCULATE TAU SWITCH POINT: SERIES TO FORMULA
        BTX = 1.D1**(-0.85533D0 + PHIM)
        BETA = PHI
      ENDIF
C
C   *********************************************
C
C   CALCULATE DRT FROM SERIES, APPROX DRT
C
      IF(ITW.GT.200.AND.(TAU.GT.BTX)) GOTO 693 
      ALP = 0.D0
      IX = 0
C
C   *********************************************
C
      DO KJ = 1,NEPS
        PHK = PHI*KJ
        SK = DSIN(PPI*PHK)
        IF(DABS(SK).GT.1.D-10) THEN 
          IX = IX + 1
          IXO = KJ
          ICN = IXO - 1
          PHN = PHI*IXO
          ISGR = (-1)**ICN
          SKN = DSIN(PPI*PHN)
          DELP = XLH(IXO)*SKN*ISGR*(TAU**PHN)
          ALP = ALP + DELP
          ARDD = DABS(DELP/ALP)
          IF(ARDD.LT.SACC) GOTO 348
        ENDIF
      ENDDO   
C
348   CONTINUE
      SUM = PIIN*ALP
      IF(SUM.LE.0.D0.OR.SUM.GT.2.D0) THEN
        SUM = 0.D0
      ENDIF
      SUMF = SUM
693   CONTINUE
C
      IF(TAU.GT.BTX) THEN
        XX = PHI*TAU
        SUMA = ST0*(XX**PHQU)*DEXP(-ONEMP*(XX**PHDIMP))*(1.D0 + ST1*
     +(XX**(-PHDIMP)) + ST2*(XX**(-2.D0*PHDIMP)))/TAU
        SUM = SUMA
      ENDIF
C
      IF(ITW.EQ.200.OR.ITW.EQ.400) THEN
        WRITE(*,*) 33,IX,TAU,SUM
      ENDIF
      RETURN
      END
