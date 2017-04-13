C
C   This is the function integrated to find general exponential response
C     SEE NELEM=34
C
      DOUBLE PRECISION FUNCTION GDAEFN2(X)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'SIZE.INC'
      REAL*8 X
      CHARACTER*1 DATTYQ
      COMMON /CM9/ TOMEGA,XIGDAE,ICHG,IWT
      COMMON /CM2/ Y(NPT2),R(NPT2),FJ(NPT2),P(NTOT),DRSS,ROE,RKE,
     + NS(NPAFR),NFREE(NTOT),NP,ICNT,MN,IRCH,IXI,DATTYQ
      COMMON /CM12/ CELCAP,ATEMP,WF,MAXFEV,ICF,MDE,JCDX
C
      XAB = DABS(X)
      XEE = DEXP(X)
      IF(ABS(MDE).EQ.8) THEN          !TRANSIENT
          IWT = 0
      ELSE
          IWT = 1         !FREQUENCY
      ENDIF
C
      XCU = X*ICHG - (XAB*XIGDAE)**P(13)
        IF(DABS(XCU).GT.3.D2) XCU = DSIGN(3.D2,XCU)
      IF(IWT.EQ.1) THEN
            GDAEFN2 = DEXP(XCU)/(1.D0 + (TOMEGA*XEE)**2)
      ELSE 
            GDAEFN2 = DEXP(XCU - TOMEGA/XEE)
      ENDIF 
        IF(MDE.LT.0) GDAEFN2 = XEE*GDAEFN2
      RETURN
      END
