C
C   This is the function integrated to find cut-off CD response
C
      DOUBLE PRECISION FUNCTION CDFG(YY)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 YY
      COMMON /CM9/ TOMEGA,PHIX,ICHG,IWT
      COMMON /CM29/ GAM,SN1,CS1,PII,CDN,QPI,MDEX
C
      XEE = 1.D0/(1.D0 + YY**QPI)
      YEE = XEE**ICHG
C
      XCU = CDN*YEE*XEE
C
      IF(IWT.EQ.1) THEN
            CDFG = XCU/(1.D0 + (TOMEGA*XEE)**2)
C
        IF(MDEX.EQ.-6) CDFG = CDFG*XEE
      ELSEIF(IWT.EQ.0.AND.ABS(MDEX).EQ.8) THEN
        TXE = TOMEGA/XEE
        IF(TXE.LE.1.D2) THEN
            CDFG = XCU*DEXP(-TXE)
        ELSE
            CDFG = XCU*DEXP(-1.D2)
        ENDIF
        IF(MDEX.EQ.-8) THEN
            CDFG = CDFG*XEE
        ENDIF
      ELSE 
        XAR = XCU - TOMEGA/XEE
        IF(XAR.LT.-5.D1) THEN
            CDFG = 0.D0
        ELSE
            CDFG = DEXP(XAR)
        ENDIF
      ENDIF 
      RETURN
      END
