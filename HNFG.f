C
C   This is the function integrated to find cut-off HN response
C
      DOUBLE PRECISION FUNCTION HNFG(YY)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 YY
      COMMON /CM9/ TOMEGA,PHIX,ICHG,IWT
      COMMON /CM29/ GAM,SN1,CS1,PII,CDN,QPI,MDEX
C
      XEE = DEXP(YY)             ! x = tau/tauo = exp(y)
C
      XGA = XEE**GAM
      THET = DATAN2(SN1,(XGA + CS1))
      SN2 = DSIN(PHIX*THET)
      XGP = XEE**(GAM*PHIX)
      X2G = XGA*XGA
      YEE = DEXP(YY*ICHG)
      XCU = PII*YEE*XGP*SN2/((X2G + 2.D0*XGA*CS1 + 1.D0)**(0.5D0*PHIX))
C
      IF(IWT.EQ.1) THEN
            HNFG = XCU/(1.D0 + (TOMEGA*XEE)**2)
        IF(MDEX.EQ.-6) HNFG = HNFG*XEE
      ELSE 
        XAR = XCU - TOMEGA/XEE
        IF(XAR.LT.-5.D1) THEN
            HNFG = 0.D0
        ELSE
            HNFG = DEXP(XAR)
        ENDIF
      ENDIF 
      RETURN
      END
