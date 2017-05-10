C
C   THIS IS THE FUNCTION INTEGRATED TO FIND GENERAL EXPONENTIAL RESPONSE
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      DOUBLE PRECISION FUNCTION GDAEFN1(X)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION X
      COMMON /CM9/ TOMEGA,PHIX,ICHG,IWT
      COMMON /CM29/ GAM,SN1,CS1,PII,CDN,QPI,MDEX
C
      SGNP = DSIGN(1.D0,PHIX)
      XABX = DABS(X*PHIX)
C                   Here, x is actually usual y
      XEE = DEXP(X)              ! y = ln(x); x=tau/tauo
C                          
      XCU = X*ICHG - SGNP*(XABX**GAM)/GAM
C
      IF(DABS(XCU).GT.3.D2) XCU = DSIGN(3.D2,XCU)
      IF(IWT.EQ.1) THEN
            GDAEFN1 = DEXP(XCU)/(1.D0 + (TOMEGA*XEE)**2)
        IF(MDEX.LT.0) GDAEFN1 = GDAEFN1*XEE
      ELSE 
        XAR = XCU - TOMEGA*DEXP(-X)
        IF(XAR.LT.-5.D1) THEN
            GDAEFN1 = 0.D0
        ELSE
            GDAEFN1 = DEXP(XAR)
        ENDIF
      ENDIF 
      RETURN
      END
