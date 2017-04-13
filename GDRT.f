C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      DOUBLE PRECISION FUNCTION GDRT(YY)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 YY
      COMMON /CM5/ TOMEGA,PHICOM,IV
      COMMON /CM9/ TCMEGA,XIGDAE,ICHG,IWT
      COMMON /CM12/ CELCAP,ATEMP,WF,MAXFEV,ICF,MDE,JCDX
C   HERE YY IS LN(TAU)
      IF(YY.GT.1.D2) THEN
          WRITE(*,*) YY
          YY = 1.D2
      ENDIF       
      XX = DEXP(YY)              !THIS IS TAU/TAUO
      XXI = XX**IV
C
      CALL KWWDRT(XX,PHICOM,DRTW,ATEMP)
C
      IF(IWT.EQ.0) THEN
        IF(MDE.lt.0) DRTW = DRTW*XX
      ELSE
        DRTW = DRTW*(XX**ICHG)
      ENDIF
C
C   TRANSIENT RESPONSE
682   IF(ABS(MDE).EQ.8) THEN
        GDRT = DRTW*DEXP(-TOMEGA/XX)
      ELSE
C   FREQUENCY RESPONSE
C
      GDRT = XXI*DRTW/(1.D0 + (TOMEGA*XX)**2)
C
      ENDIF
      RETURN
      END
