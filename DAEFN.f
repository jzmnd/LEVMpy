C
C   THIS IS THE FUNCTION INTEGRATED TO FIND THE EDAE, EDAE1 RESPONSE
C     Exponential distribution of activation energies
C
      DOUBLE PRECISION FUNCTION DAEFN(X)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 X
      COMMON /CM5/ TOMEGA,PHICOM,IV
C
      XPU = X*PHICOM
      IF(DABS(XPU).GT.3.D2) XPU = DSIGN(3.D2,XPU)
      XVU = IV*X
      IF(DABS(XVU).GT.3.D2) XVU = DSIGN(3.D2,XVU)
      DAEFN = DEXP(-XPU)/(1.D0 + (TOMEGA*DEXP(-XVU))**2)
      RETURN
      END
