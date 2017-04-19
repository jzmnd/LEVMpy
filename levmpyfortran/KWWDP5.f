C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE KWWDP5(TOMEGA,UU,SUM)
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL GAMMQ
      REAL*8 TOMEGA,UU,SUM
      DATA PISQR/1.77245385090552D0/
C
      XUU = DEXP(UU)
      IF(TOMEGA.LE.XUU) THEN
        SUM = 0.D0
        GOTO 529
      ENDIF
        TOMEG4 = TOMEGA/4.D0
        PINV = GAMMQ(0.5D0,XUU/4.D0)
        PDR = DSQRT(TOMEGA)*DEXP(-TOMEG4)
        SUM = PDR/(2.D0*PISQR*PINV)
C
529   RETURN
      END