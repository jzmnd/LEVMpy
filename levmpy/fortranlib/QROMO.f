C
C     ROMBERG INTEGRATION METHOD USING MIDPOINT RULE ON OPEN INTEGRATION
C       FUNC: function to be integrated
C       A: starting point
C       B final point
C       SS: value of integrand
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE QROMO(FUNC,A,B,SS)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER J,JMAX,JMAXP,K,KM
      REAL*8 A,B,SS
      PARAMETER (JMAX=22, JMAXP=JMAX+1, K=7, KM=K-1)
      DIMENSION S(JMAXP),H(JMAXP)
      EXTERNAL FUNC
      COMMON /CM10/ EPSG,IZR
C 
      H(1)=1.D0
      DO 11 J=1,JMAX
        CALL MIDPNT(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          CALL POLINT(H(J-KM),S(J-KM),K,0.D0,SS,DSS)
          SERR = EPSG*ABS(SS)
          IF(J.GT.10) WRITE(*,150) J,EPSG,DSS,SERR
C
150   FORMAT(I4,1P,(3D14.5))
          IF (ABS(DSS).LE.SERR) RETURN
        ENDIF
        S(J+1) = S(J)
        H(J+1) = H(J)/9.D0
11    CONTINUE
      WRITE(*,*) 'TOO MANY STEPS IN QROMO'
      END
