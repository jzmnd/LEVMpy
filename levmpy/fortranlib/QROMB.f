C
C     ROMBERG INTEGRATION METHOD USING TRAPEZIUM RULE
C       FUNC: function to be integrated
C       A: starting point
C       B final point
C       VINT: value of integrand
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE QROMB(FUNC,A,B,VINT)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER J,JMAX,JMAXP,K,KM
      REAL*8 A,B,SS,VINT
      PARAMETER (JMAX=22, JMAXP=JMAX+1, K=7, KM=K-1)
      DIMENSION S(JMAXP),H(JMAXP)
      EXTERNAL FUNC
      COMMON /CM10/ EPSG,IZR
C
      H(1)=1.D0
      DO 11 J=1,JMAX
        CALL TRAPZD(FUNC,A,B,S(J),J)
        IF (J.GE.K) THEN
          CALL POLINT(H(J-KM),S(J-KM),K,0.D0,SS,DSS)
          IF(DABS(DSS).LT.EPSG*DABS(SS).OR.SS.EQ.0.D0) THEN
            VINT = SS
            RETURN
          ENDIF
        ENDIF
        S(J+1) = S(J)
        H(J+1) = 0.25D0*H(J)
11    CONTINUE
      WRITE(*,*) 'TOO MANY STEPS IN QROMB'
      VINT = SS
      END
