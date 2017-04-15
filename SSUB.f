      SUBROUTINE SSUB(M,FREQ,P,F,NFREE)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'SIZE.INC'
      INTEGER M,NFREE
      DOUBLE PRECISION P(NTOT),F(2*M),FREQ(M),LE
      COMPLEX*16 ZT,Z1,YDL,ZE,YSD,IOMEGA
C
C   **********************     S CIRCUIT:
C
C   CALCULATES IMPEDANCE OF ELECTROCHEMICAL SUPERCAPACITOR.
C
C         M : number of data points (IN)
C      FREQ : array of frequency values (IN)
C         P : array of model parameters (IN)
C         F : model function values (OUT)
C     NFREE : free parameter array (IN)
C
C   SET PARAMETER VALUES
C
       RB = P(1)
       RE = P(2)
       RF = P(3)
       CF = P(4)
      RSD = P(5)
      CDL = P(6)
       LE = P(7)
C
C   LOOP OVER ALL FREQUENCIES
C
      DO 100 I=1,M
        OMEGA = FREQ(I)
        IOMEGA = DCMPLX(0.D0,OMEGA)
C
C   CALCULATE IMPEDANCE
C
        ZT = 1.D0/(CF*IOMEGA)
        ZT = RB + ZT
        ZT = 1.D0/ZT
        YSD = 1.D0/RSD
        ZT = ZT + YSD
        YDL = CDL*IOMEGA
        Z1 = ZT + YDL
        ZE = LE*IOMEGA
        ZE = RE + ZE
        ZT = ZE/Z1
        ZT = CDSQRT(ZT)
        Z1 = ZE*Z1
        Z1 = CDSQRT(Z1)
        Z1 = (CDEXP(Z1)+CDEXP(-Z1))/(CDEXP(Z1)-CDEXP(-Z1))
        ZT = Z1*ZT
C
C   RETURN IMPEDANCE VALUES
C
        F(I) = DREAL(ZT)
        F(I+M) = DIMAG(ZT)
100   CONTINUE
      RETURN
      END
