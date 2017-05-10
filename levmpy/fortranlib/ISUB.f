      SUBROUTINE ISUB(M,FREQ,P,F)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'SIZE.INC'
      INTEGER M
      DIMENSION AA(2,2),PP(2),QQ(2),TRSY(2)
      DOUBLE PRECISION P(NTOT),F(NPT2),FREQ(NPT2)
C
C     THIS SUBROUTINE IS FOR NLS FITTING TO Y = F(T) OR F(OMEGA)
C           WHERE Y IS REAL INPUT AND T (TIME)  IS OMEGA INPUT
C                IT CALCULATES EDAE TRANSIENT RESPONSE
C                     RUN ONLY WITH DATTYP=R CHOICE
C
C         M : number of data points (IN)
C      FREQ : array of frequency values (IN)
C         P : array of model parameters (IN)
C         F : model function values (OUT)
C
C   SET PARAMETER VALUES
C
      AY = P(1)
      BY = P(2)
      CY = P(3)
      DY = P(4)
      EY = P(5)
      FY = P(6)
C
      BYI = 1.D0/BY
C   LOOP OVER ALL FREQUENCIES (here Times)
      DO 100 I=1,M
        OMEGA = FREQ( I)
        YY = BYI*OMEGA
C
        PP(1) = 1
        PP(2) = EY
        QQ(1) = FY
        QQ(2) = 1
C   FIRST INDEX: XX ARGUMENTS; SECOND INDEX: AA AND BB PHI ARGUMENTS
        AA(1,1) = 1.D0 - CY
        AA(2,1) = 1.D0 - DY
        AA(1,2) = 1.D0 + AA(1,1)
        AA(2,2) = 1.D0 + AA(2,1)
C
C   LOOP FOR AA AND BB TRANSIENT RESPONSE
C   FIRST CALCULATE Y = 0 LIMITING VALUE
        IL = 1
        CALL TRAGA(AA,YY,PP,QQ,IL,TRSY)
C
        F(I) = AY*TRSY(1)
        F(I+M) = 0.D0
100   CONTINUE
      RETURN
      END
