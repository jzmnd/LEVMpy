      SUBROUTINE JSUB(M,FREQ,P,F)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'SIZE.INC'
      INTEGER M
      DOUBLE PRECISION P(NTOT),F(NPT2),FREQ(NPT2),L
      COMPLEX*16 YA,YC,YD,ZT,ZD,IOMEGA
      COMMON /CM47/ ICNT    
      LOGICAL NDE1N,NDE30,RC50,RP0,RP1
C              SAVE LOGICAL CHECKS FOR SUBSEQUENT CALLS
      SAVE NDE1N,NDE30,RC50,RP0,RP1
C
C   THIS CIRCUIT HAS ONE AUGMENTED DISTRIBUTED ELEMENTS AND TWO
C   REGULAR DISTRIBUTED ELEMENTS.
C
C         M : number of data points (IN)
C      FREQ : array of frequency values (IN)
C         P : array of model parameters (IN)
C         F : model function values (OUT)
C
C   SET PARAMETER VALUES
C
      R1 = P(1)
      RA = P(2)
      CA = P(3)
      R4 = P(4)
      RL = P(5)
C
      RDE1 = P(6)
      TDE1 = P(7)
      UDE1 = P(8)
      PDE1 = P(9)
      NDE1 = IDINT(P(10))
C
      RDE2 = P(11)
      TDE2 = P(12)
      UDE2 = P(13)
      PDE2 = P(14)
      NDE2 = IDINT(P(15))
C
      RDE3 = P(16)
      TDE3 = P(17)
      UDE3 = P(18)
      PDE3 = P(19)
      NDE3 = IDINT(P(20))
C
      RDE4 = P(21)
      TDE4 = P(22)
      UDE4 = P(23)
      PDE4 = P(24)
      NDE4 = IDINT(P(25))
C
      R5 = P(26)
      C5 = P(27)
      RP = P(28)
      CP = P(29)
       L = P(30)
C
C   SET VARIABLES FOR ZERO-CHECKING
C
      IF(ICNT.LE.1) THEN  
         RC50 = ((R5.EQ.0.D0).AND.(C5.GT.0.D0))
          RP0 = ((RP.EQ.0.D0).AND.(CP.NE.0.D0))
          RP1 = (RP.EQ.0)
        NDE1N = (NDE1.NE.0)
        NDE30 = (NDE3.EQ.0)
      ENDIF
C
C   LOOP OVER ALL FREQUENCIES
C
      DO 100 I=1,M
        OMEGA = FREQ(I)
        IOMEGA = DCMPLX(0.D0,OMEGA)
C
C   EVALUATE IMPEDANCE
C
        IF (RC50) THEN
          ZT = 1.D0/(C5*IOMEGA)
        ELSE
          ZT = R5/(C5*R5*IOMEGA + 1.D0)
        ENDIF
        CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE4,TDE4,UDE4,PDE4,NDE4,ZD)
        ZT = ZT + R4 + ZD
        CALL SDEA(OMEGA,IOMEGA,RA,CA,RDE3,TDE3,UDE3,PDE3,NDE3,ZD)
        IF(ZT.EQ.(0,0)) THEN
          ZT = ZD
        ELSEIF(NDE30) THEN
          ZT = ZT
        ELSE
          ZT = ZT*ZD/(ZT + ZD)
        ENDIF   
C   
        CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE2,TDE2,UDE2,PDE2,NDE2,ZD)
        ZT = ZT + ZD + R1
C
        CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE1,TDE1,UDE1,PDE1,NDE1,ZD)
C
C    CALCULATE EFFECTS OF RP AND CP AND/OR L AND RL IF PRESENT
C
        IF (RP0) THEN
          YC = CP*IOMEGA
        ELSE
        IF(RP1) THEN
          YC = (0,0)
        ELSE
          YC = (CP*RP*IOMEGA + 1.D0)/RP
        ENDIF
        ENDIF
C
        IF(NDE1N) THEN
            YD = 1.D0/ZD
        ELSE
            YD = 0.D0
        ENDIF       
        YA = YC + YD
C
        IF(ZT.EQ.(0,0).AND.YA.NE.(0,0)) THEN
            ZT = 1.0/YA
        ELSE
            ZT = (ZT/(1.D0 + ZT*YA))
        ENDIF
        ZT = ZT + L*IOMEGA + RL
        IF(I.EQ.M) THEN
C            OPEN(69,FILE='LNL')
C            WRITE(69,879) 1,I,YC,YA,ZT
            WRITE(*,879) 1,I,YC,YA,ZT
C            WRITE(69,*) I,OMEGA,ZT
C            CLOSE(69)
        ENDIF
879     FORMAT(1X,2I3,1P,(4E5.1))
        F(I) = DREAL(ZT)
        F(I+M) = DIMAG(ZT)
100   CONTINUE
      RETURN
      END
