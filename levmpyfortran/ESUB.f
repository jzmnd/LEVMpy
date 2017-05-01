      SUBROUTINE ESUB(M,FREQ,P,F)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'SIZE.INC'
      INTEGER M
      DOUBLE PRECISION P(NTOT),F(NPT2),FREQ(NPT2),L
      COMPLEX*16 ZA,ZB,ZC,ZD,ZT,IOMEGA
      COMMON /CM47/ ICNT    
      LOGICAL RD30,CP0,DE50,DE40
C
C   *********************   E CIRCUIT:
C
C   THIS SUBROUTINE CALCULATES THE IMPEDANCE OF A SET OF NESTED
C   DISTRIBUTED ELEMENT CIRCUITS.  FOR MORE DETAILED INFO,
C   CONSULT THE CIRCUIT MODEL DOCUMENTATION.
C
C         M : number of data points (IN)
C      FREQ : array of frequency values (IN)
C         P : array of model parameters (IN)
C         F : model function values (OUT)
C
C   SET PARAMETER VALUES
C
      RDE1 = P(1)
      TDE1 = P(2)
      UDE1 = P(3)
      PDE1 = P(4)
      NDE1 = IDINT(P(5))
C
      RDE2 = P(6)
      TDE2 = P(7)
      UDE2 = P(8)
      PDE2 = P(9)
      NDE2 = IDINT(P(10))
C
      RDE3 = P(11)
      TDE3 = P(12)
      UDE3 = P(13)
      PDE3 = P(14)
      NDE3 = IDINT(P(15))
C
      RDE4 = P(16)
      TDE4 = P(17)
      UDE4 = P(18)
      PDE4 = P(19)
      NDE4 = IDINT(P(20))
C
      RDE5 = P(21)
      TDE5 = P(22)
      UDE5 = P(23)
      PDE5 = P(24)
      NDE5 = IDINT(P(25))
C
      R1 = P(26)
      R2 = P(27)
      R3 = P(28)
      CP = P(29)
       L = P(30)
C
C   SET VARIABLES FOR ZERO CHECKING
C
      IF(ICNT.LE.1) THEN
         CP0 = (CP.EQ.0)   
        RD30 = ((R3.EQ.0.D0).AND.(NDE3.EQ.0))
        DE40 = (NDE4.EQ.0)
        DE50 = (NDE5.EQ.0)
      ENDIF
C
C   LOOP OVER ALL FREQUENCIES
C
      DO 100 I=1, M
        OMEGA = FREQ(I)
        IOMEGA = DCMPLX(0.D0,OMEGA)
C
C   CALCULATE IMPEDANCE OF CIRCUIT
C
        CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE5,TDE5,UDE5,PDE5,NDE5,ZD)
        ZA = ZD
        CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE3,TDE3,UDE3,PDE3,NDE3,ZD)
        ZB = ZD + R3
        IF(RD30) THEN
          ZT = ZA
        ELSEIF(DE50) THEN
          ZT = ZB
        ELSE 
          ZT = ZA*ZB/(ZA + ZB)
        ENDIF
C
        CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE2,TDE2,UDE2,PDE2,NDE2,ZD)
          ZT = ZT + ZD + R2
        CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE4,TDE4,UDE4,PDE4,NDE4,ZD)
        IF(DE40) THEN
          ZA = ZT
        ELSEIF(ZT.EQ.(0,0)) THEN
          ZA = ZD
        ELSE
          ZA = ZD*ZT/(ZD + ZT)
        ENDIF
C
        CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE1,TDE1,UDE1,PDE1,NDE1,ZD)
        ZB = ZD + R1 + ZA
C
        IF(CP0) THEN
            ZC = ZB
        ELSEIF(ZB.EQ.(0,0)) THEN
            ZC = 1.D0/(IOMEGA*CP)
        ELSEIF(ZB.NE.(0.D0,0.D0)) THEN
            ZC = ZB/(1.D0 + IOMEGA*CP*ZB)
        ENDIF
C
        ZT = L*IOMEGA + ZC
C
C   RETURN IMPEDANCE VALUES
C
        F(I) = DREAL(ZT)
        F(I+M) = DIMAG(ZT)
100   CONTINUE
      RETURN
      END
