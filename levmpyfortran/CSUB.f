      SUBROUTINE CSUB(M,FREQ,P,F)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'SIZE.INC'
      INTEGER M
      DOUBLE PRECISION P(NTOT),F(NPT2),FREQ(NPT2),L
      COMPLEX*16 YC,ZA,ZB,ZC,ZT,ZD,IOMEGA
      COMMON /CM47/ ICNT
      LOGICAL NDE20,C20,RC30,RC50,RP0,RP1
C              SAVE LOGICAL CHECKS FOR SUBSEQUENT CALLS
      SAVE NDE20,C20,RC30,RC50,RP0,RP1
C
C   ******************   C CIRCUIT:
C
C   THIS CIRCUIT HAS TWO AUGMENTED DISTRIBUTED ELEMENTS AND ONE
C   REGULAR DISTRIBUTED ELEMENT.
C
C         M : number of data points (IN)
C      FREQ : array of frequency values (IN)
C         P : array of model parameters (IN)
C         F : model function values (OUT)
C
C   SET PARAMETER VALUES
C
      R1 = P(1)
      R2 = P(2)
      C2 = P(3)
      R3 = P(4)
      C3 = P(5)
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
      RA1 = P(21)
      CA1 = P(22)
      RA2 = P(23)
      CA2 = P(24)
      R4 = P(25)
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
         C20 = (C2.EQ.0)
        RC30 = ((R3.EQ.0.D0).AND.(C3.GT.0.D0))
        RC50 = ((R5.EQ.0.D0).AND.(C5.GT.0.D0))
         RP0 = ((RP.EQ.0.D0).AND.(CP.NE.0.D0))
         RP1 = (RP.EQ.0)
        NDE20 = (NDE2.EQ.0)
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
        CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE3,TDE3,UDE3,PDE3,NDE3,ZD)
        ZT = ZT + R4 + ZD
C
        CALL SDEA(OMEGA,IOMEGA,RA2,CA2,RDE2,TDE2,UDE2,PDE2,NDE2,ZD)
        IF(ZT.EQ.(0,0)) THEN
          ZT = ZD
        ELSEIF(NDE20) THEN
          ZT = ZT
        ELSE
          ZT = ZT*ZD/(ZT + ZD)
        ENDIF   
C   
        IF(RC30) THEN
          ZA = 1.D0/(C3*IOMEGA)
        ELSE
          ZA = R3/(C3*R3*IOMEGA + 1.D0)
        ENDIF
        CALL SDEA(OMEGA,IOMEGA,RA1,CA1,RDE1,TDE1,UDE1,PDE1,NDE1,ZD)
        ZA = ZA + ZD
C
        IF(C20) THEN
          ZB = R2
        ELSE
          ZB = R2 + 1.D0/(IOMEGA*C2)
        ENDIF   
C
        IF(ZB.EQ.(0,0)) THEN
          ZC = ZA
        ELSEIF(ZA.EQ.(0,0)) THEN
          ZC = ZB
        ELSE
          ZC = ZA*ZB/(ZA + ZB)
        ENDIF
        ZT = R1 + ZC + ZT
C
C    CALCULATE EFFECTS OF RP AND CP AND/OR L IF PRESENT
C
        IF (RP0) THEN
          YC = CP*IOMEGA
        ELSE
        IF(RP1) THEN
          YC = (0,0)
        ELSE
          YC = (CP*RP*IOMEGA + 1.D0)/RP
        ENDIF
        END IF
        IF(ZT.EQ.(0,0).AND.YC.NE.(0,0)) THEN
          ZT = 1.0/YC
        ELSE
          ZT = (ZT/(1.D0 + ZT*YC))
        ENDIF
        ZT = ZT + L*IOMEGA
C
        F(I) = DREAL(ZT)
        F(I+M) = DIMAG(ZT)
100   CONTINUE
      RETURN
      END
