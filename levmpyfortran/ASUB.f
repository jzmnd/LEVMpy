      SUBROUTINE ASUB(M,FREQ,P,F)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'SIZE.INC'
      INTEGER M
      DOUBLE PRECISION P(NTOT),F(2*M),FREQ(M),L
      COMPLEX*16 ZT,ZD,Z5,YC,IOMEGA
      COMMON /CM47/ ICNT    
      LOGICAL RC10,RC20,RC30,RC40,R40,RP0,RP1
C
C   **********************     A CIRCUIT:
C
C   CALCULATES IMPEDANCE OF 6 SUBCIRCUITS IN SERIES AND
C   IN PARALLEL WITH RP AND CP AND IN SERIES WITH L.
C
C         M : number of data points (IN)
C      FREQ : array of frequency values (IN)
C         P : array of model parameters (IN)
C         F : model function values (OUT)
C
C   SET PARAMETER VALUES
C
      R1 = P(1)
      C1 = P(2)
      R2 = P(3)
      C2 = P(4)
      R3 = P(5)
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
      C3 = P(21)
      RA = P(22)
      CA = P(23)
      R4 = P(24)
      C4 = P(25)
C
      R5 = P(26)
      C5 = P(27)
      RP = P(28)
      CP = P(29)
       L = P(30)
C
C   SET LOGICAL VARIABLES USED FOR ZERO-CHECKING
C
      IF(ICNT.LE.1) THEN  
        RC10 = ((R1.EQ.0.D0).AND.(C1.NE.0.D0))
        RC20 = ((R2.EQ.0.D0).AND.(C2.NE.0.D0))
        RC30 = ((R3.EQ.0.D0).AND.(C3.NE.0.D0))
        RC40 = ((R4.EQ.0.D0).AND.(C4.NE.0.D0))
         R40 = (R4.EQ.0)
         RP0 = ((RP.EQ.0.D0).AND.(CP.NE.0.D0))
         RP1 = (RP.EQ.0)
      ENDIF
C
C   LOOP OVER ALL FREQUENCIES
C
      DO 100 I=1,M
        OMEGA = FREQ(I)
        IOMEGA = DCMPLX(0.D0,OMEGA)
C
C   CALCULATE 6 SUBCIRCUITS IN SERIES
C
        IF (RC10) THEN
          ZT = 1.D0/(C1*IOMEGA)
        ELSE
          ZT = R1/(C1*R1*IOMEGA + 1.D0)
        END IF
        IF (RC20) THEN
          ZT = ZT + 1.D0/(C2*IOMEGA)
        ELSE
          ZT = ZT + R2/(C2*R2*IOMEGA + 1.D0)
        END IF
        IF (RC30) THEN
          ZT = ZT + 1.D0/(C3*IOMEGA)
        ELSE
          ZT = ZT + R3/(C3*R3*IOMEGA + 1.D0)
        END IF
C
        CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE1,TDE1,UDE1,PDE1,NDE1,ZD)
        ZT = ZT + ZD
C
        CALL SDEA(OMEGA,IOMEGA,RA,CA,RDE2,TDE2,UDE2,PDE2,NDE2,ZD)
        ZT = ZT + ZD
C
        CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE3,TDE3,UDE3,PDE3,NDE3,ZD)
        IF(C5.EQ.0) THEN
          Z5 = R5 + ZD
        ELSE
          Z5 = R5 + ZD + 1.D0/(IOMEGA*C5)
        ENDIF
C
        IF (RC40) THEN
          YC = C4*IOMEGA
        ELSE
          IF(R40) THEN
            YC = (0,0)
          ELSE
            YC = (C4*R4*IOMEGA + 1.D0)/R4
          ENDIF
        ENDIF
        IF(Z5.EQ.(0,0)) THEN
          IF(YC.NE.(0,0)) ZT = ZT + 1.D0/YC
        ELSE
          ZT = ZT + Z5/(1.D0 + Z5*YC)
        ENDIF
C
C   CALCULATE EFFECTS OF RP AND CP AND/OR L IF PRESENT
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
        IF(ZT.EQ.(0,0).AND.YC.NE.(0,0)) THEN
          ZT = 1.0/YC
        ELSE
          ZT = (ZT/(1.D0 + ZT*YC))
        ENDIF
        ZT = ZT + L*IOMEGA
C
C   RETURN IMPEDANCE VALUES
C
        F(I) = DREAL(ZT)
        F(I+M) = DIMAG(ZT)
100   CONTINUE
      RETURN
      END
