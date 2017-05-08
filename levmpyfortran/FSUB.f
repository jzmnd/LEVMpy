      SUBROUTINE FSUB(M,FREQ,P,F)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'SIZE.INC'
      INTEGER M
      DOUBLE PRECISION P(NTOT),F(NPT2),FREQ(NPT2),L
      COMPLEX*16 Z(6),Y(6),ZA,ZB,ZD,ZT,YT,ZC1,ZC2,ZC5,IOMEGA,YC
      COMMON /CM47/ ICNT    
      LOGICAL RC60,RC70,RC90,C10,C20,C40,C50,RP0,RP1
C              SAVE LOGICAL CHECKS FOR SUBSEQUENT CALLS
      SAVE RC60,RC70,RC90,C10,C20,C40,C50,RP0,RP1
C
C     ***************       F CIRCUIT:
C
C
C   THIS SUBROUTINE CALCULATES THE IMPEDANCE OF A SET OF
C   PARALLEL SERIES ELEMENTS.  FOR MORE DETAILED INFO,
C   CONSULT THE CIRCUIT MODEL DOCUMENTATION.
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
      CA = P(5)
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
      RA = P(16)
      R3 = P(17)
      C3 = P(18)
      C4 = P(19)
      C5 = P(20)
C
      R6 = P(21)
      C6 = P(22)
      R7 = P(23)
      C7 = P(24)
      R8 = P(25)
C
      R9 = P(26)
      C9 = P(27)
      RP = P(28)
      CP = P(29)
       L = P(30)
C
C   SET VARIABLES FOR ZERO CHECKING
C
      IF(ICNT.LE.1) THEN
          RC60 = ((R6.EQ.0.D0).AND.(C6.NE.0.D0))
          RC70 = ((R2.EQ.0.D0).AND.(C7.NE.0.D0))
          RC90 = ((R3.EQ.0.D0).AND.(C9.NE.0.D0))
           C40 = (C4.EQ.0.D0)
           C20 = (C2.EQ.0.D0)
           C10 = (C1.EQ.0.D0)
           C50 = (C5.EQ.0.D0)
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
C   CALCULATE IMPEDANCE OF CIRCUIT
C
        IF(C10) THEN
          ZC1 = (0,0)
        ELSE
          ZC1 = 1.D0/(IOMEGA*C1)
        ENDIF   
        CALL SDEA(OMEGA,IOMEGA,RA,CA,RDE1,TDE1,UDE1,PDE1,NDE1,ZD)
        Z(3) = ZD + R1 + ZC1
C
        IF(C40) THEN
          ZA = R3
        ELSE 
          ZA = R3 + 1.D0/(IOMEGA*C4)
        ENDIF
        IF(ZA.EQ.(0,0).AND.C3.NE.0) THEN
          ZB = 1.D0/(IOMEGA*C3)
        ELSE 
          ZB = ZA/(1.D0 + IOMEGA*C3*ZA)
        ENDIF   
        IF(C20) THEN
          ZC2= (0,0)
        ELSE
          ZC2 = 1.D0/(IOMEGA*C2)
        ENDIF   
          Z(4) = ZB + R2 + ZC2
C
        IF (RC60) THEN
          Z(5) = 1.D0/(C6*IOMEGA)
        ELSE
          Z(5) = R6/(C6*R6*IOMEGA + 1.D0)
        END IF
        IF (RC70) THEN
          Z(5) = Z(5) + 1.D0/(C7*IOMEGA)
        ELSE
          Z(5) = Z(5) + R7/(C7*R7*IOMEGA + 1.D0)
        END IF
        IF(C50) THEN
          ZC5= (0,0)
        ELSE
          ZC5 = 1.D0/(IOMEGA*C5)
        ENDIF   
        Z(5) = Z(5) + ZC5
C
        IF(RC90) THEN
          Z(6) = 1.D0/(C9*IOMEGA)
        ELSE
          Z(6) = R9/(C9*R9*IOMEGA + 1.D0)
        END IF
        CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE2,TDE2,UDE2,PDE2,NDE2,ZD)
        Z(6) = Z(6) + ZD + R8
C
        YT = (0,0)
        DO 19 IJ = 3,6
          IF(Z(IJ).NE.(0,0)) THEN
              Y(IJ) = 1.D0/Z(IJ)
          ELSE
              Y(IJ) = (0,0)
          ENDIF
          YT = YT + Y(IJ)
19      CONTINUE
        IF(YT.NE.(0,0)) THEN
          ZT = 1.D0/YT
        ELSE
          ZT = (0,0)
        ENDIF
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
