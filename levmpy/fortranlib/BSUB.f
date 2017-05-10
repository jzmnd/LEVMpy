      SUBROUTINE BSUB(M,FREQ,P,F)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'SIZE.INC'
      INTEGER M
      DOUBLE PRECISION P(NTOT),F(NPT2),FREQ(NPT2),L
      COMPLEX*16 ZD,ZB,ZC,YC,ZT,IOMEGA
      COMMON /CM47/ ICNT    
      LOGICAL RD30,RD20,RD40,C30,C20,R230,RCT,RP0,RP1,RC3N
C              SAVE LOGICAL CHECKS FOR SUBSEQUENT CALLS
      SAVE RD30,RD20,RD40,C30,C20,R230,RCT,RP0,RP1,RC3N
C
C   ********************     B CIRCUIT:
C
C   THIS SUBROUTINE CALCULATES THE IMPEDANCE OF A SET OF NESTED
C     RC DISTRIBUTED ELEMENT CIRCUITS.  Includes a special
C     correlated calculation for some of the parameters (some redefined)
C     when NDE3=13 or 15 and UDE3 < 0.  FOR MORE DETAILED INFO,
C     CONSULT THE CIRCUIT MODEL DOCUMENTATION.
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
      R2 = P(4)
      C2 = P(5)
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
      R3 = P(26)
      C3 = P(27)
      RP = P(28)
      CP = P(29)
       L = P(30)
C
C   WHEN U < 0, CORRELATION BETWEEN SOME PARAMETERS
      IF(UDE3.LT.0.D0) THEN
        IF(UDE3.LT.-9) THEN
          R2 = P(1)*P(4)
          C2 = P(29)*P(5)
C         Here, P4 is normalized R2 and P5 is normalized C2
        ENDIF
C
        IF(NDE3.EQ.13.OR.NDE3.EQ.15) THEN
C         USES GENERAL HOMOGENEOUS DIFFUSION DCE
          IF(R2.EQ.0.D0.OR.P(27).EQ.0.D0) THEN
            RETURN
C           STOP
          ENDIF
          C3 = P(27)/R2 
          R3 = P(26)/C3
          IF(UDE3.NE.2) RDE3 = P(16)*R2
        ENDIF
C   END OF U < 0 CHOICES
      ENDIF
C
C   SET VARIABLES FOR ZERO CHECKING
C
      IF(ICNT.LE.1) THEN  
        RD30 = ((R3.EQ.0.D0).AND.(NDE3.EQ.0))
        RD20 = ((R2.EQ.0.D0).AND.(NDE2.EQ.0))
        RD40 = ((C2.EQ.0.D0).AND.(NDE4.EQ.0))  
         C30 = (C3.EQ.0.D0)
         C20 = (C2.EQ.0.D0)
        R230 = (RD20.AND.RD30)
         RCT = (R230.AND.RD40.AND.C30)
         RP0 = ((RP.EQ.0.D0).AND.(CP.NE.0.D0))
         RP1 = (RP.EQ.0)
        RC3N = ((R3.NE.0).AND.(C3.NE.0))
      ENDIF
C
C   LOOP OVER ALL FREQUENCIES
C
      DO 100 I=1,M
        OMEGA = FREQ(I)
        IOMEGA = DCMPLX(0.D0,OMEGA)
C
C   CALCULATE FIRST NESTED CIRCUIT
C
        IF(RD30) THEN
          IF(C30) THEN
              ZC = (0,0) 
          ELSE
              ZC = 1.D0/(C3*IOMEGA)
          ENDIF
        ELSE
           CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE3,TDE3,UDE3,PDE3,NDE3,ZD)
              ZC = R3 + ZD
                  ZC = ZC/(1.D0 + IOMEGA*C3*ZC)
        ENDIF
           CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE2,TDE2,UDE2,PDE2,NDE2,ZD)
          ZC = ZC + R2 + ZD
        IF(.NOT.RD40) THEN
           CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE4,TDE4,UDE4,PDE4,NDE4,ZD)
          IF(C20) THEN
              ZB = ZD
          ELSE
              ZB = (1.D0 + IOMEGA*C2*ZD)/(IOMEGA*C2)
          ENDIF
        ENDIF
        IF(RCT) THEN
          ZT = (0,0)
        ELSEIF(RD40) THEN
          ZT = ZC
        ELSEIF(R230) THEN
          ZT = ZB
        ELSE
          ZT = ZB*ZC/(ZB + ZC)
        ENDIF
C
C     CALCULATE SECOND NESTED CIRCUIT
C
        CALL SDEA(OMEGA,IOMEGA,RA,CA,RDE1,TDE1,UDE1,PDE1,NDE1,ZD)
        ZT = ZT + R1 + ZD
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
