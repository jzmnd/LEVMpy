      SUBROUTINE DSUB(M,FREQ,P,F)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER M
      EXTERNAL DAEFN,GDAEFN
      DOUBLE PRECISION P(*),F(*),FREQ(*),L
      COMPLEX*16 ZA,ZB,ZC,ZD,ZT,ZDAE,IOMEGA,YC
C              COMMON TO PASS DATA TO DAEFN FUNCTION
      COMMON /CM5/ TOMEGA,PHICOM,IV
C              COMMON TO PASS DATA TO GDAEFN FUNCTION
      COMMON /CM9/ TOMEGG,XIGDAE,CHGDAE
      COMMON /CM47/ ICNT    
      LOGICAL RD40,RD50,C40,C50,RC23,DAE0,RA0,RP0,RP1
C
C   *********************   D CIRCUIT:
C
C   THIS SUBROUTINE CALCULATES THE IMPEDANCE OF A SET OF NESTED
C   R-C/DISTRIBUTED ELEMENT CIRCUITS AND INCLUDES A POSSIBLE DAE.
C   FOR MORE DETAILED INFO, CONSULT THE CIRCUIT MODEL DOCUMENTATION. 
C   U1 AND U2 ARE LOG VARIABLES
C
C     M: number of data points (IN)
C     FREQ: array of frequency values (IN)
C     P: array of model parameters (IN)
C     F: model function values (OUT)
C
C   SET PARAMETER VALUES
C
      R1 = P(1)
      RA = P(2)
      CA = P(3)
C
      RDAE = P(4)
      TDAE = P(5)
      U1 = P(6)
      U2 = P(7)
      PHI1 = P(8)
      PHI2 = P(9)   
C
      R2 = P(10)
      C2 = P(11)
      R3 = P(12)
      R4 = P(13)
      C4 = P(14)
      DCH = P(15)
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
      R5 = P(26)
      C5 = P(27)
      RP = P(28)
      CP = P(29)
      L = P(30)
C
C   SET VARIABLES FOR ZERO CHECKING
C
      IF(ICNT.LE.1) THEN   
        RC23 = ((R3.EQ.0.D0).AND.(C2.NE.0.D0))
        RD40 = ((R4.EQ.0.D0).AND.(NDE4.EQ.0))
        RD50 = ((R5.EQ.0.D0).AND.(NDE5.EQ.0))
         C40 = (C4.EQ.0.D0)
         C50 = (C5.EQ.0.D0)
        DAE0 = (RDAE.EQ.0)
         RA0 = ((RA.EQ.0.D0).AND.(CA.NE.0))
         RP0 = ((RP.EQ.0.D0).AND.(CP.NE.0.D0))
         RP1 = (RP.EQ.0)
      ENDIF    
      IF(DAE0) GOTO 62
C
C   EXPON. DISTRIBUTION OF ACTIVATION ENERGIES, FULL EXPRESSION
C
C       EDAE PRELIMINARY CALCULATIONS
      IF(P(15).EQ.0) THEN
C
C         WHEN PHI1>800 (FIXED), EDAE1 CASE if also U2=U1
C
      IF(PHI1.GT.8.D2) PHI1 = PHI2
C
C         WHEN XU=U2 > 800(FIXED), U2=U1
C
      IF(U2.GT.8.D2) U2=U1
C
C         WHEN PHI2>800(FIXED), PHI2=-PHI1
C
      IF(PHI2.GT.8.D2) PHI2=-PHI1
C
      IF(TDAE.LT.0.D0) THEN
        TDAE = -TDAE
        TMLT = DEXP(-U1)
      ELSE
        TMLT = 1.D0
      ENDIF
C
C   CALCULATE EDAE NORMALIZATION COEFFICIENT
C
      IF (PHI1.EQ.0.D0) THEN
        RNORM = U1
      ELSE
        PHU1 = PHI1*U1
        IF(DABS(PHU1).GT.1.D2) PHU1 = DSIGN(1.D2,PHU1)
        RNORM = (1.D0 - DEXP(-PHU1))/PHI1
      ENDIF
      IF(PHI2.EQ.0.D0) THEN
        RNORM = RNORM + U2
      ELSE
        PHU2 = PHI2*U2
        IF(DABS(PHU2).GT.1.D2) PHU2 = DSIGN(1.D2,PHU2)
        RNORM = RNORM + (DEXP(PHU2) - 1.D0)/PHI2
      ENDIF
C
        RNORM = 1.D0/RNORM
      ELSE  
C
C   GAUSSIAN DISTRIBUTION OF ACTIVATION ENERGIES, FULL EXPRESSION
         RGDAE = RDAE
         TGDAE = TDAE
        XLGDAE = U1
        XUGDAE = U2
        XIGDAE = PHI1 
        THGDAE = PHI2 
C
C         CALCULATE LOWER, UPPER GDAE FINAL LIMITS
C         WHEN XL > 800, IT IS FIXED. MAKES UPPER, LOWER LIMITS SAME
          IF(XLGDAE.GT.8.D2) XLGDAE = XUGDAE
C
            SARG = 0.5D0*THGDAE*(XIGDAE**2)
            IF(DABS(SARG).GT.1.D2) SARG = DSIGN(1.D2,SARG)
                UPPER1 = XUGDAE - SARG
                BOTTM = -XLGDAE - SARG
C
C            CALCULATE GDAE NORMALIZATION COEFFICIENT
C
                CHGDAE = 0.D0
                TOMEGG = 0.D0
                    CALL QROMB(GDAEFN,BOTTM,UPPER1,VINT)
                        RNORM = 1.D0/VINT
      ENDIF 
62    CONTINUE
C
C   LOOP OVER ALL FREQUENCIES
C
      DO 100 I=1,M
        OMEGA = FREQ(I)
        IOMEGA = DCMPLX(0.D0,OMEGA)
        IF(DAE0) GOTO 38
C
        IF(P(15).EQ.0) THEN   
C        EXPON. DISTRIBUTION OF ACTIVATION ENERGIES, FULL EXPRESSION
C
C        SET UP VARIABLES TO PASS DATA TO EDAE INTEGRAL FUNCTION
C
          TOMEGA = TDAE*OMEGA*TMLT
C
C        REAL:
          IV = 1
          PHICOM = PHI1
          CALL QROMB(DAEFN,0.D0,U1,Z1)
          IF(U2.EQ.0.D0) GOTO 172
          IV = -1
          PHICOM = -PHI2
          CALL QROMB(DAEFN,0.D0,U2,OA)
          Z1 = Z1 + OA
C
C        IMAGINARY:
172       IV = 1
          PHICOM = PHI1 + 1.D0
          CALL QROMB(DAEFN,0.D0,U1,Z2)
          IF(U2.EQ.0.D0) GOTO 173
          IV = -1
          PHICOM = -PHI2 - 1.D0
          CALL QROMB(DAEFN,0.D0,U2,OB)
          Z2 = Z2 + OB
C
173       ZDAE = RDAE * RNORM * DCMPLX(DBLE(Z1), -TOMEGA*DBLE(Z2))
C
C           END OF EDAE CALCULATION
C
        ELSE
C
C      GAUSSIAN DISTRIBUTION OF ACTIVATION ENERGIES, FULL EXPRESSION
C
C    SET UP NORM. FREQ. VARIABLE TO PASS DATA TO GDAE INTEGRAL FUNCTION
C
          TOMEGG = TGDAE*OMEGA
C
C        IMAGINARY:
          CHGDAE = 1.D0
          CALL QROMB(GDAEFN,BOTTM,UPPER1,VINT)
          Z2 = VINT
C
C        REAL:
          CHGDAE = 0.D0
          CALL QROMB(GDAEFN,BOTTM,UPPER1,VINT)
          Z1 = VINT
C
          ZDAE = RGDAE * RNORM * DCMPLX(DBLE(Z1), -TOMEGG*DBLE(Z2))
C
C           END OF GDAE CALCULATION
C
        ENDIF
38    CONTINUE
C
C   CALCULATE IMPEDANCE OF REST OF CIRCUIT
C
        IF(RD50) THEN
          IF(C50) THEN
            ZC = (0,0) 
          ELSE
            ZC = 1.D0/(C5*IOMEGA)
          ENDIF
        ELSE
          CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE5,TDE5,UDE5,PDE5,NDE5,ZD)
          ZC = R5 + ZD
          ZC = ZC/(1.D0 + IOMEGA*C5*ZC)
        ENDIF
C
        IF(RD40) THEN
          IF(C40) THEN
            ZB = (0,0) 
          ELSE
            ZB = 1.D0/(C4*IOMEGA)
          ENDIF
        ELSE
          CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE4,TDE4,UDE4,PDE4,NDE4,ZD)
          ZB = R4 + ZD
          ZB = ZB/(1.D0 + IOMEGA*C4*ZB)
        ENDIF
C
        ZC = ZC + ZB
        IF(RC23) THEN
          ZA = 1.D0/(C2*IOMEGA)
        ELSE
          ZA = R3/(C2*R3*IOMEGA + 1.D0)
        ENDIF
        ZA = ZA + R2
C
        IF(ZC.EQ.(0,0)) THEN
          ZC = ZA
        ELSEIF(ZA.EQ.(0,0)) THEN
          ZC = ZC
        ELSE
          ZC = ZC*ZA/(ZC + ZA)
        ENDIF
C       
        IF(RA0) THEN
          ZA = 1.D0/(CA*IOMEGA)
        ELSE
          ZA = RA/(CA*RA*IOMEGA + 1.D0)
        ENDIF
        IF(DAE0) THEN
          ZT = ZA
        ELSEIF(ZA.EQ.(0,0)) THEN
          ZT = ZDAE
        ELSE
          ZT = ZDAE*ZA/(ZDAE + ZA)
        ENDIF   
        ZT = ZT + ZC + R1
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
C   RETURN IMPEDANCE VALUES
C
        F(I) = DREAL(ZT)
        F(I+M) = DIMAG(ZT)
  100 CONTINUE
      RETURN
      END
