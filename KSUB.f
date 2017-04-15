      SUBROUTINE KSUB(M,FREQ,P,F,NFREE)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER M,NFREE
      REAL*8 FREQ,P,F
      DOUBLE PRECISION L
      COMPLEX*16 ZT,ZD,YC,IOMEGA,YCC(20),ZX
      COMMON /CM12/ CELCAP,ATEMP,WF,MAXFEV,ICF,MDE,JCDX
      COMMON /CM47/ ICNT
      COMMON /CM73/ P39
      DIMENSION P(*),F(*),FREQ(*),NFREE(*),XTAU(20) 
C
C               **********************     K CIRCUIT:
C
C   CALCULATES IMPEDANCE OF UP TO 11 RC SUBCIRCUITS IN SERIES
C   WITH A DCE AND ALL IN PARALLEL WITH RP AND CP AND IN SERIES WITH L.
C   IF P(40) >= O, CONDUCTIVE DRT; DIELECTRIC DRT IF P(40) < 0
C   TO DO MORE THAN 11 PAIRS, MUST USE RSUB
C
C   P(30) < 0:  P(30) IS A POSITIVE SERIES RESISTANCE, NOT AN INDUCTANCE
C   ATEMP < 0 AND P(40) < 0, DRT PARAMETERS ARE DIEL. CONSTS, NOT CAPS
C
C     M: number of data points (IN)
C     FREQ: array of frequency values (IN)
C     P: array of model parameters (IN)
C     F: model function values (OUT)
C     NFREE: free parameter array (IN)
C
C   SET PARAMETER VALUES
C
      RDE1 = P(23)  
      TDE1 = P(24)
      UDE1 = P(25)
      PDE1 = P(26)
      NDE1 = IDINT(P(27)+0.01D0)
      GP = P(28)
      CP = P(29)
      L =  DABS(P(30))
C
      IF(ICNT.LE.1) THEN
        IP40A = INT(DABS(P(40)) + 0.01D0)
        IP38A = INT(DABS(P(38)) + 0.01D0)
         IP38 = NINT(SIGN(1.D0,P(38)))*IP38A
          P39 = P(39)
C
        IF(ATEMP.LT.0.D0.AND.P(40).LT.0.D0) THEN
            CELC = CELCAP
        ELSE
            CELC = 1.D0
        ENDIF
        IF(IP38.EQ.1) THEN
          JNF = 0
          DO KK = 2,ICF,2
              IF(KK.GT.22) THEN
                  KW = KK + 18
              ELSE 
                  KW = KK
              ENDIF
              JNF = JNF + NFREE(KW)
          ENDDO
          IF(JNF.NE.0) THEN
            WRITE(*,441)
441       FORMAT(1X,'**** WHEN P(38)=1, ALL TAUS MUST BE FIXED ****')
            RETURN
C          STOP
          ENDIF
        ENDIF
      ENDIF
C
C   LOOP OVER ALL FREQUENCIES
C
      DO 100 I=1,M
         OMEGA = FREQ(I)
        IOMEGA = DCMPLX(0.D0,OMEGA)
            ZT = (0.D0,0.D0)
C
C   CALCULATE ICF/2  SUBCIRCUITS IN SERIES OR PARALLEL; MAX 11
C   FOR MORE, USE R CIRCUIT
C
        JJ = 0
        DO 161 KK = 1,ICF-1
          IF(KK.GT.22) THEN
              KQ = KK + 18
          ELSE
              KQ = KK
          ENDIF
          JK = KQ + 1
          IF(MOD(JK,2).EQ.0) THEN
            JJ = JJ + 1
C
            IF(IP40A.EQ.2.OR.IP40A.EQ.3) THEN
C       HERE P(JK) IS TAU, A PARAMETER VARIABLE (OR CONSTANT)
               YCC(JJ) = P(KQ)/(1.D0 + IOMEGA*P(JK)) 
              XTAU(JJ) = P(JK)
            ELSEIF(IP40A.EQ.0.OR.IP40A.EQ.1) THEN
C       TAU IS THE RC PRODUCT HERE
              XTAU(JJ) = P(KQ)*P(JK)
               YCC(JJ) = P(KQ)/(1.D0 + IOMEGA*XTAU(JJ)) 
            ELSEIF(IP40A.EQ.4) THEN
C       TAU IS A PARAMETER VARIABLE
              XTAU(JJ) = P(JK)
               YCC(JJ) = XTAU(JJ)*P(KQ)/(1.D0 + IOMEGA*XTAU(JJ)) 
            ENDIF
          ENDIF
173     FORMAT(2X,I4,1P,(5E15.6))
C
161     CONTINUE
C
        CALL DXSPS(XTAU,YCC,ZT,WF,JJ,IP38A)
C
        ZX = ZT
C
        IF(P(40).LT.0.D0.AND.DREAL(ZT).NE.0.D0) THEN
C     **    HERE, THE DRT ZT IS THE ACTUAL COMPLEX DIELECTRIC CONSTANT
C       (WHEN ATEMP < 0), OR COMPLEX CAPACITANCE OTHERWISE
C
          ZT = ZX
          YC = IOMEGA*CELC*ZT
          ZT = 1.D0/YC
        ENDIF
C
        CALL SDEA(OMEGA,IOMEGA,0.D0,0.D0,RDE1,TDE1,UDE1,PDE1,NDE1,ZD)
C
C   IF P(36).NE.0, THEN CP IS IN PARALLEL WITH DE (ZD HERE)
        IF(P(36).GT.0) THEN
          ZD = 1.D0/(CP*IOMEGA*CELC + (1.D0/ZD))  
        ELSEIF(P(36).LT.0) THEN
          ZT = 1.D0/(CP*IOMEGA*CELC + (1.D0/ZT))  
        ENDIF
C
        ZT = ZT + ZD
C
C    CALCULATE EFFECTS OF GP AND CP AND/OR L IF PRESENT
C
        IF(P(36).EQ.0) THEN
          YC = GP + CP*IOMEGA*CELC
          IF(P(21).LT.0.D0) THEN
              P21 = DABS(P(21))
              YC = YC + IOMEGA*P21/(1.D0 + IOMEGA*P21*P(22))
          ENDIF
        ENDIF
C
        IF(ZT.EQ.(0,0).AND.YC.NE.(0,0)) THEN
          ZT = 1.0/YC
        ELSE
          ZT = (ZT/(1.D0 + ZT*YC))
        ENDIF
        IF(P(30).LT.0.D0) THEN
          ZT = ZT + L
        ELSE
          ZT = ZT + L*IOMEGA
        ENDIF
C
C   RETURN IMPEDANCE VALUES
C
        F(I) = DREAL(ZT)
        F(I+M) = DIMAG(ZT)
100   CONTINUE
      RETURN
      END
