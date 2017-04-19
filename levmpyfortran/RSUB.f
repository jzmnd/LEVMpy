      SUBROUTINE RSUB(M,FREQ,P,F,NFREE)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'SIZE.INC'
      INTEGER M,NFREE
      DOUBLE PRECISION P(NTOT),F(2*M),FREQ(M)
      COMPLEX*16 ZT,YCC(20),IOMEGA,YC,ZX
      COMMON /CM12/ CELCAP,ATEMP,WF,MAXFEV,ICF,MDE,JCDX
      COMMON /CM47/ ICNT
      COMMON /CM73/ P39
      DIMENSION XTAU(20),NFREE(NTOT)
C
C   **********************     R CIRCUIT:
C
C   CALCULATES IMPEDANCE OF UP TO 19 PARALLEL RC SUBCIRCUITS IN SERIES
C       J. R. MACDONALD 
C
C   IF P40) < 0: COMPLEX DIELECTRIC CONSTANT PARAMETERS (ATEMP < 0) OR
C       COMPLEX CAPACITANCE (ATEMP > 0); OTHERWISE CAP/IMPEDANCE    
C   IF |P(40)|.NE.0 OR 1, THE 2ND. PAR OF DRT INPUT PAIR IS TAU NOT R
C   P(37) CHOICES:  > 0: P(30) IS R IN SERIES  USE P(37) = 10
C           = 0: P(30) IS GP IN PARALLEL
C           < 0 AND >-9 : P(30) IS AN INDUCTANCE L IN SERIES
C           <-9 AND >-19 : P(30) IS A CAPACITANCE C IN SERIES
C   P(39):  IF > 0, IT IS USED IN SETTING TAU MINIMUM VALUE;
C       SEE CHOICE BELOW
C
C         M : number of data points (IN)
C      FREQ : array of frequency values (IN)
C         P : array of model parameters (IN)
C         F : model function values (OUT)
C     NFREE : free parameter array (IN)
C
C   SET PARAMETER VALUES
C
      CP = P(29)
C
      IF(ICNT.LE.1) THEN
        IP40A = INT(DABS(P(40)) + 0.01D0)
        IP38A = INT(DABS(P(38)) + 0.01D0)
         IP38 = NINT(SIGN(1.D0,P(38)))*IP38A
          P39 = P(39)
        IF(ATEMP.LT.0.D0.AND.P(40).LT.0.D0) THEN
            CELC = CELCAP
        ELSE
            CELC = 1.D0
        ENDIF
        IF(IP38.EQ.1) THEN
            JNF = 0
            DO KK = 2,ICF,2
                IF(KK.GT.28) THEN
                    KW = KK + 12
                ELSE 
                    KW = KK
                ENDIF
                JNF = JNF + NFREE(KW)
            ENDDO
            IF(JNF.NE.0) THEN
                WRITE(*,441)
441   FORMAT(1X,'**** WHEN P(38)=1, ALL TAUS MUST BE FIXED ****')
            RETURN
C            STOP
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
C       &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C        CALCULATE ICFH=ICF/2 SUBCIRCUITS IN SERIES
C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
        JJ = 0
765     DO 161 KK = 1,ICF-1
          IF(KK.GT.28) THEN
              KQ = KK + 12
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
            ELSEIF(IP40A.EQ.5) THEN
C       TAU IS A PARAMETER VARIABLE; P(KQ) IS TRANSIENT RESP.
C       OMEGA IS HERE TIME
              XTAU(JJ) = P(JK)
               YCC(JJ) = DCMPLX(P(KQ)*DEXP(-OMEGA/XTAU(JJ)),0.D0) 
            ENDIF
          IF(XTAU(JJ).EQ.0.D0) GOTO 732
          ENDIF
173     FORMAT(2X,I4,1P,(5E15.6))
C
161     CONTINUE
C
C       ***************************************************
C
        CALL DXSPS(XTAU,YCC,ZT,WF,JJ,IP38A)
        ZX = ZT
C
        IF(P(40).LT.0.D0.AND.ZT.NE.0.D0) THEN
C     ** HERE, ZT IS ACTUAL COMPLEX DIELECTRIC CONSTANT OR CAPACITANCE
            ZT = ZX
            YC = IOMEGA*CELC*ZT
            ZT = 1.D0/YC
        ENDIF
C
C    CALCULATE EFFECTS OF GP AND CP IF PRESENT
C
        IF(ATEMP.LT.0.D0.AND.P(37).GT.1.9D1) GOTO 432
C
732     IF(P(37).LE.-2.D1)  ZT = ZT + 1.D0/(IOMEGA*P(30))
C
        IF(P(37).EQ.0.D0) THEN
          YC = P(30) + CP*IOMEGA*CELC         !PARALLEL G
        ELSE
          YC = CP*IOMEGA*CELC
        ENDIF
C
        IF(ZT.EQ.(0,0).AND.YC.NE.(0,0)) THEN
          ZT = 1.0/YC
        ELSE
C       IF(P(37).GT.1.2D1) ZT = ZT + P(30)
          ZT = ZT/(1.D0 + ZT*YC)
        ENDIF
C
        IF(P(37).GT.0.D0) THEN
           ZT = ZT + P(30)            !SERIES R
        ELSEIF(P(37).LT.0.D0.AND.P(37).GT.-9.D0) THEN
           ZT = ZT + IOMEGA*P(30)         !SERIES L
        ELSEIF(P(37).LT.-9.D0.AND.P(37).GT.-19.D0) THEN
           ZT = ZT + 1.D0/(IOMEGA*P(30))      !SERIES C
        ENDIF
C
432   CONTINUE
C
C   RETURN IMPEDANCE VALUES
C
        F(I) = DREAL(ZT)
        F(I+M) = DIMAG(ZT)
100   CONTINUE
      RETURN
      END
