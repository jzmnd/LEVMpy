C                 LEVM SUBROUTINE OSUB:  LV6.FOR
C           J. Ross Macdonald
C   6/15/02; 8/4/02; 3/29/04
C   ADDED PX41 = P(41) IN /CM55/   3/14/97
C   ADDED PX45 = P(45) IN /CM55/   3/19/97
C   corrected c3,r3  4/3/99; corrected loop problems (H. Wintle) 7/01
      SUBROUTINE OSUB(M,FREQ,P,F)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER M
      REAL*8 FREQ,P,F
      DIMENSION P(*),F(*),FREQ(*),QD(8),QC(8)
      EXTERNAL DISTEL
      COMPLEX*16 ZT,YT,YZ,YE,DISTEL,IOMEGA,ZR,ZNL,ZS,ZT16,EEM,ZEL,YEL,
     + YS,YES,ZES
      LOGICAL RD30,C30
      COMMON /CM47/ ICNT
      COMMON /CM12/ CELCAP,ATEMP,WF,MAXFEV,ICF,MDE,JCDX
      COMMON /CM13/ RX,TX,UX,PHIZ,XXM1,XX1,XX2,XX3,RN,AIN,ICAV,
     +  NELEM,NCH
      COMMON /CM55/ PX1,PX41,PX45
      COMMON /CM79/ YTT(50)
C
C   ****************** O CIRCUIT: Runs as an external subroutine
C
C
C   SPECIAL FOR CALLING OTHER PROGRAMS AND DOING SEPARATE OR COMBINED
C   CONDUCTIVE AND DIELECTRIC SYSTEM DISPERSIVE (CSD,DSD) RESPONSE FITS.
C   FOUR CHOICES OF FITTING MODEL ARE AVAILABLE FOR EACH TYPE OF SYSTEM.
C   THE FIRST THREE OF THEM INVOLVE DIRECT CALCULATION USING AN
C   APPROPRIATE DRT.  THUS, THEY ALLOW CUTOFFS TO ENSURE PHYSICAL 
C   REALIZABILITY AND MAKE TRANSIENT RESPONSE FITTING POSSIBLE AS WELL. 
C
C   HERE INPUTS P(1) AND P(2) ARE INFINITE FREQUENCY R AND CG VALUES.      
C   IF ATEMP INPUT < 0, THEN P(2) AND P(15) ARE AT EPSILON LEVEL; 
C   OTHERWISE THEY ARE CAPACITANCES
C   
C
C   THE SUBROUTINE GEDAE IS USED TO CALCULATE GENERALIZED EXPONENTIAL
C   DAE RESPONSE AND HAVRILIAC-NEGAMI AND COLE-DAVIDSON RESPONSE.
C
C   SOME DISTRIBUTIONAL AND LIMITING VALUES ARE CALCULATED AND OUTPUT,
C   AT THE END OF THE ORDINARY OUTPUT FOR THE GEDAE RESPONSE MODEL 
C   AND SOME CHOICES OF DCE MODELS.  ALL DRT CALCULATIONS INVOLVE 
C   CUT-OFFS USING U1 OR U1 AND U2.
C
C   CSD RESPONSE: USE P(3) THRU P(10); DSD RESPONSE: P(13) THRU P(20).
C
C   THE FIVE AVAILABLE CHOICES FOR CSD AND FOR DSD MODELS ARE SELECTED
C   BY THE INPUT VALUES OF P(35) (CSD) AND P(40) (DSD).  THEY ARE:
C    P(35) AND/OR P(40): NCH
C   0   --: NO RESPONSE - SKIP.  ONLY P(1) AND P(2) ARE EFFECTIVE
C   1   GE: GENERALIZED EXPONENTIAL RESPONSE (THREE TYPES) 
C   2   CD: COLE-DAVIDSON RESPONSE. 
C   3   HN: HAVRILIAK-NEGAMI RESPONSE (NOT INCLUDING CD)
C   4   DC: DISTRIBUTED CIRCUIT ELEMENT RESPONSE
C
C   FOR GE, VALUES OF P(10) AND/OR P(20) HAVE NO EFFECT EXCEPT IN CASE
C   c (SEE BELOW). 
C
C   FOR CD AND HN, P(10) AND/OR P(20) MUST BE SET TO 6.  THEY MAY HAVE 
C   ANY VALID VALUES FOR DC MODELS (CHOICE 4 ABOVE).
C
C   THUS, e.g., CHOICES FOR [P(35),P(40)] YIELD FOR CSD AND DSD RESPONSE:
C   [1,1]: GE,GE.  [1,3]: GE,HN.  [2,0]: CD CSD RESPONSE ALONE.  USING
C   THE VALUE 0 ELIMINATES ALL EFFECTS, EVEN WHEN THE RELEVANT
C   PARAMETER VALUES ARE PRESENT.
C
C  **   THERE ARE THREE DIFFERENT GE FIT MODELS AVAILABLE:
C
C  a. U1 < 0, U2 = 0:  USE GAM1 ONLY.  ASYMMETRIC.  IF GAM1=1: EDAE1 CASE
C
C  b. U1 > 0, U2 = 0:  USE GAM1 ONLY.  SYMMETRIC.   IF GAM1=1: EDAE2 CASE
C   AND IF GAM1=2: GAUSSIAN DRT CASE. NOTE: 0 < GAM1 < INFINITY
C
C  c. U1 < 0, U2 > 0:  USE GAM1 AND GAM2. MOST GENERAL: POSSIBILITIY OF
C                DIFFERENT U'S, PHI'S, AND GAMMA'S.
C
C   NOTE THAT PHI1, GAM1, AND U1 AFFECT THE HIGH FREQUENCY REGION OF
C   THE RESPONSE AND PHI2,GAM2,AND U2 THE LOW-FREQUENCY REGION.  THUS,
C   -PHI1 IS THE LOG-LOG SLOPE OF A HF REGION (ESPECIALLY WHEN |PHI1|
C   IS 0.5 OR LESS) OF SAY -IM(Z) VS. OMEGA, AND PHI2 IS THE SLOPE OF
C   A CORRESPONDING LF REGION.
C    
C   FOR CSD RESPONSE: ALWAYS GAM1 = P(3) AND GAM2 = P(4).  
C   FOR DSD RESPONSE: ALWAYS GAM1 = P(13) AND GAM2 = P(14)
C
C   AS USUAL, WHEN MDE (=ORIGINAL INPUT VALUE OF MODE) >= 0, ORDINARY
C   CSD (=CSD0) OR DSD RESPONSE IS CALCULATED.  WHILE WHEN MDE <
C   0, CSD = CSD1 JRM/MOYNIHAN RESPONSE IS CALCULATED.  THE SIGN OF 
C   MDE HAS NO EFFECT ON DSD-TYPE GEDAE RESPONSE.  WHEN |MDE| = 8,
C   TRANSIENT, RATHER THAN FREQUENCY, RESPONSE IS CALCULATED AND IS 
C   AVAILABLE FOR FITTING TRANSIENT DATA.
C
C   WHEN NCH < 4, AND FOR NCH = 4 AND P(10) AND/OR P(20) = 10 OR 
C   32, MOMENTS OF THE DISTRIBUTION USED AND CERTAIN LIMITING
C   EPSILON VALUES OF THE TOTAL FITTING MODEL APPEAR AT THE END OF
C   THE ORDINARY FITTING OUTPUT.
C
C   FOR NCH = 2 OR 3, THE SIGN OF U1 DOES NOT AFFECT THE CALCULATION 
C   RESULT, BUT A NEGATIVE VALUE TURNS ON SCREEN OUTPUT SHOWING THE 
C   PROGRESSION OF THE CALCULATION: SHOWS DATA-POINT COUNT.  THESE
C   CALCULATIONS BECOME VERY SLOW FOR U1 > 15, FOR A LARGE NUMBER
C   OF DATA POINTS, AND FOR IGACC > 3.  LUCKILY, IGACC = 2 OR 3 YIELDS
C   ADEQUATE ACCURACY FOR HN AND CD FITTING.    
C
C   WHEN P(37) > 0: P30 IS A SERIES RESISTANCE, RS; 
C   WHEN P(37) =-4: P30 IS A SERIES INDUCTANCE, L.
C   WHEN P(30) =-16: P30 IS A SERIES CAPACITANCE, CS
C
C   FOR DSD ALONE: IF P(5) NOT 0, NO FINAL EXTRA OUTPUT APPEARS. BUT
C   IF P(5) = 0, FINAL OUTPUT SHOWS MOMENTS OF THE DIELECTRIC
C   DISTRIBUTION AND ALSO INCLUDES LIMITING EPSILON VALUES, ETC
C   ASSOCIATED WITH CSD (CDD0 OR CSD1) RESPONSE WHICH WOULD FOLLOW
C   IF THE DISTRIBUTION WERE USED TO CALCULATE CSD RESPONSE. 
C
C   THE QUANTITIES RN AND AIN ARE NCH < 4 DISTRIBUTION NORMALIZATION-
C   RELATED VALUES.  THEY ARE EQUAL WHEN ONLY CSD OR ONLY DSD RESPONSE
C   IS CALCULATED, BUT WHEN BOTH ARE PRESENT, AIN IS THAT FOR THE CSD 
C   RESPONSE AND RN IS THAT FOR THE DSD PART, WHICH IS CALCULATED
C   LAST OF THE TWO.  ZERO VALUES INDICATE A QUANTITY NOT CALCULATED.  
C
C   RN AND AIN ARE NOT ACTUAL NORMALIZATION VALUES BUT INVOLVE
C   APPROXIMATE OR EXACT CALCULATED NORMALIZATION VALUES.  ACTUAL, 
C   DIRECT NORMALIZATION VALUES ARE CALCULATED BY QUADRATURE FROM THE
C   DISTRIBUTION ITSELF AND THEN THEMSELVES NORMALIZED WITH THE 
C   VALUES CALCULATED DIRECTLY FROM NORMALIZATION FORMULAS.  WHEN
C   THE LATTER ARE EXACT, THE VALUES OF RN AND AIN ARE THUS UNITY.
C
C   FOR NCH = 1 AND GAM1 = 1, THE CALCULATED VALUES TAKE EXACT ACCOUNT
C   OF THE FINITE SIZE OF U1, SO ACTUAL VALUES OF RN AND AIN ARE THEN
C   ESSENTIALLY UNITY.  BUT FOR GAM1 UNEQUAL TO 1, THESE VALUES MAY
C   BE GREATER OR LESS THAN UNITY BUT APPROACH IT AS |U1| INCREASES. 
C   THUS IN THIS CASE, THE DIFFERENCE FROM UNITY IS A MEASURE OF THE
C   EFFECT OF A FINITE CUTOFF VALUE, ALLOWING ONE TO SEE HOW LARGE A
C   VALUE OF |U1| IS NECESSARY TO YIELD NEGLIGIBLE CUTOFF EFFECTS.
C
C   FOR HN/CD CASES, THE CALCULATED NORMALIZATION IS EXACT ONLY IN THE
C   NO CUTOFF LIMIT.  THIS IS OFTEN WELL APPROXIMATED BY A VALUE OF
C   U1 OF 30 OR MORE.  THUS FOR NCH = 1 OR 2, AND U1 LESS THAN
C   INFINITY, RN WILL EXCEED UNITY, AND THE DEGREE TO WHICH NO CUTOFF
C   RESPONSE IS WELL APPROXIMATED IS MEASURED BY THE EXCESS OVER UNITY.
C   A VALUE OF 1.01 THEREFORE INDICATES A PRETTY GOOD APPROXIMATION.
C   GENERALLY FOR HN/CD, WE EXPECT U1 TO FALL IN THE RANGE 5 < U < 25,
C   BUT A LIMITED DATA RANGE OFTEN PRECLUDES THE POSSIBILITY OF FINDING
C   ITS MOST APPROPRIATE VALUE BY TAKING IT FREE TO VARY IN THE FITTING.     
C
C     M: number of data points (IN)
C     FREQ: array of frequency values (IN)
C     P: array of model parameters (IN)
C     F: model function values (OUT)
C
C   SET PARAMETER VALUES
C
      P2 = P(2)
      R2 = P(26)
      C2 = P(27)
      R3 = P(28)
      C3 = P(29)
C
      MDEO = MDE
      IF(MDE.EQ.-16) MDES = MDE
      IF(ICNT.LE.1) THEN
        NCH1 = IDINT(P(35)+0.001D0)
        NCH2 = IDINT(P(40)+0.001D0)
C
C   ADDED BELOW IF ON 4/4/97  CHANGED CM13 NCH1 TO NCH
        IF(NCH1.EQ.0) THEN
            NCH = NCH2
        ELSE
            NCH = NCH1
        ENDIF
        NDE1 = IDINT(DABS(P(10))+0.001D0)  
        NDE1 = NDE1*DSIGN(1.D0,P(10))
        NDE2 = IDINT(DABS(P(20))+0.001D0)
        NDE2 = NDE2*DSIGN(1.D0,P(20))
C
      IF(NCH1.NE.3) THEN
        PX1 = P(1)
        PX41 = P(41)
        PX45 = P(45)
        JCDX = 0
      ELSE
        PX1 = 0.D0
        JCDX = 1
      ENDIF
C
C      THE FOLLOWING DCE APPEARS IN SERIES WITH THE R3 OF THE CIRCUIT
C
      NDE3 = IDINT(P(25)+.001D0)
C
      RD30 = ((NDE3.EQ.0.D0).AND.(R3.EQ.0.D0))
      C30 = C3.EQ.0.D0
      ZNL = DCMPLX(0.D0,0.D0)
C
C           CSD & DSD CASES; SET U2 = 0:
C
      IF(NCH1.GT.1.AND.NDE1.NE.9) THEN
        P(8) = 0
      ELSE
        IF(P(3).EQ.0.D0) P(3) = 1.D0
        IF(P(4).EQ.0.D0) P(4) = 1.D0
        IF(P(5).EQ.0.D0.AND.P(6).EQ.0.D0) NCH1 = 0
      ENDIF
      IF(NCH2.GT.1.AND.NDE2.NE.9) THEN
        P(18) = 0
      ELSE
        IF(P(13).EQ.0.D0) P(13) = 1.D0
        IF(P(14).EQ.0.D0) P(14) = 1.D0
        IF(P(15).EQ.0.D0.AND.P(16).EQ.0.D0) NCH2 = 0
      ENDIF
      INH = 0
      IF(NCH1.GT.0.AND.NCH2.GT.0) INH = 1
C
      IF(ATEMP.LT.0.D0.OR.ATEMP.GT.9.99D2) THEN
        CELC = CELCAP
      ELSE
        CELC = 1.D0
      ENDIF
      ENDIF
C
C   LOOP OVER ALL FREQUENCIES
C
      MDE = MDEO
      DO 100 I=1,M
        IOMEGA = DCMPLX(0.D0,FREQ(I))
        IF(NCH1.EQ.0) THEN
          ZT = ZNL
          GOTO 757        
        ELSEIF(NCH1.LE.3) THEN
          IF(I.EQ.1) THEN
            QC(1) = P(5)
            QC(2) = P(6)
            QC(3) = P(7)
            QC(4) = P(8)
            QC(5) = P(9)
            QC(6) = P(10)
            QC(7) = P(3)
            QC(8) = P(4)
          ENDIF
          CALL GEDAE(I,M,FREQ,QC,F,0,NCH1,INH)
          ZT = DCMPLX(F(I),F(I+M))
        ELSEIF(NCH1.GE.4) THEN
          JCDX = 0
          IF(MDES.EQ.-16) THEN
            MDE = -1
            ZS = DISTEL(P(5),P(6),P(7),P(9),NDE1,FREQ(I))
            ZT16 = ZS
            MDE = 1
            ZR = DISTEL(P(5),P(6),P(7),P(9),NDE1,FREQ(I))/P(5)
            ZT = ZT16/ZR
            MDE = MDES
          ELSE
            MDE = MDEO
            ZT = DISTEL(P(5),P(6),P(7),P(9),NDE1,FREQ(I))
            IF(ATEMP.EQ.-6.4D1) THEN
                YTI = 1.D0/YTT(I+M)
                GB = YTI*DEXP(GAMMLN(YTI))
                ZT = ZT*GB
            ENDIF
          ENDIF
C
C   NGAI COUPLING MODEL:
          IF((NCH1.EQ.5.OR.NCH1.EQ.6).AND.P(41).GT.0.D0) THEN
                  WT =1.D0/P(41)
            IF(FREQ(I).GE.WT) ZT = P(42)/(1.D0 + IOMEGA*P(43))
          ENDIF
C
        ELSE
          WRITE(*,*) 22,NCH1
          RETURN
C        STOP
        ENDIF
C
C   P(1) IS HF LIMIT OF RINF 
757     ZT = ZT + P(1)
        IF(ZT.NE.ZNL) THEN
              YZ = 1.D0/ZT
        ELSE
              YZ = ZNL
        ENDIF
C
C   DO DSD RESPONSE
C   USE COMPLEX EPSILON PARAMETERS HERE (P2,P15).
C
        IF(NCH2.EQ.0) THEN
          YE = ZNL
          GOTO 969
        ELSEIF(NCH2.LE.3) THEN
          IF(I.EQ.1) THEN
            QD(1) = P(15)
            QD(2) = P(16)
            QD(3) = P(17)
            QD(4) = P(18)
            QD(5) = P(19)
            QD(6) = P(20)
            QD(7) = P(13)
            QD(8) = P(14)
          ENDIF
          CALL GEDAE(I,M,FREQ,QD,F,1,NCH2,INH)
          YT = DCMPLX(F(I),F(I+M))
        ELSEIF(NCH2.GE.4) THEN
C
C   ACTUALLY, HERE YT IS A COMPLEX CAPACITANCE (ATEMP >=0) OR
C        EPSILON (ATEMP < 0)
C
          JCDX = 1
          YT = DISTEL(P(15),P(16),P(17),P(19),NDE2,FREQ(I))
C
C   NGAI COUPLING MODEL:
          IF((NCH2.EQ.5.OR.NCH2.EQ.6).AND.P(41).GT.0.D0) THEN
            WT =1.D0/P(41)
            IF(FREQ(I).GE.WT) YT = P(42)/(1.D0 + IOMEGA*P(43))
          ENDIF
        ELSE
          RETURN
C        STOP
        ENDIF
C
        IF(P(33).GT.0.D0) THEN
          YZ = YZ/(CELCAP*IOMEGA)
          EEF = 3.D0*P(33)*(YT - YZ)/(YT + 2.D0*YZ - P(33)*(YT - YZ))
          EEM = P2 + YZ*(1.D0 + EEF)
          ZT = 1.D0/(IOMEGA*CELCAP*EEM)   !CONVERT FROM E TO Z
          GOTO 197
        ENDIF
C
        F(I) = DREAL(YT)
        F(I+M) = DIMAG(YT)
        IF(F(I).EQ.0.D0) THEN
C               ! CONVERT EPS TO Y
            YE = ZNL
        ELSE
            YE = CELC*FREQ(I)*DCMPLX(-F(I+M),F(I))
        ENDIF
C
C   COMBINE CSD AND DSD RESPONSES,CINF, and parallel conductance P(12)
C
C   P2 IS HIGH-FREQUENCY VALUE OF CG=CINF OR EPSILONINF
C
969     YT = YZ + YE + P(12) + CELC*IOMEGA*P2
        ABSYT = REAL(YT)*REAL(YT) + IMAG(YT)*IMAG(YT)
        IF(ABSYT.LT.1D-150) ABSYT = 0.D0
561     IF(ABSYT.EQ.0.D0) THEN   ! FOR PS1
            ZT = ZNL
        ELSE
            ZT = 1.D0/YT
        ENDIF
C
C   CALCULATE CONTRIBUTIONS OF DCE3, L, AND R2,C2,R3,C3, AND,
C     IF P(38.NE.0, THEN ADDITIONAL DE = DE4 IN SERIES WITH 
C     P(11), NOW AS CAPACITANCE
C
        IF(P(38).EQ.6.4D1.AND.P(35).NE.0.D0) THEN
          NDE4 = IDINT(DABS(P(45))+0.001D0)
C
          ZEL = DISTEL(P(41),P(42),P(43),P(44),NDE4,FREQ(I))
          YE = ZNL ! HERE, P(11) IS A CAPACITANCE
          IF(P(11).NE.0.D0) ZEL = ZEL + 1.D0/(IOMEGA*P(11))
C
          IF(ZEL.NE.ZNL) THEN
            YEL = 1.D0/ZEL
          ELSE
            YEL = ZNL
          ENDIF
            P11 = 0.D0
C
        ELSE          ! ORIGINAL SITUATION
           YEL = 0.D0
           P11 = P(11)  ! HERE, P(11) IS A PARALLEL RESISTANCE
        ENDIF
C
        IF(RD30) THEN
          IF(C30) THEN
              ZS = ZNL
          ELSE
              YS = IOMEGA*C3
              ZS = 1.D0/YS
          ENDIF
        ELSE
          IF(P(23).LT.-9.0D2.OR.ATEMP.LT.-9.0D2) THEN
              MDE = 1
          ELSE
              MDE = -1
          ENDIF
          ZS = DISTEL(P(21),P(22),P(23),P(24),NDE3,FREQ(I))
          ZS = R3 + ZS
          IF(ZS.NE.ZNL) THEN
              YS =1.D0/ZS
          ELSE
              YS = ZNL
          ENDIF
C
          MDE = MDEO
        ENDIF
        YES = YS + YEL + IOMEGA*C3 
        IF(YES.NE.ZNL) THEN
              ZES = 1.D0/YES
        ELSE
              ZES = ZNL
        ENDIF
C
        ZR = R2 + ZES
C
C   ADD EFFECTS OF C2 AND P(11) = GE (G IN ELECTRODE PART)
        IF(ZR.EQ.ZNL) THEN
          YE = IOMEGA*C2 + P11
          IF(YE.EQ.0.D0) THEN
              ZR = ZNL
          ELSE
              ZR = 1.D0/YE
          ENDIF
        ELSE
          YE = (1.D0/ZR) + IOMEGA*C2 + P11
          ZR = 1.D0/YE
        ENDIF
C
C   IF P30<0, IT IS L; OTHERWISE IT IS R SUB S, BUT SEE BELOW:
        ZT = ZT + ZR                !NO EXTRA SERIES
        IF(P(37).GT.0.D0) THEN
              ZT = ZT + P(30)         !SERIES R
        ELSEIF(P(37).LT.0.D0.AND.P(37).GT.-9.D0) THEN
              ZT = ZT + IOMEGA*P(30)      !SERIES L
        ELSEIF(P(37).LT.-9.D0.AND.P(37).GT.-1.9D1) THEN
              ZT = ZT + 1.D0/(IOMEGA*P(30))   !SERIES C
        ENDIF
C
      CONTINUE
C
197     F(I) = DREAL(ZT)
        IF(ATEMP.EQ.-64D0) THEN
              F(I+M) = YTT(I+M)
        ELSE
              F(I+M) = DIMAG(ZT)
        ENDIF
100   CONTINUE
      RETURN
      END
