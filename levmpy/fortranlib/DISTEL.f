C     DISTRIBUTED ELEMENTS MULTIFUNCTION
C
C   This function returns impedance values from any of a variety of
C   distributed elements, and may be used in circuit models as if
C   it were any of the elements. The particular element is selected
C   by the 'NELEM' parameter. Current assignments are as follows:
C
C       1 : R-C PARALLEL CIRCUIT
C           Z = R / (1 + i*OMEGA*R*T)
C
C       2 : CPE and/or RC, MODEL 1
C           Z = 1/T*(i*OMEGA)**PHI  (CPE part ONLY)
C
C       3 : CPE and/or RC, MODEL 2
C           Z = 1/(i*T*OMEGA)**PHI (CPE part ONLY)
C
C       4 : Z-C, MODEL 1
C           Z = R / (1 + (R**U)*T*(i*OMEGA)**PHI)
C
C       5 : Z-C, MODEL 2
C           Z = R / (1 + (i*TT*OMEGA)**PHI)
C
C       6 : H-N, MODEL 1
C           Z = R / (1 + (i*TT*OMEGA)**U )**PHI
C
C       7 : H-N, MODEL 2  T=C
C           Z = R / (1 + (i*R*TT*OMEGA)**U )**PHI
C
C       8 : H-N, MODEL 3  RINPUT=R=R0; T=C0
C           Z = R*SIN(CHI) / (1 + (i*R*TT*OMEGA)**U )**PHI
C
C       9 : GENERALIZED FINITE WARBURG
C           Z = R * TANH((i*TT*OMEGA)**PHI) / (i*TT*OMEGA)**PHI
C           and other choices using P8=P(8). See Manual
C
C      10 : WILLIAMS - WATTS FRACTIONAL EXPONENT
C
C      11 : JONSCHER (GENERALIZED TO INCLUDE REAL PART)
C
C      12 : EXPONENTIAL DRT: EDAE1 DAES, ASYMMETRIC. INTEGRAL 0 TO U.
C
C      13 : EXPONENTIAL DRT: EDAE2, SYMMETRIC FORM
C
C      14 : GAUSSIAN DRT: GDAES SYMMETRIC FORM OF GAUSSIAN DAE
C
C      15: GENERAL DIFFUSION DCE: MACROSCOPIC PARAMETERS
C 
C      16: GENERAL DIFFUSION DCE: MICROSCOPIC PARAMETERS
C
C      17: PARALLEL COMBINATION OF C (TDE), R (RDE), AND A SERIES
C          COMBINATION OF ANOTHER R (UDE) AND AN L (PDE). 
C
C      18: DISSADO-HILL FITTING FUNCTION:  Z LEVEL
C
C      19: DISSADO-HILL FITTING FUNCTION:  EPSILON LEVEL
C
C      20: PLE1; Z=R*(OMEGA)**T + i*U*(OMEGA)**PHI
C
C      21: PLE2; E LEVEL. SEE LV5.FOR, Section 2100
C
C      22: PLE3; Y=R*(i*OMEGA)**T + i*U*(OMEGA)**PHI
C
C      23: LADDER; Y=(U/2)*(i*T*OMEGA)*[1+{1+4/(i*R*T*OMEGA}**PHI]
C
C      24: LADDER; Y=(U/2R)*[1 +{1+i*4*R*T*OMEGA}**PHI]
C
C      25: EFFECTIVE MEDIUM; YP=T*(i*OMEGA)**(-PHI)
C          YC=3*U*[YP - R]/[YP + 2R - U*(YP - R)]
C          E=R*(1+YC)   DIELECTRIC lEVEL
C 
C      26: SCPE; SIGMA=U*CELCAP*(i*OMEGA)**PHI=1/RESISTIVITY
C          FOR SPECIFIC DATA, USE CELCAP=EV=PERMITTIVITY OF VACUUM 
C
C      27: PCPE; EPSILON LEVEL E=U*(i*OMEGA)**(-PHI)  
C
C      28: EPSILON LEVEL IF U=PHI=0, THEN Y=R*CELCAP*(OMEGA**T), ELSE
C          E=R*(OMEGA**T)+U*(i*OMEGA)**(-PHI). USE O CIRCUIT, P(40)=4
C
C      29: MODIFIED DAVIDSON-COLE RESPONSE 
C          Z=R*(1+U)/[1+U*(1 + i*TT*OMEGA)**PHI]
C  
C      30: EXACT EDAE1 FOR PHI=0 OR 1 
C
C      31: GBEM EMA FITTING FUNCTION: EFFECTIVE MEDIUM
C
C      32: ACCURATE KWW FOR BETA=0.5; USES DRT
C
C      33: ARRHENIUS AND OTHER THERMAL FITS
C
C      34: SIMPLIFIED GENERALIZED EXPONENTIAL DAE
C
C      35: ACCURATE KWW FOR BETA=1/3; USES DRT
C
C      36: KWW FOR ARBITRARY BETA; USES DRT, SERIES CALCULATIONS, ETC.
C
C      37: CALCULATE KWW DRT FOR ARBITRARY BETA. USE IRCH=0, OR >0 FOR
C          FIT
C
C      38: NOT USED
C
C
C           R : Free parameter, usually a resistance (IN)
C           T : Free parameter, usually a time constant/relaxation time
C               Sometimes a fit model may use TT not = T, when MODE=2 (IN)
C           U : Free parameter, NOT the same as U=P(31) (IN)
C         PHI : Free parameter, typically an exponent (IN)
C       NELEM : Fixed parameter, selects DE to be used (see above) (IN)
C       OMEGA : Fixed parameter, represents frequency of AC signal
C               Sometimes time or temperature (IN)
C
C      DISTEL : Impedance of DE, COMPLEX*16 (OUT)
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      FUNCTION DISTEL(R,T,U,PHI,NELEM,OMEGA)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL DAEFN,DAEFN2,GDAEFN,GDAEFN2,WW5,GAMMLN,AIFU,GDRT
      REAL*8 R,T,U,PHI,OMEGA
      INTEGER NELEM
      LOGICAL TC,RCC,CC,RRC
      COMPLEX*16 CTEMP,EXPSAV,DISTEL,XX,XTX,CMTANH,YC,YP,YCP,
     + COMEG,YP0,YP1,CDCOTH,XP,TS
      DIMENSION GSP(43),XXJ(43),VV(5)
      DIMENSION XLL(43),XLH(43)
C     COMMON TO PASS ADDITIONAL DATA TO DAE INTEGRATION FUNCTIONS
      COMMON /CM5/ TOMEGA,PHICOM,IV
      COMMON /CM9/ TCMEGA,XIGDAE,ICHG,IWT
      COMMON /CM12/ CELCAP,ATEMP,WF,MAXFEV,ICF,MDE,JCDX
      COMMON /CM13/ RZ,TX,UU,PHIX,XXM1,XX1,XX2,XX3,RN,AIN,ICAV,
     + NELEMN,NCH
      COMMON /CM39/ P8,P18
      COMMON /CM55/ PX1,PX41,PX45
C
      DATA PI2/1.5707963267948965D0/,P0/1.D0/,SQPII/0.564189583547756D0
     + /,PE/1.16044485D4/,PI3/0.333333333333333D0/,P02/2.D0/,AKINV/
     + 1.16044485D4/,EVI/1.12940907D13/,TWPI/6.283185307179586D0/,
     + SML/1.D-10/,PII/0.318309886D0/
C
C
      CTEMP = DCMPLX(0.D0,0.D0)
      MDEX = MDE
      NELEMN = NELEM
      PHIX = PHI
C     DEFINE LOGICAL VARIABLES
      TC = (T.EQ.0.D0)
      RCC = ((R.EQ.0.D0).AND.(U.EQ.0.D0))
      CC = (U.EQ.0.D0)
      RRC = (R.EQ.0)
C   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      UU = U
      RZ = R
      TX = T  
C     CHOICE FOR EPSSUBTAU INSTEAD OF TAUSUB0 FOR INPUT T VARIABLE
C
      IF(ABS(MDE).EQ.2.AND.JCDX.EQ.0.AND.R.NE.0.D0) THEN
        RTOT = PX1 + R
          TT = T*CELCAP*RTOT*RTOT/R
      ELSE
          TT = T
      ENDIF
C
C     EXECUTE APPROPRIATE ROUTINE FOR DISTRIBUTED ELEMENT SELECTED
C
      IF (NELEM.EQ.0) THEN
        DISTEL = (0.D0,0.D0)
        MDEX = MDE
        RETURN
      ENDIF
      MDE = MDEX
      GOTO(100,200,300,400,500,600,700,800,900,1000,1100,1200,
     + 1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,
     + 2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,
     + 3700),IABS(NELEM)
C
C   INVALID FUNCTION
C
69    WRITE(*,50) NELEM,R,T,PHI
50    FORMAT(' ',I3,' IS NOT A VALID FUNCTION. PROGRAM TERMINATING',
     + '  R, T, PHI PARAMETERS: ',3D15.8)
      DISTEL = (0.D0,0.D0)
      RETURN
C      STOP
C
C--------------------------------------------------------------------
C
C   R - C PARALLEL CIRCUIT
C
100   IF ((R.EQ.0).AND.(T.GT.0)) THEN
         DISTEL = 1.D0/DCMPLX(0.D0, T*OMEGA)
      ELSE
         DISTEL = R / DCMPLX(1.D0, R*TT*OMEGA)
      ENDIF
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   CONSTANT PHASE ELEMENT, MODEL 1
C
200   IF(RCC.AND.TC) GOTO 69
      IF(RCC) THEN
        YP = DCMPLX(0.D0,0.D0)
        GOTO 212
      ELSEIF(CC) THEN
        YP = DCMPLX(1.D0/R,0.D0)
      ELSE
        YC = DCMPLX(0.D0,OMEGA*U)
        YP = YC/(1.D0 + R*YC)
      ENDIF
212   IF(TC) THEN
        YCP = DCMPLX(0.D0,0.D0)
      ELSE    
            ARG = PI2*PHI
            YCP = T*(OMEGA**PHI)*DCMPLX(DCOS(ARG),DSIN(ARG))
      ENDIF
      DISTEL = 1.D0/(YP + YCP)
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   CONSTANT PHASE ELEMENT, MODEL 2
C
300   COMEG = DCMPLX(0.D0, OMEGA)
      YCP = (R*T*COMEG)**PHI
C
      IF(RCC.OR.TC) THEN
        DISTEL = DCMPLX(0.D0,0.D0)
      ELSE
        DISTEL = R/YCP
      ENDIF
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   ZARC - COLE ELEMENT, MODEL 1
C
400   ARG = PI2*PHI
      CTEMP = T * (OMEGA**PHI) * DCMPLX(DCOS(ARG), DSIN(ARG))
      DISTEL = R / (1.D0 + (R**U)*CTEMP)
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   ZARC - COLE ELEMENT, MODEL 2
C
500   ARG = PI2*PHI
      CTEMP = (TT*OMEGA)**PHI * DCMPLX(DCOS(ARG), DSIN(ARG))
      DISTEL = R / (1.D0 + CTEMP)
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   HAVRILIAK-NEGAMI ELEMENT, MODEL 1.  USES R, U, P8, PHI, AND T
C   TT IS TAU OR EPSX, ETC. VALUE OF U: 0 OR +-1 SELECTS
C
600   XTX = DCMPLX(0.D0,OMEGA)
      XX = TT*XTX
      TOMEG = TT*OMEGA    
      TOMEGU = TOMEG**U   
C
      IF(TOMEGU.LT.1.D14) THEN    
        YC = (1.D0 + XX**U)**(-PHI)      ! IO(OMEGA)
      ELSE
        PARG = PI2*U        
        YP = DCMPLX(DCOS(PARG),DSIN(PARG))
        YCP = TOMEGU*YP
        YC = (1.D0 - PHI/YCP)*YCP**(-PHI)
      ENDIF
C
      IF(MDE.GE.0.OR.DABS(ATEMP).GT.9.D2) THEN
            DISTEL = R*YC
      ELSEIF(DABS(U).EQ.1.D0) THEN            ! CD CASE ONLY
        DISTEL = R*(1.D0 - YC)/(XX*PHI)
      ELSEIF(MDE.NE.-2) THEN              !  HN (NOT CD)
C       IF MODE < 0, FOR HN R (INPUT) BECOMES EPSSUBCINF
        DISTEL = (1.D0 - YC)/(XTX*CELCAP*R) !  R IS EPSCINF HERE
      ELSEIF(-MDE.NE.1.OR.-MDE.NE.2) THEN
        DISTEL = (0.D0,0.D0)
        RETURN
C        STOP
      ENDIF
C
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   HAVRILIAK-NEGAMI ELEMENT, MODEL 2.  USES R, U AND C: SEE NELEM=600
C
700   GOTO 600
      DISTEL = (0.D0,0.D0)
      RETURN
C      STOP
C--------------------------------------------------------------------
C
C   HAVRILIAK-NEGAMI ELEMENT, MODEL 3.
C   H-N FUNCTION USING R0*SIN; RINPUT = R0.  
C               
 800  ARG = PI2 * U
      XX = TT*DCMPLX(0.D0,OMEGA)
      DISTEL = R*DSIN(ARG)*(1.D0 + XX**U)**(-PHI)
      GO TO 5000
C
C-------------------------------------------------------------------
C
C   GENERALIZED FINITE WARBURGS (PARAMETERS R,TAU,U,P8,P18,PHI,OMEGA)
C 
C   THIS SUBROUTINE INSTANTIATES AND ACTIVATES EQS.A.10, A.11,AND A.12 
C   IN THE JRM PUBLISHED PAPER #252, IN ADDITION, IT INCORPORATES TWO
C   EXTENSIONS OF THESE EQUATIONS. FOR THE CHOICE P18=4, ITS TOTAL
C   IMPEDANCE INCLUDES TWO WEIGHTED INTERFACE ADMITTANCES IN PARALLEL, 
C     ONE INVOLVING A U**2 TERM (I*OMEGA*TAU)**PHI AND THE OTHER THE SAME 
C   BUT WITH PHI FIXED AT UNITY. FOR THE P19=8 CHOICE, A SUM OF THESE
C   TWO TERMS, EACH WEIGHTED AS BELOW, IS USED IN CALCULATING THE 
C   INTERFACE IMPEDANCE, DEFINED AS CTEMP BELOW, WEIGHTING INVOLVES THE 
C   P8=AA PARAMETER IN THE FORM AM = (1-AA) AND AA QUANTITIES. SEE PP. 
C     4-9 AND 4-10 IN THE LEVM MANUAL FOR A DISCUSSION OF THE OTHER 
C   P8-SELECTED CHOICES.
C  
900   XTX = DCMPLX(0.D0,OMEGA)  ! XTX =I*OMEGA
      A18 = DABS(P18)
      IF(P18.GT.0.D0) THEN
        JM = 0
        IF(P18.EQ.4.D0.OR.P18.EQ.8.D0) PX = 4
        IF(P18.EQ.2.D0) PX = 8
      ELSEIF(P18.LT.0.D0) THEN
        JM = 1
        IF(P18.EQ.-4.D0.OR.P18.EQ.-8.D0) PX = 3
        IF(P18.EQ.-2.D0) PX = 8
      ELSE
        JM = 0
        PX = P8                ! P18=0
      ENDIF
C
      IF(JM.EQ.0) THEN
        IF(P18.EQ.4.D0.OR.P18.EQ.8.D0) PX = 3
        AA = P8
        AM = 1.D0 - AA
      ENDIF
C      
      IF(PX.LT.6.D0) THEN 
         XX = TT*XTX       ! I*OMEGA*TAU = S**2
      ELSE
         CINF= TT
         TAU = CINF*CELCAP*R       !USES CINF FREE PAR
         XX = TAU*XTX              ! I*OMEGA*TAU 
C 
      ENDIF
      IF(P18.EQ.0.D0) THEN
        YC =  XX**PHI                ! U**2 = (I*OMEGA*TAU)**PHI
        TS = CDSQRT(YC)              
        YP = 1.D0 + YC               !P**2    
        XP = CDSQRT(YP)              !P
        EXPSAV = CDCOTH(U*XP)        !M*XP*COTH(M*XP) SPECIAL DEFINITION
        CTEMP = YP/(YC*(EXPSAV - 1.D0))  !ZSUBIN/R  EQUATION A.11
      ENDIF
      IF(A18.EQ.4.D0) THEN               
        DO KN = 1,2
          IF(KN.EQ.1) THEN
            PHIX = 1.D0
          ELSEIF(KN.EQ.2) THEN 
            PHIX = PHI
          ENDIF
        YC = XX**PHIX          !  (I*OMEGA*TAU)**PHIX             
        YP = 1.D0 + YC               !P**2  
        XP = CDSQRT(YP)              !P
C
C     CALCULATE COMPLEX HYPERBOLIC TANGENTS AND FINAL FUNCTIONS.NOTE
C     THAT THE FUNCTION CDCOTH(X) IS DEFINED AS X/TANH(X)  ******
C
      EXPSAV = CDCOTH(U*XP)            !M*XP*COTH(M*XP)
      CTEMP = YP/(YC*(EXPSAV - 1.D0))  !ZSUBSN/R
      IF(KN.EQ.1) THEN
         YP0= 1.D0/CTEMP          !  PHI=1
      ELSEIF(KN.EQ.2) THEN
         YCP = 1.D0/CTEMP         !  PHI=PHI
      ENDIF                      
      ENDDO
      CTEMP = 1.D0/(AM*YCP + AA*YP0)
      IF(JM.EQ.1) THEN 
          DISTEL = CTEMP            !INTERFACE
      ELSEIF(JM.EQ.0) THEN
          DISTEL = R/(XX + 1.D0/(1.D0 + CTEMP))
      ENDIF
C     
C
      ELSEIF(A18.EQ.8.D0) THEN
      XX = TT*XTX       ! I*OMEGA*TAU 
      YCP = XX**PHI          !S**2PHI  (I*OMEGA*TAU)**PHI
      YC =  AM*YCP + AA*XX
      YP = 1.D0 + YC               !P**2 
      XP = CDSQRT(YP)              !P
C
C     CALCULATE COMPLEX HYPERBOLIC TANGENTS AND FINAL FUNCTIONS.NOTE
C     THAT THE FUNCTION CDCOTH(X) IS DEFINED AS X/TANH(X)  ******
C
      EXPSAV = CDCOTH(U*XP)            !M*XP*COTH(M*XP)
      CTEMP = YP/(YC*(EXPSAV - 1.D0))  !ZSUBSN/R    
      IF(JM.EQ.1) THEN 
         DISTEL = CTEMP            !INTERFACE
      ELSEIF(JM.EQ.0) THEN
         DISTEL = R/(XX + 1.D0/(1.D0 + CTEMP))
      ENDIF
C
      ENDIF
C
C     SHORT-CIRCUIT ANOMALOUS DIFFUSION WHEN EXPONENT MAY NOT BE 1/2
      IF(PX.EQ.1.D0) THEN   
            DISTEL = R/CDCOTH(TS)           !R*(TANH(TS)/TS
      ELSEIF(PX.EQ.2.D0) THEN
          DISTEL = (R*CDCOTH(TS))/(TS*TS)  !R*(COTH(TS)/TS
      ENDIF  
C
C     OPEN-CIRCUIT ANOMALOUS DIFFUSION(PAPER #92, EQ.1
      IF(PX.EQ.3.D0) THEN        
          DISTEL = CTEMP                        !INTERFACE Z
      ELSEIF(PX.EQ.4.D0) THEN     
          DISTEL = R/(XX + 1.D0/(1.D0 + CTEMP)) !TOTAL Z
      ELSEIF(PX.EQ.5.D0) THEN 
          DISTEL = (R/YP)*(1.D0 + 1.D0/(YC*CDCOTH(U*XP))) !TOTAL Z 
C
C     OPEN-CIRCUIT ANOMALOUS DIFFUSION (USES EPSILON, NOT TAU)
      ELSEIF(PX.EQ.6.D0) THEN                             
          DISTEL = R*CTEMP
      ELSEIF(PX.EQ.7.D0) THEN       
          DISTEL = R/(XX + 1.D0/(1.D0 + CTEMP))
      ELSEIF(A18.GE.9.D0) THEN
C
        XX = TT*XTX       ! I*OMEGA*TAU = S**2       S
        YC = XX**PHI          !  (I*OMEGA*TAU)**PHI  U             
        YP = 1.D0 + YC               !P**2           P
        XP = CDSQRT(YP)             
C
C     CALCULATE COMPLEX HYPERBOLIC TANGENTS AND FINAL FUNCTIONS.NOTE
C     THAT THE FUNCTION CDCOTH(X) IS DEFINED AS X/TANH(X)  ******
C
        EXPSAV = CDCOTH(U*XP)            !M*XP*COTH(M*XP)   1/Q
        YP0 = 1.D0/EXPSAV                ! Q
        CTEMP = YP/(YC*(EXPSAV - 1.D0))  !ZSUBSN/R
C
      IF(P18.EQ.9.D0) THEN                 
        DISTEL = R*(YC + YP0)/(YC*(1+XX)+ YP0*(XX -YC))   !EQ.R-2 A
      ELSEIF(P18.EQ.-9.D0) THEN
        DISTEL = YP*YP0/(YC*(1.D0 - YP0))     !EQ.R-2 INTERFACE  AI
      ELSEIF(P18.EQ.1.D1) THEN                 
        DISTEL = R*(XX + YP0)/(XX*YP)          !EQ.B   A.12
      ELSEIF(P18.EQ.-1.D1) THEN
        DISTEL = (1.D0+XX)*YP0/(XX*(1.D0 - YP0))  !EQ.BI INTERFACE
      ELSEIF(P18.EQ.1.1D1) THEN                
        DISTEL = R*(YC + YP0)/(XX*YP)          !EQ.C   LEB 2009
      ELSEIF(P18.EQ.-1.1D1) THEN
        DISTEL = ((YC-XX)+YP0*(1.D0+XX))/(XX*(1.D0 - YP0)) !EQ.C
      ENDIF
C
      ENDIF
C 
      IF(MDE.EQ.-1) THEN
        DISTEL = (1.D0 - (DISTEL/R))/(XTX*CELCAP*R)
      ENDIF
C
      GO TO 5000    
C
C-------------------------------------------------------------------
C
C   WILLIAMS WATTS APPROXIMATE (FRACTIONAL EXPONENTIAL) RESPONSE
C   FOUR ACCURACY CHOICES AVAILABLE
C
1000  XX = TT*DCMPLX(0.D0,OMEGA)
      ZZ = TT*OMEGA
      ZS = ZZ*ZZ
      TOMEGA = ZZ
      PHICOM = PHI
      PP = 1.D0/PHI
C
      IF(ATEMP.GE.0.D0) THEN
        NEPS = 4
      ELSE
        IF(U.LE.0.D0) THEN
          SACC = 1.D-7
          NNP = 22
        ELSEIF(U.GT.0.D0.AND.U.LT.1.D1) THEN
          SACC = 1.D-5
          NNP = 14
        ELSEIF(U.GE.9.99D0) THEN
          SACC = 1.D-3
          NNP = 10
        ENDIF
          IF(PHI.GT.0.6D0) THEN
            NEPS = NNP
          ELSE
            NEPS = 2*NNP - 2
          ENDIF
      ENDIF
C
C   CALCULATE MOMENTS FOR FINAL OUTPUT
      DO IJ = 1,NEPS + 1
        IF(IJ.LT.3) THEN
            AGNN = 0.D0
        ELSE
            AGNN = GAMMLN(DFLOAT(IJ))
        ENDIF
            AGNL = GAMMLN(IJ*PP)
            AGNH = GAMMLN(IJ*PHI)
            XLL(IJ) = AGNL - AGNN
            XLH(IJ) = PHI*DEXP(AGNH - AGNN)
      ENDDO
C
      RN = 0.D0
      AIN = 0.D0
      DO IJ = 1,4
        XXJ(IJ) = PP*DEXP(XLL(IJ))
      ENDDO
        XX1I = 1.D0/XXJ(1)
        XDOI = XX1I/TOMEGA
      IF(MDE.GE.0) THEN
        XXM1 = 1.D100
        XX1 = XXJ(1)
        XX2 = XXJ(2)
        XX3 = XXJ(3)
      ELSE
        XXM1 = XX1I
        XX1 = XXJ(2)*XX1I
        XX2 = XXJ(3)*XX1I
        XX3 = XXJ(4)*XX1I
      ENDIF
C
C   CALCULATE RESPONSE FROM JRM APPROXIMATION METHOD IF ATEMP > 0
C
      IF(ATEMP.GT.0.D0) THEN
        CTNPHI = 1.D0/DTAN(PI2*PHI)
        PHI1 = DSIN(PI2*PHI)*(0.89879 + 0.0878113*PHI)
        PHI2 = CTNPHI*(0.627503 + 0.614423*DEXP(-3.77327*PHI))
        RFAC = 1 + (0.359585 + 34.1304*DEXP(-7.37736*PHI) + 1.86283D5*
     + DEXP(-36.8585*PHI))*CTNPHI
        TFAC = 1 + (1.3672 + 136.604*DEXP(-8.04*PHI) + 1.6615D6*
     + DEXP(-39.3333*PHI))*CTNPHI
        IF (PHI1.LT.0.D0) PHI1 = 0.D0
        IF (PHI2.LT.0.D0) PHI2 = 0.D0
C
        ARG = PI2 * PHI1
        CTEMP = (TFAC*ZZ)**PHI1 * DCMPLX(DCOS(ARG), DSIN(ARG)) + 1.D0
        ARG = DATAN(DIMAG(CTEMP)/DREAL(CTEMP))*PHI2
        CTEMP = ABS(CTEMP)**PHI2 * DCMPLX(DCOS(ARG), DSIN(ARG))
C
        XTX = RFAC/CTEMP    
        YP0 =  1.D0/(1.D0 + XX*XTX)
        IF(MDE.LT.0) THEN
          YP1 = XX1I*XTX*YP0
          DISTEL = R*YP1
        ELSE
          DISTEL = R*YP0
        ENDIF
C
        GOTO 5000
      ENDIF
C   END OF ATEMP > 0
C
C 291   FORMAT(I4,1P,(1E12.4),1P,(4E19.11))
C 167   FORMAT(I5,1P,(2E19.11),1P,(2E15.7))
C
C   CALCULATE LFREQ RESPONSE FROM SERIES, EPSILON ALGORITHM
C
469   CONTINUE
      IF(PHI.GE.0.5D0) THEN
        ZLH = 1.5D0*(PHI**5.41)
      ELSE
        ZLH = 1.653D1*(PHI**8.91)
      ENDIF
C
      IF(ZZ.LT.ZLH) THEN
      ALP = 0.D0
      ALRI = 0.D0
      ICR = 0
      ICRE = 0
      ICI = 0
      ICIE = 0
      DO KJ = 1,NEPS
      IF(MOD(KJ,2).EQ.0.AND.ICRE.EQ.0) THEN
        ICR = ICR + 1
        ISGR = (-1)**ICR
        TT = XLL(KJ) + (ICR - 1)*DLOG(ZS)
        DELP = -ISGR*PP*DEXP(TT)
        ALP = ALP + DELP
        XXJ(ICR) = ALP
        ARDD = DABS(DELP/ALP)
        IF(ARDD.LT.SACC) THEN
            ICRE = 1
        ELSE
            ICRE = 0
        ENDIF
      ELSEIF(MOD(KJ,2).NE.0.AND.ICIE.EQ.0) THEN
        ICI = ICI + 1
        ISGI = (-1)**ICI
        TI = XLL(KJ) + (ICI - 1)*DLOG(ZS)
        DELI = -ISGI*PP*DEXP(TI)
        ALRI = ALRI + DELI 
        GSP(ICI) = ALRI
        AIDD = DABS(DELI/ALRI)
        IF(AIDD.LT.SACC) THEN
            ICIE = 1
            GOTO 349
        ELSE
            ICIE = 0
        ENDIF
      ENDIF
349   ENDDO   
        IF(ICRE.EQ.0) CALL EPSALG(ICR,XXJ,ALP)
        ALRO = -ZS*ALP + 1.D0
        IF(ICIE.EQ.0) CALL EPSALG(ICI,GSP,ALRI)
        ALIO = -ZZ*ALRI
        IF(MDE.GE.0.D0) THEN
        DISTEL = R*DCMPLX(ALRO,ALIO)
        ELSEIF(MDE.GT.-3) THEN
        ALR1 = ALRI*XX1I
        ALI1 = -ZZ*ALP*XX1I
        DISTEL = R*DCMPLX(ALR1,ALI1)
        ELSEIF(MDE.EQ.-3) THEN
        DISTEL = (1.D0 - DCMPLX(ALRO,ALIO))/
     +          (R*CELCAP*DCMPLX(0.D0,OMEGA))
       ENDIF
      ENDIF
C     END OF LFREQ ZZ CHOICE
C     ********************************************
C
C   CALCULATE HFREQ RESPONSE FROM SERIES, EPSILON ALGORITHM
C
      IF(ZZ.GE.ZLH) THEN
      ICR = 0
      ICI = 0
      ICRE = 0
      ICIE = 0
      AHRO = 0.D0
      AHIO = 0.D0
      DO KJ = 1,NEPS
        ISG = (-1)**(KJ-1)
        ARG = KJ*PI2*(PHI + 1.D-7)
        OMP = ISG*(ZZ**(-KJ*PHI))
        OMPG = OMP*XLH(KJ)
        CSN = DCOS(ARG)
        SSN = DSIN(ARG)
        IF(DABS(CSN).GT.SML.AND.ICRE.EQ.0) THEN
            AHRO = AHRO + OMPG*CSN
            ICR = ICR + 1
            XXJ(ICR) = AHRO
            IF(DABS(OMPG/AHRO).LT.1.D-7) THEN
                ICRE = 1
            ELSE
                ICRE = 0
            ENDIF
        ENDIF
C 394   FORMAT(3I5,4E15.7)
        IF(DABS(SSN).GT.SML.AND.ICIE.EQ.0) THEN
            AHIO = AHIO - OMPG*SSN
            ICI = ICI + 1
            GSP(ICI) = AHIO
            IF(DABS(OMPG/AHIO).LT.1.D-7) THEN
                ICIE = 1
                GOTO 449
            ELSE
                ICIE = 0
            ENDIF
        ENDIF
449   ENDDO
C
      IF(ICRE.EQ.0) CALL EPSALG(ICR,XXJ,AHRO)
      IF(ICIE.EQ.0) CALL EPSALG(ICI,GSP,AHIO)
      IF(MDE.GE.0.D0) THEN
        DISTEL = R*DCMPLX(AHRO,AHIO)
      ELSEIF(MDE.GT.-3)THEN
        AHR1 = -AHIO*XDOI
        AHI1 = -(1.D0 -AHRO)*XDOI
        DISTEL = R*DCMPLX(AHR1,AHI1)
      ELSEIF(MDE.EQ.-3) THEN
        DISTEL = (1.D0 - DCMPLX(AHRO,AHIO))/
     *          (R*CELCAP*DCMPLX(0.D0,OMEGA))
      ENDIF
      ENDIF
C
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   JONSCHER RESPONSE (GENERALIZED TO INCLUDE REAL PART)
C
 1100 ARG1 = PI2*PHI
      ARG2 = PI2*U
      TOMEGA = T*OMEGA
      CTEMP = DCMPLX(1.D0/DTAN(ARG1), -1.D0) / TOMEGA**PHI
      DISTEL = DCMPLX(1.D0/DTAN(ARG2), -1.D0) / TOMEGA**U
      DISTEL = R*(DISTEL + CTEMP)
      GO TO 5000
C
C-------------------------------------------------------------------
C
C   EXPON. DISTRIBUTION OF ACTIVATION ENERGIES, ASYMMETRIC CASE EDAE1
C                             U = LN(R2)
C
C   SET UP VARIABLES TO PASS DATA TO DAEI FUNCTIONS THRU /CM5/
C
C   IF T<0, USE |T|= TAUA, R=E, AND ATEMP AS INPUTS: RHO,SIGMA DATA
C   XS = U; INPUT OR CALC (WHEN T<0)
C
1200  IF (TT.LT.0) THEN 
        TA = -TT
        TMP = DABS(ATEMP)
        XS = AKINV*R/TMP
        RR = DEXP(XS)
        ZT = TA*EVI     
        TAUO = TA*RR        
        R0 = ZT*(RR - 1.D0)/XS
        RINF = ZT*XS*(1.D0 - (1.D0/RR))     
        RX = R0 - RINF 
        IF(U.LT.0) THEN
            PHIN = 1.D0 - PHI/XS
        ELSE
            PHIN = PHI
        ENDIF
        TAU = DABS(U)*TAUO
        UU = XS
        TOMEGA = TAU*OMEGA
        F5 = OMEGA/TWPI
        IF(F5.GT.0.999D5.AND.F5.LT.1.001D5) THEN
            WRITE(*,39) TMP,XS,R0,RINF,TAUO,PHIN
        ENDIF
      ELSE
        TOMEGA = TT*OMEGA
        RX = R
        RINF = 0.D0
        PHIN = PHI
        UU = DABS(U)
      ENDIF
      IF (PHIN.EQ.0.D0) THEN
        RP = 1.D0 / UU
      ELSE
        PHU = PHIN*UU
        EPH = DEXP(-PHU)
        IF(DABS(PHU).GT.3.D2) PHU = DSIGN(1.D2,PHU)
        RP = PHIN / (1.D0 - EPH)
      END IF
39    FORMAT(1X,9(1PE13.4))
C
C   COMPUTE INTEGRALS
C
C     REAL:
      IV = 1
      PHICOM = PHIN
      CALL QROMB(DAEFN,0.D0,UU,Z1)
C
C     IMAGINARY:
      PHICOM = PHIN + 1.D0
      CALL QROMB(DAEFN,0.D0,UU,Z2)
C
      DISTEL = RINF + RX * RP * DCMPLX(DBLE(Z1), -TOMEGA*DBLE(Z2))
      GOTO 5000
C
C-------------------------------------------------------------------
C
C   EXPON. DISTRIBUTION OF ACTIVATION ENERGIES, SYMMETRIC CASE EDAE2
C                               U = LN(R1)
C
1300  U = DABS(U)
      IF (PHI.EQ.0.D0) THEN
        RP = 0.5D0/U
      ELSE
        PHU = PHI*U
        IF(DABS(PHU).GT.1.D2) PHU = DSIGN(1.D2,PHU)
        RP = 0.5D0*PHI / (1.D0 - DEXP(-PHU))
      END IF
C
C   SET UP VARIABLES TO PASS DATA TO EDAE FUNCTIONS THRU /CM5/
C
      TOMEGA = TT*OMEGA
C
C   COMPUTE INTEGRALS
C
C     REAL:
      IV = 0
      PHICOM = PHI
      CALL QROMB(DAEFN2,0.D0,U,Z1)
C
C     IMAGINARY:
      IV = -2
      PHICOM = PHI + 1.D0
      CALL QROMB(DAEFN2,0.D0,U,Z2)
C
      DISTEL = R * RP * DCMPLX( DBLE(Z1), -TOMEGA*DBLE(Z2) )
      GOTO 5000
C
C-------------------------------------------------------------------
C
C   GAUSSIAN DISTRIBUTION OF ACTIVATION ENERGIES, SYMMETRIC CASE
C
1400  CONTINUE
      XIGDAE = PHI              !HERE PHI IS XI
C
C      CALCULATE NORMALIZATION COEFFICIENT, RP
      ICHG = 0
      TCMEGA = 0.D0
      UPPER1 = U
      BOTTM = -U
      CALL QROMB(GDAEFN,BOTTM,UPPER1,VINT)
      RP = 1.D0/VINT
C
C   DO OMEGA VARIABLE IS T VARIABLE: |MODE| =|MDE|= 8 FOR TRANSIENT
      IF(ABS(MDE).EQ.8) THEN
        IWT = 0
      ELSE
        IWT = 1
      ENDIF
C
      IF(IWT.EQ.0) THEN 
        TCMEGA = OMEGA/T
      ELSE
        TCMEGA = TT*OMEGA   
      ENDIF
C         TOMEGA = T*OMEGA
C
C   SET UP VARIABLES TO PASS DATA TO GDAE FUNCTION THRU /CM9/
C
C   COMPUTE INTEGRALS
C
C     REAL:
      ICHG = 0
      CALL QROMB(GDAEFN,BOTTM,UPPER1,VINT)
      Z1 = VINT
C
      IF(IWT.EQ.0) GOTO 191
C     IMAGINARY:
      ICHG = 1
      CALL QROMB(GDAEFN,BOTTM,UPPER1,VINT)
      Z2 = VINT
C
191   DISTEL = R * RP * DCMPLX(DBLE(Z1), -IWT*TCMEGA*DBLE(Z2))
      GOTO 5000
C
C-------------------------------------------------------------------
C
C   GENERAL HOMOGENEOUS DIFFUSION DCE: MACROSCOPIC PARAMETERS
C
1500  TD = TT
      IF(U.LT.0.D0) U = -U
      IF(U.EQ.1.D0) THEN
        QD = R
        PD = PHI
      ELSEIF(U.EQ.2.D0) THEN
        CD = R
        IF(CD.EQ.0) THEN
          WRITE(*,40)
40        FORMAT(5X,'WRONG INPUT')
        ENDIF
        QD = TD/CD
        PD = TD*PHI
      ELSEIF(U.GE.2.999D0) THEN
        QD = R
        PD = TD*PHI
      ENDIF        
        TOMEGA = TD*OMEGA
        CTEMP = DCMPLX(0.D0,TOMEGA)
        XX = SQRT(CTEMP)
C   DEFINITION OF XTX: X*TANH(X)
        XTX = CMTANH(XX)
        DISTEL = QD*(1.D0 +(PD*XTX/CTEMP))/(XTX + PD)
      IF(MDE.NE.-1) GOTO 5000
C
      IF(U.EQ.1.AND.MDE.EQ.-1) THEN 
C       IF MODE = -1,  R (INPUT) BECOMES EPSSUBCO
        YC =  DISTEL / R
        YP = DCMPLX(0.D0,OMEGA)
C
        DISTEL = (1.D0 - YC)/(YP*CELCAP*R)  !  R IS CCO HERE
      ENDIF
C
      GO TO 5000
C
C-------------------------------------------------------------------
C
C   GENERAL HOMOGENEOUS DIFFUSION DCE: MICROSCOPIC PARAMETERS
C
1600   AA = R
       BD = T
       SD = DABS(U)
      AK1 = PHI
      RHO = AK1*SD/BD
C
      TOMEGA = OMEGA*SD*SD/BD
      CTEMP = DCMPLX(0.D0,TOMEGA)
         XX = SQRT(CTEMP)
C        DEFINITION: XTX = X*TANH(X)
         XTX = CMTANH(XX)
      DISTEL = AA*(1.D0 +(RHO*XTX/CTEMP))/((BD/SD)*XTX +AK1)
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   IDEAL ELEMENTS: R,C,R2, AND INDUCTANCE PL (SEE DESCRIPTION ABOVE)
C
1700   C = T
      R2 = U
      PL = PHI
      IF(TC) THEN
        YC = (0,0)
      ELSE
        YC = DCMPLX(0.D0,OMEGA*C)
      ENDIF
      IF(RRC) THEN
        YP = (0,0)
      ELSE
        YP = DCMPLX(1.D0/R,0.D0)
      ENDIF
        YCP = YC + YP
        XTX = DCMPLX(0.D0,OMEGA*PL) + R2
      IF(XTX.EQ.(0,0).AND.YCP.NE.(0,0)) THEN
        DISTEL = 1.D0/YCP
      ELSE
        DISTEL = XTX/(1.D0 + YCP*XTX)
      ENDIF
      GOTO 5000
C
C-------------------------------------------------------------------
C
C   DISSADO-HILL FITTING ELEMENT, Z LEVEL
C
1800  CALL DISSHILL(R,TT,U,PHI,DHRL,DHIM,OMEGA)
      CTEMP = DCMPLX(DHRL,DHIM)
      DISTEL = 1.D0/(CTEMP*DCMPLX(0.D0,OMEGA*CELCAP))
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   DISSADO-HILL FITTING ELEMENT, EPSILON LEVEL
C
1900  CALL DISSHILL(R,T,U,PHI,DHRL,DHIM,OMEGA)
      DISTEL = DCMPLX(DHRL,DHIM)
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   POWER LAW ELEMENT #1
C       Y = RDE*(OMEGA**TDE) + I UDE*(OMEGA**PDE)
C
2000  YC =  DCMPLX(R*(OMEGA**T), U*(OMEGA**PHI))
C   DISTEL = 1.D0/YC
      DISTEL = YC                !Y LEVEL
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   POWER LAW ELEMENT #2
C OLD       Y = RDE*(OMEGA**TDE) +  UDE*[(I OMEGA)**PDE]
C       Y =  UDE*[(I OMEGA)**PHI]   Y IS EPSILON'' HERE
C
C 2100   COMEG = DCMPLX(0.D0,OMEGA)
C  YC =  R*(OMEGA**T) + U*(COMEG**PHI)
C
2100  CONTINUE
      IF(T.NE.0) THEN
        PHIM = PHI - 1.D0   !PHI=N, NEARLY 1
      ELSE
        PHIM = - PHI
      ENDIF
C
      CY = DSIN(PI2*PHIM)
      DIM = CY*U*(OMEGA**PHIM)
      DISTEL = DCMPLX(1.D-90, DIM)     !EPSILON LEVEL; PHI~1
C
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   POWER LAW ELEMENT #3
C       Y = RDE*[(I OMEGA)**TDE] +  UDE*DCMPLX(0.D0,OMEGA**PDE]
C
2200  COMEG = DCMPLX(0.D0,OMEGA)
      YC =  R*(COMEG**T) + U*DCMPLX(0.D0,OMEGA**PHI)
      DISTEL = 1.D0/YC
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   LADDER NETWORK #1
C   NOW INFINITE LADDER WITH EQUAL R'S AND C'S (USE T).  USE 
C   U FOR SCALE
C
2300  CONTINUE
      COMEG = DCMPLX(0.D0,T*OMEGA)
      YC = 0.5D0*U*COMEG*(1.D0 + (1.D0 + 4.D0/(R*COMEG))**PHI) 
      DISTEL = 1.D0/YC
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   LADDER NETWORK #2
C   NOW INFINITE LADDER WITH EQUAL R'S AND C'S (USE T).  USE 
C   U FOR SCALE
C
2400  CONTINUE
      COMEG = DCMPLX(0.D0,T*OMEGA)
      YC = 0.5D0*U*(1.D0 + (1.D0 + 4.D0*R*COMEG)**PHI)/R 
      DISTEL = 1.D0/YC
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   EFFECTIVE MEDIUM RESPONSE USING PCPE
C   
2500  COMEG = DCMPLX(0.D0,OMEGA)
      YP1 = T*(COMEG**(-PHI))     !E1
      E0 = R
CC  FF = U              !XSUBC?
C
      YC = 3.D0*U*(YP1 - R)/(YP1 + 2.D0*R - U*(YP1 - R))
      DISTEL = R*(1.D0 + YC)
C
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   POWER LAW ELEMENT #7  SCPE FOR ELECTRODE EFFECTS.  Z LEVEL (SERIES)
C
2600  COMEG = DCMPLX(0.D0,OMEGA)
      DISTEL = 1.D0/((CELCAP*U)*(COMEG**(PHI)))   !Z LEVEL FOR CSD
C
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   POWER LAW # 8: PCPE FOR NEARLY CONST. LOSS. DEF. AT EPSILON LEVEL,DSD 
C   PARALLEL
C
2700  COMEG = DCMPLX(0.D0,OMEGA)
      YC = U*(COMEG**(-PHI))      !EPS LEVEL FOR DSD
      DISTEL = YC
      GO TO 5000
C
C--------------------------------------------------------------------
C   POWER LAW # 9: DISTEL IS COMPLEX DIELECTRIC CONSTANT HERE.
C   USE ONLY IN OSUB, DSD PART (E.G., P40 = 4, R=P(15),ETC)
C
C     FOR NEARLY-CONSTANT-LOSS(NCL) OR CONSTANT LOSS (IN EPS''),
C     SET U=O AND PHI=0 AND THEN USE R AND T PARAMETERS. SET T=0
C     FOR CONSTANT LOSS. ONLY APPROPRIATE FOR SIGMA' FITS
C
2800  CONTINUE
      IF(U.EQ.0.D0.AND.PHI.EQ.0.D0) THEN
        YC = DCMPLX(0.D0, - R)
        XXX = OMEGA**(T - 1.D0)     
        DISTEL = 1.D-80 + YC*XXX   !SO: Y= R*CELCAP*(OMEGA**T)
      ELSE
        COMEG = DCMPLX(0.D0,OMEGA)
        DISTEL = R*(OMEGA**T) + U*(COMEG**(-PHI))   !DEFINED AT EPS LEVEL
      ENDIF
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   GENERALIZED DAVIDSON-COLE RESPONSE 
C
2900  CONTINUE
C
      YC = (1.D0 + U)/(1 + U*((1.D0 + 
     + DCMPLX(0.D0,OMEGA*TT))**PHI))
C
      XTX = DCMPLX(0.D0,OMEGA)
C
      IF(MDE.GE.0) THEN
          DISTEL = R*YC
      ELSEIF(MDE.EQ.-1) THEN              
C       IF MODE = -1,  R (INPUT) BECOMES EPSSUBCO
        DISTEL = (1.D0 - YC)/(XTX*CELCAP*R) !  R IS CCO HERE
      ELSE
        DISTEL = (0.D0,0.D0)
        RETURN
C        STOP
      ENDIF
C
      GO TO 5000
C
C--------------------------------------------------------------------
C   EXPON. DISTRIBUTION OF ACTIVATION ENERGIES, ASYMMETRIC CASE EDAE1
C                             U =-LN(R2)
C    EXACT LOG EXPRESSIONS FOR PHI = 1 OR 0 (-.01<=PHI<=.01)
C    ******  PHI MUST BE FIXED *****
C
3000  CONTINUE
      IF(PHI.LT.0) THEN
        SUB = 1.D0
      ELSE
        SUB = 0.D0
      ENDIF
        IF(U.EQ.0.D0) THEN
          DISTEL = (0.D0,0.D0)
          RETURN
C        STOP
      ELSE
        IF(DABS(U).GT.3.0D2) U =-DSIGN(3.0D2,U)
      ENDIF
      EPH = DEXP(U)
      IF(DABS(PHI).EQ.1.D0) THEN
        RP = 1.D0 / (1.D0 - EPH)
      ELSE
        RP = 1.D0/(-U)
      ENDIF
C
      TOMEGA = TT*OMEGA
      YP = DCMPLX(0.D0,TOMEGA)
C
      XX = 1.D0 + YP
      XTX = 1.D0 + YP*EPH
C   GET DISTEL IN NORMALIZED Z-LEVEL FORM
      IF(DABS(PHI).EQ.1.D0) THEN
        DISTEL = RP*LOG(XX/XTX)/YP
      ELSE
        EPB = 1.D0/EPH
        DISTEL = RP*LOG(EPB*XTX/XX)
      ENDIF
C
      IF(PHI.LT.0.D0) THEN
        DISTEL = 1.D0/DISTEL - SUB
        DISTEL = 1.D0/DISTEL
      ENDIF
        DISTEL = R*DISTEL
      GOTO 5000
C
C--------------------------------------------------------------------
C          *************** GBEM EMA APPROXIMATE FITTING FUNCTION
C   FIT EMA DATA TO APPROX GBEM/BDM RESPONSE FUNCTION
C   JRM 11/20/93;12/21/93
C
C   PHI < 0 : USE MICRO QUANTITIES; OTHERWISE MACRO
C   PHI=2: EPCE; PHI=4 OR 8: TAU
C   HERE EPCE IS EPSILON CONTRIB FROM CONDUCTIVE SYSTEM
C   PHI MUST BE +-2, +-4, 0R +-8; |PHI|=8: SUBTRACTION OF G0N
C   IF CELCAP = EPSV, (VACUUM PERMITTIVITY) THEN WE DEAL WITH 
C       SIGMA AND RHO INSTEAD OF Y AND Z
C   R: R0 (PHI>0) OR TAUA=TA (PHI<0)
C   U: EC
C   SEE BEGINNING OF SUBROUTINE FOR DEFINITIONS OF MANY DATA QUANTITIES
C
3100  IF(U.LT.0.D0) THEN
        XC = -U
      ELSE
        XC = PE*U/DABS(ATEMP)
      ENDIF
      IF(XC.GT.3.D2) THEN
          XCL = 3.D2
      ELSE
          XCL = XC
      ENDIF
      EXC = DEXP(XCL)
      EXCI = P0/EXC
      RQ = P02*EXC/(P0 + EXCI)
      RINN = P0/RQ
      NPH = NINT(DABS(PHI))
C
      IF(NPH.NE.2.AND.NPH.NE.4.AND.NPH.NE.8) THEN
        DISTEL = (0.D0,0.D0)
        RETURN
C        STOP
      ENDIF
C
C   PHI < 0: USE TA,EC.
      IF(PHI.LT.0.D0) THEN
        SFT = 0.253D0*PI3
        EPCE = SFT*XC
        CEP0 = EPCE*CELCAP
        TA = R
        ZTA = TA/CELCAP
        R0 = ZTA*RQ
        TAU = R0*CEP0
        IF(T.LT.0.D0) TAU = TAU*DABS(T)
      ELSE
C
C   IF PHI > 0: USE R0,EC, AND T = EPCE (NPH=2) OR TAUE (NPH=4 OR 8)
C
        IF(T.LT.0.D0) THEN
          DISTEL = (0.D0,0.D0)
          RETURN
C        STOP
        ENDIF
        R0 = R
        IF(NPH.EQ.2) THEN
            EPCE = T
            CEP0 = CELCAP*EPCE
            TAU = R0*CEP0
        ELSEIF(NPH.EQ.4.OR.NPH.EQ.8) THEN
            TAU = T
            CEP0 = TAU/R0
            EPCE = CEP0/CELCAP          
        ENDIF
      ENDIF
C
C   CALCULATE AND USE NORMALIZED FREQUENCY FOR EMA
      FW = TAU*OMEGA
C
C   CALC MN(OMEGA) AND ZN (==XX) EMA RESPONSE
C
      CALL GBNM(XC,FW,RINN,NPH,XX)
C
C   RETURN UNNORMALIZED IMPEDANCE (NOT MODULUS) VALUES
C
      DISTEL = R0*XX
C
      GOTO 5000
C
C--------------------------------------------------------------------
C
C  "EXACT" WILLIAMS WATTS (FRACTIONAL EXPONENTIAL) RESPONSE: BETA = 0.5
C       RESPONSES CALCULATED FROM EXACT DRT WITH CUTOFF
C   IF MODE=MDE < 0, THEN CALC MOYNIHAN/JRM FORM DIRECTLY: CSD1 
C   IF MODE=MDE >= 0, THEN CALC GENERAL CSD=CSD0 OR DSD FORM DIRECTLY
C   IF MODE=MDE = 8, THEN CALC CSD=CSD0 OR DSD TRANSIENT RESPONSE
C   IF MODE=MDE =-8, THEN CALC CSD=CSD1 "TRANSIENT" RESPONSE
C   IF MODE=MDE = 16, THEN CALC COUPLING MODEL RESPONSE; USE U=-50 FIXED
C
C   CUTOFF DETERMINED BY VALUE OF U; *** IT IS USUALLY NEGATIVE
C   WHEN U <-40 OR SO, THERE IS NO EFFECT IN ANY PRACTICAL MEASURED
C   FREQUENCY RANGE. 
C
C   THE INPUT VALUE OF PHI IS IMMATERIAL HERE, BUT IT MUST BE FIXED
C
C   COMPUTE NORMALIZATION, INTEGRALS
C
3200  IF(MDE.EQ.16.AND.PX41.NE.0.D0) GOTO 397
      IV = 1
      IF(ZN.EQ.0.D0) THEN
        TOMEGA = 0.D0
            CALL QROMB(WW5,U,1.6D1,ZN)
        RN = 1.D0/ZN
        ANM = RN
C
C   DO XMIN AND ITS CONSEQUENCES
        XMIN = DEXP(U)
        XM4 = XMIN/4.D0
        DO KJ = 1,5
          MM = KJ - 1
          AGG = DFLOAT(MM) + 0.5D0
          GAMM = DEXP(GAMMLN(AGG))
          XXJ(KJ) = (4.D0**MM)*GAMM*SQPII*GAMMQ(AGG,XM4)/
     +               GAMMQ(0.5D0,XM4)
C
        ENDDO
C
        IF(MDE.LT.0) THEN
          XXM1 = 1.D0/XXJ(2)
          XX1 = XXJ(3)*XXM1
          XX2 = XXJ(4)*XXM1
          XX3 = XXJ(5)*XXM1
        ELSE
          XMS = DSQRT(XMIN)
          XXM1 = SQPII/XMS - 0.5D0 + PII  
          XX1 =  XXJ(2)
          XX2 =  XXJ(3)
          XX3 =  XXJ(4)
        ENDIF
          AIN = XXJ(1)
C
      ENDIF
C
      AMDE = ABS(MDE)
      IF(AMDE.NE.8) THEN
        TOMEGA = TT*OMEGA
      ELSE
        TOMEGA = OMEGA/TT
      ENDIF
C
C     REAL:
      IV = 1
      CALL QROMB(WW5,U,1.6D1,Z1)
C   TRANSIENT RESPONSE CHOICE BELOW
      IF(AMDE.EQ.8) GOTO 393
C
C     IMAGINARY:
      IV = 3
      CALL QROMB(WW5,U,1.6D1,Z2)
C
      GOTO 567
        XMIN = DEXP(U)
        XMS = DSQRT(XMIN)
        XM4 = XMIN/4.D0
        XEX = DEXP(-XM4)
        ANX = 1.D0/(1.D0 - SQPII*XMS*XEX) 
C
        GSP(1) = -2.D0
        GSP(2) = 1.D0
        GSP(3) = 0.5D0
        GSP(4) = 0.75D0
        GSP(5) = 1.5D1/8.D0
        GSP(6) = 1.05D2/1.6D1
        GSP(7) = 9.45D2/3.2D1
C
      DO IJ = 1,6
        FOM = 4.D0**MM
        AM1 = MM + 0.5D0
        OXMN = 1.D0 -(XM4**AM1)*SQPII/GSP(IJ+1)
        XXJ(IJ) = FOM*GSP(IJ)*ANX*OXMN
      ENDDO
C
      IF(MDE.LT.0) THEN
        XXM1 = 1.D0/XXJ(3)
        XX1 = XXJ(4)*XXM1
        XX2 = XXJ(5)*XXM1
        XX3 = XXJ(6)*XXM1
      ELSE
        XXM1 = XXJ(1)
        XX1 =  XXJ(3)
        XX2 =  XXJ(4)
        XX3 =  XXJ(5)
      ENDIF
        AIN = XXJ(2)
567   CONTINUE
C
      DISTEL = R*DCMPLX(DBLE(Z1),-TOMEGA*DBLE(Z2))*ANM
      GO TO 5000
C
C   DO TWO TYPES OF TRANSIENT RESPONSE
C
393   DISTEL = R*ANM*DCMPLX(DBLE(Z1),0.D0)    !USES NELEM=32 
C
      IF(NCH.EQ.5) THEN
        IF(OMEGA.EQ.PX41.AND.PX45.EQ.0.D0) THEN
          TAO = PX41/(-DLOG(DREAL(DISTEL)))   
          WRITE(*,*) NCH,MDE,TAO
          PX45 = TAO
C          READ(*,*)
          DISTEL = (0.D0,0.D0)
          RETURN
C          STOP
        ELSEIF(PX45.NE.0) THEN
          IF(OMEGA.LT.PX41) THEN
              RARG = DEXP(-OMEGA/PX45)
              DISTEL = R*DCMPLX(DBLE(RARG),0.D0)
          ENDIF
          GOTO 5000
        ENDIF
        GOTO 5000
      ELSEIF(NCH.EQ.4) THEN
        GOTO 5000
      ENDIF
C
C   COUPLING MODEL:DSD, CSDO
397   CONTINUE
        IF(OMEGA.GE.PX41) THEN      !PX41 = P(41) IS TAUC
          TOMEGA = OMEGA/TT   !T/TAUS; OMEGA IS TIME
          RARG = DEXP(-(TOMEGA**PHI))
        ELSE
          TAO = PX41*((TT/PX41)**PHI)
          IF(PX45.NE.0.D0) THEN   !PX45 = P(45) IS GAMX
              BTM = 1.D0 - PHI
              TAO = TAO*DEXP(-BTM*PX45)/PHI
          ENDIF
          RARG = DEXP(-OMEGA/TAO)
        ENDIF
      DISTEL = R*DCMPLX(DBLE(RARG),0.D0)
      GOTO 5000
C
C--------------------------------------------------------------------
C
C   FIT  ARRHENIUS TEMP DATA (Y OR Z LEVEL) OMEGA = TEMPERATURE (K)
C   R=E (EV), T,PHI,U ARE POSSIBLE INPUT PARAMETERS; ATEMP SELECTS 
C       FORM;   DATTYP MUST BE REAL ***
C   IF CELCAP=CG, THEN R0 OR Y0; IF CELCAP=EPSVAC, THEN RHOO, OR 
C   SIGMAO  FOR DISTEL OUTPUT
C   ALSO, WHEN ATEMP>16, USE Z' = U + R (TEMP)**PHI
C
3300  TEMP = OMEGA
      IF(ATEMP.GT.17.D0) GOTO 371
      TEMPU = TEMP - U                          ! (T - T0)
      IF(TEMPU.LE.0) THEN
        DISTEL = (0.D0,0.D0)
        RETURN
C        STOP
      ENDIF
C
       EK = R*AKINV                  ! E/K
       EE = EK/TEMP                  ! E/(KT)
      EEU = EK/TEMPU                 ! E/[K(T - T0)]
      TAV = T/CELCAP
C
      IF(DABS(EE).GT.2.D2) EE = DSIGN(2.D2,EE)
      IF(DABS(EEU).GT.2.D2) EEU = DSIGN(2.D2,EEU)
      EXT = DEXP(EE)
      EXTU = DEXP(EEU)
C
      IF(ATEMP.NE.-6.4D1) CTEMP = DCMPLX(0,1.D3/OMEGA) ! I1000/TEMPERATURE
C
371   CONTINUE
      IF(ATEMP.GT.1.5D3) THEN
        DISTEL = U +R*(TEMP**PHI) + CTEMP
      ELSEIF(ATEMP.EQ.-64.D0) THEN
        DISTEL = TAV*(OMEGA**PHI)*EXTU + CTEMP
      ELSEIF(ATEMP.EQ.-4.D0) THEN
        DISTEL = PHI + T*EXTU + CTEMP   
      ELSEIF(ATEMP.EQ.-2.D0) THEN
        DISTEL = PHI + TAV*EXTU + CTEMP 
      ELSEIF(ATEMP.EQ.-1.D0) THEN
        DISTEL = TAV*(OMEGA**PHI)*EXTU + CTEMP  
      ELSEIF(ATEMP.EQ.0.D0) THEN
        DISTEL = TAV*(EE**(-PHI))*EXTU + CTEMP  
      ELSEIF(ATEMP.EQ.1.D0) THEN
        DISTEL = TAV*(EE**(-PHI))*(EXT - 1.D0) + U + CTEMP  
      ELSEIF(ATEMP.EQ.2.D0) THEN
        DISTEL = 1.D0/(TAV*(OMEGA**PHI)*EXTU) + CTEMP   
      ELSEIF(ATEMP.EQ.4.D0) THEN
        DISTEL = 1.D0/(TAV*(EE**(-PHI))*EXTU) + CTEMP   
      ELSEIF(ATEMP.EQ.8.D0) THEN
        DISTEL = U + 1.D0/(TAV*(EE**(-PHI))*(EXT - 1.D0))+ CTEMP
      ELSEIF(ATEMP.EQ.16.D0) THEN
        DISTEL = -R*DLOG(TEMP) - PHI*(1.D0 - TEMP) + CTEMP
      ELSEIF(ATEMP.GT.1.61D1) THEN
        DISTEL = U + R*(TEMP**PHI) + CTEMP
      ENDIF
C
      GO TO 5000
C
C--------------------------------------------------------------------
C
C   GENERALIZED EXPONENTIAL DISTRIBUTION OF ACTIVATION ENERGIES
C
3400  XIGDAE = PHI                 ! HERE PHI IS OLD PHI
C
C      CALCULATE NORMALIZATION COEFFICIENT, RP
      ICHG = 0
      TCMEGA = 0.D0
      IF(U.GE.0) THEN                        !SYMMETRIC CASE
         UPPER1 = U
         BOTTM = -U
      ELSE                                   !ASYMMETRIC CASE, -|U| TO 0
         UPPER1 = 0
         BOTTM = U
      ENDIF
C
      CALL QROMB(GDAEFN2,BOTTM,UPPER1,VINT)
        RP = 1.D0/VINT
C
C   DO OMEGA VARIABLE IS T VARIABLE: |MODE| = |MDE|= 8 FOR TRANSIENT
      IF(ABS(MDE).EQ.8) THEN
        IWT = 0
      ELSE
        IWT = 1
      ENDIF
C
      IF(IWT.EQ.0) THEN 
        TCMEGA = OMEGA/T
      ELSE
        TCMEGA = TT*OMEGA   
      ENDIF
C
C   SET UP VARIABLES TO PASS DATA TO GDAE FUNCTION THRU /CM9/
C
C   COMPUTE INTEGRALS
C
C     REAL:
      ICHG = 0
      CALL QROMB(GDAEFN2,BOTTM,UPPER1,VINT)
      Z1 = VINT
C
      IF(IWT.EQ.0) GOTO 192
C     IMAGINARY:
      ICHG = 1
      CALL QROMB(GDAEFN2,BOTTM,UPPER1,VINT)
      Z2 = VINT
C
192   DISTEL = R * RP * DCMPLX(DBLE(Z1), -IWT*TCMEGA*DBLE(Z2))
      GOTO 5000
C
C--------------------------------------------------------------------
C
C  "EXACT" WILLIAMS WATTS (FRACTIONAL EXPONENTIAL) RESPONSE: BETA = 1/3
C       RESPONSES CALCULATED FROM EXACT DRT WITH CUTOFF
C   IF MODE=MDE < 0, THEN CALC MOYNIHAN/JRM FORM DIRECTLY: CSD1 
C   IF MODE=MDE >= 0, THEN CALC GENERAL CSD=CSD0 OR DSD FORM DIRECTLY
C   MODE = +8,-8: TRANSIENT RESPONSE
C   MODE = 16: COUPLING MODEL
C   CUTOFF DETERMINED BY VALUE OF U; *** IT WILL USUALLY BE NEGATIVE
C   WHEN -U| >= 40 OR SO, THERE IS NO EFFECT IN ANY PRACTICAL MEASURED
C   FREQUENCY RANGE.  U = - LN(TAUS/TAUC) 
C
C   THE INPUT VALUE OF PHI IS IMMATERIAL HERE
C
C   COMPUTE NORMALIZATION, INTEGRALS
C
3500  IF(MDE.EQ.16.AND.PX41.NE.0.D0) GOTO 797
      IV = 0
      IF(ZN.EQ.0.D0) THEN
        TOMEGA = 0.D0
        TOP =  16.D0
        IF(TOP.LE.U) THEN
          DISTEL = (0.D0,0.D0)
          RETURN
        ENDIF
        CALL QROMB(AIFU,U,TOP,ZN)
        RN = 1.D0/ZN
        ANM = RN
        AMDE = ABS(MDE)
        WRITE(*,501) 35,U,TOP,ZN,RN,MDE
501     FORMAT(I4,1P,4E15.6,1P,I4)
      ENDIF
C
      IF(AMDE.NE.8) THEN
        TOMEGA = TT*OMEGA
      ELSE
        TOMEGA = OMEGA/TT
      ENDIF
C
C     REAL:
      IV = 0
      CALL QROMB(AIFU,U,TOP,Z1)
C
C   TRANSIENT RESPONSE CHOICE BELOW
      IF(AMDE.EQ.8) GOTO 793
C
C     IMAGINARY:
      IV = 1
      CALL QROMB(AIFU,U,TOP,Z2)
C
      DISTEL = R*DCMPLX(DBLE(Z1),-TOMEGA*DBLE(Z2))*ANM
      GO TO 5000
C
C   DO TWO TYPES OF TRANSIENT RESPONSE
C
793   DISTEL = R*ANM*DCMPLX(DBLE(Z1),0.D0)    !USES NELEM=32 
C
      IF(NCH.EQ.5) THEN
        IF(OMEGA.EQ.PX41.AND.PX45.EQ.0.D0) THEN
          TAO = PX41/(-DLOG(DREAL(DISTEL)))   
          WRITE(*,*) NCH,MDE,TAO
          PX45 = TAO
C          READ(*,*)
          DISTEL = (0.D0,0.D0)
          RETURN
C          STOP
        ELSEIF(PX45.NE.0) THEN
          IF(OMEGA.LT.PX41) THEN
              RARG = DEXP(-OMEGA/PX45)
              DISTEL = R*DCMPLX(DBLE(RARG),0.D0)
          ENDIF
          GOTO 5000
        ENDIF
        GOTO 5000
      ELSEIF(NCH.EQ.4) THEN
        GOTO 5000
      ENDIF
C
C   COUPLING MODEL:DSD, CSDO
797   CONTINUE
        IF(OMEGA.GE.PX41) THEN      !PX41 = P(41) IS TAUC
            TOMEGA = OMEGA/TT   !T/TAUS; OMEGA IS TIME
            RARG = DEXP(-(TOMEGA**PHI))
        ELSE
            TAO = PX41*((TT/PX41)**PHI)
            IF(PX45.NE.0.D0) THEN   !PX45 = P(45) IS GAMX
                BTM = 1.D0 - PHI
                TAO = TAO*DEXP(-BTM*PX45)/PHI
            ENDIF
C
            RARG = DEXP(-OMEGA/TAO)
        ENDIF
      DISTEL = R*DCMPLX(DBLE(RARG),0.D0)
      GOTO 5000
C
C---------------------------------------------------------------------
C
C   WILLIAMS WATTS (FRACTIONAL EXPONENTIAL) RESPONSE USING
C       SERIES AND APPROX FORMULA FOR DRT; CUTOFF POSSIBLE
C       ARBITRARY BETA
C   IF MODE=MDE < 0, THEN CALC CSD1 FORM DIRECTLY 
C   IF MODE=MDE >= 0, THEN CALC GENERAL CSD=CSD0 OR DSD FORM DIRECTLY
C   MODE = +8,-8: TRANSIENT RESPONSE
C   CUTOFF DETERMINED BY VALUE OF U; *** IT WILL USUALLY BE NEGATIVE
C   WHEN -U >= 40 OR SO, THERE IS NO EFFECT IN ANY PRACTICAL MEASURED
C   FREQUENCY RANGE UNLESS BETA IS VERY SMALL.  U =  -LN(TAUS/TAUC) 
C
C
C   COMPUTE NORMALIZATION, INTEGRALS
C
3600  IV = 0
      PHICOM = PHI
      IF(ZN.EQ.0.D0) THEN
        AMDE = ABS(MDE)
        IWR = 1
        ITW = IDINT(DABS(ATEMP) + 1.D-3)
        IF(U.GT.-1.D3) THEN
            BOT = U
        ELSEIF(U.EQ.-1.D3) THEN
            BOT = - 6.D0/PHI
        ELSEIF(U.LT.-1.D3) THEN
            BOT = - 1.2D1/PHI
        ENDIF
        TOP = 1.6D1
        WRITE(*,601) 36,PHI,BOT,TOP,MDE
601     FORMAT(I4,1P,(3E16.6),1P,I4)
        U = BOT
C
      ENDIF
C
      IF(ICAV.EQ.1) THEN
        TOMEGA = 0.D0
        IWT = 1
        DO IJ = 1,5
            ICHG = IJ - 2
                CALL QROMB(GDRT,BOT,TOP,ZN)
            VV(IJ) = ZN
        IF(ZN.EQ.0.D0) THEN
          WRITE(*,61)
61        FORMAT(3X,'U IS TOO LARGE !!!!')
          DISTEL = (0.D0,0.D0)
          RETURN
C          STOP
        ENDIF
        ENDDO
        AIN = 1.D0/VV(2)
        XX1 = AIN*VV(3)
        XXM1 = AIN*VV(1)
        XX2 = AIN*VV(4)
        XX3 = AIN*VV(5)
C
        IF(MDE.LT.0) THEN
          XXM1 = 1.D0/XX1
          XX1 = VV(4)/VV(3)
          XX2 = VV(5)/VV(3)
          XX3 = 0.D0
        ENDIF
C
        RN = PHI
        ICAV = 0
        DISTEL = (0.D0,0.D0)
        RETURN
      ENDIF
C
      IWT = 0
      TOMEGA = 0.D0
      IF(BOT.GT.TOP) THEN
        WRITE(*,71)
71    FORMAT(6X,'****  INTEGRATION LIMITS:  U = BOT > TOP! '/6X,
     +'####  TRY, E.G., U = .99*TOP,  WITH U FIXED'/)
        DISTEL = (0.D0,0.D0)
        RETURN
C      STOP
      ENDIF
C
      CALL QROMB(GDRT,BOT,TOP,ZN)
      ANM = 1.D0/ZN
      IWR = 0
      BETA = PHI
      TTN = TT
      UUN = U
      BOT = U
      IF(AMDE.NE.8) THEN
        TOMEGA = TT*OMEGA
      ELSE
        TOMEGA = OMEGA/TT
      ENDIF
C
C     REAL:
      IV = 0
      CALL QROMB(GDRT,BOT,TOP,Z1)
C
C   TRANSIENT RESPONSE CHOICE BELOW
      IF(AMDE.EQ.8) GOTO 996
C
C     IMAGINARY:
      IV = 1
      CALL QROMB(GDRT,BOT,TOP,Z2)
C
      DISTEL = R*DCMPLX(DBLE(Z1),-TOMEGA*DBLE(Z2))*ANM
C
      GOTO 5000
C
C   DO TRANSIENT RESPONSE
C
996   CONTINUE
C
      Z2 = Z1*ANM
      DISTEL = R*DCMPLX(DBLE(Z2),0.D0)
C
      GOTO 5000
C
C--------------------------------------------------------------------
C   FOR CALCULATING OR FITTING KWW DRT FOR ARBITRARY BETA
C      WILLIAMS WATTS (FRACTIONAL EXPONENTIAL) RESPONSE USING
C      SERIES AND APPROX FORMULA FOR DRT; CUTOFF POSSIBLE
C   USE DATTYP=R.  "FREQUENCY COLUMN OF DATA IS TAU, NOT FREQUENCY;
C   KWW DRT APPEARS IN REAL COLUMN OF DATA
C   INPUT PARAMETERS: R, TAU, BETA, ATEMP (USUALLY <=0)
C
3700  CONTINUE
C
      TOMEGA = OMEGA/TT
      IF(PHI.EQ.0.5D0.AND.ATEMP.LT.-9.D2) THEN
        CALL KWWDP5(TOMEGA,U,SUM)
      ELSE
        CALL KWWDRT(TOMEGA,PHI,SUM,ATEMP)
      ENDIF
C
      IF(MDE.EQ.-1.D0) THEN      !GET K1 FROM K0
          SUM = SUM*OMEGA
      ENDIF
C
      DISTEL = R*SUM
      GOTO 5000
C
C--------------------------------------------------------------------
C
C   END OF DISTRIBUTED ELEMENTS
C   CHECK FOR EPSILON OPTION AND DO TRANSFORM IF APPROPRIATE
C
5000  IF(NELEM.LT.0) DISTEL = 1.D0/(DISTEL*DCMPLX(0.D0,OMEGA*CELCAP))
      RETURN
      END
