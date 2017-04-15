      SUBROUTINE HSUB(M,FREQ,P,F)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER M
      DOUBLE PRECISION P(*),F(*),FREQ(*)
      COMPLEX*16 QSN,ZBC,HI,FI,ZA,ZB,ZC,ZT,YC,DISTEL,COMEGA,ZA3,
     +  THSQ(2),PSI,AX,BX,RAD,RADS,IOMEGA,AMT(2),GAM(2),CDCOTH,A11,
     +  TT(2),THP(2),ANUM,EG(2),DNM,PTHI,THDI,YMN0,YMN1,A22,ZBCSR,
     +  P11,P12,P21,P22,YMN,E1THS,E2THS,TSM1,TSM2,ANS,YABN,YCDN,
     +  RH(2),GSI,HS,HGS,CXMEGA,ZBCS,ZBCF,EBC,RCC,QX,COOMEGA,EDCD   
      DIMENSION XIA(2),RH0(2),RHIN(2)
      COMMON /CM12/ CLCAP,ATEMP,WF,MAXFEV,ICF,MDE,JCDX
      COMMON /CM47/ ICNT    
      LOGICAL RD30,RD20,C30,C20,RA0,ZA0,RP0,RP1,CA0
C
C   H CIRCUIT:    1/17/11 3/31/11  7/9/11  11/8/11   2/23/12  3/7/12  1/18/13 6/13/13  8/6/13
C
C   *****************    BLOCKING/CONDUCTIVE MODEL
C
C   THIS SUBROUTINE CALCULATES THE IMPEDANCE OF A SET OF NESTED
C   R-C/DISTRIBUTED ELEMENT CIRCUITS.  P11 - P29 ARE INPUTS ASSOCIATED
C   WITH THE BLOCKING/CONDUCTIVE DCE MODEL
C
C     M: number of data points (IN)
C     FREQ: array of frequency values (IN)
C     P: array of model parameters (IN)
C     F: model function values (OUT)
C
C   SET PARAMETER VALUES
C
      RDE2 = P(1)
      TDE2 = P(2)
      UDE2 = P(3)
      PDE2 = P(4)
      NDE2 = IDINT(P(5))
C
      RDE3 = P(6)
      TDE3 = P(7)
      UDE3 = P(8)
      PDE3 = P(9)
      NDE3 = IDINT(P(10))
C
      RA = P(17)
      IRA = 1
      CA = P(18)
      R2 = P(21)
      C2 = P(22)
      R3 = P(23)
      C3 = P(24)
      RP = P(28)
      CP = P(29)
C
C       SET VARIABLES FOR ZERO CHECKING
C
      IF(ICNT.LE.1) THEN
         RP0 = ((RP.EQ.0.D0).AND.(CP.NE.0.D0))
         RP1 = (RP.EQ.0)
        RD30 = ((R3.EQ.0.D0).AND.(NDE3.EQ.0))
        RD20 = ((R2.EQ.0.D0).AND.(NDE2.EQ.0))
         C30 = (C3.EQ.0.D0)
         C20 = (C2.EQ.0.D0)
         RA0 = (RA.EQ.0.D0)
         CA0 = (CA.EQ.0.D0)
         ZA0 = (RD30.AND.C30)
C         ICH0 = (ICH.EQ.0)
      ENDIF
      ICH = IDINT(P(25))
      IF(ICH.EQ.1.OR.ICH.EQ.3) IBC = 1
      IF(ICH.EQ.2.OR.ICH.EQ.4) IBC = 2      
C
      IF(P(11).GE.0) THEN
          INN = 0
          ELL = P(27)
      ELSE
          INN = 1
          ELL= P(12)
      ENDIF
C
      R1 = P(19)
      C1 = P(20)
      ALLA = 8.854187817D-14/CLCAP     !Usual definition of CLCAP
      AREA = ELL/ALLA     
C    
      DNUM = 6.900902285D0
C
      IF(ATEMP.LT.0) THEN
        TEMP = -0.1D0*ATEMP
      ELSE
        TEMP = ATEMP
      ENDIF
      IF(ICH.EQ.0) THEN
        ALINDUC = P(30)
      ELSEIF(ICH.LT.3.OR.ICH.GE.5) THEN
        PHIN = P(30)
        IF(PHIN.LT.0) THEN
            PHI = -PHIN
        ELSE
            PHI = PHIN
        ENDIF
      ENDIF
C
      IF(ICH.EQ.5) THEN
         RH0(1) = P(26)
         RH0(2) = P(27)
        RHIN(1) = P(17)
        RHIN(2) = P(18)
            ELL = P(12)
            PIM = P(14)
             XI = P(16)
             RA = 0.D0
             CA = 0.D0
C           XIA(1) = P(11)
         XIA(2) = P(29)
             CP = 0.D0
            AN0 = -P(11)
          AKRAT = P(15)
        IF(P(13).GE.0.D0) THEN
          PIZ= P(13)
        ELSE
          PIZ= 1.D0
          XIA(1)=-P(13)
        ENDIF
C
      ELSE
        IF(IBC.EQ.1) THEN
            RA = P(17)
            CA = P(18)
            R1 = P(19)
            C1 = P(20)
C           EEPS = P(26)
            ELL = P(27)
        ELSEIF(IBC.EQ.2) THEN
            ELL = P(19)
            EEPS = P(20)
        ENDIF
C   
        IF(ICH.EQ.1.OR.ICH.EQ.2) AN0 = DABS(P(11))
        IF(ICH.EQ.3) THEN
            ENO = P(11)
            ALAM = P(15)
            AN0 = P(30)
        ENDIF
          IF(ICH.EQ.4) THEN
              ENO = P(12)
              AN0 = P(30)
          ENDIF
      ENDIF
C
      IF(ICH.EQ.1.OR.ICH.EQ.3.OR.ICH.GE.5) THEN
C
            XI = P(16)
      ELSEIF(IBC.EQ.2) THEN
            EZN = P(13) 
            EZP = P(14) 
            EKG = P(15)
            EKR = P(16)  
            EMUN = P(17)  
            EMUP = P(18) 
      ENDIF
C
          IF(PIM.LE.1.0D38.OR.PIM.GE.1.0D-38) THEN
                DM = 2.0D0
             ELSE
                DM = 1.0D0
          ENDIF
      IF(ICH.EQ.1.OR.ICH.EQ.2.OR.ICH.GE.5) THEN
         IF(INN.EQ.0) THEN
            BLAM = P(15)
            R1 = P(19)
            TAUD = R1*C1
            AM = P(12)
            TAUGR = XI*TAUD
            EEPS = P(26)
         ELSE
            ELL = P(12)
            AKRAT = P(15)
            BLAM = AKRAT/AN0
C       WRITE(*,*) INN,AKRAT,AN0,BLAM
C       PAUSE 210
            SLAM = 0.5D0*BLAM
            RLAM = SLAM + DSQRT(SLAM*SLAM + BLAM)
            RATnDN = BLAM/RLAM
            ENO = RATnDN*AN0      !n0 CONCENTRATION
            ALAM = RLAM
          IF(ICH.EQ.1) THEN
            EMUN = P(19)
            R1 = 6.24150965D18*ALLA/(EMUN*(1.D0+(1.D0/PIM))*ENO)
          ELSEIF(ICH.EQ.5) THEN
            R1=P(19)
          ENDIF
C           R1 = 6.242D18*ELL/(EMUN*(1.D0+(1.D0/PIM))*ENO)
            TAUD = R1*C1
            EEPS = C1/CLCAP
            ADL = DNUM*DSQRT((TEMP*EEPS)/(DM*ENO))
            AM = ELL/(2.D0*ADL)
            TAUGR = XI*TAUD
        ENDIF
C       
        IF(ICH.EQ.2) THEN
            BLAM = EKG/(EKR*AN0)
        ENDIF       
            SLAM = 0.5D0*BLAM
            RLAM = SLAM + DSQRT(SLAM*SLAM + BLAM)
            RATnDN = BLAM/RLAM
            ENO = RATnDN*AN0      !n0 CONCENTRATION
            ALAM = RLAM
      ENDIF
C
C
      IF(IBC.EQ.2) THEN
        IF(EKG.LT.1.D36) EZP = EZN
        EP0 = EZN*ENO/EZP
        C1 = 8.85419D-14*EEPS/ELL
        EDA = EZP*EMUP*EP0 + EZN*EMUN*ENO
        R1 = 6.242D18*ELL/EDA
        TAUD = R1*C1
        TAUGR = 1.D0/(EKR*ENO)
        EDB = EZP*EZP*EP0 + EZN*EZN*ENO
        DEBL = DNUM*DSQRT(TEMP*EEPS/EDB)     ! CHECK
C
        AM = ELL/(2.D0*DEBL)
        PIZ = EZN/EZP
        PIM = EMUN/EMUP
        EREC = 2.D0*EZN*EZN/EDB         ! ???
        XI = EREC/(TAUD*EKR)
C       ALAM = TAUD*EKG*XI                                    
        IF(ICH.EQ.4) ALAM = EKG/(EKR*ENO)   
      ENDIF
C       
C   
      IF(PIZ.EQ.1.D0) THEN
            EZP = 1.D0
            EZN = 1.D0
      ENDIF
      EKR = 1.D0/(ENO*TAUGR)
      EKG = ALAM/TAUGR
      IF(INN.GT.0) THEN
        EMUN = 6.24150965D18*ALLA/((1.D0+(1.D0/PIM))*R1*ENO) 
        EMUP = EMUN/PIM
        EKT = 1.16045D4/TEMP
        DN = EMUN/EKT
        DP = EMUP/EKT
        IF(ICH.EQ.1) THEN
            R2I = P(17)                         !R2I= RH02  R1I=RHO1 (Dimensionless reaction rates)
        ELSEIF(ICH.EQ.5) THEN
            R2I = P(27)
            R1I = P(26)
        ENDIF
        AK2 = 2.D0*DN*R2I/P(12)
        AK1 = 2.D0*DN*R1I/P(12)

C   WRITE(*,*) INN, ICH, DN, ELL, R2I, AK2
C   PAUSE 212
      ELSE
        EMUN = 1.D30 
        EMUP = EMUN/PIM
      ENDIF
C
C
59    FORMAT(4(1PD8.4)) 
49    FORMAT(4(1PD12.5)) 
51    FORMAT(3X,'EEPS',10X,'EKG',8X,'TAUD',7X,'R1')
52    FORMAT(3X,'BLAM',9X,'EKR',7X,'TAUGR',7X,'C1')
53    FORMAT(3X,'DM',11X,'ADL',8X,'ELL',7X,'AN0') 
54    FORMAT(3X,'ALAM',9X,'PIM',9X,'AK2',8X,'ENO')
55    FORMAT(3X,'AK1',10X,'CLCAP',10X,'XI',9X,'EDS')
56    FORMAT(4X,'EMUN',11X,'EMUP',9X,'DN',10X,'DP')
61    FORMAT(3X,'REZ',10X,'REE',10X,'RCM',9X,'MRAT')
60    FORMAT(3X,'PIM',10X,'ADL1',10X,'RCM1',9X,'TEMP')  
C
C 
      IF(ICH.NE.0.AND.PIZ.NE.0.D0) THEN
        PIZI = 1.D0/PIZ
        DEL1 = 1.D0/(1.D0 + PIZ)
        DEL2 = 1.D0/(1.D0 + PIZI)
      ENDIF
C
       PIMI = 1.D0/PIM
       EPS1 = 1.D0/(1.D0 + PIM)
       EPS2 = 1.D0/(1.D0 + PIMI)
      BLAM1 = DEL1/EPS1
      BLAM2 = DEL2/EPS2
         BB = BLAM1*BLAM2
         DD = BLAM1 + BLAM2
C
C    CONTINUE
C    LOOP OVER ALL FREQUENCIES
C
C   WRITE(*,*) EEPS,TAUD,R1,C1
C   WRITE(*,*) PIM,ALAM,XI,PHI
C   WRITE(*,*) AKRAT,EMUN,AM,ICH
C   PAUSE 'BEFORE FIT'
C
C
      DO 100 I=1,M
        OMEGA = FREQ(I)
        IOMEGA = DCMPLX(0.D0,OMEGA)
C
C   NO BCD ELEMENT IF ICH=0
C        IF(ICH.EQ.0) GO TO 38
C
        COOMEGA = IOMEGA*TAUD
C   PAUSE 396
C
        COMEGA = (IOMEGA*TAUD)**PHI
        PSI = 1.D0 + COMEGA
        HI = 1.D0/(ALAM + COMEGA*XI)
        FI = 2.D0*BB*HI
        AX = 4.D0*COMEGA*PSI*(BB + FI)
        BX = 1.D0 + COMEGA*(DD + FI)
C   IF(I.EQ.1) THEN
C   WRITE(*,*) I,DM,PIM,AM,ICH
C   WRITE(*,*) PSI,HI,FI,BB,DD
C   PAUSE 397
C   ENDIF
C   
C     IF(DM.EQ.1.D0.AND.PIM.GT.1.D49) THEN
      IF(PIM.GE.1.D40) THEN
        QX = 0.5D0*(1.D0 + 2.D0*HI)/(1.D0 + HI)
        THSQ(2) = PSI*QX
        AMT(2) = AM*CDSQRT(THSQ(2))
        GAM(2) = CDCOTH(AMT(2))
        TT(2) = GAM(2) - 1.D0
C       IF(RA.GT.0.D0) THEN
        RHON = RA
        IRA = 0
        YC = ((COMEGA*TT(2))/PSI) + RHON      !RHON = RH02 
        ZBCS = 1.D0/YC
C   WRITE(*,*) OMEGA, TAUD, COMEGA
C   WRITE(*,*) RA, AM, AMT(2), TT(2)
C   WRITE(*,*) ICH,I,PIM,ZBCS
C   PAUSE 599
C
      ELSE
        RADS = BX*BX -AX
        RAD = CDSQRT(RADS)
        THSQ(2) = (AX/(4.D0*BX))*(1.D0 + (AX/((BX + RAD))**2))
        THSQ(1) = BX - THSQ(2)
C   WRITE(*,*) ICH,AX,BX,RAD,THSQ(2),THSQ(1),AM
C   PAUSE 555
        DO 87 IK=1,2
            THP(IK) = THSQ(IK) - PSI
            AMT(IK) = AM*CDSQRT(THSQ(IK))
            GAM(IK) = CDCOTH(AMT(IK))
             TT(IK) = GAM(IK) - 1.D0
             EG(IK) = THP(IK)*GAM(IK)
C   WRITE(*,*) I,IK,AM,AMT(IK),TT(IK),EG(IK)
87      CONTINUE
        DNM= EG(1)*TT(2) - EG(2)*TT(1)
C   WRITE(*,*) ICH, OMEGA, TAUD, COOMEGA, DNM
C   PAUSE 556 
C
C   IF(I.LE.2) THEN
C       WRITE(*,*) BB,DD,HI,FI
C       WRITE(*,*) AX,BX,THSQ(1),THSQ(2)
C       WRITE(*,*) EG(1),EG(2),TT(1),TT(2),BX
C       PAUSE 400
C   ENDIF

C
C   PARTIAL BLOCKING CALCULATIONS BELOW
      IF(ICH.EQ.5) THEN
        DO 413 IV = 1,2
           CXMEGA = COMEGA*XIA(IV)
           RH(IV) = (RH0(IV) + CXMEGA*RHIN(IV))/(1.D0 +CXMEGA)
C       IF(I.EQ.1) THEN
C      WRITE(*,*) I,IV,RH(IV)
C       PAUSE 444
C      ENDIF
C
413     CONTINUE
        GSI = 1.D0/(1.D0 + EPS1*RH(2) + EPS2*RH(1))
        HS = RH(1)*RH(2) + EPS1*RH(1) + EPS2*RH(2)
        HGS = HS*GSI
C
        THDI = 1.D0/(THSQ(1) - THSQ(2))     
        PTHI = THDI/PSI     
        ANS = THDI*(THP(1)*TT(1) - THP(2)*TT(2))
        YMN0 = PTHI*COMEGA*DNM
C       
        A11 = DEL1 + COMEGA*(BLAM1 + EPS2*FI)
        A22 = DEL2 + COMEGA*(BLAM2 + EPS1*FI)
        TSM1 = -THSQ(2) + A11
        TSM2 = -THSQ(2) + A22
        E1THS = -EPS1*THSQ(2)
        E2THS = -EPS2*THSQ(2)
C
        P11 = E1THS + COMEGA*(EPS1*(BLAM2 + FI) + TSM2)
        P21 = E1THS + EPS1 + COMEGA*(DEL1 + TSM1)
        P12 = E2THS + COMEGA*(EPS2*(BLAM1 + FI) + TSM1)
        P22 = E2THS + EPS2 + COMEGA*(DEL2 + TSM2)
        YMN1 = PTHI*((RH(1)*P11 + RH(2)*P12)*TT(1) + (RH(1)*P21
     + +RH(2)*P22)*TT(2))
        YMN = YMN0 + YMN1
        YABN = GSI*YMN + HGS
        YCDN = (YMN + HS)/ANS
C       ZBC = R1*(1.D0/YABN + 1.D0/YCDN)
        ZBCS = (1.D0/YABN + 1.D0/YCDN)

C   WRITE(*,*) ICH,I,GSI,YMN,YMN1,HS,ZBCS
C   PAUSE 624
      ELSE  
C   
C   LINES BELOW FOR COMPLETE BLOCKING
      ANUM = EG(1) - EG(2)
C   WRITE(*,*) ICH,I,ANUM,DNM
C   PAUSE 5561
      QSN=PSI*ANUM/DNM
C           ZBC = R1*QSN/COMEGA
            ZBCS = QSN/COMEGA
C   WRITE(*,*) ICH,I,ZBCS
C   PAUSE 557
      ENDIF
      ENDIF
C
      ZBCSR = ZBCS*R1      !ITERFACE IMPEDANCE
      ZBCF = R1*(1.D0+ZBCS)/(1.D0+COOMEGA*(1.D0+ZBCS))    !FULL IMP
      ZBC = ZBCF
C   WRITE(*,*) R1,ZBCF
C
      IF(ICH.EQ.1.OR.ICH.EQ.2.OR.ICH.EQ.5) THEN   !PHIN = P(30) ALREADY SET
        IF(PHIN.LT.0.D0) THEN
            ZBC = ZBCSR !USE NEGATIVE PHI TO OUTPUT INTERFACE Z ONLY        
         ELSE
            ZBC = ZBCF
        ENDIF
      ENDIF

C   PAUSE 600
C
      IF(I.EQ.1) THEN
C     IF(MAXFEV.EQ.0.AND.IBC.EQ.1) THEN   !1&3 TO 2 MICRO
        IF(MAXFEV.EQ.0) THEN    !1&5 TO OUTPUT
          IF(PIZ.EQ.1.D0) THEN
            EZP = 1.D0
            EZN = 1.D0
          ENDIF
        EKR = 1.D0/(ENO*TAUGR)
C   WRITE(*,*) EKR
C           PAUSE 123
        EKG = ALAM/TAUGR
      IF(INN.GT.0) THEN
        EMUN = 6.241509645D18*ALLA/((1.D0+(1.D0/PIM))*R1*ENO) 
        EMUP = EMUN/PIM
      ELSE
        EMUN = 1.D30 
        EMUP = EMUN/PIM
      ENDIF
C
      EBC =1.D0/(IOMEGA*CLCAP*ZBC)
      REE = DREAL(EBC)
      REZ = DREAL(ZBC)
      RCC = DCMPLX(AM,0.D0)
      RCM = DREAL(CDCOTH(RCC))
C       AMRAT = (EPS(0)/(EEPS)/(MCTNH(M)    !
      AMRAT = (REE/EEPS)/RCM 
C
C       WRITE(*,*) DNUM,TEMP,ENO,AKRAT,EEPS
C       WRITE(*,*) IKX,IBC,MAXFEV,TAUD,R1
C           PAUSE 124
      IF(IKX.GT.0) GOTO 663
      WRITE(*,51)
C   PAUSE 125
      WRITE(*,49) EEPS,EKG,TAUD,R1
C   PAUSE 126
      WRITE(*,52)  
      WRITE(*,49) BLAM,EKR,TAUGR,C1
      WRITE(*,53) 
      WRITE(*,49) DM,ADL,ELL,AN0
C   PAUSE 127
      WRITE(*,54) 
      WRITE(*,49) ALAM,PIM,AK2,ENO
      WRITE(*,56)
      WRITE(*,49) EMUN,EMUP,DN,DP
      WRITE(*,55)
      
C        IF(PIM.GE.1.D38.OR.PIM.LE.1.D-38) THEN
      SQRT2 = DSQRT(2.D0)
      ADL1 = SQRT2*ADL
C           AM1= ELL/(2.D0*ADL1)
      EDS = P(24)/CLCAP
      RCM1 = RCM/SQRT2        
C          ENDIF
C
      WRITE(*,49) AK1,CLCAP,XI,EDS
C       WRITE(*,*) AK1,CLCAP,AKRAT,XI,EDS,P(24)
C       PAUSE 777
C
      WRITE(*,61) 
      WRITE(*,49) REZ,REE,RCM,AMRAT
      WRITE(*,60) 
      WRITE(*,49) PIM,ADL1,RCM1,TEMP
C       ELSE
C       ENDIF                                                     
C       PAUSE 128
          IKX = IKX +1
      ENDIF     
663   CONTINUE
C           IKX = IKX + 1
C       WRITE(*,*) IKX,IBC,MAXFEV
C           PAUSE 129 
C
C       WRITE(*,*) I,OMEGA,CLCAP
C       WRITE(*,*) INN,ZBC,REE
C       PAUSE 333
      ENDIF
C
C     MNUN = DINT(DABS(P(29)))  XXXX
C     MM = MOD(I,MNUM)
C     MM = MOD(I,10)        
C       IF(MM.EQ.0) THEN
C       WRITE(*,*) I,ICH,ZBCF,ZBCSR
C       PAUSE 'AFTER FIT'
C       ENDIF
C
C   END OF BLOCKING/CONDUCTIVE PART OF CALCULATION
C
38    CONTINUE
C
C   WRITE(*,*) P(40)
C   PAUSE 40
C
C   CALCULATE FIRST NESTED CIRCUIT
C
      IF(P(40).LT.4) THEN                                           !PNP P25=5  P40=4
C       GO TO 936
C   ELSE     
      IF (RD30) THEN
        IF (C30) THEN
          ZA = (0.D0,0.D0)
        ELSE
          ZA = 1.D0 / (C3*IOMEGA)
        END IF
      ELSE
        ZA = 1.D0/ (C3*IOMEGA + 1.D0/(R3 + DISTEL(RDE3,TDE3,UDE3,
     +       PDE3,NDE3,OMEGA)))
      ENDIF
C
C   CALCULATE SECOND NESTED CIRCUIT
C
      IF (RD20.AND.ZA0) THEN
        IF (C20) THEN
          ZB = (0.D0,0.D0)
        ELSE
          ZB = 1.D0/(C2*IOMEGA)
        END IF
      ELSE
        ZB = 1.D0 / (C2*IOMEGA + 1.D0/(R2 + ZA + DISTEL(RDE2,TDE2,
     +       UDE2,PDE2,NDE2,OMEGA)))
      END IF
C
C   CALCULATE THIRD NESTED CIRCUIT
C
      RA1=RA*IRA
      IF (RA0) THEN
      IF(CA0) THEN
        ZA = (0,0)
      ELSE
            ZA = 1.D0/(CA*IOMEGA)
      ENDIF
      ELSE
         ZA = RA1/(CA*RA1*IOMEGA + 1.D0)
      ENDIF
      IF(ICH.EQ.0) THEN
        ZC = ZA 
      ELSEIF(ZA.EQ.(0,0)) THEN
        ZC = ZBC
      ELSE
        ZC = ZBC*ZA/(ZBC +ZA)
      ENDIF
C
C     ZC = ZC + ZB + R1 
      ZC = ZC + ZB
C
C    Calculate effects of RP AND CP and/or P33 if present.  NO LONGER USED
C
      IF(RP0) THEN
        YC = CP*IOMEGA
      ELSE
        IF(RP1) THEN
        YC = (0,0)
            YC = (CP*RP*IOMEGA + 1.D0)/RP
        ENDIF
      END IF
      IF(ZC.EQ.(0,0).AND.YC.NE.(0,0)) THEN
        ZT = 1.0/YC
      ELSE
        ZT = (ZC/(1.D0 + ZC*YC))
      ENDIF
      ZT = ZT + ALINDUC*IOMEGA
C
C
      ELSEIF(P(40).GE.4) THEN
C
      IF(RD30) THEN
C   IF(P(40).GE.4) THEN
C       WRITE(*,*) P(40)
C   PAUSE 340

      IF (C30) THEN
          ZA3 = (0.D0,0.D0)
        ELSE
          ZA3 = 1.D0 / (C3*IOMEGA)
        END IF
      ELSE
        ZA3 = 1.D0/ (C3*IOMEGA + 1.D0/(R3 + DISTEL(RDE3,TDE3,UDE3,
     +       PDE3,NDE3,OMEGA)))
      ENDIF
C
      YMN0 =IOMEGA*CLCAP*DISTEL(RDE2,TDE2,UDE2,PDE2,NDE2,OMEGA)  
C                          ! USE DIELECTRIC PARAMETERS in distel:      HND MODEL
      IF(ICH.EQ.1)THEN
      RA1 = RA*IRA
      IF (RA0) THEN
      IF(CA0) THEN
        ZA = (0,0)
      ELSE
        ZA = 1.D0/(CA*IOMEGA)
      ENDIF
      ELSE
        ZA = RA1/(CA*RA1*IOMEGA + 1.D0)
      ENDIF
      IF(ICH.EQ.0) THEN
        ZC = ZA 
      ELSEIF(ZA.EQ.(0,0)) THEN
        ZC = ZBC
      ELSE
        ZC = ZBC*ZA/(ZBC +ZA)
      ENDIF
C
      ELSEIF(ICH.EQ.5) THEN
        ZC = ZBC
      ENDIF
C
      IF(P(19).GE.1.D20) THEN 
        ZC = (0.D0,0.D0)
        EDCD = DISTEL(RDE2,TDE2,UDE2,PDE2,NDE2,OMEGA)   
        YMN0 = IOMEGA*C1 + YMN0
      ENDIF
C   WRITE(*,*) P(19), ICH, NDE2,ZA3, ZBC,ZC, YMN0, EDCD
C   PAUSE 444
      IF(NDE2.EQ.0) THEN
        YMN0 = (0.D0,0.D0)
        ZC = ZC + ZA3
      ELSEIF(ZC.EQ.0.D0.AND.YMN0.NE.(0.D0,0.D0)) THEN
        ZC = (1.D0/YMN0) + ZA3
      ELSEIF(P(19).LT.1.D20) THEN 
        ZC = 1.D0/(YMN0 + (1.D0/ZC)) + ZA3
      ENDIF
C   WRITE(*,*) P(19), NDE2,YMN0,ZBC, ZC
C   PAUSE 447
C   WRITE(*,*) P(40)
C   PAUSE 640
C      
      ENDIF
C
C       WRITE(*,*) P(40)
C    PAUSE 640
C    EBCC =1.D0/(IOMEGA*CLCAP*ZT)
C    ETFUL = 1.D0/(IOMEGA*CLCAP*ZT)
C       WRITE(*,*) P(40),ZA,ZA3,ZBC,YMN0,ZC,ZT,EDCD,EBCC,ETFUL
C    PAUSE 641
C       IF P37=0, P33 IS  R SUB S, BUT SEE BELOW:   Various final series elements   SET FINAL ZT
C               
      IF(P(37).GE.0.D0) THEN
          ZT = ZC + P(33)         !SERIES R               !  P37=0 P33=RsubS
      ELSEIF(P(37).LT.0.D0.AND.P(37).GT.-9.D0) THEN
          ZT = ZC + IOMEGA*P(33)      !SERIES L            ! P37 = -4
      ELSEIF(P(37).LT.-9.D0.AND.P(37).GT.-1.9D1) THEN
          ZT = ZC + 1.D0/(IOMEGA*P(33))   !SERIES C         !P37 = -16
      ENDIF   
C
C         WRITE(*,*) ZC, ZT, P(33),P(37)
C         PAUSE 66
C   
C   RETURN IMPEDANCE VALUES
C   WRITE(*,*) ICH,I,PIM,ZT
C   PAUSE 667
C
      F(I) = DREAL(ZT)
      F(I+M) = DIMAG(ZT)
100   CONTINUE
      RETURN
      END