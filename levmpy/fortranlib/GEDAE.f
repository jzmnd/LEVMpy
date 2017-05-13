C
C   CALCULATE GENERALIZED EDAE, HN, OR CD RESPONSE
C      Generalized exponential distribution of relaxation times model (GED)
C      Cole-Davidson (CD) model
C      Havriliak-Negami (HN) model
C
C       II : ith element of array (IN)
C        M : number of data points (IN)
C     FREQ : frequency array (IN)
C       QX : model parameters, subset of 8 (IN)
C        F : output array (IN,OUT)
C      JCD : 0 or 1: conductive-system vs dielectric-system dispersion (IN)
C      NCH : 0 to 6: selects model, see O-circuit (IN)
C      INH : 0 or 1: 1 if NCH1 and NCH2 are > 0  (IN)
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE GEDAE(II,M,FREQ,QX,F,JCD,NCH,INH)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL GDAEFN1,CDFG,HNFG
      INTEGER II,M,JCD,NCH,INH
      REAL*8 FREQ,QX,F
      COMPLEX*16 ZI0C,TCOMEGA
      INCLUDE 'SIZE.INC'
      DIMENSION QX(8),F(*),FREQ(*),VV(5),PHI(2),GAM(2),UU(2)
C      COMMON TO PASS DATA TO DAEFN FUNCTION
      COMMON /CM9/ TOMEGA,PHIX,ICHG,IWT
      COMMON /CM13/ RX,TX,UX,PHIZ,XXM1,XX1,XX2,XX3,RNX,AIN,ICAV,
     + NELEM,NCHA
      COMMON /CM29/ GAMX,SN1,CS1,PII,CDN,QPI,MDEX
      COMMON /CM12/ CELCAP,ATEMP,WF,MAXFEV,ICF,MDE,JCDX
      COMMON /CM55/ PX1,PX41,PX45
      DATA PI/3.1415926535898D0/,PII/0.3183098861838D0/
C
      SAVE RDAE,TDAE,UU,U1A,PHI,GAM,US,RNE
C
      ISW = 1
      IF(II.EQ.1.OR.INH.EQ.1) THEN
        MDEXA = IABS(MDE)
        IF(JCD.EQ.0) THEN           ! CSD ONLY
            MDEX = MDE
        ELSE
            MDEX = MDEXA            ! DSD
        ENDIF
C
        RDAE = QX(1)
C
        IF(MDEXA.EQ.2.AND.JCD.EQ.0) THEN    
            RTOT = PX1 + RDAE
            TDAE = QX(2)*RTOT*RTOT*CELCAP/RDAE   !EPSTAU CSD
        ELSE
            TDAE = QX(2)            ! TAU
        ENDIF
C
        UU(1) = QX(3)
        UU(2) = QX(4)
          U1A = -DABS(UU(1))
        PHI(1) = QX(5)         ! HERE PHI1 IS OLD PHI
C
        PHI(2) = QX(6)
        GAM(1) = QX(7)
        GAM(2) = QX(8)
C
        IF(MDEXA.EQ.8) THEN
          IWT = 0
        ELSE
          IWT = 1
        ENDIF
C
        IF(NCH.EQ.2) THEN
          SN1 = DSIN(PHI(1)*PI)
          QPI = 1.D0/(1.D0 - PHI(1))
          CDN = PII*SN1*QPI
        ELSEIF(NCH.EQ.3) THEN
          SN1 = DSIN(GAM(1)*PI)
          CS1 = DCOS(GAM(1)*PI)
        ENDIF
C
C      CALCULATE NORMALIZATION COEFFICIENT, RN
        IF(NCH.EQ.1) THEN
          IF(UU(1).GT.0) THEN
              US = UU(1)
              SYM = 0.5D0
          ELSE
              US = 0
              SYM = 1.D0
          ENDIF
          IF(UU(2).EQ.0) THEN
            MC = 1
          ELSE
            MC = 2
          ENDIF
C
          RNI = 0.D0
          DO IJ = 1,MC
            KK = 3 - 2*IJ
            PP = 1.D0/GAM(IJ)
            GAMP = GAMMLN(PP)
            XGI = DEXP(GAMP)
            XGP = GAM(IJ)**((1.D0 - GAM(IJ))/GAM(IJ))
            IF(PHI(IJ).NE.0.D0) THEN
                RNU = (1.D0 - DEXP(-PHI(IJ)*DABS(UU(IJ))))/PHI(IJ)
                RNI = RNI + XGI*XGP*RNU
            ELSE
                RNI = RNI + DABS(UU(IJ))*XGI*XGP
            ENDIF
          ENDDO
          RNE = SYM/RNI
        ELSE
          RNE = 1.D0
        ENDIF
C
        ICHG = 0
        TOMEGA = 0.D0
        ISW = 0
      ENDIF
      IF(ISW.EQ.0) GOTO 199
189   CONTINUE
      IF(ISW.EQ.0) RN = 1.D0/ZRE
      ISW = 1
C
C   CALCULATE FINAL LIMITS FOR CSD QUANTITIES
C
      IF(ICAV.EQ.1) THEN
        TOMEGA = 0.D0
        DO IJ = 1,5
          IM = IJ - 2
          ICHG = IM
          PHIX = PHI(1)
          GAMX = GAM(1)
C
          VINT2 = 0.D0
C
          IF(NCH.EQ.1) THEN
            CALL QROMB(GDAEFN1,U1A,US,VINT1)
          ELSEIF(NCH.EQ.2) THEN
            UTOP = (DEXP(-U1A) - 1.D0)**(1.D0 - PHI(1))
            CALL QROMO(CDFG,0.D0,UTOP,VINT1)
          ELSEIF(NCH.EQ.3) THEN
            CALL QROMB(HNFG,U1A,-U1A,VINT1)
          ENDIF
C
          IF(UU(2).GT.0) THEN
C     REAL:
            PHIX = PHI(2)
            GAMX = GAM(2)       
            CALL QROMB(GDAEFN1,0.D0,UU(2),VINT2)
C       VINT2 = VINT2*RNE(2)
          ENDIF
          ZRE = (VINT1 + VINT2)*RNE
          VV(IJ) = ZRE
        ENDDO
C
        AIN = 1.D0/VV(2)
        XX1 = VV(3)/VV(2)
        XC0 = 1.D0/XX1
        XXM1 = AIN*VV(1)
        XX2 = AIN*VV(4)
        XX3 = AIN*VV(5)
C
        RNX = RN
        IF(MDE.LT.0) THEN
          IF(NCH.EQ.2.OR.NCH.EQ.3) THEN
              XXM1 = 1.D0/XX1
              XX1 = VV(4)/VV(3)
              XX2 = VV(5)/VV(3)
              XX3 = 0
          ELSEIF(NCH.EQ.1.AND.JCD.EQ.0) THEN
              AIN = VV(1)
              RNX = AIN
          ENDIF
        ENDIF
        IF(M.EQ.-1) RETURN
C
      ENDIF
      ICAV = 0
C
C   DO OMEGA VARIABLE OR T VARIABLE: |MDEX| = |MDE|= 8 FOR TRANSIENT
C
      IF(IWT.EQ.0) THEN 
        TOMEGA = FREQ(II)/QX(2)
      ELSE
        TOMEGA = TDAE*FREQ(II)  
        TCOMEGA = DCMPLX(0,FREQ(II))
      ENDIF
C
C   SET UP VARIABLES TO PASS DATA TO GDAE FUNCTION THRU /CM9/
C   COMPUTE INTEGRALS
C
199   CONTINUE
C     REAL:
      PHIX = PHI(1)
      GAMX = GAM(1)
      IF(ICAV.EQ.0) ICHG = 0
C
      VINT2 = 0.D0
C
      IF(NCH.EQ.1) THEN
            CALL QROMB(GDAEFN1,U1A,US,VINT1)
      ELSEIF(NCH.EQ.2) THEN
        UTOP = (DEXP(-U1A) - 1.D0)**(1.D0 - PHI(1))
            CALL QROMO(CDFG,0.D0,UTOP,VINT1)
      ELSEIF(NCH.EQ.3) THEN
            CALL QROMB(HNFG,U1A,-U1A,VINT1)
      ENDIF
C
      IF(UU(2).GT.0) THEN
C     REAL:
        PHIX = PHI(2)
        GAMX = GAM(2)       
            CALL QROMB(GDAEFN1,0.D0,UU(2),VINT2)
C       VINT2 = VINT2*RNE(2)
      ENDIF
      ZRE = (VINT1 + VINT2)*RNE   
C
      IF(ISW.EQ.0) GOTO 189
      IF(IWT.EQ.0.D0) GOTO 192
C
C     IMAGINARY:
      PHIX = PHI(1)
      GAMX = GAM(1)
          ICHG = 1
      VINT2 = 0.D0
      IF(NCH.EQ.1) THEN
              CALL QROMB(GDAEFN1,U1A,US,VINT1)
      ELSEIF(NCH.EQ.2) THEN
              CALL QROMO(CDFG,0.D0,UTOP,VINT1)
      ELSEIF(NCH.EQ.3) THEN
              CALL QROMB(HNFG,U1A,-U1A,VINT1)
      ENDIF
      IF(UU(2).GT.0) THEN
C     IMAGINARY:
        PHIX = PHI(2)
        GAMX = GAM(2)       
              CALL QROMB(GDAEFN1,0.D0,UU(2),VINT2)
      ENDIF
      ZIM = (VINT1 + VINT2)*RNE
C
      IF(NCH.GT.1.AND.UU(1).LT.0) WRITE(*,*) II
192   CONTINUE
C
      ZI0C = RN*DCMPLX(ZRE, -IWT*TOMEGA*ZIM)
      IF(NCH.EQ.2) THEN
        IF(MDE.EQ.-1) ZI0C = (1.D0 - ZI0C)/(TDAE*TCOMEGA*PHI(1)) 
      ELSEIF(NCH.EQ.3) THEN                    ! R IS EPSCINF HERE
        IF(MDE.EQ.-1) THEN
            ZI0C = (1.D0 - ZI0C)/(TCOMEGA*QX(1))   !QX TCOME
            RDAE = 1.D0/CELCAP
        ENDIF
      ENDIF
C
C   IF MDE = -6, CD AND HN GIVE CSD1a.
C
      F(II) = RDAE*DREAL(ZI0C)
      F(II+M) = RDAE*DIMAG(ZI0C)
      RETURN
      END
