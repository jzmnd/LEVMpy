      SUBROUTINE GSUB(M,FREQ,P,F)
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'SIZE.INC'
      INTEGER M
      DOUBLE PRECISION P(NTOT),F(NPT2),FREQ(NPT2)
      CHARACTER*1 DATTYP
      EXTERNAL GDAETN
      COMMON /CM79/ OMEGAT,PSI,GAM,DELI
      COMMON /CM78/ DATTYP
      DATA SML /1.D-50/
C
C     THIS SUBROUTINE IS FOR NLS FITTING TO Y = F(X) OR F(OMEGA)
C           WHERE Y IS REAL INPUT AND X (OR T, possibly TEMPERATURE)
C           IS THE OMEGA INPUT VARIABLE 
C           FITTING FUNCTION POSSIBILITIES DESCRIBED IN DOCUMENTATION
C           %%%%%%%%%%%%%% RUN ONLY WITH DATTYP=R OR I CHOICE. FITTING
C                   USES WHICHEVER DATA IS THUS DESIGNATED &&&&&&&&&&&
C
C         M : number of data points (IN)
C      FREQ : array of frequency values (IN)
C         P : array of model parameters (IN)
C         F : model function values (OUT)
C
C   SET PARAMETER VALUES
C
      IF(DATTYP.EQ.'C') THEN
        RETURN
C        STOP
      ENDIF
      AY = P(1)
      BY = P(2)
      CY = P(3)
      DY = P(4) + SML
      EY = P(5)
C
      FY = P(6)
      GY = P(7)
      HY = P(8)
C
      PXY = P(9)
      QY = P(10)
C
      RY = P(11)
      SY = P(12)
      TY = P(13)
C
      WY = P(14)
      YY = P(15)
      ZY = P(16)
C
      AX = P(17)
      BX = P(18) + SML
C
      CX = P(19)
      DX = P(20) + SML
      EX = P(21)
C
      FX = P(22)
      GX = P(23)
      HX = P(24) + SML
C
      XI = P(25)
      GAM = P(26)
      DEL = P(27)
      XLG = P(28)
      XUG = P(29)
      XTA = P(30)
C
      IF(WY.NE.0.D0.AND.XTA.GT.0.D0) THEN
        PSI = XI - GAM
        DELI = 1.D0/DEL
        XTI = 1.D0/XTA
C       CALCULATE GDAE NORMALIZING VALUE
        OMEGAT = 0.D0
                CALL QROMB(GDAETN,-XLG,XUG,ZTO)
        ZTI = 1.D0/ZTO
      ENDIF
C   LOOP OVER ALL FREQUENCIES
      DO 100 I=1, M
        OMEGA = FREQ( I)
        IF(OMEGA.EQ.0.D0) OMEGA = 1.D-60
C
C   START FUNCTION CALCULATION
C
      IF(AY.NE.0) THEN
        IF(EY.NE.0) THEN
          IF(DY.NE.0.D0) THEN
                   YT0 = CY*(DY-OMEGA)/(DY*(OMEGA-EY))   
                   IF(DABS(YT0).GT.1.D2) YT0 = DSIGN(1.D2,YT0) 
                   YT1 = AY*(OMEGA**BY)*DEXP(YT0)
          ENDIF
        ELSE
        IF(XTA.LT.0) THEN
            YT1 = AY*DEXP(-OMEGA*DY)
        ELSE
            YT1 = AY*DEXP(-OMEGA/DY)
        ENDIF
        ENDIF
        ELSE
           YT1 = 0.D0
        ENDIF
      IF(FY.NE.0.D0) THEN
        IF(GY.NE.0.D0) THEN
               YT2 = FY*(GY-OMEGA)/(GY*(OMEGA-HY))
        ENDIF
      ELSE
          YT2 = 0.D0
      ENDIF
      IF(PXY.NE.0.D0) THEN
          YT3 = PXY*(OMEGA**QY) 
      ELSE
          YT3 = 0.D0
      ENDIF
      IF(TY.NE.0.D0 ) YT4 = TY*(RY + SY*OMEGA)
C
C       IF(WY.NE.0.D0.AND.XTA.EQ.0.D0) THEN
      IF(YY.NE.0.D0.AND.XTA.EQ.0.D0) THEN
C             YT5 = WY/(1.D0 + YY*(OMEGA**ZY))
      IF(P(40).EQ.0.D0) THEN
             YT5 = WY/(1.D0 + ((OMEGA/YY)**ZY))
      ELSEIF(P(40).GT.0.D0) THEN
             YT5 = WY/((1.D0 + (OMEGA/YY))**ZY)
      ELSEIF(P(40).EQ.-1.D0) THEN
             YT5 = WY/((1.D0 + (OMEGA/YY)**ZY)**AX)
      ELSEIF(P(40).EQ.-2.D0) THEN
             YT5 = WY/((1.D0 + (OMEGA/BX) + (OMEGA/YY)**ZY)**AX)
      ENDIF
      ELSEIF(WY.NE.0.D0.AND.XTA.GT.0.D0.AND.YY.EQ.0.D0) THEN
C
C   PROGRAM TO CALCULATE GDAE TRANSIENT RESPONSE INVOLVING
C                  GAUSSIAN DAE
C
        OMEGAT = XTI*OMEGA
        CALL QROMB(GDAETN,-XLG,XUG,ZT)
        YT5 =  WY*ZT*ZTI
      ELSE
        YT5 = 0.D0
      ENDIF
C
      IF(AX.NE.0.D0) THEN
        IF(XTA.LT.0) THEN
                    YT6 = AX*DEXP(-OMEGA*BX)                                                                                                                                                                     
        ELSE
                    YT6 = AX*DEXP(-OMEGA/BX)
        ENDIF
      ELSE
           YT6 = 0.D0
      ENDIF
      IF(CX.NE.0.D0) THEN
        IF(XTA.LT.0) THEN
                   YT7 = CX*DEXP(-(OMEGA*DX)**EX)
        ELSE
                       YT7 = CX*DEXP(-(OMEGA/DX)**EX)
        ENDIF
      ELSE 
               YT7 = 0.D0
      ENDIF
      IF(FX.NE.0.D0) THEN
        YT8 = FX*GX*(OMEGA**(GX-1.D0))*(HX**(-GX))*DEXP(-(OMEGA/HX)**GX)
      ELSE
               YT8 = 0.D0
      ENDIF
C
      YT = YT1 + YT2 + YT3 + YT4 + YT5 + YT6 + YT7 + YT8
C
C     POLYNOMIAL FIT AND OTHER FIT FUNCTIONS
C
      IF(P(40).EQ.-4.D0) THEN
          YT = P(41) + OMEGA*(P(42) + OMEGA*(P(43) + OMEGA*P(44))) + YT
      ELSEIF(P(40).EQ.-8.D0) THEN
        OMEGI = 1.D0/OMEGA
      YT = P(41) + OMEGI*(P(42) + OMEGI*(P(43) + OMEGI*P(44))) + YT
C
      ELSEIF(P(40).EQ.-1.6D1) THEN
        YT = P(41)*DSIN(P(42)*OMEGA**P(43)) + YT
      ELSEIF(P(40).EQ.-3.2D1) THEN
        YT = P(41)*DATAN(P(42)*OMEGA**P(43)) + YT
      ELSEIF(P(40).EQ.-6.4D1) THEN
        YT = P(41)*DTANH(P(42)*OMEGA**P(43)) + YT
      ELSEIF(P(40).EQ.-1.28D2) THEN
        YT = P(46)*(OMEGA**P(47))*DEXP(-((P(48)*OMEGA)**P(49)))
C
      ENDIF
C
      IF(DATTYP.EQ.'R') THEN
            F(I) = YT
            F(I+M) = 0.D0
      ELSEIF(DATTYP.EQ.'I') THEN
            F(I) = 0.D0
            F(I+M) = YT
      ENDIF
100   CONTINUE
      RETURN
      END