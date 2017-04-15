C
C   JRM 12/21/93
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE GBNM(XC,FW,RINN,NPH,ZM)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER NPH
      REAL*8 XC,FW,RINN
      COMPLEX*16 ZM,CFWI,YM,AMM
      DATA P0/1.D0/,P1/0.99938487D+00/,P2/-0.18911536D+00/,
     +P3/-0.17480496D+01/,P4/0.30317138D+00/,P6/0.43155789D+01/,P7/
     +0.18521559D+02/,P8/0.70162962D+01/,P9/0.10471997D+02/,P11L/-0.2616
     +6399D-01/,P12/0.27391968D-00/,P13/0.76160040D+00/,P14/0.81790503D-
     +01/,P17/0.56807327D+00/,P18/0.34347179D+00/,
     +P19/-0.36166320D+05/,P20/0.16284434D+07/,P27/0.62794829D-01/,
     +P24/3.6083435D0/,P5/0.35977899D+00/,P10/0.39121674D+00/,
     +P15/0.65292569D0/,P21/-0.23328332D+00/,P22/0.40533549D+01/,
     +P23/0.50536137D0/   
C
      IF(XC.LT.5.D1) THEN
        TADD = DEXP(-XC)
      ELSE
        TADD = 0.D0
      ENDIF
      IF(XC.LT.3.D0) THEN
        WRITE(*,*) 'XC TOO SMALL'
        RETURN
C        STOP
      ENDIF
C
      TC1 = P0 - (1.5D0/XC) + TADD 
C   GOOD APPROX DOWN TO XC=5
      FW = FW*TC1
      CFWI = 1.D0/DCMPLX(0.D0,FW)
      CALL LOGAC(FW,FX)
      FXS = FX*FX
      IF(FW.GT.P0) THEN
        CALL LOGAC(FW**P12,ALP)
        ALB = P11L + ALP
        CALL LOGAC(ALB,YL)
        ACR = (FXS + P6*FX +P7)/(FXS + P8*FX + P9)
        AMR = P1*FX*ACR*(P0 + (P2*(YL**P3))) 
C
        CA = (P17 + P4*FX + P27*FXS)/(P0 + P18*FX + P27*FXS)
        CB = (P0 + P20*FW*FW)/(P0 + FW*(P19 + P20*FW))
        AMI = CA*CB*DATAN(FX/P17)
      ELSE
C       LOW OMEGA PART
        FS = FW*FW
        CALL LOGAC(FS,FNU)
        ALI = P24*FS
        CALL LOGAC(ALI,FNC)
        AMR = FNU/(P0 + FNC/(P0 + P13*(FNC**P14)))
C
        CA = (P10 + P5*FX + P23*FXS)/(P0 + P15*FX + P23*FXS)
        CB = (P0 + P22*FW*FW)/(P0 + FW*P21 + P22*FS)
        AMI = CA*CB*DATAN(FX/P10)
      ENDIF
C
      AMM = DCMPLX(AMR,AMI)
C       SUBTRACT G0N = 1 (AT YN LEVEL) HERE 
      AMM = AMM/(P0 - AMM*CFWI)
C       ADD BACK G0N IF PHI.NE.8.D0
      IF(NPH.EQ.8) THEN
        ADD = 0.D0
      ELSE
        ADD = P0
      ENDIF
      YM = P0/(CFWI*AMM) + ADD
C      EPS = YM*CFWI
      ZM = (1.D0 - RINN)*(P0/YM) + RINN
C
      RETURN
      END
