C
C   SPECIAL FOR DISSADO-HILL FITTING :   USES SERIES
C
C     USES SUBROUTINES: SET, HYPER, EPSALG, AND GAMMA 
C
C     R1: scale
C     TA: tau
C     AM: m
C     AGN: 1 - n
C     WLL is middle barrier
C
C     FIXED VALUES: EPS = 1.D-8; NMAX = 11; WLL = 0.69
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE DISSHILL(R1,TA,AM,AGN,DHRL,DHIM,OMEGA)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 R1,TA,AM,AGN,DHRL,DHIM,OMEGA
      DATA PI2/1.570796326795D0/
C
      IF(AM.GT.9.D0) AM = AGN
C
C   ASYMMETRIC CASE: LIKE EDAE1
      IF(AGN.GT.9.D0.AND.AM.LT.0.D0) AGN = 1.D-6 - AM
C
      AN = 1.D0 - AGN
C
      IF(IQZ.EQ.0) THEN
        OMEGO = OMEGA
        IQZ = 1
      ENDIF
C
      DOM = DABS(OMEGA - OMEGO)
      IF(DOM.LT.1.D-12) DOM = 0.D0
C
C   NORMALIZE AMPLITUDE; LATER MULT BY R1
C
      IF(DOM.EQ.0) THEN
C
        TCS = DCOS(PI2*AM)
        TSS = DSIN(PI2*AM)
        AN1 = AN - 1.D0
        TNC = DCOS(PI2*AN1)
        TNS = DSIN(PI2*AN1)
        ANM = AN - AM
C
C   CALCULATE NORMALIZATION CONSTANT, F0
C
        S1 = 1.D0
        S2 = AN
        S3 = 1.D0 + AM
        IF(AM.GE.0.D0) GOTO 500
        VV = -ANM
        IF(VV.LT.0.D0) GOTO 390
        ZZ = 1.D0 - AN
        CALL GAMMA(ZZ,G1,HH)
        H1 = HH
        ZZ = 1.D0 + AM
        CALL GAMMA(ZZ,G1,HH)
        H2 = HH
        ZZ = 1.D0 - ANM
        CALL GAMMA(ZZ,G1,HH)
        H3 = HH
C        IF(H3*AM.EQ.0.D0) PAUSE 1
        FO = H1*H2/(H3*AM)
        GOTO 650
390   CONTINUE
        ZZ = 1.D0 - AN
        CALL GAMMA(ZZ,G1,HH)
        H1 = G1
        ZZ = 1.D0 + AM
        CALL GAMMA(ZZ,G1,HH)
        H2 = HH
        ZZ = 2.D0 - ANM
        CALL GAMMA(ZZ,G1,HH)
        H3 = HH
C        IF(H3*AM.EQ.0.D0) PAUSE 2
        FO = H1*H2*(1.D0 - ANM)/(H3*AM)
        GOTO 650
500   CONTINUE
        ZZ = AM
        CALL GAMMA(ZZ,G1,HH)
        H1 = HH
        ZZ = 1 - AN
        CALL GAMMA(ZZ,G1,HH)
        H2 = G1
        ZZ = -ANM
        IF(ZZ.GE.0.D0.AND.ZZ.LT.1.D-11) ZZ = 1.D-11
        IF(ZZ.LT.-1.D-11) GOTO 610
        CALL GAMMA(ZZ,G1,HH)
        H3 = G1
        GOTO 640
610     ZZ = 1.D0 - ANM
        CALL GAMMA(ZZ,G1,HH)
        H3 = HH
C        IF(H3.EQ.0.D0) PAUSE 3
C
640     FO = (H1*H2)/H3
650   CONTINUE
C        IF(FO.EQ.0.D0) PAUSE 4
        FOI = DABS(R1/FO)
      ENDIF
C
      XW = OMEGA*TA
      IF(XW.GE.1.000001D0) GOTO 147
C
      CALL SET(S1,S2,S3,S4,S5,D6,DD,LL,AA,BB,CC,IC)
      CALL HYPER(-XW,LL,S4,S5,AA,BB,CC,DD,D6,IC)
      XP = XW**AM
      T1 = XP*TCS
      T2 = XP*TSS
      T5 = DATAN(XW)
      XP2 = (1.D0 + XW**2)**(ANM/2)
      U1 = XP2*DCOS(T5*ANM)
      U2 = XP2*DSIN(T5*ANM)
      GR = FO + AN1*(U1*(S4*T1 - S5*T2) - U2*(S5*T1 + S4*T2))/AM
      GI = -(U1*(S5*T1 + S4*T2) + U2*(S4*T1 -S5*T2))*(1.D0 - AN)/AM
      GR = GR*FOI
      GI = GI*FOI
C
      GOTO 333
147   XI = 1.D0/XW
      S1 = 1.D0 - AN
      S2 = 1.D0 - ANM
      S3 = 2.D0 - AN
C
      CALL SET(S1,S2,S3,S4,S5,D6,DD,LL,AA,BB,CC,IC)
      CALL HYPER(XI,LL,S4,S5,AA,BB,CC,DD,D6,IC)
      ANP = XI**(-AN1)
      T1 = ANP*TNC
      T2 = ANP*TNS
      GR = (S4*T1 - S5*T2)*FOI
      GI = (S4*T2 + S5*T1)*FOI
333   CONTINUE
C
      DHRL = GR
      DHIM = GI
C
      RETURN
      END
