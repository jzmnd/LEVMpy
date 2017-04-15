C
C     MAIN LEAST SQUARES FIT CALCULATIONS
C       CALLS LMDER, SVDM, MODEL
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE FITCALC(K,FTOL,GTOL,XTOL,X,MAXFEV,NPRINT,NFEV,PEX,
     + NFREI,SDWR,SDWI,IOCNT,SGMAF,FV1,JIT,SGMSQM)
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL FCN1
      INTEGER K,MAXFEV,NPRINT,NFEV,NFREI,IOCNT,JIT
      REAL*8 FTOL,GTOL,XTOL,X,PEX,SDWR,SDWI,SGMAF,FV1,SGMSQM
      CHARACTER*1 DATTYP,DFIT,PFIT,FUN,DATTYQ
      INCLUDE 'SIZE.INC'
      DIMENSION X(*),FVEC(NPT2),FJAC(NPT2,NPAFR),WA(801),
     + IPVT(NPAFR),FV1(*),YTEMP(NPTS),PEX(NTOT)
      COMMON /CM1/ FREQ(NPTS),M,DATTYP
      COMMON /CM2/ Y(NPT2),R(NPT2),FJ(NPT2),P(NTOT),DRSS,ROE,RKE,
     + NS(NPAFR),NFREE(NTOT),N,ICNT,MN,IRCH,IXI,DATTYQ
      COMMON /CM3/ CELCAP,FUN,DFIT,PFIT
      COMMON /CM11/ MQY,ISPR,ICX,NDF,FQQ
      COMMON /CM16/ IOP,IORIG,NYC,J,IPRINT,LDFJAC,MODE,IFP,IRE,
     + ISTP,JFP,NPH,INE
      COMMON /CM18/ SDWC,SDRC,DIAG(NPAFR),IPAR,IOPR
      COMMON /CM34/ MDA,IWT,IXW,INFP,IPL
      DATA SML/1.D-80/,TPII/1.5915494309189534D-1/
C
C     FOR IMAGINARY DATA FIT, MOVE DATA AND R TO FIRST PART OF VECTORS
C
      IF(DATTYP.EQ.'I') THEN
        DO 9998 IW=1,M
          YTEMP(IW) = Y(IW)
          Y(IW) = Y(IW+M)
          R(IW) = R(IW+M)
9998      CONTINUE
          K = M
      ENDIF
      ICNT = 0
C
C     EXECUTE THE NON-LINEAR LEAST SQUARES FITTING ROUTINE
C
      FACTOR = 1.D2
C
C     TO ADD TIKHONOV REGULARIZATION POSSIBILITES:
C     CHANGED FCN TO FCN1, K IN CALL TO KPL: =MPN
      IF(P(33).EQ.0.D0) THEN
        KPL = K
        LDFX = LDFJAC
      ELSE
        KPL = K + NFREI
        LDFX = LDFJAC + NFREI
      ENDIF
C
      IF(ABS(MODE).EQ.4) THEN
        MDX = 2
      ELSE
        MDX = 0
      ENDIF
C
      WRITE(*,10)
10    FORMAT(3X,'** RUNNING LEVENBERG–MARQUARDT ALGORITHM **')
C
      CALL LMDER(FCN1,KPL,NFREI,X,FVEC,FJAC,LDFX,FTOL,XTOL,GTOL,
     + MAXFEV,DIAG,MDX,FACTOR,NPRINT,INFO,NFEV,NJEV,
     + IPVT,WA(N+1),WA(2*N+1),WA(3*N+1),WA(4*N+1),WA(5*N+1))
C
      FNORM = ENORM(K,FVEC)
C
C     FOR IMAGINARY DATA FIT, MOVE DATA AND R BACK TO LAST PART OF VECTORS
C
      IF(DATTYP.EQ.'I') THEN
        DO 9999 IW=1,M
          Y(IW+M) = Y(IW)
          Y(IW) = YTEMP(IW)
          FVEC(IW+M) = FVEC(IW)
          R(IW+M) = R(IW)
9999      CONTINUE
          K = 2*M
      ENDIF        
C
      IF(ISPR.EQ.0) GOTO 947
C
      INFP = INFO
      IF(ISTP.GT.0) THEN
        IF(INFO.LT.0) THEN
          WRITE(*,*) INFO,ICNT,666
          WRITE(*,843)
843       FORMAT(4X,'U & XI OR XI CONVERGENCE PROBLEM. TRY NEW INI
     +TIAL VALUES, (e.g., U = 0.1);'
     +4X,'OR IF FIRST NO. IS -13 OR -14, TRY UNITY WEIGHTING, BUT THEN'
     +/4X,'MAKE SURE U (P31) AND XI (P32) ARE FIXED, NOT FREE'/)
        ENDIF   
      IC=INFO+1
      GO TO(751,752,753,754,755,756,757,758,755),IC
751   WRITE(*,761)
C      WRITE(4,761)
      RETURN
752   WRITE(*,762)
C      WRITE(4,762)
      GO TO 759
753   WRITE(*,763)
C      WRITE(4,763)
      GO TO 759
754   WRITE(*,764)
C      WRITE(4,764)
      GO TO 759
755   WRITE(*,765)
C      WRITE(4,765)
      GO TO 759
756   WRITE(*,766)
C      WRITE(4,766)
      GO TO 759
757   WRITE(*,767)
C      WRITE(4,767)
      GO TO 759
758   WRITE(*,768)
C      WRITE(4,768)
759   CONTINUE
761   FORMAT(/10X,'IMPROPER INPUT PARAMETERS'/)
762   FORMAT(/5X,'ALGORITHM TERMINATES HAVING ESTIMATED THAT THE RELATIV
     +E ERROR'/5X,'IN THE SUM OF SQUARES IS AT MOST FTOL'/)
763   FORMAT(/5X,'ALGORITHM TERMINATES HAVING ESTIMATED THAT THE RELATIV
     +E ERROR'/5X,'BETWEEN THE VECTOR OF FREE PARAMETERS AND THE DESIRED
     + SOLUTION'/5X,'IS AT MOST XTOL'/)
764   FORMAT(/5X,'ALGORITHM TERMINATES HAVING ESTIMATED THAT THE RELATIV
     +E ERROR'/5X,'IN THE SUM OF SQUARES IS LESS THAN FTOL AND THE L2 NO
     +RM OF'/5X,'THE VECTOR OF FREE PARAMETERS AND THE DESIRED SOLUTION'
     +/5X,'IS LESS THAN XTOL.'/)
765   FORMAT(/5X,'ALGORITHM TERMINATES HAVING ESTIMATED THAT FVEC IS ORT
     +HOGONAL TO THE COLUMNS OF THE JACOBIAN'/5X,'TO PRECISION GTOL, THI
     +S IS NOT A GOOD INDICATOR OF A CONVERGED MINIMUM.'/)
766   FORMAT(/5X,'TERMINATION: NUMBER OF CALLS TO FCN WITH IFLAG=1 EXCEE
     +DS MAXFEV.'/)
767   FORMAT(/5X,'TERMINATION: FTOL IS TOO SMALL.NO FURTHER REDUCTION IN
     +THE SUM OF SQUARES IS POSSIBLE.'/)
768   FORMAT(/5X,'TERMINATION: XTOL IS TOO SMALL. NO FURTHER IMPROVEMENT
     + IN THE FREE'/5X,'PARAMETERS IS POSSIBLE.'/)
        IF(IPRINT.GT.0) THEN
C           WRITE(3,750) FNORM, NFEV, NJEV
           WRITE(*,750) FNORM, NFEV, NJEV
C           WRITE(4,750) FNORM,NFEV,NJEV
        ENDIF
750   FORMAT(15X,'FINAL L2 NORM OF THE RESIDUALS  ',D15.7 /15X,'NUM
     +BER OF FUNCTION EVALUATIONS  ',I10 /15X,'NUMBER OF JACOBIAN EVAL
     +UATIONS  ',I10 /)
      ENDIF
C
      DO 760 IW = 1, NFREI
         P(NS(IW)) = X(IW)
760   CONTINUE
      IF((IPRINT.EQ.2.AND.ISTP.EQ.2).OR.(IPRINT.EQ.3.AND.
     +ISTP.GT.0)) THEN
            WRITE(*,315)
C            WRITE(3,315)
            DO 770 IW = 1, N
              IF(NFREE(IW).NE.0) THEN
                  WRITE(*,310) IW,P(IW)
C                  WRITE(3,310) IW,P(IW)
              ENDIF
770         CONTINUE
        ENDIF
310     FORMAT(' ',12X,'P(',I2,') = ',1PD20.13)
315     FORMAT(15X,'FINAL PARAMETER ESTIMATES')
C
                IF(MN.EQ.0) MN = 1
                KJN = K - J + 1
                   AN1 = -1.D0/FLOAT(KJN)
            PUW = P(31)*P(31)
            PXI = 2.D0*P(32)
        IF(IXW.EQ.1) THEN
          IF(IWT.EQ.0) THEN
            YDDS = 1.D0
            DO 487 I = J,K
                FTD = PUW + DABS(Y(I))**PXI
                YDDS = YDDS*(FTD**AN1)
487         CONTINUE
            YPXS = YDDS
        ELSEIF(IWT.EQ.1) THEN
            CALL MODEL(NTOT,P,FV1)
            YDMS = 1.D0
            DO 488 I = J,K
              IF(FV1(I).EQ.0.D0) THEN
                FTD = PUW
              ELSE
                FTD = PUW + DABS(FV1(I))**PXI
              ENDIF
              IF(FTD.EQ.0.D0) THEN
                YDMS = 0.D0
              ELSE
                YDMS = YDMS*(FTD**AN1)
              ENDIF
488         CONTINUE
            YPXS = YDMS + SML
        ENDIF
      ELSE
        YPXS = 1.D0
      ENDIF
      IF(IRCH.GE.5.AND.IWT.EQ.1) YPXS = 1.D0
C
        IF(IXW.EQ.1.AND.JIT.LT.3) THEN
            YPSRI= 1.D0/DSQRT(YPXS)
            DO 489 I = J,K
                R(I) = R(I)*YPSRI
489         CONTINUE
        ENDIF
C
        IF(JIT.EQ.4) YPXS = 1.D0
C
      IF(P(33).EQ.0.D0) THEN
        KJN = KPL
      ELSE
        KJN = KPL - NFREI
      ENDIF
      NDF = KJN - NFREI
C
        ANEFF = 1.D0/FLOAT(NDF)
                SGMUN = FNORM*FNORM
                SPSQ = YPXS*SGMUN
      IF(ROE.LT.0.D0) THEN
        ICX = 1
        ROE = 0.D0
      ELSE
        ICX = 0
      ENDIF
C
      IF(RKE.EQ.0.D0) THEN
        SUMSQ = SPSQ
      ELSE
                SPSR = SPSQ*(KJN/FLOAT(MN))**2
                SUMSQ = (1.D0 + 0.4D0/ROE)*SPSR
      ENDIF
C     CALCULATE AKAI FIT QUALITY FACTOR, FQQ, WHICH BECOMES FQF
      IF(DATTYP.EQ.'I') THEN
        KQ = M
      ELSE
        KQ = K
      ENDIF
        FQQ = KQ*DLOG(SUMSQ+SML) + 2.D0*NFREI
        WRITE(*,8001) KQ,NFREI,FQQ
8001    FORMAT(39X,2I6,1PD12.3)
        SGMASQ = SUMSQ*ANEFF
      IF(JIT.NE.3) THEN
        SGMASV = SGMASQ
      ELSE
        SGMASV = 0.6D0/SGMSQM
      ENDIF
      SGMFUN = ANEFF*SGMUN
      SGMAF = DSQRT(SGMASQ)
C
      IF(ISTP.GT.0) THEN
             IF(RKE.NE.0.D0) WRITE(*,2103) RKE,MN
              WRITE(*,853) SUMSQ, SGMASQ, SGMAF, FQQ, NDF
C             WRITE(3,853) SUMSQ, SGMASQ, SGMAF, FQQ, NDF
C             WRITE(4,853) SUMSQ, SGMASQ, SGMAF, FQQ, NDF
      ENDIF
2103  FORMAT(/' OBJECTIVE FUNCTION, S, INVOLVES LS/AB SWITCH VALUE, RKE=
     +',D10.3, ' AND MN =',I4)
853   FORMAT (10X,'WEIGHTED SUM OF SQUARES, S:',T43,1PD20.13/10X,
     +'SIGMAF SQUARED ESTIMATE, XS:',T43,D20.13/10X,
     +'SIGMAF ESTIMATE:',T43,D20.13/10X,
     +'FIT QUALITY FACTOR:',T46,1PD12.3/10X,
     +'NO. DEGREES OF FREEDOM:',T50,1I5)
856   FORMAT(//10X,'OBSERVED VARIABLES'/5X,'M',4X,'MEASURED',6X,'ESTIMAT
     +ED',3X,'UNCERTANTIES',4X,'RESIDUALS',4X,'RESID./UNCTY')
8561  FORMAT(//10X,'OBSERVED VARIABLES'/5X,'M',4X,'MEASURED',6X,'ESTIMAT
     +ED',3X,'UNCERTANTIES',4X,'RESIDUALS',4X,'RESID./MODEL')
C
      IF(ISTP.GT.0.AND.IPRINT.GT.1) THEN
        WRITE(*,673) DFIT,DFIT,DFIT,DFIT
C        WRITE(3,673) DFIT,DFIT,DFIT,DFIT
673     FORMAT(/4X,'I',6X,'OMEGA',9X,A1,'(IW)',6X,'SIGMA(',A1,
     +'(IW))',4X,A1,'(IW+M)',3X,'SIGMA(',A1,'(IW+M))')
          DO 672 IW=1,M
672       WRITE(*,674)IW,FREQ(IW),Y(IW),R(IW),Y(IW+M),R(IW+M)
674        FORMAT(I5,5(1PD14.5))
C
      END IF
C
C     FINAL OUTPUT OF LEV FOLLOWS
C
      WRITE(*,20)
20    FORMAT(5X,'** RUNNING SVDM ALGORITHM **')
C
      CALL SVDM(K,NFREI,SGMFUN,X,PEX,NS,IORIG,ISTP)
C
      IF(IPRINT.EQ.0) GOTO 947
      CALL MODEL(NTOT,P,FV1)
      YPX = DSQRT(YPXS)
      IF(IPRINT.GE.1.AND.ISTP.GT.0) THEN
        IF(IPAR.LE.0) THEN
            WRITE(*,8561)
C            WRITE(4,8561)
        ELSE
            WRITE(*,8561)
C            WRITE(4,856)
        ENDIF
      ENDIF 
        V2S = 0.D0
        V4S = 0.D0
        V2R = 0.D0
       V2RA = 0.D0
        V4R = 0.D0
       V4RA = 0.D0
        V5R = 0.D0
        V5W = 0.D0
        V6R = 0.D0
        V6W = 0.D0
C      OPEN(7,FILE='AUXPNTL1')
      DO 859 IW=1,M
      IF(IOPR.EQ.1) THEN
        FREQX = TPII*FREQ(IW) 
      ELSE
        FREQX = FREQ(IW) 
      ENDIF
      IF(IPL.EQ.1) THEN
        FREQXX = DLOG10(FREQX + 1.D-20)
      ELSE
        FREQXX = FREQX
      ENDIF
          IF(DATTYP.EQ.'I') GOTO 10000
           V1 = -FV1(IW) + Y(IW)
           V2 = V1/Y(IW)
          V2R = V2R + V2
         V2RA = V2RA + DABS(V2)
          V2S = V2S + V2*V2
           V5 = V1/R(IW)
      IF(IPAR.GT.0) THEN
        V7 = V5
      ELSE
        V7 = V1/FV1(IW)
      ENDIF
          V5R = V5R + V5
          V5W = V5W + V5*V5
          IF(IPRINT.GE.2.AND.ISTP.GT.0) WRITE(*,857) IW,Y(IW),FV1(IW),
     +R(IW),V1,V2
          IF(IPRINT.GE.1.AND.ISTP.GT.0) THEN
        WRITE(*,857) IW,Y(IW),FV1(IW),R(IW),V1,V7
        IF(DATTYP.EQ.'R') THEN
            V7L = DLOG10(DABS(V7)+SML) + MODE
            WRITE(*,865)IW,FREQXX,V1,V5,V7,V7L
        ENDIF
      ENDIF
10000   CONTINUE
857      FORMAT(1X,I5,5(1PE14.5))
         IF(DATTYP.EQ.'R') GOTO 860
           V3 = -FV1(IW+M) + Y(IW+M)
           V4 = V3/Y(IW+M)
          V4R = V4R + V4
         V4RA = V4RA + DABS(V4)
          V4S = V4S + V4*V4
           V6 = V3/R(IW+M)
      IF(IPAR.GT.0) THEN
        V8 = V6
      ELSE
        V8 = V3/FV1(IW+M)
      ENDIF
          V6R = V6R + V6
          V6W = V6W + V6*V6
        IF(IPRINT.GE.2.AND.ISTP.GT.0) WRITE(*,857)IW,Y(IW+M),FV1(IW+M),
     +R(IW+M),V3,V4
        IF(IPRINT.GE.1.AND.ISTP.GT.0) THEN
          WRITE(*,857)IW,Y(IW+M),FV1(IW+M),R(IW+M),V3,V8
          WRITE(*,865)IW,FREQXX,V1,V3,V7,V8
      ENDIF
860   CONTINUE
859   CONTINUE
865      FORMAT(1X,I5,1PE14.5,4(1PE13.4))
c      CLOSE(7)
C
      IF(DATTYP.EQ.'I') THEN
        KX = M
        KY = 2*M
      ELSE
        KX = K
        KY = K
      ENDIF
      SUMSQW = V5W + V6W + 1.D-30
         PDK = 1.D0/ALOG(FLOAT(KX))
        SUMP = 0.D0
       SUMFR = 0.D0
        DO 104 IKK = J,KY
                SUMFR = SUMFR + FVEC(IKK)
        VRW = (FV1(IKK) - Y(IKK))/R(IKK)
        QQK = VRW*VRW/SUMSQW
        IF(QQK.EQ.0.D0) THEN
            PPK = 0.D0
        ELSE
            PPK = QQK*DLOG(QQK+SML)
        ENDIF
           SUMP = SUMP + PPK
104        CONTINUE
        SSE = 1.D0 + PDK*SUMP
      IF(NFREI.GE.M) THEN
        ANPEF = 1.D0
      ELSE
        ANPEF = 1.D0/FLOAT(M - NFREI)
      ENDIF
         ANK = 1.D0/FLOAT(KX)
       SUMRR = -(V2R - V4R)
        AVRR = SUMRR*ANK
        AVFR = SUMFR*ANK*YPX
        SDRR = DSQRT(V2S*ANPEF)
        SDRI = DSQRT(V4S*ANPEF)
        SDRC = DSQRT((V2S + V4S)*ANEFF)
        SDWR = DSQRT(V5W*ANPEF)
        SDWI = DSQRT(V6W*ANPEF)
        SDWC = DSQRT(SUMSQW*ANEFF)
       AVARR = V2RA/M
       AVARI = V4RA/M
C
      IF(ISTP.GT.0) THEN 
          WRITE(*,107)
C          WRITE(3,107)
C          WRITE(3,109) AVARR,SDWR,SDWI,SDWC,AVARI,NYC
          WRITE(*,109) AVARR,SDWR,SDWI,SDWC,AVARI,NYC
          WRITE(*,108)
C          WRITE(3,108)
C          WRITE(3,105) AVRR,SDRR,SDRI,SDRC,SSE,IOCNT
          WRITE(*,105) AVRR,SDRR,SDRI,SDRC,SSE,IOCNT
          IF(NFREI.GE.M) WRITE(*,812)
      ENDIF
107   FORMAT(/7X,'AV ARR',6X,'SD WR',7X,'SD WI',7X,'SD WC',7X,'AV ARI',
     +5X,'# IT')
108   FORMAT(/7X,'AV RR',6X,'SD RR',7X,'SD RI',7X,'SD RC',8X,'SSE',
     +4X,'# OPIT')
105   FORMAT(2X,5(1PD12.3),I5/)
109   FORMAT(2X,5(1PD12.3),I5)
812   FORMAT(/1X,'**!!! FORMULA for RE and IM SD CALCULATION CHANGED **
     +*!!!  NFREE > M !!',/)
C
947   CONTINUE
      RETURN
      END
