C
C       LEVM PROGRAM: MAIN CALCULATIONS -- LEV2A.FOR: LV1.FOR
C        CHANGES FOR STAGED WEIGHTING: FREE/FIXED ITER; FFI
C        CALLS SPFIT
C
C              K : number of functions (2*M for complex) (IN)
C
C           FTOL : termination occurs when both the actual and predicted
C                  relative reductions in the sum of squares are at most
C                  ftol (IN)
C           GTOL : termination occurs when the cosine of the angle
C                  between fvec and any column of the jacobian is at
C                  most gtol in absolute value (IN)
C           XTOL : termination occurs when the relative error between
C                  two consecutive iterates is at most xtol (IN)
C
C              X : array containing an initial estimate of the solution
C                  vector (IN,OUT)
C
C         MAXFEV : termination occurs when the number of calls to fcn
C                  with iflag = 1 has reached maxfev (IN)
C         NPRINT : enables controlled printing of iterates (IN)
C           NFEV : number of calls output (IN,OUT)
C            PEX : parameter list (IN)
C          NFREI : number of free parameters (IN)
C            FV1 : output values array (IN,OUT)
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE MAINCLC(K,FTOL,GTOL,XTOL,X,MAXFEV,NPRINT,NFEV,PEX,
     + NFREI,FV1)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 DATTYP,DFIT,PFIT,FUN,DATTYQ
      INTEGER K,MAXFEV,NPRINT,NFEV,NFREI
      REAL*8 FTOL,GTOL,XTOL,X,PEX,FV1
      INCLUDE  'SIZE.INC'
      DIMENSION X(K),PEX(NTOT),FV1(K),NFREO(40),
     + NFREF(40),NFREW(40)
      COMMON /CM1/ FREQ(NPTS),M,DATTYP
      COMMON /CM2/ Y(NPT2),R(NPT2),FJ(NPT2),P(NTOT),DRSS,ROE,RKE,
     + NS(NPAFR),NFREE(NTOT),NP,ICNT,MN,IRCH,IXI,DATTYQ
      COMMON /CM3/ CELCAP,FUN,DFIT,PFIT
      COMMON /CM11/ MQY,ISPR,ICX,NDF,FQQ
      COMMON /CM16/ IOP,IORIG,NYC,J,IPRINT,LDFJAC,MODE,IFP,IRE,ISTP,JFP,
     + NPH,INE
      COMMON /CM34/ MD,IWT,IXW,INFP,IPL
C
      ISTP = 1
      IOCNT = 0
      RKE = 0.D0
C
      WRITE(*,10)
10    FORMAT(2X,'==========  RUNNING MAIN CALCULATION  =========='/)
C
C     HERE K (OLD KY) IS MD OR 2*MD
C     SET UP FOR POSSIBLE STAGED WEIGHTING WHEN XI FREE
      ISPR = 1
      DO 126 JWX = 1,2
        IF(CELCAP.LT.0.AND.JWX.EQ.2.AND.IRCH.EQ.1) NFREE(32) = 2
        IF(JWX.EQ.1) THEN
          CALL SPFIT(K,FTOL,GTOL,XTOL,X,MAXFEV,NPRINT,NFEV,PEX,
     +    NFREI,FV1)
C
          IF(DATTYP.NE.'C') GOTO 984
C
        ELSEIF(CELCAP.LT.0.AND.NFREE(32).NE.0) THEN
C     SAVE ORIG NFREE AND DEFINE NEW SETS
          ISPR = 0
          DO 127 I = 1, 32
            NFREO(I) = NFREE(I)
            NFREF(I) = NFREE(I)
            IF(I.LT.31) NFREW(I) = 0
127       CONTINUE
          NFREF(31) = 0
          NFREF(32) = 0
          NFREW(31) = NFREO(31)
          NFREW(32) = NFREO(32)
C
          IF(IRCH.EQ.1) THEN
            IRCH = 2
            IWT = 1
C     ABOVE LINES EFFECTIVELY SET IRCH = -2
            IF(CELCAP.EQ.-2.D0.OR.CELCAP.EQ.-4.D0) NFREW(31) = 2
            NFREW(32) = 2
          ENDIF
          IWC = 0
          DIFFXI = 1.D0
116       DO 130 IE = 1,75
            IF(DIFFXI.LT.1.D-5) GOTO 131    
            DO 129 IWX = 1,2
              LLL = 0
              DO 128 I = 1, 32
              IF(IWX.EQ.1) THEN
                  NFREE(I) = NFREW(I)
              ELSE
                  NFREE(I) = NFREF(I)
              ENDIF
              IF (NFREE(I).EQ.0) GOTO 128
              LLL = LLL + 1
              NS(LLL) = I
               X(LLL) = P(I)
128         CONTINUE
            NFREI = LLL
            IF(IWX.EQ.1) THEN
              X32O = X(NFREI)
            ENDIF
C
            NFIX = NP - NFREI
            INS = 0
            IXW = 0
            IXI = 0
            IF(NFREE(32).GT.0.AND.IRCH.GT.1) IXI = 1
            IF(IXI.EQ.1.OR.IWT.EQ.1) IXW = 1
            IF(NFREE(31).EQ.0.AND.NFREE(32).NE.0) INS = 1
            IF(NFREE(31).NE.0.AND.NFREE(32).NE.0) INS = 2
            MQY = NFREI - INS
C
            IF(CELCAP.LE.-3D0.AND.IWX.EQ.2) THEN
              IRCH = 2
               IXI = 0
               IXW = 0
               IWT = 0
C     ABOVE LINES EFFECTIVELY SET IRCH = 2
            ENDIF
C
            CALL SPFIT(K,FTOL,GTOL,XTOL,X,MAXFEV,NPRINT,NFEV,PEX,
     +      NFREI,FV1)
C
            IF(DATTYP.NE.'C') GOTO 984
C
            IF(JWX.EQ.1) GOTO 120
            IF(IWX.EQ.1) THEN
              IWC = IWC + 1
              IF(IWC.EQ.1) WRITE(*,125)
              X32 = X(NFREI)
              DIFFXI = DABS(X32 - X32O)
              WRITE(*,118) IWC,X32O,X32,DIFFXI
C              WRITE(3,118) IWC,X32O,X32,DIFFXI
            ENDIF
            IF(IWX.EQ.2.AND.DIFFXI.LT.1.D-5) THEN
              IF(ISPR.EQ.0) THEN
                  DIFFXI = 1.01D-5
                  IF(IWC.GT.1) ISPR = 1
                  GOTO 116
              ENDIF
            ENDIF
C       
129         CONTINUE
130       CONTINUE
131       CONTINUE
120     ENDIF
126   CONTINUE
984   CONTINUE
125   FORMAT(15X,'******** STAGED FREE/FIXED WEIGHTING ********'/
     + 3X,'ITERATION COUNT      OLD XI        NEW XI       ABS. DIFF')   
118   FORMAT(4X,I5,8X,3(1PD15.6))
      RETURN
      END
