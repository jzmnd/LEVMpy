C
C     APPLIES RESIDUAL WEIGHTING
C
C           KY : number of functions (2*M for complex) (IN)
C       DATTYP : data type C, R or I (IN)
C           YQ : data input (IN)
C           FT : weights on function (IN,OUT)
C
C       MODIFIED JEREMY SMITH 3/31/2017
C
      SUBROUTINE RWTS(KY,DATTYP,YQ,FT)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*1 DATTYP,DATTYQ
      INTEGER KY
      REAL*8 YQ,FT
      INCLUDE 'SIZE.INC'
      COMMON /CM2/ Y(NPT2),R(NPT2),FJT(NPT2),P(NTOT),DRSS,ROE,RKE,
     + NS(NPAFR),NFREE(NTOT),NP,ICNT,MN,IRCH,IXI,DATTYQ
      COMMON /CM34/ MD,IWT,IXW,INFP,IPL
      COMMON /CM35/ JIT,IPF,NPRIN
      COMMON /CM36/ SHW(NPT2),ISW
      DIMENSION YQ(NPT2),FT(NPT2),RT(NPT2)
      DATA AK1S,AK2S/1.D-6,1.D-8/,OAA5,OAB5,OAG5/8.12D-4,9.33D-4,
     + 2.31D-4/,OAA6,OAB6,OAG6/1.0D-25,6.9966D-3,5.3903D-6/
C
      IF(IRCH.EQ.5) THEN
        OAA = OAA5
        OAB = OAB5
        OAG = OAG5
      ELSEIF(IRCH.EQ.6) THEN
        OAA = OAA6
        OAB = OAB6
        OAG = OAG6
      ENDIF
        XI = P(32)
        XT2 = 2.D0*XI
      GOTO (100,200,300,400,500,600) IRCH
C
100   CONTINUE
C       UNITY WEIGHTING
      IF(JIT.GE.2) GOTO 200
      DO 650 IJ = 1,KY
        FT(IJ) = 1.D0
650   CONTINUE
      GOTO 1677
C
200   CONTINUE
C       P0WER LAW (PROPORTIONAL) WEIGHTING
        DO 670 I = 1,KY
            RT(I) = DABS(YQ(I))
            FT(I) = RT(I)**XT2
670     CONTINUE
        GOTO 675
C
300   CONTINUE
C       MODULUS PL WEIGHTING
      IF(DATTYP.NE.'C') GOTO 871
      DO 678 I=1,MD
        FT(I) = YQ(I)*YQ(I) + YQ(I+MD)*YQ(I+MD)
        RT(I) = DSQRT(FT(I))
        RT(I+MD) = RT(I)
        FT(I) = FT(I)**XI
        FT(I+MD) = FT(I)
678   CONTINUE
      GOTO 675
C
400   CONTINUE
C     SOLATRON 1174 FRA PL WEIGHTING FOLLOWS (SEE SPINOLO REF.)
      IF(DATTYP.NE.'C') GOTO 871
      DO 778 I=1,MD
         YQS1 = YQ(I)*YQ(I)
         YQS2 = YQ(I+MD)*YQ(I+MD)
        FT(I) = AK1S*(2.D0*YQS1 + YQS2) + AK2S
        RT(I) = DSQRT(FT(I))
        FT(I+MD) = AK1S*(2.D0*YQS2 + YQS1) + AK2S
        RT(I+MD) = DSQRT(FT(I+MD))
        FT(I) = FT(I)**XI
        FT(I+MD) = FT(I+MD)**XI
778   CONTINUE
      GOTO 675
871   WRITE(*,*) 'DATTYP MUST BE C FOR THIS IRCH VALUE'
      RETURN
C      STOP
C
500   CONTINUE
C     SOLATRON #1250,1286 FRA WEIGHTING FOLLOWS (SEE ORAZEM REF.)
C
      IF(DATTYP.NE.'C'.AND.IWT.EQ.0) GOTO 872
      I1 = IDINT(DABS(P(34)) + 1.D-3)
      IF(I1.EQ.0) I1 = MD 
      I2 = IDINT(P(36) + 1.D-3)
      IF(I2.EQ.0) I2 = MD 
      I3 = IDINT(P(38) + 1.D-3)
      IF(I3.EQ.0) I3 = MD
C
      DO 878 I=1,MD
      IF(I.LE.I1) THEN
        IF(P(33).LE.0.D0) THEN
            WRITE(*,*) 'P(33) MUST BE POSITIVE FOR THIS IRCH VALUE'
            RETURN
C            STOP
        ELSE
            RFI = 1.D0/P(33)
        ENDIF
      ELSEIF(I.LE.I2) THEN
        IF(P(35).LE.0.D0) THEN
            WRITE(*,*) 'P(35) MUST BE POSITIVE FOR THIS IRCH VALUE'
            RETURN
C            STOP
        ELSE
            RFI = 1.D0/P(35)
        ENDIF
      ELSEIF(I.LE.I3) THEN
        IF(P(37).LE.0.D0) THEN
            WRITE(*,*) 'P(37) MUST BE POSITIVE FOR THIS IRCH VALUE'
            RETURN
C            STOP
        ELSE
            RFI = 1.D0/P(37)
        ENDIF
      ELSE
        IF(P(39).LE.0.D0) THEN
            WRITE(*,*) 'P(39) MUST BE POSITIVE FOR THIS IRCH VALUE'
            RETURN
C            STOP
        ELSE
            RFI = 1.D0/P(39)
        ENDIF
      ENDIF
C
C     DO RECTANGULAR WEIGHTING OR MAGNITUDE AND PHASE WEIGHTING
      IF(IPF.EQ.0) THEN
                FMS = YQ(I)*YQ(I) + YQ(I+MD)*YQ(I+MD)
        YQRA = DABS(YQ(I))
        YQIA = DABS(YQ(I + MD))
        RT(I) = OAA*YQIA + OAB*YQRA + OAG*FMS*RFI
        RT(I + MD) = RT(I)
C     CALCULATE RECT. REAL AND IMAG. CORRESP. TO MAG. AND PHASE
C     ***  MAGNITUDE AND PHASE WEIGHTING NOT YET IMPLEMENTED  ***
      ENDIF
C
C     IPF= 0, RECT; =1: P; =2: D; =3: L.  SEE PFIT IN MANUAL
      IF(IPF.NE.0) THEN
        WRITE(*,*) 'THIS INPUT OPTION IS NOT YET IMPLEMENTED'
        RETURN
C        STOP
      ENDIF
C
      FT(I) =  RT(I)*RT(I)
      FT(I+MD) = RT(I + MD)*RT(I + MD)
878   CONTINUE
      GOTO 1677
872   WRITE(*,*) 'IRCH MUST BE -5 FOR DATTYP UNEQUAL TO C'
      RETURN
C      STOP
C
600   CONTINUE
C     SOLATRON #1250, PAR 273 WEIGHTING FOLLOWS
      GOTO 500
C
675   CONTINUE
C
      PUW = P(31)*P(31)
      IF(IXW.EQ.1) THEN
        AN1 = -1.D0/FLOAT(KY)
        YPD2 = 1.D0
        DO 8001 IB = 1,KY
          FT(IB) = PUW + FT(IB)
        IF(FT(IB).NE.0.D0) THEN
          YPD2 = YPD2*(FT(IB)**AN1)
        ELSE
            YPD2 = 0.D0
        ENDIF
8001    CONTINUE
C
        DO 717 IC = 1,KY
          FT(IC) = (DABS(FT(IC)))*YPD2
717     CONTINUE
C           DO 720 IB = 1,KY
        GOTO 1677
C
      ELSE
        DO 8002 IB = 1,KY
          FT(IB) = PUW + FT(IB)
C          R(IB) = DSQRT(FT(IB))
8002    CONTINUE
        GOTO 1677
      ENDIF
C
1677  CONTINUE
C
C     IF REQUIRED, DO WEIGHT SHAPING AT EXTREMES OF FREQUENCY
C
      IF(ISW.EQ.1) THEN
        DO 81 IX = 1,KY 
            RT(IX) = RT(IX)*SHW(IX)
            FT(IX) = RT(IX)*RT(IX)
81      CONTINUE
      ENDIF
C
      DO 8007 IB = 1,KY
        FT(IB) = DSQRT(FT(IB))
8007  CONTINUE
      RETURN
      END
