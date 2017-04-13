C      FIRST PART: LEVM CNLS PROGRAM TO RUN ON MICROSOFT FORTRAN
C
C                J.R. Macdonald
C            LEVM PROGRAM: MAIN -- LEV0.FOR
C
C  COPY COMMENTED BY JS
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 BLANK,DD,EE,FF,LM,MM,NN,PP,YY,KZ,SS,
     + DATTYP,DINP,DFIT,PINP,PFIT,FREEQ,NEG,FUN,DATTYQ,QQ*4,ALPHA*80,
     + RR,ZZ,DATTYC,OO
      INCLUDE 'SIZE.INC'
C
C    *****  NOW HANDLES UP TO NPTS DATA POINTS AND NPAFR FREE PARAM *****
C    !!!!!  IF ANY VALUES IN SIZE.INC ARE CHANGED, YOU MUST RECOMPILE
C        LV0.FOR, LV1.FOR, LV2.FOR, AND LV5.FOR   !!!!!
C
      EXTERNAL DISTEL
      DIMENSION X(NPAFR),FV1(NPT2),PEX(NTOT),TFREQ(NPTS),PSX(NTOT),
     +  PS(NTOT),NFRES(NTOT),NFREX(NTOT),QX(7),FV2(NPT2),OFREQ(NPTS)
      COMMON/CM1/FREQ(NPTS),MD,DATTYP
      COMMON /CM2/ Y(NPT2),R(NPT2),FJ(NPT2),P(NTOT),DRSS,ROE,RKE,
     +NS(NPAFR),NFREE(NTOT),NP,ICNT,MN,IRCH,IXI,DATTYQ
      COMMON/CM3/CELCAP,FUN,DFIT,PFIT
      COMMON /CM10/EPSG,IZR
      COMMON /CM11/ MQY,ISPR,ICX,NDF,FQQ
      COMMON /CM12/ CLCAP,ATEMP,WF,MAXFEV,ICF,MDE,JCDX
      COMMON /CM13/ RX,TX,UU,PHI,XXM1,XX1,XX2,XX3,RN,AIN,ICAV,NELEM,NCH
      COMMON/CM16/IOP,IORIG,NYC,JY,IPRINT,LDFJAC,MODE,IFP,IRE,ISTP,JFP,
     +NPH,INE
      COMMON /CM18/ SDWC,SDRC,DIAG(NPAFR),IPAR,IOPR
      COMMON /CM34/ MDA,IWT,IXW,INFO,IPL
      COMMON /CM36/ SHW(NPT2),ISW
      COMMON /CM78/ DATTYC
      COMMON /CM79/ YTT(NPT2)
      COMMON /CM35/ JIT,IPF,NPRINT
      COMMON /CM39/ P8,P18
      COMMON /CM53/ IRG,IPRIN,STRN(20),AI
      COMPLEX*16 FREQIM,TRANS,DISTEL
      DATA BLANK/' '/,DD/'D'/,EE/'E'/,FF/'F'/,KZ/'K'/,
     + MM/'M'/,NN/'N'/,PP/'P'/,SS/'S'/,TT/'T'/,
     + YY/'Y'/,LM/'L'/,OO/'O'/,RR/'R'/,ONENEG/-1.D0/,
     + ZZ/'Z'/,SML/1.D-40/
      DATA  TUPI/6.283185307179586D0/,TPDEG/1.7453292519943296D-2/,
     + TPII/1.5915494309189534D-1/,EPV/8.85418782D-14/
C 
C     OPEN THE PRINT OUTPUT FILES
C     SET UNIT 3 TO BE THE MAIN PRINTOUT AND
C     UNIT 4 TO BE THE AUXILIARY PRINTOUT
C       MD IS NUMBER OF DATA POINTS
C
      OPEN(3,FILE='PNTOUTL')
      OPEN(4,FILE='AUXPNTL')
C
C     FILE 1 IS THE INPUT DATA FILE (INFL)
      OPEN(1,FILE='INFL')
C
                IOCNT = 0
                ICAV = 0
                JCDX = 0
                INDE = 2
      WRITE(3,20)
      WRITE(4,20)
C
   20 FORMAT(' LEVM : COMPLEX NONLINEAR LEAST SQUARES IMMITTANCE ',
     1       'DATA FITTING PROGRAM',//15X,'VERSION 8.06 - 2/05'//)
C
C     READ AND LIST INPUT INFORMATION
C
C   LINE 1: CHARACTER DESCRIPTION OF RUN
C
      READ(1,50) ALPHA
   50 FORMAT(A80)
      WRITE(3,55) ALPHA
      WRITE(4,55) ALPHA
   55 FORMAT(A80)
C
C   LINE 2: DATA AND FIT TYPES
C
      READ(1,51) IOP,DINP,DFIT,PINP,PFIT,FREEQ,NEG,FUN,CELCAP,DATTYP,
     + IPAR,ROE,IFP,IRE
   51 FORMAT(I4,2X,7A1,D10.4,A1,I10,1D10.4,2I6,1X,D10.4,2X,D10.4)
                DATTYQ = DATTYP
                DATTYC = DATTYP
        IF(PFIT.EQ.'R') THEN
                IPF = 0
        ELSEIF(PFIT.EQ.'P') THEN
                IPF = 1
        ELSEIF(PFIT.EQ.'D') THEN
                IPF = 2
        ELSEIF(PFIT.EQ.'L') THEN
                IPF = 3
        ENDIF
C
         WRITE(3,53) DINP,PINP,DFIT,PFIT,FUN
         WRITE(4,53) DINP,PINP,DFIT,PFIT,FUN
   53    FORMAT(/' DATA ENTERED IN ',2A1,' FORMAT TO BE USED IN ',
     1         2A1,' FIT;      CIRCUIT MODEL : ',A1)
C
C       SET FOR "EXACT" INPUT PARAMETER REL SD CALCULATIONS
        IF(IPAR.LT.0) THEN
                IORIG = 1
        ELSE
                IORIG = 0
        ENDIF
C 
C   LINE 3: DATA, WEIGHTING, AND FITTING SPECIFICATIONS
C
      READ(1,60) MD,NP,MAXFEV,NPRINT,IRCH,MODE,ICP,IPRINT,IGACC,ATEMP
60    FORMAT(9I5,D10.4)
C
        IF(MD.LT.0) THEN
                MD = IABS(MD)
                IMD = -1
        ELSE
                IMD = 1
        ENDIF   
C
        MDE = MODE
        IF(NP.LT.0) THEN
                NP = -NP
                INP = 1
        ELSE
                INP = 0
        ENDIF
        IF(NP.LT.32) NP = 32
        IF(NP.GE.NTOT) NP = NTOT
        IF(NPAFR.GT.32) THEN
                NZO = NPAFR
        ELSE
                NZO = 40
        ENDIF
C
                IPRIN = IPRINT
        IF(IPRINT.LT.0) THEN
                IPRINT = -IPRINT
C               IPL = 1
                IPL = IPRINT
        ELSE
                IPL = 0
        ENDIF
C
        IZR = 0
        IF(DATTYP.EQ.'C'.AND.DINP.EQ.ZZ.AND.DFIT.EQ.ZZ.AND.PINP.EQ.RR.
     +  AND.PFIT.EQ.RR.AND.INP.EQ.1) IZR = 1
C
        MDA = MD
        NRCH = IRCH
C
                IFPO = IFP
        IF(IFP.LT.0) THEN
                JFP = 0
                IFP = -IFP
        ELSE
                JFP = 1
        ENDIF
C
        IF(IRCH.GE.0) THEN
                IWT = 0
                IF(IRCH.LT.2) THEN
                        QQ = 'UNIT'
                ELSE
                        QQ = 'DATA'
                ENDIF
        ELSE
                IWT = 1
                IRCH = -IRCH
                QQ = 'FUNC'
        ENDIF
C
C        IF DATTYP IS C, FIT COMPLEX DATA; IF IT IS R, FIT REAL DATA;
C        IF IT IS I, FIT IMAGINARY DATA
C
                KY = 2*MD 
        IF(IOP.LT.0) THEN
                IACY = -IOP
                IOPO = IOP
                IOP = 0
        ELSE
                IACY = 0
        ENDIF
       IF(DATTYP.NE.'C') IOP = 0
       IF(DATTYP.EQ.'C') THEN
                JY = 1
                WRITE(3,61)
                WRITE(4,61)
61              FORMAT(/12X,'********  FIT OF COMPLEX DATA  ********')
        ELSE IF (DATTYP.EQ.'R') THEN
                JY = 1
                KY = MD
                WRITE(3,62)
                WRITE(4,62)
62              FORMAT(/12X,' *****  FIT OF REAL PART OF DATA  *****')
        ELSE IF (DATTYP.EQ.'I') THEN
                JY = 1 + MD
                WRITE(3,63)
                WRITE(4,63)
63              FORMAT(/12X,'*****  FIT IMAGINARY PART OF DATA  *****')
        ELSE
                WRITE(3,64)
                WRITE(4,64)
64              FORMAT(/' ERROR: IMPROPER CHOICE OF DATTYP',
     +                  ' FITTING PARAMETER; USE C,R, or I')
        CALL STOPX
        ENDIF
C
      NYC=0
C
C     LINES 4: MODEL PARAMETERS
C
      READ(1,110) ( P(I), I = 1, NP )
110   FORMAT(5D15.8)
C
	    P8 = P(8)
	    P18 = P(18)
  
        IRG = 0
        IF(P(31).LT.0.D0) THEN
                IRG = INT(DABS(P(31))+.01)
                P(31) = 0.D0
        ENDIF
C
C       SET PEX VALUES AS "EXACT" INPUTS **************
        IF(IORIG.EQ.1) THEN
                DO 666 ILM = 1,NTOT
                PEX(ILM) = P(ILM)
666             CONTINUE
        ENDIF
C
C   LINE 5: PARAMETER TYPES - FREE OR FIXED
C
      READ(1,120) ( NFREE(I), I = 1, NTOT )
  120 FORMAT(50I1)
C
C       MAKE A VECTOR OF THE FREE PARAMETERS (X) AND INDICES (NS)
C
        IF(NFREE(31).NE.0.AND.NFREE(32).EQ.0) THEN
                PAUSE 'ILLEGAL CHOICE OF FREE/FIX FOR U AND XI; STOP!'
                STOP
        ENDIF
        IF(IRCH.LT.2.AND.((NFREE(31) + NFREE(32)).NE.0)) THEN
                PAUSE 'BAD CHOICE OF FREE/FIX FOR U AND XI; DO ''go''
     1 to RESET AND CONTINUE'
                NFREE(31) = 0
                NFREE(32) = 0
        ENDIF
        IXI = 0
        IF(NFREE(32).GT.0.AND.IRCH.GT.1) IXI = 1
        LLL = 0
        LS = 0
      DO 128 I = 1, NTOT
         IF ( NFREE(I) .EQ. 0 ) GOTO 128
         LLL = LLL + 1
         NS(LLL) = I
         X(LLL) = P(I)
128   CONTINUE
C
      IF(FUN.EQ.RR.OR.FUN.EQ.KZ) THEN
        INE = 0
        ING = 0
        IF(NFREE(29).GT.0) INE = 1
        IF(NFREE(30).GT.0) ING = 1
        INE = INE + ING
        ENDIF
        NFREI = LLL
C
        NFIX=NP-NFREI
        INS = 0
        IXW = 0
        IF(IXI.EQ.1.OR.IWT.EQ.1) IXW = 1
        IF(NFREE(31).EQ.0.AND.NFREE(32).NE.0) INS = 1
        IF(NFREE(31).NE.0.AND.NFREE(32).NE.0) INS = 2
        MQY = NFREI - INS
C
        IF(MAXFEV.LE.2) GOTO 457
C
C       ENTER NEW INITIAL VALUES OF FREE PARAMETERS IF ICP .NE. 0
C
        IF(ICP.EQ.0) GOTO 452
        DO 455 IJ = 1,NFREI
                WRITE(*,458) IJ,NS(IJ),P(NS(IJ))
458             FORMAT(2X, 'INDEX #:',I3,4X,'PARAMETER #:',I3,
     14X,'PARAMETER=',1P,1D16.7/)                       
        WRITE(*,464)
464     FORMAT(2X,'IF NO CHANGE IN PARAMETER VALUE DESIRED, ENTER 99',/
     12X,'OTHERWISE, ENTER NEW VALUE OF PARAMETER BUT DO NOT USE D DESCR
     2IPTOR')
                READ(*,*) X(IJ)
                IF(IDINT(X(IJ)+.01).EQ.99) X(IJ) = P(NS(IJ))
455     CONTINUE
C
C         FOR |MODE|=4  USE INPUT MODEL PARAMETER SCALING
C
452     IF(ABS(MODE).NE.4) GO TO 457
        DO 456 IJ = 1,NFREI
        WRITE(*,'(2X,''ENTER VALUE OF ISC=LOG10 OF SCALE FACTOR:'')')
                WRITE(*,*)
                WRITE(*,458) IJ,NS(IJ),P(NS(IJ))
                READ(*,*) ISC
                DIAG(IJ) = 1.D1**(-ISC)                                 
456     CONTINUE
457     CONTINUE        
C
C   SET DEFAULT VALUES
C
        IF (CELCAP .EQ. 0.D0) CELCAP = 1.D0
        IF(IMD.EQ.-1) CELCAP = EPV
        IF(CELCAP.GT.0) THEN
                CLCAP = CELCAP
        ELSE
                CLCAP = 1.D0
        ENDIF
C
C       SET CHOICE TO PRINT CONDUCTIVITY (IN LINOTCON.DAT) AS WELL
C       AS ADMITTANCE
                ICOND = 0
          IPATM = INT(DABS(ATEMP) + 0.01D0)
          IPATA = SIGN(1.D0,ATEMP)
          IPTI = IPATA*IPATM
        IF(IPTI.EQ.-2.AND.CELCAP.NE.1.D0.AND.DFIT.EQ.YY) THEN
                ICOND = 1
                FC = EPV/CELCAP
        ELSEIF(IPTI.EQ.-3.AND.DATTYP.EQ.'R') THEN
                ICOND = 2
        ENDIF
C
        DRSS = 1.0D-8
        GTOL = 0.D0
        IF(IACY.EQ.0) THEN      
                FTOL = 1.D-30
                XTOL = 1.D-48
C               XTOL = 0.D0
        ELSE
                FTOL = 1.D1**(-IACY)
                XTOL = FTOL
        ENDIF
C               
        LDFJAC = KY
C
C   OUTPUT STATISTICS FOR RUN
C
      WRITE(3,200) MD, NRCH, NFREI
  200 FORMAT('- # OF DATA POINTS =',I3,5X,'WEIGHT: IRCH =',I3,5X,
     +       '# OF FREE PARAMETERS =',I2)
      WRITE(3,205) NPRINT, MAXFEV,QQ
  205 FORMAT ('- PRINTS EVERY',I3,' ITERATIONS',5X,'MAX # ITERATIONS=',
     +       I3,5X,'MAIN WT USES: ',A4,/)
C
C   OUTPUT WEIGHTING STATISTICS
C
      IF(IXI.EQ.0) THEN 
        IF(IRCH.EQ.1) WRITE(3,259)
        IF(IRCH.EQ.1) WRITE(4,259)
        IF(IRCH.EQ.2) WRITE(3,260) P(32)
        IF(IRCH.EQ.2) WRITE(4,260) P(32)
        IF(IRCH.EQ.0) WRITE(3,261)
        IF(IRCH.EQ.0) WRITE(4,261)
        IF(IRCH.EQ.3) WRITE(3,264) P(32)
        IF(IRCH.EQ.3) WRITE(4,264) P(32)
        IF(IRCH.EQ.4) WRITE(3,2644)
        IF(IRCH.EQ.4) WRITE(4,2644)
        IF(IRCH.EQ.5) WRITE(3,2645)
        IF(IRCH.EQ.5) WRITE(4,2645)
        IF(IRCH.EQ.6) WRITE(3,2646)
        IF(IRCH.EQ.6) WRITE(4,2646)
      ELSE
        WRITE(3,2641)
        WRITE(4,2641)
      ENDIF
259   FORMAT(15X,'UNIT WEIGHTING ASSIGNED TO EACH POINT'/)
260   FORMAT(5X,'WEIGHTS INVOLVE THE MAGNITUDES OF DATA OR FUNCTION VALU
     1ES RAISED'/,5X,'TO THE POWER',1P,(D12.5)/)
261   FORMAT (15X,'WEIGHTS (STANDARD DEVIATIONS) ARE READ IN'/) 
264   FORMAT(5X,' MODULUS WEIGHTING: RESULTS RAISED TO THE POWER',
     21P,(D12.5)/)
2641  FORMAT(5X,'WEIGHTS USE XI OR U**2 AND XI: BOTH MAY BE FREE PARAMET
     1ERS'/)
2644  FORMAT(5X,'#1174 SPINOLO FRA WEIGHTING'/)
2645  FORMAT(5X,'#1250 & 1286 ORAZEM-AGARWAL FRA WEIGHTING'/)
2646  FORMAT(5X,'#1250 & PAR 273 ORAZEM-AGARWAL FRA WEIGHTING'/)
C
      WRITE(3,265)
      WRITE(4,265)
265   FORMAT(' INITIAL PARAMETER GUESSES AND FIXED (0) OR FREE (1 OR 2)
     1STATUS')
C       OUTPUT OF INPUT FOR I = 1 TO 16, 17 T0 34, and 41 TO 50 ******
      DO 330 I = 1, 16
                JO = 16 + I
         WRITE(3,310) I,P(I),NFREE(I),JO,P(JO),NFREE(JO)
         WRITE(4,310) I,P(I),NFREE(I),JO,P(JO),NFREE(JO)
330   CONTINUE
        IF(NP.GT.40) THEN
                WRITE(3,*)
                WRITE(4,*)
           DO 338 I = 41,45
                JO = 5 + I
                WRITE(3,310) I,P(I),NFREE(I),JO,P(JO),NFREE(JO)
                WRITE(4,310) I,P(I),NFREE(I),JO,P(JO),NFREE(JO)
338        CONTINUE
        ENDIF
C
      IF(NP.GT.32) THEN
        WRITE(3,313)
        WRITE(4,313)
        DO 397 I = 33, 36
                 JO = 4 + I
                 WRITE(3,312) I,P(I),JO,P(JO)
                 WRITE(4,312) I,P(I),JO,P(JO)
397     CONTINUE
      ENDIF
C
  310    FORMAT(' ',4X,'P(',I2,') = ',1PD15.6,I5,
     18X,'P(',I2,') = ',1PD15.6,I5)
  312    FORMAT(' ',4X,'P(',I2,') = ',1PD15.6,13X,
     1'P(',I2,') = ',1PD15.6)
313     FORMAT(/15X,'THE FOLLOWING PARAMETERS ARE ALWAYS FIXED')
C
C   LINES 6: DATA TO BE FIT
C
      DO 380 I = 1, MD
         IM = I + MD
         READ(1,370) IC, FREQ(I), Y(I), Y(IM)
C
	 IF(ROE.LT.0.D0) THEN	  !TO CHANGE TO SPECIFIC DATA
		CLC= DABS(ROE)
C	 WRITE(*,*) CLC,Y(I),Y(IM)
		Y(I)=CLC*Y(I)
		Y(IM)=CLC*Y(IM)
C	 WRITE(*,*) MD,Y(I),Y(IM)
	 ENDIF


C	 PAUSE 	 	
C		
C    CHANGE SIGN OF IMAGINARY PART OF RECTANGULAR DATA IF SPECIFIED
C
        IF (NEG.EQ.NN) Y(IM) = -Y(IM)
C               YTT(I) = Y(I)
                YTT(IM) = Y(IM)
C
CC      DO I=1,MD
CCC     IF(I.LT.4.OR.I.GT.146) THEN
CCC             IF(I.EQ.2) WRITE(*,*) MD,MDS,KY,MAXFEV,NPT2,NPTS
CCC                     WRITE(*,370) I, FREQ(I), Y(I), Y(IM)
CCC     ENDIF
CC      ENDDO
C XXX    IF(I.LT.3.OR.I.EQ.701.OR.I.EQ.702.OR.I.GT.702) THEN 
C XXX            IF(I.EQ.2) WRITE(*,*) MD,MDS,KY,MAXFEV,NPT2,NPTS
C XXX                    WRITE(*,370) I, FREQ(I), YTT(I), YTT(IM)
C XXX    ENDIF
C
380   CONTINUE
CC
	IF(ROE.LT.0.D0) ROE=0.D0
C XXX    READ(*,*)
CC
370   FORMAT(I5,3(5X,D20.13))
                NEG = BLANK
3702     FORMAT(10X,E20.13)
C
C   *** END OF INPUT LINES (except LINES 7, read later in PGM) ***
C
C               &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C
        IF(FUN.EQ.OO.AND.P(10).EQ.33D0) THEN
                IP10T = 1
           ELSE
                IP10T = 0
        ENDIF
C
      IF(FUN.EQ.KZ.OR.FUN.EQ.RR.OR.FUN.EQ.SS) THEN
          IF(FUN.EQ.RR.OR.FUN.EQ.SS) THEN
                JIA = 28
                JIB = 12
          ELSE
                JIA = 22
                JIB = 18
          ENDIF
C
          IP34A = INT(DABS(P(34)) + 0.01D0)
          IP34 = SIGN(1.D0,P(34))*IP34A
          IP35A = INT(DABS(P(35)) + 0.01D0)
          IP38A = INT(DABS(P(38)) + 0.01D0)
          IP40A = INT(DABS(P(40)) + 0.01D0)
          IP40S = SIGN(1.D0,P(40))
C
        IF(P(34).GT.7.1D0) THEN
                PAUSE "P(34) TOO LARGE"
                STOP
        ENDIF
C
        IF(IRE.GT.-1.D1) IRE = -10
        IF(IP35A.GT.0.D0) THEN
                IF(IRCH.GT.3) THEN
                        PAUSE 'IRCH MUST BE < 4 HERE'
                        STOP
                ENDIF
                ICF = 2.D0*P(35)
C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       ELSE
                ICF = 0
        DO KK = 1, NTOT
          IF((KK.LE.JIA.OR.KK.GT.40).AND.P(KK).NE.0.D0) ICF = ICF + 1
        ENDDO
       ENDIF    
C
        NPH = ICF/2
        IKW = 1
        IF(FUN.EQ.KZ.AND.NPH.GT.11) THEN
                PAUSE "TOO MANY DRT PAIRS FOR KSUB; USE RSUB"
                STOP
        ENDIF
C
        IF(IP34.GT.0) IP38A = 1
        IF(IP38A.EQ.1) THEN
                  IF(P(36).EQ.0) THEN
                        IF(FREQ(MD).GT.FREQ(1)) THEN
                                TAUMIN = 1.D0/FREQ(MD)
                                TAUMAX = 1.D0/FREQ(1)
                        ELSE
                                TAUMIN = 1.D0/FREQ(1)
                                TAUMAX = FREQ(MD)
                        ENDIF
                                TAUM = TAUMIN
                                FRAT = TAUMAX/TAUMIN
                  ELSE
                                TAUM = P(36)
                                FRAT = P(39)/P(36)
                  ENDIF

                        WF = (1.D0/DFLOAT(NPH - 1))*DLOG(FRAT)

        ELSE
                IF(IP40A.LE.1.AND.IP34.LE.0)  THEN
                        GOTO 783
                ENDIF
        ENDIF
       IF(IP34.GT.0.AND.IP34.LT.8) THEN
                MAXFEV = 1
                IRCH = 2
                IKW = 0
C
                IF(IP34.LT.3.OR.IP34.EQ.5) THEN
                                P(38) = 1
                        IF(IP40A.NE.4)  P(40) = 3*IP40S
                ELSEIF(IP34.EQ.3.OR.IP34.EQ.4) THEN
                        IF(IP40A.NE.4) P(40) = 2*IP40S
                        P(38) = 2
C       CHANGE R,TAU TO R,C
                ELSEIF(IP34.EQ.6) THEN
C                       P(40) = IP40S
                        P(38) = 2
C       CHANGE R,C TO R,TAU
                ELSEIF(IP34.EQ.7) THEN
                        P(40) = 2*IP40S
                        P(38) = 2
                ENDIF
C               
       ENDIF
781             IF(MOD(ICF,2).NE.0.AND.IP35A.GE.0) THEN
                        WRITE(*,*) ICF,IP35A
                        WRITE(*,*)
                        PAUSE "# OF DRT-PAIR PARAMETERS MUST BE EVEN"
                        PAUSE "SET P(35) TO TOTAL NUMBER OF PAIRS"
                        STOP
                ENDIF
C
        IF(FUN.NE.SS) THEN 
        CALL RESORT(P,NFREE,TAUM,WF,X,NTOT,NS,JIA,JIB)
        ENDIF
C
783             IF(ICF.LT.2) THEN
                        PAUSE "# OF FREE DRT PARAMETERS MUST BE >= 2"
                        STOP
                ENDIF
      ENDIF
C
C                       *************************************
C
                IOPR = 0
        IF(IRE.LT.0) THEN
                IDMOD = 1
                OPEN(9,FILE='OUTIN')
        IF(IRE.LT.-9) THEN
                        IDMOD = 2
C           OPEN FILES FOR PLOT INPUT: DATA AND FIT RESULTS
                        OPEN(2,FILE='LINOUT')
                        OPEN(22,FILE='LVOUT')
                        OPEN(17,FILE='LINOTR')
                        OPEN(18,FILE='LINOTL')
                        OPEN(21,FILE='LINOTLN')
                        IF(ICOND.NE.0) OPEN(19,FILE='LINOTCON')
                        OPEN(20,FILE='LINOTRL.DAT')
        ENDIF
C       IRE=-11, OUTPUT COLUMN IS FREQUENCY, NOT ANGULAR FREQUENCY
                IF(IRE.LE.-11) IOPR = 1
        ELSE
                IDMOD = 0
        ENDIF
C
C       &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
C       PREPARE WEIGHT SHAPING IF REQUESTED
C
        IF(IGACC.LT.0) THEN
                ISW = 1
C
                WRITE(*,'(2X,''ENTER WT SHAPING FREQUENCIES: FL1, FL2, F
     1H1, FH2, ALL ON ONE LINE:'')')
                READ(*,*) FL1,FL2,FH1,FH2
                WRITE(3,764) FL1,FL2,FH1,FH2
764             FORMAT(/1X,'WEIGHT SHAPING FREQS:',1P,4E12.4)
                UKL2 = (FL2)**0.25
                UKH1 = 1.D0/(FH1)**0.25
                USL12 = (FL1*FL2)**0.25
                UH12 = DSQRT(FH1*FH2)
                USH12 = 1.D0/DSQRT(UH12)
C
                DO 1007 IH = 1,MD
                        UIX = DSQRT(FREQ(IH))
                        UIXS = DSQRT(UIX)
                        UIXI = 1.D0/UIX
                        UIXSI = 1.D0/UIXS
                   IF(FREQ(IH).LE.FL1) THEN
                                SHW(IH) = USL12*UIXI
                   ELSEIF(FREQ(IH).LE.FL2) THEN
                                SHW(IH) = UKL2*UIXSI
                   ELSEIF(FREQ(IH).LE.FH1) THEN
                                SHW(IH) = 1.D0
                   ELSEIF(FREQ(IH).LE.FH2) THEN
                                SHW(IH) = UKH1*UIXS
                   ELSEIF(FREQ(IH).GT.FH2) THEN
                                SHW(IH) = USH12*UIX                             
                   ENDIF        
                   SHW(IH+MD) = SHW(IH)
1007            CONTINUE
        ELSE
                ISW = 0
        ENDIF
                EPSG = 1.D1**(-ABS(IGACC))
C
C   IF FREQUENCY IS IN HERTZ, CONVERT TO ANGULAR FREQUENCY
C
      IF (FREEQ.EQ.FF) THEN
         IF (IPRINT.GE.2) WRITE(3,410) DINP,DINP
  410    FORMAT ('-',10X,'FREQUENCY',17X,A1,'(I)',20X,A1,'(I+MD)')
         IF (IPRINT.GE.2) WRITE(3,420) (FREQ(I),Y(I),Y(I+MD), I=1,MD)
420      FORMAT (3(5X,1PD20.13))
         DO 430 I=1,MD
          FREQ(I) = TUPI*FREQ(I)        ! FREQ IN OMEGA r/s
430      CONTINUE
      ELSE
         IF (IPRINT.GT.2) THEN 
             WRITE(3,460) DINP,DINP
             WRITE(3,420) (FREQ(I),Y(I),Y(I+MD), I=1,MD)
         ENDIF
      END IF
  460    FORMAT(12X,'OMEGA',17X,A1,'(I)',17X,A1,'(I+MD)')
C
C       CONTINUE CONSTRUCTING FILE OUTIN
        IF(IRE.LT.0.AND.FUN.NE.OO.AND.MAXFEV.GT.2) THEN
          IF(IACY.NE.0) IOP = IOPO
C
        ENDIF
C
C   TRANSFORM INPUT DATA TO FIT FORMAT IF NECESSARY
C
      IF ((DINP.EQ.DFIT) .AND. (PINP.EQ.PFIT)) THEN
         WRITE(3,500) DINP,PINP
  500    FORMAT(/'-NO TRANSFORMATION MADE: DATA IN ',2A1,' FORMAT')
      ELSE
C
C   CHANGE ALL DATA TO RECTANGULAR FORMAT
C
        IF((PINP.EQ.PP).OR.(PINP.EQ.DD).OR.(PINP.EQ.LM)) THEN
            DO 510 I=1,MD
              IM = I + MD
              ARG = Y(IM)
              IF(PINP.EQ.LM) Y(I) = DEXP(Y(I))
              IF (PINP.EQ.DD) ARG = ARG*TPDEG
              RX = Y(I)*DCOS(ARG)
              RY = Y(I)*DSIN(ARG)
              Y(I) = RX
              Y(IM) = RY
  510       CONTINUE
        END IF
C
C   CHANGE ALL DATA TO Z FORMAT
C
C      CHANGE  Y -> Z, E -> M
C
         IF ((DINP.EQ.YY) .OR. (DINP.EQ.EE)) THEN
            DO 520 I=1, MD
              IM = I + MD
              TRANS = DCMPLX(Y(I), Y(IM))
              TRANS = 1.D0/TRANS
              Y(I) = DREAL(TRANS)
              Y(IM) = DIMAG(TRANS)
  520       CONTINUE
         END IF
C
C      CHANGE  M (OR E) -> Z
C
         IF ((DINP.EQ.EE) .OR. (DINP.EQ.MM)) THEN
            DO 515 I=1, MD
              IM = I + MD
              FREQIM = DCMPLX(0.D0, FREQ(I)*CLCAP)
              TRANS = DCMPLX(Y(I), Y(IM))
              TRANS = TRANS/FREQIM
              Y(I) = DREAL(TRANS)
              Y(IM) = DIMAG(TRANS)
  515       CONTINUE
         END IF
C
C   NOW THAT DATA IS Z, CHANGE IT TO FITTING FORMAT
C
C      CHANGE  Z -> M  IF APPROPRIATE
C
         IF ((DFIT.EQ.EE) .OR. (DFIT.EQ.MM)) THEN
            DO 540 I=1, MD
              IM = I + MD
              FREQIM = DCMPLX(0.D0,FREQ(I)*CLCAP)
        FREQIM =FREQIM
              TRANS = DCMPLX(Y(I), Y(IM))*FREQIM
              Y(I) = DREAL(TRANS)
              Y(IM) = DIMAG(TRANS)
540         CONTINUE
         END IF
C
C      CHANGE  Z -> Y, M -> E  IF APPROPRIATE
C
         IF ((DFIT.EQ.YY) .OR. (DFIT.EQ.EE)) THEN
            DO 550 I=1, MD
              IM = I + MD
              TRANS = DCMPLX(Y(I), Y(IM))
              TRANS = 1.D0/TRANS
              Y(I) = DREAL(TRANS)
              Y(IM) = DIMAG(TRANS)
 550       CONTINUE
         END IF
C
C   CHANGE TO POLAR DATA IF APPROPRIATE
C
         IF ((PFIT.EQ.PP).OR.(PFIT.EQ.DD).OR.(PFIT.EQ.LM)) THEN
            DO 560 I=1,MD
              IM = I + MD
              YR = Y(I)
              YI = Y(IM)
              Y(I) = DSQRT(YR*YR + YI*YI)
              Y(IM) = DATAN(YI/YR)
              IF(PFIT.EQ.LM) Y(I) = DLOG(Y(I))
              IF (PFIT.EQ.DD) Y(IM) = Y(IM)/TPDEG
  560       CONTINUE
         END IF
         WRITE(3,570) DINP,PINP,DFIT,PFIT,CLCAP
  570    FORMAT(/'-DATA TRANSFORMED: ',2A1,' --> ',2A1,
     1          '      CELL CAPACITANCE:',D10.4,/)
      END IF
C
        DO KD = 1,MD
                TFREQ(KD) = TPII*FREQ(KD)
                OFREQ(KD) = FREQ(KD)
           IF(IDMOD.EQ.2) THEN
                IF(IOPR.EQ.1) THEN
                        WRITE(2,1020) TFREQ(KD),Y(KD),Y(KD+MD)
                ELSE
                        WRITE(2,1020) FREQ(KD),Y(KD),Y(KD+MD)
                ENDIF
           ENDIF
        ENDDO
C
                IF(IDMOD.EQ.2) WRITE(2,3702) ONENEG
C
C   *** TRANSFORMATION COMPLETE ***
C
        NFEV = 2
C
        IF(FUN.EQ.OO) THEN
                IF(P(35).EQ.4.D0.AND.P(40).NE.4.D0) THEN
                        NELEM = P(10)
                ELSEIF(P(35).NE.4.D0.AND.P(40).EQ.4.D0) THEN
                        NELEM = P(20)
                ELSE
                        NELEM = P(10)
                ENDIF
        ENDIF
C
        NELEM0 = NELEM
        IF(FUN.EQ.OO.AND.(NELEM.EQ.7.OR.NELEM.EQ.10.OR.NELEM.EQ.32
     *  .OR.NELEM.EQ.36)) THEN
                INDE = 0
          ELSEIF(IRE.LT.0) THEN
                INDE = 1
          ELSE
                INDE = 2
        ENDIF
C
        ECINF = 0.D0
        EC0 = 0.D0
C
        IF(NELEM.EQ.7.AND.P(7).NE.1.D0) THEN
                P(10) = 6.D0
                NELEM = 6
                INDE = 1
        ENDIF
C
C
        MDS = MD*IMD
C
C   DETERMINE THE (FIXED) WEIGHTING FACTORS FOR IRCH = 0
C    HERE R IS A VECTOR(I=1,MD) OF THE STANDARD DEVIATIONS ASSOCIATED
C    WITH EACH DATA POINT
C
        IF(IRCH.EQ.0) THEN
                 DO 653 I=1,MD
C    READ LINES 7 (if they exist)
         READ(1,*) IC,QRXX,PRYY
C                       
        R(I) = QRXX
        R(I+MD) = PRYY
653     CONTINUE
652    FORMAT(I5,25X,2(5X,D20.13))
                IF(R(1).EQ.0) THEN
                        WRITE(*,736)
736     FORMAT(3X,'FOR IRCH=0, YOU MUST INCLUDE INPUTS IN LINES 7')
                        CALL STOPX
                ENDIF   
        ENDIF
C
C   DON'T CALL FITTING ROUTINE IF < 3 ITERATIONS ARE SPECIFIED
      IF (MAXFEV.LE.2) GO TO 855
C
C
C           ***********   ################   *********
C
C               CALL MAIN FITTING ROUTINES
C
C           ***********   ################   *********
C
CCC     PAUSE "PREMAIN"
      CALL MAINCLC(KY,FTOL,GTOL,XTOL,X,MAXFEV,NPRINT,NFEV,PEX,
     + NFREI,FV1)
C
CCC     PAUSE "MAIN"
      IF(FUN.EQ.OO.AND.P(35).EQ.4) THEN
              ICAV = 1
              CALL MODEL(NZO,P,FV2)
              ICAV = 0
      ENDIF
C
C       Make output-input file OUTIN if IRE < 0.
C
C       CONTINUE CONSTRUCTING FILE OUTIN
C
      IF(FUN.EQ.OO.AND.IRE.LT.0) THEN
         IF(P(5).GT.0.D0) THEN
              ETAU = P(6)/(P(5)*CELCAP)
           IF(XXM1.EQ.0.D0) THEN
              ECINF = 0.D0
           ELSE
              ECINF = ETAU/XXM1
           ENDIF
              EC0 = ETAU*XX1
         ELSE
              ECINF = 0.D0
              EC0 = 0.D0
         ENDIF
      ENDIF
C
      IF(INDE.EQ.0) THEN
              CONTINUE
          ELSEIF((NELEM.EQ.7.AND.P(6).NE.1.D0).OR.P(40).EQ.4.D0) THEN                  
              ECINF = 0.D0
              EC0 = 0.D0
      ENDIF
C
      IF(IDMOD.GE.1) THEN
              DO 3755 I = 1,NFREI
                      P(NS(I)) = X(I)
3755          CONTINUE
C
C     RESET FREQUENCY COLUMN TO Hz IF FREEQ=FF ORIGINALLY
      ENDIF
C
      IF(FREEQ.EQ.FF) THEN
         DO II =1, MD
              FREQ(II) = TFREQ(II)
         END DO
      ENDIF
C
69    FORMAT(9I5,D10.4,2X,D10.4,2X,D10.4)
C
      IF(NELEM.LT.7) GOTO 973
      IF(NELEM.EQ.7.AND.(P(7).LT.1.D0.AND.P(9).LT.1.D0)) GO TO 973
C
c       @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
855   CONTINUE       
C
C   ENTRY POINT FOR MAXFEV = 0,1,2
C       FOR THESE CHOICES, A NEW INPUT FILE IS PRODUCED
C  0    IF MAXFEV=0, NO FIT BUT CALCULATES NEW DATA BASED ON PARAMETER VALUES
C
C  1    IF MAXFEV=1, NO FIT BUT IGNORES PARAMETER VALUES AND CALCULATES DATA
C       AT NEW LEVEL, BASED ON VALUE OF DFIT INPUT CHOICE.  USE TO CONVERT TO
C       DIFFERENT DATA LEVEL (EPSILON TO Z, ETC.  SEE CELCAP INPUT DISCUSSION
C
C  2    IF MAXFEV=2, NO FIT BUT CALCULATES NEW Y DATA WITHOUT PARAMETRS WHOSE
C       NFREE VALUES WERE SET TO 3
C
        IF(MAXFEV.LT.2.AND.FUN.EQ.OO.AND.P(35).EQ.4) THEN
                        ICAV = 1
                        CALL MODEL(NZO,P,FV1)
                        ICAV = 0
        ENDIF
C
        IF(MAXFEV.LE.2) THEN
                IF(MAXFEV.EQ.0) THEN
                        CALL MODEL(NZO,P,FV1)
                WRITE(3,1010)
C
                ELSEIF(MAXFEV.EQ.2) THEN             !SUBTRACTION
                   DO 129 IV = 1,NP
                        IF(IABS(NFREE(IV)).EQ.3) THEN
                                PS(IV) = P(IV)
                                PSX(IV) = 0.D0  
                           IF(IV.LE.NZO) THEN
                                NFRES(IV) = NFREE(IV)
                                NFREX(IV) = 0
                           ENDIF
                        ELSE
                        IF(IV.GT.30.AND.IV.LT.41) THEN
                                        PS(IV) = P(IV)
                                ELSE
                                        PS(IV) = 0.D0
                                ENDIF
                                PSX(IV) = P(IV)
                           IF(IV.LE.NZO) THEN
                                NFRES(IV) = 0
                                NFREX(IV) = NFREE(IV)
                           ENDIF
                         ENDIF
129                CONTINUE
                        DO KP = 1,NZO
                                NFREE(KP) = NFRES(KP)
                        ENDDO
                        CALL MODEL(NZO,PS,FV1)
                        DO KR = 1,MD
                                TEMP = Y(KR)
                                TEMP2 = Y(KR+MD)
                                Y(KR) = Y(KR) - FV1(KR)
                                Y(KR+MD) = Y(KR+MD) - FV1(KR+MD)
                        ENDDO
                ENDIF
C
         IF(IDMOD.GE.1) THEN
                MAXFEVN = 91
                IF(IKW.EQ.0) THEN
                        KRCH = IRCH
                ELSE
                        KRCH = NRCH
                ENDIF
C
                IF(MAXFEV.NE.2) CALL MODEL(NZO,P,FV1)
C
        IF(INDE.NE.0) GOTO 506
           IF(P(5).GT.0.D0.AND.MAXFEV.NE.2) THEN
                ETAU = P(6)/(P(5)*CELCAP)
                ECINF = ETAU/XXM1
                EC0 = ETAU*XX1
           ELSE
                ECINF = 0.D0
                EC0 = 0.D0
           ENDIF
C
506     CONTINUE
C
        IF(INDE.NE.0) THEN
                ECINF = 0.D0
                EC0 = 0.D0
        ENDIF
C
C       IF(FUN.EQ.OO.AND.MAXFEV.EQ.0.AND.NELEM0.EQ.7) THEN
C               CALL GEDAE(1,-1,FREQ,QX,F,0,2,0)
C
        IF(FUN.EQ.OO.AND.NELEM0.EQ.7.AND.P(35).EQ.4.D0) THEN
                IF(NELEM0.EQ.P(10).AND.P(9).GT.0.D0) THEN       
                        XXM1 = 1.D0/P(9)
                ENDIF
                XXM1 = XXM1 + SML
C       
                ETAU = P(6)/(P(5)*CELCAP)
                ECINF = ETAU/XXM1
CCC             EC0 = ETAU*XX1
                EC0 = 0.D0
C
        ENDIF
C
        IF(MAXFEV.LE.2) THEN                    ! #66
C
      WRITE(9,55) ALPHA
      WRITE(9,51) IOP,DFIT,DFIT,PFIT,PFIT,FREEQ,BLANK,FUN,CELCAP,DATTYP,
     1IPAR,ROE,IFPO,IRE,ECINF,EC0
C
      WRITE(9,69) MDS,NP,MAXFEVN,NPRINT,KRCH,MODE,ICP,IPRIN,IGACC,
     1ATEMP,SDWC,SDRC
C
        IF(MAXFEV.EQ.0.AND.IP40A.EQ.5.AND.MDE.EQ.-16) THEN
         CALL DSORT(P,NFREE,JIA,JIB)
          DO 199 IX = 1, NFREI
                        IF(IX.LE.14) THEN
                                IY = IX
                        ELSE
                                IY = IX + 6
                        ENDIF
                        P(2*IY-1) = STRN(IX)
199       ENDDO
        ENDIF
C
           IF(MAXFEV.LE.1) THEN
                WRITE(9,110) ( P(I), I = 1, NP )
                WRITE(9,120) ( NFREE(I), I = 1, NZO )
           ELSEIF(MAXFEV.EQ.2) THEN
                WRITE(9,110) ( PSX(I), I = 1, NP )
                WRITE(9,120) ( NFREX(I), I = 1, NZO )
           ENDIF
C
C       RESET FREQUENCY COLUMN TO Hz IF FREEQ=FF ORIGINALLY
C
        IF(FREEQ.EQ.FF) THEN
           DO II =1, MD
                FREQ(II) = TFREQ(II)
           END DO
        ENDIF
C
       IF(MAXFEV.EQ.0) THEN
                DO 3807 I = 1, MD
                        IM = I + MD
                IF(IP10T.EQ.1) FV1(IM) = 1.D3/FREQ(I)
C
                        WRITE(9,370) I, FREQ(I), FV1(I), FV1(IM)
3807            CONTINUE                
       ELSEIF(MAXFEV.EQ.1.OR.MAXFEV.EQ.2) THEN
                DO 3808 I = 1,MD
                        IM = I + MD
                        DP39 = DABS(P(39))
                IF(P(39).LT.0.D0.AND.DP39.LE.1.D17) THEN
                                Y(I) = DP39*Y(I)
                                Y(IM) = DP39*Y(IM)
          ELSEIF(DP39.GT.1.D17.AND.DFIT.EQ.YY.AND.FUN.EQ.'O') THEN
                                Y(I) = P(5)*Y(I)
                                Y(IM) =P(5)*Y(IM)
                ENDIF
C
                IF(IP10T.EQ.1) Y(IM) = 1.D3/FREQ(I)
C
           IF(FUN.NE.'R') THEN
                IF(P(38).GT.0.D0.AND.P(38).LE.1.D17) THEN 
                        FREQ(I) = P(38)*FREQ(I)
                ELSEIF(P(38).GT.1.D17.AND.FUN.EQ.'O') THEN 
                        FREQ(I) = P(6)*FREQ(I)
                ENDIF
           ENDIF
C
        IF(FUN.EQ.'O'.AND.P(36).NE.0.D0.AND.DFIT.EQ.YY) THEN
                Y(I) = Y(I) - 1.D0/(P(5))
                Y(IM) = -Y(I)/(TUPI*FREQ(I)*EPV)
CC
        ENDIF
C
                        WRITE(9,370) I, FREQ(I), Y(I), Y(IM)
CCC
3808    CONTINUE                
C
       ENDIF
C
C
CC
           IF(IRCH.EQ.0) THEN
                 DO 6533 I=1,MD
C          WRITE LINES 7 (if they exist)
                         WRITE(9,652) I,R(I),R(I+MD)
6533            CONTINUE
                IF(R(1).EQ.0) THEN
                        WRITE(*,736)
                        CALL STOPX
                ENDIF   
            ENDIF
        ENDIF                                   ! #66
CCC
        CLOSE(9)                ! CLOSE FOR MAXFEV <= 2
C
          ENDIF
C
       ENDIF
C
C       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
973         WRITE(4,1010)
1010    FORMAT(//'MODEL DATA -- FREQUENCY, REAL, IMAGINARY')
C
        DO 9893 KXW=1,MD
                KXM = KXW + MD
         IF(MAXFEV.LE.2) WRITE(3,1020) FREQ(KXW), FV1(KXW), FV1(KXM)
C
           IF(MAXFEV.EQ.0) THEN
                IF(DATTYP.EQ.'R') THEN
                        WRITE(4,1020) FREQ(KXW),FV1(KXW),0.D0
                ELSEIF(DATTYP.EQ.'I') THEN
                        WRITE(4,1020) FREQ(KXW),0.D0,FV1(KXM)
                ELSE
                        WRITE(4,1020) FREQ(KXW),FV1(KXW),FV1(KXM)
                ENDIF   
           ELSEIF(MAXFEV.EQ.1.OR.MAXFEV.EQ.2) THEN
                IF(DATTYP.EQ.'R') THEN
                        WRITE(4,1020) FREQ(KXW),Y(KXW),0.D0
                ELSEIF(DATTYP.EQ.'I') THEN
                        WRITE(4,1020) FREQ(KXW),0.D0,Y(KXM)
                ELSE
                        WRITE(4,1020) FREQ(KXW),Y(KXW),Y(KXM)
                ENDIF   
           ELSE
                IF(DATTYP.EQ.'R') THEN
                        WRITE(4,1020) FREQ(KXW),FV1(KXW)
                ELSEIF(DATTYP.EQ.'I') THEN
                        WRITE(4,1020) FREQ(KXW),FV1(KXM)
                ELSE
                        WRITE(4,1020) FREQ(KXW),FV1(KXW),FV1(KXM)
                ENDIF   
           ENDIF
C
          IF(IDMOD.EQ.2) THEN
                IF(IOPR.EQ.1) THEN        ! FF=-11 OR SMALLER
                   IF(FREEQ.EQ.' '.AND.IRE.LE.-11) THEN
                        FREQX = TFREQ(KXW)
                   ELSE
                        FREQX = FREQ(KXW)
                   ENDIF
                ELSE
                        FREQX = OFREQ(KXW)
                ENDIF
CC
                ALF = DLOG10(FREQX+SML)
                ALRI = DLOG10(DABS(Y(KXW)+SML))
                ALII = DLOG10(DABS(Y(KXM)+SML))
                ALRO = DLOG10(DABS(FV1(KXW)+SML))
                ALIO = DLOG10(DABS(FV1(KXM)+SML))
C
                WRITE(2,1020) FREQX,FV1(KXW),FV1(KXM)
        WRITE(17,1028) FREQX,Y(KXW),DABS(Y(KXM)),FV1(KXW),
     1                 DABS(FV1(KXM))
C
                IF(P(35).EQ.4.D0.AND.P(10).NE.33.D0) THEN
                        IF(P(5).LE.0.D0) P(5) = SML - P(5)
                        ALFN = ALF + DLOG10(P(6)+SML)
                        ALRON = ALRO + DLOG10(P(5)+SML)
                        ALION = ALIO + DLOG10(P(5)+SML)
                WRITE(21,1028) ALFN,ALRON,ALION
                ENDIF
C
                WRITE(18,1028) ALF,ALRI,ALII,ALRO,ALIO
                WRITE(22,1020) FREQX,FV1(KXW),FV1(KXM)
                WRITE(20,1029) FREQX,ALF,Y(KXW),ALRI,DABS(Y(KXM)),ALII,
     *          FV1(KXW),ALRO,DABS(FV1(KXM)),ALIO
            IF(ICOND.EQ.1) THEN
                CONR = FC*Y(KXW)
                CONI = FC*Y(KXM)
                FCNR = FC*FV1(KXW)
                FCNI = FC*FV1(KXM)
                WRITE(19,1028) FREQX,CONR,CONI,FCNR,FCNI
            ELSEIF(ICOND.EQ.2) THEN
                TIN = 1.D3/(FREQX + SML)
                RELRR = (FV1(KXW) - Y(KXW))/(Y(KXW) + SML)
                WRITE(19,1028) FREQ(KXW),TIN,ALRI,ALRO,RELRR
            ENDIF
          ENDIF
9893    CONTINUE
C
        IF(MAXFEV.LE.2) GOTO 317
CCC
C
1020    FORMAT(5X,3E25.13)
1028   FORMAT(1P,(1E12.4),1P,(1E20.12),1P,(1E14.6),1P,(1E20.12),1P,
     *       (1E14.6))
1029   FORMAT(1X,1P,(10E14.6))
C
        IF(INFO.EQ.5) WRITE(*,937)
937     FORMAT(5X, '!!! INCOMPLETE CONVERGENCE:  MAXFEV TOO SMALL !!!') 
          IF(IDMOD.EQ.2) THEN
                WRITE(2,3702) ONENEG
                CLOSE(2)
                CLOSE(22)
                CLOSE(17)
                CLOSE(18)
                CLOSE(19)
                CLOSE(20)
                CLOSE(21)
          ENDIF
C
      CONTINUE
C
C       **************************************************
C
        IF(NELEM.LT.7) GOTO 669
        IF(NELEM.EQ.7.AND.(P(7).LT.1.D0.AND.P(9).LT.1.D0)) GOTO 669
C
C       &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C       FINAL MOMENT AND LIMIT CALCULATIONS
C
        IF(NCH.GT.6) GOTO 317
C
C       COUPLING MODEL:
C       MUST SET TC IN P(41) BEFORE RUNNING!!!!
C
        IF(NCH.EQ.5.OR.NCH.EQ.6) THEN
                TC = P(41)
                TSIG = P(6)
                TSIGN = TSIG/TC
                BT = P(9)
                BTI = 1.D0/BT
                BT2 = 0.5D0*BT
                BTM = 1.D0 - BT
                TSNB= TSIGN**BT
                TAOFCM = TC*TSNB
                IF(ABS(MDE).NE.8) THEN
                        WC = 1.D0/TC
                ELSE
                        WC = TC                 !TRANSIENT RESPONSE
                        GAMX = P(45)                    
                        WRITE(*,*) 123,TAOFCM,P(45)
                        GOTO 317
                ENDIF
                AAO = P(42)
                AAON = AAO/P(5)
                TAOO = P(43)
                EGAM = 0.42D0
                QZ= BTI*DEXP(-BTM*EGAM)
                QZB = QZ**(-BTI)
C
                YS = WC*TSIG
                YSI = 1.D0/YS
                YSQ1 = 1.D0 + YS*YS
                YSP = YSQ1**(-BT2)
                THT = BT*DATAN(YS)
                YSPY = YSP*YSI
C
                IF(NELEM.EQ.10.OR.NELEM.EQ.32) THEN
                        AGNL = GAMMLN(BTI)
                        XXB = BTI*DEXP(AGNL)
                        TONC = QZ*(TSIGN**BT)
C
                ELSEIF(NELEM.EQ.6) THEN
                        XXB = BT
                ENDIF
C
                XX1I = 1.D0/XXB
                DSN = DSIN(THT)
                DCS = DCOS(THT)
                        GAM1C = XX1I*YSPY*DSN 
                        GAM2C = XX1I*(YSI - YSPY*DCS) 
C
                DO IJ = 1,4
                        QX(IJ) = FV1(IJ)
                ENDDO
                        FR1 = FREQ(1)
                        FR2 = FREQ(2)
C
                FREQ(1) = 0.9999999D0*WC
                FREQ(2) = 1.00000001D0*WC
                MD = 2
C
                CALL MODEL(NZO,P,FV1)
C
                GAM1 = FV1(1)
                GAM2 = -FV1(3)          
                GAM1P = FV1(2)
                GAM2P = -FV1(4)
                WRITE(*,*)
                WRITE(*,*) FREQ(1),GAM1,GAM2
                WRITE(*,*) FREQ(1),GAM1C,GAM2C
                WRITE(*,*) FREQ(2),GAM1P,GAM2P
                WRITE(*,*) 100,NELEM,GAM1C/GAM1,GAM2C/GAM2
                WRITE(*,*)
C
                TRANS = DCMPLX(FV1(1),-FV1(3))
                TRANS = 1.D0/TRANS
                AAC = 1.D0/DREAL(TRANS)
                TAOEPS = -TC*AAC*DIMAG(TRANS)
                TAOC1 =  TAOEPS - AAC*P(2)*CELCAP
CCC     WRITE(*,*) AAC,TAOEPS,TAOC1 
CC      READ(*,*)
                VV = GAM2/GAM1
                VSQ1 = 1.D0 + VV*VV
                TAOC1N = VV
                TON = VV                        !FROM CONTINUITY
                TSNC = QZB*(TON**BTI)
                TOR = TONC/TON
                TSR = TSNC/TSIGN
                WRITE(*,212) 333,BT,TOR,TSR
C
                TCM = TC**(-BTM)
                AAC = GAM1*VSQ1
                AACN = AAC/P(5)
                TSIGC = (TCM*TAOC1)**(1.D0/BT)
                TAOC = (TSIG**BT)*(1.D0/TCM)
                TSIGCSN = (TAOC1N*DSN)**(1.D0/BT)
                TAOCSN = ((TSIGN)**BT)/DSN
                UU = DLOG(TSIGN)
                QJ = VV/(TSIGN**BT)
                GAMX = -(1.D0/BTM)*DLOG(BT*QJ)
                P(44) = QJ
                TAOCK = TC*BTI*DEXP(-GAMX*BTM)*(TSIGN**BT)
                WRITE(*,212) 111,GAMX,TAOCK,TAOFCM
C
                TSNL = DLOG10(TSIGN)
                AN = AACN*(TSIGN**BTM)
                GAMXX = GAMX*1.732464D0
                WRITE(73,211) TSNL,AN,GAMXX,BT
C       
                WRITE(*,211) BT,TC,TAOO,XXB
                WRITE(*,*) AACN,TSIG,TAOC1
                WRITE(*,*) AAON,TSIGC,TAOC
                WRITE(*,*) 200,MDE,TSIGC/TSIG, TAOC/TAOC1
                WRITE(*,*) 300,MDE,TSIGCSN/TSIGN, TAOCSN/TAOC1N
                WRITE(*,*) 400,TSIG,TC,UU
C
CC              OPEN(9,FILE='OUTIN')
                P(42) = AAC
                P(43) = TAOC1
                DO IJ = 1,4
                        FV1(IJ) = QX(IJ)
                ENDDO
                        FREQ(1) = FR1
                        FREQ(2) = FR2
C
        IF(IRE.LT.0) THEN
          IF(IACY.NE.0) IOP = IOPO
        ENDIF
        ENDIF
C

C       *********************************************************
C
        OPEN(9,FILE='OUTIN')
        IF(INDE.GE.1) GOTO 669
C
        IF(FUN.EQ.OO.AND.P(5).EQ.0.D0.AND.P(6).EQ.0.D0.AND.P(40).EQ.
     *  0.D0) GOTO 317
C
                MDEX = IABS(MDE)
        IF(MDEX.LE.3.OR.MDEX.EQ.6) THEN 
                CONTINUE
C --------------- BY AG ---------------
	  ELSEIF (MDEX.EQ.8  .AND.  INDE.EQ.0) THEN
			  GOTO 669
C --------------- BY AG ---------------
        ELSE
                GOTO 317
        ENDIF
C
                ICAV = 1
                INL = 0
                NCHK = NCH
C
        IF(FUN.EQ.OO.AND.P(10).GT.0.D0) NELEM = P(10)
CCC
        IF(NELEM.EQ.7.OR.NELEM.EQ.10.OR.NELEM.EQ.32.OR.NELEM.EQ.36.OR.
     *  MDE.EQ.-6) INL = 1
        IF(NELEM.NE.NELEM0) INL = 1
C
        IF(NELEM.EQ.7) GOTO 587
C
        IF(NELEM.EQ.10.OR.NELEM.EQ.32.OR.NELEM.EQ.36) GOTO 609
        IF(FUN.EQ.OO) THEN
                IF(NCH.NE.4) GOTO 609
CCC
                IF(NELEM.EQ.0.D0.OR.INL.EQ.0) GOTO 317
        ENDIF
        IF(INL.EQ.0) GOTO 317
C
587             IF(MDE.GE.0.D0) THEN 
                        NCH = 2
                        NCHK = 3
                        IF(FUN.EQ.OO) RX = P(5)
                        QX(1) = RX              !  HERE RX IS RHOO
                        RHOO = RX
                ELSE
                        NCH = 3
                ENDIF   
C
                IF(FUN.EQ.OO) THEN
                TAO = P(6)
                PHI = P(9)
                IF(NCHK.EQ.3.OR.NCHK.EQ.2) THEN
CCC                     QX(7) = P(3)            ! GAMMA
CCC                     QX(3) = P(7)            ! UU UPPER LIMIT
                        QX(7) = P(7)            ! GAMMA
                        QX(3) = P(3)            ! UU UPPER LIMIT

                ELSE
                        QX(7) = P(7)            ! GAMMA
                        QX(3) = 1.5D1           ! INTERNAL UU LIMIT
                ENDIF
        ELSE
                TAO = TX
                QX(7) = UU                      !  HERE GAMMA
                QX(3) = 1.5D1                   ! INTERNAL UU LIMIT
        ENDIF
                GAMM = QX(7)
                QX(2) = TAO                     !  HERE T
                QX(5) = PHI
C
CCC     WRITE(*,*) '    ******   PLEASE WAIT FOR FINAL OUTPUT   ******'
C
                IF(GAMM.EQ.1.D0) THEN
                        NCH = 2
                ELSE
                        NCH = 3
                ENDIF
C
        IF(P(10).EQ.7.D0.AND.(P(7).EQ.1.D0.AND.P(9).EQ.1.D0).OR.P(7).
     *  LT.0.08D0) THEN
                JCDX = 0
                GOTO 469
        ENDIF
C
        IF(NELEM0.EQ.7.AND.P(7).EQ.1.D0.AND.MDE.LT.0) THEN
                NCH = 2
        WRITE(*,*) '******  CD RESPONSE FROM DRT INTEGRATION  ******'
CCC     ENDIF
                CALL GEDAE(1,-1,FREQ,QX,F,0,NCH,0)
        ENDIF
CC
                XXM1 = XXM1 + SML
C
                IF(NCH.EQ.3) THEN
                        IF(MDE.LT.-1) THEN
                                GOTO 317
                        ELSEIF(MDE.EQ.-1) THEN
                                IF(GAMM.NE.1.D0) THEN
                                        RHOO = TX/(CELCAP*XXM1*RX)
                                ELSE            
                                        RHOO = RX
                                ENDIF
                        ENDIF                   
                ELSE
                                RHOO = RX
                ENDIF
                GOTO 613
C
609             NCHK = NCH
C
        IF(FUN.EQ.OO) THEN
           IF(P(40).EQ.0.D0) THEN
                RHOO = P(5)
                TAO = P(6)
           ELSE                 
                RHOO = P(1)
                TAO = P(16)
           ENDIF                
        IF(P(10).EQ.3.6D1) FREQIM = DISTEL(P(5),P(6),P(7),P(9),P(10),
     *  0.D0)
        ELSE
                RHOO = RX
                TAO = TX
        ENDIF
C
613     CONTINUE
C
        IF(FUN.EQ.OO.AND.MDE.EQ.-6) JCDX = 0
           IF(FUN.EQ.OO.AND.NCHK.NE.4) THEN
                IF(NCH.NE.3) THEN
                        RHOO = P(5)
                ELSEIF(NCH.EQ.3.AND.MDE.EQ.-1) THEN
                        XXM1 = XXM1 + SML
                        RHOO = TAO/(CELCAP*XXM1*P(5))
                ENDIF
           ENDIF
C
469        IF(FUN.EQ.OO.AND.JCDX.EQ.0) THEN
                IF(NCH.NE.3.OR.MDE.EQ.-6) THEN
                        RTOT = P(1) + RHOO                              
                        REFI = RHOO/(RTOT*RTOT)
                ENDIF
           ELSE
                RTOT = RHOO
                REFI = 1.D0/RHOO
           ENDIF
C
          IF(MDEX.LE.1.OR.MDE.EQ.-6) THEN
                IF(RHOO.NE.0) THEN
                        ETAU = TAO*REFI/CELCAP
                ELSE
                        ETAU = 0.D0
                ENDIF
          ELSEIF(MDEX.EQ.2) THEN
                ETAU = TAO                      !HERE TAO IS ETAU
                TAO = CELCAP*ETAU/REFI
          ELSE
                TAO = P(6)
                ETAU = 0.D0
          ENDIF
C
        IF(P(7).EQ.1.D0.AND.P(9).EQ.1.D0) THEN
                WRITE(*,*) 'WRONG INPUT: NO DEBYE HERE'
                        GOTO 317
        ENDIF
                ECINF = SML + ETAU/XXM1
                IF(FUN.EQ.OO) THEN
                        IF(ATEMP.LT.0) THEN     
                                EDINF = P(2)
                        ELSE
                                EDINF = P(2)/CELCAP
                        ENDIF
                ELSE
                        IF(P(29).NE.0.D0) THEN
                                WRITE(*,239)    
239   FORMAT(7X,'????  CORRECT EDINF CALC NEEDS CELCAP = EPSVACUUM ???')
                                WRITE(*,*)
C
                        ENDIF
                        EDINF = P(29)/CELCAP
                ENDIF
                EC0 = ETAU*XX1
                XRAT =  EC0/ECINF
CCC                     E0 = FV1(2)
CCC                     E0M = E0 - EDINF
CCC     IF(E0.GT.4.D4.OR.E0M.LT.1.D-3.OR.NELEM.EQ.7.OR.NELEM.EQ.36) 
CCC     *                  E0 = EDINF + EC0
                E0 = EDINF + EC0
C
CCC                     EINF = FV1(3)
CCC             IF(NELEM.EQ.7.OR.NELEM.EQ.36) EINF = EDINF + ECINF
                IF(INDE.EQ.0) EINF = EDINF + ECINF
C
CC      WRITE(*,*) MDE,E0,EC0
C
                AVTAU = TAO*XX1
                IF(MDE.EQ.-6) THEN 
                        AIN = AIN/XXM1
                        RN = RN/XXM1
                ENDIF
C               
        IF(P(40).EQ.4.D0.AND.(NELEM.EQ.10.OR.NELEM.EQ.36)) THEN
                EC0 = 0.D0
                ECINF = 0.D0
                EINF = EDINF
                E0 = EDINF + P(15)
          ELSEIF(P(35).EQ.4.D0.AND.NELEM.EQ.10.AND.MDE.EQ.-3) THEN
                        XXM1 = XXM1 + SML
                        RHOO = TAO/(CELCAP*XXM1*P(5))
                        ECINF = P(5)    
                        RTOT = P(1) + RHOO                              
                        REFI = RHOO/(RTOT*RTOT)
                        ETAU = TAO*REFI/CELCAP
                        EC0 = ETAU*XX1
                        E0 = EDINF + EC0
        ENDIF
C
        WRITE(*,213)
        WRITE(3,213)
213     FORMAT(8X,'XXM1           XX1            XX2          XX3')
        WRITE(*,214)
        WRITE(3,214)
214     FORMAT(8X,'AVTAU          XRAT           RN           AIN')
C
        WRITE(*,211) XXM1,XX1,XX2,XX3
        WRITE(*,211) AVTAU,XRAT,RN,AIN
        WRITE(3,211) XXM1,XX1,XX2,XX3
        WRITE(3,211) AVTAU,XRAT,RN,AIN
C
        WRITE(*,215)
        WRITE(3,215)
215     FORMAT(/,8X,'RHOO           EDINF          EINF         E0')
        WRITE(*,216)
        WRITE(3,216)
216     FORMAT(8X,'TAO            ETAU           ECINF        EC0')
        WRITE(*,211) RHOO,EDINF,EINF,E0
        WRITE(*,211) TAO,ETAU,ECINF,EC0
        WRITE(3,211) RHOO,EDINF,EINF,E0
        WRITE(3,211) TAO,ETAU,ECINF,EC0
C
        IF(P(1).NE.0.D0) THEN
         WRITE(*,*) '     *****  SINCE P(1) NOT ZERO, FINAL EINF= EDINF
     #  ******'
         WRITE(3,*) '     *****  SINCE P(1) NOT ZERO, FINAL EINF= EDINF
     #  ******'
        ENDIF
CC
669     CONTINUE
        IF(MAXFEV.LE.2) GOTO 317
C
           IF(NELEM.EQ.10.OR.NELEM.EQ.32.OR.NELEM.EQ.36.OR.NELEM.EQ.7)
     *     THEN
                CONTINUE
            ELSEIF((NELEM.EQ.7.AND.P(6).NE.1.D0).OR.P(40).EQ.4.D0) THEN
                ECINF = 0.D0
                EC0 = 0.D0
            ENDIF
C
      IF(IDMOD.GE.1) THEN
                WRITE(9,55) ALPHA
C
      WRITE(9,51) IOP,DFIT,DFIT,PFIT,PFIT,FREEQ,BLANK,FUN,CELCAP,DATTYP,
     1IPAR,ROE,IFPO,IRE,ECINF,EC0
C
      WRITE(9,69) MDS,NP,MAXFEV,NPRINT,NRCH,MODE,ICP,
     1IPRIN,IGACC,ATEMP,SDWC,SDRC
C
        IF(DABS(P(40)).GE.1.D0.AND.(FUN.EQ.KZ.OR.FUN.EQ.RR)) CALL
     *  DSORT(P,NFREE,JIA,JIB)
C
                WRITE(9,110) ( P(I), I = 1, NP )
                WRITE(9,120) ( NFREE(I), I = 1, NTOT )
                DO 678 I = 1, MDA
                        IM = I + MDA
C
                IF(IP10T.EQ.1) Y(IM) = 1.D3/FREQ(I)
C
                        WRITE(9,370) I, FREQ(I), Y(I), Y(IM)
678             CONTINUE
      ENDIF
C
        IF(IRCH.EQ.0) THEN
                 DO 6531 I=1,MD
                         WRITE(9,652) I,R(I),R(I+MD)
6531            CONTINUE
                IF(R(1).EQ.0) THEN
                        WRITE(*,736)
                        CALL STOPX
                ENDIF   
        ENDIF
C
211     FORMAT(1X,1P,(4E15.5))
212     FORMAT(1X,I5,1P,(4E15.6))
C
317   CLOSE(1)
      CLOSE(3)
      CLOSE(4)
      CLOSE(9)  
C     END OF MAIN ROUTINE
      CALL STOPX
      END
C
C       PROGRAM TO SORT PARAMETER DATA
C           J.R. Macdonald
C           
        SUBROUTINE DSORT(P,NFREE,JIA,JIB)
        IMPLICIT REAL*8 (A-H,O-Z)
        COMPLEX*16 ZT,YCC(20)
        EXTERNAL DXSPS
      DIMENSION P(*),PR(20),PC(20),TAU(20),ISP(20),RTSD(20),TM(6),
     1NFREE(*),IORIG(20),STR(20),PSX(40),STRS(20),TAUS(20)
        COMMON /CM12/ CELCAP,ATEMP,WF,MAXFEV,ICF,MDE,JCDX
        COMMON /CM20/ RXSD(40)
        COMMON /CM53/ IRG,IPRIN,STRN(20),AI
C
                OPEN(72,FILE='OUTPARM')
                OPEN(43,FILE='OUTDRPT')
C
        SCALE = 0.D0
        IP40A = INT(DABS(P(40)) + 0.01D0)
        IP38A = INT(DABS(P(38)) + 0.01D0)
        IP38 = SIGN(1.D0,P(38))*IP38A
        SUME = 0.D0
C
C               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                NPH = ICF/2
                ICR = 0
C
                OPEN(41,FILE='OUTDRT')
C
C       SEPARATE EVEN AND ODD PARAMETER VALUES
C
        JJ = 0
        IA = 0
        DO IK = 1,ICF
                IF(IK.GT.JIA) THEN
                        KJ = IK + JIB
                ELSE
                        KJ = IK
                ENDIF
C
                IF(NFREE(KJ).EQ.0) THEN
                        IA = IA + 1
                        PSX(IK) = 0.D0
                ELSE
                        PSX(IK) = RXSD(IK-IA)
                ENDIF
        ENDDO
C
        DO KK = 1,ICF-1
                IF(KK.GT.JIA) THEN
                        KQ = KK + JIB
                ELSE
                        KQ = KK
                ENDIF
                JK = KQ + 1
                IF(MOD(JK,2).EQ.0) THEN
                        JJ = JJ + 1
C       P(KQ) IS R OR C(DIEL)
                        PR(JJ) = P(KQ)
                        PC(JJ) = P(JK)
                        RXSD(JJ) = PSX(KK)
                        RTSD(JJ) = PSX(KK+1)
C
                   IF(IP40A.EQ.1) THEN
C       PC IS CAPACITANCE (OR R DIEL HERE)
                        TAU(JJ) = DABS(PR(JJ)*PC(JJ))
                        STR(JJ) = PR(JJ)
                   ELSEIF(IP40A.EQ.2.OR.IP40A.EQ.3.OR.IP40A.EQ.10) THEN
C       Here, PC is initially TAU
                        TAU(JJ) = PC(JJ)
                        STR(JJ) = PR(JJ)
                        PC(JJ) = TAU(JJ)/PR(JJ)
                   ELSEIF(IP40A.EQ.4) THEN
                        TAU(JJ) = PC(JJ)
C                       STR(JJ) = PR(JJ)*PC(JJ)
                        STR(JJ) = PR(JJ)
                        PC(JJ) = TAU(JJ)/PR(JJ)
                   ELSEIF(IP40A.EQ.5) THEN
                        TAU(JJ) = PC(JJ)
                        STR(JJ) = PR(JJ)
                   ENDIF
                        SCALE = SCALE + STR(JJ)
                        SUME = SUME + 1.D0/PC(JJ)
                ENDIF
        ENDDO
C
209     FORMAT(2I4,1X,1P,(4E15.6))
C
                IF(MDE.EQ.8) THEN       
                        NPN = NPH - 1
                ELSE
                        NPN = NPH
                ENDIF
C
        CALL SORT(TAU,ISP,NPH)
        DO KM = 1,NPH
                TAUS(KM) = TAU(ISP(KM))
                STRS(KM) = STR(ISP(KM))
        ENDDO
C
        IF(IP38A.EQ.0) THEN
                AR = 1.D0
                SCLI = 1.D0/SCALE
                AI = SCLI
        ELSEIF(IP38A.EQ.1) THEN
                SCALE = SCALE*WF
                SCLI = 1.D0/SCALE
                AR = WF
                AI = SCLI
        ELSEIF(IP38A.GT.1.AND.IP38.LT.5) THEN
           DO KJ = 1,6
C               GET NORMALIZATION VALUE
             IF(KJ.EQ.1) THEN
                DO IJ = 1, NPH
                        YCC(IJ) = DCMPLX(STRS(IJ),0.D0)
                ENDDO
C               GET AV TAU
             ELSEIF(KJ.EQ.2) THEN       
                DO IJ = 1, NPH
                   YCC(IJ) = TAUS(IJ)*DCMPLX(STRS(IJ),0.D0)
                ENDDO
C               GET AV TAU*TAU
             ELSEIF(KJ.EQ.3) THEN       
                DO IJ = 1, NPH
                   YCC(IJ) = (TAUS(IJ)**2)*DCMPLX(STRS(IJ),0.D0)
                ENDDO
C               GET AV 1/TAU
             ELSEIF(KJ.EQ.4) THEN       
                DO IJ = 1, NPH
                   YCC(IJ) = DCMPLX(STRS(IJ),0.D0)/TAUS(IJ)
                ENDDO
C               GET AV TAU*TAU*TAU
             ELSEIF(KJ.EQ.5) THEN       
                DO IJ = 1, NPH
                   YCC(IJ) = DCMPLX(STRS(IJ),0.D0)*(TAUS(IJ)**3)
                ENDDO
C               GET VALUE OF f(t) at t=P(36)
             ELSEIF(KJ.EQ.6) THEN       
                DO IJ = 1, NPH
                   YCC(IJ) = DEXP(-P(36)/TAUS(IJ))*DCMPLX(STRS(IJ),0.D0)
                ENDDO
             ENDIF
C
C C               SUBROUTINE DXSPS IS IN RSUB
C
                CALL DXSPS(TAU,YCC,ZT,WF,NPN,IP38A)
C
                TM(KJ) = DREAL(ZT)              
                IF(TM(KJ).EQ.0.D0) VTM = VTM + TM(KJ)
C
           ENDDO
C
                SCALEN = TM(1)
                AI = 1.D0/SCALEN
                SCALE = SCALEN
C
                WRITE(72,*)
C
        ELSE
                WRITE(*,824)
824             FORMAT(3X,'WRONG VALUE OF P(38)')
                CALL STOPX
        ENDIF
632     CONTINUE
C
        SUMS = 0.D0             !
        SUMRR = 0.D0            !
        DO IK = 1,NPH
                IORIG(IK) = 2*ISP(IK) - 1
                IF(IORIG(IK).GT.JIA) IORIG(IK) = IORIG(IK) + JIB
C
                STS = AI*STRS(IK)
                STRN(IK) = STS
                SUMS = SUMS +STS                !
                PCS = PC(ISP(IK))
                PRR = TAUS(IK)*PR(ISP(IK))
                SUMRR = SUMRR + PRR             !
                PCC = 1.D0/PR(ISP(IK))
C
        IF(IP40A.NE.4) THEN 
      WRITE(72,183) IORIG(IK),TAUS(IK),RTSD(ISP(IK)),STS,
     *RXSD(ISP(IK)),PR(ISP(IK)),PCS
        ELSE 
      WRITE(72,183) IORIG(IK),TAUS(IK),RTSD(ISP(IK)),STS,
     *RXSD(ISP(IK)),PR(ISP(IK)),PCS,PRR,PCC
        ENDIF
       WRITE(*,183) IORIG(IK),TAUS(IK),RTSD(ISP(IK)),STS,
     *RXSD(ISP(IK)),PR(ISP(IK)),PCS
C
                ALT = DLOG10(TAUS(IK))
                ALD = DLOG10(DABS(STS))
                IF(RTSD(ISP(IK)).LE.0.D0) THEN
                        ALRTSD = 0.D0
                ELSE
                        ALRTSD = DLOG10(RTSD(ISP(IK)))
                ENDIF
                IF(RXSD(ISP(IK)).LE.0.D0) THEN
                        ALRXSD = 0.D0
                ELSE
                        ALRXSD = DLOG10(RXSD(ISP(IK)))
                ENDIF
C
                IF(IPRIN.GE.0) THEN
                  WRITE(41,180) TAUS(IK),ALT,STS,ALD,RXSD(ISP(IK)),
     1            RTSD(ISP(IK)),ALRXSD,ALRTSD
                ELSE
                        ALR = DLOG10(DABS(PR(ISP(IK))))
                        ALI = DLOG10(DABS(PCS))
                        WRITE(41,180) TAUS(IK),ALT,ALD,ALR,ALI
                ENDIF
                WRITE(43,179) TAUS(IK),STS,-1.E0
C
C               *************************************
           IF(IRG.GT.0) THEN
C               OPEN DRTLSD FOR VARIOUS TYPES OF RECT. PLOTTING
C
                OPEN(42,FILE="DRTLSD.DAT")
C
                PMN = 1.D0 - RXSD(ISP(IK))
                PLM = 1.D0 - RTSD(ISP(IK))
                IF(PMN.LE.0.D0) PMN = 1.D-3
                IF(PLM.LE.0.D0) PLM = 1.D-3
                ALTM = TAUS(IK)*PLM
                DMN = STS*PMN
C
                ALTP = TAUS(IK)*(1.D0 + RTSD(ISP(IK)))
                DMP = STS*(1.D0 + RXSD(ISP(IK)))
             IF(IRG.EQ.2) THEN
                ALTM = DLOG10(ALTM)
                ALTP = DLOG10(ALTP)
                DMN = DLOG10(DABS(DMN))
                DMP = DLOG10(DABS(DMP))
             ELSEIF(IRG.EQ.3) THEN
                ALTM = DLOG10(ALTM)
                ALTP = DLOG10(ALTP)
             ENDIF
                WRITE(42,180) ALTM,DMN
                WRITE(42,180) ALTM,DMP
                WRITE(42,180) ALTP,DMP
                WRITE(42,180) ALTP,DMN
                WRITE(42,180) ALTM,DMN
        WRITE(42,180) 1.D5,1.D5
           ENDIF
179     FORMAT(5X,3(5X,E20.13))
180     FORMAT(1P,(8E12.4))
        ENDDO
C               ******************************
C                
        EPC = 1.D0/SUME
                IF(IP40A.EQ.4) THEN
                        IF(IP38A.EQ.2) SUMRR = TM(2)
                        WRITE(72,193) SUMS,SUMRR,EPC,SCALE              !
                        WRITE(*,193) SUMS,SUMRR,EPC,SCALE
                ELSE
                        WRITE(72,193) SUMS,SCALE,EPC
                        WRITE(*,193) SUMS,SCALE,EPC,AI
                ENDIF
                IF(IP38A.GT.1) THEN
                WRITE(72,193) AI,TM(1),TM(2),TM(3),TM(4),TM(5),TM(6)
                        TMI = TM(4)*AI
                        TM2 = TM(3)*AI
                        TM1 = TM(2)*AI
                        TM3 = TM(5)*AI
                        TMT = TM(6)*AI
                   IF(IP40A.EQ.4) THEN
                        TMI = 1.D0/TM1
                        TM1 = TMI*TM2
                        TM2 = TMI*TM3
                   ENDIF
                        WRITE(*,193) TMI,TM1,TM2,TMT
                        WRITE(72,193) TMI,TM1,TM2,TMT
                        ATMA = DABS(ATEMP)
                  IF(ATMA.LT.1.D80.AND.ATMA.NE.0.D0) THEN
                        ATINV = 1.D0/ATMA
C
                        TMI = TMI*ATMA                  !MOMENT X**-1
                        TM2 = TM2*ATINV*ATINV           !MOMENT X**2
                        TM1 = TM1*ATINV                 !AV X**1
                        WRITE(*,193) TMI,TM1,TM2,ATMA
                        WRITE(72,193) TMI,TM1,TM2,ATMA
                  ENDIF
                ENDIF
C
183     FORMAT(I3,1X,1P,(1E12.6),1P,(7E12.4))
193     FORMAT(2X,1P,(6E15.6))
194     FORMAT(/,2X,1P,(3E12.3),1P,(1E16.7))
C
                CLOSE(41)
                CLOSE(42)
                CLOSE(43)
                CLOSE(72)
C
        RETURN
        END
C
        SUBROUTINE SORT(TAU,ISP,NN)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION TAU(20),ISP(20)
C 
        DO KK = 1,NN
                ISP(KK) = KK
        ENDDO
        DO 20 I = 1, NN-1
                ASM = TAU(ISP(I))
                IND = I
                DO 15 J = I + 1,NN
                        IF(TAU(ISP(J)).LT.ASM) THEN
                                ASM = TAU(ISP(J))
                                IND = J
                        ENDIF
15              CONTINUE
                IF(IND.NE.I) THEN
                        ITMP = ISP(I)
                        ISP(I) = ISP(IND)
                        ISP(IND) = ITMP
                ENDIF
20      CONTINUE
        RETURN
        END
C
        SUBROUTINE RESORT(P,NFREE,TAUM,WF,X,NTOT,NS,JIA,JIB)
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION ISQ(20),P(*),NFREE(*),NFRE(20),NFRO(20),
     *  PE(20),PO(20),NFRET(40),PET(40),X(*),NS(*)
C
        IP34A = INT(DABS(P(34)) + 0.01D0)
        IP34 = SIGN(1.D0,P(34))*IP34A
        IP40A = INT(DABS(P(40)) + 0.01D0)
        IP40S = SIGN(1.D0,P(40))
C
        JX = 0
        ICO = 0
        DO KK = 1, NTOT
          IF((KK.LE.JIA.OR.KK.GT.40).AND.P(KK).NE.0.D0) THEN
                ICO = ICO + 1
                PET(ICO) = P(KK)
                NFRET(ICO) = NFREE(KK)
          ENDIF
        ENDDO
C
        IF(ICO.EQ.0.AND.P(35).GT.0) ICO = 2*P(35)
C
        IF(IP34.EQ.1) THEN
                NFRX = 1
                NFRT = 0
                IFX = 1
        ELSEIF(IP34.EQ.2) THEN
                NFRX = 2
                NFRT = 0
                IFX = 1 
        ELSEIF(IP34.EQ.3) THEN
                NFRX = 1
                NFRT = 1
                IFX = 0
        ELSEIF(IP34.EQ.4) THEN
                NFRX = 2
                NFRT = 1
                IFX = 0
        ELSEIF(IP34.EQ.5) THEN
                NFRX = 1
                NFRT = 0
                IFX = 2
        ELSE
                IFX = 3
        ENDIF
C
        IP35A = INT(DABS(P(35)) + 0.01D0)
        IP35 = SIGN(1.D0,P(35))*IP35A
C       
C       SEPARATE ORIGINAL VALUES INTO EVEN/ODD PARTS
        INFE = 0
        DO KK = 1, ICO
                KP = KK + 1
                IF(MOD(KP,2).EQ.0) THEN
                   JX = JX + 1
C               ONLY EVEN P TO PE 
                   PE(JX) = PET(KP)
                   NFRE(JX) = NFRET(KP)
                   IF(NFRE(JX).EQ.0) INFE = 1
C               ONLY ODD P TO PO
                   PO(JX) = PET(KK)
                   NFRO(JX) = NFRET(KK)
                ENDIF
        ENDDO
393     FORMAT(1X,2I4,1P,(2E15.6),2I4)
C
        IF(P(35).GT.0) THEN
                JKX = IP35A
        ELSE
                JKX = JX
        ENDIF
C
        IF(IFX.EQ.1.OR.IFX.EQ.2) THEN
                DO JJ = 1, JKX
                        PET(JJ) = TAUM
                        TAUM = TAUM*DEXP(WF)
                        IF(IFX.EQ.2) ISQ(JJ) = JJ
                ENDDO
        ELSE
                CALL SORT(PE,ISQ,JX)
        ENDIF
C
        DO KK = 1, JKX
                KK2 = 2*KK
           IF(KK2.LE.JIA) THEN
                KE = KK2
           ELSE
                KE = KK2 + JIB
           ENDIF
                KO = KE - 1
C
              IF(IFX.EQ.1) THEN
                NFREE(KO) = NFRX
                NFREE(KE) = NFRT
                P(KO) = 1.D0
                P(KE) = PET(KK)
              ELSEIF(IFX.EQ.0) THEN
                NFREE(KO) = NFRX
                NFREE(KE) = NFRT
                P(KO) = PO(ISQ(KK))
                P(KE) = PE(ISQ(KK))
              ELSEIF(IFX.EQ.2) THEN
                NFREE(KO) = NFRX
                NFREE(KE) = NFRT
                P(KO) = PO(ISQ(KK))
                P(KE) = PET(KK)
              ELSE
                NFREE(KE) = NFRE(ISQ(KK))
                NFREE(KO) = NFRO(ISQ(KK))
C
                IF(IP34.EQ.6) THEN
                        IF(IP40A.LT.4) THEN
                                P(KE) = PE(ISQ(KK))/PO(ISQ(KK))
                        ELSE
                                P(KO) = PE(ISQ(KK))*PO(ISQ(KK))
                                P(KE) = PE(ISQ(KK))/P(KO)
                        ENDIF
                ELSEIF(IP34.EQ.7) THEN
                        P(KE) = PE(ISQ(KK))*PO(ISQ(KK))
                ELSE
                        P(KE) = PE(ISQ(KK))
                ENDIF
                        IF(IP40A.NE.4) P(KO) = PO(ISQ(KK))
              ENDIF
C
        ENDDO   
C
        IF(IP34A.EQ.6) P(40) = IP40S
C
392     FORMAT(1X,2I4,1P,(2E15.6),2I4)
348     FORMAT(1X,I4,1P,(2E15.5),3I4)
C
        IF(IFX.LT.3) THEN
          IPD = 2*JKX
          IF(IPD.GT.JIA) IPD = IPD + JIB
          DO KK = IPD,NTOT
          IF(KK.LE.JIA.OR.KK.GT.40) THEN
                IF(KK.GT.IPD) THEN
                        P(KK) = 0
                        NFREE(KK) = 0
                ENDIF
          ENDIF
          ENDDO
        ENDIF
c
      DO 628 I = 1, NTOT
        IF (NFREE(I) .EQ. 0 ) GO TO 628
        LLL = LLL + 1
        NS(LLL) = I
        X(LLL) = P(I)
628   CONTINUE
C
        IF(INFE.EQ.1) P(35) = JKX
        IF(P(34).GT.0) P(34) = - 1.D0
C
        RETURN
        END

C       ********************************************
C       *
C       *       Added by Andriy Gorkovenko.
C       *
C       ********************************************

        SUBROUTINE STOPX
        CHARACTER*80 BUF
        INTEGER*2 RES
C       Get command line argument
        CALL GETARG(1, BUF, RES)
C       If argument is 'S' then pause the program
        IF ( BUF .EQ. 'S' ) THEN
        PAUSE 'PRESS ENTER TO CLOSE CONSOLE...'
        ENDIF
        STOP
        RETURN
        END
