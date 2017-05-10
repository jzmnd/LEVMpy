C
C       MODIFIED JEREMY SMITH 3/31/2017
C
      SUBROUTINE REGMAX(MP,IOR,ING,IFLAG,NFREI,MPN,MK,XR,FVEC,
     + FJAC,P33A)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 DATTYP
      INTEGER MP,IOR,ING,IFLAG,NFREI,MPN,MK
      REAL*8 XR,FVEC,FJAC,P33A
      INCLUDE 'SIZE.INC'
      DIMENSION IREG(40,40),MREG(40,40),
     + XR(*),XREG(40),FVEC(MPN),FJAC(MPN,NFREI)
      COMMON /CM2/ Y(NPT2),R(NPT2),FJ(NPT2),P(NTOT),DRSS,ROE,RKE,
     + NS(NPAFR),NFREE(NTOT),N,ICNT,MN,IRCH,IXI,DATTYP
      COMMON /CM16/ IOP,IORIG,NYC,J,IPRINT,LDFJYC,MODE,IFP,IRE,ISTP,JFP,
     + NPH,INE
C
      MPS = MP - IOR
C
      DO II = 1,MP
        DO JJ = 1,MP
          MREG(II,JJ) = 0
        ENDDO
        XREG(II) = 0
      ENDDO
C
      IF(IOR.EQ.0) THEN 
        DO II = 1,MPS
          DO JJ = 1,MP
          IF(II.EQ.JJ) THEN
            IREG(II,JJ) = 1
          ELSE
            IREG(II,JJ) = 0
          ENDIF
        ENDDO
      ENDDO
      ELSEIF(IOR.EQ.1) THEN
        DO II = 1,MPS
          DO JJ = 1,MP
          IF(JJ.LT.II.OR.JJ.GT.(II+1)) IREG(II,JJ) = 0
          IF(JJ.EQ.II) IREG(II,JJ) = -1
          IF(JJ.EQ.(II+1)) IREG(II,JJ) = 1
          ENDDO
        ENDDO
      ELSEIF(IOR.EQ.2) THEN
        DO II = 1,MPS
          DO JJ = 1,MP
            IF(JJ.LT.II.OR.JJ.GT.(II+3)) IREG(II,JJ) = 0
            IF(JJ.EQ.II) IREG(II,JJ) = -1
            IF(JJ.EQ.(II+1)) IREG(II,JJ) = 2
            IF(JJ.EQ.(II+2)) IREG(II,JJ) = -1
            IF(II.GE.(MP-1)) IREG(II,JJ) = 0
          ENDDO
        ENDDO
      ELSEIF(IOR.EQ.3) THEN
        DO II = 1,MPS
          DO JJ = 1,MP
            IF(JJ.LT.II.OR.JJ.GT.(II+3)) IREG(II,JJ) = 0
            IF(JJ.EQ.II) IREG(II,JJ) = -1
            IF(JJ.EQ.(II+1)) IREG(II,JJ) = 3
            IF(JJ.EQ.(II+2)) IREG(II,JJ) = -3
            IF(JJ.EQ.(II+3)) IREG(II,JJ) = 1
            IF(II.GE.(MP-2)) IREG(II,JJ) = 0
          ENDDO
        ENDDO
      ENDIF
C
C     CALCULATE H = BTRANSPOSE*B
      DO IJ = 1,MP
        DO JK = 1,MP
          DO MI = 1,MPS 
        MREG(IJ,JK) = MREG(IJ,JK) + IREG(MI,IJ)*IREG(MI,JK)
          ENDDO
        ENDDO
      ENDDO
C
      DP3 = P33A
C
C     CALCULATE TRANS(BP)*(BP)
      XRES = 0
      REGS = 0
      DO IJ = 1,MPS
        DO JK = 1,MP
          XRES = XRES + IREG(IJ,JK)*XR(JK)
        ENDDO
          REGS = REGS + XRES*XRES
      ENDDO
C
      DO IJ = 1,MP
        DO JK = 1,MP
           IF(ING.EQ.0) THEN
          XREG(IJ) = XREG(IJ) + MREG(IJ,JK)*XR(JK)
           ELSE
          XREG(IJ) = XREG(IJ) + MREG(IJ,JK)*XR(JK)*DEXP(-XR(JK))
           ENDIF
        ENDDO
      ENDDO
C
47    FORMAT(1X,19I4)
C
      IF(IFLAG.EQ.1) THEN
        DO 10 I = 1,MP
            FVEC(I+MK) = DP3*XREG(I)
10      CONTINUE
      ELSE
        DO 30 I = 1,MP
          DO 20 J = 1,MP
            FJAC(I+MK,J) = DP3*MREG(I,J)
20        CONTINUE
30      CONTINUE
      ENDIF
C
876   CONTINUE
      RETURN
      END
