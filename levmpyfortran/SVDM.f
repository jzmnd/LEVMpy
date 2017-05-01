C
C     DECOMPOSITION OF MATRIX
C
C           ND : number of functions (2*M for complex) (IN)
C           MP : number of free parameters (IN)
C       SGMASQ : 
C            X : solution vector (IN)
C          PEX : parameter list (IN)
C           NS : 
C        IORIG : 
C         ISTP : 
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE SVDM(ND,MP,SGMASQ,X,PEX,NS,IORIG,ISTP)
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'SIZE.INC'
      INTEGER ND,MP,NS,IORIG,ISTP
      REAL*8 SGMASQ,X,PEX
      DIMENSION CVM(NPAFR,NPAFR),CINV(NPAFR,NPAFR),CORR(NPAFR,NPAFR),
     + X(*),NS(*),INDX(NPAFR),XSD(NPAFR),PEX(*),PRELE(NPAFR)
      COMMON /CM11/ MQY,ISPR,ICX,NDF,FQQ
      COMMON /CM14/ V(NPT2,NPAFR)
      COMMON /CM20/ RXSD(NPAFR)
C
      DO 14 I=1,MP
        DO 13 J=1,I
          SUM=0.D0
          DO 12 KX=1,ND
            SUM=SUM+V(KX,I)*V(KX,J)
12        CONTINUE
          CVM(I,J)=SUM
          CVM(J,I)=SUM
13      CONTINUE
14    CONTINUE
        DO 121 I = 1,MP
          DO 111 J = 1,MP
            CINV(I,J) = 0.D0
111     CONTINUE
          CINV(I,I) = 1.D0
121     CONTINUE
        CALL LUDCMP(CVM,MP,NPAFR,INDX,D,ISD)
        DO 131 J = 1,MP
          CALL LUBKSB(CVM,MP,NPAFR,INDX,CINV(1,J))
131     CONTINUE
C
      IF(ISTP.GT.0) THEN  
C        OUTPUT PARAMETER NUMBERS AND ESTIMATED VALUES
C         WRITE(3,16)
        WRITE(*,16)
        IF(ISD.EQ.1) THEN
          WRITE(*,*) '!!!!!**  SINGULAR MATRIX - BEWARE  **!!!!!'
          WRITE(*,*)
        ENDIF
        WRITE(*,41)
C         WRITE(3,41)  
41      FORMAT(3X,'ORIGINAL PARAMETER NUMBERS',/,
     +3X,'PARAMETER ESTIMATES',/,3X,'ESTIMATED STD DEV OF PARAMETERS'
     +,/,3X,'ESTIMATED RELATIVE STD DEV OF PARAMETERS')
        IF(IORIG.GT.0) THEN
          WRITE(*,48)
C        WRITE(3,48)
      ENDIF
48    FORMAT(3X,'RELATIVE ERRORS OF PARAMETERS')
      WRITE(*,16)
C        WRITE(3,16)
C
        PDAV = 0.D0
        PDRMS = 0.D0
        PSABS = 0.D0
        PSRMS = 0.D0
C
      DO 31 JJ =1,MP
        IF(CINV(JJ,JJ).LE.0.D0) CINV(JJ,JJ) = DABS(CINV(JJ,JJ))
          XSD(JJ) = DSQRT(CINV(JJ,JJ)*SGMASQ)
        IF(X(JJ).NE.0.D0) THEN
          RXSD(JJ) = DABS(XSD(JJ)/X(JJ))
        ELSE
          RXSD(JJ) = XSD(JJ)
        ENDIF
        IF(JJ.LE.MQY) THEN
          PDAV = PDAV + RXSD(JJ)
          PDRMS = PDRMS + RXSD(JJ)**2
        ENDIF
        IF(IORIG.GT.0) THEN
          IF(PEX(NS(JJ)).NE.0.D0) THEN    
            PRELE(JJ) = (X(JJ) - PEX(NS(JJ)))/PEX(NS(JJ))
          ELSE
            PRELE(JJ) = 0
          ENDIF
          IF(JJ.LE.MQY) THEN
            PSABS = PSABS + DABS(PRELE(JJ))
            PSRMS = PSRMS + PRELE(JJ)**2
          ENDIF
        ENDIF
31    CONTINUE
        PDAV = PDAV/MQY
        PDRMS = DSQRT(PDRMS/MQY)
C     OUTPUT PARAMETERS, ESTIMATED STD DEV, AND REL STD DEV
C
      MQP = INT(MP/6) + 1
      DO 762 JQ = 1,MQP
        JQP = 6*JQ
        IXS = JQP - 5
        IF(MP.LE.JQP) THEN
            IXL = MP
        ELSE
            IXL = JQP
        ENDIF
        IF(IXS.LE.MP) THEN
          WRITE(*,'(2X,I7,6I12)') (NS(I),I=IXS,IXL)
C        WRITE(3,'(2X,I7,6I12)') (NS(I),I=IXS,IXL)
          WRITE(*,'(2X,1P,6D12.4)') (X(I),I=IXS,IXL)
C        WRITE(3,'(2X,1P,6D12.4)') (X(I),I=IXS,IXL)
          WRITE(*,'(2X,1P,6D12.4)') (XSD(I),I=IXS,IXL)
C        WRITE(3,'(2X,1P,6D12.4)') (XSD(I),I=IXS,IXL)
          WRITE(*,'(2X,1P,6D12.4)') (RXSD(I),I=IXS,IXL)
C        WRITE(3,'(2X,1P,6D12.4)') (RXSD(I),I=IXS,IXL)
          IF(IORIG.GT.0) THEN
            WRITE(*,'(2X,1P,6D12.4)') (PRELE(I),I=IXS,IXL)
C          WRITE(3,'(2X,1P,6D12.4)') (PRELE(I),I=IXS,IXL)
          ENDIF
C        WRITE(3,16)
        WRITE(*,16)
        ENDIF
C
762   CONTINUE
C
      IF(MQY.GT.0) THEN
        WRITE(*,39) PDAV,PDRMS
C        WRITE(3,39) PDAV,PDRMS
39    FORMAT(6X,'PDAV=',2X,1P,1D12.4,6X,'PDRMS=',2X,1P,1D12.4)
C
C        WRITE(*,173) NDF,FQQ 
173   FORMAT(7X,'NDF=',I4,7X,'FQF=',1P,1D15.6) 
C     temporarily changed 1.d3 to 1.d40  8/17/94 !!!!!!!!!!!!!
        IF(PDAV.GT.1.D40) THEN
          WRITE(*,749)
749   FORMAT(/2X,'****  NOTE LARGE VALUE OF ONE OR MORE RELATIVE
     + STANDARD DEVIATIONS ****'/,
     + 14X,'SIGN OF WRONG OR POOR INITIAL CHOICE OF A FREE PARAMETER,'
     + /,21X,'ALSO LACK OF CONVERGENCE AND/OR SINGULAR MATRIX'//)
C     CALL NOTESUB
        ENDIF
C
        IF(IORIG.GT.0) THEN
          PSABS = PSABS/MQY
          PSRMS = DSQRT(PSRMS/MQY)
          WRITE(*,43) PSABS,PSRMS
C          WRITE(3,43) PSABS,PSRMS
43    FORMAT(6X,'PRSDAV=',1P,1D12.4,6X,'PRSDRMS=',1P,1D12.4)
        ENDIF
      ENDIF
16    FORMAT(/)
C
C            OUTPUT VARIANCE-COVARIANCE MATRIX: CINV
C            OUTPUT PARAMETER CORRELATION MATRIX: CORR
C
      IF(ICX.EQ.1) WRITE(*,42)
      WRITE(*,42)  
42    FORMAT(/3X,'ESTIMATED PARAMETER CORRELATION MATRIX')
C
      IBAD = 0
      DO 763 JQ = 1,MQP
        JQP = 6*JQ
        IXS = JQP - 5
        IF(MP.LE.JQP) THEN
            IXL = MP
        ELSE
            IXL = JQP
        ENDIF
C
        DO 17 I = IXS,IXL
          DO 18 J = 1,I
            CORP = CINV(I,I)*CINV(J,J)
            IF(CORP.LE.0.D0) THEN
                CORP = 1.D0
                IBAD = IBAD + 1
            ENDIF
            CORR(I,J) = CINV(I,J)/DSQRT(CORP)
18        CONTINUE
          IF(ICX.EQ.1) WRITE(*,'(2X,6D12.4)') (CORR(I,J),J=1,I)
          WRITE(*,'(2X,6D12.4)') (CORR(I,J),J=1,I)
17      CONTINUE
769   FORMAT(/2X,'  #!!!!  BAD CORRELATION MATRIX  !!!!#'
     + /5X,'PROBABLY AT LEAST ONE FREE PARAMETER IS NOT USED IN 
     + FITTING FUNCTION'/5X,'OR INITIAL VALUE OF A PARAMETER IS TOO FAR 
     + FROM THE CORRECT VALUE'//)
C
763   CONTINUE
        IF(IBAD.GE.1) THEN
          WRITE(*,769)
        ENDIF
      ENDIF
      RETURN
      END
