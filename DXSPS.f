C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE DXSPS(XA,YCC,ARR,WF,NN,IP38A)
      IMPLICIT REAL*8(A-H,O-Z)
      COMPLEX*16 YCC(*),ARC(20),ARR,ZT,ERROR
      INTEGER NN,IP38A,ISP(20)
      DOUBLE PRECISION WF,XB(101),Q(101),XA(*)
      COMMON /CM73/ P39
      COMMON /CM47/ ICNT
      DATA PT5/0.5D0/
C
C   JRM 6/22/94;12/23/94
C
      ARR = DCMPLX(0.D0,0.D0)
C
C       DO SIMPLE SUM FOR DV OR FOR EQUAL-SPACING QUADRATURE
      IF(IP38A.EQ.0.OR.IP38A.EQ.1) THEN
        DO IU = 1, NN
            ARR = ARR + YCC(IU)
        ENDDO
C   
        IF(IP38A.EQ.1) ARR = ARR*WF
C
        GOTO 339
      ENDIF
C
C       ######################################
C
      IF(IP38A.LE.4) THEN
        CALL SORT(XA,ISP,NN)
        DO IP = 1,NN
            XB(IP) = XA(ISP(IP))
            ARC(IP) = YCC(ISP(IP))
            Q(IP) = DLOG(XB(IP))
        ENDDO
C
        DO IP = 1,NN
          IF(IP38A.EQ.2.OR.IP38A.EQ.4) THEN
            IF(IP.EQ.1.OR.IP.EQ.NN) THEN
              YCC(IP) = 2.D0*ARC(IP)
            ELSE
              YCC(IP) = ARC(IP)
            ENDIF
          ELSE
              YCC(IP) = ARC(IP)
          ENDIF
        ENDDO
C
      ELSE
        RETURN
C        STOP
      ENDIF
C
C   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C   VARIABLE SPACING TRAP QUADRATURE
C
      IF(IP38A.EQ.2) THEN
        DO 118 K = 1,NN-1
            J = K + 1
          DXV = Q(J) - Q(K)
          ARR = ARR + DXV*(YCC(J) + YCC(K))
118    CONTINUE
        ARR = ARR*PT5
C
        GOTO 339
      ENDIF
C
      IF(IP38A.EQ.3) THEN
        CALL CUBINT(Q,YCC,NN,1,NN,ARR,ERROR,ZT,IND)
        ELSEIF(IP38A.EQ.4) THEN
        CALL CUBINT(Q,YCC,NN,1,NN,ZT,ERROR,ARR,IND)
      ENDIF   
        IF(IND.NE.1) THEN
          RETURN
        ENDIF
        DXL = Q(2) - Q(1)
        DXH = Q(NN) - Q(NN-1)
C
      ARR = ARR + P39*(DXL*ARC(1) + DXH*ARC(NN))
339   RETURN
      END
