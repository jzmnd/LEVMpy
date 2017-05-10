C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE CUBINT(X,F,N,IA,IB,RESULT,ERROR,ZT,IND)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 X
      INTEGER N,IA,IB,IND
      COMPLEX*16 F(N),D1,D2,D3,R1,R2,R3,R4,RESULT,ERROR,ZT,RELERR
      DIMENSION X(N)
C
      IND = 0
      IF(N.LT.4.OR.IA.LT.1.OR.IB.GT.N) RETURN
      IND = 1
      RESULT = 0.D0
      ERROR = 0.D0
      IF(IA.EQ.IB) RETURN
      IF(IA.LT.IB) GOTO 2
      IND = 1
      IT = IB
      IB = IA
      IA = IT
2     S = 0.D0
      C = 0.D0
      R4 = 0.D0
      J = N - 2
      IF(IA.LT.N-1.OR.N.EQ.4) J = MAX0(3,IA)
      K = 4
      IF(IB.GT.2.OR.N.EQ.4) K = MIN0(N,IB+2) - 1
      DO 1 I = J,K
        IF(I.GT.J) GOTO 5
        H2 = X(J-1) - X(J-2)
        D3 = (F(J-1) -F(J-2))/H2
        H3 = X(J) - X(J-1)
        D1 = (F(J) - F(J-1))/H3
        H1 = H2 + H3
        D2 = (D1 - D3)/H1
        H4 = X(J+1) - X(J)
        R1 = (F(J+1) - F(J))/H4
        R2 = (R1 - D1)/(H4 + H3)
        H1 = H1 + H4
        R3 = (R2 - D2)/H1
        IF(IA.GT.1) GOTO 8
        RESULT = H2*(F(1) + H2*(0.5D0*D3 - H2*(D2/6.D0 - (H2 +
     +          H3 + H3)*R3/12.D0)))
        S = -H2**3*(H2*(3.D0*H2 + 5.D0*H4) + 1.D1*H3*H1)/6.D1
        GOTO 8
5       H4 = X(I+1) - X(I)
        R1 = (F(I+1) -F(I))/H4
        R4 = H4 + H3
        R2 = (R1 - D1)/R4
        R4 = R4 + H2
        R3 = (R2 - D2)/R4
        R4 = (R3 - D3)/(R4 + H1)
8       IF(I.GT.IB.OR.I.LE.IA) GOTO 11
        RESULT = RESULT + H3*((F(I) + F(I-1))*0.5D0 - H3*H3*(
     +          D2 + R2 + (H2 -H4)*R3)/12.D0)
        C = H3**3*(2.D0*H3*H3 + 5.D0*(H3*(H4 + H2) + 2.D0*H2*H4)
     +          )/1.2D2 
        ERROR = ERROR + (C + S)*R4
        IF(I.EQ.J) GOTO 14
        S = C
        GOTO 15 
14      S = S + C + C
        GOTO 15
11      ERROR = ERROR + R4*S
15      IF(I.LT.K) GOTO 20
        IF(IB.LT.N) GOTO 22
        RESULT = RESULT + H4*(F(N) - H4*(0.5D0*R1 + H4*(R2/6.D0
     +           + (H3 + H3 + H4)*R3/12.D0)))
        ERROR = ERROR - H4**3*R4*(H4*(3.D0*H4 + 5.D0*H2) + 1.D1*
     +          *H3*(H2 + H3 + H4))/6.D1
22      IF(IB.GE.N-1) ERROR = ERROR + S*R4
        GOTO 1
20      H1 = H2
        H2 = H3
        H3 = H4
        D1 = R1
        D2 = R2
        D3 = R3
1     CONTINUE
      IF(IND.NE.1) THEN
        RETURN
      ENDIF
      ZT = RESULT + ERROR
      RELERR = ERROR/(RESULT + (1.D-30,1.D-29))
      RETURN
      END
