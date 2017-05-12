C
C     Given a matrix A(1:N,1:N), with physical dimension NP by NP, this
C       routine replaces it by the LU decomposition of a rowwise
C       permutation of itself.
C         A and N are input
C         A is output
C         INDX(1:N) is an output vector that records the row permutation
C           effected by the partial pivoting
C         D is output as ±1 depending on whether the number of row
C           interchanges was even or odd, respectively
C
C     This routine is used in combination with lubksb to solve linear
C       equations or invert a matrix.
C
      SUBROUTINE LUDCMP(A,N,NP,INDX,D,ISD)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N,NP,INDX,IMAX,ISD
      REAL*8 A,D,TINY
      PARAMETER (TINY=1.0D-20)
      INCLUDE 'SIZE.INC'
      DIMENSION A(NP,NP),INDX(N),VV(NPAFR)
C
      D = 1.D0
      ISD = 0
      IMAX = 1
      DO 12 I=1,N
        AAMAX=0.D0
        DO 11 J=1,N
          IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
11      CONTINUE
        IF (AAMAX.EQ.0.D0) THEN
C          WRITE(*,*) '*** SINGULAR MATRIX ***'
          AAMAX = 1.D0
          ISD = 1
        ENDIF
        VV(I)=1.D0/AAMAX
12    CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
              A(I,J)=SUM
            ENDIF
14        CONTINUE
        ENDIF
        AAMAX=0.D0
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*DABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
16      CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
17        CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.D0) A(J,J)=TINY
          DUM=1.D0/A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
18        CONTINUE
        ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.D0) A(N,N)=TINY
      RETURN
      END
