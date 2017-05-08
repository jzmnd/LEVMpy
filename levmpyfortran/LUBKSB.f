C
C     Solves the set of N linear equations A · X = B. Here A is input,
C       not as the matrix A but rather as its LU decomposition,
C       determined by the routine ludcmp.
C         INDX is input as the permutation vector returned by ludcmp
C         B(1:N) is input as the right-hand side vector B, and returns 
C           with the solution vector X.
C
C     A, N, NP, and INDX are not modified by this routine and can be
C       left in place for successive calls with different right-hand
C       sides B. This routine takes into account the possibility that 
C       B will begin with many zero elements, so it is efficient for
C       use in matrix inversion.
C
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N,NP,INDX
      REAL*8 A,B
      DIMENSION A(NP,NP),INDX(*),B(NP)
C
      II=0
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0) THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
11        CONTINUE
        ELSE IF (DABS(SUM).GT.1.D-200) THEN
          II=I
        ENDIF
        B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
13        CONTINUE
        ENDIF
        IF(A(I,I).EQ.0.D0) THEN
          A(I,I) = 1.D-10
        ENDIF
        B(I)=SUM/A(I,I)
14    CONTINUE
      RETURN
      END
