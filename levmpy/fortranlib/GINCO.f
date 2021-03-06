C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE GINCO(A,X,FINCGM,SINCGM)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A,X,FINCGM,SINCGM
      DATA B0,B1,B2,B3,B4,B5/-.5772156649D0,.99999193D0,-.24991055D0,
     + 5.519968D-2,-9.76004D-3,1.07857D-3/
C
      IF(X.LT.0.) THEN
        WRITE(*,*) 'X IS LESS THAN ZERO'
        RETURN
C        STOP
      ENDIF
      IFLG = 0
      IF(A.GT.0.D0.AND.X.LT.A+1.D0) THEN
        CALL GSER(GAMSER,A,X,GLN)
C
C     FOLLOWING CHOICE GIVES CORRECT DIFFERENCES OF INC GAMMA FUNCTIONS
C     FOR SMALL X ARGUMENTS
C
        IFLG = 1
        GAMMQ = 1.D0 - GAMSER
        GAMMQD = - GAMSER
        DGLN = DEXP(GLN)
        FINCGM = DGLN*GAMMQ
        SINCGM = DGLN*GAMMQD        
      ELSEIF(A.EQ.0.D0.AND.X.LT.0.51) THEN
        FINCGM = B0+X*(B1+X*(B2+X*(B3+X*(B4+X*B5))))-DLOG(X)
        N = 0
        GAMMQ = 0.D0
        SINCGM = FINCGM
      ELSE
        CALL GCF(GAMMQ,A,X,GLN)
        IFLG = 2
        FINCGM = GAMMQ
        SINCGM = FINCGM
      ENDIF
      RETURN
      END
