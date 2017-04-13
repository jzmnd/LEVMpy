C
C     Calls either GSER or GCF
C       If X < A+1 then call GSER
C       If X >= A+1 then call GCF
C
      FUNCTION GAMMQ(A,X)
      DOUBLE PRECISION A,GAMMQ,X
      DOUBLE PRECISION GAMMCF,GAMSER,GLN
      IF(X.LT.0.D0.OR.A.LE.0.D0) THEN
        WRITE(*,*) 'BAD ARGUMENTS IN GAMMQ'
        GAMMQ = 0.D0
        RETURN
C        STOP
      ENDIF
      IF(X.LT.A+1.D0) THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ = 1.D0-GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMQ = GAMMCF
      ENDIF
      RETURN
      END
