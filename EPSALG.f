C
C     EPS ALGORITHM
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE EPSALG(NMAX,SS,SF)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER NMAX
      REAL*8 SF,SS
      DIMENSION EX(114,114)
      DIMENSION SS(113)
C
      IF(NMAX.GT.113) THEN
        WRITE(*,*) 'USE SMALLER ODD NMAX'
C        READ(*,*)
        RETURN
C        STOP
      ENDIF
C
      IF(MOD(NMAX,2).EQ.0) NMAX = NMAX - 1
C
11    DO II = 1,2
        DO NN = 1,NMAX
            IF(II.EQ.1) EX(NN,II) = 0.D0
            IF(II.EQ.2) THEN
                EX(NN,II) = SS(NN)  
            ENDIF
        ENDDO
      ENDDO
C
      DO II = 2,NMAX
        DO NN = 1,(NMAX - II + 1)
          IF(EX(NN+1,II).EQ.EX(NN,II)) THEN
            NMAX = NMAX - 2
            IF(NMAX.GE.7) THEN
              GOTO 11
            ELSE
              WRITE(*,*) 'EPS ALGORITHM FAILS'
C              READ(*,*)
              RETURN
C              STOP
            ENDIF
          ELSE
          EX(NN,II+1) = EX(NN+1,II-1) + 1.D0/(EX(NN+1,II) - EX(NN,II))        
          ENDIF
        ENDDO
      ENDDO
      SF = EX(1,NMAX+1)
558   RETURN
      END
