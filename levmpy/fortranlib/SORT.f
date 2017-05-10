C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE SORT(TAU,ISP,NN)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER ISP,NN
      REAL*8 TAU
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
15      CONTINUE
        IF(IND.NE.I) THEN
              ITMP = ISP(I)
            ISP(I) = ISP(IND)
          ISP(IND) = ITMP
        ENDIF
20    CONTINUE
      RETURN
      END
