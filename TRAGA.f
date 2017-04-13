C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      SUBROUTINE TRAGA(AA,YY,PP,QQ,IL,TRSY)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER IL
      REAL*8 AA,YY,PP,QQ,TRSY
      DIMENSION XX(2,2),AA(2,2),PG(2),GC(2),PP(2),QQ(2),
     +  TRSY(2),SC(2)
C
C   Index IA: Make both PGs, phis; IY: two parts of PG, x args.
      DO 111 IA = 1,2
        IF(YY.GT.0.D0) THEN  
          IF(YY.LT.2.0D-12*PP(2)) THEN
            TRSY(IL) = DNORMI
            GOTO 956
              ENDIF 
                XX(IA,1) = PP(IA)*YY
                XX(IA,2) = QQ(IA)*YY
            DO 222 IY = 1,2
              CALL GINCO(AA(IA,IL),XX(IA,IY),GC(IY),SC(IY))
222         CONTINUE            
            DGC = GC(1) - GC(2)
            DSC = SC(1) - SC(2)
            IF(DGC.GT.0.D0) THEN
                GINC = DGC
            ELSE
                GINC = DSC
            ENDIF
            PG(IA) = YY**(-AA(IA,IL))*GINC
            ELSE
            IF(AA(IA,IL).EQ.0.D0) THEN
              PG(IA) = DLOG(QQ(IA)/PP(IA))                
            ELSE
              PG(IA) = (QQ(IA)**AA(IA,IL) - PP(IA)**AA(IA,IL))/AA(IA,IL)
            ENDIF
          ENDIF
111     CONTINUE
        TRSY(IL) = PG(1) + PG(2)
        IF(YY.EQ.0.D0) DNORMI = TRSY(IL)
956    CONTINUE
      RETURN
      END
