C
C     Selects the required distributed element and calls DISTEL
C
      SUBROUTINE SDEA(OMEGA,IOMEGA,RA,CA,RDE,TDE,UDE,PDE,NDE,ZX)
      IMPLICIT REAL*8(A-H,O-Z)
      EXTERNAL DISTEL
      COMPLEX*16 ZX,IOMEGA,DISTEL,ZA
      INTEGER NDE
      REAL*8 OMEGA,RA,CA,RDE,TDE,UDE,PDE
C
      IF(CA.EQ.0) THEN
            ZA = RA
      ELSEIF(RA.EQ.0) THEN
            ZA = 1.D0/(IOMEGA*CA)
      ELSE
            ZA = RA/(1.D0 + IOMEGA*CA*RA)
      ENDIF
C
      IF(NDE.NE.0) THEN
        ZX = DISTEL(RDE,TDE,UDE,PDE,NDE,OMEGA)
        IF(ZX.NE.(0,0)) THEN
            IF(ZA.NE.(0,0)) ZX = ZA*ZX/(ZX+ ZA)
        ELSE
            ZX = ZA
        ENDIF
      ELSE
        ZX = ZA
      ENDIF
C
      RETURN
      END
