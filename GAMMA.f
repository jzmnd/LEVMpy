      SUBROUTINE GAMMA(ZZ,G1,HH)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 ZZ,G1,HH
      DATA AC/-0.57719163D0/, AD/0.98820589D0/,AE/-0.897056937D0/,
     +  AF/0.918206857D0/,AG/-0.75670408D0/,AH/0.48219939D0/,
     +  AI/-0.193527818D0/,AJ/0.035868343D0/
C
      IF(ZZ.NE.0.D0) THEN
C           HERE G1(Z) IS GAMMA(1 + Z); HH = GAMMA(Z)
      G1 = 1.D0 + ZZ*(AC + ZZ*(AD + ZZ*(AE + ZZ*(AF + ZZ*(AG +
     +  ZZ*(AH + ZZ*(AI + ZZ*AJ)))))))
            HH = G1/ZZ
      ELSE
        G1 = 1.D0
        HH = 666
C        PAUSE 66
      ENDIF
C
      RETURN
      END
