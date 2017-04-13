      SUBROUTINE LOGAC(AL,ACL)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 AL,ACL
C
      IF(AL.GT.1.D-8) THEN
        ACL = DLOG(1.D0 + AL)
      ELSE
        ACL = AL*(1.D0 - 0.5D0*AL)
      ENDIF
      RETURN
      END
