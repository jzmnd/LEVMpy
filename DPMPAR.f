C     DPMPAR from MINPACK
C
C     This function provides double precision machine parameters
C     when the appropriate set of data statements is activated (by
C     removing the c from column 1) and all other data statements are
C     rendered inactive. Most of the parameter values were obtained
C     from the corresponding Bell Laboratories Port Library function.
C
C       I is an integer input variable set to 1, 2, or 3 which
C         selects the desired machine parameter. If the machine has
C         t base b digits and its smallest and largest exponents are
C         emin and emax, respectively, then these parameters are
C
C         dpmpar(1) = b**(1 - t), the machine precision,
C
C         dpmpar(2) = b**(emin - 1), the smallest magnitude,
C
C         dpmpar(3) = b**emax*(1 - b**(-t)), the largest magnitude.
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      DOUBLE PRECISION FUNCTION DPMPAR(I)
      INTEGER I
      INTEGER MCHEPS(4), MINMAG(4), MAXMAG(4)
      DOUBLE PRECISION DMACH(3)
      EQUIVALENCE (DMACH(1),MCHEPS(1))
      EQUIVALENCE (DMACH(2),MINMAG(1))
      EQUIVALENCE (DMACH(3),MAXMAG(1))
C
C     Machine constants for IEEE machines
C
      DATA DMACH(1) /2.22044604926D-16/
      DATA DMACH(2) /2.22507385852D-308/
      DATA DMACH(3) /1.79769313485D+308/
C
C     machine constants for 8088/8086 MACHINES
C
C      DATA MCHEPS(1)/0/,MCHEPS(2)/1018167296/
C      DATA MINMAG(1)/0/,MINMAG(2)/1048576/
C      DATA MAXMAG(1)/-1/,MAXMAG(2)/2146435071/
        DPMPAR = DMACH(I)
      RETURN
      END
