C
C     ALL COMMON BLOCKS FROM LV0.FOR
C
C       MODIFIED FOR LEVMpy JEREMY SMITH 3/31/2017
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*1 DATTYP,DFIT,PFIT,FUN,DATTYQ,DATTYC
      INCLUDE 'SIZE.INC'
C
      COMMON /CM1/ FREQ(NPTS),M,DATTYP
      COMMON /CM2/ Y(NPT2),R(NPT2),FJ(NPT2),P(NTOT),DRSS,ROE,RKE,
     + NS(NPAFR),NFREE(NTOT),N,ICNT,MN,IRCH,IXI,DATTYQ
      COMMON /CM3/ CELCAP,FUN,DFIT,PFIT
      COMMON /CM4/ FNORM
      COMMON /CM5/ TOMEGA,PHICOM,IV
      COMMON /CM9/ TOMEGAX,PHIX,ICHG,IWTX
      COMMON /CM10/ EPSG,IZR
      COMMON /CM11/ MQY,ISPR,ICX,NDF,FQQ
      COMMON /CM12/ CLCAP,ATEMP,WF,MAXFEV,ICF,MDE,JCDX
      COMMON /CM13/ RX,TX,UX,PHI,XXM1,XX1,XX2,XX3,RN,AIN,ICAV,NELEM,NCH
      COMMON /CM14/ FJACC(NPT2,NPAFR)
      COMMON /CM16/ IOP,IORIG,NYC,JY,IPRINT,LDFJAC,MODE,IFP,IRE,ISTP,
     + JFP,NPH,INE
      COMMON /CM18/ SDWC,SDRC,DIAG(NPAFR),IPAR,IOPR
      COMMON /CM20/ RXSD(NPAFR)
      COMMON /CM29/ GAM,SN1,CS1,PII,CDN,QPI,MDEX
      COMMON /CM34/ MDA,IWT,IXW,INFP,IPL
      COMMON /CM35/ JIT,IPF,NPRINT
      COMMON /CM36/ SHW(NPT2),ISW
      COMMON /CM39/ P8,P18
      COMMON /CM47/ ICOUNT
      COMMON /CM55/ PX1,PX41,PX45
      COMMON /CM73/ P39
      COMMON /CM78/ DATTYC
      COMMON /CM79/ YTT(NPT2)
      COMMON /CM80/ OMEGAT,PSI,GAMA,DELI
C
      END