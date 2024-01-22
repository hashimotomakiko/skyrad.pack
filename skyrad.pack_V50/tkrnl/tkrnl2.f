C$OUTF tkrnl2.f
C$LIBF /Users/makibo/Fortran/LIB/BASLIB/A
C$LIBF /Users/makibo/Fortran/LIB_Nakajima/LBC1/A
C$LIBF /Users/makibo/Fortran/LIB_Nakajima/LBM1/A
C$LIBF /Users/makibo/Fortran/LIB_Nakajima/LBM2/A
C$LIBF /Users/makibo/Fortran/LIB_Nakajima/LBR1/A
C
C --- Note
C   Make the table of kernel matrix
C   Different read-format for refractive indices
C
C --- History
C 89. 8.17  CREATED FROM PHS:CKRNL (IBM TO MAC)
C 90.11. 6  INITIALIZE ERR=' ', GOTO 1
C 97. 2.15  Modified by Rao
C 98. 6.15  VCI= -x.xx
C *******************************************************************
C 01.01.29 Overhauled by M.Yamano for analyzing system Ver.3 & Ver.4.
C *******************************************************************
C 01.09.04  Changed input parameters for refractive index
C 02.05.10  Renewed by M.Yamano for Ver.5.
C 03.10.25  Renewed.
c 04.12.06  Start on Mac.
C
C --- I/O files
C  tkrnl.par  I    input parameters for creating kernels
C  MIEKER     O    output kernel datafile
C ----------------------------------------------------------------------
      PARAMETER (KNANG =200)
      PARAMETER (KSIZE =100)
      PARAMETER (KCR   =20)
      PARAMETER (KCI   =20)
      PARAMETER (KNANG1=KNANG+1,KSIZE1=KSIZE+1)
      PARAMETER (KNN=10)
C
C
      CHARACTER ERC*64
      CHARACTER OFIL*80
C-----AREA FOR -MIES-
      DIMENSION Q(KNANG1,4,KSIZE),THETA(KNANG),XX(KSIZE1)
C-----WORKING
      DIMENSION XB(KNN),NB(KNN),VCR(KCR),VCI(KCI)
C
      DATA IFORM /2/
C
      ERC=' '
      IUI=1
      OPEN (IUI,FILE='tkrnl.par',STATUS='OLD')
      READ(IUI,*) OFIL
      CALL UTLSPC(OFIL,NL)
      IF (NL.LE.0) THEN
         ERC='illegal kernel filename'
         GOTO 90
      ENDIF
C
C-----SET ANGLES
      READ(IUI,*) IJOB,NN
      READ(IUI,*) (XB(I),I=1,NN)
      READ(IUI,*) (NB(I),I=1,NN-1)
      CALL XSET(IJOB,NN,XB,NB,KNANG,NANG,THETA,ERC)
      IF(ERC.NE.' ') GOTO 90
C
C-----SET SIZE PARAMETER
      READ(IUI,*) IJOB,NN
      READ(IUI,*) (XB(I),I=1,NN)
      READ(IUI,*) (NB(I),I=1,NN-1)
      CALL XSET(IJOB,NN,XB,NB,KSIZE1,NSIZE,XX,ERC)
      IF(ERC.NE.' ') GOTO 90
      INTVL=NSIZE-1
      DEL=LOG(XB(2)/XB(1))/INTVL
C
C-----READ POLARIZATION STATES, AVERAGING NUMBER, NUMBER OF INDEX
      READ(IUI,*) IPOL
      READ(IUI,*) NXI
C
C --- Changed.(01/09/04)
      READ(IUI,*)
      READ(IUI,*) NCR,IFGR,CRMIN,CRMAX
      IF (NCR.GT.KCR) THEN
         ERC='NCR.GT.KCR !'
         GOTO 90
      ENDIF
      IF (IFGR.EQ.0) THEN
         READ(IUI,*) (VCR(I),I=1,NCR)
         CRMIN=VCR(1)
         CRMAX=VCR(NCR)
      ENDIF
      READ(IUI,*) NCI,IFGI,CIMIN,CIMAX
      IF (NCI.GT.KCI) THEN
         ERC='NCI.GT.KCI !'
         GOTO 90
      ENDIF
      IF (IFGI.EQ.0) THEN
         READ(IUI,*) (VCI(I),I=1,NCI)
         CIMIN=VCI(1)
         CIMAX=VCI(NCI)
         DO I=1,NCI
           VCI(I)=-VCI(I)
         ENDDO
      ENDIF
      IF ((NCR.LE.3).OR.(NCI.LE.4)) THEN
         ERC='illegal setting (NCR or NCI) !'
         GOTO 90
      ENDIF
      IF ((CRMIN.GE.CRMAX).OR.(CIMIN.GE.CIMAX)) THEN
         ERC='illegal setting (CRMIN or CRMAX or CIMIN or CIMAX) !'
         GOTO 90
      ENDIF
      IF (IFGR.NE.0) THEN
         IF (IFGR.GT.0) THEN
            DL=(CRMAX-CRMIN)/FLOAT(NCR-1)
            CRL=CRMIN-DL
            DO I=1,NCR
              CRL=CRL+DL
              VCR(I)=CRL
            ENDDO
         ELSE
            DL=LOG(CRMAX/CRMIN)/FLOAT(NCR-1)
            CRL=LOG(CRMIN)-DL
            DO I=1,NCR
              CRL=CRL+DL
              VCR(I)=EXP(CRL)
            ENDDO
         ENDIF
      ENDIF
      IF (IFGI.NE.0) THEN
         IF (IFGI.GT.0) THEN
            DL=(CIMAX-CIMIN)/FLOAT(NCI-1)
            CIL=CIMIN-DL
            DO I=1,NCI
              CIL=CIL+DL
              VCI(I)=-CIL
            ENDDO
         ELSE
            DL=LOG(CIMAX/CIMIN)/FLOAT(NCI-2)
            CIL=LOG(CIMIN)-DL
            VCI(1)=0.0
            DO I=2,NCI
              CIL=CIL+DL
              VCI(I)=-EXP(CIL)
            ENDDO
         ENDIF
      ENDIF
C
C --- Output different from the version 2
C
      IUO=2
      OPEN (IUO,FILE=OFIL(1:NL),STATUS='UNKNOWN')
      WRITE(IUO,10) IFORM
   10 FORMAT(I3,' : IFORM')
      WRITE(IUO,11) NCR,NCI
   11 FORMAT(2I3,' : NCR NCI')
      WRITE(IUO,12) (VCR(I),I=1,NCR)
   12 FORMAT(7F11.3)
      WRITE(IUO,13) (VCI(I),I=1,NCI)
   13 FORMAT(1P7E11.2)
C
C-----COMPLEX REFRACTIVE INDEX LOOP
C
      IPOL0=4
      IREF=0
      DO ICI=1,NCI
        DO ICR=1,NCR
          IREF=IREF+1
          CRD=VCR(ICR)
          CID=VCI(ICI)
C
C-----MIE CALCULATION WITH INTERVAL AVERAGING
          CALL MIES(IPOL0,CRD,CID,NSIZE,NXI,XX,NANG,THETA,KNANG1,Q,ERC)
          IF(ERC.NE.' ') GOTO 90
C
          DO I=1,INTVL
            CO=3.0/4.0/SQRT(XX(I)*XX(I+1))*DEL
            IF (IPOL0.GE.1) THEN
               DO K=1,IPOL0
                 DO J=1,NANG
                   Q(J,K,I)=Q(J,K,I)*CO
                 ENDDO
               ENDDO
            ENDIF
            DO K=1,2
              Q(NANG+1,K,I)=Q(NANG+1,K,I)*CO
            ENDDO
          ENDDO
          IF (IPOL.EQ.1) THEN
             DO I=1,INTVL
               DO J=1,NANG
                 Q(J,1,I)=(Q(J,1,I)+Q(J,2,I))/2.0
               ENDDO
             ENDDO
          ENDIF
C
C-----WRITE HEADING DATA
          IF (IREF.EQ.1) THEN
             WRITE(IUO,20) INTVL,NANG,IPOL
   20        FORMAT(3I4)
             WRITE(IUO,30) (SQRT(XX(I)*XX(I+1)),I=1,INTVL)
             WRITE(IUO,30) (THETA(I),I=1,NANG)
   30        FORMAT(1P7E11.3)
          ENDIF
C
C-----WRITE RESULTS
          WRITE(IUO,40) IREF
   40     FORMAT('CM',I4)
          WRITE(IUO,50) CRD,CID
   50     FORMAT(1P2E11.3)
          DO I=1,INTVL
            DO K=1,IPOL
              WRITE(IUO,30) (Q(J,K,I),J=1,NANG)
            ENDDO
          ENDDO
          DO K=1,2
            WRITE(IUO,30) (Q(NANG+1,K,I),I=1,INTVL)
          ENDDO
        ENDDO
      ENDDO
      GOTO 100
C
   90 WRITE(IUO,95) ERC
   95 FORMAT(' ERROR: ',A64)
C
 100  STOP
      END
C
      SUBROUTINE MIES(IP,CR,CI,NSIZE,NXI,X,NANG,THETA,KNA1,Q,ERR)
C TABLE OF EFFICIENCY FACTOR FOR EXTINCTION, SCATTERING AND ITS ANGULAR
C  DISTRIBUTIONS. AVERAGING OVER SIZE INTERVALS TO SUPPRESS
C  THE FLUCTUATION.
C--- HISTORY
C 89. 8.10  USE MMIES (MAC)
C--- INPUT
C IP      I         IP=0  QEXT, QABS
C                   IP=2  QEXT, QABS, Q1, Q2
C                   IP=4  QEXT, QABS, Q1, Q2, Q3, Q4.
C CR      I         REFRACTIVE INDEX = CR + CI I
C CI      I           CI: GENERALLY NEGATIVE.
C NSIZE   I         NO. OF SIZE PARAMETERS AT DIVISION POINTS.
C NXI     I         NO. OF SUBDIVISIONS   ABOUT 10.
C X     R(NSIZE)    X = KA SIZE PARAMETER
C                   I-TH INTERVAL   ( X I, X I+1 )
C NANG    I         NUMBER OF SCATTERING ANGLES.
C THETA R(NANG)     SCATTERING ANGLES IN DEGREE,  I=1, NANG.
C KNA1    I         ARGUMENT SIZE OF NANG+1.
C--- OUTPUT
C Q     R(KNA1,     QLIJ  GENERALIZED EFFICIENCY FACTORS
C       4,INTVL)    ANGULAR EFFICIECY  I=1,NANG, L=1,4,   J=1,INTVL
C                   EXTINCTION EFF.    I=NANG+1, L=1
C                   ABSORPTION                   L=2
C ERR     C*16      ERROR INDICATER
C--- PARAMTERS
C NMMAX     I       1.2*XMAX+5
C NMMAX2    I       1.2*REAL(CM)*XMAX+20
C KNANG     I       DECLARED NUMBER OF SCATTERING ANGLES
C
      CHARACTER ERR*(*)
      DIMENSION X(*),THETA(*),Q(KNA1,4,*)
      PARAMETER (NMMAX =1500)
      PARAMETER (NMMAX2=2500)
      PARAMETER (KNANG =200)
      PARAMETER (KNANG1=KNANG+1)
C WORKING AREA
      DIMENSION Q1(KNANG1,4)
C
      ERR=' '
      INIT=1
      XMAX=X(NSIZE)
      IF(IP.LE.0) THEN
        NAS=NANG+1
        NAF=NAS
        IPP=2
      ELSE
        NAS=1
        NAF=NANG+1
        IPP=IP
      ENDIF
C
      DO 1 I=1,NSIZE-1
      DO 2 J=1,IPP
      DO 2 NA=NAS,NAF
    2 Q(NA,J,I)=0
      X1=X(I)
      X2=X(I+1)
      XM=0.5*(X1+X2)
      H=(X2-X1)/REAL(NXI)
      XX=X1-0.5*H
C MEAN OVER INTERVAL (X(I),X(I+1))
      DO 3 L=1,NXI
      XX=XX+H
      CALL MMIES(INIT,XMAX,IP,CR,CI,XX,NANG,THETA,KNANG1,Q1,ERR)
      IF(ERR.NE.' ') RETURN
        DO 4 J=1,IPP
        DO 4 NA=NAS,NAF
    4   Q(NA,J,I)=Q(NA,J,I)+Q1(NA,J)
    3 CONTINUE
C
      DO 5 J=1,IPP
      DO 5 NA=NAS,NAF
    5 Q(NA,J,I)=Q(NA,J,I) /REAL(NXI)
C
    1 CONTINUE
      RETURN
      END
      SUBROUTINE XSET(IJOB,NN,X,N,KNA,NA,XA,ERR)
C SET X-VALUES
C                   EXAMPLE
C  I-----------I---------------I-------------I    X-AXIS.
C X(1)        X(2)            X(3)          X(4)  KNOT VALUE.
C      N(1)=4      N(2)=3            N(3)=7       NUMBER OF SUBINTERVAL.
C
C  I--+--+--+--I---+---+---I-+-+-+-+-+-+-I        RESULTS.
C XA(1),(2), (3), ....                (NA=15)
C--- HISTORY
C 89. 8.19  MODIFIED
C--- INPUT
C IJOB      I         DIVIDE BY EQUAL SPACING ON LINEAR- (0) OR
C                       LOG- (1) AXIS.
C NN        I         NUMBER OF INTERVAL + 1.
C X       R(NN)       X-VALUE AT KNOT.
C N       I(NN)       NUMBER OF SUBINTERVAL IN EACH INTERVAL
C                      (X(I), X(I+1)).
C                      N(NN): DUMMY
C KNA       I         DECLARED SIZE OF XA .GE. NA.
C--- OUTPUT
C NA        I         TOTAL POINT NUMBER.
C XA      R(NA)       X-VALUES.
C ERR     C*16        IF ' ' THEN NO ERROR.
C
      CHARACTER ERR*(*)
      DIMENSION X(NN),N(NN),XA(KNA)
      NA=1
      XA(1)=X(1)
      IF(NN.LE.1) GOTO 3
C LINEAR
      IF(IJOB.LE.0) THEN
        DO 1 I=1,NN-1
        IF(N(I).LE.0) GOTO 3
        DX=(X(I+1)-X(I))/N(I)
        X1=X(I)
        DO 1 J=1,N(I)
        X1=X1+DX
        NA=NA+1
        IF(NA.GT.KNA) GOTO 3
    1   XA(NA)=X1
C LOG
      ELSE
        IF(X(1).LE.0.0) GOTO 3
        DO 4 I=1,NN-1
        IF(N(I).LE.0) GOTO 3
        IF(X(I+1).LE.0.0) GOTO 3
        DX=EXP(LOG(X(I+1)/X(I))/N(I))
        X1=X(I)
        DO 4 J=1,N(I)
        X1=X1*DX
        NA=NA+1
        IF(NA.GT.KNA) GOTO 3
    4   XA(NA)=X1
      ENDIF
C
      ERR=' '
      RETURN
    3 ERR='ERROR IN XSET'
      RETURN
      END
      SUBROUTINE UTLSPC(A,IL)
C
C   Get first string delimited by space.
C
C --- History
C   2000.11.22 Created by M.Yamano
C
C --- Input
C A     C       string(maximum length of A = 256)
C
C --- Output
C A     C       (non-space part + spaces) string
C IL    I       length of non-space part in string A
C ----------------------------------------------------------------------
      CHARACTER A*(*),B*256
C
      NL=LEN(A)
      IF (NL.GT.256) NL=256
      N=0
      IL=0
      DO WHILE((N.LT.NL).AND.(IL.GE.0))
        N=N+1
        IF (A(N:N).NE.' ') THEN
           IL=IL+1
           B(IL:IL)=A(N:N)
        ELSE
           IL=-IL
        ENDIF
      ENDDO
      IL=ABS(IL)
      IF (IL.GT.0) A=B(1:IL)
      RETURN
      END
C
      SUBROUTINE MMIES(INIT,XMAX,IP,CR,CI,X,NANG,THETA,KNA1,Q,ERR)
C FOR ONE PARTICLE
C EFFICIENCY FACTOR FOR EXTINCTION, SCATTERING AND ITS ANGULAR
C  DISTRIBUTIONS. MIE THEORY.
C  MAC VERSION WITH -AIMAG-
C--- HISTORY
C 87.11.12  CREATED FROM MIES
C 89. 8.10  RE-ORGANIZED WITH NEW CSJBES/CSNBES (MAC)
C 02. 3. 1  Modified by Yamano.(Qabs(Ci=0) is forced to be zero.)
C--- INPUT
C INIT    I         IF 1 THEN CALCULATE ANGULAR DEPENDENT PARTS.
C XMAX    R         MAXIMUM SIZE PARAMETER FOR SERIES OF CALCULATION.
C IP      I         IP=0  QEXT, QABS
C                   IP=2  QEXT, QABS, Q1, Q2
C                   IP=4  QEXT, QABS, Q1, Q2, Q3, Q4.
C CR      I         REFRACTIVE INDEX = CR + CI I.
C CI      I            CI: GENERALLY NEGATIVE.
C X       R         SIZE PARAMETER.
C NANG    I         NUMBER OF SCATTERING ANGLES.
C THETA R(NANG)     SCATTERING ANGLES IN DEGREE,  I=1, NANG.
C KA1     I         SIZE OF FIRST ARGUMENT OF ARRAY Q   .LE. NANG+1.
C--- OUTPUT
C INIT    I         0
C Q   R(KNA1,*)     QLIJ  GENERALIZED EFFICIENCY FACTORS
C                    =I(L,THETA(I))/(PI*X(J)**2)
C                   IF I.LE.NANG THEN
C                     I(L=1, THETA(I))=ABS(S1(THETA(I)))**2
C                     I(L=2, THETA(I))=ABS(S2(THETA(I)))**2
C                     I(L=3, THETA(I))=REAL(CONJG(S1)*S2)
C                     I(L=4, THETA(I))=IMAG(CONJG(S1)*S2),
C                   IF I=NANG+1 THEN
C                     EXTINCTION EFFICIECY  IF L=1,
C                     ABSORPTION EFFICIECY  IF L=2.
C ERR    C*16       ERROR INDICATER. IF ' ' THEN NORMAL
C--- PARAMETER
C KNANG     I       ARGUMENT SIZE CORRESPONDING TO NANG.
C NMMAX     I       1.2*XMAX+5.
C NMMAX2    I       1.2*REAL(CM)*XMAX+20.
C
C AREA FOR THIS ROUTINE
      SAVE PAI,TAU
      COMPLEX  A,B,S1,S2,CM,S12
      CHARACTER ERR*(*)
      DIMENSION THETA(*),Q(KNA1,*)
      PARAMETER (NMMAX =1500)
      PARAMETER (NMMAX2=2500)
      PARAMETER (KNANG =200)
      PARAMETER (PI=3.14159265358979, R=PI/180.0)
C WORKING AREAS
      DIMENSION PAI(NMMAX,KNANG),TAU(NMMAX,KNANG),A(NMMAX),B(NMMAX)
C
      ERR=' '
C CONSTANTS
      CM = CMPLX( CR, CI)
      NANG3=(NANG+1)/2
      NANG1=NANG-1
      NANG2=NANG-2
      NMX=1.2*XMAX+5
      IF(NMX.GT.NMMAX) THEN
        ERR='ILL NMX'
        RETURN
      ENDIF
C
      IF(IP.GT.0 .AND. INIT.GT.0) THEN
        INIT=0
        DO 1 NA=1,NANG
        SC=COS(THETA(NA)*R)
        PAI(1,NA)=1
        PAI(2,NA)=3*SC
        TAU(1,NA)=SC
        TAU(2,NA)=6*SC*SC-3
        DO 1 N=3,NMX
        PAI(N,NA)=REAL(N+N-1)/REAL(N-1)*SC*PAI(N-1,NA)
     &         -REAL(N)    /REAL(N-1)   *PAI(N-2,NA)
    1   TAU(N,NA)=REAL(N)*SC*PAI(N,NA)-REAL(N+1)*PAI(N-1,NA)
      ENDIF
C
      SUME=0
      SUMS=0
      XX=X
      CALL SIZCO(XX,CM,M3,A,B,ERR)
      IF(ERR.NE.' ') RETURN
      QE1=0
      QS1=0
      DO 3 N=1,M3
      CO=REAL(N+N+1)
      QE1=QE1+CO*REAL(A(N)+B(N))
    3 QS1=QS1+CO*(ABS(A(N))**2+ABS(B(N))**2)
      CO=2.0/(XX*XX)
      SUME=SUME+CO*QE1
      SUMS=SUMS+CO*QS1
      IF(IP.GT.0) THEN
        DO 2 J=1,IP
        DO 2 NA=1,NANG
    2   Q(NA,J)=0
        CO=PI*XX*XX
        DO 4 NA=1,NANG
        S1=(0.0, 0.0)
        S2=(0.0, 0.0)
        DO 5 N=1,M3
        CO3=REAL(N+N+1)/REAL(N*N+N)
        S1=S1+CO3*(A(N)*PAI(N,NA)+B(N)*TAU(N,NA))
    5   S2=S2+CO3*(B(N)*PAI(N,NA)+A(N)*TAU(N,NA))
        Q(NA,1)=Q(NA,1)+ABS(S1)**2/CO
        Q(NA,2)=Q(NA,2)+ABS(S2)**2/CO
        IF(IP.GT.2) THEN
          S12=CONJG(S1)*S2
          Q(NA,3)=Q(NA,3)+REAL(S12)/CO
          Q(NA,4)=Q(NA,4)+AIMAG(S12)/CO
        ENDIF
    4   CONTINUE
      ENDIF
      Q(NANG+1,1)=SUME
c
c --- 02/03/01 by Yamano.(Qabs(Ci=0) is forced to be zero.) --
c
c     Q(NANG+1,2)=SUME-SUMS
      IF (ABS(SUME/SUMS-1).LE.1.0E-7) THEN
         Q(NANG+1,2)=0.
      ELSE
         Q(NANG+1,2)=SUME-SUMS
      ENDIF
C ------------------------------------------------------------
      RETURN
      END
      SUBROUTINE SIZCO(X,CM,M3,A,B,ERR)
C CALCULATION OF  A&B FOR A GIVEN SIZE PARAMETER.
C--- HISTORY
C 89. 8.10  REGISTERED
C--- INPUT
C X         R          SIZE PARAMETER
C CM        CR         COMPLEX REFRACTIVE INDEX
C--- OUTPUT
C M3        I          MAX. ORDER OF A(N) AND B(N)
C                        M3=1.2*X+5
C A       CR(M3)       A(N), I=1,M3
C B       CR(M3)       B(N), I=1,M3
C ERR      C*10        ERROR INDICATER
C--- PARAMETER
C NMMAX     I          1.2*XMAX+5
C NMMAX2    I          1.2*REAL(CM)*XMAX+20
C
      COMPLEX CM,DY,Y,GX,A,B,ZETA,RT
C AREA FOR THIS SUBROUTINE
      CHARACTER ERR*(*)
      DIMENSION A(*),B(*)
      PARAMETER (NMMAX =1500)
      PARAMETER (NMMAX2=2500)
      PARAMETER (NMMAX1=NMMAX+15)
C WORKING AREAS
      DIMENSION SJBES(NMMAX),SNBES(NMMAX),DX(NMMAX1),DY(NMMAX2)
C
      ERR=' '
      IF(X .LE. 1.0) THEN
        M1=1.2*X+7
      ELSE
        IF(X .LE. 20.0) THEN
          M1=1.2*X+10
        ELSE
          IF(X .LE. 50.0) THEN
            M1=1.2*X+12
          ELSE
            M1=1.2*X+18
          ENDIF
        ENDIF
      ENDIF
C
      XM=X*REAL(CM)
      IF(XM .LE. 1.0) THEN
        M2=1.2*XM+7
      ELSE
        IF(XM .LE. 20.0) THEN
          M2=1.2*XM+10
        ELSE
          IF(XM .LE. 50.0) THEN
            M2=1.2*XM+12
          ELSE
            M2=1.2*XM+18
          ENDIF
        ENDIF
      ENDIF
C
      DX(M1+1)=0
      DO 1 N=1,M1
      NN=M1+1-N
    1 DX(NN)=REAL(NN+1)/X-1.0/(DX(NN+1)+REAL(NN+1)/X)
      DY(M2+1)=(0.0,0.0)
      Y=CM*X
      DO 2 N=1,M2
      NN=M2+1-N
      DY(NN)=REAL(NN+1)/Y-1.0/(DY(NN+1)+REAL(NN+1)/Y)
    2 CONTINUE
      GX=(0.0,-1.0)
      M3=1.2*X+5
      CALL CSJBES(M3,X,SJBES,ERR)
      IF(ERR.NE.' ') RETURN
      CALL CSNBES(M3,X,SNBES,ERR)
      IF(ERR.NE.' ') RETURN
      DO 3 N=1,M3
C PSAI=X*SJBES(N)
      PSAI=SJBES(N)
C ZETA=PSAI-(0.0,1.0)*X*SNBES(N)
      ZETA=PSAI-(0.0,1.0) *SNBES(N)
      RT=PSAI/ZETA
      GX=1.0/(REAL(N)/X-GX)-REAL(N)/X
      A(N)=RT*(DY(N)-CM*DX(N))/(DY(N)-CM*GX)
      B(N)=RT*(CM*DY(N)-DX(N))/(CM*DY(N)-GX)
    3 CONTINUE
      RETURN
      END
      SUBROUTINE CSJBES(NMAX,X,SJBES,ERR)
C   SPHERICAL BESSEL FUNCTION OF FIRST KIND
C--- HISTORY
C 89. 8.10  MODIFIED AND REGISTERED (MAC)
C--- INPUT
C NMAX       I        MAX. ORDER
C X          R        INDEPENDENT VARIABLE
C--- OUTPUT
C SJBES    R(NMAX)    SPHERICAL BESSEL FUNCTION OF FIRST KIND
C                      N=1,NMAX
C ERR       C*10      ERROR INDICATER.  IF ' ' THEN NORMAL
C
      CHARACTER ERR*(*)
      DIMENSION SJBES(NMAX)
C
      ERR=' '
      IF(NMAX.GE.30000) THEN
        ERR='LARGE NMAX'
        RETURN
      ENDIF
      IF(X .GE. 3.0E4) THEN
        ERR='LARGE X'
        RETURN
      ENDIF
      IF(X.LT.0.0) THEN
        ERR='ILL X'
        RETURN
      ENDIF
      IF(NMAX.LE.0) THEN
        ERR='ILL NMAX'
        RETURN
      ENDIF
C SMALL X
      IF(X.LE.7.0E-4) THEN
        T1=3
        T2=1
        DO 1 N=1,NMAX
        IF(N.LE.10) THEN
          T3=T2*X/T1
          T1=T1+2
          T2=T3
          SJBES(N)=T3
        ELSE
          SJBES(N)=0
        ENDIF
    1   CONTINUE
        RETURN
      ENDIF
C LARGE X
      IF(X.LT.0.2) THEN
        Y=X*X
        W=X*(1.0-0.1*Y+Y*Y/280.0)/3.0
      ELSE
        W=(SIN(X)-X*COS(X))/(X*X)
      ENDIF
CC
      IF(X .GE. 100.0) THEN
        L=0.02*X+18
      ELSE
        IF(X .GE. 10.0) THEN
          L=0.1*X+10
        ELSE
          IF(X .GT. 1.0) THEN
            L=0.5*X+5
          ELSE
            L=5
          ENDIF
        ENDIF
      ENDIF
CC
      NM=MAX0(NMAX,INT(X))+L-1
      Z=1.0/X
      T3=0
      T2=1.0E-35
      DO 2 I=1,NM
      K=NM-I+1
      T1=REAL(K+K+3)*Z*T2-T3
      IF(K.GT.NMAX) THEN
        IF(ABS(T1) .GE. 1.0E25) THEN
          T1=T1*1.0E-25
          T2=T2*1.0E-25
        ENDIF
      ELSE
        SJBES(K)=T1
        IF(ABS(T1) .GE. 1.0E25) THEN
          T1=T1*1.0E-25
          T2=T2*1.0E-25
          DO 3 J=K,NMAX
    3     SJBES(J)=SJBES(J)*1.0E-25
        ENDIF
      ENDIF
      T3=T2
      T2=T1
    2 CONTINUE
C
      T1=W/T1
      DO 4 N=1,NMAX
    4 SJBES(N)=T1*SJBES(N)
      RETURN
      END
      SUBROUTINE CSNBES(NMAX,X,SNBES,ERR)
C SPHERICAL BESSEL FUNCTION OF SECOND KIND
C--- HISTORY
C 89. 8.10  MODIFIED AND REGISTERED (MAC)
C 92. 3.21  INTRODUCE EMAX
C--- INPUT
C NMAX       I        MAX. ORDER
C X          R        INDEPENDENT VARIABLE
C ERR       C*10      IF ERR.NE.' ' THEN INITIALIZATION.
C--- OUTPUT
C SNBES    R(NMAX)    SPHERICAL BESSEL FUNCTION OF SECOND KIND
C                      N=1,NMAX
C ERR       C*10      ERROR INDICATER.  IF ' ' THEN NORMAL
C
      SAVE EMAX,INIT
      CHARACTER ERR*(*)
      DIMENSION SNBES(NMAX),CP(3)
      DATA INIT/1/
C
      IF(ERR.NE.' ' .OR. INIT.NE.0) THEN
        INIT=0
        CALL CPCON(CP)
        EMAX=10.0**(CP(3)-5)
      ENDIF
      ERR=' '
      IF(NMAX.GE.30000) THEN
        ERR='LARGE NMAX'
        RETURN
      ENDIF
      IF(NMAX.LE.0) THEN
        ERR='ILL NMAX'
        RETURN
      ENDIF
      IF(X.LE.0.0) THEN
        ERR='ILL X'
        RETURN
      ENDIF
C SMALL X
      Z=1.0/X
      IF(X.LE.7.0E-4) THEN
        M=30.0/LOG10(Z)-1
        IF(NMAX.GT.M) THEN
          ERR='SMALL X'
          DO 4 I=1,NMAX
    4     SNBES(I)=-EMAX
          RETURN
        ENDIF
        QN0=Z
        QN1=1
        DO 1 N=1,NMAX
        QN2=QN0*QN1*Z
        QN1=QN1+2
        QN0=QN2
    1   SNBES(N)=-QN2
        RETURN
      ENDIF
C LARGE X
      QN0=-Z*COS(X)
      QN1=Z*(QN0-SIN(X))
      SNBES(1)=QN1
      IF(NMAX.LE.1) RETURN
      DO 2 N=2,NMAX
      QN2=REAL(N+N-1)*Z*QN1-QN0
      IF(QN2.LE.-EMAX) THEN
        ERR='UNDERFLOW'
        DO 3 I=N,NMAX
    3   SNBES(I)=-EMAX
        RETURN
      ENDIF
      QN0=QN1
      QN1=QN2
    2 SNBES(N)=QN2
      RETURN
      END
      SUBROUTINE CPCON(C)
C MACHINE CONSTANTS OF COMPUTER
C--- OUTPUT
C C     R(3)      (1)  MINIMUM POSITIVE X FOR  1+X      .NE. 1
C                      AVERAGE AS A RESULT OF COMPLEX ARITHMETIC
C                      OPERATIONS.
C                 (2)  MINIMUM EXPONENT Y FOR  10.0**Y  .NE. 0
C                 (3)  MAXIMUM EXPONENT Z FOR  10.0**Z  IS MAX. VALUE
C                  IF INIT=1 (DATA STATEMENT) THEN  SET AS Z=Y
C                  IF INIT=2 THEN THIS ROUTINE GETS ACTUAL VALUE OF Z.
C                  - SEE NOTE -
C--- HISTORY
C 90. 1.20  CREATED
C     6.27  CHANGE THE ALGORITHM TO GET X AND Y TAKING INTO ACCOUNT THE
C           HIGH ACCURACY CO-PROCESSOR AND GRACEFUL UNDERFLOW.
C 92. 3.21  N=1000 from N=200
C     7. 9  BUG IN X-DEFINITION
C--- NOTE
C THIS PROGRAM WILL GENERATE -UNDERFLOW ERROR- AND -OVERFLOW ERROR-
C  MESSAGES.  ON SOME COMPUTER -OVERFLOW ERROR- MESSAGE MAY BE
C  FATAL ERROR.  IN THAT CASE, PLEASE SET INIT = 1 IN THE DATA
C  SATEMENT FOR SUPPRESSING THE PROCEDURE OF GETTING C(3).
      DIMENSION C(*)
      CHARACTER CH*80
C RESOLUTION OF COMPUTATION
      SAVE INIT,X,Y,Z
      DATA INIT/1/
      IF(INIT.LE.0) THEN
        C(1)=X
        C(2)=Y
        C(3)=Z
        RETURN
      ENDIF
CC TEST SUM(K=1,M) COS((2K-1)*PI/(2M+1)) = 0.5
CC SIMPLE CHECK, X = X + E, IS NOT VALIDE WHEN THE COMPUTER
CC  USE A HIGH ACCURATE CO-PROCESSOR.
      N=500
      PI=ASIN(1.0)*2
      M0=10
      X=0
      DO 1 M=1,M0
      Y=0
      DO 2 K=1,M
    2 Y=Y+COS((2*K-1)*PI/(2*M+1))
      Y=ABS(2*Y-1)
      X=X+Y
    1 CONTINUE
      X=X/M0
      C(1)=X
C EXPONENT FOR MINIMUM POSITIVE VALUE
C  THIS PROCEDURE WILL GENERATE -UNDERFLOW ERROR MESSAGE-
      Y2=1
      N=1000
      DO 3 I=1,N
      Y1=Y2
      Y3=Y1/10
CC FOR GRACEFUL UNDERFLOW
CC EVEN Y2 BECOMES 0 AS OUTPUT, Y2 IS NOT 0 INSIDE
CC COMPUTER WHEN GRACEFUL UNDERFLOW IS APPLIED.
CC SO WE REPLACE THE VALUE OF Y2 BY OUTPUT.
      CH='0'
      WRITE(CH,7) Y3
    7 FORMAT(1P,E12.5)
      Y2=0
      READ(CH,*) Y2
      IF(ABS(10*Y2/Y1-1) .GT. 5*X) GOTO 4
    3 CONTINUE
      I=N+1
    4 Y=1-I
      C(2)=Y
C EXPONENT FOR MAXIMUM POSITIVE VALUE
C THIS PROCEDURE WILL GENERATE -OVERFLOW MESSAGE-
      IF(INIT.LE.1) THEN
        Z=-Y
       ELSE
        Z2=1
        DO 5 I=1,N
        Z1=Z2
        Z2=Z1*10
        IF(ABS(Z2/Z1/10-1) .GT. 5*X) GOTO 6
    5   CONTINUE
        I=N+1
    6   Z=I-1
      ENDIF
      C(3)=Z
C
      INIT=0
      RETURN
      END
