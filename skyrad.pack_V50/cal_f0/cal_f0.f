C$OUTF cal_f0.f
C$LIBF /Users/makibo/Fortran/LIB/SKYLIB/A
C$LIBF /Users/makibo/Fortran/LIB/POMLIB/A
C$LIBF /Users/makibo/Fortran/LIB/MATLIB/A
C$LIBF /Users/makibo/Fortran/LIB/BASLIB/A
C$LIBF /Users/makibo/Fortran/LIB_Nakajima/LBC1/A
C$LIBF /Users/makibo/Fortran/LIB_Nakajima/LBM1/A
C$LIBF /Users/makibo/Fortran/LIB_Nakajima/LBM2/A
C$LIBF /Users/makibo/Fortran/LIB_Nakajima/LBR1/A
C$LIBF /Users/makibo/Fortran/LIB_Nakajima/LBR2/A
C
C --- Note
C   Skyradiometer data analyzing program
C   For PREDE instruments POM-01L/POM-01MKII
C
C  (1) Determination of calibration constants(F0)
C
C --- History
C   2002.11.15 Created from 'fproc?' by M.Yamano
C   2003.02.25 Added parameters TMAX, WMIN.
C   2003.03.05 Added 'ins.para' replacement option.
C   2003.09.16 Debugged.
C   2004.12.06 Start on Mac.
C   2006.04.05 Added output option IPLT.
C
C --- I/O files
C  cal_f0.par          I    processing control file
C  fname               I    processing dates file
C  F0d/(yy)yymmdd.w??  I    measured datafile
C  f0_w??.out          O    results files for calculation of F0
C  F0.out              O    Determined F0 values
C  F0Nd.w??            O   (normal)  Langley plot for a data set (IPLT>0)
C  F0Id.w??            O    Improved Langley plot for a data set (IPLT>0)
C  F0Ad.w??            O    Improved Langley plot for all data sets(IPLT>0)
C ----------------------------------------------------------------------
      PARAMETER (KNN   =3000)
      PARAMETER (KNDY  =500)
      PARAMETER (KNW   =20)
      PARAMETER (KDAY  =30)
C
      PARAMETER (PI=3.141592653589793,RAD=PI/180.0)
C
      CHARACTER ERC*64
      CHARACTER FLI*80,FLD*80,FLO*80,LB*2,ATRF*4,ATR*4,DIR*4
      CHARACTER PFIL*80,FLI0*80
C
      DIMENSION IYD(KNN),IMD(KNN),IDD(KNN),TMD(KNN),EMS(KNN)
     &         ,FOBS(KNN),TAR(KNN),TAO(KNN),TAUE(KNN),TAUS(KNN)
     &         ,F0N(KNDY),SGN(KNDY),TTN(KNDY),NON(KNDY)
     &         ,F0I(KNDY),SGI(KNDY),W0I(KNDY),NOI(KNDY)
     &         ,XXX(KNDY),IXXX(KNDY),IDDAY(KNDY),IAMPM(KNDY)
      DIMENSION WVDAT(KNW),F0DAT(3,KNW)
C
      DATA NMIN,SGMX,TMIN,TMAX,WMIN,WMAX /10, 0.1, 0.01, 0.3, 0., 2./
      DATA DIR,ATR /'F0d/','.w00'/
C
      ERC=' '
      IUI=1
      OPEN (IUI,FILE='cal_f0.par',STATUS='OLD')
      CALL RDCLFP(IUI,IPLT,IRPL,PFIL,NFL,ERC)
      CLOSE (IUI)
      IF (ERC.NE.' ') GOTO 900
C
      NFNM0=0
      IUF=3
      IUO=9
      IF (IPLT.GT.0) THEN
         IUP=10
      ELSE
         IUP=0
      ENDIF
      JEND=0
      NW=0
      IW=0
      DO WHILE ((IW.LT.KNW).AND.(JEND.EQ.0))
        IW=IW+1
        WRITE(LB,5) IW
    5   FORMAT(I2)
        IF (IW.LT.10) THEN
           IST=2
        ELSE
           IST=1
        ENDIF
        ATRF=ATR(1:1+IST)//LB(IST:2)
        WL00=-999.
        NDT=0
        OPEN (IUF,FILE='fname',STATUS='OLD')
C
   10   READ(IUF,15,END=100,ERR=100) FLI
   15   FORMAT(A80)
        CALL GETPOS(FLI,'.',80,IL)
        IF (IL.EQ.1) GOTO 100
        IF (IL.GT.0) FLI=FLI(1:IL-1)
        CALL UTLSPC(FLI,NFNM)
        FLD=DIR//FLI(1:NFNM)//ATRF
        IF ((IRPL.GT.0).AND.(NFNM0.EQ.0)) THEN
           NFNM0=NFNM
           FLI0=FLI
        ENDIF
        NTT=NDT
        OPEN (IUI,FILE=FLD,STATUS='OLD',ERR=50)
        CALL RDDATF(IUI,WL0,JWV,NDT,IYD,IMD,IDD,TMD,EMS
     &             ,FOBS,TAR,TAO,TAUE,TAUS)
   50   CLOSE (IUI)
        IF (NTT.EQ.0) WL00=WL0
        IF (ABS(WL0-WL00).GT.1.E-9) THEN
           ERC='Mismatch of input data files !'
           GOTO 900
        ENDIF
        IF (NDT.LT.KNN)  GOTO 10
C
  100   IF (NDT.GT.0) THEN
           FLO='f0_'//ATRF(2:4)//'.out'
           OPEN (IUO,FILE=FLO,STATUS='UNKNOWN')
           IF (IPLT.GT.0) THEN
              OPEN (IUP  ,FILE='F0Nd'//ATRF,STATUS='UNKNOWN')
              OPEN (IUP+1,FILE='F0Id'//ATRF,STATUS='UNKNOWN')
              OPEN (IUP+2,FILE='F0Ad'//ATRF,STATUS='UNKNOWN')
           ENDIF
           CALL PROCF0(IUP,NDT,JWV,IYD,IMD,IDD,EMS,FOBS,TAR,TAO,TAUS
     &          ,NDAY,IDDAY,IAMPM,F0N,SGN,TTN,NON,F0I,SGI,W0I,NOI
     &          ,F0II,SGII,W0II,NOII,ERC)
           IF (ERC.NE.' ') GOTO 190
C
           CALL F0MEAN(NDAY,F0N,SGN,TTN,NON,F0I,SGI,W0I,NOI,NN,FN00,SGN0
     &                ,NI,FI00,SGI0,NMIN,SGMX,TMIN,TMAX,WMIN,WMAX)
           CALL WTCALF(IUO,WL0,NDAY,IYD,IMD,IDD,IAMPM,IDDAY
     &                ,F0N,SGN,TTN,NON,F0I,SGI,W0I,NOI
     &                ,F0II,SGII,W0II,NOII,NN,FN00,SGN0
     &                ,NI,FI00,SGI0,NMIN,SGMX,TMIN,TMAX,WMIN,WMAX)
           DO I=1,NDAY
             IF (TTN(I).GT.0.0) THEN
                XXX(I)=SGN(I)*TTN(I)
             ELSE
                XXX(I)=999
             ENDIF
           ENDDO
           CALL SORT1(NDAY,XXX,IXXX)
           CALL WTSRTF(IUO,WL0,NDAY,IYD,IMD,IDD,IAMPM,IDDAY,IXXX
     &                      ,F0N,SGN,TTN,NON,F0I,SGI,W0I,NOI)
           NW=NW+1
           WVDAT(NW)=WL00
           F0DAT(1,NW)=FN00
           F0DAT(2,NW)=FI00
           F0DAT(3,NW)=F0II
  190      CLOSE (IUO)
           CLOSE (IUP)
           CLOSE (IUP+1)
           CLOSE (IUP+2)
        ELSE
           JEND=1
        ENDIF
        CLOSE (IUF)
      ENDDO
      IF (NW.GT.0) THEN
         OPEN (IUO,FILE='F0.out',STATUS='UNKNOWN')
         WRITE(IUO,200) (WVDAT(I),I=1,NW)
  200    FORMAT('WL :',1P7E10.3)
         WRITE(IUO,210) (F0DAT(1,I),I=1,NW)
  210    FORMAT('F0N:',1P7E10.3)
         WRITE(IUO,220) (F0DAT(2,I),I=1,NW)
  220    FORMAT('F0I:',1P7E10.3)
         WRITE(IUO,230) (F0DAT(3,I),I=1,NW)
  230    FORMAT('F0A:',1P7E10.3)
         CLOSE (IUO)
         IF (IRPL.GT.0) CALL RPINSP(IUO,IRPL,PFIL,NFL,FLI0,NFNM0
     &                             ,NW,WVDAT,F0DAT,ERC)
      ENDIF
C
  900 IF (ERC.NE.' ') print*, ERC
      STOP
      END
C
      SUBROUTINE RPINSP(IU,IRPL,PFIL,NFL,FLI,NFNM,NWI,WVDAT,F0DAT,ERC)
      PARAMETER (KNW   =20)
      PARAMETER (KDAY  =30)
C
      CHARACTER ERC*(*)
      CHARACTER SRNO*20
      CHARACTER PFIL*80,FLI*80
      DIMENSION WL(KNW),SOLID(KNW),F0D(KNW,KDAY),TJLC(KDAY)
      DIMENSION WVDAT(KNW),F0DAT(3,KNW)
C
      OPEN (IU,FILE=PFIL(1:NFL),STATUS='OLD')
      CALL RDINSP(IU,SRNO,ITYP,NW,WL,SOLID,NDY,F0D,TJLC,ERC)
      CLOSE (IU)
      IF (ERC.NE.' ') GOTO 900
C
      IF (NFNM.EQ.6) THEN
         READ(FLI,10) IY,IM,ID
   10    FORMAT(3I2)
         IF (IY.GT.95) THEN
            IY=IY+1900
         ELSE
            IY=IY+2000
         ENDIF
      ELSE
         READ(FLI,20) IY,IM,ID
   20    FORMAT(I4,2I2)
      ENDIF
      TM=0.0
      NDY=1
      DO IW=1,NW
        IP=0
        I=0
        DO WHILE((I.LT.NWI).AND.(IP.EQ.0))
          I=I+1
          IF (ABS(WVDAT(I)-WL(IW)).LT.0.001E-4) IP=I
        ENDDO
        IF (IP.GT.0) THEN
           F0D(IW,1)=F0DAT(IRPL,IP)
        ELSE
           F0D(IW,1)=-9.999E-4
        ENDIF
      ENDDO
C
      OPEN (IU,FILE=PFIL(1:NFL),STATUS='UNKNOWN')
      WRITE(IU,50) SRNO,ITYP,NW
   50 FORMAT('- instrument parameters -'/'"',A20,'"  : instrument S/N'
     &   /I2,'  : instrument type (POM-01L:10,11/ POM-01MKII:20,21,22)'
     &   /I2,'  : NW(number of wavelengths)/ WL[cm]/ SVA[sr]')
      WRITE(IU,60) (WL(I),I=1,NW)
      WRITE(IU,60) (SOLID(I),I=1,NW)
   60 FORMAT(1P20E11.3)
      WRITE(IU,70) NDY
   70 FORMAT(I2,'  : NDAY(number of calib. constant data)/ '
     &         ,'YYYY MM DD HR(GMT) F0(IW=1,NW)')
      WRITE(IU,80) IY,IM,ID,TM,(F0D(IW,1),IW=1,NW)
   80 FORMAT(I4,2I3,0PF5.1,1P20E11.3)
      CLOSE (IU)
C
  900 RETURN
      END
C
      SUBROUTINE RDDATF(IU,WL0,JWV,N,IYD,IMD,IDD,TMD,EMS
     &                        ,FOBS,TAR,TAO,TAUE,TAUS)
C
C   Read data from 'f0_wl?.dat' for determination F0. (old version)
C   Read data from '(yy)yymmdd.w??' for determination F0.
C
C --- history
C   2000.12.08 Created by M.Yamano
C   2002.11.15 Renewed.
C
C --- Input
C IU          I         device No. for reading
C N           I         data No.
C
C --- Output
C N           I         data No.
C ----------------------------------------------------------------------
      PARAMETER (KNN   =3000)
C
      CHARACTER CH*1
      DIMENSION IYD(KNN),IMD(KNN),IDD(KNN),TMD(KNN),EMS(KNN)
     &         ,FOBS(KNN),TAR(KNN),TAO(KNN),TAUE(KNN),TAUS(KNN)
C
      READ(IU,*,END=900,ERR=900) WL0,IW,IO3,IWV
      IF (IW.EQ.IWV) THEN
         JWV=1
      ELSE
         JWV=0
      ENDIF
      READ(IU,10,END=900,ERR=900) CH
   10 FORMAT(A1)
C
   20 N=N+1
      READ(IU,*,END=890,ERR=890) IYD(N),IMD(N),IDD(N),TMD(N),EMS(N)
     &                       ,FOBS(N),TAR(N),TAO(N),TAUE(N),TAUS(N)
      IF (N.GE.KNN) GOTO 900
      GOTO 20
C
  890 N=N-1
  900 RETURN
      END
C
      SUBROUTINE PROCF0(IUP,NN,JWV,IYD,IMD,IDD,EMS,FOBS,TAR,TAO,TAUS
     &                 ,NDY,IDDAY,IAMPM,F0N,SGN,TTN,NON,F0I,SGI,W0I,NOI
     &                 ,F0II,SGII,W0II,NOII,ERC)
C
C   Main processing routine for determination of F0.
C
C --- history
C   2001.01.24 Created by M.Yamano
C   2001.03.14 Added F0II,SGII,W0II,NOII.
C   2006.04.06 Added argument IUP.
C
C --- Input
C
C --- Output
C ERC           C*64      ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNN   =3000)
      PARAMETER (KNDY  =500)
C
      CHARACTER ERC*(*)
      CHARACTER HD*11,AP*2
      DIMENSION IYD(KNN),IMD(KNN),IDD(KNN),EMS(KNN)
     &         ,FOBS(KNN),TAR(KNN),TAO(KNN),TAUS(KNN)
     &         ,F0N(KNDY),SGN(KNDY),TTN(KNDY),NON(KNDY)
     &         ,F0I(KNDY),SGI(KNDY),W0I(KNDY),NOI(KNDY)
     &         ,IDDAY(KNDY),IAMPM(KNDY)
      DIMENSION X(KNN),Y(KNN),YI(KNN),Z(KNN)
     &         ,YII(KNN),ZZ(KNN)
C
      DATA EMMX, NDLIM, SRG/4.5, 4, 0.5/
C
      LDY=0
      N=0
      K=0
      KK=0
      NDY=0
      DO WHILE(N.LT.NN)
        N=N+1
        LDAY=IYD(N)*10000+IMD(N)*100+IDD(N)
        IF (LDY.EQ.0) THEN
           WRITE(HD,10) IYD(N),IMD(N),IDD(N)
   10      FORMAT(I5,2I3)
           LDY=LDAY
        ENDIF
        IF (LDAY.EQ.LDY) THEN
           IF ((FOBS(N).GT.0.).AND.(EMS(N).LE.EMMX)) THEN
              K=K+1
              KK=KK+1
              IF (KK.GT.KNN) THEN
                 ERC='PROCF0: KK.GT.KNN'
                 GOTO 900
              ENDIF
              IF (JWV.EQ.0) THEN
                 X(K)=EMS(N)
              ELSE
                 X(K)=SQRT(EMS(N))
              ENDIF
              Y(K)=LOG(FOBS(N))
              YI(K)=LOG(FOBS(N))+EMS(N)*(TAR(N)+TAO(N))
              IF (TAUS(N).GT.0.) THEN
                 Z(K)=EMS(N)*TAUS(N)
              ELSE
                 Z(K)=0.
              ENDIF
              ZZ(KK)=Z(K)
              YII(KK)=YI(K)
           ENDIF
        ELSE
           N=N-1
           IF (K.GE.NDLIM) CALL CALCDY(IUP,HD,N,K,NDLIM,SRG,
     &                                 NDY,X,Y,YI,Z,IDDAY,IAMPM,
     &                          F0N,SGN,TTN,NON,F0I,SGI,W0I,NOI)
           LDY=0
           K=0
           IF (N.GT.NN-NDLIM) N=NN
        ENDIF
      ENDDO
      IF (K.GE.NDLIM) CALL CALCDY(IUP,HD,N,K,NDLIM,SRG,
     &                            NDY,X,Y,YI,Z,IDDAY,IAMPM,
     &                     F0N,SGN,TTN,NON,F0I,SGI,W0I,NOI)
C
      IF (IUP.GT.0) THEN
         IUPP=IUP+2
      ELSE
         IUPP=IUP
      ENDIF
      HD=' yyyy mm dd'
      AP=' *'
      CALL CALCFI(IUPP,HD,AP,1,KK,NDLIM,SRG,ZZ,YII,NOII,F0II,SGII,W0II)
  900 RETURN
      END
C
      SUBROUTINE CALCDY(IUP,HD,ND,KD,NDLIM,SRG,NDY,X,Y,YI,Z
     &                 ,IDDAY,IAMPM,F0N,SGN,TTN,NON,F0I,SGI,W0I,NOI)
      PARAMETER (KNN   =3000)
      PARAMETER (KNDY  =500)
C
      CHARACTER HD*11,AP*2
      DIMENSION X(KNN),Y(KNN),YI(KNN),Z(KNN)
     &         ,F0N(KNDY),SGN(KNDY),TTN(KNDY),NON(KNDY)
     &         ,F0I(KNDY),SGI(KNDY),W0I(KNDY),NOI(KNDY)
     &         ,IDDAY(KNDY),IAMPM(KNDY)
C
      IF (IUP.GT.0) THEN
         IUPI=IUP+1
      ELSE
         IUPI=IUP
      ENDIF
      K=1
      DO WHILE((X(K+1).LT.X(K)).AND.(K.LT.KD-1))
        K=K+1
      ENDDO
      IF (K.GE.NDLIM) THEN
         IF (ABS(X(1)-X(K)).GE.SRG) THEN
            NDY=NDY+1
            IDDAY(NDY)=ND
            IAMPM(NDY)=1
            AP=' 1'
            CALL CALCFN(IUP ,HD,AP,1,K,NDLIM,SRG,X,Y,N1,F01,SG1,TT1)
            CALL CALCFI(IUPI,HD,AP,1,K,NDLIM,SRG,Z,YI,N2,F02,SG2,W02)
            NON(NDY)=N1
            F0N(NDY)=F01
            SGN(NDY)=SG1
            TTN(NDY)=TT1
            NOI(NDY)=N2
            F0I(NDY)=F02
            SGI(NDY)=SG2
            W0I(NDY)=W02
         ENDIF
      ENDIF
      IF (KD-K+1.GE.NDLIM) THEN
         IF (ABS(X(K)-X(KD)).GE.SRG) THEN
            NDY=NDY+1
            IDDAY(NDY)=ND
            IAMPM(NDY)=2
            AP=' 2'
            CALL CALCFN(IUP ,HD,AP,K,KD,NDLIM,SRG,X,Y,N1,F01,SG1,TT1)
            CALL CALCFI(IUPI,HD,AP,K,KD,NDLIM,SRG,Z,YI,N2,F02,SG2,W02)
            NON(NDY)=N1
            F0N(NDY)=F01
            SGN(NDY)=SG1
            TTN(NDY)=TT1
            NOI(NDY)=N2
            F0I(NDY)=F02
            SGI(NDY)=SG2
            W0I(NDY)=W02
         ENDIF
      ENDIF
  900 RETURN
      END
C
      SUBROUTINE CALCFN(IU,HD,AP,K1,K2,NDLIM,SRG,X,Y,NON,F0N,SGN,TTN)
      PARAMETER (KNN   =3000)
C
      CHARACTER HD*11,AP*2
      DIMENSION X(KNN),Y(KNN),XX(KNN),YY(KNN),IFG(KNN)
C
      DATA GMN /0.99/
C
      N=0
      DO I=K1,K2
        N=N+1
        IFG(N)=1
        XX(N)=X(I)
        YY(N)=Y(I)
      ENDDO
      CALL CLLNGE(IU,HD,AP,N,XX,YY,IFG,NON,F0N,SGN1,SGN2,TTN,GMN,SRG)
      SGN=SGN1
      RETURN
      END
C
      SUBROUTINE CALCFI(IU,HD,AP,K1,K2,NDLIM,SRG,Z,YI,NOI,F0I,SGI,W0I)
      PARAMETER (KNN   =3000)
C
      CHARACTER HD*11,AP*2
      DIMENSION YI(KNN),Z(KNN),YY(KNN),ZZ(KNN),IFG(KNN)
C
      DATA GMI /0.95/
C
      N=0
      DO I=K1,K2
        IF (Z(I).GT.0.) THEN
           N=N+1
           IFG(N)=1
           ZZ(N)=Z(I)
           YY(N)=YI(I)
        ENDIF
      ENDDO
      IF (N.GE.NDLIM) THEN
         CALL CLLNGE(IU,HD,AP,N,ZZ,YY,IFG,NOI,F0I,SGI1,SGI2,TTI,GMI,SRG)
         SGI=SGI1
         IF (TTI.LT.99) THEN
            W0I=1/TTI
         ELSE
            W0I=999
         ENDIF
      ELSE
         NOI=N
         F0I=999
         SGI=999
         W0I=999
      ENDIF
      RETURN
      END
C
      SUBROUTINE CLLNGE(IU,HD,AP,NN,XX,YY,IFG,NE,F0,SG,SGR,TT,GMM,SRG)
      PARAMETER (KNN   =3000)
C
      CHARACTER HD*11,AP*2
      DIMENSION XX(KNN),YY(KNN),IFG(KNN),JFG(KNN)
     &         ,X(KNN),Y(KNN),CC(3),WK1(3),WK2(3,3)
      DATA RAT,MLP,CO /0.5, 10, 1./
C
      DO I=1,NN
        JFG(I)=1
      ENDDO
      G=0.
      RT=1
      IRP=1
      LP=-1
      DO WHILE((RT.GE.RAT).AND.(ABS(G).LT.GMM)
     &                    .AND.(LP.LT.MLP).AND.(IRP.EQ.1))
        LP=LP+1
        DO I=1,NN
          IFG(I)=JFG(I)
        ENDDO
        NE=0
        IMX=0
        IMN=0
        DO I=1,NN
          IF (IFG(I).GT.0) THEN
             NE=NE+1
             X(NE)=XX(I)
             Y(NE)=YY(I)
             IF (NE.EQ.1) THEN
                IMX=1
                IMN=1
             ENDIF
             IF (X(NE).GT.X(IMX)) IMX=NE
             IF (X(NE).LT.X(IMN)) IMN=NE
          ENDIF
        ENDDO
        CALL REGRES(NE,X,Y,A,B,G)
        XRG=ABS(LOG(X(IMX))-LOG(X(IMN)))
        IF (XRG.GE.SRG) THEN
           CALL PLYFIT(2,NE,X,Y,CC,3,WK1,WK2,IER)
           IF (IER.EQ.0) THEN
              SG1=0
              DO K=1,NE
                SG1=SG1+(Y(K)-(CC(1)+CC(2)*X(K)))**2
              ENDDO
              CALL CLMEAN(NE,X,XA,STD,SGX)
              F0=EXP(CC(1))
              SG=SQRT(SG1/NE)
              SGR=SG/SGX
              TT=-CC(2)
C
              NJ=0
              DO I=1,NN
                DY=ABS(YY(I)-(CC(1)+CC(2)*XX(I)))
                IF (DY.LE.CO*SG) THEN
                   NJ=NJ+1
                   JFG(I)=1
                ELSE
                   JFG(I)=0
                ENDIF
              ENDDO
              N=0
              IRP=0
              DO WHILE((N.LT.NN).AND.(IRP.EQ.0))
                N=N+1
                IF (JFG(N).NE.IFG(N)) IRP=1
              ENDDO
           ELSE
              NJ=0
              F0=999
              SG=999
              SGR=999
              TT=999
           ENDIF
        ELSE
           NJ=0
           F0=999
           SG=999
           SGR=999
           TT=999
        ENDIF
        RT=FLOAT(NJ)/NN
c       print*,LP,NE,NN,G,RT,XRG
      ENDDO
      IF (ABS(G).LT.GMM) NE=-NE
      IF (IU.GT.0) THEN
         WRITE(IU,10) HD,AP,NE,NN,ABS(G)
   10    FORMAT(A11,A2,2I5,0PF10.4)
         IF (ABS(NE).GT.0) THEN
            DO K=1,ABS(NE)
              WRITE(IU,20) K,X(K),Y(K)
   20         FORMAT(I5,0P2F10.4)
            ENDDO
         ENDIF
      ENDIF
      RETURN
      END
C
      SUBROUTINE CLLNG0(N,XX,YY,F0,SG,SGR,TT,SRG)
      PARAMETER (KNN   =3000)
C
      DIMENSION XX(KNN),YY(KNN),CC(3),WK1(3),WK2(3,3)
C
      IMN=1
      IMX=1
      DO I=1,N
        IF (XX(I).LT.XX(IMN)) IMN=I
        IF (XX(I).GT.XX(IMX)) IMX=I
      ENDDO
      XRG=ABS(LOG(XX(IMX))-LOG(XX(IMN)))
      IF (XRG.GE.SRG) THEN
         CALL PLYFIT(2,N,XX,YY,CC,3,WK1,WK2,IER)
         IF (IER.EQ.0) THEN
            SG1=0
            DO K=1,N
              SG1=SG1+(YY(K)-(CC(1)+CC(2)*XX(K)))**2
            ENDDO
            CALL CLMEAN(N,XX,XA,STD,SGX)
            F0=EXP(CC(1))
            SG=SQRT(SG1/N)
            SGR=SG/SGX
            TT=-CC(2)
         ELSE
            F0=999
            SG=999
            SGR=999
            TT=999
         ENDIF
      ELSE
         F0=999
         SG=999
         SGR=999
         TT=999
      ENDIF
      RETURN
      END
C
      SUBROUTINE F0MEAN(NDAY,F0N,SGN,TTN,NON,F0I,SGI,W0I,NOI,NN,FN00
     &              ,SGN0,NI,FI00,SGI0,NMIN,SGMX,TMIN,TMAX,WMIN,WMAX)
C
C   Get averaged F0.
C
C --- history
C   2001.01.24 Created by M.Yamano
C   2003.02.25 Added arguments TMAX, WMIN.
C
C --- Input
C
C --- Output
C ----------------------------------------------------------------------
      PARAMETER (KNDY  =500)
C
      DIMENSION F0N(KNDY),SGN(KNDY),TTN(KNDY),NON(KNDY)
     &         ,F0I(KNDY),SGI(KNDY),W0I(KNDY),NOI(KNDY)
      DIMENSION DD(KNDY)
C
      NN=0
      DO I=1,NDAY
        IF ((NON(I).GE.NMIN).AND.(SGN(I).LT.SGMX)
     &             .AND.(TTN(I).GT.TMIN).AND.(TTN(I).LT.TMAX)) THEN
           NN=NN+1
           DD(NN)=F0N(I)
        ENDIF
      ENDDO
      IF (NN.GT.0) THEN
         CALL CLMEAN(NN,DD,FN00,STD,SGN0)
         SGN0=SGN0/FN00
      ELSE
         FN00=-99.
         SGN0=-99.
      ENDIF
C
      NI=0
      DO I=1,NDAY
        IF ((NOI(I).GE.NMIN).AND.(SGI(I).LT.SGMX)
     &             .AND.(W0I(I).GT.WMIN).AND.(W0I(I).LT.WMAX)) THEN
           NI=NI+1
           DD(NI)=F0I(I)
        ENDIF
      ENDDO
      IF (NI.GT.0) THEN
         CALL CLMEAN(NI,DD,FI00,STD,SGI0)
         SGI0=SGI0/FI00
      ELSE
         FI00=-99.
         SGI0=-99.
      ENDIF
      RETURN
      END
C
      SUBROUTINE WTCALF(IU,WL0,NDAY,IYD,IMD,IDD,IAMPM,IDDAY
     &                 ,F0N,SGN,TTN,NON,F0I,SGI,W0I,NOI
     &                 ,F0II,SGII,W0II,NOII,NN,FN00,SGN0
     &                 ,NI,FI00,SGI0,NMIN,SGMX,TMIN,TMAX,WMIN,WMAX)
C
C   Write results for determination of F0.
C
C --- history
C   2000.12.08 Created by M.Yamano
C   2001.03.14 Added F0II,SGII,W0II,NOII.
C   2002.11.15 Modified.
C   2003.02.25 Added arguments TMAX, WMIN.
C
C --- Input
C
C --- Output
C ----------------------------------------------------------------------
      PARAMETER (KNN   =3000)
      PARAMETER (KNDY  =500)
C
      DIMENSION IYD(KNN),IMD(KNN),IDD(KNN),IDDAY(KNDY),IAMPM(KNDY)
     &         ,F0N(KNDY),SGN(KNDY),TTN(KNDY),NON(KNDY)
     &         ,F0I(KNDY),SGI(KNDY),W0I(KNDY),NOI(KNDY)
C
      WRITE(IU,10) WL0,NDAY
   10 FORMAT(/1PE12.3,I5,' : WL  NDAY')
      WRITE(IU,12) NMIN,TMAX,WMAX,SGMX,TMIN,WMIN
   12 FORMAT(/' Averaging condition : '
     &           ,'NMIN =',I6,'    TMAX =',0PF6.3,'    WMAX =',F5.1
     &       /23X,'SGMX =',F6.3,'    TMIN =',F6.3,'    WMIN =',F5.1)
      WRITE(IU,15) FN00,SGN0,NN
   15 FORMAT(/'   F0N(averaged) =',1PE11.3,'   SGN/F0N ='
     &                            ,0PF8.3,' ( NN =',I3,')')
      WRITE(IU,16) FI00,SGI0,NI
   16 FORMAT( '   F0I(averaged) =',1PE11.3,'   SGI/F0I ='
     &                            ,0PF8.3,' ( NI =',I3,')')
      WRITE(IU,20)
   20 FORMAT(/'  Y   M  D AP    F0N       SGN      TTN    N'
     &       ,'      F0I       SGI      W0I    N')
      DO I=1,NDAY
        J=IDDAY(I)
        WRITE(IU,30) IYD(J),IMD(J),IDD(J),IAMPM(I)
     &     ,F0N(I),SGN(I),TTN(I),NON(I),F0I(I),SGI(I),W0I(I),NOI(I)
   30   FORMAT(I4,2I3,I2,1X,1P2E10.3,0PF8.3,I4,1X,1P2E10.3,0PF8.3,I4)
      ENDDO
      WRITE(IU,40) F0II,SGII,W0II,NOII
   40 FORMAT(/46X,1P2E10.3,0PF8.3,I5)
      RETURN
      END
C
      SUBROUTINE WTSRTF(IU,WL0,NDAY,IYD,IMD,IDD,IAMPM,IDDAY,IXXX
     &                 ,F0N,SGN,TTN,NON,F0I,SGI,W0I,NOI)
C
C   Write results for determination of F0(after sorting).
C
C --- history
C   2000.12.08 Created by M.Yamano
C   2002.11.15 Modified.
C
C --- Input
C
C --- Output
C ----------------------------------------------------------------------
      PARAMETER (KNN   =3000)
      PARAMETER (KNDY  =500)
C
      DIMENSION IYD(KNN),IMD(KNN),IDD(KNN),IDDAY(KNDY),IAMPM(KNDY)
     &         ,F0N(KNDY),SGN(KNDY),TTN(KNDY),NON(KNDY)
     &         ,F0I(KNDY),SGI(KNDY),W0I(KNDY),NOI(KNDY)
     &         ,IXXX(KNDY)
C
      WRITE(IU,10) WL0,NDAY
   10 FORMAT(/1PE12.3,I5,' : WL  NDAY(SORTED)')
      WRITE(IU,20)
   20 FORMAT(/'  Y   M  D AP    F0N       SGN      TTN    N'
     &       ,'      F0I       SGI      W0I    N')
      DO K=1,NDAY
        I=IXXX(K)
        J=IDDAY(I)
        WRITE(IU,30) IYD(J),IMD(J),IDD(J),IAMPM(I)
     &     ,F0N(I),SGN(I),TTN(I),NON(I),F0I(I),SGI(I),W0I(I),NOI(I)
   30   FORMAT(I4,2I3,I2,1X,1P2E10.3,0PF8.3,I4,1X,1P2E10.3,0PF8.3,I4)
      ENDDO
      RETURN
      END
C
      SUBROUTINE RDCLFP(IU,IPLT,IRPL,PFIL,NFL,ERC)
C
C   Read parameters from file('cal_f0.par')
C
C --- history
C   2003.03.05 Created by M.Yamano
C   2006.04.03 Modified. (If IRPL.GT.3 then IRPL=0)
C   2006.04.05 Added argument IPLT.
C
C --- Input
C IU       I       device No. for reading
C
C --- Output
C IPLT     I       output option for plot files
C IRPL     I       'ins.para' replace option
C PFIL     C*80    'ins.para' file name
C NFL      I       length of 'ins.para' file name -> PFIL(1:NFL)
C ERC      C*64    ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
C
      CHARACTER ERC*(*)
      CHARACTER PFIL*80
C
      ERC=' '
      READ(IU,*,ERR=800,END=800) IPLT
      READ(IU,*,ERR=800,END=800) IRPL
      IF (IRPL.GT.0) THEN
C
C --- read 'ins.para' parameters
C
         READ(IU,*,ERR=800,END=800) PFIL
         CALL UTLSPC(PFIL,NFL)
         IF (NFL.LE.0) THEN
            ERC='RDCLFP: illegal ins.para filename'
            GOTO 900
         ENDIF
         OPEN (IU+1,FILE=PFIL(1:NFL),STATUS='OLD',ERR=850)
         CLOSE (IU+1)
         IF (IRPL.GT.3) IRPL=0
      ELSE
         IRPL=0
      ENDIF
      GOTO 900
C
  800 ERC='RDCLFP: Read Error !'
      GOTO 900
  850 ERC='RDCLFP: No ins.para file !'
C
  900 RETURN
      END
C
      SUBROUTINE GETPOS(A,CH,NL,IL)
C
C   Get the position of a specified character in string.
C
C --- History
C   2002.11.15 Created by M.Yamano
C
C --- Input
C A     C       string (LEN(A)<= 256)
C CH    C*1     the specified character
C NL    I       length of searched part in string A (NL<=256)
C
C --- Output
C IL    I       the position of a specified character in string A
C ----------------------------------------------------------------------
      CHARACTER A*(*),CH*1
C
      NLA=LEN(A)
      NL0=MIN(NLA,NL)
      IF (NL0.GT.256) NL0=256
      IL=0
      N=0
      DO WHILE((N.LT.NL0).AND.(IL.EQ.0))
        N=N+1
        IF (A(N:N).EQ.CH) IL=N
      ENDDO
      RETURN
      END
C
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
      SUBROUTINE SORT1(NX,X,IX)
C SORTING OF X1,X2,... -> Y1<=Y2<=Y3...
C Y(I)=X(IX(I))
C IX: LABEL(ORDER)
C$ENDI
      DIMENSION X(NX),IX(NX)
      DO 1 I=1,NX
    1 IX(I)=I
      IF(NX.LE.1) RETURN
      DO 2 I=1,NX-1
      DO 3 J=I+1,NX
      IF(X(IX(J))-X(IX(I))) 4,3,3
    4 IX1=IX(I)
      IX(I)=IX(J)
      IX(J)=IX1
    3 CONTINUE
    2 CONTINUE
      RETURN
      END
      SUBROUTINE RDINSP(IU,SRNO,ITYP,NW,WL,SOLID,NDY,F0D,TJLC,ERC)
C
C   Read instrumental parameters.(from 'ins.para' file)
C
C --- History
C   2001.02.14 Created by M.Yamano
C   2002.05.11 Renewed by M.Yamano for Ver.5.
C   2002.11.20 Modified.
C   2003.11.05 Added processing for ITYP=30 (new format) data.
C
C --- Input
C IU     I           device No. for reading
C
C --- Output
C SRNO   C*20        instrument S/N
C ITYP   I           data file format type No.
C                     10,11   : PREDE - skyradiometer POM-01L
C                     20,21,22: PREDE - skyradiometer POM-01MKII
C                     30      : PREDE - skyradiometer POM-02
C NW     I           number of wavelengths
C WL     R(NW)       wavelengths[cm]
C SOLID  R(NW)       solid view angles[sr]
C NDY    I           number of F0D data
C F0D    R(NW,NDY)   calibration constants
C TJLC   R(NDY)      F0D time data(Julian day)
C ERC    C*64        ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KDAY  =30)
C
      CHARACTER ERC*(*)
      CHARACTER SRNO*(*)
      DIMENSION WL(KNW),SOLID(KNW),F0D(KNW,KDAY),TJLC(KDAY)
C
      ERC=' '
      READ(IU,*,ERR=800,END=800)
      READ(IU,*,ERR=800,END=800) SRNO
      READ(IU,*,ERR=800,END=800) ITYP
      IF ((ITYP.NE.10).AND.(ITYP.NE.11).AND.
     &    (ITYP.NE.20).AND.(ITYP.NE.21).AND.(ITYP.NE.22).AND.
     &                                      (ITYP.NE.30)) THEN
         ERC='RDINSP: illegal setting of skyradiometer type(ITYP).'
         GOTO 900
      ENDIF
      READ(IU,*,ERR=800,END=800) NW
      IF (NW.GT.KNW) THEN
         ERC='RDINSP: NW.GT.KNW'
         GOTO 900
      ENDIF
      READ(IU,*,ERR=800,END=800) (WL(I),I=1,NW)
      READ(IU,*,ERR=800,END=800) (SOLID(I),I=1,NW)
      READ(IU,*,ERR=800,END=800) NDY
      IF (NDY.GT.KDAY) THEN
         ERC='RDINSP: NDY.GT.KDAY'
         GOTO 900
      ENDIF
      IF (NDY.GT.0) THEN
         DO I=1,NDY
           READ(IU,*,ERR=800,END=800) IYC,IMC,IDC,TMC
     &                               ,(F0D(IW,I),IW=1,NW)
           TJLC(I)=TJUL1(IYC,IMC,IDC,TMC,0.)
         ENDDO
      ENDIF
      GOTO 900
C
  800 ERC='RDINSP: Read error.'
C
  900 RETURN
      END
C
      SUBROUTINE PLYFIT(M1,N,X,Y,A,KM1,WK1,WK2,IER)
C FITTING BY POLYNOMIALS
C--- HISTORY
C 88. 2.25 CREATED
C--- INPUT
C M1     I      ORDER OF POLYNOMIALS + 1.
C N      I      NBR OF DATA.
C X    R(N)     DATA-X.
C Y    R(N)     DATA-Y.
C KM1    I      FIRST DIMENSION SIZE FOR WK2.
C WK1  R(M1)    WORKING AREA.
C WK2 R(KM1,M1) WORKING AREA.
C--- OUTPUT
C A    R(M1)    Y = SUM(I=1,M1) ( A(I)*X**(I-1)).
C IER    I      ERROR INDICATOR.
C$ENDI
C
      DIMENSION X(N),Y(N),A(M1),WK1(M1),WK2(KM1,M1)
C
      DO 1 I=1,M1
      DO 1 J=1,I
      S=0
      DO 2 K=1,N
    2 S=S+X(K)**(I+J-2)
      WK2(I,J)=S
    1 WK2(J,I)=S
      CALL TINVCH(M1,WK2,DT,KM1,IER)
      IF(IER.NE.0) RETURN
C
      DO 3 I=1,M1
      S=0
      DO 4 J=1,N
    4 S=S+Y(J)*X(J)**(I-1)
    3 WK1(I)=S
      DO 5 I=1,M1
      S=0
      DO 6 J=1,M1
    6 S=S+WK2(I,J)*WK1(J)
    5 A(I)=S
C
      RETURN
      END
      SUBROUTINE REGRES(N,DX,DY,A,B,R)
C
C   Get coefficients of regression and correlation.
C
C --- history
C   2000.11.28 Created by M.Yamano
C
C --- Input
C N          I         number of data
C DX         R(N)      X-data array
C DY         R(N)      Y-data array
C
C --- Output
C A,B        R         regression coefficients
C                        Y = A*X + B
C R          R         correlation coefficient of X-Y data
C ----------------------------------------------------------------------
C
      DIMENSION DX(N),DY(N)
C
      X=0.
      Y=0.
      XY=0.
      X2=0.
      Y2=0.
      DO I =1,N
        X=X+DX(I)
        Y=Y+DY(I)
        XY=XY+DX(I)*DY(I)
        X2=X2+DX(I)*DX(I)
        Y2=Y2+DY(I)*DY(I)
      ENDDO
      FN=FLOAT(N)
      CX=X2-X*X/FN
      CY=Y2-Y*Y/FN
      CXY=XY-X*Y/FN
      IF ((CX.GT.0.).AND.(CY.GT.0.)) THEN
         A=CXY/CX
         B=Y/FN-A*X/FN
         R=CXY/SQRT(CX*CY)
      ELSE
         A=-999.
         B=-999.
         R=-999.
      ENDIF
      RETURN
      END
C
      SUBROUTINE CLMEAN(N,DT,DAVE,DSTD,DRMS)
C
C   Get average, STD and RMSD of data array.
C
C --- history
C   2000.12.07 Created by M.Yamano
C   2001.10.09 Modified.
C
C --- Input
C N          I         number of data
C DT         R(N)      data array
C
C --- Output
C DAVE       R         average
C DSTD       R         standard deviation
C DRMS       R         root mean square deviation
C ----------------------------------------------------------------------
C
      DIMENSION DT(N)
C
      S=0.
      S2=0.
      DAVE=0.
      DSTD=0.
      DRMS=0.
      IF (N.GT.0) THEN
         DO I = 1,N
           S=S+DT(I)
           S2=S2+DT(I)*DT(I)
         ENDDO
         FN=FLOAT(N)
         DAVE=S/FN
         IF ((N.GT.1).AND.(S2.GT.S*S/FN)) THEN
            DSTD=SQRT((S2-S*S/FN)/(FN-1.))
            DRMS=SQRT((S2-S*S/FN)/FN)
         ENDIF
      ENDIF
      RETURN
      END
C
      FUNCTION TJUL1(IY,IM,ID,TM,ALNGS)
C Julian day
C--- HISTORY
C 98.  6. 15  Created
C--- INPUT
C IY      I      YEAR YYYY
C IM      I      MONTH
C ID      I      DAY
C TM      R      TIME (HOUR, I.E., 10:30 = 10.5)
C ALNGS   R      LONGITUDE FOR LOCAL STANDARD TIME (DEGREE)
C--- OUTPUT
C TJUL    R      Julian day from 0:00 of Jun. 1, 1900
      DT=51.0/3600.0
      U=TM-ALNGS/15
      E=U+DT
      XI=INT((14-IM)/12)
      XJ=-XI+IY+4800
      X=INT(1461*XJ/4.0)+INT((12*XI+IM-2)*367/12.0)
     & -INT(INT(XJ/100.0+1)*3/4.0)+ID-2447095.5+E/24
C     TJUL1=X/36525
      TJUL1=X-0.5
      return
      end
      SUBROUTINE      TINVCH(N,A,DT,NN,INDER)
C INVERSION OF SYMMETRIC POSITIVE DEFINITE MATRIX
C     SQUARE ROOT METHOD
C     MATRIX A IS SYMMETRIC AND POSITIVE DEFINITE
C--- HISTORY
C 90. 1.17  REGISTERED
C--- INPUT
C N        I        DIMENSION OF THE MATRIX
C A      R(NN,N)    MATRIX
C NN       I        SIZE OF THE FIRST ARGUMENT OF A
C--- OUTPUT
C A                 INVERSE OF THE MATRIX
C DT       R        DETERMINATION OF THE MATRIX
C INDER    I        0: NORMAL TERMINATION
C$ENDI
      DIMENSION    A(NN,N)
      INDER=0
      IF(N-1)   1,2,3
    1 INDER=-1
      WRITE(*,602)    N
  602 FORMAT(1H0,8H *** N =  ,I5,5X,30HN SHOULD BE POSITIVE IN TINVCH  )
      RETURN
    2 DT=A(1,1)
      A(1,1)=1.0/A(1,1)
      IF(DT.LE.0.0)   GO TO  60
      RETURN
    4 INDER=-1
      WRITE(*,603)   N,NN
  603 FORMAT(1H0,8H *** N =  ,I5,5X,4HNN =,I5/
     & 60H    N SHOULD BE LESS THAN OR EQUAL TO NN IN TINVCH          )
      RETURN
    3 IF(N.GT.NN)    GO TO  4
      DT=A(1,1)
      IF(DT.LE.0.0)   GO TO  60
      A(1,1)=SQRT(DT)
      DO  5  J=2,N
    5 A(1,J)=A(1,J)/A(1,1)
      DO  10  K=2,N
      Z=A(K,K)
      K1=K-1
      DO  11   J=1,K1
   11 Z=Z-A(J,K)*A(J,K)
      IF(Z.LE.0.0)   GO TO  60
      DT=DT*Z
      A(K,K)=SQRT(Z)
      IF(K.EQ.N) GO TO 10
      Y=1.0/A(K,K)
      J1=K+1
      DO 13   J=J1,N
      Z=A(K,J)
      DO  14  I=1,K1
   14 Z=Z-A(I,K)*A(I,J)
   13 A(K,J)=Z*Y
   10 CONTINUE
      A(N,N)=1.0/A(N,N)
      DO   40   IA=2,N
      I=N-IA+1
      Y=1.0/A(I,I)
      A(I,I)=Y
      I1=I+1
      DO  41    JA=I1,N
      J=N-JA+I1
      Z=0.0
      DO  42    K=I1,J
   42 Z=Z-A(I,K)*A(K,J)
   41 A(I,J)=Z*Y
   40 CONTINUE
      DO  50  J=1,N
      DO  51   I=J,N
      Z=0.0
      DO  52   K=I,N
   52 Z=Z+A(I,K)*A(J,K)
   51 A(I,J)=Z
      DO  53  I=1,J
   53 A(I,J)=A(J,I)
   50 CONTINUE
      RETURN
   60 INDER=1
      WRITE(*,601)
  601 FORMAT(1H0,52H*** GIVEN MATRIX TO TINVCH IS NOT POSITIVE DEFINITE.
     &,40H   RETURN WITH NO FURTHER CALCULATION     )
      RETURN
      END
