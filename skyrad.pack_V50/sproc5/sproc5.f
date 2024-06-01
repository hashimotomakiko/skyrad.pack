C$OUTF sproc5.f
C$LIBF /Users/yamanomaki/Fortran/LIB/SKYLIB/A
C$LIBF /Users/yamanomaki/Fortran/LIB/POMLIB/A
C$LIBF /Users/yamanomaki/Fortran/LIB/MATLIB/A
C$LIBF /Users/yamanomaki/Fortran/LIB/BASLIB/A
C$LIBF /Users/yamanomaki/Fortran/LIB_Nakajima/LBC1/A
C$LIBF /Users/yamanomaki/Fortran/LIB_Nakajima/LBM1/A
C$LIBF /Users/yamanomaki/Fortran/LIB_Nakajima/LBM2/A
C$LIBF /Users/yamanomaki/Fortran/LIB_Nakajima/LBR1/A
C$LIBF /Users/yamanomaki/Fortran/LIB_Nakajima/LBR2/A
C
C --- Note
C   Skyradiometer data analyzing program
C   For PREDE instruments POM-01L/POM-01MKII/POM-01/POM-02
C
C  (2) Sky radiance analyses - Ver.5.0
C
C --- History
C   2007.10.31 Made from Ver.4.2 by M.Yamano (sub REDML5)
C   2008.01.18 Modified.
C   2008.10.17 Modified.
C   2009.01.20 Start on Mac(g95).
C   2011.08.31 Modified by M.Hashimoto.
C   2024.01.17 Modified.
C   2024.05.29 Modified. (KCR:15 => 20, KCI:10 => 20) : Info.from Momoi-san, 2018
C   2024.05.29 Modified. (KNMEAS: 1240 => KNANG*KNW+2*KNW) * Info.from Momoi-san, 2018
C   2024.05.29 Debug (sub. RDDAT5) : from Info.from Momoi-san, 2018
C
C --- I/O files
C  sproc.par            I    processing control file
C  fname                I    processing dates file
C  ins.para             I    instrumental parameters file(any name is OK)
C  METEO.DAT            I    meteorological data file(any name is OK)
C  MIEKER               I    kernel datafile
C  DT5/(yy)yymmdd.DT5   I    measured datafile
C  Tag/(yy)yymmdd.tag   I    tagfile for input data
C  Out/(yy)yymmdd.out   O    full results file(if IOUT>0)
C  Par/(yy)yymmdd.par   O    optical parameters results file(if IPAR>0)
C  Vol/(yy)yymmdd.vol   O    volume spectrum results file(if IVOL>0)
C  Aur/(yy)yymmdd.aur   O    sky radiance results file(if IAUR>0)
C  Phs/(yy)yymmdd.phs   O    phase function results file(if IPHS>0)
C  F0d/(yy)yymmdd.w??   O    input data file for F0 determination(if IPF0>0)
C  sproc.tag            O    processing tagfile
C  sproc5.log           O    processing logfile
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
      PARAMETER (KNANG =60)
      PARAMETER (KNSIZE=60)
      PARAMETER (KINTVL=60)
      PARAMETER (KKDT=KDT*KNW)
      PARAMETER (KCR   =20)
      PARAMETER (KCI   =20)
      PARAMETER (KNA1  =50)
      PARAMETER (KNFI  =1000)
      PARAMETER (KNLN  =4)
      PARAMETER (KNMEAS=KNANG*KNW+2*KNW,KNA1U=KNA1)
      PARAMETER (KDAY  =30)
C
      PARAMETER (PI=3.141592653589793,RAD=PI/180.0)
C
      CHARACTER ERC*64
      CHARACTER FLI*80,FLD*80,FLO*80,FLT*80,LB*2,ATR(6)*4,DIR(6)*4
      CHARACTER ATAG*60,HEAD*80,ASTR(4)*4
C
      DIMENSION CONCA(KNLN),CONCM(KNLN)
      DIMENSION TJLM(KDAY),PRSD(KDAY),O3D(KDAY)
     &         ,TJLC(KDAY),F0D(KNW,KDAY)
      DIMENSION WL(KNW),CR(KNW),CI(KNW),GA(KNW),SOLID(KNW)
     &         ,IANL(KNW)
      DIMENSION WLD(KNW),AUR(KNW,KDT),OPT(KNW)
      DIMENSION TH(KNW,KDT),FI(KNW,KDT),SCA(KNW,KDT),IAUR(KNW,KDT)
      DIMENSION TAUR(KNW),TAUO3(KNW)
      DIMENSION SIZE(KNSIZE),VOL(KNSIZE),CMX(KNSIZE)
     &         ,RMODE(KNSIZE),CMODE(KNSIZE),SSL(KNSIZE)
      DIMENSION TA(KNW),WA(KNW),AURC(KNW,KDT),OPTC(KNW)
      DIMENSION THETA(KNANG),PHRTR(KNMEAS),WCR(KNW),WCI(KNW)
      DIMENSION SRETV(KNSIZE),SRETR(KNW),SRETI(KNW)
      DIMENSION IOUT(6)
C
      DATA DIR /'Out/','Par/','Vol/','Aur/','Phs/','F0d/'/
      DATA ATR /'.out','.par','.vol','.aur','.phs','.w00'/
      DATA ASTR /'    ',' *  ', ' ** ','    '/
C
      ERC=' '
      IUI=1
      OPEN (IUI,FILE='sproc.par',STATUS='OLD')
      CALL RDPARS(IUI,IOUT,NW,NO3,NWV,WL,IANL,CR,CI,GA,SOLID,ITYP
     &           ,NLN,CONCA,CONCM,MDY,TJLM,PRSD,O3D,NDY,TJLC,F0D
     &           ,JTAU,IDCR,IDCI,INVM,IPLC,IPHS
     &           ,ANGMN,ANGMX,LOOPMX,EPSA,EPSA1,EPSA2
     &           ,NSIZE,RMIN,RMAX,SIZE,NMODE,RMODE,SSL,ERC)
      CLOSE(IUI)
      IF (ERC.NE.' ') GOTO 890
      IF ((JTAU.LT.0).AND.(IOUT(4).GT.0)) IOUT(4)=0
C
      IUK=2
      OPEN (IUK,FILE='MIEKER',STATUS='OLD')
      IUF=3
      OPEN (IUF,FILE='fname',STATUS='OLD')
      IUT=4
      IUTG=7
      OPEN (IUTG,FILE='sproc.tag',STATUS='UNKNOWN')
      IUL=9
      OPEN (IUL,FILE='sproc5.log',STATUS='UNKNOWN')
      IUO=10
      IUW=20
      NTG=0
C
   10 READ(IUF,15,END=900,ERR=900) FLI
   15 FORMAT(A80)
      CALL GETPOS(FLI,'.',80,IL)
      IF (IL.EQ.1) GOTO 900
      IF (IL.GT.0) FLI=FLI(1:IL-1)
      CALL UTLSPC(FLI,NFNM)
      CLOSE(IUT)
      FLT='Tag/'//FLI(1:NFNM)//'.tag'
      OPEN (IUT,FILE=FLT,STATUS='OLD')
      READ(IUT,25,END=10,ERR=10) ATAG
      DO I=1,5
        IF (IOUT(I).GT.0) THEN
           CLOSE(IUO+I)
           FLO=DIR(I)//FLI(1:NFNM)//ATR(I)
           OPEN (IUO+I,FILE=FLO,STATUS='UNKNOWN')
        ENDIF
      ENDDO
      IF (IOUT(6).GT.0) THEN
         DO IW=1,NW
           WRITE(LB,17) IW
   17      FORMAT(I2)
           IF (IW.LT.10) THEN
              IST=2
           ELSE
              IST=1
           ENDIF
           CLOSE(IUW+IW)
           FLO=DIR(6)//FLI(1:NFNM)//ATR(6)(1:1+IST)//LB(IST:2)
           OPEN (IUW+IW,FILE=FLO,STATUS='UNKNOWN')
         ENDDO
         NN=0
      ENDIF
      CLOSE(IUI)
      FLD='DT5/'//FLI(1:NFNM)//'.DT5'
      OPEN (IUI,FILE=FLD,STATUS='OLD')
C
   20 CALL RDDAT5(IUI,IY,IM,ID,JH,JM,JS,TIM,TJL,ALNG,ALAT,ALT
     &               ,NA,TH0,FI0,DST,NWD,WLD,TH,FI,AUR,IFLG,ERC)
      IF (IFLG.NE.0) THEN
         IF (IFLG.LT.0) THEN
            WRITE(IUL,105) IY,IM,ID,TIM,ERC
            write(*,105) IY,IM,ID,TIM,ERC
         ENDIF
         GOTO 10
      ENDIF
      READ(IUT,25) ATAG
   25 FORMAT(A60)
C
      CALL CKDTS5(NW,WL,CR,CI,GA,SOLID,NDY,F0D
     &           ,NWD,WLD,NA,TH,FI,TH0,AUR,ITYP,ERC)
      IF (ERC.NE.' ') GOTO 100
C
      CALL STDTS5(TJL,TH0,NA,NW,WL,TH,FI,SCA,AUR,OPT,ANGMN,ANGMX
     &           ,TAUR,TAUO3,MDY,TJLM,PRSD,O3D,NDY,TJLC,F0D,ERC)
      IF (ERC.NE.' ') GOTO 100
C
C Get refractive index & volume spectrum
C
      CALL REDML5(JTAU,IPLC,INVM,IDCR,IDCI,ANGMN,ANGMX,LOOPMX
     &     ,IUK,EPSA,EPSA1,EPSA2,GA,NLN,CONCA,CONCM,CR,CI,NW,NA
     &     ,NMODE,RMIN,RMAX,WL,TAUR,TAUO3,TH0,TH,FI,SCA,OPT,AUR
     &     ,SOLID,LOOP,EROR,SRETV,SRETR,SRETI,MRK
     &     ,TA,WA,OPTC,AURC,IAUR,SCAMX,JDCR,JDCI
     &     ,NANG,THETA,PHRTR,RMODE,SSL,CMODE,WCR,WCI,IFLAG,ERC)
      IF (IFLAG.EQ.-9) GOTO 890
C
C Write results
C
  100 WRITE(IUL,105) IY,IM,ID,TIM,ERC
      write(*,105) IY,IM,ID,TIM,ERC
  105 FORMAT(I4,2I3,F6.2,1X,A63)
      IF (ERC.EQ.' ') THEN
         NTG=NTG+1
         NN=NN+1
         IF (NTG.EQ.1) WRITE(IUTG,110)
  110    FORMAT('  TNo yyyy mm dd Hour   Long    Lat  '
     &         ,'  Hs    SA(max) SA  r i  error  (LP)')
         WRITE(HEAD,115) NTG,ATAG(5:49)
     &                  ,SCAMX,JDCR,JDCI,EROR,LOOP,ASTR(MRK+1)
  115    FORMAT(I5,A45,F6.1,2I2,F8.4,' (',I2,')',A4)
         WRITE(IUTG,120) HEAD
  120    FORMAT(A80)
         IF (INVM.GE.20) THEN
            CALL VOLLOG(NMODE,RMODE,SSL,CMODE,NSIZE,SIZE,VOL)
            CALL INTERV(NMODE,RMODE,SRETV,NSIZE,SIZE,CMX)
         ELSE
            DO I=1,NSIZE
              VOL(I)=CMODE(I)
              CMX(I)=SRETV(I)
            ENDDO
         ENDIF
         IF (IOUT(1).GT.0)
     &      CALL WTOUT5(IUO+1,HEAD,INVM,NSIZE,SIZE,VOL
     &            ,NMODE,RMODE,CMODE,WCR,WCI,NW,WL,OPT,TA,WA
     &            ,NA,TH,FI,SCA,IAUR,AUR,AURC,NANG,THETA,PHRTR,JTAU
     &            ,CMX,SRETR,SRETI)
         IF (IOUT(2).GT.0)
     &      CALL WTPAR5(IUO+2,HEAD,NW,WL,OPT,TA,WA,WCR,WCI,SRETR,SRETI)
         IF (IOUT(3).GT.0) CALL WTVOL5(IUO+3,HEAD,NSIZE,SIZE,VOL,CMX)
         IF (IOUT(4).GT.0)
     &      CALL WTAUR4(IUO+4,HEAD,NW,NA,TH,FI,SCA,IAUR,AUR,AURC)
         IF (IOUT(5).GT.0)
     &      CALL WTPHS0(IUO+5,HEAD,NW,NANG,THETA,PHRTR)
         IF (IOUT(6).GT.0) THEN
            EM=1/COS(TH0*RAD)
            CALL WTDATF(IUW,NN,IY,IM,ID,TIM,EM,EROR,NO3,NWV
     &                 ,NW,WL,TAUR,TAUO3,AUR,TA,WA)
         ENDIF
      ENDIF
      GOTO 20
C
  890 WRITE(*,*) 'ERROR CODE: ',ERC
  900 STOP
      END
C
      SUBROUTINE RDPARS(IU,IOUT,NW,NO3,NWV,WL,IANL,CR,CI,GA,SOLID,ITYP
     &                 ,NLN,CONCA,CONCM,MDY,TJLM,PRSD,O3D,NDY,TJLC,F0D
     &                 ,JTAU,IDCR,IDCI,INVM,IPLC,IPHS
     &                 ,ANGMN,ANGMX,LOOPMX,EPSA,EPSA1,EPSA2
     &                 ,NSIZE,RMIN,RMAX,SIZE,NMODE,RMODE,SSL,ERC)
C
C   Read parameters from file('sproc.par')
C
C --- history
C   2001.01.25 Created by M.Yamano
C   2002.02.28 Bug corrected.
C   2002.05.15 Renewed by M.Yamano for Ver.5.
C   2002.11.15 Added IOUT(6)
C   2002.11.20 Modified.
C   2004.01.21 Added ANGMN,ANGMX check
C   2006.05.25 Added output argument ITYP.
C
C --- Input
C IU          I         device No. for reading
C
C --- Output
C IOUT        I(6)      (full/par/vol/aur/phs/w00)
C                        output option (1: create/0: not)
C NW          I         number of wavelengths
C NO3         I         no. of wavelength for ozone absorption
C NWV         I         no. of wavelength for water vapor absorption
C WL          R(NW)     wavelengths[cm]
C IANL        I(NW)     processing option(positive: processing/ 0: not)
C CR          R(NW)     real part of refractive index
C CI          R(NW)     imaginary part of refractive index
C GA          R(NW)     assumed ground albedo
C SOLID       R(NW)     solid view angles[sr]
C NLN         I         number of layers
C CONCA       R(NLN)    fraction of aerosol   optical thickness
C CONCM       R(NLN)    fraction of molecular optical thickness
C MDY         I         number of meteorological data(pressure and ozone)
C TJLM        R(MDY)    time data of pressure and ozone
C PRSD        R(MDY)    pressure data[atm]
C O3          R(MDY)    ozone data[cm,STP]
C NDY         I         number of calibration constant data
C TJLC        R(NDY)    time data of calibration constants(F0)
C F0D         R(NW,NDY) calibration constants data
C JTAU        I         controlling input data
C                        2: ver.5 analysis(in case of no direct data)
C                        1: OPT is observed and AUR+OPT for inversion.
C                        0: only AUR for inversion (no use of OPT)
C                       -1: only OPT for inversion
C IDCR        I         0: give CR; 1: determine CR by this routine.
C IDCI        I         0: give CI; 1: determine CI by this routine.
C INVM        I          =<19: box kernel
C                        >=20: log-normal multi mode kernel
C IPLC        I         polarization correction control
C                        1: with polarization correction
C                        0: without correction
C IPHS        I         phase function control used only for ver3
C                        0: use Mie-inverted phase func.for radiative transfer
C                        1: use extrapolated observed phase fuction
C ANGMN       R         minimum scattering angle[deg] for useful sky radiance
C ANGMX       R         maximum scattering angle[deg] for useful sky radiance
C LOOPMX      I          >0  Maximum number of iteration
C EPSA        R         ABSOLUTE CONVERGENCE CRITERION.
C EPSA1       R         RELATIVE CONVERGENCE CRITERION.
C EPSA2       R         Give-up convergence criterion after 2nd-loop.
C NSIZE       I         number of particle radii
C RMIN,RMAX   R         Minimum and maximum radii[cm] for size integration.
C SIZE        R(NSIZE)  particle radii[cm]
C NMODE       I         number of mode
C RMODE       R(NSIZE)  mode radius[cm]
C SSL         R(NSIZE)  log-dispersion
C                dV/dlnr = sum(n=1,N) Cn * exp(-(ln(r/SIZE(n))/SSL(n))**2/2)
C ERC         C*64      ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KNLN  =4)
      PARAMETER (KNSIZE=60)
      PARAMETER (KDAY  =30)
C
      CHARACTER ERC*(*)
      CHARACTER PFIL*80
      CHARACTER CMNT*40,SITE*20,SRNO*20
      DIMENSION WL(KNW),IANL(KNW),CR(KNW),CI(KNW),GA(KNW),SOLID(KNW)
      DIMENSION CONCA(KNLN),CONCM(KNLN)
      DIMENSION TJLM(KDAY),PRSD(KDAY),O3D(KDAY)
     &         ,TJLC(KDAY),F0D(KNW,KDAY)
      DIMENSION SIZE(KNSIZE),RMODE(KNSIZE),SSL(KNSIZE)
      DIMENSION WLI(KNW),SOLIDI(KNW),F0DI(KNW,KDAY),IWANL(KNW)
      DIMENSION IOUT(6)
C
      ERC=' '
      READ(IU,*,ERR=800,END=800) (IOUT(I),I=1,6)
C
C --- read 'ins.para' parameters
C
      READ(IU,*,ERR=800,END=800) PFIL
      CALL UTLSPC(PFIL,NL)
      IF (NL.LE.0) THEN
         ERC='RDPARS: illegal ins.para filename'
         GOTO 900
      ENDIF
      OPEN (IU+1,FILE=PFIL(1:NL),STATUS='OLD')
      CALL RDINSP(IU+1,SRNO,ITYP,NWI,WLI,SOLIDI,NDY,F0DI,TJLC,ERC)
      CLOSE (IU+1)
      IF (ERC.NE.' ') GOTO 900
C
C --- read 'METEO.DAT' filename
C
      READ(IU,*,ERR=800,END=800) PFIL
      CALL UTLSPC(PFIL,NL)
      IF (NL.LE.0) THEN
         ERC='RDPARS: illegal METEO.DAT filename'
         GOTO 900
      ENDIF
      OPEN (IU+1,FILE=PFIL(1:NL),STATUS='OLD')
      CALL RDMETP(IU+1,CMNT,SITE,MDY,PRSD,O3D,TJLM,ERC)
      CLOSE(IU+1)
      IF (ERC.NE.' ') GOTO 900
C
      READ(IU,*)
      READ(IU,*,ERR=800,END=800) NLN
      IF (NLN.GT.KNLN) THEN
         ERC='RDPARS: NLN.GT.KNLN !'
         GOTO 900
      ENDIF
      READ(IU,*,ERR=800,END=800) (CONCA(L),L=1,NLN)
      READ(IU,*,ERR=800,END=800) (CONCM(L),L=1,NLN)
C
      READ(IU,*)
      READ(IU,*,ERR=800,END=800) NW,NO3,NWV
      IF (NW.GT.KNW) THEN
         ERC='RDPARS: NW.GT.KNW !'
         GOTO 900
      ENDIF
      READ(IU,*,ERR=800,END=800) (WL(I),I=1,NW)
      READ(IU,*,ERR=800,END=800) (IWANL(I),I=1,NW)
      READ(IU,*,ERR=800,END=800) (CR(I),I=1,NW)
      READ(IU,*,ERR=800,END=800) (CI(I),I=1,NW)
      READ(IU,*,ERR=800,END=800) (GA(I),I=1,NW)
C
      IF (NW.GT.NWI) THEN
         ERC='RDPARS: NW mismatching !'
         GOTO 900
      ELSE
         N=0
         DO WHILE(N.LT.NW)
           N=N+1
           K=0
           DO WHILE((K.LT.NWI).AND.(K.GE.0))
             K=K+1
             IF (ABS(WLI(K)-WL(N)).LT.0.001E-4) K=-K
           ENDDO
           IF (K.LT.0) THEN
              K=ABS(K)
              WL(N)=WLI(K)
              SOLID(N)=SOLIDI(K)
              IF (NDY.GT.0) THEN
                 DO J=1,NDY
                   F0D(N,J)=F0DI(K,J)
                 ENDDO
              ENDIF
           ELSE
              SOLID(N)=-999.
              IWANL(N)=0
              IF (NDY.GT.0) THEN
                 DO J=1,NDY
                   F0D(N,J)=-999.
                 ENDDO
              ENDIF
           ENDIF
         ENDDO
      ENDIF
C
      N=0
      DO IW=1,NW
        IF (IWANL(IW).GT.0) THEN
           N=N+1
           IANL(N)=IWANL(IW)
           WL(N)=WL(IW)
           CR(N)=CR(IW)
           CI(N)=CI(IW)
           GA(N)=GA(IW)
           SOLID(N)=SOLID(IW)
           IF (NDY.GT.0) THEN
              DO J=1,NDY
                F0D(N,J)=F0D(IW,J)
              ENDDO
           ENDIF
           IF (IW.EQ.NO3) NO3=N
           IF (IW.EQ.NWV) NWV=N
        ELSE
           IF (IW.EQ.NO3) NO3=0
           IF (IW.EQ.NWV) NWV=0
        ENDIF
      ENDDO
      NW=N
      IF (NW.LE.0) THEN
         ERC='RDPARS: No available wavelength !'
         GOTO 900
      ENDIF
C
      READ(IU,*)
      READ(IU,*,ERR=800,END=800) JTAU,IDCR,IDCI,INVM,IPLC,IPHS
      READ(IU,*,ERR=800,END=800) ANGMN,ANGMX
      IF ((JTAU.GE.0).AND.(ANGMN.GE.ANGMX)) THEN
         ERC='RDPARS: ANGMN.GE.ANGMX !'
         GOTO 900
      ENDIF
      READ(IU,*,ERR=800,END=800) LOOPMX,EPSA,EPSA1,EPSA2
      READ(IU,*,ERR=800,END=800) NSIZE,RMIN,RMAX
      IF (NSIZE.GT.KNSIZE) THEN
         ERC='RDPARS: NSIZE.GT.KNSIZE !'
         GOTO 900
      ENDIF
      IF (RMAX.GT.RMIN) THEN
         DL=LOG(RMAX/RMIN)/FLOAT(NSIZE)
         SZL=LOG(RMIN)-DL
         DO I=1,NSIZE
           SZL=SZL+DL
           SIZE(I)=EXP((SZL*2+DL)/2)
         ENDDO
      ELSE
         ERC='RDPARS: RMIN.GE.RMAX !'
         GOTO 900
      ENDIF
C
      IF (INVM.GE.20) THEN
         READ(IU,*,ERR=800,END=800) NMODE
         READ(IU,*,ERR=800,END=800) (RMODE(I),I=1,NMODE)
         READ(IU,*,ERR=800,END=800) (SSL(I),I=1,NMODE)
      ELSE
         NMODE=NSIZE
         DO I=1,NMODE
           RMODE(I)=SIZE(I)
         ENDDO
      ENDIF
      GOTO 900
C
  800 ERC='RDPARS: Read Error !'
C
  900 RETURN
      END
C
      SUBROUTINE RDDAT5(IU,IY,IM,ID,JH,JM,JS,TIM,TJL,ALNG,ALAT,ALT
     &                 ,NA,TH0,FI0,DST,NW,WL,TH,FI,AUR,IFLG,ERC)
C
C   Read datafile(yymmdd.DT5) for analysis Ver.5.0.
C
C --- History
C   2001.01.26 Created by M.Yamano
C   2007.05.12 Debugged.
C   2024.05.29 Debugged. (ref Momoi-san debug)
C
C --- Input
C IU          I         device No. for reading
C
C --- Output
C IY,IM,ID    I         measured date IY/IM/ID
C JH,JM,JS    I         measured time JH:JM:JS
C TIM         R         measured time in hour
C TJL         R         Julian day from 0:00 of Jun. 1, 1900
C ALNG,ALAT   R         longitude,latitude in degree of measurement site
C ALT         R         altitude in meter of measurement site
C NA          I         number of scattering angles
C TH0,FI0     R         solar zenith and azimuth angles in degree
C DST         I         distance between the sun and the earth
C NW          I         number of wavelength
C WL          R(NW)     wavelength
C TH,FI       R(NW,NA)  zenith and azimuth angles of measured direction
C AUR         R(NW,NA)  intensity F(sca=0) or I/EM/F/(solid view angle)
C IFLG        I         return flag(0:normal /9:end /-9:error)
C ERC         C*64      ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
C
      CHARACTER ERC*(*)
      CHARACTER CH*1
      DIMENSION WL(KNW),TH(KNW,KDT),FI(KNW,KDT),AUR(KNW,KDT)
C
      ERC=' '
      IFLG=0
      READ(IU,*,ERR=800,END=890) IY,IM,ID,JH,JM,JS,TM
      READ(IU,*,ERR=800,END=890) ALNGS,ALNG,ALAT,ALT
      TIM=JH+JM/60.0+JS/3600.0
      TJL=TJUL1(IY,IM,ID,TIM,ALNGS)
C
      READ(IU,*,ERR=800,END=890) NA,TH0,FI0,DST
      READ(IU,*,ERR=800,END=890) NW
      IF (NW.GT.KNW) THEN
         ERC='RDDAT5: NW.GT.KNW !'
         GOTO 900
      ENDIF
      READ(IU,*,ERR=800,END=890) (WL(IW),IW=1,NW)
      READ(IU,10) CH
   10 FORMAT(A1)
      DO I=1,NA
        READ(IU,*,ERR=800,END=890) (TH(IW,I),FI(IW,I),AUR(IW,I),IW=1,NW)

C 20240529 modified : from Momoi-san
C Modified by Momoi 2018.07.02: Coordinate in principal plane
C -90<th<0, 0<fi<90 => 0<th<90, 90<fi<180
        IF( TH(1,I).LT.0.0 ) THEN
          TH(1:NW,I) = - TH(1:NW,I)
          FI(1:NW,I) = 180.0 + FI(1:NW,I)
        ENDIF
      
      ENDDO
      GOTO 900
C
  800 ERC='RDDAT5: Read error.'
      IFLG=-9
      GOTO 900
  890 ERC='RDDAT5: End of file.'
      IFLG=9
C
  900 RETURN
      END
C
      SUBROUTINE CKDTS5(NW,WL,CR,CI,GA,SOLID,NDY,F0D
     &                 ,NWD,WLD,NA,TH,FI,TH0,AUR,ITYP,ERC)
C
C   Check input data for analysis Ver.5.0.
C
C --- history
C   2004.01.21 Created from old version for Ver.4.2.
C   2004.01.30 Debugged.
C   2006.03.13 Added argument ITYP.
C   2006.12.16 Added check for the last (all positive) diffuse data.
C   2007.04.30 Debugged.
C   2024.01.17 Edited for Ver.5.0.
C
C --- Input
C
C --- Output
C ERC           C*64      ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
      PARAMETER (KDAY  =30)
C
      CHARACTER ERC*(*)
      DIMENSION F0D(KNW,KDAY)
      DIMENSION WL(KNW),CR(KNW),CI(KNW),GA(KNW),SOLID(KNW)
      DIMENSION WLD(KNW),TH(KNW,KDT),FI(KNW,KDT),AUR(KNW,KDT)
      DIMENSION THC(KNW,KDT),FIC(KNW,KDT),AURC(KNW,KDT)
C
      DATA AURMX,RLMT /10.0, 0.2/
      DATA DANG/0.2/
C
      ERC=' '
      IW=0
      N=0
      DO WHILE(N.LT.NW)
        N=N+1
        K=0
        DO WHILE((K.LT.NWD).AND.(K.GE.0))
          K=K+1
          IF (ABS(WLD(K)-WL(N)).LT.0.001E-4) K=-K
        ENDDO
        IF (K.LT.0) THEN
           IW=IW+1
           WL(IW)=WL(N)
           CR(IW)=CR(N)
           CI(IW)=CI(N)
           GA(IW)=GA(N)
           SOLID(IW)=SOLID(N)
           IF (NDY.GT.0) THEN
              DO I=1,NDY
                F0D(IW,I)=F0D(N,I)
              ENDDO
           ENDIF
           K=ABS(K)
           DO J=1,NA
             THC(IW,J)=TH(K,J)
             FIC(IW,J)=FI(K,J)
             AURC(IW,J)=AUR(K,J)
           ENDDO
        ENDIF
      ENDDO
      NW=IW
      IF (NW.LE.0) THEN
         ERC='CKDTS5: available NW.LE.0 !'
         GOTO 900
      ENDIF
C
      J=NA+1
      DO WHILE(J.GT.1)
        J=J-1
        N=0
        DO IW=1,NW
          IF (AURC(IW,J).GT.0.0) N=N+1
        ENDDO
        IF (N.EQ.NW) J=-J
      ENDDO
      NA=ABS(J)
C
      NE=0
      NO=0
      DO IW=1,NW
        DO J=1,NA
          TH(IW,J)=THC(IW,J)
          FI(IW,J)=FIC(IW,J)
          AUR(IW,J)=AURC(IW,J)
          IF (AUR(IW,J).GT.0.) THEN
             NE=NE+1
             IF (J.GT.1) THEN
                IF (AUR(IW,J).GT.AURMX) NO=NO+1
             ENDIF
          ENDIF
        ENDDO
C
        IF ((ITYP.GE.20).AND.(ITYP.LT.30)) THEN
           IZERO=1
        ELSE
           IZERO=0
        ENDIF
        N=0
        DO WHILE((IZERO.EQ.0).AND.(N.LT.NA))
          N=N+1
          IF (AUR(IW,N).GT.0.) THEN
             IF ((ABS(TH(IW,N)-TH0).LE.DANG)
     &      .AND.(ABS(FI(IW,N)).LE.DANG)) IZERO=N
          ENDIF
        ENDDO
        IF (IZERO.NE.1) THEN
           ERC='CKDTS5: First datum in not direct'
           GOTO 900
        ENDIF
      ENDDO
      IF (NE.GT.0) THEN
         RAT=FLOAT(NO)/FLOAT(NE)
         IF (RAT.GT.RLMT) ERC='CKDTS5: abnormal obs. data'
      ELSE
         ERC='CKDTS5: No available data.'
      ENDIF
C
  900 RETURN
      END
C
      SUBROUTINE STDTS5(TJL,TH0,NA,NW,WL,TH,FI,SCA,AUR,OPT,ANGMN,ANGMX
     &            ,TAUR,TAUO3,MDAY,TJLM,PRSD,O3D,NDAY,TJLC,F0D,ERC)
C
C   Set input data for analysis Ver.5.0.
C
C --- history
C   2004.01.21 Created from old version for Ver.4.2.
C   2024.01.17 Edited for Ver.5.0.
C
C --- Input
C TJL      R           Julian day from 0:00 of Jun. 1, 1900
C TH0      R           solar zenith angle in degree
C NA       I           number of scattering angles
C NW       I           number of wavelength
C WL       R(NW)       wavelength
C TH,FI    R(NW,NA)    zenith and azimuth angles of measured direction
C AUR      R(NW,NA)    intensity F(sca=0) or I/EM/F/(solid view angle)
C ANGMN    R           Minimum scattering angle for useful sky radiance
C ANGMX    R           Maximum scattering angle for useful sky radiance
C MDAY     I           number of meteorological data
C TJLM     R(MDAY)     time data(Julian day) for meteorological data
C PRSD     R(MDAY)     pressure data[atm]
C O3D      R(MDAY)     ozone data[cm,STP]
C NDAY     I           number of F0D data
C TJLC     R(NDAY)     F0D time data(Julian day)
C F0D      R(NW,NDAY)  calibration constants
C
C --- Output
C NA       I           number of available scattering angles
C TH,FI    R(NW,NA)    zenith and azimuth angles of measured direction
C SCA      R(NW,NA)    Scattering angles of observed data
C AUR      R(NW,NA)    intensity F(sca=0) or I/EM/F/(solid view angle)
C OPT      R(NW)       Observed aerosol optical thickness
C TAUR     R(NW)       Rayleigh optical thickness
C TAUO3    R(NW)       Ozone optical thickness
C ERC      C*64        ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
      PARAMETER (KDAY  =30)
C
      PARAMETER (PI=3.141592653589793,RAD=PI/180.0)
C
      CHARACTER ERC*(*)
      DIMENSION TJLM(KDAY),PRSD(KDAY),O3D(KDAY)
     &         ,TJLC(KDAY),F0D(KNW,KDAY),F00(KDAY)
      DIMENSION WL(KNW),AUR(KNW,KDT),OPT(KNW)
      DIMENSION TH(KNW,KDT),FI(KNW,KDT),SCA(KNW,KDT)
      DIMENSION TAUR(KNW),TAUO3(KNW)
C
      DATA DANG /0.5/
      DATA TAUMX /10./
C
      AM00=COS(TH0*RAD)
      EM=1./COS(TH0*RAD)
      DO IW=1,NW
        DO I=1,NA
          AM=COS(TH(IW,I)*RAD)
          SCA(IW,I)=ACOS(AM00*AM+SQRT((1.-AM00**2)*(1.-AM**2))
     &             *COS(FI(IW,I)*RAD))/RAD
        ENDDO
      ENDDO
C
      IF (NA.GT.1) THEN
         L=1
         DO I=2,NA
           IFG=0
           N=0
           DO WHILE ((N.LT.NW).AND.(IFG.EQ.0))
             N=N+1
             IF ((SCA(N,I)+DANG.GE.ANGMN)
     &                        .AND.(SCA(N,I)-DANG.LE.ANGMX)) IFG=1
           ENDDO
           IF (IFG.GT.0) THEN
              L=L+1
              DO IW=1,NW
                TH(IW,L)=TH(IW,I)
                FI(IW,L)=FI(IW,I)
                SCA(IW,L)=SCA(IW,I)
                AUR(IW,L)=AUR(IW,I)
              ENDDO
           ENDIF
         ENDDO
         NA=L
      ENDIF
C
      O3=AINTP(TJL,MDAY,TJLM,O3D)
      PRS=AINTP(TJL,MDAY,TJLM,PRSD)
      DO IW=1,NW
        IF (NDAY.GT.0) THEN
           DO I=1,NDAY
             F00(I)=F0D(IW,I)
           ENDDO
           F0=AINTP(TJL,NDAY,TJLC,F00)
        ELSE
           F0=-999.
        ENDIF
        WL0=WL(IW)*1.0E4
        TAUR(IW)=RAYLS(WL0,PRS)
        TAUO3(IW)=AO3(WL0)*O3
        IF (F0.GT.0) THEN
           FDW=AUR(IW,1)/F0
           OPT(IW)=-LOG(FDW)/EM-TAUR(IW)-TAUO3(IW)
           IF (OPT(IW).LE.0.) OPT(IW)=0.
           IF (OPT(IW).GE.TAUMX) THEN
              ERC='STDTS5: illegal measured data(Too lagre OPT)'
              GOTO 900
           ENDIF
        ELSE
           OPT(IW)=-1.
        ENDIF
      ENDDO
C
  900 RETURN
      END
C
      SUBROUTINE REDML5(JTAU,IPLC,INVM,IDCR,IDCI,ANGMN,ANGMX,LOOPMX
     &     ,IUK,EPSA,EPSA1,EPSA2,GA,NLN,CONCA,CONCM,CR,CI,NW,NA
     &     ,NSIZE,RMIN,RMAX,WL,TAUR,TAUO3,TH0,TH,FI,SCA,OPT,AUR
     &     ,SOLID,LOOP,EROR,SRETV,SRETR,SRETI,IFLG
     &     ,TA,WA,OPTC,AURC,IAUR,SCAMX,JDCR,JDCI
     &     ,NANG,THETA,PHRTR,SIZE,SSL,VOL,WCR,WCI,IFLAG,ERC)
C
C --- Note
C   Skyradiance Analysis - Ver.5.0
C
C --- History
C   2007.10.05 Made by M.Yamano (from REDMLZ)
C   2008.01.18 Modified.
C   2008.10.17 Modified.
C
C --- Input
C JTAU     I           controlling optical thickness constraint
C                      1: AUR+OPT for inversion step by step
C                      0: only AUR for inversion (no use of OPT)
C                     -1: only OPT for inversion
C                      2: AUR+OPT for inversion all at once
C IPLC     I           polarization correction control
C                        1: WITH POLARIZATION CORRECTION
C                        0: WITHOUT CORRECTION
C INVM     I           =<19: box kernel
C                      >=20: log-normal multi mode kernel
C IDCR     I           0: give CR; 1: determine CR by this routine.
C IDCI     I           0: give CI; 1: determine CI by this routine.
C ANGMN    R           Minimum scattering angle for useful sky radiance
C ANGMX    R           Maximum scattering angle for useful sky radiance
C LOOPMX   I           >0  Maximum number of iteration
C                      otherwise : return after getting kernel
C IUK      I           DEVICE NO. FROM WHICH KERNELS ARE READ.
C EPSA     R           ABSOLUTE CONVERGENCE CRITERION.
C EPSA1    R           RELATIVE CONVERGENCE CRITERION.
C EPSA2    R           Give-up convergence criterion after 2nd-loop.
C GA       R(KNW)      ASSUMED GROUND ALBEDO FOR EACH CHANNEL.
C NLN      I           Number of Layer
C CONCA    R(KNLN)     Fraction of Aerosol   Optical Thickness
C CONCM    R(KNLN)     Fraction of Molecular Optical Thickness
C CR,CI    R(KNW)      COMPLEX REFRACTIVE INDEX OF AEROSOLS.
C                      ASSUME THAT SIZE DISTRIBUTION IS INDEPENDENT
C                      OF WAVELENGTHS
C                      Meaningful if IDCR=0 or IDCI=0
C NW       I           NBR OF WAVELENGTHS.
C NA       I           NBR OF ANGLES for AUR
C NSIZE    I           NBR OF PARTICLE RADII.
C RMIN     R           MINIMUM PARTICLE RADIUS IN CM
C RMAX     R           MAXIMUM PARTICLE RADIUS IN CM
C WL       R(KNW)      WAVELENGTH IN CM.
C TAUR     R(KNW)      Rayleigh optical thickness
C TAUO3    R(KNW)      Ozone optical thickness
C TH0      R           SOLAR ZENITH ANGLE IN DEGREE (0,90)
C TH       R(KNW,KDT)  EMERGENT ZENITH ANGLES IN DEGREE (0,90)
C FI       R(KNW,KDT)  AZIMUTHAL ANGLES IN DEGREE  (0,180)
C SCA      R(KNW,KDT)  Scattering angles of observed data
C OPT      R(KNW)      Observed aerosol optical thickness
C AUR      R(KNW,KDT)  I/EM/F/SOLID VIEW ANGLE.
C SOLID    R(KNW)      Solid view angles
C SIZE     R(KNSIZE)   mode radius (give when INVM>=20)
C SSL      R(KNSIZE)   log-dispersion
C                       dV/dlnr
C                         =sum(n=1,N) Cn*exp(-(ln(r/SIZE(n))/SSL(n))**2/2)
C --- Output
C LOOP     I           Number of iteration
C EROR     R           Reconstruction error for AURC
C SRETV    R(KNSIZE)   Retrieval reliability for dV/dlnr
C SRETR    R(KNW)      Retrieval reliability for Cr
C SRETI    R(KNW)      Retrieval reliability for Ci
C IFLG     I           1: convergence with criterion EPSA
C                      2: convergence with criterion EPSA1
C                      3: no convergence with criterion EPSA nor EPSA1
C                      0: convergence with criterion LOOPMX
C                      negative: no convergence with criterion EPSA2
C TA       R(KNW)      Retrieved optical thickness.
C WA       R(KNW)      Retrieved single scattering albedo.
C OPTC     R(KNW)      Reconstructed aerosol optical thickness
C AURC     R(KNW,KDT)  RECONSTRUCTED AUREOLE.
C IAUR     I(KNW,KDT)  availability flag for AUR(>0 then available)
C SCAMX    R           Upper limit of available scattering angles.
C NANG     I           NBR OF SCATTERING ANGLE ON THE KERNEL TABLE.
C THETA    R(KNANG)    SCATTERING ANGLES(IN DEGREE) ON THE KERNEL TABLE.
C PHRTR    R(KNMEAS)   Retrieved phase function
C VOL      R(KNSIZE)   Optimized volume spectrum(dV/dlnR)
C WCR      R(KNW)      Optimized real part of refractive index
C WCI      R(KNW)      Optimized imaginary part of refractive index
C IFLAG    I           Return flag( =  0: normal
C                                   = -1: illegal data(can't be analized)
C                                   = -9: fatal error
C ERC      C*64        ERROR CODE. IF ' ' THEN NORMAL.
C
C------------------------------------------- AREAS FOR THIS ROUTINE.
      SAVE
      PARAMETER (KCR   =20)
      PARAMETER (KCI   =20)
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
      PARAMETER (KNANG =60)
      PARAMETER (KNSIZE=60)
      PARAMETER (KINTVL=60)
      PARAMETER (KKSZ=KNSIZE+KNW*2)
      PARAMETER (KKDT=KDT*KNW+KNW)
      PARAMETER (KNA1  =50)
      PARAMETER (KNFI  =50)
      PARAMETER (KNLN  =4)
      PARAMETER (KNMEAS=KNANG*KNW+2*KNW,KNA1U=KNA1)
C
      PARAMETER (PI=3.141592653589793,RAD=PI/180.0)
      CHARACTER ERC*(*)
      DIMENSION CONCA(KNLN),CONCM(KNLN)
      DIMENSION WL(KNW),CR(KNW),CI(KNW),GA(KNW),SOLID(KNW)
     &         ,TAUR(KNW),TAUO3(KNW),OPT(KNW),AUR(KNW,KDT)
      DIMENSION WM(KNW),TA(KNW),WA(KNW),PHRTR(KNMEAS)
      DIMENSION SIZE(KNSIZE),SSL(KNSIZE),VOL(KNSIZE)
      DIMENSION VCR(KCR),VCI(KCI),AK(KNMEAS,KNSIZE,KCR,KCI)
     &         ,THETA(KNANG),PPM(KNANG),PRTR(KNMEAS)
      DIMENSION TH(KNW,KDT),FI(KNW,KDT),SCA(KNW,KDT),IAUR(KNW,KDT)
     &         ,ISTH(KNW,KDT),ISFI(KNW,KDT)
      DIMENSION AM1U(KNW,KNA1),FAI(KNW,KNFI),NA1U0(KNW),NFI0(KNW)
      DIMENSION AURC(KNW,KDT),OPTC(KNW),ASGM(KNW,KDT),FSGM(KNW)
      DIMENSION ANGTOP(KNW),ANGBTM(KNW),ARNG(3)
      DIMENSION JCR(KNW),JCI(KNW),JCR0(KNW),JCI0(KNW)
     &         ,WCR(KNW),WCI(KNW),WCR0(KNW),WCI0(KNW)
     &         ,AKNEL(KNMEAS,KNSIZE)
      DIMENSION SRETV(KNSIZE),SRETR(KNW),SRETI(KNW)
      DOUBLE PRECISION XX(KNANG)
C
      DATA INIT /0/
      DATA RANG /0.1/
      DATA ARNG / 30., 20., 70./
C
      ERC=' '
      IFLAG=0
C
C --- Initial call loop.
C
      IF (INIT.LE.0) THEN
C Read kernel
         CALL CKRNL8(IUK,INVM,SSL,NCR,NCI,VCR,VCI,NW
     &              ,WL,RMIN,RMAX,NSIZE,NANG,THETA,SIZE,AK,ERC)
         IF (ERC.NE.' ') GOTO 890
C
         IPOL=1
         DPF=0.03
         DPL=DPF/(2-DPF)
         DO J=1,NANG
           XX(J)=COS(THETA(J)*RAD)
           PPM(J)=3./(8.*PI*(2.+DPL))*(1.+DPL+(1.-DPL)*XX(J)**2)
         ENDDO
         INIT=1
      ENDIF
C
C --- Initial call loop end.
C
      CALL ININDX(IDCR,IDCI,NCR,NCI,VCR,VCI,NW,CR,CI,JCR0,JCI0)
c     WRITE(6,10)
c  10 FORMAT('Initial refractive indices :')
c     DO IW=1,NW
c       JCR0(IW)=ABS(JCR0(IW))
c       JCI0(IW)=ABS(JCI0(IW))
c       WRITE(6,15) WL(IW)*1.E4,VCR(JCR0(IW)),VCI(JCI0(IW))
c  15   FORMAT(' WL=',0PF6.3,' um: [ Cr= ',F6.3,' , Ci= ',F6.3,' ]')
c     ENDDO
C
      CALL ANGODR(NW,NA,TH,FI,ISTH,ISFI,NA1U0,AM1U,NFI0,FAI,ERC)
      IF (ERC.NE.' ') GOTO 800
C
      DO IW=1,NW
        WM(IW)=1.
      ENDDO
C
      SCAMAX=0.
      DO IW=1,NW
        IF (JTAU.GE.0) THEN
           IAUR(IW,1)=-9
           DO I=2,NA
             IF (AUR(IW,I).GT.0.) THEN
                IF (SCA(IW,I).GT.SCAMAX) SCAMAX=SCA(IW,I)
                IAUR(IW,I)=0
             ELSE
                IAUR(IW,I)=-9
             ENDIF
           ENDDO
        ENDIF
      ENDDO
C
C Set initial volume spectrum
C
      AMIN=ANGMN
      AMAX=ARNG(1)
      IF (AMIN.GT.AMAX) THEN
         JTAU0=-1
      ELSE
         JTAU0=JTAU
      ENDIF
      CALL ANGRNG(JTAU0,NW,NA,SCA,IAUR,ANGTOP,ANGBTM
     &           ,NANG,THETA,AMIN,AMAX,RANG,ERC)
      IF (ERC.NE.' ') GOTO 800
C
      DO IW=1,NW
        JCR(IW)=JCR0(IW)
        JCI(IW)=JCI0(IW)
        WCR(IW)=VCR(JCR(IW))
        WCI(IW)=VCI(JCI(IW))
      ENDDO
      CALL CLKRNL(NW,IPOL,NANG,NSIZE,NCR,NCI,VCR,VCI
     &                            ,JCR,JCI,WCR,WCI,AK,AKNEL)
      CALL INVOLM(JTAU0,IPLC,NLN,CONCA,CONCM,NW,WL,GA,SOLID
     &           ,TAUO3,TAUR,WM,TH0,NA,TH,FI,ISTH,ISFI
     &           ,OPT,AUR,IAUR,NA1U0,AM1U,NFI0,FAI,NSIZE,SIZE
     &           ,RMIN,RMAX,IPOL,NANG,THETA,AKNEL,PPM
     &           ,VOL,EPSA,EPSA1,IFLG,ERC)
      IF (IFLG.EQ.-9) GOTO 890
      IF (IFLG.LT.0) GOTO 800
C
C Loop for inversion
C
      JDCR=0
      JDCI=0
      IF (JTAU.LT.2) THEN
         DO IW=1,NW
           WCR0(IW)=WCR(IW)
           WCI0(IW)=WCI(IW)
         ENDDO
         IF (JTAU0.GE.0) JTAU0=0
         IDX=0
         CALL STSGMA(JTAU,NW,NA,IAUR,OPT,ASGM,FSGM)
         CALL INVLP5(JTAU0,IPLC,LOOPMX,EPSA,EPSA1,EPSA2
     &              ,NLN,CONCA,CONCM,NW,WL,GA,SOLID,TAUO3,TAUR,WM
     &              ,TH0,NA,TH,FI,ISTH,ISFI,OPT,AUR,IAUR,FSGM,ASGM
     &              ,NA1U0,AM1U,NFI0,FAI,NSIZE
     &              ,IPOL,NANG,THETA,AK,AKNEL,PPM,IDX
     &              ,WCR0,WCI0,NCR,NCI,VCR,VCI,JCR,JCI,WCR,WCI
     &              ,VOL,OPTC,AURC,PRTR,LOOP
     &              ,EROR,SRETV,SRETR,SRETI,IFLG,ERC)
         IF (IFLG.LT.0) THEN
            ERC=ERC(1:6)//'(0)'//ERC(7:61)
            GOTO 800
         ENDIF
      ENDIF
C
      IF (JTAU.GT.0) THEN
         CALL ANGRNG(JTAU,NW,NA,SCA,IAUR,ANGTOP,ANGBTM
     &              ,NANG,THETA,ANGMN,ANGMX,RANG,ERC)
         IF (ERC.NE.' ') GOTO 800
         DO IW=1,NW
           WCR0(IW)=WCR(IW)
           WCI0(IW)=WCI(IW)
         ENDDO
         IDX=2
         CALL STSGMA(JTAU,NW,NA,IAUR,OPT,ASGM,FSGM)
         CALL INVLP5(JTAU,IPLC,LOOPMX,EPSA,EPSA1,EPSA2
     &              ,NLN,CONCA,CONCM,NW,WL,GA,SOLID,TAUO3,TAUR,WM
     &              ,TH0,NA,TH,FI,ISTH,ISFI,OPT,AUR,IAUR,FSGM,ASGM
     &              ,NA1U0,AM1U,NFI0,FAI,NSIZE
     &              ,IPOL,NANG,THETA,AK,AKNEL,PPM,IDX
     &              ,WCR0,WCI0,NCR,NCI,VCR,VCI,JCR,JCI,WCR,WCI
     &              ,VOL,OPTC,AURC,PRTR,LOOP
     &              ,EROR,SRETV,SRETR,SRETI,IFLG,ERC)
         IF (IFLG.LT.0) THEN
            ERC=ERC(1:6)//'(1)'//ERC(7:61)
            GOTO 800
         ENDIF
         JDCR=JTAU
         JDCI=JTAU
      ENDIF
C
  500 NWLP=NW
      NWLE=NW
      NANG2=NANG*IPOL
      NMEAS2=NANG2*NWLP
      NMEAS1=NMEAS2+NWLE
      NMEAS=NMEAS1+NWLE
      SCAMX=0.
      DO IW=1,NW
        IF (ANGBTM(IW).GT.SCAMX) SCAMX=ANGBTM(IW)
        TA(IW)=PRTR(NMEAS2+IW)
        WA(IW)=PRTR(NMEAS1+IW)/PRTR(NMEAS2+IW)
        II=NANG2*(IW-1)
        DO I=1,NANG
          PHRTR(I+II)=PRTR(I+II)/PRTR(NMEAS1+IW)
        ENDDO
      ENDDO
      GOTO 900
C
  800 IFLAG=-1
      GOTO 900
C
  890 IFLAG=-9
  900 RETURN
      END
C
C ----------------------------------------------------------------------
      SUBROUTINE ANGODR(NW,NA,TH,FI,ISTH,ISFI,NA1U,AM1U,NFI,FAI,ERC)
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
      PARAMETER (KNA1  =50)
      PARAMETER (KNFI  =50)
C
      PARAMETER (PI=3.141592653589793,RAD=PI/180.0)
      CHARACTER ERC*(*)
      DIMENSION AM1U(KNW,KNA1),FAI(KNW,KNFI),NA1U(KNW),NFI(KNW)
      DIMENSION TH(KNW,KDT),FI(KNW,KDT),ISTH(KNW,KDT),ISFI(KNW,KDT)
      DIMENSION ANG(KDT),IBF(KDT)
C
C ORDERING zenith ANGLES.
C
      ERC=' '
      DO IW=1,NW
        DO J=1,NA
          ANG(J)=TH(IW,J)
        ENDDO
        CALL SORT1(NA,ANG,IBF)
        NA1U(IW)=1
        A1=ANG(IBF(1))
        A2=A1
        AM1U(IW,1)=COS(A1*RAD)
        ISTH(IW,IBF(1))=1
        DO I=2,NA
          A0=(A1+A2)*.5
          IF (A0.GT.20.) THEN
             DANG=0.5
          ELSE
             DANG=0.1
          ENDIF
          A2=ANG(IBF(I))
          IF (ABS(A1-A2).GT.DANG) THEN
             NA1U(IW)=NA1U(IW)+1
             IF (NA1U(IW).GT.KNA1) THEN
                ERC='ANGODR: EKNA1'
                GOTO 900
             ENDIF
             AM1U(IW,NA1U(IW))=COS(A2*RAD)
             A1=A2
          ENDIF
          ISTH(IW,IBF(I))=NA1U(IW)
        ENDDO
      ENDDO
C
C ORDERING azimuthal ANGLES.
C
      DO IW=1,NW
        DO J=1,NA
          ANG(J)=FI(IW,J)
        ENDDO
        CALL SORT1(NA,ANG,IBF)
        NFI(IW)=1
        A1=ANG(IBF(1))
        A2=A1
        FAI(IW,1)=A1
        ISFI(IW,IBF(1))=1
        DO I=2,NA
          A0=(A1+A2)*.5
          IF (A0.GT.20.) THEN
             DANG=0.5
          ELSE
             DANG=0.1
          ENDIF
          A2=ANG(IBF(I))
          IF (ABS(A1-A2).GT.DANG) THEN
             NFI(IW)=NFI(IW)+1
             IF (NFI(IW).GT.KNFI) THEN
                ERC='ANGODR: EKNFI'
                GOTO 900
             ENDIF
             FAI(IW,NFI(IW))=A2
             A1=A2
          ENDIF
          ISFI(IW,IBF(I))=NFI(IW)
        ENDDO
      ENDDO
C
  900 RETURN
      END
C
      SUBROUTINE ANGRNG(JTAU,NW,NA,SCA,IAUR,ANGTOP,ANGBTM
     &                 ,NANG,THETA,AMIN,AMAX,RANG,ERC)
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
      PARAMETER (KNANG =60)
C
      CHARACTER ERC*(*)
      DIMENSION IAUR(KNW,KDT),SCA(KNW,KDT),THETA(KNANG)
      DIMENSION ANGTOP(KNW),ANGBTM(KNW)
C
C Angular range of observation data
C
      ERC=' '
      DO IW=1,NW
        IF (JTAU.GE.0) THEN
           ANGTOP(IW)=180.
           ANGBTM(IW)=0.
           DO I=1,NA
             IF (ABS(IAUR(IW,I)).LE.1.) THEN
                IF ((SCA(IW,I).GE.AMIN-RANG)
     &                       .AND.(SCA(IW,I).LE.AMAX+RANG)) THEN
                   IAUR(IW,I)=1
                   IF (SCA(IW,I).LT.ANGTOP(IW)) ANGTOP(IW)=SCA(IW,I)
                   IF (SCA(IW,I).GT.ANGBTM(IW)) ANGBTM(IW)=SCA(IW,I)
                ELSE
                   IAUR(IW,I)=0
                ENDIF
             ENDIF
           ENDDO
        ELSE
           ANGTOP(IW)=0.
           ANGBTM(IW)=0.
        ENDIF
C
        J=0
        JTOP=0
        DO WHILE((J.LT.NANG).AND.(JTOP.EQ.0))
          J=J+1
          IF (THETA(J).GE.ANGTOP(IW)-RANG) JTOP=J
        ENDDO
        IF (JTOP.EQ.0) JTOP=NANG
        J=NANG+1
        JBTM=0
        DO WHILE((J.GT.1).AND.(JBTM.EQ.0))
          J=J-1
          IF (THETA(J).LE.ANGBTM(IW)+RANG) JBTM=J
        ENDDO
        IF (JBTM.EQ.0) JBTM=1
C
        IF (JTOP.GT.JBTM) THEN
           ERC='ANGRNG: illegal angular setting'
           GOTO 900
        ENDIF
      ENDDO
C
  900 RETURN
      END
C
      SUBROUTINE ININDX(IDCR,IDCI,NCR,NCI,VCR,VCI,NW,CR,CI,JCR,JCI)
      PARAMETER (KNW   =20)
      PARAMETER (KCR   =20)
      PARAMETER (KCI   =20)
C
      DIMENSION VCR(KCR),VCI(KCI),CR(KNW),CI(KNW),JCR(KNW),JCI(KNW)
C
      DO IW=1,NW
        ICR0=0
        ICR=0
        DO WHILE((ICR0.EQ.0).AND.(ICR.LT.NCR))
          ICR=ICR+1
          IF (ABS(VCR(ICR)/CR(IW)-1.).LE.1.E-5) THEN
             ICR0=ICR
          ELSE
             IF (VCR(ICR).GT.CR(IW)) ICR0=-ICR
          ENDIF
        ENDDO
        IF (ICR0.LT.0) THEN
           ICR0=ABS(ICR0)
           IF (ICR0.GT.1) THEN
              DR1=CR(IW)-VCR(ICR0-1)
              DR2=VCR(ICR0)-CR(IW)
              IF (DR1.LE.DR2) ICR0=ICR0-1
           ENDIF
        ELSE IF (ICR0.EQ.0) THEN
           ICR0=NCR
        ENDIF
        JCR(IW)=ICR0
C
        ICI0=0
        ICI=0
        DO WHILE((ICI0.EQ.0).AND.(ICI.LT.NCI))
          ICI=ICI+1
          IF (CI(IW).EQ.0.) THEN
             IF (ABS(VCI(ICI)).EQ.0.) ICI0=ICI
          ELSE
             IF (ABS(VCI(ICI)/CI(IW)-1.).LE.1.E-5) THEN
                ICI0=ICI
             ELSE
                IF (ABS(VCI(ICI)).GT.ABS(CI(IW))) ICI0=-ICI
             ENDIF
          ENDIF
        ENDDO
        IF (ICI0.LT.0) THEN
           ICI0=ABS(ICI0)
           IF (ICI0.GT.1) THEN
              DI1=ABS(CI(IW))-ABS(VCI(ICI0-1))
              DI2=ABS(VCI(ICI0))-ABS(CI(IW))
              IF (DI1.LE.DI2) ICI0=ICI0-1
           ENDIF
        ELSE IF (ICI0.EQ.0) THEN
           ICI0=NCI
        ENDIF
        JCI(IW)=ICI0
      ENDDO
C
      RETURN
      END
C
      SUBROUTINE STSGMA(JTAU,NW,NA,IAUR,OPT,ASGM,FSGM)
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
C
      DIMENSION IAUR(KNW,KDT),OPT(KNW),ASGM(KNW,KDT),FSGM(KNW)
C
      DATA SGA,SGF1,SGF0,SGV /0.05, 0.01, 0.05, 0.04/
C
      DO IW=1,NW
        IF (OPT(IW).GT.0.) THEN
           SGF=SGF1/OPT(IW)
        ELSE
           SGF=SGF0
        ENDIF
        IF (JTAU.GE.0) THEN
           DO I=1,NA
             ASGM(IW,I)=SGA+SGF+SGV
           ENDDO
        ELSE
           DO I=1,NA
             ASGM(IW,I)=1.
           ENDDO
        ENDIF
        FSGM(IW)=1.
        IF (ABS(JTAU).GT.0) THEN
           IF (OPT(IW).GT.0.) FSGM(IW)=SGF
        ENDIF
      ENDDO
C
      RETURN
      END
C
      SUBROUTINE VOLLOG(NMODE,RMODE,SSL,CMODE,NSIZE,SIZE,VOL)
C
C   Calculation of log-normal multi mode size distribution.
C
C --- history
C   2000.12.20 Created by M.Yamano
C
C --- Input
C NMODE    I         number of mode radii
C RMODE    R(NMODE)  mode radius
C SSL      R(NMODE)  log-dispersion
C CMODE    R(NMODE)  peak volume Cn in the definition
C                    dV/dlnR=sum(1,n) VOL(n)*exp(-(ln(R/size(n))/SSL(n))**2/2)
C NSIZE    I         number of particle radii.
C SIZE     R(NSIZE)  particle radii
C
C --- Output
C VOL      R(NSIZE)  dV/dlnR
C
C ----------------------------------------------------------------------
C
      DIMENSION RMODE(NMODE),SSL(NMODE),CMODE(NMODE)
     &         ,SIZE(NSIZE),VOL(NSIZE)
C
      DO I=1,NSIZE
        S=0.
        DO J=1,NMODE
          S=S+CMODE(J)*EXP(-(LOG(SIZE(I)/RMODE(J))/SSL(J))**2/2.0)
        ENDDO
        VOL(I)=S
      ENDDO
      RETURN
      END
C
      SUBROUTINE INTERV(INTVL,RDS,VOL,NSIZE,SIZE,VOLS)
C
C   Linear interpolation of volume spectrum.
C
C --- history
C   2001.01.17 Created by M.Yamano
C
C --- Input
C INTVL  I             number of radii for given volume spectrum.
C RDS    R(INTVL)      radii data for given volume spectrum.
C VOL    R(INTVL)      volume data for given volume spectrum.
C NSIZE  I             number of radii for interpolated volume spectrum.
C SIZE   R(NSIZE)      radii data for interpolated volume spectrum.
C
C --- Output
C VOLS   R(NSIZE)      volume data for interpolated volume spectrum.
C ----------------------------------------------------------------------
      PARAMETER (KINTVL=60)
C
      PARAMETER (PI=3.141592653589793,RAD=PI/180.0)
      DIMENSION RDS(KINTVL),VOL(KINTVL),SIZE(KINTVL),VOLS(KINTVL)
      DIMENSION XX(KINTVL),YY(KINTVL)
C
      DO I=1,INTVL
        XX(I)=LOG(RDS(I))
        IF (VOL(I).LE.0.) VOL(I)=1.E-20
        YY(I)=LOG(VOL(I))
      ENDDO
C
      I2=1
      N=0
      DO WHILE(N.LT.NSIZE)
        N=N+1
        X1=LOG(SIZE(N))
        I=I2-1
        I2=0
        DO WHILE((I.LT.INTVL).AND.(I2.LE.0))
          I=I+1
          IF (XX(I).GE.X1) I2=I
        ENDDO
        IF (I2.GT.1) THEN
           Y1=(YY(I2)-YY(I2-1))*(X1-XX(I2-1))/(XX(I2)-XX(I2-1))+YY(I2-1)
           VOLS(N)=EXP(Y1)
        ELSE
           IF (ABS(X1-XX(I2)).LT.1.E-8) THEN
              VOLS(N)=EXP(YY(I2))
           ELSE
              VOLS(N)=0.
           ENDIF
        ENDIF
      ENDDO
      RETURN
      END
C
      SUBROUTINE WTOUT5(IU,HEAD,INVM,NSIZE,SIZE,VOL,NMODE,RMODE,CMODE
     &                 ,WCR,WCI,NW,WL,OPT,TA,WA,NA,TH,FI,SCA,IAUR
     &                 ,AUR,AURC,NANG,THETA,POBSN,JTAU
     &                 ,SRETV,SRETR,SRETI)
C
C   Output results(full) of analysis Ver.5.0.
C
C --- history
C   2007.10.31 Created by M.Yamano
C
C --- Input
C
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
      PARAMETER (KNANG =60)
      PARAMETER (KNSIZE=60)
C 20240529 Modified: from Momoi-san, 2018.07.16
C      PARAMETER (KNMEAS=1240)
      PARAMETER (KNMEAS=KNANG*KNW+2*KNW)
C
      CHARACTER HEAD*80
      DIMENSION WCR(KNW),WCI(KNW)
      DIMENSION SIZE(KNSIZE),VOL(KNSIZE)
     &         ,RMODE(KNSIZE),CMODE(KNSIZE)
      DIMENSION WL(KNW),OPT(KNW),TA(KNW),WA(KNW)
      DIMENSION TH(KNW,KDT),FI(KNW,KDT),SCA(KNW,KDT),IAUR(KNW,KDT)
     &         ,AUR(KNW,KDT),AURC(KNW,KDT)
      DIMENSION THETA(KNANG),POBSN(KNMEAS)
      DIMENSION SRETV(KNSIZE),SRETR(KNW),SRETI(KNW)
C
      WRITE(IU,10) HEAD
   10 FORMAT(A80//' Refractive Indices')
      WRITE(IU,11) (WL(IW)*1.E4, IW=1,NW)
   11 FORMAT('  WL ',0P7F10.4)
      WRITE(IU,12) (WCR(IW),IW=1,NW)
   12 FORMAT('  Cr ',0P7F10.4)
C      WRITE(IU,13) (SRETR(IW),IW=1,NW)
C   13 FORMAT('  R(r)',0P7F10.3)
      WRITE(IU,14) (WCI(IW),IW=1,NW)
   14 FORMAT('  Ci ',0P7F10.5)
C      WRITE(IU,15) (SRETI(IW),IW=1,NW)
C   15 FORMAT('  R(i)',0P7F10.3)
C
CCC 20110831 hashimoto
       WRITE(IU,20)
   20  FORMAT(/ '     Radius     Volume')
       DO J=1,NSIZE
         WRITE(IU,25) SIZE(J),VOL(J)
   25    FORMAT(1P2E12.3)
       ENDDO
CCC 20110831 hashimoto end
C      WRITE(IU,20)
C   20 FORMAT(/ '  No.     Radius         Volume        R(v)')
C      DO J=1,NSIZE
C        WRITE(IU,25) J,SIZE(J),VOL(J),SRETV(J)
C   25   FORMAT(I4,1P2E15.4,0PF10.3)
C      ENDDO
CCC -- hashimoto --
C
      WRITE(IU,40)
   40 FORMAT(/' Cross sections')
      WRITE(IU,42) (WL(IW)*1.E4, IW=1,NW)
   42 FORMAT('  WL ',0P7F10.4)
      WRITE(IU,44) (OPT(IW),IW=1,NW)
   44 FORMAT('  OPT',0P7F10.4)
      WRITE(IU,46) (TA(IW),IW=1,NW)
   46 FORMAT('  TA ',0P7F10.4)
      WRITE(IU,48) (WA(IW),IW=1,NW)
   48 FORMAT('  WA ',0P7F10.4)
C
      IF (JTAU.GE.0) THEN
         WRITE(IU,50)
   50    FORMAT(/' fg    TH     FI     SCA     AUR        AURC')
         DO IW=1,NW
           DO I=1,NA
             WRITE(IU,55) IAUR(IW,I),TH(IW,I),FI(IW,I),SCA(IW,I)
     &                    ,AUR(IW,I),AURC(IW,I)
   55        FORMAT(I3,0P3F7.1,1P2E12.3)
           ENDDO
         ENDDO
      ENDIF
C
      WRITE(IU,60)
   60 FORMAT(/' THETA POBSN(IW=1,NW)')
      DO I=1,NANG
        WRITE(IU,65) THETA(I),(POBSN((IW-1)*NANG+I),IW=1,NW)
   65   FORMAT(F6.1,1P7E10.3)
      ENDDO
      RETURN
      END
C
      SUBROUTINE WTPAR5(IU,HEAD,NW,WL,OPT,TA,WA,WCR,WCI,SRETR,SRETI)
C
C   Output results(optical parameters) of analysis Ver.5.0.
C
C --- history
C   2007.10.31 Created by M.Yamano
C
C --- Input
C
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
C
      CHARACTER HEAD*80
      DIMENSION WCR(KNW),WCI(KNW)
      DIMENSION WL(KNW),OPT(KNW),TA(KNW),WA(KNW)
      DIMENSION SRETR(KNW),SRETI(KNW)
C
      WRITE(IU,10) HEAD
   10 FORMAT(A80)
CCC 20110831 hashimoto
      DO IW=1,NW
        WRITE(IU,20) WL(IW)*1.E4,OPT(IW),TA(IW),WA(IW)
     &              ,WCR(IW),WCI(IW)
   20   FORMAT(0P4F10.4,F10.4,F10.5)
      ENDDO
CCC 20110831 hashimoto end
C      DO IW=1,NW
C        WRITE(IU,20) WL(IW)*1.E4,OPT(IW),TA(IW),WA(IW)
C     &              ,WCR(IW),SRETR(IW),WCI(IW),SRETI(IW)
C   20   FORMAT(0P4F10.4,F10.4,F8.3,F10.5,F8.3)
C      ENDDO
CCC -- hashimoto --
      RETURN
      END
C
      SUBROUTINE WTVOL5(IU,HEAD,NSIZE,SIZE,VOL,SRETV)
C
C   Output results(volume spectrum) of analysis Ver.5.0.
C
C --- history
C   2007.10.31 Created by M.Yamano
C
C --- Input
C
C ----------------------------------------------------------------------
      PARAMETER (KNSIZE=60)
C
      CHARACTER HEAD*80
      DIMENSION SIZE(KNSIZE),VOL(KNSIZE),SRETV(KNSIZE)
C
      WRITE(IU,10) HEAD
   10 FORMAT(A80)
      DO J=1,NSIZE
CCC 20110831 hashimoto 
        WRITE(IU,20) SIZE(J),VOL(J)
   20   FORMAT(1P2E12.3)
CCC 20110831 hashimoto end
CCC        WRITE(IU,20) J,SIZE(J),VOL(J),SRETV(J)
CCC   20   FORMAT(I5,1P2E15.4,0PF10.3)
      ENDDO
      RETURN
      END
C
      SUBROUTINE WTAUR4(IU,HEAD,NW,NA,TH,FI,SCA,IAUR,AUR,AURC)
C
C   Output results(sky radiance) of analysis Ver.4.0.
C
C --- history
C   2001.01.19 Created by M.Yamano
C   2001.08.22 Modified.
C   2003.10.29 Changed.: IAUR(I) -> IAUR(IW,I)
C
C --- Input
C
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
C
      CHARACTER HEAD*80
      DIMENSION TH(KNW,KDT),FI(KNW,KDT),SCA(KNW,KDT),IAUR(KNW,KDT)
     &         ,AUR(KNW,KDT),AURC(KNW,KDT)
C
      WRITE(IU,10) HEAD
   10 FORMAT(A80)
      DO IW=1,NW
        DO I=1,NA
          WRITE(IU,20) IAUR(IW,I),TH(IW,I),FI(IW,I),SCA(IW,I)
     &                 ,AUR(IW,I),AURC(IW,I)
   20     FORMAT(I3,0P3F7.1,1P2E12.3)
        ENDDO
      ENDDO
      RETURN
      END
C
      SUBROUTINE WTPHS0(IU,HEAD,NW,NANG,THETA,POBSN)
C
C   Output results(phase function) of analysis Ver.3.0/4.0.
C
C --- history
C   2001.01.19 Created by M.Yamano
C   2001.08.22 Modified.
C   2001.09.05 Renamed. (WTPHS4 -> WTPHS0)
C
C --- Input
C
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KNANG =60)
C 20240529 Modified: from Momoi-san, 2018
C      PARAMETER (KNMEAS=1240)
      PARAMETER (KNMEAS=KNANG*KNW+2*KNW)
C
      CHARACTER HEAD*80
      DIMENSION THETA(KNANG),POBSN(KNMEAS)
C
      WRITE(IU,10) HEAD
   10 FORMAT(A80)
      DO I=1,NANG
        WRITE(IU,20) THETA(I),(POBSN((IW-1)*NANG+I),IW=1,NW)
   20   FORMAT(F6.1,1P7E10.3)
      ENDDO
      RETURN
      END
C
      SUBROUTINE WTDATF(IU,NN,IY,IM,ID,TIM,EM,EROR,NO3,NWV
     &                 ,NW,WL,TAUR,TAUO3,AUR,TA,WA)
C
C   Create 'f0_wl?.dat' for analysis Ver.5.0.
C
C --- history
C   2004.01.21 Created from old version for Ver.4.2.
C   2006.04.25 Changed writing format for air mass (F7.4 -> F8.4)
C   2024.01.17 Edited for Ver.5.0.
C
C --- Input
C
C --- Output
C ERC           C*64      ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
C
      DIMENSION WL(KNW),TAUR(KNW),TAUO3(KNW),AUR(KNW,KDT)
      DIMENSION TA(KNW),WA(KNW)
C
      IF (NN.EQ.1) THEN
         DO IW=1,NW
           WRITE(IU+IW,10) WL(IW),IW,NO3,NWV
   10      FORMAT(1PE11.3,3I4,' : WL(IW)  IW  IO3  IWV')
           WRITE(IU+IW,15)
   15      FORMAT('   Y   M  D   HR      M      FOBS   TAUR   TAUO3'
     &           ,'    TAUE      TAUS      ERR')
         ENDDO
      ENDIF
C
      DO IW=1,NW
        WRITE(IU+IW,20) IY,IM,ID,TIM,EM,AUR(IW,1)
     &                 ,TAUR(IW),TAUO3(IW),TA(IW),WA(IW)*TA(IW),EROR
   20   FORMAT(I5,2I3,F7.3,F8.4,1PE10.3,0P2F7.4,1P3E10.3)
      ENDDO
      RETURN
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
      SUBROUTINE RDINSP(IU,SRNO,ITYP,NW,WL,SOLID,NDY,F0D,TJLC,ERC)
C
C   Read instrumental parameters.(from 'ins.para' file)
C
C --- History
C   2001.02.14 Created by M.Yamano
C   2002.05.11 Renewed by M.Yamano for Ver.5.
C   2002.11.20 Modified.
C   2003.11.05 Added processing for ITYP=30 data.
C   2007.04.29 Added processing for ITYP=31 data.
C
C --- Input
C IU     I           device No. for reading
C
C --- Output
C SRNO   C*20        instrument S/N
C ITYP   I           data file format type No.
C                     10,11   : PREDE - skyradiometer POM-01L
C                     20,21,22: PREDE - skyradiometer POM-01MKII
C                     30,31   : PREDE - skyradiometer POM-02
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
     &                     (ITYP.NE.30).AND.(ITYP.NE.31)) THEN
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
      SUBROUTINE RDMETP(IU,CMNT,SITE,NDY,PRSD,O3D,TJLC,ERC)
C
C   Read meteorological parameters.(from 'METEO.DAT' file)
C
C --- History
C   2000.11.29 Created by M.Yamano
C   2002.05.15 Renewed by M.Yamano for Ver.5.
C   2006.04.04 Added check of NDY (if NDY<=0 -> error)
C
C --- Input
C IU     I           device No. for reading
C
C --- Output
C CMNT   C*40        comment
C SITE   C*20        observation site name
C NDY    I           number of meteorological data
C PRSD   R(NDY)      pressure data[atm]
C O3D    R(NDY)      ozone data[cm,STP]
C TJLC   R(NDY)      time data(Julian day)
C ERC    C*64        ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KDAY  =30)
C
      CHARACTER ERC*(*)
      CHARACTER CMNT*(*),SITE*(*)
      DIMENSION PRSD(KDAY),O3D(KDAY),TJLC(KDAY)
C
      ERC=' '
      READ(IU,*)
      READ(IU,*,ERR=800,END=800) CMNT
      READ(IU,*,ERR=800,END=800) SITE
      READ(IU,*,ERR=800,END=800) ALNGS
      READ(IU,*,ERR=800,END=800) ALNG,ALAT,ALT
      READ(IU,*,ERR=800,END=800) NDY
      IF (NDY.GT.KDAY) THEN
         ERC='RDMETP: NDY.GT.KDAY !'
         GOTO 900
      ENDIF
      IF (NDY.LE.0) THEN
         ERC='RDMETP: illegal NDY !'
         GOTO 900
      ENDIF
      DO I=1,NDY
        READ(IU,*,ERR=800,END=800) IYC,IMC,IDC,TMC,PRSD(I),O3D(I)
        TJLC(I)=TJUL1(IYC,IMC,IDC,TMC,ALNGS)
      ENDDO
      GOTO 900
C
  800 ERC='RDMETP: Read error.'
C
  900 RETURN
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
      FUNCTION AINTP(XX,N,X,Y)
C LINEAR INTERPOLATION
C--- HISTORY
C 90. 3. 6 CREATED
C--- INPUT
C XX    R     INTERPOLATION POINT
C N     R     NBR OF DATA
C X   R(NN)   INDEPENDENT VARIABLE,  X(I) .LT. X(I+1)
C Y   R(NN)   DEPENDENT   VARIABLE
C--- OUTPUT
C AINTP R     OUTPUT
C$ENDI
C
      DIMENSION X(N),Y(N)
C 1 POINT
      IF(N.LE.1) THEN
      AINTP=Y(1)
      RETURN
      ENDIF
C 2 POINTS OR XX.LE.X(1)
      IF(XX.LE.X(1) .OR. N.EQ.2) THEN
      AINTP=Y(1)+(Y(2)-Y(1))*(XX-X(1))/(X(2)-X(1))
      RETURN
      ENDIF
C XX.GE.XX(N)
      IF(XX.GE.X(N)) THEN
      AINTP=Y(N-1)+(Y(N)-Y(N-1))*(XX-X(N-1))/(X(N)-X(N-1))
      RETURN
      ENDIF
C MORE THAN 2 POINTS
      DO 1 I=1,N
      IF(XX.LE.X(I)) GOTO 2
    1 CONTINUE
      I=N
    2 K1=MAX(1,I-1)
      K2=K1+1
      IF(K2.GT.N) THEN
        K2=N
        K1=K2-1
      ENDIF
C
      AINTP=Y(K1)+(Y(K2)-Y(K1))*(XX-X(K1))/(X(K2)-X(K1))
      RETURN
      END
      FUNCTION RAYLS(WL,PRS)
C RAYLEIGH SCATTERING OPTICAL THICKNESS
C--- HISTORY
C 89.12. 7 CREATED
C--- INPUT
C WL     R      WAVELENGTH IN MICRON
C PRS    R      PRESSURE IN ATM
C--- OUTPUT
C RAYLS  R      OPTICAL THICKNESS
C$ENDI
      RAYLS=0.00864/WL**(3.916+0.074*WL+0.05/WL)*PRS
      RETURN
      END

      function AO3(WL)
C Absorption coefficient of ozone
C--- history
C 98. 6.21 Created
C--- input
C WL      R     Wavelength (micron) (0.3, 4) micron
C--- output
C AO3     R     Absorption coefficient of ozone (cm-1 STP)
C
      parameter (KD=22)
      dimension WLD(KD),AO3D(KD)
      data WLD/0.30,0.31,0.32,0.34,0.36,0.44,0.46,0.48,0.50,0.52
     & ,0.54,0.56,0.58,0.60,0.62,0.64,0.66,0.68,0.70,0.72,0.75,4.00/
      data AO3D/10.00,2.70,0.80,0.04,0.00,0.00,0.01,0.01,0.03,0.05
     & ,0.08,0.10,0.12,0.13,0.11,0.08,0.06,0.04,0.02,0.01,0.00,0.00/
C
      do 1 K=1,KD-1
      if((WL-WLD(K))*(WL-WLD(K+1)) .LE. 0.0) goto 2
    1 continue
      AO3=-1.0E10
      return
    2 if(WL.GE.0.36) goto 3
      if(WL.GT.0.34) K=3
      AO3=EXP(LOG(AO3D(K))+LOG(AO3D(K+1)/AO3D(K))/(WLD(K+1)-WLD(K))
     &  *(WL-WLD(K)))
      return
    3 AO3=AO3D(K)+(AO3D(K+1)-AO3D(K))/(WLD(K+1)-WLD(K))*(WL-WLD(K))
      return
      end
      SUBROUTINE INVLP5(JTAU,IPLC,LOOPMX,EPSA,EPSA1,EPSA2
     &          ,NLN,CONCA,CONCM,NW,WL,GA,SOLID,TAUO3,TM,WM
     &          ,TH0,NA,TH,FI,ISTH,ISFI,OPT,AUR,IAUR,FSGM,ASGM
     &          ,NA1U0,AM1U,NFI0,FAI,NSIZE
     &          ,IPOL,NANG,THETA,AK,AKNEL,PPM,IDX
     &          ,WCR,WCI,NCR,NCI,VCR,VCI,JCR,JCI,WCR0,WCI0
     &          ,VOL0,OPTC0,AURC0,PRTR0,LOOP
     &          ,EROR,SRETV,SRETR,SRETI,IFLG,ERC)
C
C --- Note
C   Determination of volume spectrum for skyradiance analysis.
C
C --- History
C   2007.10.05 Made by M.Yamano (from INVLOP)
C   2008.01.18 Modified.
C   2008.10.17 Modified.
C
C --- Input
C JTAU     I            controlling optical thickness constraint
C                        1: OPT is observed and AUR+OPT for inversion.
C                        0: only AUR for inversion (no use of OPT)
C                       -1: only OPT for inversion
C IPLC     I            polarization correction control
C                        1: WITH POLARIZATION CORRECTION
C                        0: WITHOUT CORRECTION
C LOOPMX   I            Maximum number of iteration
C EPSA     R            ABSOLUTE CONVERGENCE CRITERION.
C EPSA1    R            RELATIVE CONVERGENCE CRITERION.
C EPSA2    R            Give-up convergence criterion after 2nd-loop.
C NLN      I            Number of Layer
C CONCA    R(KNLN)      Fraction of Aerosol   Optical Thickness
C CONCM    R(KNLN)      Fraction of Molecular Optical Thickness
C NW       I            NBR OF WAVELENGTHS.
C WL       R(KNW)       WAVELENGTH IN CM.
C GA       R(KNW)       ASSUMED GROUND ALBEDO FOR EACH CHANNEL.
C SOLID    R(KNW)       Solid view angles
C TAUO3    R(KNW)       Ozone optical thickness
C TM       R(KNW)       Molecular Optical Thickness
C WM       R(KNW)       Molecular single scattering albedo
C TH0      R            SOLAR ZENITH ANGLE IN DEGREE (0,90)
C NA       I            NBR OF ANGLES for AUR
C TH       R(KNW,KDT)   EMERGENT ZENITH ANGLES IN DEGREE (0,90)
C FI       R(KNW,KDT)   AZIMUTHAL ANGLES IN DEGREE  (0,180)
C ISTH     I(KNW,KDT)   Ordering No. of zenith angles
C ISFI     I(KNW,KDT)   Ordering No. of azimuthal angles
C OPT      R(KNW)       Observed aerosol optical thickness
C AUR      R(KNW,KDT)   I/EM/F/SOLID VIEW ANGLE.
C IAUR     I(KNW,KDT)   availability flag for AUR(>0 then available)
C FSGM     R(KNW)       sigma of OPT for inversion
C ASGM     R(KNW,KDT)   sigma of AUR for inversion
C NA1U0    I(KNW)       NUMBER OF EMERGENT NADIR ANGLES IN SPHERE.
C AM1U     R(KNW,KNA1U) CONSINE OF THE EMERGENT NADIR ANGLES.
C                          - FOR UPWARD, + FOR DOWNWARD.
C NFI0     I(KNW)       NUMBER OF AZIMUTHAL ANGLES.
C FAI      R(KNW,KNFI)  AZIMUTHAL ANGLES IN DEGREES.
C NSIZE    I            NBR OF PARTICLE RADII.
C IPOL     I            fixed(=1) for CKRNL8
C NANG     I            NBR OF SCATTERING ANGLE ON THE KERNEL TABLE.
C THETA    R(KNANG)     SCATTERING ANGLES(IN DEGREE) ON THE KERNEL TABLE.
C AK       R(KNMEAS,KNSIZE,KCR,KCI)
C AKNEL    R(KNMEAS,KNSIZE)
C                       KERNEL MATRIX:  (I,J) ,I=1,NMEAS, J=1,NSIZE
C PPM      R(KNANG)     Phase function of Rayleigh scattering
C IDX      I            inversion flag
C                        0: dV/dlnr
C                        1: dV/dlnr + Cr
C                       -1: dV/dlnr + Ci
C                        2: dV/dlnr + Cr + Ci
C WCR      R(KNW)       real part value of current refractive index
C WCI      R(KNW)       imaginary part value of current refractive index
C NCR      I            real part number of kernel matrix
C NCI      I            imaginary part number of kernel matrix
C VCR      R(KCR)       real part data array of kernel matrix
C VCI      R(KCI)       imaginary part data array of kernel matrix
C JCR      I(KNW)       real part index of current refractive index
C JCI      I(KNW)       imaginary part index od current refractive index
C --- Output
C LOOP     I            Number of iteration
C EROR     R            Reconstruction error in AURC
C IFLG     I            1: convergence with criterion EPSA.
C                       2: convergence with criterion EPSA1.
C                       not positive: no convergence
C OPTC0    R(KNW)       Retrieved optical thickness
C AURC0    R(KNW,KDT)   RECONSTRUCTED AUREOLE.
C IAUR     I(KNW,KDT)   availability flag for AUR(>0 then available)
C PRTR0    R(KNMEAS)    RETRIEVED MIE PHASE FUNCTION AFTER REDUCING MULTI-SCAT.
C VOL0     R(KNSIZE)    Optimized volume spectrum(dV/dlnR)
C WCR0     R(KNW)       Optimized real part of refractive index
C WCI0     R(KNW)       Optimized imaginary part of refractive index
C ERC      C*64         ERROR CODE. IF ' ' THEN NORMAL.
C
C------------------------------------------- AREAS FOR THIS ROUTINE.
      SAVE
      PARAMETER (KCR   =20)
      PARAMETER (KCI   =20)
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
      PARAMETER (KNANG =60)
      PARAMETER (KNSIZE=60)
      PARAMETER (KKSZ=KNSIZE+KNW*2)
      PARAMETER (KKDT=KDT*KNW+KNW)
      PARAMETER (KNA1  =50)
      PARAMETER (KNFI  =50)
      PARAMETER (KNLN  =4)
      PARAMETER (KNMEAS=KNANG*KNW+2*KNW,KNA1U=KNA1)
C
      PARAMETER (PI=3.141592653589793,RAD=PI/180.0)
      CHARACTER ERC*(*)
      DIMENSION CONCA(KNLN),CONCM(KNLN)
      DIMENSION WL(KNW),GA(KNW),SOLID(KNW),TAUO3(KNW)
     &         ,OPT(KNW),AUR(KNW,KDT),IAUR(KNW,KDT)
      DIMENSION TM(KNW),WM(KNW)
      DIMENSION VOL0(KNSIZE),VOL(KNSIZE)
      DIMENSION VCR(KCR),VCI(KCI),AK(KNMEAS,KNSIZE,KCR,KCI)
     &         ,THETA(KNANG),PPM(KNANG),PRTR(KNMEAS),PRTR0(KNMEAS)
      DIMENSION TH(KNW,KDT),FI(KNW,KDT),ISTH(KNW,KDT),ISFI(KNW,KDT)
      DIMENSION AM1U(KNW,KNA1),FAI(KNW,KNFI),NA1U0(KNW),NFI0(KNW)
      DIMENSION AURC(KNW,KDT),AURC0(KNW,KDT),OPTC(KNW),OPTC0(KNW)
     &         ,ASGM(KNW,KDT),FSGM(KNW)
      DIMENSION ERRV(KNW)
      DIMENSION DLA(KKDT),DLAC(KKDT,KKSZ),DEV(KKSZ),DEV0(KKSZ),IEF(KKDT)
     &         ,DDV(KKSZ),SIG(KKSZ)
      DIMENSION SRET(KKSZ),SRETV(KNSIZE),SRETR(KNW),SRETI(KNW)
      DIMENSION JCR(KNW),JCI(KNW),WCR(KNW),WCI(KNW),WCR0(KNW),WCI0(KNW)
     &         ,AKNEL(KNMEAS,KNSIZE)
      DIMENSION WGTM(KKDT),WGTP(KKSZ),DEVA(KKSZ),DAS(KKSZ)
C
      DATA COEF,VRATE /10., 0.01/
C
      ERC=' '
      DO IW=1,NW
        WCR0(IW)=WCR(IW)
        WCI0(IW)=WCI(IW)
      ENDDO
      CALL VLRTRN(IPLC,NLN,CONCA,CONCM,NW,WL,GA,SOLID
     &           ,TAUO3,TM,WM,TH0,NA,TH,FI,ISTH,ISFI,OPT,AUR,IAUR
     &           ,NA1U0,AM1U,NFI0,FAI,NSIZE,VOL0
     &           ,IPOL,NANG,THETA,AKNEL,PPM,OPTC0,AURC0,PRTR,ERC)
      IF (ERC.NE.' ') GOTO 890
C
      CALL SETWGT(JTAU,NW,NA,ASGM,FSGM,MA,WGTM,IDX,NSIZE,NINV,WGTP)
      CALL SETDEV(NSIZE,VOL0,NW,WCR0,WCI0,IDX,DEVA,DDV,VRATE)
C
      ERRB=1.E10
      XJ0=ERRB
      IFLG=0
      LOOP=0
      DO WHILE((IFLG.EQ.0).AND.(LOOP.LT.LOOPMX))
        LOOP=LOOP+1
        CALL SETDEV(NSIZE,VOL0,NW,WCR0,WCI0,IDX,DEV0,DDV,VRATE)
        CALL CALDLT(JTAU,NW,NA,NINV,WGTP,DEV0,DEVA,DAS
     &             ,MA,WGTM,IAUR,AUR,AURC0,OPT,OPTC0,DLA,IEF,XJ1)
C
        DO IR=1,NINV
          DO J=1,NINV
            IF (J.EQ.IR) THEN
               DEV(J)=DEV0(J)+DDV(J)
            ELSE
               DEV(J)=DEV0(J)
            ENDIF
          ENDDO
          CALL GETDEV(DEV,NSIZE,VOL,NW,WCR,WCI
     &                         ,NCR,NCI,VCR,VCI,JCR,JCI,IDX)
          CALL CHKDEV(IR,NSIZE,NW,IDX,WCR,WCR0,WCI,WCI0
     &                         ,NCR,NCI,VCR,VCI,JCR,JCI)
          CALL CLKRNL(NW,IPOL,NANG,NSIZE,NCR,NCI,VCR,VCI
     &                         ,JCR,JCI,WCR,WCI,AK,AKNEL)
          CALL VLRTRN(IPLC,NLN,CONCA,CONCM,NW,WL,GA,SOLID
     &            ,TAUO3,TM,WM,TH0,NA,TH,FI,ISTH,ISFI,OPT,AUR,IAUR
     &            ,NA1U0,AM1U,NFI0,FAI,NSIZE,VOL
     &            ,IPOL,NANG,THETA,AKNEL,PPM,OPTC,AURC,PRTR,ERC)
          IF (ERC.NE.' ') GOTO 890
C
          DO IW=1,NW
            IF (JTAU.GE.0) THEN
               K=NA*(IW-1)
               DO I=1,NA
                 DLAC(K+I,IR)=(LOG(AURC(IW,I))-LOG(AURC0(IW,I)))
     &                                                   /DDV(IR)
               ENDDO
               K=NA*NW
            ELSE
               K=0
            ENDIF
            IF (ABS(JTAU).GT.0)
     &         DLAC(K+IW,IR)=(LOG(OPTC(IW))-LOG(OPTC0(IW)))/DDV(IR)
          ENDDO
        ENDDO
C
        CALL INVMAP(MA,NINV,WGTM,WGTP,DLA,DAS,DLAC
     &                             ,DEV,SIG,SRET,IEF,ERC)
        IF (ERC.NE.' ') GOTO 890
C
        CALL NXTDEV(DEV,NSIZE,VOL0,NW,WCR0,WCI0,IDX,DEV0)
        CALL GETDEV(DEV0,NSIZE,VOL0,NW,WCR0,WCI0
     &               ,NCR,NCI,VCR,VCI,JCR,JCI,IDX)
        CALL CLKRNL(NW,IPOL,NANG,NSIZE,NCR,NCI,VCR,VCI
     &                       ,JCR,JCI,WCR0,WCI0,AK,AKNEL)
        CALL VLRTRN(IPLC,NLN,CONCA,CONCM,NW,WL,GA,SOLID
     &       ,TAUO3,TM,WM,TH0,NA,TH,FI,ISTH,ISFI,OPT,AUR,IAUR
     &       ,NA1U0,AM1U,NFI0,FAI,NSIZE,VOL0
     &       ,IPOL,NANG,THETA,AKNEL,PPM,OPTC0,AURC0,PRTR0,ERC)
        IF (ERC.NE.' ') GOTO 890
C
        CALL CALERR(JTAU,NW,NA,OPT,AUR,OPTC0,AURC0,IAUR,ERRV)
        EROR=0.
        N=0
        DO IW=1,NW
          IF (ERRV(IW).LT.1.E10) THEN
             N=N+1
             EROR=EROR+ERRV(IW)
          ENDIF
        ENDDO
        IF (N.GT.0) THEN
           EROR=EROR/FLOAT(N)
        ELSE
           EROR=1.E10
           IFLG=-1
        ENDIF
        IF ((XJ1.LE.XJ0).OR.(EROR.LE.ERRB)) THEN
           IF (ABS(EROR-ERRB).LE.EPSA1) IFLG=2
           IF (EROR.LE.EPSA) IFLG=1
           IF (LOOP.GE.3) THEN
              IF ((EROR.GT.EPSA2)) IFLG=-1
           ENDIF
        ELSE
           IFLG=3
        ENDIF
c       print*, LOOP,XJ1,EROR,IFLG
        IF ((LOOP.GE.3).AND.(IFLG.EQ.0)) THEN
           CALL CHKERR(NW,NA,AUR,AURC0,IAUR,COEF,EROR)
           ERRB=EROR
           XJ0=XJ1
        ENDIF
      ENDDO
      IF (IFLG.GE.0) THEN
         CALL STSRET(IDX,NW,NSIZE,SIG,SRETV,SRETR,SRETI)
      ELSE
         ERC='INVLP5: EROR.GT.EPSA2'
      ENDIF
      GOTO 900
C
  890 IFLG=-9
  900 RETURN
      END
C
      SUBROUTINE SETWGT(JTAU,NW,NA,ASGM,FSGM,MA,WGTM
     &                          ,IDX,NSIZE,NINV,WGTP)
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
      PARAMETER (KNSIZE=60)
      PARAMETER (KKDT=KDT*KNW+KNW)
      PARAMETER (KKSZ=KNSIZE+KNW*2)
C
      DIMENSION ASGM(KNW,KDT),FSGM(KNW),WGTM(KKDT),WGTP(KKSZ)
C
c     DATA SGV,SGR,SGI / 2.5, 0.05, 1./
      DATA SGV,SGR,SGI / 2.5, 0.1,  1./
C
      DO I=1,NSIZE
        WGTP(I)=SGV**2
      ENDDO
      NINV=NSIZE
      IF (ABS(IDX).GT.0) THEN
         IF (IDX.GT.0) THEN
            DO I=1,NW
              WGTP(NINV+I)=SGR**2
            ENDDO
         ELSE
            DO I=1,NW
              WGTP(NINV+I)=SGI**2
            ENDDO
         ENDIF
         NINV=NINV+NW
         IF (IDX.EQ.2) THEN
            DO I=1,NW
              WGTP(NINV+I)=SGI**2
            ENDDO
            NINV=NINV+NW
         ENDIF
      ENDIF
C
      IF (JTAU.GE.0) THEN
         DO IW=1,NW
           K=NA*(IW-1)
           DO I=1,NA
             WGTM(K+I)=ASGM(IW,I)**2
           ENDDO
         ENDDO
         MA=NA*NW
      ELSE
         MA=0
      ENDIF
      IF (ABS(JTAU).GT.0) THEN
         DO IW=1,NW
           WGTM(MA+IW)=FSGM(IW)**2
         ENDDO
         MA=MA+NW
      ENDIF
      RETURN
      END
C
      SUBROUTINE CALDLT(JTAU,NW,NA,NINV,WGTP,DEV0,DEVA,DAS
     &                 ,MA,WGTM,IAUR,AUR,AURC0,OPT,OPTC0,DLA,IEF,XJ)
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
      PARAMETER (KNSIZE=60)
      PARAMETER (KKDT=KDT*KNW+KNW)
      PARAMETER (KKSZ=KNSIZE+KNW*2)
C
      DIMENSION AUR(KNW,KDT),AURC0(KNW,KDT),OPT(KNW),OPTC0(KNW)
      DIMENSION IAUR(KNW,KDT),DLA(KKDT),IEF(KKDT),WGTM(KKDT)
      DIMENSION DAS(KKSZ),DEV0(KKSZ),DEVA(KKSZ),WGTP(KKSZ)

C
      DO J=1,NINV
        DAS(J)=DEV0(J)-DEVA(J)
      ENDDO
C
      DO IW=1,NW
        IF (JTAU.GE.0) THEN
           K=NA*(IW-1)
           DO I=1,NA
             IF (IAUR(IW,I).GT.0) THEN
                DLA(K+I)=LOG(AUR(IW,I))-LOG(AURC0(IW,I))
                IEF(K+I)=1
             ELSE
                DLA(K+I)=0.
                IEF(K+I)=0
             ENDIF
           ENDDO
           K=NA*NW
        ELSE
           K=0
        ENDIF
        IF (ABS(JTAU).GT.0) THEN
           IF (OPT(IW).GT.0) THEN
              DLA(K+IW)=LOG(OPT(IW))-LOG(OPTC0(IW))
              IEF(K+IW)=1
           ELSE
              DLA(K+IW)=0.
              IEF(K+IW)=0
           ENDIF
        ENDIF
      ENDDO
      XJ=CLJFNC(MA,NINV,WGTM,WGTP,DLA,DAS)
      RETURN
      END
C
      FUNCTION CLJFNC(MA,NINV,WGTM,WGTP,DLA,DAS)
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
      PARAMETER (KNSIZE=60)
      PARAMETER (KKDT=KDT*KNW+KNW)
      PARAMETER (KKSZ=KNSIZE+KNW*2)
C
      DIMENSION DLA(KKDT),DAS(KKSZ),WGTM(KKDT),WGTP(KKSZ)
C
      XBX=0.0
      DO J=1,NINV
        XBX=XBX+0.5*DAS(J)**2/WGTP(J)
      ENDDO
      YRY=0.0
      DO I=1,MA
        YRY=YRY+0.5*DLA(I)**2/WGTM(I)
      ENDDO
      CLJFNC=XBX+YRY
      RETURN
      END
C
C
      SUBROUTINE STSRET(IDX,NW,NSIZE,SRET,SRETV,SRETR,SRETI)
      PARAMETER (KNW   =20)
      PARAMETER (KNSIZE=60)
      PARAMETER (KKSZ=KNSIZE+KNW*2)
C
      DIMENSION SRET(KKSZ),SRETV(KNSIZE),SRETR(KNW),SRETI(KNW)
C
      DO I=1,NW
        SRETR(I)=0.
        SRETI(I)=0.
      ENDDO
C
      DO I=1,NSIZE
        SRETV(I)=SRET(I)
      ENDDO
      IF (ABS(IDX).GT.0) THEN
         IF (IDX.GT.0) THEN
            DO I=1,NW
              SRETR(I)=SRET(NSIZE+I)
            ENDDO
         ELSE
            DO I=1,NW
              SRETI(I)=SRET(NSIZE+I)
            ENDDO
         ENDIF
         IF (IDX.EQ.2) THEN
            DO I=1,NW
              SRETI(I)=SRET(NSIZE+NW+I)
            ENDDO
         ENDIF
      ENDIF
      RETURN
      END
C
      SUBROUTINE SETDEV(NSIZE,VOL,NW,WCR,WCI,IDX
     &                                  ,DEV,DDV,VRATE)
      PARAMETER (KNW   =20)
      PARAMETER (KNSIZE=60)
      PARAMETER (KKSZ=KNSIZE+KNW*2)
C
      DIMENSION VOL(KNSIZE),WCR(KNW),WCI(KNW),DEV(KKSZ),DDV(KKSZ)
C
      DO J=1,NSIZE
        DEV(J)=LOG(VOL(J))
        DDV(J)=VRATE
      ENDDO
      NT=NSIZE
      IF (ABS(IDX).GT.0) THEN
         IF (IDX.GT.0) THEN
            DO IW=1,NW
              DEV(NT+IW)=LOG(WCR(IW))
              DDV(NT+IW)=VRATE
            ENDDO
            IF (IDX.GT.1) NT=NT+NW
         ENDIF
         IF (IDX.NE.1) THEN
            DO IW=1,NW
              DEV(NT+IW)=LOG(CIFUNC(WCI(IW)))
              DDV(NT+IW)=VRATE
            ENDDO
         ENDIF
      ENDIF
      RETURN
      END
C
      SUBROUTINE NXTDEV(DEV,NSIZE,VOL,NW,WCR,WCI,IDX,DEV0)
      PARAMETER (KNW   =20)
      PARAMETER (KNSIZE=60)
      PARAMETER (KKSZ=KNSIZE+KNW*2)
C
      DIMENSION VOL(KNSIZE),WCR(KNW),WCI(KNW),DEV(KKSZ),DEV0(KKSZ)
C
      DO J=1,NSIZE
        V=-VOL(J)*DEV(J)
        IF (V.GE.VOL(J)) V=VOL(J)*.5
        VOL(J)=VOL(J)-V
        DEV0(J)=LOG(VOL(J))
      ENDDO
      NT=NSIZE
      IF (ABS(IDX).GT.0) THEN
         IF (IDX.GT.0) THEN
            DO IW=1,NW
              V=-WCR(IW)*DEV(NT+IW)
              IF (V.GE.WCR(IW)) V=WCR(IW)*.5
              WCR(IW)=WCR(IW)-V
              DEV0(NT+IW)=LOG(WCR(IW))
            ENDDO
            IF (IDX.GT.1) NT=NT+NW
         ENDIF
         IF (IDX.NE.1) THEN
            DO IW=1,NW
              V=-CIFUNC(WCI(IW))*DEV(NT+IW)
              IF (V.GE.ABS(WCI(IW))) V=ABS(WCI(IW))*.5
              WCI(IW)=ABS(WCI(IW))-V
              DEV0(NT+IW)=LOG(CIFUNC(WCI(IW)))
            ENDDO
         ENDIF
      ENDIF
      RETURN
      END
C
      SUBROUTINE GETDEV(DEV,NSIZE,VOL,NW,WCR,WCI
     &                     ,NCR,NCI,VCR,VCI,JCR,JCI,IDX)
      PARAMETER (KCR   =20)
      PARAMETER (KCI   =20)
      PARAMETER (KNW   =20)
      PARAMETER (KNSIZE=60)
      PARAMETER (KKSZ=KNSIZE+KNW*2)
C
      DIMENSION VOL(KNSIZE),WCR(KNW),WCI(KNW),DEV(KKSZ)
     &         ,JCR(KNW),JCI(KNW),VCR(KCR),VCI(KCI)
C
      DO J=1,NSIZE
        VOL(J)=EXP(DEV(J))
      ENDDO
      IF (ABS(IDX).GT.0) THEN
         NT=NSIZE
         IF (IDX.GT.0) THEN
            DO IW=1,NW
              WCR(IW)=EXP(DEV(NT+IW))
              IF (WCR(IW).GT.VCR(NCR-1)) WCR(IW)=VCR(NCR-1)
              IF (WCR(IW).LT.VCR(1)) WCR(IW)=VCR(1)
              JCR(IW)=0
              N=NCR+1
              DO WHILE((JCR(IW).EQ.0).AND.(N.GT.1))
                N=N-1
                IF (VCR(N).LE.WCR(IW)) JCR(IW)=N
              ENDDO
              IF (JCR(IW).EQ.0) JCR(IW)=1
            ENDDO
            IF (IDX.GT.1) NT=NT+NW
         ENDIF
         IF (IDX.NE.1) THEN
            DO IW=1,NW
              WCI(IW)=ZIFUNC(EXP(DEV(NT+IW)))
              IF (WCI(IW).GT.ABS(VCI(NCI-1))) WCI(IW)=ABS(VCI(NCI-1))
              JCI(IW)=0
              N=NCI+1
              DO WHILE((JCI(IW).EQ.0).AND.(N.GT.1))
                N=N-1
                IF (ABS(VCI(N)).LE.WCI(IW)) JCI(IW)=N
              ENDDO
              IF (JCI(IW).EQ.0) JCI(IW)=1
              WCI(IW)=-WCI(IW)
            ENDDO
         ENDIF
      ENDIF
      RETURN
      END
C
      SUBROUTINE CHKDEV(IR,NSIZE,NW,IDX,WCR,WCR0,WCI,WCI0
     &                                 ,NCR,NCI,VCR,VCI,JCR,JCI)
      PARAMETER (KCR   =20)
      PARAMETER (KCI   =20)
      PARAMETER (KNW   =20)
C
      DIMENSION WCR(KNW),WCR0(KNW),WCI(KNW),WCI0(KNW)
     &         ,JCR(KNW),JCI(KNW),VCR(KCR),VCI(KCI)
C
      IF (IR.GT.NSIZE) THEN
         K=IR-NSIZE
         IF (IDX.GT.0) THEN
            IF (K.LE.NW) THEN
               IF (ABS(WCR(K)-WCR0(K)).LT.0.01) THEN
                  IF (JCR(K).EQ.NCR-1) THEN
                     JCR(K)=JCR(K)-1
                  ELSE
                     JCR(K)=JCR(K)+1
                  ENDIF
                  WCR(K)=VCR(JCR(K))
               ENDIF
            ELSE
               K=K-NW
               IF (ABS(WCI(K)-WCI0(K)).LT.0.0001) THEN
                  IF (JCI(K).EQ.NCI-1) THEN
                     JCI(K)=JCI(K)-1
                  ELSE
                     JCI(K)=JCI(K)+1
                  ENDIF
                  WCI(K)=-ABS(VCI(JCI(K)))
               ENDIF
            ENDIF
         ELSE
            IF (ABS(WCI(K)-WCI0(K)).LT.0.0001) THEN
               IF (JCI(K).EQ.NCI-1) THEN
                  JCI(K)=JCI(K)-1
               ELSE
                  JCI(K)=JCI(K)+1
               ENDIF
               WCI(K)=-ABS(VCI(JCI(K)))
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END
C
      FUNCTION CIFUNC(CI)
C
      CIFUNC=0.05+ABS(CI)
      RETURN
      END
C
      FUNCTION ZIFUNC(ZI)
C
      ZIFUNC=ZI-0.05
      IF (ZIFUNC.LT.0.) ZIFUNC=0.
      RETURN
      END
C
      SUBROUTINE CALERR(JTAU,NW,NA,OPT,AUR,OPTC,AURC,IAUR,EROR)
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
C
      DIMENSION AUR(KNW,KDT),OPT(KNW),AURC(KNW,KDT),OPTC(KNW)
     &         ,IAUR(KNW,KDT),EROR(KNW)
C
      DO IW=1,NW
        SUM=0.
        N=0
        IF (JTAU.GE.0) THEN
           DO I=1,NA
             IF (IAUR(IW,I).GT.0) THEN
                N=N+1
                SUM=SUM+(AURC(IW,I)/AUR(IW,I)-1.)**2
             ENDIF
           ENDDO
        ENDIF
        IF (ABS(JTAU).GT.0) THEN
           IF (OPT(IW).GT.0.) THEN
              N=N+1
              SUM=SUM+(OPTC(IW)/OPT(IW)-1.)**2
           ENDIF
        ENDIF
        IF (N.GT.0) THEN
           EROR(IW)=SQRT(SUM/FLOAT(N))
        ELSE
           EROR(IW)=1.E10
        ENDIF
      ENDDO
      RETURN
      END
C
      SUBROUTINE CHKERR(NW,NA,AUR,AURC,IAUR,COEF,ERR)
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
C
      DIMENSION AUR(KNW,KDT),AURC(KNW,KDT),IAUR(KNW,KDT)
C
      DO IW=1,NW
        DO I=1,NA
          IF (ABS(IAUR(IW,I)).EQ.1) THEN
             A0=ABS(AURC(IW,I)/AUR(IW,I)-1.)
             IF (A0.LT.ERR*COEF) THEN
                IAUR(IW,I)=1
             ELSE
                IAUR(IW,I)=-1
c               print*,IW,I,A0,ERR,ERR*COEF
             ENDIF
          ENDIF
        ENDDO
      ENDDO
      RETURN
      END
C
      SUBROUTINE CKRNL8(IUK,INVM,SSL,NCR,NCI,VCR,VCI,NWL
     &                 ,WL,RMIN,RMAX,NSIZE,NANG,THETA,SIZE,AK,ERC)
C
C GET KERNEL MATRIX FOR MIE CALCULATION.
C with new type kernel files (IFORM=2) and matrix reading
C for rectangular CR and CI
C Caution : Total scattering cross section output difference from CKRNL5
C
C --- notes
C FOR THE KERNEL-TABLE PRIOR TO SEPT,1987: CEXT and CSCA
C After that                             : CEXT and CABS
C IPOL=1 ONLY
C
C --- history
C 85. 4.24
C     6.15  MISTAKE IN SZIMUL   ->CALCULATION OF SIZE
C 86. 6. 5  ADD 179 IF , 201  FOR LARGE DELTA
C 86.11.28 CHANGE IF STATEMENT (IF 174,176,177->174,174,177)
C            AND ADD CORRESPONDING IF FOR TOTAL CROSS SECTIONS.
C 88. 4. 5 FOR IBM. FOR CEXT, CABS
C 93. 6.19 Introduce ERC
C 98. 6. 6 Matrix reading, add molti-mode weighting
C          Absorption to scattering output
C 00.11.29 SIZE(NSIZE) is changed to be set outside of this subroutine
C
C --- Input
C IUK     I         READ DEVICE NO. OF KERNEL FILE FROM CONKRNL
C INVM    I         =<19: box kernel
C                     SIZE are given by this program
C                   >=20: log-normal multi mode kernel
C                   give mode sized by SIZE and mode dispersion by SSL
C NWL     I         NO. OF WAVELENGTHS
C WL      R(NWL)    WAVELENGTHS IN CM FOR TOTAL CROSS SECTIONS.
C RMIN    R         MINIMUM PARTICLE RADII IN CM
C RMAX    R         MINIMUM PARTICLE RADII IN CM
C NSIZE   I         NO. OF SIZE INTERVALS
C SIZE    R(NSIZE)  PARTICLE MODE RADIUS IN CM AT LOG.CENTER
C                   OF EACH SIZE INTERVAL
C SSL     R(NSIZE)  log-dispersion of log-normal multi mode kernel
C                     (only for INVM>19)
C
C--- Output
C NCR,NCI I         number of gird points for real/imaginary part data
C VCR     R(KCR)    real part data array of kernel matrix
C VCI     R(KCI)    imaginary part data array of kernel matrix
C NANG    I         NO. OF SCATTERING ANGLES
C THETA   R(KNANG)  SCATTERING ANGLES IN DEGREE
C AK      R(KNMEAS,KNSIZE,KCR,KCI)
C                   KERNEL MATRIX:  (I,J) ,I=1,NMEAS, J=1,NSIZE
CC  DEFINITION OF KERNEL GIVEN BY THIS ROUTINE
CC   OBSERVATION VECTOR  P=
CC   (((PHASE F. FOR S=S(1) TO S(NANG)) FOR IP=1 TO IPOL)
CC     FOR WL=WL(1) TO WL(NWL))
CC   (E.CROSS S. FOR WL=WL(1) TO WL(NWL))
CC   (S.CROSS S. FOR WL=WL(1) TO WL(NWL))
CC   P(I)= SUM( AK(I,J)*V(J)) FOR J=1 TO NSIZE
CC   V=  DV/DLNR
C
C--- CHECK PARAMETERS
C KINTVL     NO. OF SUBINTERVALS IN THE KERNEL TABLE.
C KSIZE1     NSIZE+1
C KBUF       MAX(NANG, INTVL)
C$ENDI
      PARAMETER (PAI=3.141592653589793)
C 20240529 Modified : from Momoi-san, 2018
C      PARAMETER (KNMEAS=1240)
      PARAMETER (KINTVL=60)
      PARAMETER (KNSIZE=60)
      PARAMETER (KNANG =60)
      PARAMETER (KCR   =20)
      PARAMETER (KCI   =20)
      PARAMETER (KNW   =20)
C 20240529 Modified : from Momoi-san, 2018
      PARAMETER (KNMEAS=KNANG*KNW+2*KNW)
      CHARACTER ERC*(*)
      DIMENSION AK(KNMEAS,KNSIZE,KCR,KCI),THETA(KNANG),SIZE(KNSIZE)
     &,WL(KNW),VCR(KCR),VCI(KCI),SSL(KNSIZE)
C WORKING AREA
      PARAMETER (KSIZE1=KNSIZE+1,KFN=10,KINT=KFN*KINTVL)
      CHARACTER FNAME2*1
      DIMENSION BUF(KNANG,KINTVL),SLBND(KSIZE1),SZPARA(KINTVL)
     &  ,SZD(KNW,KNSIZE,KINTVL)
C
      ERC=' '
      REWIND IUK
      READ(IUK,*) IFORM
      IF(IFORM.NE.2) then
        ERC='CKRNL8: Kernel file format not suitable'
        return
      endif
      READ(IUK,*) NCR, NCI
      READ(IUK,*) (VCR(I),I=1,NCR)
      READ(IUK,*) (VCI(I),I=1,NCI)
      READ(IUK,*) INTVL,NANG,IPOL
      IF(IPOL.NE.1) THEN
        ERC='CKRNL8: only for IPOL=1'
        return
      endif
      IF(INTVL.GT.KINTVL) GOTO 3
      IF(NSIZE.GT.KNSIZE) goto 3
      IF(NANG.GT.KNANG) GOTO 3
      NANG2=IPOL*NANG
      NMEAS2=NANG2*NWL
      NMEAS1=NMEAS2+NWL
      NMEAS=NMEAS1+NWL
      IF(NMEAS.GT.KNMEAS) GOTO 3
      READ(IUK,*) (SZPARA(I),I=1,INTVL)
      DELTA=LOG(SZPARA(2)/SZPARA(1))
      IF(NANG.GT.0) READ(IUK,*) (THETA(I),I=1,NANG)
      WLMIN=1.E10
      WLMAX=-1.0
      IF(NWL.GT.0) then
        WLMIN=AMIN1(WLMIN,WL(1))
        WLMAX=AMAX1(WLMAX,WL(NWL))
      endif
      RMAX1=WLMIN*SZPARA(INTVL)/2./PAI
      RMAX=AMIN1(RMAX,RMAX1)
      RMIN1=WLMAX*SZPARA(1)/2./PAI
      RMIN=AMAX1(RMIN,RMIN1)
C
      IF(INVM.LE.19) then
        DLSIZE=LOG(RMAX/RMIN)/NSIZE
        SZL=LOG(RMIN)-DLSIZE
        DO 171 I=1,NSIZE+1
        SZL=SZL+DLSIZE
  171   SLBND(I)=SZL
C --- SIZE(NSIZE) is changed to be set outside of this subroutine(00/1/29)
c       DO 172 I=1,NSIZE
c 172   SIZE(I)=EXP((SLBND(I)+SLBND(I+1))/2)
      ENDIF
C --------------------------------------------------------------------------
      DO 12 IW=1,NWL
      WNUM=2*PAI/WL(IW)
      DO 12 I=1,INTVL
      R1=SZPARA(I)/WNUM
      RL1=LOG(R1)-0.5*DELTA
      RL2=RL1+DELTA
      IF(INVM.LE.19) THEN
        DO 175 IS=1,NSIZE
        SZD(IW,IS,I)=0
        IF(RL1-SLBND(IS)) 174,174,177
  174   IF(RL2-SLBND(IS)) 175,175,179
  179   IF(RL2-SLBND(IS+1)) 200,200,201
  201   DL1=(SLBND(IS+1)-SLBND(IS))/DELTA
        GOTO 180
  200   DL1=(RL2-SLBND(IS))/DELTA
        GO TO 180
  177   IF(RL1-SLBND(IS+1)) 181,175,175
  181   IF(RL2-SLBND(IS+1)) 176,176,182
  182   DL1=(SLBND(IS+1)-RL1)/DELTA
        GO TO 180
  176   DL1=1
  180   SZD(IW,IS,I)=DL1*WNUM
  175   CONTINUE
       else
        DRL=(RL2-RL1)/KFN
        DL1=DRL/DELTA
        DO 15 IS=1,NSIZE
        SUM=0
        RLL=RL1-DRL/2
        do 13 L=1,KFN
        RLL=RLL+DRL
        R3=EXP(RLL)
        IF(R3.LT.RMIN .or. R3.GT.RMAX) goto 13
        EXP1=0.5*(LOG(R3/SIZE(IS))/SSL(IS))**2
        IF(EXP1.LT.100.0) SUM=SUM+EXP(-EXP1)
   13   continue
   15   SZD(IW,IS,I)=SUM*WNUM*DL1
      endif
   12 continue
C ------------------------------------------------------------
      DO 60 ICI=1,NCI
      DO 60 ICR=1,NCR
      READ(IUK,27,END=29) FNAME2
   27 FORMAT(A1)
      READ(IUK,*) CR1,CI1
C read differential cross section
      IF(NANG.GT.0) then
        DO 14 I=1,INTVL
   14   READ(IUK,*) (BUF(K,I),K=1,NANG)
      endif
C
      DO 4 IW=1,NWL
      DO 4 IS=1,NSIZE
      DO 4 K=1,NANG
      M=(IW-1)*NANG2+K
      SUM=0
      DO 16 I=1,INTVL
   16 SUM=SUM+SZD(IW,IS,I)*BUF(K,I)
    4 AK(M,IS,ICR,ICI)=SUM
C read extinction and absorption cross section
      do 5 K=1,2
    5 read(IUK,*) (BUF(K,I),I=1,INTVL)
      do 17 K=1,2
      DO 17 IW=1,NWL
      M=NMEAS2+(K-1)*NWL+IW
      DO 17 IS=1,NSIZE
      SUM=0
      DO 18 I=1,INTVL
   18 SUM=SUM+SZD(IW,IS,I)*BUF(K,I)
   17 AK(M,IS,ICR,ICI)=SUM
C Absorption cross section to scattering cross section !!!
      DO  50 I=1,NWL
      DO  50 J=1,NSIZE
   50 AK(NMEAS1+I,J,ICR,ICI)=AK(NMEAS2+I,J,ICR,ICI)
     &                      -AK(NMEAS1+I,J,ICR,ICI)
C
   60 continue
      RETURN
   29 ERC='CKRNL8: ILLEGAL COMPLEX REFRACTIVE INDEX.'
      RETURN
    3 ERC='CKRNL8: ILLEGAL PARAMETER SIZE.'
      RETURN
      END
      SUBROUTINE CLKRNL(NW,IPOL,NANG,NSIZE,NCR,NCI,VCR,VCI
     &                                    ,JCR,JCI,WCR,WCI,AK,AKC)
C
C  Calculation of kernel for current refractive index.
C
C --- history
C   2000.11.17 Created by M.Yamano
C
C --- Input
C NW     I         number of wavelength
C IPOL   I         fixed(=1) for CKRNL8
C NANG   I         number of scattering angles on the kernel table
C NSIZE  I         number of particle radii
C NCR    I         real part number of kernel matrix
C NCI    I         imaginary part number of kernel matrix
C VCR    R(KCR)    real part data array of kernel matrix
C VCI    R(KCI)    imaginary part data array of kernel matrix
C JCR    I(KNW)    real part index of current refractive index
C JCI    I(KNW)    imaginary part index od current refractive index
C WCR    R(KNW)    real part value of current refractive index
C WCI    R(KNW)    imaginary part value of current refractive index
C AK     R(KNMEAS,KNSIZE,KCR,KCI)
C                  KERNEL MATRIX:  (I,J) ,I=1,NMEAS, J=1,NSIZE
C
C --- Output
C AKC    R(KNMEAS,KNSIZE)
C                  KERNEL MATRIX for current index
C ERC    C*64      ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
C
      PARAMETER (KNW   =20)
      PARAMETER (KNANG =60)
      PARAMETER (KNSIZE=60)
      PARAMETER (KNMEAS=KNANG*KNW+2*KNW)
      PARAMETER (KCR   =20)
      PARAMETER (KCI   =20)
C
      DIMENSION AK(KNMEAS,KNSIZE,KCR,KCI),AKC(KNMEAS,KNSIZE)
      DIMENSION JCR(KNW),JCI(KNW),WCR(KNW),WCI(KNW)
     &         ,VCR(KCR),VCI(KCI)
C
      NWLE=NW
      NWLP=NW
      NANG2=NANG*IPOL
      NMEAS2=NANG2*NWLP
      NMEAS1=NMEAS2+NWLE
      NMEAS=NMEAS1+NWLE
C
      DO IW=1,NW
        ICR0=JCR(IW)
        IF (ICR0.LT.NCR) THEN
           ICR1=ICR0+1
           DCR=LOG(WCR(IW)/VCR(ICR0))/LOG(VCR(ICR1)/VCR(ICR0))
        ELSE
           ICR1=ICR0
           DCR=0
        ENDIF
        ICI0=JCI(IW)
        IF (ICI0.LT.NCI) THEN
           ICI1=ICI0+1
           DCI=LOG(CIFUNC(WCI(IW))/CIFUNC(VCI(ICI0)))
     &        /LOG(CIFUNC(VCI(ICI1))/CIFUNC(VCI(ICI0)))
        ELSE
           ICI1=ICI0
           DCI=0
        ENDIF
        II=NANG2*(IW-1)
        DO I=1,NANG2
          DO J=1,NSIZE
            A00=LOG(AK(I+II,J,ICR0,ICI0))
            A01=LOG(AK(I+II,J,ICR0,ICI1))
            A10=LOG(AK(I+II,J,ICR1,ICI0))
            A11=LOG(AK(I+II,J,ICR1,ICI1))
            AM0=A00+(A10-A00)*DCR
            AM1=A01+(A11-A01)*DCR
            AKC(I+II,J)=EXP(AM0+(AM1-AM0)*DCI)
          ENDDO
        ENDDO
        DO J=1,NSIZE
          A00=LOG(AK(NMEAS2+IW,J,ICR0,ICI0))
          A01=LOG(AK(NMEAS2+IW,J,ICR0,ICI1))
          A10=LOG(AK(NMEAS2+IW,J,ICR1,ICI0))
          A11=LOG(AK(NMEAS2+IW,J,ICR1,ICI1))
          AM0=A00+(A10-A00)*DCR
          AM1=A01+(A11-A01)*DCR
          AKC(NMEAS2+IW,J)=EXP(AM0+(AM1-AM0)*DCI)
        ENDDO
        DO J=1,NSIZE
          A00=LOG(AK(NMEAS1+IW,J,ICR0,ICI0))
          A01=LOG(AK(NMEAS1+IW,J,ICR0,ICI1))
          A10=LOG(AK(NMEAS1+IW,J,ICR1,ICI0))
          A11=LOG(AK(NMEAS1+IW,J,ICR1,ICI1))
          AM0=A00+(A10-A00)*DCR
          AM1=A01+(A11-A01)*DCR
          AKC(NMEAS1+IW,J)=EXP(AM0+(AM1-AM0)*DCI)
        ENDDO
      ENDDO
      RETURN
      END
C
      SUBROUTINE INVOLM(JTAU,IPLC,NLN,CONCA,CONCM,NW,WL,GA,SOLID
     &                 ,TAUO3,TM,WM,TH0,NA,TH,FI,ISTH,ISFI
     &                 ,OPT,AUR,IAUR,NA1U0,AM1U,NFI0,FAI,NSIZE,SIZE
     &                 ,RMIN,RMAX,IPOL,NANG,THETA,AKNEL,PPM
     &                 ,VOLI,EPSA,EPSA1,IFLG,ERC)
C--- Note
C   Set initial volume spectrum for skyradiance analysis.
C
C--- History
C  2000.11.17 Created by M.Yamano
C  2002.05.22 VOLSPC -> VLSPCN
C  2003.12.10 Added OPTC(IW) into sub VLRTRN arguments
C  2004.01.23 Added arguments JTAU,OPT(IW)/ Deleted FDW(IW)
C  2004.01.24 Added IAUR(IW,I) instead of SCA,ANGTOP,ANGBTM
C  2008.10.17 Changed initial volume spectrum (Junge -> 2-modal log-normal).
CC 2010.02.15 Modified by M.Hashimoto (VLSPCN -> VLSPC2; S -> ln(S))
C
C--- Input
C JTAU     I           controlling optical thickness constraint
C                      1: OPT is observed and AUR+OPT for inversion.
C                      0: only AUR for inversion (no use of OPT)
C                     -1: only OPT for inversion
C IPLC     I            polarization correction control
C                         1: WITH POLARIZATION CORRECTION
C                         0: WITHOUT CORRECTION
C NLN      I            Number of Layer
C CONCA    R(KNLN)      Fraction of Aerosol   Optical Thickness
C CONCM    R(KNLN)      Fraction of Molecular Optical Thickness
C NW       I            NBR OF WAVELENGTHS.
C WL       R(KNW)       WAVELENGTH IN CM.
C GA       R(KNW)       ASSUMED GROUND ALBEDO FOR EACH CHANNEL.
C SOLID    R(KNW)       Solid view angles
C TAUO3    R(KNW)       Ozone optical thickness
C TM       R(KNW)       Rayleigh optical thickness
C WM       R(KNW)       Rayleigh single scattering albedo.
C TH0      R            SOLAR ZENITH ANGLE IN DEGREE (0,90)
C NA       I            NBR OF ANGLES for AUR
C TH       R(KNW,KDT)   EMERGENT ZENITH ANGLES IN DEGREE (0,90)
C FI       R(KNW,KDT)   AZIMUTHAL ANGLES IN DEGREE  (0,180)
C ISTH     I(KNW,KDT)   Ordering No. of zenith angles
C ISFI     I(KNW,KDT)   Ordering No. of azimuthal angles
C OPT      R(KNW)       Observed aerosol optical thickness
C AUR      R(KNW,KDT)   I/EM/F/SOLID VIEW ANGLE.
C IAUR     I(KNW,KDT)   availability flag for AUR(>0 available)
C NA1U0    I(KNW)       NUMBER OF EMERGENT NADIR ANGLES IN SPHERE.
C AM1U     R(KNW,KNA1U) CONSINE OF THE EMERGENT NADIR ANGLES.
C                          - FOR UPWARD, + FOR DOWNWARD.
C NFI0     I(KNW)       NUMBER OF AZIMUTHAL ANGLES.
C FAI      R(KNW,KNFI)  AZIMUTHAL ANGLES IN DEGREES.
C NSIZE    I            NBR OF PARTICLE RADII.
C SIZE     R(KNSIZE)    mode radius of size distribution
C RMIN     R            MINIMUM PARTICLE RADIUS IN CM
C RMAX     R            MAXIMUM PARTICLE RADIUS IN CM
C IPOL     I            fixed(=1) for CKRNL8
C NANG     I            NBR OF SCATTERING ANGLE ON THE KERNEL TABLE.
C THETA    R(KNANG)     SCATTERING ANGLES(IN DEGREE) ON THE KERNEL TABLE.
C AKNEL    R(KNMEAS,KNSIZE)
C                       KERNEL MATRIX:  (I,J) ,I=1,NMEAS, J=1,NSIZE
C PPM      R(KNANG)     Phase function of Rayleigh scattering
C EPSA     R            ABSOLUTE CONVERGENCE CRITERION.
C EPSA1    R            RELATIVE CONVERGENCE CRITERION.
C
C--- Output
C VOLI     R(KNSIZE)    Optimized volume spectrum(dV/dlnR)
C IFLG     I            Return flag(>0 : available)
C ERC      C*64         ERROR CODE. IF ' ' THEN NORMAL.
C
C--- REPLACE PARAMETERS
C
C KNW      SIZE FOR NO. OF WAVELENGTHS.
C KDT      SIZE FOR NO. OF DATA.
C KNANG    SIZE FOR NO. OF SCATTERING ANGLES.
C KNSIZE   SIZE FOR NO. OF PARTICLE RADII.
C KNMEAS   SIZE FOR LENGTH OF OBSERVED DATA VECTOR,
C           I.E.  = IPOL*KNANG*KNW+2*KNW
C KNA1     SIZE FOR EMERGENT ZENITH ANGLES.
C KNLN     SIZE FOR NO. OF SUBLAYERS.
C KNFI     SIZE FOR NO. OF AZIMUTHAL ANGLES.
C------------------------------------------- AREAS FOR THIS ROUTINE.
      SAVE
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
      PARAMETER (KNANG =60)
      PARAMETER (KNSIZE=60)
      PARAMETER (KNA1  =50)
      PARAMETER (KNFI  =50)
      PARAMETER (KNLN  =4)
      PARAMETER (KNMEAS=KNANG*KNW+2*KNW, KNA1U=KNA1)
C
      PARAMETER (PI=3.141592653589793,RAD=PI/180.0)
      CHARACTER ERC*(*)
C-----FOR VLRTRN
      DIMENSION CONCA(KNLN),CONCM(KNLN)
      DIMENSION WL(KNW),GA(KNW),SOLID(KNW),TAUO3(KNW),TM(KNW),WM(KNW)
     &         ,OPT(KNW),AUR(KNW,KDT),IAUR(KNW,KDT)
      DIMENSION TH(KNW,KDT),FI(KNW,KDT),ISTH(KNW,KDT),ISFI(KNW,KDT)
      DIMENSION AM1U(KNW,KNA1U),FAI(KNW,KNFI),NA1U0(KNW),NFI0(KNW)
      DIMENSION OPTC(KNW),AURC(KNW,KDT),PRTR(KNMEAS)
      DIMENSION THETA(KNANG),AKNEL(KNMEAS,KNSIZE),PPM(KNANG)
      DIMENSION SIZE(KNSIZE),VOLI(KNSIZE),VOLS(KNSIZE)
C 10.02.15 PR(6,4) added by M.Hashimoto
      DIMENSION PR(6,4)
CC      COMMON /VSPC/RMN,RMX,PR(5,5),NMODE
      DATA INIT/0/
C
      ERC=' '
      IF ((NSIZE.GT.KNSIZE).OR.(NA.GT.KDT)) THEN
         ERC='INVOLM: Dimension size overflow'
         GOTO 900
      ENDIF
C
C Initial volume spectrum (Junge size distribution)
C
      IF (INIT.LE.0) THEN
c        NMODE=1
c        RMN=RMIN
c        RMX=RMAX
c        PR(1,1)=1
c        PR(2,1)=1.E-12
c        PR(3,1)=0.1E-4
c        PR(4,1)=4
CC 10.02.15 M.Hashimoto added 'CC' and code between CC== and CC==
CC         NMODE=2
CC         RMN=RMIN
CC         RMX=RMAX
CC         PR(1,1)=2
CC         PR(2,1)=1.E-12
CC         PR(3,1)=0.4
CC         PR(4,1)=0.1E-4
CC         PR(1,2)=2
CC         PR(2,2)=1.E-12
CC         PR(3,2)=0.8
CC         PR(4,2)=2.E-4
CC===
CC 10.02.15 M.Hashimoto
        NMODE  =2
        PR(1,2)=NMODE
        PR(1,3)=RMIN
        PR(1,4)=RMAX
        PR(2,1)=2
        PR(3,1)=1.E-12
        PR(4,1)=1.492
        PR(5,1)=0.1E-4
        PR(2,2)=2
        PR(3,2)=1.E-12
        PR(4,2)=2.225
        PR(5,2)=2.E-4
CC 10.02.15 M.Hashimoto SIZE=>PR(1,1), VLSPCN(R)->VLSPC2(PR)
         DO J=1,NSIZE
           PR(1,1)=SIZE(J)
           VOLS(J)=VLSPC2(PR)
CC           R=SIZE(J)
CC           VOLS(J)=VLSPCN(R)
         ENDDO
         INIT=1
      ENDIF
C
      DO J=1,NSIZE
        VOLI(J)=VOLS(J)
      ENDDO
      ERR=1.E10
      C=1.
      IFLG=0
      LOOPMX=15
      LOOP=-1
      DO WHILE((IFLG.EQ.0).AND.(LOOP.LT.LOOPMX))
        LOOP=LOOP+1
        CVL=1./C
        DO J=1,NSIZE
          VOLI(J)=CVL*VOLI(J)
        ENDDO
        CALL VLRTRN(IPLC,NLN,CONCA,CONCM,NW,WL,GA,SOLID
     &             ,TAUO3,TM,WM,TH0,NA,TH,FI,ISTH,ISFI,OPT,AUR,IAUR
     &             ,NA1U0,AM1U,NFI0,FAI,NSIZE,VOLI
     &             ,IPOL,NANG,THETA,AKNEL,PPM,OPTC,AURC,PRTR,ERC)
        IF (ERC.NE.' ') GOTO 900
C
        C=0.
        S=0.
        N=0
        DO IW=1,NW
          IF (JTAU.GE.0) THEN
             DO I=1,NA
               IF (IAUR(IW,I).GT.0) THEN
                  N=N+1
                  C=C+AURC(IW,I)/AUR(IW,I)
                  S=S+(AURC(IW,I)/AUR(IW,I)-1.)**2
               ENDIF
             ENDDO
           ENDIF
           IF (ABS(JTAU).GT.0) THEN
              IF (OPT(IW).GT.0.) THEN
                 N=N+1
                 C=C+OPTC(IW)/OPT(IW)
                 S=S+(OPTC(IW)/OPT(IW)-1.)**2
              ENDIF
           ENDIF
        ENDDO
        ERRB=ERR
        IF (N.GT.0) THEN
           C=C/FLOAT(N)
           ERR=SQRT(S/FLOAT(N))
        ELSE
           ERC='INVOLM: No available data.'
           IFLG=-1
        ENDIF
c       IF (ABS(ERR-ERRB).LE.EPSA1) IFLG=2
        IF (ABS(C-1.).LE.0.0001) IFLG=2
        IF (ERR.LE.EPSA) IFLG=1
c       print*, LOOP,C,CVL,ERR,IFLG
      ENDDO
c     IF (IFLG.EQ.0) THEN
c          ERC='INVOLM: Not converged data.'
c     ENDIF
C
  900 RETURN
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
      SUBROUTINE VLRTRN(IPLC,NLN,CONCA,CONCM,NW,WL,GA,SOLID
     &                 ,TAUO3,TM,WM,TH0,NA,TH,FI,ISTH,ISFI,OPT,AUR,IAUR
     &                 ,NA1U0,AM1U,NFI0,FAI,NSIZE,VOL
     &                 ,IPOL,NANG,THETA,AKNEL,PPM,OPTC,AURC,PRTR,ERC)
C--- Note
C   Calculation of skyradiance(L/EM/F)
C   for known refractive index and volume spectrum.
C
C--- History
C  2001.01.17 Created by M.Yamano
C  2002.07.03 Debugged.(XX is set)
C  2003.05.09 Debugged.(TAUMX check is added)
C  2003.11.26 Changed.(FDWC(IW)=FDW(IW)*FDW1 -> FDWC(IW)=FD)
C  2003.12.10 Added output argument OPTC(IW)
C  2004.01.21 Changed FDW(IW) -> OPT(IW)/ Deleted FDWC(IW)
C  2004.01.24 Added IAUR(IW,I) instead of AUR,SCA,ANGTOP,ANGBTM
C
C--- Input
C IPLC     I           polarization correction control
C                        1: WITH POLARIZATION CORRECTION
C                        0: WITHOUT CORRECTION
C NLN      I           Number of Layer
C CONCA    R(KNLN)     Fraction of Aerosol   Optical Thickness
C CONCM    R(KNLN)     Fraction of Molecular Optical Thickness
C NW       I           NBR OF WAVELENGTHS.
C WL       R(KNW)      WAVELENGTH IN CM.
C GA       R(KNW)      ASSUMED GROUND ALBEDO FOR EACH CHANNEL.
C SOLID    R(KNW)      Solid view angles
C TAUO3    R(KNW)      Ozone optical thickness
C TM       R(KNW)      Rayleigh optical thickness
C WM       R(KNW)      Rayleigh single scattering albedo.
C TH0      R           SOLAR ZENITH ANGLE IN DEGREE (0,90)
C NA       I           NBR OF ANGLES for AUR
C TH       R(KNW,KDT)  EMERGENT ZENITH ANGLES IN DEGREE (0,90)
C FI       R(KNW,KDT)  AZIMUTHAL ANGLES IN DEGREE  (0,180)
C ISTH     I(KNW,KDT)  Ordering No. of zenith angles
C ISFI     I(KNW,KDT)  Ordering No. of azimuthal angles
C OPT      R(KNW)      Observed aerosol optical thickness
C AUR      R(KNW,KDT)   I/EM/F/SOLID VIEW ANGLE. CC: 20220404 modified by hashimoto 
C IAUR     I(KNW,KDT)  availability flag for AUR(ABS(IAUR)=1 then available)
C NA1U0    I(KNW)      NUMBER OF EMERGENT NADIR ANGLES IN SPHERE.
C AM1U     R(KNA1U)    CONSINE OF THE EMERGENT NADIR ANGLES.
C                         - FOR UPWARD, + FOR DOWNWARD.
C NFI0     I(KNW)      NUMBER OF AZIMUTHAL ANGLES.
C FAI      R(KNFI)     AZIMUTHAL ANGLES IN DEGREES.
C NSIZE    I           NBR OF PARTICLE RADII.
C VOL      R(KNSIZE)   DV/DLNR
C IPOL     I           fixed(=1) for CKRNL8
C NANG     I           NBR OF SCATTERING ANGLE ON THE KERNEL TABLE.
C THETA    R(KNANG)    SCATTERING ANGLES(IN DEGREE) ON THE KERNEL TABLE.
C AKNEL    R(KNMEAS,KNSIZE)
C                      KERNEL MATRIX:  (I,J) ,I=1,NMEAS, J=1,NSIZE
C PPM      R(KNANG)    Phase function of Rayleigh scattering
C
C--- Output
C OPTC     R(KNW)      Reconstructed optical thickness
C AURC     R(KNW,KDT)  RECONSTRUCTED AUREOLE.
C PRTR     R(KNMEAS)   RETRIEVED MIE PHASE FUNCTION AFTER REDUCING MULTI-SCAT.
C ERC      C*64        ERROR CODE. IF ' ' THEN NORMAL.
C
C--- REPLACE PARAMETERS
C
C KLGN1    MAXIMUM ORDER FOR LEGENDRE EXPANSION OF PHASE FUNCTION.
C KNW      SIZE FOR NO. OF WAVELENGTHS.
C KDT      SIZE FOR NO. OF DATA.
C KNANG    SIZE FOR NO. OF SCATTERING ANGLES.
C KNSIZE   SIZE FOR NO. OF PARTICLE RADII.
C KNMEAS   SIZE FOR LENGTH OF OBSERVED DATA VECTOR,
C           I.E.  = IPOL*KNANG*KNW+2*KNW
C KNA1     SIZE FOR EMERGENT ZENITH ANGLES.
C KNA0     SIZE FOR INCIDENT ZENITH ANGLES.
C KNDM     SIZE FOR NO. OF DISCRETE PATHS FOR TRANSFER.
C KNLN     SIZE FOR NO. OF SUBLAYERS.
C KNTAU    SIZE FOR NO. OF USER-DEFINED OPTICAL DEPTHS.
C KNFI     SIZE FOR NO. OF AZIMUTHAL ANGLES.
C------------------------------------------- AREAS FOR THIS ROUTINE.
      PARAMETER (KNW   =20)
      PARAMETER (KDT   =50)
      PARAMETER (KNANG =60)
      PARAMETER (KNSIZE=60)
      PARAMETER (KNA1  =50)
      PARAMETER (KNFI  =50)
      PARAMETER (KNLN  =4)
      PARAMETER (KNMEAS=KNANG*KNW+2*KNW,KNA1U=KNA1)
      PARAMETER (KLGN1 =400)
      PARAMETER (KNA0  =1)
      PARAMETER (KNTAU =1)
      PARAMETER (KNDM  =16)
      PARAMETER (KPLK1 =2)
C
      PARAMETER (PI=3.141592653589793,RAD=PI/180.0)
      CHARACTER ERC*(*)
      DIMENSION CONCA(KNLN),CONCM(KNLN)
      DIMENSION WL(KNW),GA(KNW),SOLID(KNW),TAUO3(KNW),TM(KNW),WM(KNW)
     &         ,OPT(KNW),AUR(KNW,KDT),IAUR(KNW,KDT)
C 20220404 hashimoto  &  ,OPT(KNW),IAUR(KNW,KDT)
      DIMENSION TH(KNW,KDT),FI(KNW,KDT),ISTH(KNW,KDT),ISFI(KNW,KDT)
      DIMENSION AM1U(KNW,KNA1U),FAI(KNW,KNFI),NA1U0(KNW),NFI0(KNW)
      DIMENSION OPTC(KNW),AURC(KNW,KDT)
      DIMENSION AKNEL(KNMEAS,KNSIZE),VOL(KNSIZE),PRTR(KNMEAS)
C-----FOR RADIATIVE TRANSFER.
      DIMENSION PPA(KNANG),PPM(KNANG)
C-----FOR RTRN1
      DIMENSION AM1UR(KNA1U),FAIR(KNFI),THK(KNLN),OMG(KNLN)
     &         ,NLGN1(KNLN),THETA(KNANG),PHSF(KNANG,KNLN)
     &         ,AI(KNA1U,KNA0,KNFI,KNTAU)
C-----OTHER WORKING AREAS.
      DIMENSION COR(2),ALF(4)
      DOUBLE PRECISION XX(KNANG)
C
      DATA TAUMX /100./
C
      ERC=' '
      IF ((NSIZE.GT.KNSIZE).OR.(NA.GT.KDT)) THEN
         ERC='VLRTRN: Dimension size overflow'
         GOTO 900
      ENDIF
C
      NWLE=NW
      NWLP=NW
      NANG2=NANG*IPOL
      NMEAS2=NANG2*NWLP
      NMEAS1=NMEAS2+NWLE
      NMEAS=NMEAS1+NWLE
      DO I=1,NMEAS
        PRTR(I)=0.
        DO J=1,NSIZE
          PRTR(I)=PRTR(I)+AKNEL(I,J)*VOL(J)
        ENDDO
      ENDDO
      IW=0
      DO WHILE((IW.LT.NW).AND.(ERC.EQ.' '))
        IW=IW+1
        IF (PRTR(NMEAS2+IW).GE.TAUMX)
     &     ERC='VLRTRN: Give up - Too large TAU'
      ENDDO
      IF (ERC.NE.' ') GOTO 900
C
      DO J=1,NANG
        XX(J)=COS(THETA(J)*RAD)
      ENDDO
C
      EM=1.0/COS(TH0*RAD)
      IMTHD=3
      DO IW=1,NW
C 98.0528 for diffuse correction
        NA1U=NA1U0(IW)
        DO I=1,NA1U
          AM1UR(I)=AM1U(IW,I)
        ENDDO
        NFI=NFI0(IW)
        DO I=1,NFI
          FAIR(I)=FAI(IW,I)
        ENDDO
C
        II=NANG2*(IW-1)
        TA1=PRTR(NMEAS2+IW)
        WA1=PRTR(NMEAS1+IW)/PRTR(NMEAS2+IW)
        DO I=1,NANG
          PPA(I)=PRTR(I+II)/PRTR(NMEAS1+IW)
        ENDDO
        TM1=TM(IW)
        WM1=WM(IW)
        CALL CLUTAU(NLN,CONCA,CONCM,TA1,WA1,TM1,WM1,NANG,PPA,PPM
     &             ,THK,OMG,NLGN1,PHSF,TAU)
        SOLID1=SOLID(IW)
        TRO3=EXP(-EM*TAUO3(IW))
        GALB=GA(IW)
        CALL PRTRN1(IMTHD,NLN,THK,OMG,NLGN1,TAU,NANG,THETA,PHSF,TH0
     &             ,SOLID1,TRO3,GALB,NA1U,AM1UR,NFI,FAIR,AI,FD,FDW1,ERC)
        IF (ERC.NE.' ') GOTO 900
C
        IF (OPT(IW).GT.0.) THEN
           OPTC(IW)=TAU-TM(IW)
        ELSE
           OPTC(IW)=-1.
        ENDIF
C
        LMAX9=3
        CALL LGNDF1(LMAX9,NANG,XX,PPA,ALF)
        P112=ALF(3)/ALF(1)
        P122=0
        DO I=1,NA
          IF (ABS(IAUR(IW,I)).EQ.1) THEN
             COR(1)=0
             IF (IPLC.GT.0) THEN
                THW=TH(IW,I)
                FIW=FI(IW,I)
                CALL PCOR(COR,TH0,THW,FIW
     &                   ,TM1,WM1,DPF,TM1,TA1,WA1,P112,P122)
             ENDIF
             ISTHW=ISTH(IW,I)
             ISFIW=ISFI(IW,I)
             AURC(IW,I)=(AI(ISTHW,1,ISFIW,1)+COR(1))*TRO3/EM/FD
          ELSE
             AURC(IW,I)=-1.
          ENDIF
        ENDDO
      ENDDO
C
  900 RETURN
      END
C
      SUBROUTINE INVMAP(NM,NP,WGTM,WGTP,DM,DP,DK,DX,SG,SRET,IFG,ERC)
C
C --- Note
C   Maximum a posteriori (or maximum likelihood) retrieval.
C
C --- history
C   2007.10.31 Made by M.Yamano (from VOLINV2 made by E.Kobayashi.)
C   2008.09.09 The definition of SG is changed.
C
C --- Input
C NM       I           number of measurements(y)
C NP       I           number of parameters(X)
C WGTM     R(KDT)      weight for measurements
C WGTP     R(KSZ)      weight for parameters
C DM       R(KDT)      y-f(Xk)
C DP       R(KSZ)      Xk-Xa
C DK       R(KDT,KSZ)  KERNEL MATRIX: (f(Xk+dX)-f(Xk))/dX
C IFG      I(KDT)      availability flag for measurements(>0: available)
C
C --- Output
C DX       R(KSZ)      retrieved dX=Xk+1 - Xk
C SG       R(KSZ)      RMSD of parameters
C SRET     R(KSZ)      reliability for parameters
C ERC      C*64        ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KDT   =1020)
      PARAMETER (KSZ   =100)
C
      CHARACTER ERC*(*)
      DIMENSION DM(KDT),DP(KSZ),WGTM(KDT),WGTP(KSZ),IFG(KDT)
      DIMENSION DK(KDT,KSZ),DX(KSZ),SG(KSZ),SRET(KSZ)
      DIMENSION A1(KSZ),A2(KSZ,KSZ),B(KSZ),C(KSZ,KSZ),SM(KSZ,KSZ)
C
      ERC=' '
      IF ((NP.GT.KSZ).OR.(NM.GT.KDT)) THEN
         ERC='INVMAP: Dimension size overflow'
         GOTO 900
      ENDIF
C
      DO J=1,NP
        S=0.
        DO K=1,NM
          IF (IFG(K).GT.0) S=S+DK(K,J)*(1/WGTM(K))*DM(K)
        ENDDO
        A1(J)=S-(1/WGTP(J))*DP(J)
        DO I=1,J
          S=0.
          DO K=1,NM
            IF (IFG(K).GT.0) S=S+DK(K,I)*(1/WGTM(K))*DK(K,J)
          ENDDO
          A2(I,J)=S
        ENDDO
        A2(J,J)=A2(J,J)+(1/WGTP(J))
      ENDDO
C
      G1=0
      DO I=1,NP
        IF (ABS(A2(I,I)).LE.0.0) THEN
           ERC='INVMAP: illegal kernel matrices'
           GOTO 900
        ENDIF
        G1=G1+LOG(ABS(A2(I,I)))
      ENDDO
      G1=EXP(G1/REAL(NP))
      DO J=1,NP
        B(J)=A1(J)/G1
        DO I=1,J
          C(I,J)=A2(I,J)/G1
        ENDDO
      ENDDO
C
      CALL TNVCH2(NP,C,DTT,KSZ,ERC)
      IF(ERC.NE.' ') GOTO 900
C ---
      DO J=1,NP
        DO I=1,J
          S=0.
          DO K=1,NP
            S=S+C(I,K)*(1/WGTP(K))*C(K,J)
          ENDDO
          SM(I,J)=S
        ENDDO
      ENDDO
C ---
      DO I=1,NP
        S=0.
        DO J=1,NP
          S=S+C(I,J)*B(J)
        ENDDO
        DX(I)=S
        SG(I)=SQRT(C(I,I)/G1/WGTP(I))
        SRET(I)=SQRT(SM(I,I)/C(I,I)/G1)
      ENDDO
C
  900 RETURN
      END
CC 10.02.15 M.Hashimoto      
      FUNCTION VLSPC2(PR)
CC      FUNCTION VLSPCN(R)
C
C Volume spectrum of partile polydisperison: v(x) = dV / d ln r
C
C--- HISTORY
C 88. 5. 6  REGISTERED BY T.NAKAJIMA
C 94.12.19  IF(R.LE.RMIN .OR. R.GE.RMAX) -> IF(R.LT.RMIN .OR. R.GT.RMAX)
C 02.05.15  new function for volume spectrum by M.Yamano.
C    *** the same contents as VOLSPC in /home/takka/Fortran/LBR1 ***
C    *** only the alignment of data in COMMON was changed        ***
C
C 08.09.05  The definition of S for log-normal(ITP=2) is changed.
CC 10.02.15 Modified by M.Hashimoto
C
C--- INPUT
C R        R        Particle radius in cm
C
C--- COMMON  /VSPC/
C RMIN     R        Minimum particle radius (cm)
C RMAX     R        Maximum particle radius (cm)
C PR    R(5,5)      Parameter packet
C                   PR(I,J)  I: parameters,  J: mode number.
C  PR(1,j): Type of function (ITP) for the mode.
C   ITP=1: power law
C     PR(2,j)=C, PR(3,j)=R0 (cm),  PR(4,j)=P
C     vj = C * (R/R0)**(4-P) if R>R0; = C * (R/R0)**4 if R<R0
C   ITP=2: log-normal
C     PR(2,j)=C, PR(3,j)=S ,  PR(4,j)=RM
C     vj = C * exp(-(ln(R/RM)/S)**2 / 2)
C   ITP=3: modified gamma
C     PR(2,j)=C, PR(3,j)=ALFA, PR(4,j)=BETA, PR(5,j)=GAMMA
C     vj = C * (R1)**(ALFA+4) exp (-BETA*R1**GAMMA) where R1=R*1.0E4
C NMODE    I        Number of modes
C                   v(x) = sum (n=1, NMODE) vn(x)
C
CC===
CC 10.02.15 M.Hashimoto
C$NAME  VLSPC2
C--- INPUT
C PR        R(6,4)      Parameter paket
C
C    PR(1,1)=R          Particle radius in cm
C    PR(1,2)=NMODE      Number of mode radius
C    PR(1,3)=RMIN       Minimum particle radius in cm
C    PR(1,4)=RMAX       Maximum particle radius in cm
C
C    For each j-th mode (<= 4)
C
C  PR(2,j): Type of function (ITP) for the mode.
C   ITP=1: power law
C     PR(3,j)=C, PR(4,j)=R0 (cm),  PR(5,j)=P
C     vj = C * (R/R0)**(4-P) if R>R0; = C * (R/R0)**4 if R<R0
C   ITP=2: log-normal
C     PR(3,j)=C, PR(4,j)=S ,  PR(5,j)=RM
C     vj = C * exp(-(ln(R/RM)/ln(S))**2 / 2)
C   ITP=3: modified gamma
C     PR(3,j)=C, PR(4,j)=ALFA, PR(5,j)=BETA, PR(6,j)=GAMMA
C     vj = C * (R1)**(ALFA+4) exp (-BETA*R1**GAMMA) where R1=R*1.0E4
C
C
C--- OUTPUT
C VLSPC2   RF       dV/d ln r
CC===
CC VLSPCN   RF       dV/d ln r
C$ENDI
C
      DIMENSION PR(6,4)
CC      COMMON /VSPC/RMIN,RMAX,PR(5,5),NMODE
C
CC==
      R     =PR(1,1)
      NMODE =PR(1,2)+0.001
      RMIN  =PR(1,3)
      RMAX  =PR(1,4)
      VLSPC2=0.0
      IF(R.LT.RMIN .OR. R.GT.RMAX) RETURN
      DO M=1,NMODE
        ITP=PR(2,M)+0.001
        IF(ITP.EQ.2) GOTO 4
        IF(ITP.EQ.3) GOTO 5
C POWER LAW
        C    =PR(3,M)
        RC   =PR(4,M)
        PDNDR=PR(5,M)
        IF(R.LE.RC) THEN
          PN=4.0
        ELSE
          PN=4.0-PDNDR
        ENDIF
        E1=PN*LOG(R/RC)
        GOTO 100
C LOG-NORMAL
    4   C =PR(3,M)
        S =PR(4,M)
        RM=PR(5,M)
        E1=-0.5*(LOG(R/RM)/LOG(S))**2
        GOTO 100
C MODIFIED GAMMA
    5   R1=R*1.0E4
        S  =PR(3,M)
        ALF=PR(4,M)
        BET=PR(5,M)
        GAM=PR(6,M)
        E1=(ALF+4)*LOG(R1)-BET*R1**GAM
  100   IF(E1.GT.-100.0) VLSPC2=VLSPC2+C*EXP(E1)
      ENDDO
      RETURN
      END
CC===
CC      VLSPCN=0.0
CC      IF ((R.LT.RMIN).OR.(R.GT.RMAX)) GOTO 900
C
CC      DO M=1,NMODE
CC        ITP=INT(PR(2,M)+0.1)
CC        ITP=INT(PR(1,M)+0.1)
C
C --- Power Law.
C
CC        IF (ITP.EQ.1) THEN
CC           C=PR(2,M)
CC           RC=PR(3,M)
CC           PDNDR=PR(4,M)
CC           IF (R.LE.RC) THEN
CC              PN=4.0
CC           ELSE
CC              PN=4.0-PDNDR
CC           ENDIF
CC           E1=PN*LOG(R/RC)
C
C --- Log - Normal.
C
CC        ELSE IF (ITP.EQ.2) THEN
CC           C=PR(2,M)
CC           S=PR(3,M)
CC           RM=PR(4,M)
c          E1=-0.5*(LOG(R/RM)/LOG(S))**2
CC           E1=-0.5*(LOG(R/RM)/S)**2
C
C --- Modified Gamma.
C
CC        ELSE IF (ITP.EQ.3) THEN
CC           R1=R*1.0E4
CC           C=PR(2,M)
CC           ALF=PR(3,M)
CC           BET=PR(4,M)
CC           GAM=PR(5,M)
CC           E1=(ALF+4)*LOG(R1)-BET*R1**GAM
CC        ENDIF
C
CC        IF (E1.GT.-100.0) VLSPCN=VLSPCN+C*EXP(E1)
CC      ENDDO
C
CC  900 RETURN
CC      END
C
      SUBROUTINE CLUTAU(NLN,CONCA,CONCM,TA1,WA1,TM1,WM1,NANG,PPA,PPM
     &                 ,THK,OMG,NLGN1,PHSF,TAU)
C
C   Calculation of total optical thickness.
C
C --- history
C   2000.12.15 Created by M.Yamano
C
C --- Input
C NLN     I              Number of Layer
C CONCA   R(KNLN)        Fraction of Aerosol   Optical Thickness
C CONCM   R(KNLN)        Fraction of Molecular Optical Thickness
C TA1     R              Aerosol Optical Thickness
C WA1     R              Aerosol single scattering albedo
C TM1     R              Molecular Optical Thickness
C WM1     R              Molecular single scattering albedo
C NANG    I              Number of scattering angles
C PPA     R(KNANG)       Normalized phase function of aerosols
C PPM     R(KNANG)       Normalized phase function of molecules
C
C --- Output
C THK     R(KNLN)        OPTICAL THICKNESS OF SUBLAYERS FROM TOP TO BOTTOM.
C OMG     R(KNLN)        SINGLE SCATTERING ALBEDO.
C NLGN1   I(KNLN)        MAXIMUM ORDER OF MOMENTS + 1.
C PHSF    R(KNANG,KNLN)  PHASE FUNCTION. GIVE WHEN INDP=1.
C TAU     R              total optical thickness
C ERC     C*64           ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNLN  =4)
      PARAMETER (KNANG =60)
C
      DIMENSION CONCA(KNLN),CONCM(KNLN),PPA(KNANG),PPM(KNANG)
      DIMENSION THK(KNLN),OMG(KNLN),NLGN1(KNLN),PHSF(KNANG,KNLN)
C
      TAU=0.
      DO L=1,NLN
        TAUA=TA1*CONCA(L)
        TAUM=TM1*CONCM(L)
        THK(L)=TAUA+TAUM
        TAU=TAU+THK(L)
        S1=WM1*TAUM+WA1*TAUA
        IF (THK(L).GT.0.) THEN
           OMG(L)=S1/THK(L)
        ELSE
           OMG(L)=1.
        ENDIF
        IF (S1.GT.0.) THEN
           DO J=1,NANG
             PHSF(J,L)=(WM1*TAUM*PPM(J)+WA1*TAUA*PPA(J))/S1
           ENDDO
        ELSE
           DO J=1,NANG
             PHSF(J,L)=-1.
           ENDDO
        ENDIF
        NLGN1(L)=150
      ENDDO
      RETURN
      END
C
      SUBROUTINE PRTRN1(IMTHD,NLN,THK,OMG,NLGN1,TAU,NANG,THETA,PHSF,TH0
     &              ,SOLID,TRO3,GALB,NA1U0,AM1U,NFI0,FAI,AI,FD,FDW1,ERC)
C
C
C
C --- history
C   2000.12.20 Created by M.Yamano
C
C --- Input
C
C --- Output
C ERC           C*64      ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNLN  =4)
      PARAMETER (KNANG =60)
      PARAMETER (KNA1  =50)
      PARAMETER (KNFI  =50)
      PARAMETER (KLGN1 =400)
      PARAMETER (KNA0  =1)
      PARAMETER (KNTAU =1)
      PARAMETER (KNDM  =16)
      PARAMETER (KPLK1 =2)
      PARAMETER (KNA1U=KNA1)
C
      PARAMETER (PI=3.141592653589793,RAD=PI/180.0)
C
      CHARACTER ERC*(*)
      DIMENSION THETA(KNANG)
      DIMENSION THK(KNLN),OMG(KNLN),NLGN1(KNLN),PHSF(KNANG,KNLN)
      DIMENSION AM1U(KNA1),FAI(KNFI)
      DIMENSION AM0(KNA0),G(KLGN1,KNLN),CPLK(KPLK1,KNLN),UTAU(KNTAU)
     &  ,FLXD(KNA0,KNTAU),FLXU(KNA0,KNTAU),AI(KNA1U,KNA0,KNFI,KNTAU)
C
      ERC=' '
      INDA=1
      INDT=1
      INDP=1
      NDA=8
      EPSP=0
      EPSU=0
      FSOL=1
      NPLK1=0
      BGND=0
      AM00=COS(TH0*RAD)
      NA0=1
      AM0(NA0)=AM00
      EM=1.0/AM00
      NTAU=1
      UTAU(NTAU)=TAU*0.9999
C
      RDRD=SQRT(SOLID/PI)/RAD
      NRDRD=5
      NA1U=NA1U0
      NFI=NFI0
      CALL DFCRCT(TH0,NRDRD,RDRD,NA1U,AM1U,NFI,FAI,ERC)
      IF (ERC.NE.' ') GOTO 900
C
      CALL RTRN1(INDA,INDT,INDP,IMTHD,NDA,NA1U,AM1U,NA0,AM0
     &          ,NFI,FAI,NLN,THK,OMG,NLGN1,G,NANG,THETA,PHSF
     &          ,EPSP,EPSU,GALB,FSOL,NPLK1,CPLK,BGND,NTAU,UTAU
     &          ,FLXD,FLXU,AI,ERC)
      IF (ERC.NE.' ') GOTO 900
C
      CALL CORRR(TH0,AI,NA1U,NFI,NRDRD,RDRD,ES)
      FD=(EXP(-EM*UTAU(NTAU))+ES)*TRO3
      FDW1=(FD-ES*TRO3)/FD
C
  900 RETURN
      END
C
      SUBROUTINE LGNDF1(LMAX1,N,X,Y,ALF)
C------
C LEGENDRE EXPANSION
C------ HISTORY
C 86.12.16 CREATED
C 87.11. 5 FOR IBM/MVS
C 89. 8.20 FOR MAC WITH SINGLE PRECISION
C 90.11. 5 X IS DOUBLE PRECISION AGAIN, FOR
C          SIGLE PRECISION, PLEASE USE LGNDF2.
C------ VARIABLES
C  TYPE = I(INTEGER), C(CHARACTER), R(REAL)
C  TYPE ALSO SHOWS THE ARRAY STRUCTURE AS I(N)
C------ INPUT VARIABLES
C VARIABLE  TYPE    INTERPRETATION
C LMAX1      I      MAXIMUM ORDER + 1.
C N          I      NUMBER OF DATA ON (-1, 1).
C X        D(NA)    INDEPENDENT VARIABLES ON (-1, 1)
C Y        R(NA)    Y(X).
C------ OUTPUT VARIABLES
C ALF    R(LMAX1)   COEFFICIENTS OF LEGENDRE EXPANTION.
C$ENDI
C------ VARIABLES FOR THE ROUTINE.
      DOUBLE PRECISION X(N)
      DIMENSION Y(N),ALF(LMAX1)
C------ WORKING AREAS.
      PARAMETER (PI=3.141592653589793, RAD=PI/180.0)
      PARAMETER (EXMX=100.0)
      PARAMETER (NG=5)
      DOUBLE PRECISION X1,X2,X3,X4,XX
      DIMENSION GX(NG),GW(NG)
      DATA INIT/1/
C------ SHIFTED GAUSSIAN QUADRATURE.
      IF(INIT.GE.1) THEN
      INIT=0
      CALL QGAUSN(GW,GX,NG)
      ENDIF
C------- EXPANSION.
      LMAX=LMAX1-1
      DO 15 L1=1,LMAX1
   15 ALF(L1)=0.0
      DO 17 I=1,N-1
      IF(I-2) 18,18,19
   18 I1=1
      I4=4
      GO TO 20
C *
   19 IF(I-N+3) 21,22,22
   21 I1=I-1
      I4=I+2
      GO TO 20
C *
   22 I1=N-3
      I4=N
   20 I2=I1+1
      I3=I2+1
      X1=X(I1)
      X2=X(I2)
      X3=X(I3)
      X4=X(I4)
      IF(Y(I1).LE.0..OR.Y(I2).LE.0..OR.Y(I3).LE.0..OR.
     & Y(I4).LE.0.) GO TO 91
      ISIGN=1
      ALP1=LOG(Y(I1))
      ALP2=LOG(Y(I2))
      ALP3=LOG(Y(I3))
      ALP4=LOG(Y(I4))
      GO TO 92
   91 ALP1=Y(I1)
      ALP2=Y(I2)
      ALP3=Y(I3)
      ALP4=Y(I4)
      ISIGN=-1
   92 CONTINUE
      DO 23 J=1,NG
      XX=X(I+1)+GX(J)*(X(I)-X(I+1))
      WW=GW(J)*(X(I)-X(I+1))
      PP=(XX-X2)*(XX-X3)*(XX-X4)/(X1-X2)/(X1-X3)/(X1-X4)*ALP1
     &  +(XX-X1)*(XX-X3)*(XX-X4)/(X2-X1)/(X2-X3)/(X2-X4)*ALP2
     &  +(XX-X1)*(XX-X2)*(XX-X4)/(X3-X1)/(X3-X2)/(X3-X4)*ALP3
     &  +(XX-X1)*(XX-X2)*(XX-X3)/(X4-X1)/(X4-X2)/(X4-X3)*ALP4
      IF(ISIGN.EQ.-1) GO TO 109
      IF(PP.GT.-EXMX) THEN
      PP=EXP(PP)
      ELSE
      PP=0.
      ENDIF
  109 DO 24 L1=1,LMAX1
      L=L1-1
      IF(L-1) 26,27,28
   26 PL=1.
      PL2=PL
      GO TO 29
C *
   27 PL=XX
      PL1=PL
      GO TO 29
C *
   28 PL=(REAL(2*L-1)*XX*PL1-REAL(L-1)*PL2)/REAL(L)
      PL2=PL1
      PL1=PL
   29 ALF(L1)=ALF(L1)+PP*PL*WW
   24 CONTINUE
   23 CONTINUE
   17 CONTINUE
      DO 41 L1=1,LMAX1
   41 ALF(L1)=ALF(L1)*REAL(2*L1-1)/2.0
      RETURN
      END
      SUBROUTINE PCOR(COR,THI,THO,PHI,TAUR,SSAR,DPF,
     & TAURO,TAUA,SSAA,P112,P122)
C CORRECTING POLARIZATION EFFECT TO THE TOTAL INTENSITY.
C--- HISTORY
C 87. 7. 9  FROM COR87 BY OGAWA, COR/PI
C 88. 4. 6  ADD IF-STATEMENT FOR ZNTI,ZNTO=60
C           INTRODUCE ZZ FOR OPT(3).
C--- INPUT
C THI   R   INCIDENT ZENITH ANGLE (DEGREE)
C THO   R   EMERGENT ZENITH ANGLE (DEGREE)
C PHI    R   AZIMUTHAL ANGLE (DEGREE).
C TAUR   R   RAYLEIGH OPTICAL THICKNESS OF THE LAYER.
C SSAR   R   RAYLEIGH SINGLE SCATTERING ALBEDO.
C DPF    R   DEPOLARIZATION FACTOR (RHO).
C TAURO R    RAYLEIGH OPTICAL DEPTH AT OBSERVING LEVEL.
C TAUA   R   AEROSOL OPTICAL THICKNESS OF THE LAYER.
C SSAA   R   AEROSOL SINGLE SCATTERING ALBEDO.
C P112   R   BETA2 = ALFA2 = 5*G2 (LEGENDRE COEFF. OF PHASE FUNCTION)
C P122   R   = -GAMMA2
C--- OUTPUT
C COR  R(2)  CORRECTION FOR INTENSITY   1(DOWN), 2(UP)
C            SPECIFIC INTENSITY (F0=1) + COR
C--- CPU TIME
C 0.40 MSEC PER CALL ON IBM/3081
C$ENDI
      PARAMETER (PI=3.141592653, RD=PI/180.0)
      DIMENSION DA(3),DB(3),CF(3),ETA(4),EF(4,2),COR(2)
C ANGLES <> 60
      IF(ABS(THI-60.0).GE.0.5) THEN
      ZNTI=THI
      ELSE
      ZNTI=59.5
      ENDIF
      IF(ABS(THO-60.0).GE.0.5) THEN
      ZNTO=THO
      ELSE
      ZNTO=59.5
      ENDIF
C
      F=(2.0*P112-1.0)/9.0
      SA=(1.0-F)*SSAA*TAUA
      TA=(1-F*SSAA)*TAUA
      SR=SSAR*TAUR
      TR=TAUR
      TAU=TR+TA
      TAUSCA=SR+SA
      SSA=TAUSCA/TAU
C
      C=(1.0-DPF)/(1.0+DPF/2.0)
      A=(C*SR+SA)/TAUSCA
      B=(C*SR+SQRT(2.0/3.0)*P122/(1.0-F)*SA)/TAUSCA/A
      TAU1=TAU*TAURO/TAUR
C
      C=SQRT(TAU)
      DA(1)=0.6*C
      DA(2)=0.5*C
      DA(3)=0.4*C
      DB(1)=0.1*TAU
      DB(2)=0.1*TAU
      DB(3)=0.1*TAU
      DO 1 M1=1,3
    1 CF(M1)=(SSA*A*B)**2*(1.0-DA(M1))/(1.0-SSA*A*DA(M1))
     &  *(1.0-DB(M1))/(1.0-B**2*DB(M1))
C
      X000= 0.0085*TAU
      X002=-0.0415*TAU**0.52
      X022=0.25
      X122=0.36
      X222=0.195
      ETA(1)=0.5
      ETA(2)=0.514-0.373*EXP(-0.826*TAU)
      ETA(3)=0.752-0.581*EXP(-0.584*TAU)
      ETA(4)=1.248-1.160*EXP(-0.843*TAU)
C
      AM =COS(ZNTI*RD)
      P20E=(3.0*AM**2-1.0)/2.0
      P21E=SQRT(3.0/2.0)*AM*SQRT(1.0-AM**2)
      P22E=SQRT(3.0/8.0)*(1.0-AM**2)
C
      AM0=COS(ZNTO*RD)
      P20S=(3.0*AM0**2-1.0)/2.0
      P21S=SQRT(3.0/2.0)*AM0*SQRT(1.0-AM0**2)
      P22S=SQRT(3.0/8.0)*(1.0-AM0**2)
C
      DO 3 IE=1,4
      E =ETA(IE)
      EN=-E
      EXP1=EXP(-TAU/AM0-TAU/E)
      CALL GFCD(AM,AM0,TAU1,TAU,C1,D1)
      CALL GFCD(AM,E,TAU1,TAU,C2,D2)
      CALL GFCD(AM,EN,TAU1,TAU,C3,D3)
      EF(IE,1)=AM0*E*((C1-C2)/(AM0-E)+
     &                (C1-C3*EXP1)/(AM0+E))
      EF(IE,2)=AM0*E*((D1-D2)/(AM0-E)
     &               +(D1-D3*EXP1)/(AM0+E))
    3 CONTINUE
C
      DO 20 IUD=1,2
      DI0= ((SSA**2*X000+SSA*X002*(P20E+P20S))*EF(1,IUD)
     &      +X022*P20E*P20S*EF(2,IUD))*CF(1)
      DI1= X122*P21E*P21S*EF(3,IUD)*CF(2)
      IF(IUD.EQ.2) DI1=-DI1
      DI2= X222*P22E*P22S*EF(4,IUD)*CF(3)
C /PI -> FOR UNIT OF F0=1.
      COR(IUD)=(DI0+DI1*COS(PHI*RD)+DI2*COS(2.0*PHI*RD))/PI
   20 CONTINUE
C
      RETURN
      END
      SUBROUTINE GFCD(AM,AM0,TAU1,TAU0,C,D)
C GEOMETRICAL FACTORS C AND D
C AM, AM0 > 0
      EXP1=EXP(-TAU1/AM0)
      D=AM0/(AM0+AM)*(EXP1-EXP(TAU1/AM-TAU0/AM0-TAU0/AM))
      IF(ABS(AM-AM0).LE.1.0E-5) THEN
      C=TAU1/AM0*EXP1
      RETURN
      ELSE
      C=AM0/(AM0-AM)*(EXP1-EXP(-TAU1/AM))
      RETURN
      ENDIF
      END
      SUBROUTINE  TNVCH2(N,A,DT,NN,ERR)
C INVERSION OF SYMMETRIC POSITIVE DEFINITE MATRIX
C     SQUARE ROOT METHOD
C     MATRIX A IS SYMMETRIC AND POSITIVE DEFINITE
C--- HISTORY
C 71. 4.30 CREATED BY SAKATA MASATO
C 89.11.10 ADDED ERR.
C 90. 1.17  REGISTERED
C--- INPUT
C N        I        DIMENSION OF THE MATRIX
C A      R(NN,N)    MATRIX
C NN       I        SIZE OF THE FIRST ARGUMENT OF A
C--- OUTPUT
C A                 INVERSE OF THE MATRIX
C DT       R        DETERMINATION OF THE MATRIX
C ERR    C*64       Error code if ' ' then NORMAL TERMINATION
C$ENDI
      CHARACTER ERR*(*)
      DIMENSION    A(NN,N)
      ERR=' '
      IF(N-1)   1,2,3
    1 ERR='ERROR IN TINVCH: N .LE.0'
      RETURN
    2 DT=A(1,1)
      A(1,1)=1.0/A(1,1)
      IF(DT.LE.0.0)   GO TO  60
      RETURN
    4 ERR='ERROR IN TINVCH: N.GT.NN'
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
   60 ERR='ERROR IN TINVCH: GIVEN MATRIX IS NOT POSITIVE DEFINITE'
      RETURN
      END
      SUBROUTINE DFCRCT(TH0,NRDRD,RDRD,NA1U,AM1U,NFI,FAI,ERC)
C
C   Diffuse correction.
C
C --- history
C   2000.12.15 Created by M.Yamano
C
C --- Input
C TH0        R           SOLAR ZENITH ANGLE IN DEGREE (0,90)
C NRDRD      I
C RDRD       R
C NA1U       I           NUMBER OF EMERGENT NADIR ANGLES IN SPHERE.
C AM1U       R(KNA1U)    CONSINE OF THE EMERGENT NADIR ANGLES.
C NFI        I           NUMBER OF AZIMUTHAL ANGLES.
C FAI        R(KNFI)     AZIMUTHAL ANGLES IN DEGREES.
C
C --- Output
C NA1U       I           NUMBER OF EMERGENT NADIR ANGLES IN SPHERE.
C AM1U       R(KNA1U)    CONSINE OF THE EMERGENT NADIR ANGLES.
C NFI        I           NUMBER OF AZIMUTHAL ANGLES.
C FAI        R(KNFI)     AZIMUTHAL ANGLES IN DEGREES.
C ERC        C*64        ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNA1  =50)
      PARAMETER (KNFI  =50)
C
      PARAMETER (PI=3.141592653589793,RAD=PI/180.0)
      CHARACTER ERC*(*)
      DIMENSION AM1U(KNA1),FAI(KNFI)
C
      ERC=' '
      AM00=COS(TH0*RAD)
      IF (TH0.GT.RDRD/2) THEN
         NA1U=NA1U+1
         IF (NA1U.GT.KNA1) THEN
            ERC='DFCRCT: EKNA1'
            GOTO 900
         ENDIF
         AM1U(NA1U)=AM00
         DO I=1,NRDRD
           NFI=NFI+1
           IF (NFI.GT.KNFI) THEN
              ERC='DFCRCT: EKNFI'
              GOTO 900
           ENDIF
           A2=(COS(RDRD*(I-0.5)/REAL(NRDRD)*RAD)-AM00**2)/(1-AM00**2)
           FAI(NFI)=ACOS(A2)/RAD
         ENDDO
      ELSE
         NFI=NFI+1
         IF (NFI.GT.KNFI) THEN
            ERC='DFCRCT: EKNFI'
            GOTO 900
         ENDIF
         FAI(NFI)=0
         DO I=1,NRDRD
           NA1U=NA1U+1
           IF (NA1U.GT.KNA1) THEN
              ERC='DFCRCT: EKNA1'
              GOTO 900
           ENDIF
           AM1U(NA1U)=COS((TH0+RDRD*(I-0.5)/REAL(NRDRD))*RAD)
         ENDDO
      ENDIF
C
  900 RETURN
      END
C
      SUBROUTINE RTRN1(INDA,INDT,INDP,IMTHD,NDA,NA1U,AM1U,NA0,AM0
     &,NFI,FI,NLN,THK,OMG,NLGN1,G,NANG,ANG,PHSF,EPSP,EPSU,GALB,FSOL
     &,NPLK1,CPLK,BGND,NTAU,UTAU,FLXD,FLXU,AI,ERR)
C
C SOLVE THE RADIATIVE TRANSFER IN ATMOSPHERE SYSTEM.
C  -USE FTRN1-
C BY TERUYUKI NAKAJIMA
C
C--- HISTORY
C 89.11. 7 CREATED
C 90.11.22 IF(INDP.LE.0 .AND. (INDA.GE.1 .AND. IMTHD.GE.1)) INDPL=1 ->
C                  EQ
C          CHANGE CODES BETWEEN -DO 4 - AND -DO 8 -;
C          IF(ABS(CS1).GT.1.0) CS1=SIGN(1.0,CS1); CHANGE LOOP-13,14
C          ADDING LOOP-47.
C          (Effected when IMTHD=3 and INDP=1)
C 95. 4. 8 Change note for CPLK
C--- INPUT
C INDA       I       0: FLUX ONLY.
C                    1: INTENSITY USING -AM1U-.
C                    2: INTENSITY AT GAUSSIAN QUADRATURE POINTS.
C                       -NA1U- AND -AM1U- ARE SET AUTOMATICALLY AS
C                       2*NDA AND (-AMUA, +AMUA).
C INDT       I       1: INTENSITY/FLUX AT USER DEFINED DEPTH -UTAU-.
C                    2: INTENSITY/FLUX AT SUBLAYER INTERFACES.
C                       -NTAU- AND -UTAU- ARE SET AUTOMATICAALY AS
C                       -NLN1- AND -DPT-.
C INDP       I      -1: GIVE -G- AND USE -G- FOR INTENSITY INTERPOLATION
C                       FROM GAUSSIAN POINTS TO USER POINTS.
C                    0: GIVE -G- AND CONSTRUCT PHASE FUNCTION FOR
C                       INTENSITY INTERPOLATION.
C                    1: GIVE PHASE FUNCTION -PHSF-.
C IMTHD      I      -1: NT,  0: DMS-METHOD  FOR INTENSITY/FLUX
C                    1: MS,  2:TMS,  3:IMS-METHOD FOR INTENSITY.
C NDA        I       NUNBER OF NADIR-QUADRATURE ANGLES IN THE
C                    HEMISPHERE.
C NA1U       I       NUMBER OF EMERGENT NADIR ANGLES IN SPHERE.
C AM1U    R(KNA1U)   CONSINE OF THE EMERGENT NADIR ANGLES.
C                      - FOR UPWARD, + FOR DOWNWARD.
C NA0        I       NUMBER OF SOLAR INCIDENCES.
C AM0     R(NA0)     CONSINE OF SOLAR ZENITH ANGLES .GT. 0.
C NFI        I       NUMBER OF AZIMUTHAL ANGLES.
C FI     R(KNFI)     AZIMUTHAL ANGLES IN DEGREES.
C NLN        I       NUMBER OF ATMOSPHERIC SUBLAYERS.
C THK     R(KNLN)    OPTICAL THICKNESS OF SUBLAYERS FROM TOP TO BOTTOM.
C OMG     R(KNLN)    SINGLE SCATTERING ALBEDO.
C NLGN1   I(KNLN)    MAXIMUM ORDER OF MOMENTS + 1.
C                      GIVE WHEN -INDP- .LE. 0.
C                      GIVE WHEN IMTHD=3 REGARDLESS OF -INDP-.
C                      OTHERWISE, A VALUE .LE. 2*NDA+1 IS GIVEN
C                      AUTOMATICALLY BY THE ROUTINE.
C G   R(KLGN1,KNLN)  LEGENDRE MOMENTS OF PHASE FUNCTION.
C                      GIVE WHEN INDP .LE. 0.
C NANG       I       NUMBER OF SCATTERING ANGLES FOR PHASE FUNCTIONS.
C                      GIVE WHEN INDP=1.
C                      GIVE WHEN INDP=0 .AND. (INDA.GE.1) .AND.
C                                             (IMTHD.GE.1).
C ANG     R(KNANG)   SCATTERING ANGLES IN DEGREES.
C                      GIVE WHEN INDP=1.
C                      GIVE WHEN INDP=0 .AND. (INDA.GE.1) .AND.
C                                             (IMTHD.GE.1).
C                      ANGLES SHOULD HAVE ENOUGH RESOLUTION FOR LEGENDRE
C                      EXPANSION (INDP=1) AND INTERPOLATION BY CUBIC
C                      POLYNOMIAL.
C PHSF    R(KNANG,   PHASE FUNCTION. GIVE WHEN INDP=1.
C EPSP       R       TRUNCATION CRITERION OF PHASE FUNCTION MOMENTS.
C EPSU       R       CONVERGENCE CRITERION OF INTENSITY (INDA.GT.0)
C GALB       R       GROUND ALBEDO (LAMBERTIAN).
C FSOL       R       SOLAR IRRADIANCE AT THE SYSTEM TOP.
C NPLK1      I       NUMBER OF ORDER TO APPROXIMATE PLANK + 1.
C                      IF 0 THEN NO THERMAL.
C CPLK R(KPLK1,KNLN) Plank function
C                      = SUM(K1=1,NPLK1) CPLK(K1,L) * TAU**(K1-1)
C                      TAU IS THE OPTICAL DEPTH MEASURED FROM THE TOP
C                      OF THE SUBLAYER (SAME UNIT AS FSOL).
C BGND       R       THERMAL EMISSION FROM THE GROUND
C                    = (1-GALB)*Plank function
C                      (SAME UNIT AS FSOL).
C NTAU       I       NUMBER OF USER DEFINED OPTICAL DEPTHS.
C UTAU    R(KNTAU)   OPTICAL DEPTHS WHERE THE FIELD IS CALCULTED.
C                      TOP TO BOTTOM.
C--- OUTPUT
C NA1U, AM1U         SAME AS 2*NDA, (-AUMA, +AMUA) WHEN INDA=2.
C NLGN1, G           NORMALIZED BY G(1), AND CUT OFF BY THE CRITERION
C                     -EPSP-.
C NTAU, UTAU         SAME AS NLN+1 AND -DPT- WHEN INDT=1.
C PHSF               NORMALIZED PHASE FUNCTION (1 OVER UNIT SPHERE).
C                    RECONSTRUCTED WHEN (INDA.GE.1) .AND. (IMTHD.GE.1)
C                    .AND. (INDP=0).
C FLXD    R(KNA0,    DOWNWARD FLUX AT -UTAU-.
C           KNTAU)
C FLXU               UPWARD   FLUX AT -UTAU-.
C AI      R(KNA1U,   I(MU1U(I), MU0(J), FI(K), UTAU(IT))
C           KNA0,    INTENSITY AT -UTAU- (WATT/M2/MICRON).
C           KNFI,
C           KNTAU)
C ERR      C*64      ERROR INDICATER. IF " " THEN NORMAL END.
C--- CORRESPONDENCE BETWEEN VARIABLES AND PARAMETERS
C KNA1U       NA1U
C KNA1        NA1
C KNA0        NA0
C KNDM        NDA
C KNFI        NFI
C KNLN        NLN
C KNLN1       KNLN+1
C KNTAU       NTAU
C KLGN1       MXLGN1 = MAX(NLGN1)
C              THIS SHOULD BE VERY LARGE WHEN IMS-METHD IS USED
C              SO AS TO APPROXIMATE P**2, OR WHEN INDP=1 IS ASSIGNED
C              TO RECONSTRUCT THE ORIGINAL PHASE FUNCTION FROM -G-.
C KLGT1       MXLGT1 = SAME AS MXLGN1 BUT FOR THAT OF TRUNCATED
C              VARIABLES
C KPLK1       NPLK1
C$ENDI
C PARAMETERS
      PARAMETER (KNA1U =50)
      PARAMETER (KNA1  =50)
      PARAMETER (KNA0  =1)
      PARAMETER (KNDM  =16)
      PARAMETER (KNFI  =50)
      PARAMETER (KLGN1 =400)
      PARAMETER (KNLN  =4)
      PARAMETER (KNTAU =1)
      PARAMETER (KPLK1 =2)
      PARAMETER (KNANG =60)
      PARAMETER (NCHK1 =3)
      PARAMETER (NCHK2 =3)
C
      PARAMETER (KNLN1=KNLN+1,KLGT1=2*KNDM)
      PARAMETER (PI=3.141592654,RAD=PI/180.0)
C AREAS FOR THIS ROUTINE
      CHARACTER ERR*(*)
      DIMENSION AM1U(KNA1U),AM0(KNA0),FI(KNFI),THK(KNLN),OMG(KNLN)
     &,NLGN1(KNLN),G(KLGN1,KNLN),ANG(KNANG),PHSF(KNANG,KNLN)
     &,CPLK(KPLK1,KNLN),UTAU(KNTAU),FLXD(KNA0,KNTAU)
     &,FLXU(KNA0,KNTAU),AI(KNA1U,KNA0,KNFI,KNTAU)
C AREAS FOR FTRN1
      DIMENSION AMUA(KNDM),WA(KNDM),GT(KLGT1,KNLN),NLGT1(KNLN)
     &,CPLKT(KPLK1,KNLN),OMGT(KNLN),THKT(KNLN),UTAUT(KNTAU)
     &,AIF(KNA1U,KNA0,KNTAU)
C WORKING AREAS
      DIMENSION CS(KNANG),YY(KNANG),GG(KLGN1),COSM(KNFI)
     &,PHS(KNLN),FF(KNLN),DPT(KNLN1),DPTT(KNLN1)
     &,PHST(KNLN),SGL2(KNTAU),PHSB(KNLN),COR(KNTAU)
       DIMENSION IM1U(KNA1U),JM0(KNA0),BM1U(KNA1U),BM0(KNA0)
     &,EU1(KNA1U),EU0(KNA0),IICHK0(KNA0),IICHK1(KNA1U)
C CHECK AND SET OF VARIABLES
      CALL CHKRT(INDA,INDT,INDP,IMTHD,NDA,NA1U,AM1U,NA0,AM0
     &,NFI,FI,NLN,THK,OMG,NLGN1,G,NANG,ANG,PHSF,EPSP,EPSU,GALB,FSOL
     &,NPLK1,CPLK,BGND,NTAU,UTAU,AMUA,WA,DPT,MXLGN2,ERR)
      IF(ERR.NE.' ') RETURN
C FOURIER EXPANSION PARAMETERS
      NDA2=2*NDA
      NDA21=NDA2+1
      MXLGT1=NDA2
      IF(MXLGT1.GT.KLGT1) THEN
        ERR='-KLGT1- IS TOO SMALL'
        RETURN
      ENDIF
CC FLUX
      IF(INDA.LE.0) THEN
        IF(IMTHD.LT.0) THEN
          MXLGN1=MIN(NDA2,MXLGN2)
         ELSE
          MXLGN1=MIN(NDA21,MXLGN2)
        ENDIF
CC INTENSITY
       ELSE
        IF(INDP.LE.0) THEN
          MXLGN1=MXLGN2
         ELSE
          IF(IMTHD.LE.-1) THEN
            MXLGN1=NDA2
           ELSE
            IF(IMTHD.LE.2) THEN
              MXLGN1=NDA21
             ELSE
              MXLGN1=MXLGN2
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      IF(MXLGN1.GT.KLGN1) THEN
        ERR='-KLGN1- IS TOO SMALL'
        RETURN
      ENDIF
C INDICATER FOR RE-CALCULATING PHASE FUNCTION USING -G-
      INDPL=0
      IF(INDP.EQ.0 .AND. (INDA.GE.1 .AND. IMTHD.GE.1)) INDPL=1
C CLEAR VARIABLES
      IF(INDA.GE.1) THEN
        DO 1 IT=1,NTAU
        DO 1 J=1,NA0
        DO 1 K=1,NFI
        DO 1 I=1,NA1U
    1   AI(I,J,K,IT)=0
      ENDIF
C SET COS(SCATTERING ANGLE)
      IF((INDP.EQ.1) .OR. (INDPL.GE.1)) THEN
        DO 2 I=1,NANG
    2   CS(I)=COS(ANG(I)*RAD)
      ENDIF
C LOOP FOR SUBLAYERS
      DO 3 L=1,NLN
CC LEGENDRE EXPANSION OF THE PHASE FUNCTIONS
CC  WE TRUNCATE THE SERIES WHEN ABS(G) BECOMES SMALLER THAN -EPSP-
CC   FOR SUCCESSIVE THREE TIMES.
        IF(INDP.EQ.1) THEN
          DO 4 I=1,NANG
    4     YY(I)=PHSF(I,L)
          CALL LGNDF3(MXLGN1,NANG,CS,YY,GG)
          GG0=GG(1)
          ICHK=0
          DO 5 K1=1,MXLGN1
          G(K1,L)=GG(K1)/GG0
          IF(ABS(G(K1,L)).LE.EPSP) THEN
            ICHK=ICHK+1
            IF(ICHK.GE.NCHK1) GOTO 6
          ENDIF
    5     CONTINUE
          K1=MXLGN1
    6     NLGN1(L)=K1
C CS IS FROM 1 TO -1 FOR LGNDF3, SO THAT GG0<0 (90.11.22)
          GG0=ABS(GG0)
          DO 10 I=1,NANG
   10     PHSF(I,L)=PHSF(I,L)/GG0/4/PI
CC CONSTRUCTION OF PHASE FUNCTION FROM -G- WHEN INDP=0 AND IMTHD.GE.1
         ELSE
          GG0=G(1,L)
          ICHK=0
          K2=NLGN1(L)
          DO 9 K1=1,K2
          G(K1,L)=G(K1,L)/GG0
          IF(ABS(G(K1,L)).LE.EPSP) THEN
            ICHK=ICHK+1
            IF(ICHK.GE.NCHK1) GOTO 46
          ENDIF
    9     CONTINUE
          K1=MXLGN1
   46     NLGN1(L)=K1
          IF(INDPL.GE.1) THEN
            DO 7 I=1,NANG
            CS1=CS(I)
            INIT=1
            SUM=0
            DO 8 K1=1,NLGN1(L)
    8       SUM=SUM+(2*K1-1)*G(K1,L)*PLGD(INIT,CS1)
    7       PHSF(I,L)=SUM/4/PI
          ENDIF
        ENDIF
CC DELTA-M MOMENTS (ICHKT=1 MEANS WE DO TRUNCATION AT LEAST ONE TIME)
        ICHKT=0
        NLGT1(L)=MIN(MXLGT1,NLGN1(L))
        IF(IMTHD.GE.0 .AND. NLGN1(L).GT.MXLGT1) THEN
          FF(L)=G(MXLGT1+1,L)
          ICHKT=1
         ELSE
          FF(L)=0
        ENDIF
        DO 11 K1=1,NLGT1(L)
   11   GT(K1,L)=(G(K1,L)-FF(L))/(1-FF(L))
        OMGT(L)=(1-FF(L))*OMG(L)/(1-FF(L)*OMG(L))
        THKT(L)=(1-FF(L)*OMG(L))*THK(L)
    3 CONTINUE
C RESET -MXLGN1- AND -MXLGT1-
      MXLGN1=1
      MXLGT1=1
      DO 12 L=1,NLN
      MXLGN1=MAX(MXLGN1,NLGN1(L))
   12 MXLGT1=MAX(MXLGT1,NLGT1(L))
C SET TRUNCATED USER DEFINED DEPTHS
      DPTT(1)=DPT(1)
      DO 13 L=1,NLN
   13 DPTT(L+1)=DPTT(L)+THKT(L)
      DO 14 IT=1,NTAU
      DO 47 L=1,NLN
      IF((UTAU(IT)-DPT(L))*(UTAU(IT)-DPT(L+1)).LE.0) THEN
        DTAU=UTAU(IT)-DPT(L)
        DTAUT=(1-FF(L)*OMG(L))*DTAU
        UTAUT(IT)=DPTT(L)+DTAUT
      ENDIF
   47 CONTINUE
      IF(UTAUT(IT).GT.DPTT(NLN+1)) UTAUT(IT)=DPTT(NLN+1)
   14 CONTINUE
C SCALING COEFFICIENTS FOR THERMAL EMISSION
      IF(NPLK1.GT.0) THEN
        DO 15 L=1,NLN
        DO 15 K1=1,NPLK1
   15   CPLKT(K1,L)=CPLK(K1,L)/(1-OMG(L)*FF(L))**(K1-1)
      ENDIF
C SET -MMAX1-
      MMAX1=1
      IF(INDA.GE.1) MMAX1=MXLGT1
C INITIALIZE ANGLE INDEX FOR CHECKING CONVERGENCE
        NB0=NA0
        DO 24 J=1,NA0
        IICHK0(J)=0
        JM0(J)=J
   24   BM0(J)=AM0(J)
      IF(INDA.GT.0) THEN
        NB1U=NA1U
        DO 25 I=1,NA1U
        IICHK1(I)=0
        IM1U(I)=I
   25   BM1U(I)=AM1U(I)
      ENDIF
C FOURIER SUM
      INIT=1
      DO 16 M1=1,MMAX1
        M=M1-1
        FM=PI
        IF(M.EQ.0) FM=2*PI
        CALL FTRN1(INIT,M,INDA,INDT,NB1U,BM1U,NB0,BM0
     &  ,NDA,AMUA,WA,NLN,NLGT1,GT,DPTT,OMGT,NPLK1,CPLKT,GALB,BGND
     &  ,FSOL,NTAU,UTAUT,FLXD,FLXU,AIF,ERR)
        IF(ERR.NE.' ') RETURN
        IF(INDA.GT.0) THEN
          DO 44 JB=1,NB0
   44     EU0(JB)=0
          DO 45 IB=1,NB1U
   45     EU1(IB)=0
          DO 17 K=1,NFI
          FI1=M*FI(K)*RAD
   17     COSM(K)=COS(FI1)/FM
          DO 18 JB=1,NB0
          J=JM0(JB)
          DO 18 IB=1,NB1U
          I=IM1U(IB)
          DO 18 K=1,NFI
          DO 18 IT=1,NTAU
          DAI=AIF(IB,JB,IT)*COSM(K)
          AI(I,J,K,IT)=AI(I,J,K,IT)+DAI
          IF(ABS(AI(I,J,K,IT)).GT.0) THEN
            EAI=ABS(DAI/AI(I,J,K,IT))
           ELSE
            IF(ABS(DAI).LE.0) THEN
              EAI=0
             ELSE
              EAI=100
            ENDIF
          ENDIF
          EU0(JB)=MAX(EU0(JB),EAI)
          EU1(IB)=MAX(EU1(IB),EAI)
   18     CONTINUE
C CHECK CONVERGENCE FOR AM0
          CALL CONVU(NB0,BM0,JM0,EU0,IICHK0,EPSU,NCHK2)
          IF(NB0.LE.0) GOTO 30
          IF(INDA.NE.2) THEN
            CALL CONVU(NB1U,BM1U,IM1U,EU1,IICHK1,EPSU,NCHK2)
            IF(NB1U.LE.0) GOTO 30
          ENDIF
C
        ENDIF
   16 CONTINUE
   30 IF(IMTHD.LE.0 .OR. INDA.LE.0 .OR. ICHKT.LE.0) RETURN
C--- INTENSITY CORRECTION.
        DO 19 IS=1,NA0
        DO 19 IU=1,NA1U
        DO 19 IZ =1,NFI
CC COS(SCATTERING ANGLE)
        CS1=AM0(IS)*AM1U(IU)+SQRT((1-AM0(IS)**2)*(1-AM1U(IU)**2))
     &    *COS(FI(IZ)*RAD)
        IF(ABS(CS1).GT.1.0) CS1=SIGN(1.0,CS1)
        ANG1=ACOS(CS1)/RAD
        DO 20 L=1,NLN
        L9=L
        INIT=1
CC INTERPOLATION OF ORIGINAL PHASE FUNCTION FROM GIVEN PHASE FUNCTION.
        IF(INDP.GE.0) THEN
          PHS(L)=PINT4(INIT,ANG1,KNANG,NANG,ANG,PHSF,L9)
CC INTERPOLATION OF ORIGINAL PHASE FUNCTION FROM -G-.
         ELSE
          SUM=0
          DO 21 K1=1,NLGN1(L)
   21     SUM=SUM+(2*K1-1)*G(K1,L)*PLGD(INIT,CS1)
          PHS(L)=SUM/4/PI
        ENDIF
CC INTERPOLATION OF TRUNCATED PHASE FUNCTION FROM -GT-.
        INIT=1
        SUM=0
        DO 22 K1=1,NLGT1(L)
   22   SUM=SUM+(2*K1-1)*GT(K1,L)*PLGD(INIT,CS1)
        PHST(L)=SUM/4/PI
   20   CONTINUE
        CALL INTCR1(IMTHD,AM0(IS),AM1U(IU),CS1,NTAU,UTAU,UTAUT
     &   ,NLN,THK,THKT,OMG,OMGT,PHS,PHST,FF,KLGN1,MXLGN1,NLGN1
     &   ,NLGT1,G,EPSP,NCHK1,COR,SGL2,PHSB,ERR)
        DO 23 IT=1,NTAU
   23   AI(IU,IS,IZ,IT)=AI(IU,IS,IZ,IT)+COR(IT)*FSOL
C ANGULAR LOOP END
   19 CONTINUE
C NORMAL END
      ERR=' '
      RETURN
      END
      SUBROUTINE CHKRT(INDA,INDT,INDP,IMTHD,NDA,NA1U,AM1U,NA0,AM0
     &,NFI,FI,NLN,THK,OMG,NLGN1,G,NANG,ANG,PHSF,EPSP,EPSU,GALB,FSOL
     &,NPLK1,CPLK,BGND,NTAU,UTAU,AMUA,WA,DPT,MXLGN2,ERR)
C CHECK-IN VARIABLES FOR -RTRN1- AND SET SOME VARIABLES.
C--- HISTORY
C 90. 1.20   CREATED
C    11.22   ADD: ELSE MXLGN2=2*NDA+1
C    12. 1   ADD LOOP-19.
C INPUT/OUTPUT SEE MAIN ROUTINE.
      PARAMETER (KNA1U =50)
      PARAMETER (KNA0  =1)
      PARAMETER (KNDM  =16)
      PARAMETER (KNFI  =50)
      PARAMETER (KLGN1 =400)
      PARAMETER (KNLN  =4)
      PARAMETER (KNTAU =1)
      PARAMETER (KPLK1 =2)
      PARAMETER (KNANG =60)
C
      PARAMETER (KNLN1=KNLN+1,KLGT1=2*KNDM)
      PARAMETER (PI=3.141592654)
C AREAS FOR THIS ROUTINE
      CHARACTER ERR*(*)
      DIMENSION AM1U(KNA1U),AM0(KNA0),FI(KNFI),THK(KNLN),OMG(KNLN)
     &,NLGN1(KNLN),G(KLGN1,KNLN),ANG(KNANG),PHSF(KNANG,KNLN)
     &,CPLK(KPLK1,KNLN),UTAU(KNTAU)
      DIMENSION AMUA(KNDM),WA(KNDM),DPT(KNLN1)
      DIMENSION CCP(3)
      CALL CPCON(CCP)
      EPS=CCP(1)*10
      ERR=' '
C EPSP, EPSU, FSOL
      IF(EPSP.LT.0) THEN
        ERR='ILLEGAL VALUE OF EPSP'
        RETURN
      ENDIF
      IF(EPSU.LT.0) THEN
        ERR='ILLEGAL VALUE OF EPSU'
        RETURN
      ENDIF
      IF(FSOL.LT.0) THEN
        ERR='ILLEGAL VALUE OF FSOL'
        RETURN
      ENDIF
C INDA
      IF(INDA.LT.0 .OR. INDA.GT.2) THEN
        ERR='ILLEGAL VALUE OF INDA'
        RETURN
      ENDIF
C INDT
      IF(INDT.LT.0 .OR. INDT.GT.2) THEN
        ERR='ILLEGAL VALUE OF INDT'
        RETURN
      ENDIF
C INDP
      IF(INDP.LT.-1 .OR. INDP.GT.1) THEN
        ERR='ILLEGAL VALUE OF INDP'
        RETURN
       ENDIF
C IMTHD
      IF(IMTHD.GT.3) THEN
        ERR='ILLEGAL VALUE OF IMTHD'
        RETURN
      ENDIF
C NDA
      IF(NDA.LE.0 .OR. NDA.GT.KNDM) THEN
        ERR='ILLEGAL VALUE OF NDA'
        RETURN
      ENDIF
C SET QUADRATURE
      CALL QUADA(NDA,AMUA,WA)
C NA0, AM0
      IF(NA0.LE.0 .OR. NA0.GT.KNA0) THEN
        ERR='ILLEGAL VALUE OF NA0'
        RETURN
      ENDIF
      DO 20 I=1,NA0
      IF(AM0(I).LE.0.0 .OR. AM0(I).GT.1.0) THEN
        ERR='ILLEGAL VALUE OF AM0'
        RETURN
      ENDIF
   20 CONTINUE
C INDA, AM1U, FI
      IF(INDA.GT.0) THEN
        IF(NFI.LE.0 .OR. NFI.GT.KNFI) THEN
          ERR='ILLEGAL VALUE OF NFI'
          RETURN
        ENDIF
        IF(INDA.EQ.2) NA1U=2*NDA
        IF(NA1U.LE.0 .OR. NA1U.GT.KNA1U) THEN
          ERR='ILLEGAL VALUE OF NA1U'
          RETURN
        ENDIF
        IF(INDA.EQ.2) THEN
          DO 1 I=1,NDA
          AM1U(I)        = -AMUA(I)
    1     AM1U(NA1U+1-I) =  AMUA(I)
        ENDIF
      ENDIF
C NLN
      IF(NLN.LE.0 .OR. NLN.GT.KNLN) THEN
        ERR='ILLEGAL VALUE OF NLN'
        RETURN
      ENDIF
      DPT(1)=0
      DO 17 L=1,NLN
   17 DPT(L+1)=DPT(L)+THK(L)
C NTAU
      IF(INDT.EQ.2) NTAU=NLN+1
      IF(NTAU.LE.0 .OR. NTAU.GT.KNTAU) THEN
        ERR='ILLEGAL VALUE OF NTAU'
        RETURN
      ENDIF
C UTAU
      IF(INDT.EQ.2) THEN
        DO 3 IT=1,NTAU
    3   UTAU(IT)=DPT(IT)
      ENDIF
C NLGN1,  MXLGN2=MAX(NLGN1)
      IF(INDP.LE.0 .OR. (INDA.GT.0 .AND. IMTHD.EQ.3)) THEN
        MXLGN2=1
        DO 4 L=1,NLN
        N1=NLGN1(L)
        IF(N1.LE.0 .OR. N1.GT.KLGN1) THEN
          ERR='ILLEGAL VALUE OF NLGN1'
          RETURN
        ENDIF
        MXLGN2=MAX(MXLGN2,N1)
    4   CONTINUE
       ELSE
        MXLGN2=2*NDA+1
      ENDIF
C NANG
      IF(INDP.GT.0 .OR. (IMTHD.GE.1 .AND. INDA.GE.1)) THEN
        IF(NANG.LE.3 .OR. NANG.GT.KNANG) THEN
          ERR='YOU SHOULD SET AT LEAST FOUR ANGLES FOR THIS CONDITION'
          RETURN
        ENDIF
        DO 12 I=2,NANG
          IF(ANG(I).LE.ANG(I-1)) THEN
            ERR='YOUR SHOULD SET SCATTERING ANGLE FROM 0 TO 180 DEGREES'
            RETURN
          ENDIF
   12   CONTINUE
      ENDIF
C CHECK ORDER OF DPT AND UTAU.
      DO 5 L=1,NLN
      IF(DPT(L+1).LT.DPT(L)) THEN
        ERR='DPT SHOULD BE SET FROM TOP TO BOTTOM'
        RETURN
      ENDIF
    5 CONTINUE
      IF(NTAU.GE.2) THEN
        DO 6 IT=2,NTAU
        IF(UTAU(IT).LT.UTAU(IT-1)) THEN
          ERR='UTAU SHOULD BE SET FROM TOP TO BOTTOM'
          RETURN
        ENDIF
    6   CONTINUE
      ENDIF
      DO 19 IT=1,NTAU
      IF(UTAU(IT).GT.DPT(NLN+1)) THEN
        IF(ABS(UTAU(IT)-DPT(NLN+1)).LE.EPS) THEN
          UTAU(IT)=DPT(NLN+1)
         ELSE
          ERR='UTAU IS OUT OF BOUNDS'
          RETURN
        ENDIF
      ENDIF
   19 CONTINUE
C RESET CPLK AND BGND, GALB
      IF(GALB.LT.0.0 .OR. GALB.GT.1.0) THEN
        ERR='ILLEGAL VALUE OF GALB'
        RETURN
      ENDIF
      IF(NPLK1.GT.KPLK1) THEN
        ERR='TOO LARGE -NPLK1-, CHANGE -KPLK1-'
        RETURN
       ELSE
        IF(NPLK1.LE.0) THEN
          BGND=0
          DO 18 L=1,NLN
          DO 18 K1=1,NPLK1
   18     CPLK(K1,L)=0
        ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE CONVU(N,AM,JJ,ER,IC,EPSU,NCHK)
C CHECK CONVERGENCE OF INTENSITY IN THE DIRECTION OF MU
C--- HISTORY
C 90. 1.20 CREATED
C--- INPUT
C N      I        NUMBER OF STREAMS TO BE CHECKED
C AM    R(N)      DIRECTION OF STREAM
C JJ    I(N)      STREAM NUMBER IN THE ORIGINAL ORDER
C ER    R(N)      MAXIMUM ERROR FOR THE STREAM REGARDLESS OF OTHER
C                  ANGLES AND LAYERS
C IC    I(N)      NUMBER OF CONSECTIVE SERIES WITH -ER- LESS THAN -EPSU-
C EPSU    R       CONVERGENCE CRITERION
C NCHK    I       MAX. NUMBER FOR THE VALUE OF -IC- BY WHICH THE
C                  ROUTINE CONFIRMS CONVERGENCE
C--- OUTPUT
C N, AM, JJ, IC   UPDATED AFTER DROPPING CONVERGENT STREAMS
C$ENDI
      DIMENSION AM(*),JJ(*),ER(*),IC(*)
      J1=0
      DO 1 J=1,N
      IF(ABS(ABS(AM(J))-1).LE.0) THEN
        IC(J)=999
        ER(J)=0
      ENDIF
      IF(ER(J).LT.EPSU) THEN
        IC(J)=IC(J)+1
        IF(IC(J).GE.NCHK) GOTO 1
       ELSE
        IC(J)=0
      ENDIF
      J1=J1+1
      IC(J1)=IC(J)
      JJ(J1)=JJ(J)
      AM(J1)=AM(J)
      ER(J1)=ER(J)
    1 CONTINUE
      N=J1
      RETURN
      END
      SUBROUTINE INTCR1(IMTHD,AM0,AM1U,CS1,NTAU,UTAU,UTAUT
     &,NLN,THK,THKT,OMG,OMGT,PHS,PHST,FF,KLGN1,MXLGN1,NLGN1
     &,NLGT1,G,EPSP,NCHK1,COR,SGL2,PHSB,ERR)
C SIGLE AND SECOND SCATTERING CORRECTION
C ASSUME DELTA-M METHOD
C--- HISTORY
C 90. 1.28 CREATED
C--- INPUT
C IMTHD      I         1: MS, 2: TMS, 3: IMS
C AM0        R         COS(SOLAR ZENITH ANGLE)
C AM1U       R         COS(EMERGENT NADIR ANGLE)
C                     .LT.O: UPWARD,  .GT.0: DOWNWARD
C CS1        R         COS(SCATTERING ANGLE)
C NTAU       I         NUMBER OF USER DEFINED ANGLE
C UTAU     R(NTAU)     USER DEFINED OPTICAL DEPTH
C UTAUT    R(NTAU)     TRUNCATED USER DEFINED DEPTH
C NLN        I         NUMBER OF LAYER
C THK      R(NLN)      OPTICAL THICKNESS OF SUBLAYER
C THKT     R(NLN)      TRUNCATED THICKNESS
C OMG      R(NLN)      SINGLE SCATTERING ALEBEDO
C OMGT     R(NLN)      TRUNCATED SINGLE SCATTERING ALBEDO
C PHS      R(NLN)      PHASE FUNCTION  (INTEGRAL = 1 OVER UNIT SPHERE)
C PHST     R(NLN)      TRUNCATED PHASE FUNCTION (1 OVER UNIT SPHERE)
C FF       R(NLN)      TRUNCATION FRACTION
C KLGN1      I         FIRST ARGUMENT SIZE OF G
C MXLGN1     I         MAX(NLGN1) (ONLY FOR IMTHD=3)
C NLGN1    I(NLN)      MAX ORDER OF LEGENDER SERIES +1 IN EACH SUBLAYER
C                        (ONLY FOR IMTHD=3)
C NLGT1    I(NLN)      SAME AS -NLGN1- BUT FOR TRUNCATION
C G      R(KLGN1,NLN)  PHASE FUNCTION MOMENT (ONLY FOR IMTHD=3)
C EPSP       R         CONVERGENCE CRITERION OF LEGENDRE SERIS OF
C                        PHASE FUNCTION ** 2
C NCHK1      I         NUMBER OF CONSECTIVE CONVERGENCES BEFORE
C                        FINAL DECISION
C--- OUTPUT
C COR      R(NTAU)     AI= AI+COR*FSOL IS CORRECTED INTENSITY
C                        AT EACH USER
C ERR      C*64        ERROR INDICATER
C                       DEFINED DEPTH (UTAU)
C--- WORK
C SGL2     R(NTAU)
C PHSB     R(NLN)
C--- AREA FOR THIS ROUTINE
      PARAMETER (PI=3.141592654, RAD=PI/180.0)
      CHARACTER ERR*(*)
      DIMENSION UTAU(NTAU),UTAUT(NTAU),THK(NLN),THKT(NLN)
     &,OMG(NLN),OMGT(NLN),PHS(NLN),PHST(NLN),FF(NLN),NLGN1(NLN)
     &,NLGT1(NLN),G(KLGN1,NLN),COR(NTAU),SGL2(NTAU),PHSB(NLN)
C
      ERR=' '
      DO 1 IT=1,NTAU
    1 COR(IT)=0
      IF(IMTHD.LE.0 .OR. IMTHD.GE.4) RETURN
C--- MS-METHOD IN REFERENCE-NT
      IF(IMTHD.EQ.1) THEN
C SGL2 = U1WAVE OF EQ.(14)
        DO 2 L=1,NLN
    2   PHSB(L)=(1-FF(L))*PHST(L)
        DO 3 IT=1,NTAU
    3   SGL2(IT)=SGLR(AM1U,AM0,NLN,THK,OMG,PHSB,UTAU(IT))
C CORRECTION BY MS-METHOD.  SEE EQ.(14)
        DO 4 IT=1,NTAU
CC SGL1 = U1 OF EQ.(14)
        SGL1=SGLR(AM1U,AM0,NLN,THK,OMG,PHS,UTAU(IT))
    4   COR(IT)=-SGL2(IT)+SGL1
        RETURN
      ENDIF
C--- TMS AND IMS METHODS
C  SGL2 = U1* OF EQ.(15)
      DO 5 IT=1,NTAU
    5 SGL2(IT)=SGLR(AM1U,AM0,NLN,THKT,OMGT,PHST,UTAUT(IT))
      DO 6 L=1,NLN
    6 PHSB(L)=PHS(L)/(1-FF(L))
C CORRECTION OF INTENSITY BY TMS-METHOD.  SEE EQ.(15)
      DO 7 IT=1,NTAU
CC  SGL1 = U1WAVE* OF EQ.(15)
      SGL1=SGLR(AM1U,AM0,NLN,THKT,OMGT,PHSB,UTAUT(IT))
    7 COR(IT)=-SGL2(IT)+SGL1
      IF(IMTHD.EQ.2 .OR. AM1U.LE.0.0) RETURN
C--- SECONDARY SCATTERING CORRECTION FOR IMS-METHOD
      DO 8 IT=1,NTAU
        UTAUS=UTAU(IT)
C GETTING MEAN OPTICAL CONSTANTS ABOVE THE USER DEFINED LEVEL-UTAUS
        EH=0
        SH=0
        SHH=0
        PHSPK=0
        DPT2=0
        DO 9 L=1,NLN
          DPT1=DPT2
          DPT2=DPT1+THK(L)
          IF(UTAUS.LE.DPT1) GOTO 10
          IF(UTAUS.LT.DPT2) THEN
            TAU=UTAUS-DPT1
           ELSE
            TAU=THK(L)
          ENDIF
          EH= EH     +       TAU
          SH= SH     +OMG(L)*TAU
          SHH=SHH    +OMG(L)*TAU*FF(L)
    9     PHSPK=PHSPK+OMG(L)*TAU*(PHS(L)-(1-FF(L))*PHST(L))
   10   IF(ABS(EH).LE.0) GOTO 8
C WH: MEAN SINGLE SCATTERING ALBEDO
        WH=SH/EH
C FH: MEAN TRUNCATED FRACTION
        FH=SHH/SH
        IF(FH.LE.0.0) GOTO 8
C PHSPK: MEAN TRUNCETED PEAK OF PHASE FUNCTION
        PHSPK=PHSPK/SHH
C AM3: VARIABLE APPEARING IN EQ.(23)
        AM3=AM0/(1-FH*WH)
C LEGENDRE SUM FOR EQ.(23)
C  WE TRUNCATE THE SERIES WHEN ABS(GPK**2) BECOMES SMALLER THAN -EPSP-
C   FOR SUCCESSIVE THREE TIMES.
        PHSPK2=0
        INIT=1
        ICHK=0
        DO 11 K1=1,MXLGN1
CC PHSPK: 2PHAT-PHAT2 OF EQ.(23)
CC  MEAN VALUE FOR THE LAYER ABOVE THE USER-DEFINED LEVEL
          GPK=0
          DPT2=0
          DO 12 L=1,NLN
          DPT1=DPT2
          DPT2=DPT1+THK(L)
          IF(K1.LE.NLGN1(L) .AND. OMG(L).GT.0) THEN
            IF(UTAUS.LE.DPT1) GOTO 13
            IF(UTAUS.LT.DPT2) THEN
              TAU=UTAUS-DPT1
             ELSE
              TAU=THK(L)
            ENDIF
            IF(K1.LE.NLGT1(L)) THEN
              GP=FF(L)
             ELSE
              GP=G(K1,L)
            ENDIF
            GPK=GPK+GP*OMG(L)*TAU
          ENDIF
   12     CONTINUE
   13     GPK=GPK/SHH
          GPK2=GPK**2
          IF(GPK2.LE.EPSP) THEN
            ICHK=ICHK+1
            IF(ICHK.GE.NCHK1) GOTO 14
          ENDIF
   11     PHSPK2=PHSPK2+(2*K1-1)*GPK2*PLGD(INIT,CS1)
   14   PHSPK2=PHSPK2/4/PI
CC CORRECTION BY UU = UU + UHAT IN IMS-METHOD
CC   SGL3: UHAT OF EQ.(23)
        SGL3=(FH*WH)**2/(1-FH*WH)*(2*PHSPK-PHSPK2)
     &    * HF(UTAUS, AM1U, AM3, AM3)
        COR(IT)=COR(IT)-SGL3
    8 CONTINUE
      RETURN
      END
      FUNCTION PINT4(INIT,ANG1,KNA,NA,ANG,P,L)
C INTERPOLATION OF THE PHASE FUNCTION.
C--- HISTORY
C 89.11. 8 CREATED FROM PINT3. CHANGE INTERPOLATION X-> ANG
C 90. 1.23 DEBUG
C--- INPUT
C INIT       I      1 THEN SEARCH ANG1-INTERVAL ELSE NOT SEARCH.
C ANG1       R      SCATTERING ANGLE IN DEGREE FOR INTERPOLATION.
C NA         I      NO. OF SCATTERING ANGLES.
C ANG      R(NA)    SCATTERING ANGLES IN DEGREE.
C P     R(KNA,L)    PHASE FUNCTION
C L          I      LAYER NUMBER
C--- OUTPUT VARIABLES
C INIT       I      0
C PINT4      R      INTERPOLATED VALUE OF P AT ANG1.
C$ENDI
C--- VARIABLES FOR THE ROUTINE.
      SAVE I1,I2,I3
      DIMENSION ANG(NA),P(KNA,L)
C--- WORKING AREAS.
      PARAMETER (PI=3.141592654, RAD=PI/180.0)
C
      IF(INIT.GE.1) THEN
        INIT=0
        DO 1 I=1,NA-1
        IF((ANG1-ANG(I))*(ANG1-ANG(I+1)).LE.0.0) GOTO 2
    1   CONTINUE
        I=NA
    2   IF(I-1) 3,3,4
    3   I1=1
        I3=3
        GO TO 5
C *
    4   IF(I-NA+1) 6,7,7
    6   I1=I-1
        I3=I+1
        GO TO 5
C *
    7   I1=NA-2
        I3=NA
    5   I2=I1+1
      ENDIF
      XX=ANG1
      X1=ANG(I1)
      X2=ANG(I2)
      X3=ANG(I3)
      ALP1=P(I1,L)
      ALP2=P(I2,L)
      ALP3=P(I3,L)
      ISIGN=-1
      IF(ALP1.GT.0.0 .AND. ALP2.GT.0.0 .AND. ALP3.GT.0.0) THEN
        ISIGN=1
        ALP1=LOG(ALP1)
        ALP2=LOG(ALP2)
        ALP3=LOG(ALP3)
      ENDIF
      PP=(XX-X2)*(XX-X3)/(X1-X2)/(X1-X3)*ALP1
     &  +(XX-X1)*(XX-X3)/(X2-X1)/(X2-X3)*ALP2
     &  +(XX-X1)*(XX-X2)/(X3-X1)/(X3-X2)*ALP3
      IF(ISIGN.GE.1) PP=EXPFN(PP)
      PINT4=PP
      RETURN
      END
      FUNCTION HF(TAU,AM1,AM2,AM3)
C GEOMETRICAL FACTOR FOR THE SCONDERY SCATTERING EQ.(24) OF NT.
C    HF=INTEG(0,TAU)DT*INTEG(0,T)DT1*EXP(T*(1/MU1-1/MU2)
C             + T1*(1/MU2-1/MU3))*EXP(-TAU/AM1)/AM1/AM2
C--- REFERENCE
C NT:  T. NAKAJIMA AND M. TANAKA, 1988, JQSRT, 40, 51-69
C--- HISTORY
C 88. 9.22  CREATED BY T. NAKAJIMA
C 89. 5. 4  USE EXPFN
C--- INPUT
C TAU      R         OPTICAL THICKNESS OF THE LAYER.
C AM1      R         COS(ZENITH ANGLE-1).
C AM2      R         COS(ZENITH ANGLE-2).
C AM3      R         COS(ZENITH ANGLE-3).
C--- OUTPUT
C HF       F         GEOMETRICAL FACTOR.
C
      SAVE INIT,EPS
      DIMENSION CCP(3)
      DATA INIT/1/
C
C SET EPS: IF ABS(1/AM1 - 1/AM0)*TAU .LE. EPS THEN
C                      THE ROUTINE SETS ALMUCANTAR CONDITION-IALM.
      IF(INIT.GT.0) THEN
        INIT=0
        CALL CPCON(CCP)
        EPS=CCP(1)*30
      ENDIF
C
      X1=1/AM1-1/AM2
      X2=1/AM2-1/AM3
      X3=1/AM1-1/AM3
      EX1=-TAU/AM1
      EX2=-TAU/AM2
      EX3=-TAU/AM3
      EX1=EXPFN(EX1)
      EX2=EXPFN(EX2)
      EX3=EXPFN(EX3)
C
      IF(ABS(X2*TAU).LE.EPS) GOTO 1
C X2 <> 0
CC I1
      IF(ABS(X3*TAU).LE.EPS) THEN
        AI1=EX1*(TAU+X3*TAU*TAU/2)
       ELSE
        AI1=(EX3-EX1)/X3
      ENDIF
CC I2
      IF(ABS(X1*TAU).LE.EPS) THEN
        AI2=EX1*(TAU+X1*TAU*TAU/2)
       ELSE
        AI2=(EX2-EX1)/X1
      ENDIF
      HF=(AI1-AI2)/AM1/AM2/X2
      RETURN
C X2 =  0
    1 IF(ABS(X1*TAU).LE.EPS) THEN
        HF=TAU**2*(0.5-X1*TAU/3)*EX1/AM1/AM2
       ELSE
        HF=((TAU-1/X1)*EX2+EX1/X1)/AM1/AM2/X1
      ENDIF
      RETURN
      END
      SUBROUTINE CORRR(TH0,AI,NA1U,NFAI,NRDRD,RDRD,ES)
C
C  Integration of radiance within IFOV.
C
C --- history
C 98. 5.28 Created
C     6.28 debugged
      PARAMETER (KNA1U =50)
      PARAMETER (KNA0  =1)
      PARAMETER (KNFI  =50)
      PARAMETER (PI=3.141592653,RAD=PI/180.0)
      dimension AI(KNA1U,KNA0,KNFI,*)
      DIMENSION X(10),Y(10),A(10),B(10),C(10),D(10)
C
      do 1 I=1,NRDRD
      X(I)=RDRD*(I-0.5)/REAL(RDRD)
      IF(TH0.LT.RDRD/2) then
        Y(I)=AI(NA1U-NRDRD+I,1,NFAI,1)
       else
        Y(I)=AI(NA1U,1,NFAI-NRDRD+I,1)
      endif
    1 continue
      CALL CSPL1(NRDRD,X,Y,A,B,C,D)
      N=4*NRDRD
      DR=RDRD/N
      ES=0
      do 2 I=1,N
      X1=DR*(I-0.5)
      ES1=CSPLI(X1,NRDRD,X,A,B,C,D)*2*PI*DR*X1*RAD**2
    2 ES=ES+ES1
      return
      end
      SUBROUTINE  QGAUSN( GWT, GMU, M )
C  COMPUTE WEIGHTS AND ABSCISSAE FOR ORDINARY GAUSSIAN QUADRATURE
C   (NO WEIGHT FUNCTION INSIDE INTEGRAL) ON THE INTERVAL (0,1)
C--- HISTORY
C 90. 1.17  REGISTERED
C--- INPUT
C M        I       ORDER OF QUADRATURE RULE
C--- OUTPUT
C GMU    R(M)      ARRAY OF ABSCISSAE (0, 1)
C GWT    R(M)      ARRAY OF WEIGHTS   SUM=1
C--- NOTES
C REFERENCE:  DAVIS,P.J. AND P. RABINOWITZ, METHODS OF NUMERICAL
C             INTEGRATION,ACADEMIC PRESS, NEW YORK, 1975, PP. 87
C METHOD:     COMPUTE THE ABSCISSAE AS ROOTS OF THE LEGENDRE
C             POLYNOMIAL P-SUB-N USING A CUBICALLY CONVERGENT
C             REFINEMENT OF NEWTON'S METHOD.  COMPUTE THE
C             WEIGHTS FROM EQ. 2.7.3.8 OF DAVIS/RABINOWITZ.
C             ACCURACY:  AT LEAST 13 SIGNIFICANT DIGITS.
C--- INTERNAL VARIABLES
C PM2,PM1,P : 3 SUCCESSIVE LEGENDRE POLYNOMIALS
C PPR       : DERIVATIVE OF LEGENDRE POLYNOMIAL
C P2PRI     : 2ND DERIVATIVE OF LEGENDRE POLYNOMIAL
C TOL       : CONVERGENCE CRITERION
C X,XI      : SUCCESSIVE ITERATES IN CUBICALLY-
C             CONVERGENT VERSION OF NEWTON'S METHOD
C            ( SEEKING ROOTS OF LEGENDRE POLYNOMIAL )
C$ENDI
      REAL     CONA, GMU( M ), GWT( M ), PI, T
      INTEGER  LIM, M, NP1
      DOUBLE   PRECISION  EN, NNP1, P, PM1, PM2, PPR, P2PRI, PROD,
     &                    TMP, TOL, X, XI
      DATA     TOL / 1.0D-13 /
      DATA     PI  / 3.1415926535898 /
C
      IF ( M.LE.1 )  THEN
         M = 1
         GMU( 1 ) = 0.5
         GWT( 1 ) = 1.0
         RETURN
      END IF
C
      EN   = M
      NP1  = M + 1
      NNP1 = M * NP1
      CONA = REAL( M-1 ) / ( 8 * M**3 )
C+---------------------------------------------------------------------+
C|         INITIAL GUESS FOR K-TH ROOT OF LEGENDRE POLYNOMIAL,         |
C|         FROM DAVIS/RABINOWITZ  EQ. (2.7.3.3A)                       |
C+---------------------------------------------------------------------+
      LIM  = M / 2
      DO 30  K = 1, LIM
         T = ( 4*K - 1 ) * PI / ( 4*M + 2 )
         X = COS ( T + CONA / TAN( T ) )
C
C+---------------------------------------------------------------------+
C|             RECURSION RELATION FOR LEGENDRE POLYNOMIALS             |
C|       INITIALIZE LEGENDRE POLYNOMIALS: (PM2) P-SUB-0, (PM1) P-SUB-1 |
C+---------------------------------------------------------------------+
 10       PM2 = 1.D0
         PM1 = X
         DO 20 NN = 2, M
            P   = ( ( 2*NN - 1 ) * X * PM1 - ( NN-1 ) * PM2 ) / NN
            PM2 = PM1
            PM1 = P
 20       CONTINUE
C
         TMP   = 1.D0 / ( 1.D0 - X**2 )
         PPR   = EN * ( PM2 - X * P ) * TMP
         P2PRI = ( 2.D0 * X * PPR - NNP1 * P ) * TMP
         XI    = X - ( P / PPR ) * ( 1.D0 +
     &               ( P / PPR ) * P2PRI / ( 2.D0 * PPR ) )
C
         IF ( DABS(XI-X) .GT. TOL ) THEN
C|          CHECK FOR CONVERGENCE
            X = XI
            GO TO 10
         END IF
C
C       ** ITERATION FINISHED--CALC. WEIGHTS, ABSCISSAE FOR (-1,1)
         GMU( K ) = - X
         GWT( K ) = 2.D0 / ( TMP * ( EN * PM2 )**2 )
         GMU( NP1 - K ) = - GMU( K )
         GWT( NP1 - K ) =   GWT( K )
 30    CONTINUE
C
      IF ( MOD( M,2 ) .NE. 0 )  THEN
C|       SET MIDDLE ABSCISSA AND WEIGHT FOR RULES OF ODD ORDER
         GMU( LIM + 1 ) = 0.0
         PROD = 1.D0
         DO 40 K = 3, M, 2
            PROD = PROD * K / ( K-1 )
 40       CONTINUE
         GWT( LIM + 1 ) = 2.D0 / PROD**2
      END IF
C
      DO 50  K = 1, M
C|       CONVERT FROM (-1,1) TO (0,1)
         GMU( K ) = 0.5 * GMU( K ) + 0.5
         GWT( K ) = 0.5 * GWT( K )
 50    CONTINUE
C
      RETURN
      END
      SUBROUTINE FTRN1(INIT,M,INDA,INDT,NA1U,AM1U,NA0,AM0
     &,NDA,AMUA,WA,NLN,NLGN1,G,DPT,OMG,NPLK1,CPLK,GALB,BGND
     &,FSOL,NTAU,UTAU,FLXD,FLXU,AI,ERR)
C
C SOLVE THE RADIATIVE TRANSFER IN ATMOSPHERE SYSTEM FOR EACH FOURIER
C COMPONENT.
C  -DOM AND ADDING METHOD-
C BY TERUYUKI NAKAJIMA
C SINCE THIS SYSTEM DOES NOT INCLUDE A OCEAN SURFACE,
C KNLT=KNLN IN TRN1.
C
C--- HISTORY
C 89.11. 2 CREATED FROM STRN7 WITH THERMAL EMISSION
C 90.12. 1 CALL AINTEG(M,L,  ->  CALL AINTEG(M,LT,
C 93. 5. 4 Split GRNDL
C 95. 5.25 change GRNDL by GRNDL3
C--- INPUT
C INIT       I       1: INITIALIZE THE PART (DEPENDENT AMUA, UTAU)
C                    0:  BYPASS THE M-INDEPENDENT PART.
C M          I       FORIER ORDER.
C INDA       I       0: FLUX ONLY.
C                    1: INTENSITY USING AM1U.
C                    2: NA1U AND AM1U ARE SUPPOSED TO BE
C                       2*NDA AND (-AMUA, +AMUA).
C INDT       I       0: SET USER DEFINED DEPTH.
C                    1: SAME AS ABOVE.
C                    2: SET NTAU AND UTAU AS NLN1 AND DPT.
C NA1U       I       NUMBER OF EMERGENT ZENITH ANGLES IN THE SPHERE.
C                      NA1U=2*NDA WHEN INDA=2.
C AM1U    R(KNA1U)   CONSINE OF THE EMERGENT ZENITH ANGLES.
C                      + FOR DOWNWARD, - FOR UPWARD.
C                      AM1U = (-AMUA, +AMUA) WHEN INDA=2.
C NA0        I       NUMBER OF THE SOLAR INCIDENCES.
C AM0     R(NA0)     CONSINE OF SOLAR ZENITH ANGLES .GT.0.
C NDA        I       NO. OF ZENITH-QUADRATURE ANGLES IN THE HEMISPHERE.
C AMUA    R(KNDM)    COSINE OF THE ZENITH-QUADRATURE ANGLES. 1 TO 0.
C WA      R(KNDM)    CORRESPONDING QUADRATURE WEIGHTS.
C NLN        I       NUMBER OF ATMOSPHERIC SUBLAYERS.
C NLGN1   I(KNLN)    MAXIMUM ORDER OF THE LEGENDRE SESIES OF THE PHASE
C                      FUNCTION + 1
C G     R(KLGN1,     LEGENDRE MOMENTS OF PHASE FUNCTION.
C           KNLN)      G(1,L)=1
C DPT     R(KNLN1)   OPTICAL DEPTH AT THE INTERFACES BETWEEN SUBLAYERS.
C                      TOP TO BOTTOM (TOP = 0 FOR NORMAL APPLICATION).
C OMG     R(KNLN)    SINGLE SCATTERING ALBEDO.
C NPLK1      I       NUMBER OF ORDER TO APPROXIMATE PLANK + 1.
C                      IF 0 THEN NO THERMAL.
C CPLK    R(KPLK1    PLANK FUNCTION (B) =
C         ,KNLN)       SUM(IB=1,NPLK1) CPLK(IB,L) TAU**(IB-1).
C                      TAU IS OPTICAL DEPTH MEASURED FROM
C                      THE TOP OF THE SUBSURFACE.
C GALB       R       GROUND ALBEDO (LAMBERTIAN).
C BGND       R       THERMAL EMISSION FROM THE GROUND
C                      BGND=0 WHEN NPLK1=0.
C FSOL       R       SOLAR IRRADIANCE (W/M2/MICRON).
C NTAU       I       NUMBER OF USER DEFINED OPTICAL DEPTHS.
C                      NTAU=NLN+1 WHEN INDT=2.
C UTAU    R(KNTAU)   OPTICAL DEPTHS WHERE THE FIELD IS CALCULTED.
C                      TOP TO BOTTOM.
C                      UTAU=DPT WHEN INDT=2.
C--- OUTPUT
C INIT               0
C FLXD  R(KNA0,KNTAU) DOWNWARD FLUX AT UTAU.
C FLXU               SAME AS FLXD BUT FOR UPWARD FLUX.
C AI      R(KNA1U,   I(MU1U(I), MU0(J), L)
C           KNA0,    INTENSITY AT UTAU.
C           KNTAU)
C ERR      C*64      ERROR INDICATER.
C--- CORRESPONDENCE BETWEEN VARIABLES AND PARAMETERS
C KNA1U       NA1U
C KNA1        NA1
C KNA0        NA0
C KNDM        NDA
C KNLN        NLN
C KNLN1       KNLN+1
C KNTAU       NTAU
C KLGN1       NLGN1
C KPLK1       NPLK1
C--- NOTES FOR LOCAL VARIABLES
C NA1        I       NUMBER OF STREAMS FOR AM1.
C AM1     R(KNA1)    ABS(AM1U).
C IIAM1   I(KNA1U)   DIRECTION NUMBER OF AM1 FOR EACH AM1U.
C                      + FOR DOWNWARD, - FOR UPWARD.
C IITAU   I(KNTAU)   SUBLAYER NUMBER FOR USER DEFINED DEPTHS.
C$ENDI
      SAVE EPS,WMP,WMM,IITAU
C PARAMETERS
      PARAMETER (KNA1U =50)
      PARAMETER (KNA1  =50)
      PARAMETER (KNA0  =1)
      PARAMETER (KNDM  =16)
      PARAMETER (KNLN  =4)
      PARAMETER (KNTAU =1)
      PARAMETER (KPLK1 =2)
C
      PARAMETER (KNLN1=KNLN+1,KLGN1=2*KNDM)
      PARAMETER (PI=3.141592654,RAD=PI/180.0)
C AREAS FOR THIS ROUTINE
      CHARACTER ERR*(*)
      DIMENSION AM1U(KNA1U),AM0(KNA0),AMUA(KNDM),WA(KNDM)
     & ,NLGN1(KNLN),G(KLGN1,KNLN),DPT(KNLN1),OMG(KNLN)
     &,CPLK(KPLK1,KNLN),UTAU(KNTAU),FLXD(KNA0,KNTAU),FLXU(KNA0,KNTAU)
     &,AI(KNA1U,KNA0,KNTAU)
C WORKING AREAS
      DIMENSION AM1(KNA1),WMP(KNDM),WMM(KNDM)
     &,PL1(KNA1,KLGN1),PL0(KNA0,KLGN1),PLA(KNDM,KLGN1),GBUF(KLGN1)
     &,PT(KNDM,KNDM),PR(KNDM,KNDM),PT0(KNDM,KNA0),PR0(KNDM,KNA0)
     &,PR1(KNA1,KNDM),PT1(KNA1,KNDM),PR10(KNA1,KNA0),PT10(KNA1,KNA0)
     &,ALFA(KNDM,KNA0),BETA(KNDM,KNA0),UUP(KNDM,KNA0),UDN(KNDM,KNA0)
     &,AII(KNA1U,KNA0,KNLN1),AIB(KNA0),IIAM1(KNA1U),IITAU(KNTAU),CCP(3)
C FOR HOMOG2
      DIMENSION CPL(KPLK1),R(KNDM,KNDM),T(KNDM,KNDM)
     &,ER(KNDM,KNA0),ET(KNDM,KNA0),ZEIG(KNDM)
     &,Q(KNDM,KNDM),QI(KNDM,KNDM),C11(KNDM,KNDM),C22(KNDM,KNDM)
     &,VP(KNDM,KNA0),VM(KNDM,KNA0),DP(KNDM,KPLK1),DM(KNDM,KPLK1)
C FOR TRN1 (INCLUDE GROUND AS A SUBLAYER)
      PARAMETER (KNLNM=KNLN1,KNLNM1=KNLNM+1,KNLTM=KNLN1)
      DIMENSION IUP(KNLNM),IDN(KNLNM),NDD(KNLNM1)
     &, RE(KNDM,KNDM,KNLTM),  TE(KNDM,KNDM,KNLTM)
     &,SER(KNDM,KNA0,KNLNM),SET(KNDM,KNA0,KNLNM)
     &,RUP(KNDM,KNA0,KNLNM1),RDN(KNDM,KNA0,KNLNM1)
C FOR INTENSITY INTERPOLATION IN MULTI-LAYER SYSTEM
      DIMENSION QE(KNDM,KNDM,KNLN),QIE(KNDM,KNDM,KNLN)
     &,C1E(KNDM,KNDM,KNLN),C2E(KNDM,KNDM,KNLN),ZEE(KNDM,KNLN)
     &,DPE(KNDM,KPLK1,KNLN),DME(KNDM,KPLK1,KNLN)
     &,VPE(KNDM,KNA0,KNLN),VME(KNDM,KNA0,KNLN)
C SET AND CLEAR VARIABLES
      ERR=' '
      IF(INIT.GT.0) THEN
        INIT=0
        CALL CPCON(CCP)
        EPS=CCP(1)*10
        DO 1 I=1,NDA
        WMP(I)=SQRT(WA(I)*AMUA(I))
    1   WMM(I)=SQRT(WA(I)/AMUA(I))
C SET SUBLAYER WHERE USER-DEFINED DEPTH RESIDES -IITAU-.
        DO 40 IT=1,NTAU
        ODP=UTAU(IT)
        EPS1=ODP*EPS
        DO 41 L=1,NLN
          IF((ODP-DPT(L))*(ODP-DPT(L+1)).LE.0.0) GOTO 42
   41   CONTINUE
        ERR='OUT OF BOUNDS IN LAYER SETTING'
        RETURN
   42   IF(ABS(ODP-DPT(L)).LE.EPS1) THEN
          IITAU(IT)=-L
         ELSE
          IF(ABS(ODP-DPT(L+1)).LE.EPS1) THEN
            IITAU(IT)=-(L+1)
           ELSE
            IITAU(IT)=L
          ENDIF
        ENDIF
   40   CONTINUE
      ENDIF
C SET AM1, IIAM1
      IF(INDA.GT.0) THEN
        NA1=1
        AM1(1)=ABS(AM1U(1))
        IIAM1(1)= 1
        IF(NA1U.GE.2) THEN
          DO 43 IU=2,NA1U
          X=ABS(AM1U(IU))
          DO 44 J=1,NA1
          IF(ABS(X-AM1(J)).LE.EPS) THEN
            IIAM1(IU)= J
            GOTO 43
          ENDIF
   44     CONTINUE
          NA1=NA1+1
          IF(NA1.GT.KNA1) THEN
            ERR='SETTING ERROR OF AM1U'
            RETURN
          ENDIF
          AM1(NA1)=X
          IIAM1(IU)= NA1
   43     CONTINUE
        ENDIF
      ENDIF
C CHECK MAXIMUM NLGN1
      MXLGN1=0
      DO 2 L=1,NLN
    2 MXLGN1=MAX(MXLGN1,NLGN1(L))
      IF(MXLGN1.LE.0) THEN
        ERR='NLGN1 ARE ALL ZERO'
        RETURN
      ENDIF
      M1=M+1
C SET LEGENDRE POLYNOMIALS.
      CALL PLGND(M1,MXLGN1,NA0,KNA0,AM0 ,PL0)
      CALL PLGND(M1,MXLGN1,NDA,KNDM,AMUA,PLA)
C SCALING FOR SYMMETLICITY
        DO 5 I=1,NDA
        DO 5 J=1,MXLGN1
    5   PLA(I,J)=WMM(I)*PLA(I,J)
C SOLVE THE EIGENVALUE PROBLEM OF ATMOSPHERIC SUBLAYERS.
      DO 6 L=1,NLN
        LB=L
        NAL1=NLGN1(L)
        NDD(L)=NDA
        IUP(L)=L
        IDN(L)=L
        T1=DPT(L)
        T2=DPT(L+1)
        W=OMG(L)
        IF(NPLK1.GT.0) THEN
          DO 7 I=1,NPLK1
    7     CPL(I)=2*PI*(1-W)*CPLK(I,L)
        ENDIF
CC SCATTERING MEDIA
CC PT, PR
        CALL EQ12(GBUF,G,NAL1,LB,KLGN1)
        CALL PHAS2(M1,NAL1,NDA,NDA,KNDM,KNDM,KNDM,PT,PR
     &     ,GBUF,PLA,PLA)
CC PT0, PR0
        CALL PHAS2(M1,NAL1,NDA,NA0,KNDM,KNDM,KNA0,PT0,PR0
     &     ,GBUF,PLA,PL0)
CC EIGENVALUE PROBLEM
        CALL HOMOG2(M,T1,T2,W,NDA,AMUA,WMM,NA0,AM0
     &      ,PR,PT,PR0,PT0,FSOL,NPLK1,CPL,R,T,ER,ET,ZEIG
     &      ,Q,QI,C11,C22,VP,VM,DP,DM,ERR)
        IF(ERR.NE.' ') RETURN
        CALL EQ32(RE ,R  ,NDA,NDA,LB,KNDM,KNDM,KNLTM,KNDM,KNDM)
        CALL EQ32(TE ,T  ,NDA,NDA,LB,KNDM,KNDM,KNLTM,KNDM,KNDM)
        CALL EQ32(SER,ER ,NDA,NA0,LB,KNDM,KNA0,KNLNM,KNDM,KNA0)
        CALL EQ32(SET,ET ,NDA,NA0,LB,KNDM,KNA0,KNLNM,KNDM,KNA0)
C STORE EXTRA FOR INTERNAL FIELD AT ARBITRARY DEPTH.
        IF(INDT.NE.2 .OR. INDA.NE.2) THEN
          CALL EQ32(QE ,Q  ,NDA,NDA,LB,KNDM,KNDM,KNLN,KNDM,KNDM)
          CALL EQ32(QIE,QI ,NDA,NDA,LB,KNDM,KNDM,KNLN,KNDM,KNDM)
          CALL EQ32(C1E,C11,NDA,NDA,LB,KNDM,KNDM,KNLN,KNDM,KNDM)
          CALL EQ32(C2E,C22,NDA,NDA,LB,KNDM,KNDM,KNLN,KNDM,KNDM)
          CALL EQ32(VPE,VP ,NDA,NA0,LB,KNDM,KNA0,KNLN,KNDM,KNA0)
          CALL EQ32(VME,VM ,NDA,NA0,LB,KNDM,KNA0,KNLN,KNDM,KNA0)
          CALL EQ21(ZEE,ZEIG,NDA,LB,KNDM)
          IF(NPLK1.GT.0) THEN
            CALL EQ32(DPE,DP,NDA,NPLK1,LB,KNDM,KPLK1,KNLN,KNDM,KPLK1)
            CALL EQ32(DME,DM,NDA,NPLK1,LB,KNDM,KPLK1,KNLN,KNDM,KPLK1)
          ENDIF
        ENDIF
    6 CONTINUE
C LAMBERT SURFACE
      NLN1=NLN+1
      IF(GALB.LE.0.0 .AND. BGND.LE.0.0) THEN
        NLT=NLN
       ELSE
        NLT=NLN1
        NDD(NLT)=NDA
        IUP(NLT)=NLT
        IDN(NLT)=NLT
        T1=DPT(NLT)
        CALL GRNDL3(FSOL,GALB,BGND,T1,M,NDA,AMUA,WA,NA0,AM0,R,T,ER,ET)
        CALL EQ32(RE , R,NDA,NDA,NLT,KNDM,KNDM,KNLTM,KNDM,KNDM)
        CALL EQ32(TE , T,NDA,NDA,NLT,KNDM,KNDM,KNLTM,KNDM,KNDM)
        CALL EQ32(SER,ER,NDA,NA0,NLT,KNDM,KNA0,KNLNM,KNDM,KNA0)
        CALL EQ32(SET,ET,NDA,NA0,NLT,KNDM,KNA0,KNLNM,KNDM,KNA0)
      ENDIF
      NLT1=NLT+1
      NDD(NLT1)=NDA
C ADDING OF THE SUBLAYERS.
      CALL TRN1(NLT,NDD,NA0,IUP,IDN,RE,TE,SER,SET,RUP,RDN,ERR)
      IF(ERR.NE.' ') RETURN
      IF(INDT.EQ.2 .AND. INDA.NE.1) THEN
        IF(INDA.EQ.2) THEN
          DO 8 L=1,NLN1
          DO 8 J=1,NA0
          DO 9 I=1,NDA
    9     AI(I,J,L)=RUP(I,J,L)/WMP(I)
          DO 10 I=1,NDA
   10     AI(NA1U+1-I,J,L)=RDN(I,J,L)/WMP(I)
    8     CONTINUE
        ENDIF
        IF(M.EQ.0) THEN
          DO 11 L=1,NLN1
          DO 11 J=1,NA0
          FU1=0
          FD1=0
          EX1=-UTAU(L)/AM0(J)
          TRNS0=EXPFN(EX1)*FSOL
CC FOR SCALED INTENSITY
          DO 12 I=1,NDA
          FU1=FU1+WMP(I)*RUP(I,J,L)
   12     FD1=FD1+WMP(I)*RDN(I,J,L)
          FLXU(J,L)=FU1
          FLXD(J,L)=FD1+TRNS0*AM0(J)
   11     CONTINUE
        ENDIF
        RETURN
      ENDIF
C--- INTERPOLATION OF THE FIELD.
      IF(INDA.EQ.1) CALL PLGND(M1,MXLGN1,NA1,KNA1,AM1 ,PL1)
      DO 13 L=1,NLN
        LT=L
        LB=LT+1
        DPTH=DPT(LT)
        TAU=DPT(LB)-DPTH
        T1=DPT(LT)
        T2=DPT(LB)
        W=OMG(L)
C INTEGRAL CONSTANTS:  ALFA, BETA
        CALL CINGR(NDA,NA0,LT,AM0,TAU,RDN,RUP,VPE,VME,C1E,C2E
     &   ,NPLK1,DPE,DME,ALFA,BETA)
C FLUXES OR INTENSITY WHEN INDA=2
        IF(M.EQ.0 .OR. INDA.EQ.2) THEN
          DO 14 IT=1,NTAU
          L1=IABS(IITAU(IT))
          IF(L1.EQ.NLN1) L1=NLN
          IF(L1.NE.L) GOTO 14
          DTAU=UTAU(IT)-DPTH
          IF(IITAU(IT).GT.0) THEN
CC IN THE SUBLAYER
            CALL ADISC(M,L1,NDA,NA0,AM0,WMP,TAU,DTAU,ZEE,QE,QIE
     &      ,VPE,VME,ALFA,BETA,NPLK1,DPE,DME,UDN,UUP)
           ELSE
CC JUST INTERFACE
            LU=IABS(IITAU(IT))
            DO 15 J=1,NA0
            DO 15 I=1,NDA
            UDN(I,J)=RDN(I,J,LU)/WMP(I)
   15       UUP(I,J)=RUP(I,J,LU)/WMP(I)
          ENDIF
CC INTENSITY WHEN INDA=2
          IF(INDA.EQ.2) THEN
            DO 18 J=1,NA0
            DO 16 I=1,NDA
   16       AI(I,       J,IT)=UUP(I,J)
            DO 17 I=1,NDA
   17       AI(NA1U+1-I,J,IT)=UDN(I,J)
   18       CONTINUE
          ENDIF
CC FLUX
          IF(M.EQ.0) THEN
            DO 19 J=1,NA0
            FU1=0
            FD1=0
            EX1=-UTAU(IT)/AM0(J)
            TRNS0=EXPFN(EX1)*FSOL
            DO 20 I=1,NDA
            FU1=FU1+AMUA(I)*WA(I)*UUP(I,J)
   20       FD1=FD1+AMUA(I)*WA(I)*UDN(I,J)
            FLXU(J,IT)=FU1
   19       FLXD(J,IT)=FD1+TRNS0*AM0(J)
          ENDIF
   14     CONTINUE
        ENDIF
C--- INTENSITIES IN THE USER DEFINED DIRECTIONS (WITHOUT CONTRIBUTION
C     FROM THE INTERFACES).
C PHASE FUNCTIONS FOR ANGULAR INTERPOLATION
        IF(INDA.EQ.1) THEN
          NAL1=NLGN1(L)
          CALL EQ12(GBUF,G,NAL1,LT,KLGN1)
C PT10, PR10
          CALL PHAS2(M1,NAL1,NA1,NA0,KNA1,KNA1,KNA0,PT10,PR10
     &     ,GBUF,PL1,PL0)
C PT1, PR1
          CALL PHAS2(M1,NAL1,NA1,NDA,KNA1,KNA1,KNDM,PT1,PR1
     &     ,GBUF,PL1,PLA)
          DO 21 IU=1,NA1U
C AT INTERFACE.
          I=IIAM1(IU)
          IF(AM1U(IU).GT.0.0) THEN
CC DOWNWARD INCREMENT AT BOTTOM OF SUBLAYER (STAMNES INTEGRATION)
            CALL AINTEG(M,LT,NDA,NA0,I,FSOL,AM1U(IU),AM0,W,T2,T1,T2
     &       ,ZEE,QE,QIE,VPE,VME,PT10,PR10,PT1,PR1,ALFA,BETA
     &       ,NPLK1,DPE,DME,CPLK,AIB)
            DO 22 J=1,NA0
   22       AII(IU,J,LB)=AIB(J)
           ELSE
CC UPWARD INCREMENT AT TOP OF SUBLAYER (STAMNES INTEGRATION)
            CALL AINTEG(M,LT,NDA,NA0,I,FSOL,AM1U(IU),AM0,W,T1,T1,T2
     &       ,ZEE,QE,QIE,VPE,VME,PT10,PR10,PT1,PR1,ALFA,BETA
     &       ,NPLK1,DPE,DME,CPLK,AIB)
            DO 23 J=1,NA0
   23       AII(IU,J,L)=AIB(J)
          ENDIF
C AT USER DEFINED DEPTH
          DO 24 IT=1,NTAU
          L1=IITAU(IT)
CC  WE DO NOT CALCULATE THE FIELD AT THIS STAGE IF THE USER DEFINED
CC     DEPTH IS EXACTLY ON THE INTERFACES (L1.LT.0).
          IF(L1.EQ.L) THEN
           CALL AINTEG(M,LT,NDA,NA0,I,FSOL,AM1U(IU),AM0,W,UTAU(IT),T1,T2
     &      ,ZEE,QE,QIE,VPE,VME,PT10,PR10,PT1,PR1,ALFA,BETA
     &      ,NPLK1,DPE,DME,CPLK,AIB)
            DO 25 J=1,NA0
   25       AI(IU,J,IT)=AIB(J)
          ENDIF
   24     CONTINUE
   21     CONTINUE
        ENDIF
   13 CONTINUE
      IF(INDA.LE.0 .OR. INDA.EQ.2) RETURN
C--- ADDING INTENSITIES (INDA=1)
      DO 26 IU=1,NA1U
      I=IIAM1(IU)
      IF(AM1U(IU).GT.0.0) THEN
CC DOWNWARD
        DO 27 J=1,NA0
   27   AIB(J)=0
        DO 28 L=1,NLN1
        IF(L.GT.1) THEN
           EX1=-(DPT(L)-DPT(L-1))/AM1(I)
           EX1=EXPFN(EX1)
          DO 29 J=1,NA0
   29     AIB(J)=AIB(J)*EX1+AII(IU,J,L)
        ENDIF
        DO 30 IT=1,NTAU
        LU=IITAU(IT)
        IF(IABS(LU).EQ.L) THEN
          IF(LU.LT.0) THEN
            DO 31 J=1,NA0
   31       AI(IU,J,IT)=AIB(J)
           ELSE
            DTAU=UTAU(IT)-DPT(L)
            EX1=-DTAU/AM1(I)
            EX1=EXPFN(EX1)
            DO 32 J=1,NA0
   32       AI(IU,J,IT)=AI(IU,J,IT)+AIB(J)*EX1
          ENDIF
        ENDIF
   30   CONTINUE
   28   CONTINUE
CC UPWARD ADDING
       ELSE
        DO 33 J=1,NA0
        IF(NLN.LT.NLT) THEN
          AIB(J)=RUP(1,J,NLN1)/WMP(1)
         ELSE
          AIB(J)=0
        ENDIF
   33   CONTINUE
        DO 34 L=NLN1,1,-1
        IF(L.LE.NLN) THEN
           EX1=-(DPT(L+1)-DPT(L))/AM1(I)
           EX1=EXPFN(EX1)
          DO 35 J=1,NA0
   35     AIB(J)=AIB(J)*EX1+AII(IU,J,L)
        ENDIF
        DO 36 IT=1,NTAU
        LU=IITAU(IT)
        IF(LU.LT.0) THEN
          IF(IABS(LU).EQ.L) THEN
            DO 37 J=1,NA0
   37       AI(IU,J,IT)=AIB(J)
          ENDIF
         ELSE
          IF(LU.EQ.L-1) THEN
            DTAU=DPT(L)-UTAU(IT)
            EX1=-DTAU/AM1(I)
            EX1=EXPFN(EX1)
            DO 38 J=1,NA0
   38       AI(IU,J,IT)=AI(IU,J,IT)+AIB(J)*EX1
          ENDIF
        ENDIF
   36   CONTINUE
   34   CONTINUE
      ENDIF
   26 CONTINUE
      RETURN
      END
      SUBROUTINE HOMOG2(M,T1,T2,OMG,N1,AM1,WMM,N2,AM2
     &    ,PR,PT,PR0,PT0,FSOL,NPLK1,CPLK,R,T,ER,ET,ZEIG
     &    ,Q,QI,C11,C22,VP,VM,DP,DM,ERR)
C SOLVE THE TRANSFER IN A HOMOGENEOUS SCATTERING AND EMITTING
C  MEDIUM BY THE DISCRETE ORDINATE METHOD.
C--- HISTORY
C 89.10.31 CREATED FROM HOMOG1 INCLUDING THERMAL RADIATION.
C 95.11.22 EPS=CCP(1)*10 -> EPS=EXP(LOG(CCP(1))*0.8)
C--- INPUT
C T1         R       OPTICAL DEPTH AT THE LAYER TOP.
C T2         R       OPTICAL DEPTH AT THE LAYER BOTTOM.
C OMG        R       SINGLE SCATTERING ALBEDO.
C N1         I       NO. OF THE QUADRATURE STREAMS.
C AM1     R(KNDM)    MU (I), I=1,N1.
C WMM     R(KNDM)    SQRT(W1/M1)
C N2         I       NO. OF THE SOLAR DIRECTIONS.
C AM2     R(KNA0)    MU0(I), I=1,N2.
C PR      R(KNDM,    SCALED P+-(I,J) I,J=1,N1
C           KNDM)
C PT                 SCALED P++(I,J)
C PR0     R(KNDM,    SCALED P0+-(I,J)  I=1,N1; J=1,N2.
C PT0       KNDM)    SCALED P0++(I,J)  I=1,N1; J=1,N2.
C FSOL       R       SOLAR IRRADIANCE AT THE TOP OF THE SYSTEM.
C NPLK1      I       MAX ORDER OF PLANK FUNCTION EXPANSION BY TAU + 1.
C                      IF NPLK1=0 THEN NO THERMAL.
C CPLK    R(NPLK1)   2*PI*(1-W)*B(N)
C--- OUTPUT
C R       R(KNDM,    REFLECTION   MATRIX   RIJ, I,J=1,N1.
C T         KNDM)    TRANSMISSION MATRIX   TIJ.
C ER      R(KNDM,    UPGOING   SOURCE MATRIX EU(I,J), I=1,N1; J=1,N2.
C ET        KNA0)    DOWNGOING SOURCE MATRIX ED(I,J).
C ZEIG    R(KNDM)    ROOT OF THE EIGENVALUES OF Z.
C Q                  Q-MATRIX
C QI                 INVERSE OF Q
C C11                C1 = INVERSE OF 2A-
C C22                C2 = INVERSE OF 2B-
C VP                 VS+
C VM                 VS-
C DP   R(KNDM,KPLK1) THERMAL EMISSION INTENSITY EXAPNSION (DNWARD)
C DM   R(KNDM,KPLK1) THERMAL EMISSION INTENSITY EXAPNSION (UPWARD)
C ERR      C*64      ERROR INDEX
C--- PARAMETER
C KNDM       I       DECLARED SIZE FOR NDA
C KNA0       I       DECLARED SIZE FOR NA0
C KPLK1      I       DECLARED SIZE FOR NPLK1
      PARAMETER (KNDM  =16)
      PARAMETER (KNA0  =1)
      PARAMETER (KPLK1 =2)
      PARAMETER (PI=3.141592653)
C--- AREAS FOR THE ROUTINE.
      CHARACTER ERR*(*)
      DIMENSION AM1(KNDM),WMM(KNDM),AM2(KNA0)
     &,PR(KNDM,KNDM),PT(KNDM,KNDM),PR0(KNDM,KNA0),PT0(KNDM,KNA0)
     &,CPLK(KPLK1),ZEIG(KNDM),R(KNDM,KNDM),T(KNDM,KNDM)
     &,ER(KNDM,KNA0),ET(KNDM,KNA0),Q(KNDM,KNDM),QI(KNDM,KNDM)
     &,C11(KNDM,KNDM),C22(KNDM,KNDM),VP(KNDM,KNA0),VM(KNDM,KNA0)
     &,DP(KNDM,KPLK1),DM(KNDM,KPLK1)
C--- WORKING AREAS
      PARAMETER (KROWIJ=KNDM*(KNDM+1)/2,KNDM2=2*KNDM)
      DIMENSION CCP(3),X(KNDM,KNDM),Y(KNDM,KNDM),XI(KNDM,KNDM)
     &,SP(KNDM,KNA0),SM(KNDM,KNA0),G1(KNDM,KNA0),GAM(KNDM,KNA0)
     &,E0(KNA0),C(KNDM),SL(KNDM),SL1(KNDM)
     &,AP(KNDM,KNDM),AM(KNDM,KNDM),BP(KNDM,KNDM),BM(KNDM,KNDM)
     &,IW(KNDM2),DP0(KNDM),DP1(KNDM),DM0(KNDM),DM1(KNDM)
C PRECISION
      CALL CPCON(CCP)
C     EPS=CCP(1)*10
      EPS=EXP(LOG(CCP(1))*0.8)
C
      ERR=' '
      TAU=T2-T1
C X, Y MATRICES
      DO 1 I=1,N1
        DO 2 J=1,N1
        X(I,J)=-OMG*(PT(I,J)-PR(I,J))
    2   Y(I,J)=-OMG*(PT(I,J)+PR(I,J))
        X(I,I)=1.0/AM1(I)+X(I,I)
    1   Y(I,I)=1.0/AM1(I)+Y(I,I)
C DECOMPOSITION OF XY
      IF(M.EQ.0 .AND. 1.0-OMG.LE.EPS) THEN
        IW0=1
       ELSE
        IW0=0
      ENDIF
      IF(M.EQ.0 .AND. IW0.EQ.0 .AND. NPLK1.GT.0) THEN
        IPK=1
       ELSE
        IPK=0
      ENDIF
      CALL GETQM(IW0,N1,X,Y,ZEIG,Q,QI,XI,IMN,ERR)
      IF(ERR.NE.' ') RETURN
C THERMAL SOURCE
CC C
      IF(IPK.EQ.0) THEN
        IF(NPLK1.GT.0) THEN
          DO 25 J=1,NPLK1
          DO 25 I=1,N1
          DP(I,J)=0
   25     DM(I,J)=0
        ENDIF
       ELSE
        DO 3 I=1,N1
          SUM1=0
          DO 4 J=1,N1
    4     SUM1=SUM1+Q(J,I)*WMM(J)
          DO 5 J=NPLK1,1,-1
          IF(J+2.GT.NPLK1) THEN
            CPLK1=0
            ELSE
            CPLK1=DP(I,J+2)
          ENDIF
    5     DP(I,J)=((J+1)*J*CPLK1+SUM1*CPLK(J))/ZEIG(I)**2
    3   CONTINUE
CC D
        CALL AXB(DM,Q,DP,N1,N1,NPLK1,KNDM,KNDM,KNDM)
        DO 6 J=1,NPLK1
          DO 7 I=1,N1
          SUM=0
          IF(J+1.LE.NPLK1) THEN
            DO 8 K=1,N1
    8       SUM=SUM+QI(K,I)*DP(K,J+1)
          ENDIF
          DP(I,J)=DM(I,J)-J*SUM
    7     DM(I,J)=DM(I,J)+J*SUM
    6   CONTINUE
      ENDIF
C SIGMA+ - (FOR SINGLE SCATTERING)
      DO 9 I=1,N1
      DO 9 J=1,N2
      SP(I,J)=OMG*(PT0(I,J)+PR0(I,J))
    9 SM(I,J)=OMG*(PT0(I,J)-PR0(I,J))
C LOWER G
      DO 10 I=1,N1
      DO 10 J=1,N2
      SUM=0
      DO 11 K=1,N1
   11 SUM=SUM+X(I,K)*SP(K,J)
   10 G1(I,J)=-SUM-SM(I,J)/AM2(J)
C GAMMA
      CALL AXB(GAM,QI,G1,N1,N1,N2,KNDM,KNDM,KNDM)
      DO 12 I=1,N1
      DO 12 J=1,N2
   12 GAM(I,J)=GAM(I,J)/(1.0/AM2(J)**2-ZEIG(I)**2)
C VS+ AND -
      DO 14 J=1,N2
      TRNS0=-T1/AM2(J)
      TRNS0=EXPFN(TRNS0)*FSOL
      DO 14 I=1,N1
      SUM1=0.0
      SUM2=0.0
      DO 13 K=1,N1
      SUM1=SUM1+Q (I,K)*GAM(K,J)
   13 SUM2=SUM2+QI(K,I)*GAM(K,J)/AM2(J)+XI(I,K)*SM(K,J)
      VP(I,J)=(SUM1+SUM2)/2.0*TRNS0
   14 VM(I,J)=(SUM1-SUM2)/2.0*TRNS0
C E0
      DO 15 I=1,N2
      EX1=-TAU/AM2(I)
   15 E0(I)=EXPFN(EX1)
C BASE FUNCTION  C(TAU) AND S(TAU).
      DO 16 I=1,N1
   16 CALL CSFN(TAU,TAU,ZEIG(I),C(I),SL(I),SL1(I))
C A+-, B+-
      DO 17 I=1,N1
      DO 17 J=1,N1
      SUM1=Q (I,J)*C  (J)
      SUM2=QI(J,I)*SL (J)
      SUM3=Q (I,J)*SL1(J)
      SUM4=QI(J,I)*C  (J)
      AP(I,J)=SUM1-SUM2
      AM(I,J)=SUM1+SUM2
      BP(I,J)=SUM3-SUM4
   17 BM(I,J)=SUM3+SUM4
C C11 AND C22 -> THEIR INVERSION.
      DO 18 I=1,N1
      DO 18 J=1,N1
      C11(I,J)=2.0*AM(I,J)
   18 C22(I,J)=2.0*BM(I,J)
      CALL TNVSS2(N1,C11,DT,0.0,KNDM,IW,ERR)
      IF(ERR.NE.' ') THEN
        ERR='ERROR TO GET -C11- (HOMOG2)'
        RETURN
      ENDIF
      CALL TNVSS2(N1,C22,DT,0.0,KNDM,IW,ERR)
      IF(ERR.NE.' ') THEN
        ERR='ERROR TO GET -C22- (HOMOG2)'
        RETURN
      ENDIF
C R, T MATRICES
      DO 19 I=1,N1
      DO 19 J=1,N1
      SUM1=0
      SUM2=0
      DO 20 K=1,N1
      SUM1=SUM1+AP(I,K)*C11(K,J)
   20 SUM2=SUM2+BP(I,K)*C22(K,J)
      R(I,J)=SUM1+SUM2
   19 T(I,J)=SUM1-SUM2
C ER, ET MATRICES
      DO 22 I=1,N1
      DP1(I)=0
      DM1(I)=0
      IF(IPK.EQ.1) THEN
        DO 21 J=1,NPLK1
        DP1(I)=DP1(I)+DP(I,J)*TAU**(J-1)
   21   DM1(I)=DM1(I)+DM(I,J)*TAU**(J-1)
        DP0(I)=DP(I,1)
        DM0(I)=DM(I,1)
       ELSE
        DP0(I)=0
        DM0(I)=0
      ENDIF
   22 CONTINUE
      DO 23 J=1,N2
      DO 23 I=1,N1
      SUM1=0
      SUM2=0
      DO 24 K=1,N1
      VP0=VP(K,J)      +DP0(K)
      VM1=VM(K,J)*E0(J)+DM1(K)
      SUM1=SUM1+R(I,K)*VP0+T(I,K)*VM1
   24 SUM2=SUM2+T(I,K)*VP0+R(I,K)*VM1
      ER(I,J)=VM(I,J)      +DM0(I)-SUM1
   23 ET(I,J)=VP(I,J)*E0(J)+DP1(I)-SUM2
      RETURN
      END
      SUBROUTINE GETQM(IW0,N,X,Y,ZEIG,Q,QI,XI,IMN,ERR)
C SOLVE    XY = Q ZEIG**2 INVERSE(Q)
C ROOT DECOMPOSITION METHOD
C--- HISTORY
C 89. 8. 4 CREATED
C--- INPUT
C IW0       I        IF 1 THEN RENOMALIZATION (FOR M=0 AND W0=1)
C N         I        ORDER OF MATRICES
C X     R(KNDM,KNDM) SYMMETRIC MATRIX
C Y     R(KNDM,KNDM) SYMMETRIC MATRIX
C--- OUTPUT
C ZEIG    R(KNDM)    SQRT(EIGENVALUE)
C Q     R(KNDM,KNDM) ROTATION MATIX
C QI    R(KNDM,KNDM) INVERSE OF Q
C XI    R(KNDM,KNDM) INVERSE OF X
C IMN       I        LOCATION OF MINIMUM EIGENVALUE
C ERR     C*64       ERROR INDICATER
C--- PRPC-PARAMETER
C KNDM      I        DECLARED SIZE OF MATRICES
C
      PARAMETER (KNDM  =16)
      CHARACTER ERR*(*)
      DIMENSION X(KNDM,KNDM),Y(KNDM,KNDM),ZEIG(KNDM),Q(KNDM,KNDM)
     & ,QI(KNDM,KNDM),XI(KNDM,KNDM)
C WORKING AREA
      PARAMETER (KROWIJ=(KNDM*(KNDM+1))/2)
      DIMENSION V(KNDM,KNDM),SQX(KNDM,KNDM),SQXI(KNDM,KNDM)
     & ,ROWIJ(KROWIJ),WK(KNDM)
C ROOT DECOMPOSITION OF X (USE V AND ZEING FOR U AND XEIG).
      ERR=' '
      K=0
      DO 1 I=1,N
      DO 1 J=1,I
      K=K+1
    1 ROWIJ(K)=X(I,J)
      CALL SYMTRX(ROWIJ,N,ZEIG,V,KNDM,WK,IERR)
      IF(IERR.GT.128) THEN
        ERR='ERROR IN DECOMPOSITION OF X IN GETQM'
        RETURN
      ENDIF
      DO 2 I=1,N
      IF(ZEIG(I).LE.0.0) THEN
        ERR='NON-POSITIVE EIGENVALUE OF X'
        RETURN
      ENDIF
    2 ZEIG(I)=SQRT(ZEIG(I))
      DO 5 I=1,N
      DO 5 J=1,N
      SUM1=0
      SUM2=0
      SUM3=0
      DO 6 K=1,N
      SUM1=SUM1+V(I,K)*ZEIG(K)*   V(J,K)
      SUM2=SUM2+V(I,K)/ZEIG(K)*   V(J,K)
    6 SUM3=SUM3+V(I,K)/ZEIG(K)**2*V(J,K)
      SQX (I,J)=SUM1
      SQXI(I,J)=SUM2
    5 XI  (I,J)=SUM3
      CALL AXB(V,SQX,Y,N,N,N,KNDM,KNDM,KNDM)
      CALL AXB(Q,V,SQX,N,N,N,KNDM,KNDM,KNDM)
C ROOT DECOMPOSITION OF Z (USE Q FOR Z).
      K=0
      DO 7 I=1,N
      DO 7 J=1,I
      K=K+1
    7 ROWIJ(K)=Q(I,J)
      CALL SYMTRX(ROWIJ,N,ZEIG,V,KNDM,WK,IERR)
      IF(IERR.GT.128) THEN
        ERR='ERROR IN DECOMPOSITION OF Z IN GETQM'
        RETURN
      ENDIF
C CHECK MINIMUM EIGENVALUE
      IMN=1
      ZMN=ZEIG(1)
      IF(N.GE.2) THEN
        DO 8 J=2,N
        IF(ZEIG(J).LT.ZMN) THEN
          IMN=J
          ZMN=ZEIG(J)
        ENDIF
    8   CONTINUE
      ENDIF
C RENORMALIZATION
      IF(IW0.EQ.1) ZEIG(IMN)=0
C
      DO 9 I=1,N
      IF(ZEIG(I).LT.0.0) THEN
        ERR='NON-POSITIVE EIGENVALUE OF Z'
        RETURN
      ENDIF
    9 ZEIG(I)=SQRT(ZEIG(I))
C Q-MATRICES
      DO 10 I=1,N
      DO 10 J=1,N
      SUM1=0
      SUM2=0
      DO 11 K=1,N
      SUM1=SUM1+SQX(I,K)*V(K,J)
   11 SUM2=SUM2+V(K,I)*SQXI(K,J)
      Q (I,J)=SUM1
   10 QI(I,J)=SUM2
      RETURN
      END
      SUBROUTINE CINGR(NDA,NA0,L,AM0,TAU,RDN,RUP,VPE,VME,C1E,C2E
     &  ,NPLK1,DPE,DME,ALFA,BETA)
C GET ALFA AND BETA (INTEGRAL CONSTANTS).
C PARAMETERS
      PARAMETER (KNA0  =1)
      PARAMETER (KNDM  =16)
      PARAMETER (KNLN  =4)
      PARAMETER (KPLK1 =2)
C
      PARAMETER (KNLN1=KNLN+1)
      PARAMETER (KNLNM=KNLN1,KNLNM1=KNLNM+1)
      PARAMETER (PI=3.141592653,RAD=PI/180.0)
      DIMENSION AM0(KNA0),ALFA(KNDM,KNA0),BETA(KNDM,KNA0)
     &,RUP(KNDM,KNA0,KNLNM1),RDN(KNDM,KNA0,KNLNM1)
     &,DPE(KNDM,KPLK1,KNLN),DME(KNDM,KPLK1,KNLN)
     &,VPE(KNDM,KNA0,KNLN),VME(KNDM,KNA0,KNLN)
     &,C1E(KNDM,KNDM,KNLN),C2E(KNDM,KNDM,KNLN)
C WORK AREAS
      DIMENSION BUF1(KNDM),BUF2(KNDM)
C
      L1=L+1
      DO 1 J=1,NA0
        EX1=-TAU/AM0(J)
        EX1=EXPFN(EX1)
        DO 2 I=1,NDA
        SUM1=RDN(I,J,L)-VPE(I,J,L)
        SUM2=RUP(I,J,L1)-VME(I,J,L)*EX1
        IF(NPLK1.GT.0) THEN
          SUM1=SUM1-DPE(I,1,L)
          SUM=0
          DO 3 K=1,NPLK1
    3     SUM=SUM+DME(I,K,L)*TAU**(K-1)
          SUM2=SUM2-SUM
        ENDIF
        BUF1(I)=SUM2+SUM1
    2   BUF2(I)=SUM2-SUM1
        DO 4 I=1,NDA
        SUM1=0
        SUM2=0
        DO 5 K=1,NDA
        SUM1=SUM1+C1E(I,K,L)*BUF1(K)
    5   SUM2=SUM2+C2E(I,K,L)*BUF2(K)
        ALFA(I,J)=SUM1
    4   BETA(I,J)=SUM2
    1 CONTINUE
      RETURN
      END
      SUBROUTINE ADISC(M,L,NDA,NA0,AM0,WMP,TC,T1,ZEE,QE,QIE
     &  ,VPE,VME,ALFA,BETA,NPLK1,DPE,DME,UDN,UUP)
C INTENSITY AT A USER DEFINED DEPTH.
C--- HISTORY
C 87. 3. 9
C 89.11. 2 MODIFIED
C--- INPUT
C M        I          FOURIER ORDER
C L        I          LAYER NUMBER.
C NDA      I          NUMBER OF QUADRATURE POINTS.
C NA0      I          NUMBER OF SOLAR ZENITH ANGLES.
C AM0    R(KNA0)      COS(SOLAR ZENITH ANGLES).
C WMP    R(KNDM)      SQRT(W*M)
C TC       R          OPTICAL THICKNESS OF THE LAYER.
C T1       R          OPITCAL DEPTH OF INTERPOLATION MEASURED FROM TOP.
C ZEE    R(KNDM,KNLN)  EIGENVALUES
C QE    R(KNDM,KNDM,KNLN)
C QIE   R(KNDM,KNDM,KNLN)
C VPE   R(KNDM,KNA0,KNLN)
C VME   R(KNDM,KNA0,KNLN)
C DPE   R(KNDM,KPLK1,KNLN)
C DME   R(KNDM,KPLK1,KNLN)
C ALFA  R(KNDM,KNA0)  INTEGRAL CONSTANT ALFA.
C BETA  R(KNDM,KNA0)  INTEGRAL CONSTATN BETA.
C--- OUTPUT
C UUP   R(KNDM,KNA0)  UPWARD INTENSITY (UNSCALE)
C UDN   R(KNDM,KNA0)  DNWARD INTENSITY (UNSCALE)
C
      PARAMETER (KNLN  =4)
      PARAMETER (KNA0  =1)
      PARAMETER (KNDM  =16)
      PARAMETER (KPLK1 =2)
C--- AREAS FOR THE ROUTINE
      DIMENSION AM0(NA0),ZEE(KNDM,KNLN),WMP(KNDM)
     &,QE(KNDM,KNDM,KNLN),QIE(KNDM,KNDM,KNLN)
     &,VPE(KNDM,KNA0,KNLN),VME(KNDM,KNA0,KNLN)
     &,ALFA(KNDM,KNA0),BETA(KNDM,KNA0)
     &,DPE(KNDM,KPLK1,KNLN),DME(KNDM,KPLK1,KNLN)
     &,UUP(KNDM,KNA0),UDN(KNDM,KNA0)
C--- WORK AREAS
      DIMENSION C(KNDM),SL(KNDM),SL1(KNDM),E0(KNA0)
     & ,AP(KNDM),AM(KNDM),BP(KNDM),BM(KNDM)
C A+-, B+-
      DO 1 I=1,NDA
    1 CALL CSFN(TC,T1,ZEE(I,L),C(I),SL(I),SL1(I))
      DO 2 J=1,NA0
      EX1=-T1/AM0(J)
    2 E0(J)=EXPFN(EX1)
C U+ = S1D,  U- = S1U
      DO 3 I=1,NDA
        DO 4 J=1,NDA
        AP(J)=QE(I,J,L)*C  (J) - QIE(J,I,L)*SL(J)
        AM(J)=QE(I,J,L)*C  (J) + QIE(J,I,L)*SL(J)
        BP(J)=QE(I,J,L)*SL1(J) - QIE(J,I,L)*C (J)
    4   BM(J)=QE(I,J,L)*SL1(J) + QIE(J,I,L)*C (J)
        DP=0
        DM=0
        IF(M.EQ.0 .AND. NPLK1.GT.0) THEN
          TAUN=1
          DO 5 J=1,NPLK1
          DP=DP+DPE(I,J,L)*TAUN
          DM=DM+DME(I,J,L)*TAUN
    5     TAUN=TAUN*T1
        ENDIF
        DO 3 J=1,NA0
        SUM1=0
        SUM2=0
        DO 6 K=1,NDA
        SUM1=SUM1+AP(K)*ALFA(K,J)+BP(K)*BETA(K,J)
    6   SUM2=SUM2+AM(K)*ALFA(K,J)+BM(K)*BETA(K,J)
        UDN(I,J)=(SUM1+VPE(I,J,L)*E0(J)+DP)/WMP(I)
        UUP(I,J)=(SUM2+VME(I,J,L)*E0(J)+DM)/WMP(I)
    3 CONTINUE
      RETURN
      END
      SUBROUTINE AINTEG(M,L,NDA,NA0,I,FSOL,AM1U,AM0,W,T,T1,T2
     &     ,ZEE,QE,QIE,VPE,VME,PT10,PR10,PT1,PR1,ALFA,BETA
     &     ,NPLK1,DPE,DME,CPLK,AI)
C--- INPUT
C M       I         FOURIER ORDER
C L       I         LAYER NUMBER
C NDA
C NA0
C FSOL
C I       I         STREAM NUMBER FOR AM1(I)
C AM1U    R         USER DEFINED DIRECTION
C AM0   R(KNA0)
C W       R         SINGLE SCATTERING ALBEDO
C T       R         OPTICAL DEPTH FOR INTERPOLATION
C T1      R         OPTICAL DEPTH AT THE TOP OF THE SUBLAYER
C T2      R         OPTICAL DEPTH AT THE BOTTOM OF THE SUBLAYER.
      PARAMETER (KNDM  =16)
      PARAMETER (KNA1U =50)
      PARAMETER (KNA1  =50)
      PARAMETER (KNA0  =1)
      PARAMETER (KNLN  =4)
      PARAMETER (KPLK1 =2)
C
      PARAMETER (PI=3.141592654)
      DIMENSION AM0(KNA0),ZEE(KNDM,KNLN),QE(KNDM,KNDM,KNLN)
     & ,QIE(KNDM,KNDM,KNLN),VPE(KNDM,KNA0,KNLN),VME(KNDM,KNA0,KNLN)
     & ,PT10(KNA1,KNA0),PR10(KNA1,KNA0),PT1(KNA1,KNDM),PR1(KNA1,KNDM)
     & ,DPE(KNDM,KPLK1,KNLN),DME(KNDM,KPLK1,KNLN),AI(KNA0)
     & ,ALFA(KNDM,KNA0),BETA(KNDM,KNA0),CPLK(KPLK1,KNLN)
C WORKING AREAS
      DIMENSION H1(KNDM),H2(KNDM),C(KNDM),SL(KNDM),SL1(KNDM),EI(KNA0)
     & ,PKI(KPLK1)
      TT1=0
      TT2=T2-T1
      TT =T -T1
      AM1=ABS(AM1U)
      DO 1 J=1,NDA
        SUM1=0
        SUM2=0
        DO 2 K=1,NDA
        SUM1=SUM1+(PT1(I,K)+PR1(I,K))*QE (K,J,L)
    2   SUM2=SUM2+(PT1(I,K)-PR1(I,K))*QIE(J,K,L)
        H1(J)=W*SUM1
    1 H2(J)=W*SUM2
C DOWNWARD INTENSITY
      IF(AM1U.GT.0.0) THEN
        DO 3 J=1,NA0
    3   CALL EXINT(AM1U,AM0(J),TT,TT1,TT,EI(J))
        DO 4 J=1,NDA
    4   CALL CSINT(AM1U,ZEE(J,L),TT2,TT,TT1,TT,C(J),SL(J),SL1(J))
        PDU=0
        IF(M.EQ.0 .AND. NPLK1.GT.0) THEN
          CALL PKINT(AM1U,TT,TT1,TT,NPLK1,PKI)
          DO 6 J=1,NPLK1
          SUM1=0
          DO 5 K=1,NDA
    5     SUM1=SUM1+PT1(I,K)*DPE(K,J,L)+PR1(I,K)*DME(K,J,L)
    6     PDU=PDU+(W*SUM1+2*PI*(1-W)*CPLK(J,L))*PKI(J)
        ENDIF
        DO 7 J=1,NA0
        TRNS0=-T1/AM0(J)
        TRNS0=EXPFN(TRNS0)*FSOL
        PVU=PT10(I,J)*TRNS0
        HU=0
        DO 8 K=1,NDA
        HU=HU+(H1(K)*C  (K)-H2(K)*SL(K))*ALFA(K,J)
     &       +(H1(K)*SL1(K)-H2(K)*C (K))*BETA(K,J)
    8   PVU=PVU+PT1(I,K)*VPE(K,J,L)+PR1(I,K)*VME(K,J,L)
        PVU=W*PVU*EI (J)
        AI(J)=(HU+PVU+PDU)/AM1
    7   CONTINUE
       ELSE
C UPWARD INTENSITY
        DO 13 J=1,NA0
   13   CALL EXINT(AM1U,AM0(J),TT,TT,TT2,EI(J))
        DO 14 J=1,NDA
   14   CALL CSINT(AM1U,ZEE(J,L),TT2,TT,TT,TT2,C(J),SL(J),SL1(J))
        PDU=0
        IF(M.EQ.0 .AND. NPLK1.GT.0) THEN
          CALL PKINT(AM1U,TT,TT,TT2,NPLK1,PKI)
          DO 16 J=1,NPLK1
          SUM1=0
          DO 15 K=1,NDA
   15     SUM1=SUM1+PR1(I,K)*DPE(K,J,L)+PT1(I,K)*DME(K,J,L)
   16     PDU=PDU+(W*SUM1+2*PI*(1-W)*CPLK(J,L))*PKI(J)
        ENDIF
        DO 17 J=1,NA0
        TRNS0=-T1/AM0(J)
        TRNS0=EXPFN(TRNS0)*FSOL
        PVU=PR10(I,J)*TRNS0
        HU=0
        DO 18 K=1,NDA
        HU=HU+(H1(K)*C  (K)+H2(K)*SL(K))*ALFA(K,J)
     &       +(H1(K)*SL1(K)+H2(K)*C (K))*BETA(K,J)
   18   PVU=PVU+PR1(I,K)*VPE(K,J,L)+PT1(I,K)*VME(K,J,L)
        PVU=W*PVU*EI (J)
        AI(J)=(HU+PVU+PDU)/AM1
   17   CONTINUE
      ENDIF
      RETURN
      END
      FUNCTION EXX(EX,T1,T2,X)
C EXP(EX)*INTEG(FROM T1 TO T2) DT EXP(XT)
C--- HISTORY
C 89.11. 3 MODIFIED
      SAVE EPS,EXMN,INIT
      DIMENSION CCP(3)
      DATA INIT/1/
      IF(INIT.EQ.1) THEN
        INIT=0
        CALL CPCON(CCP)
        EPS=CCP(1)*100
        EXMN=CCP(2)*0.8*2.3
      ENDIF
      EX1=EX+X*T1
      EX2=EXPFN(EX1)
      T22=T2-T1
      IF(ABS(T22*X).LE.EPS) THEN
         EXX=EX2*(T22+X*T22**2/2.0)
        ELSE
         EX3=EX+X*T2
         IF(EX3.LE.EXMN) THEN
           EXX=-EX2/X
          ELSE
           EXX=(EXP(EX3)-EX2)/X
         ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE EXINT(AM1U,AM0,DPT,T1,T2,C)
C INTEGRAL(FROM T1 TO T2) DT EXP(-(THK-T)/AM1U-T/AM0)
C--- HISTORY
C 90. 1.13 CREATED
C--- INPUT
C AM1U     R      USERDEFINED EMERGENT MU.
C                 IF .GT. 0 THEN DOWN,  IF .LT. 0 THEN UPWARD.
C AM0      R      SOLAR DIRECTION.
C DPT      R      OPTICAL DEPTH FOR INTERPOLATION.
C T1       R      LOWER LIMIT OF INTEGRATION.
C T2       R      UPPER LIMIT OF INTEGRATION.
C--- OUTPUT
C C        R      INTEGRATION.
      EX=-DPT/AM1U
      X=1/AM1U-1/AM0
      C=EXX(EX,T1,T2,X)
      RETURN
      END
      SUBROUTINE CSINT(AM1U,ZE,THK,DPT,T1,T2,C,SL,SL1)
C INTEGRAL(FROM T1 TO T2) DT EXP(-(DPT-T)/AM1U) C(T, THK)
C   WHERE C(T, THK) = (EXP(-Z*(THK-T)+EXP(-Z*T))/2
C WE ALSO DEFINE THE INTEGRAL FOR THE FUNCTIONS: S*Z, S/Z.
C--- HISTORY
C 90. 1.13 CREATED
C--- INPUT
C AM1U     R      USERDEFINED EMERGENT MU.
C                 IF .GT. 0 THEN DOWN,  IF .LT. 0 THEN UPWARD.
C ZE       R      EIGENVALUE (L).
C THK      R      OPTICAL THICKNESS OF LAYER.
C DPT      R      OPTICAL DEPTH FOR INTERPOLATION.
C T1       R      LOWER LIMIT OF INTEGRATION.
C T2       R      UPPER LIMIT OF INTEGRATION.
C--- OUTPUT
C C        R      INTEGRATION OF C.
C SL       R      INTEGRATION OF S*L.
C SL1      R      INTEGRATION OF S/L.
      SAVE EPS,INIT
      DIMENSION CCP(3)
      DATA INIT/1/
      IF(INIT.EQ.1) THEN
        INIT=0
        CALL CPCON(CCP)
        EPS=CCP(1)*100
      ENDIF
C C
      EX=-DPT/AM1U-ZE*THK
      X=1/AM1U+ZE
      EX1=EXX(EX,T1,T2,X)
      EX=-DPT/AM1U
      X=1/AM1U-ZE
      EX2=EXX(EX,T1,T2,X)
      C =(EX1+EX2)/2
      IF(ABS(ZE*T1).GT.EPS .OR. ABS(ZE*T2).GT.EPS) THEN
        SL1=(EX1-EX2)/2.0/ZE
       ELSE
        EX=(T1-DPT)/AM1U
        EX1=EXPFN(EX)
        EX=(T2-DPT)/AM1U
        EX2=EXPFN(EX)
        EX=-DPT/AM1U
        X=1/AM1U
        EX3=EXX(EX,T1,T2,X)
        SL1=AM1U*(T2*EX2-T1*EX1)-(AM1U+THK/2)*EX3
      ENDIF
      SL =ZE**2*SL1
      RETURN
      END
      SUBROUTINE PKINT(AM1U,DPT,T1,T2,NPLK1,C)
C INTEGRAL(FROM T1 TO T2) DT EXP(-(DPT-T)/AM1U) T**N
C--- HISTORY
C 90. 1.13 CREATED
C--- INPUT
C AM1U     R      USERDEFINED EMERGENT MU.
C                 IF .GT. 0 THEN DOWN,  IF .LT. 0 THEN UPWARD.
C DPT      R      OPTICAL DEPTH FOR INTERPOLATION.
C T1       R      LOWER LIMIT OF INTEGRATION.
C T2       R      UPPER LIMIT OF INTEGRATION.
C NPLK1    I      .GE. 1
C--- OUTPUT
C C     R(NPLK1)  INTEGRATIONS
C
      DIMENSION C(*)
      EX=-DPT/AM1U
      X=1/AM1U
      C(1)=EXX(EX,T1,T2,X)
      IF(NPLK1.GE.2) THEN
        EX=(T1-DPT)/AM1U
        EX1=EXPFN(EX)
        EX=(T2-DPT)/AM1U
        EX2=EXPFN(EX)
        DO 1 N1=2,NPLK1
        EX1=T1*EX1
        EX2=T2*EX2
    1   C(N1)=AM1U*(EX2-EX1-(N1-1)*C(N1-1))
      ENDIF
      RETURN
      END
      SUBROUTINE CSFN(TC,T,ZE,C,SL,SL1)
C GET C SL AND SL1 FUNCTIONS.
C--- HISTORY
C 89.11. 2  CREATED
C---
C 2 C = E(TC-T) + E(T),   2S = E(TC-T) - E(T)
C   E(T) = EXP(-ZE*T)
C SL = S L,   SL1 = S L**-1
C
      SAVE INIT,EPS
      DIMENSION CCP(3)
      DATA INIT/1/
      IF(INIT.EQ.1) THEN
        INIT=0
        CALL CPCON(CCP)
        EPS=CCP(1)*10
      ENDIF
      EXP1=-ZE*(TC-T)
      EXP2=-ZE*T
      EXP3=EXPFN(EXP1)
      EXP4=EXPFN(EXP2)
      C  =(EXP3+EXP4)/2
      S  =(EXP3-EXP4)/2
      SL =ZE*S
      IF(ABS(EXP1).LE.EPS .AND. ABS(EXP2).LE.EPS) THEN
        SL1 =T-TC/2
        ELSE
        SL1 =S/ZE
      ENDIF
      RETURN
      END
      SUBROUTINE LGNDF3(LMAX1,N,X,Y,G)
C
C LEGENDRE EXPANSION (SAME AS LGNDF2 BUT GENERATING G-MOMENTS).
C--- HISTORY
C 90. 1.20 CREATED FROM LGNDF2, USE EXPFN.
C          DIRECTION OF INTEGRATION FROM X(1) TO X(N).
C 91. 2.16 STRIP NG FROM SAVE STATEMENT.
C 92. 4. 3 KILL THE STATEMENT OF GW=GW/2 AFTER QGAUSN
C     6.22 ADD GW/2 AGAIN BECAUSE THE ABOVE CHANGE IS MISTAKE
C--- INPUT
C LMAX1      I      MAXIMUM ORDER + 1.
C N          I      NUMBER OF DATA ON (-1, 1). .GE. 4.
C X        R(NA)    INDEPENDENT VARIABLES ON (-1, 1)
C Y        R(NA)    Y(X).
C--- OUTPUT
C G    R(LMAX1)     Y = SUM(L1=1,LMAX1) (2*L1-1)*G(L1)*PL(L1)
C                      WHERE PL(L1) IS (L1-1)TH ORDER LEGENDRE
C                      POLYNOMIAL.
C$ENDI
      SAVE INIT,GW,GX
C VARIABLES FOR THE ROUTINE.
      DIMENSION X(N),Y(N),G(LMAX1)
C WORKING AREAS.
      PARAMETER (PI=3.141592653589793, RAD=PI/180.0)
      PARAMETER (NG=5)
      DIMENSION GX(NG),GW(NG)
      DATA INIT/1/
C SHIFTED GAUSSIAN QUADRATURE.
      IF(INIT.GE.1) THEN
        INIT=0
        CALL QGAUSN(GW,GX,NG)
        DO 1 I=1,NG
    1   GW(I)=GW(I)/2
      ENDIF
C CLEAR
      DO 2 L1=1,LMAX1
    2 G(L1)=0
C LOOP FOR ANGLE
      DO 3 I=1,N-1
C CUBIC INTERPOLATION
      IF(I .LE. 2) THEN
        I1=1
        I4=4
       ELSE
        IF(I .LE. N-2) THEN
          I1=I-1
          I4=I+2
         ELSE
          I1=N-3
          I4=N
        ENDIF
      ENDIF
CC
      I2=I1+1
      I3=I2+1
      X1=X(I1)
      X2=X(I2)
      X3=X(I3)
      X4=X(I4)
      IF(   (Y(I1) .LE. 0) .OR. (Y(I2) .LE. 0) .OR. (Y(I3) .LE. 0)
     & .OR. (Y(I4) .LE. 0) ) THEN
        ISIGN=-1
        ALP1=Y(I1)
        ALP2=Y(I2)
        ALP3=Y(I3)
        ALP4=Y(I4)
       ELSE
        ISIGN=1
        ALP1=LOG(Y(I1))
        ALP2=LOG(Y(I2))
        ALP3=LOG(Y(I3))
        ALP4=LOG(Y(I4))
      ENDIF
C LOOP FOR GAUSSIAN INTEGRATION
      DO 4 J=1,NG
CC INTERPOLATED VALUE OF Y
      XX=X(I)+GX(J)*(X(I+1)-X(I))
      WW=GW(J)*(X(I+1)-X(I))
      PP=(XX-X2)*(XX-X3)*(XX-X4)/(X1-X2)/(X1-X3)/(X1-X4)*ALP1
     &  +(XX-X1)*(XX-X3)*(XX-X4)/(X2-X1)/(X2-X3)/(X2-X4)*ALP2
     &  +(XX-X1)*(XX-X2)*(XX-X4)/(X3-X1)/(X3-X2)/(X3-X4)*ALP3
     &  +(XX-X1)*(XX-X2)*(XX-X3)/(X4-X1)/(X4-X2)/(X4-X3)*ALP4
      IF(ISIGN .EQ. 1) PP=EXPFN(PP)
C LEGENDRE SUM
      PL1=0
      PL=1
      G(1)=G(1)+PP*WW
      IF(LMAX1.GE.2) THEN
        DO 5 L1=2,LMAX1
        PL2=PL1
        PL1=PL
        PL=((2*L1-3)*XX*PL1-(L1-2)*PL2)/(L1-1)
    5   G(L1)=G(L1)+PP*PL*WW
      ENDIF
    4 CONTINUE
    3 CONTINUE
      RETURN
      END
      FUNCTION PLGD(INIT,X)
C LEGENDRE POLYNOMIALS
C--- HISTORY
C 87.11.12 CREATED
C 90. 1.16 SAVE STATEMENT
C--- INPUT
C INIT   I    IF 1 THEN L=0
C             IF 0 THEN L=L+1
C X      R    (-1,1)
C--- OUT
C INIT=0
C$ENDI
C
      SAVE L,PL,PL1,PL2
      IF(INIT.GT.0) THEN
        INIT=0
        L=-1
      ENDIF
      L=L+1
      IF(L.EQ.0) THEN
        PL=1
       ELSE
        IF(L.EQ.1) THEN
          PL1=PL
          PL=X
         ELSE
          PL2=PL1
          PL1=PL
          PL=((2*L-1)*X*PL1-(L-1)*PL2)/REAL(L)
        ENDIF
      ENDIF
      PLGD=PL
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
      SUBROUTINE QUADA(NDA,AMUA,WA)
C CONSTRUCTION OF QUADRATURE IN THE ATMOSPHERE.
C------ HISTORY
C 1986.10.2  USE QGAUSN FROM STAMNES AND TSAY.
C--- INPUT
C NDA        I      NO. OF QUADRATURE STREAMS IN THE HEMISPHERE OF ATMOS
C--- OUTPUT
C AMUA    R(KNDA)   MUA(I), I=1,N1    MUA(1) > MUA(2) > ...
C WA      R(KNDA)   CORRESPONDING WEIGHTS.
C$ENDI
C--- VARIABLES FOR THE ROUTINE
      DIMENSION AMUA(NDA),WA(NDA)
C SHIFTED GAUSSIAN QUADRATURE ON (0, 1).
      CALL QGAUSN(WA,AMUA,NDA)
C REORDERING.
      DO 3 I=1,NDA/2
      I1=NDA+1-I
      X=AMUA(I)
      AMUA(I)=AMUA(I1)
      AMUA(I1)=X
      X=WA(I)
      WA(I)=WA(I1)
    3 WA(I1)=X
      RETURN
      END
      FUNCTION SGLR (AM1, AM0, NL, TC, W, P, UT)
C SINGLY SCATTERED INTENSITY (MULTI-LAYER) BY EQ. (12) OF NT.
C--- REFERENCE
C NT:  T. NAKAJIMA AND M. TANAKA, 1988, JQSRT, 40, 51-69
C--- HISTORY
C 89. 2.22   CREATED BY T. NAKAJIMA
C              FROM AIDS BY CHANGING THE MEANING OF T, W, P,
C              AND UT.  UT IS A LEVEL INSIDE THE MULTI-LAYER SYSTEM.
C     5. 4   PATCH FOR EXPONENTIAL OVERFLOW.
C     6.20   REFORM STYLE.
C    12.20   WITH CPEX(3)
C 90.12. 1   PUT IF(IALM... OUT OF DO-LOOP (FOR OPTIMIZATION)
C FOLLOWINGS ARE DEFINITIONS OF SYMBOLS IN THE ROUTINE, WITH THE FORMAT:
C NAME    TYPE       CONTENT
C WHERE TYPE = S (SUBROUTINE), F (FUNCTION), I (INTEGER), R (REAL).
C FOR THE TYPE = I AND R, THE SIZE OF ARRAY IS SHOWN BY ( ).
C--- INPUT
C AM1      R         COS(EMERGENT ZENITH ANGLE) WITH THE SIGN OF
C                    TRANSMITTED( .GT.0),  REFLECTED( .LT.0).
C AM0      R         COS(INCIDENT ZENITH ANGLE) .GT. 0.
C NL       I         NO. OF SUBLAYERS.
C TC     R(NL)       OPTICAL THICKNESS OF SUBLAYER.
C W      R(NL)       SINGLE SCATTERING ALBEDO OF SUBLAYER.
C P      R(NL)       PHASE FUNCTION OF SUBLAYER.
C                    INTEGRATION OVER UNIT SPHERE SHOULD BE 1.
C UT       R         USER DEFINED OPTICAL DEPTH FOR THE INTENSITY
C                    WITHIN THE MULTI-LAYERED SYSTEM.
C--- OUTPUT
C SGLR     F         SINGLY SCATTERED INTENSITY.
      SAVE INIT,EPS
      DIMENSION TC(NL),W(NL),P(NL),CCP(3)
      LOGICAL IALM,IST
      DATA INIT/1/
      IF(INIT.GT.0) THEN
        INIT=0
        CALL CPCON(CCP)
        EPS=CCP(1)*30
      ENDIF
C SET EPS: IF ABS(1/AM1 - 1/AM0) .LE. EPS THEN
C                      THE ROUTINE SETS ALMUCANTAR CONDITION-IALM.
      SGLR=0
      X=1/AM1-1/AM0
      IF(ABS(X).LE.EPS) THEN
        IALM=.TRUE.
        ELSE
        IALM=.FALSE.
      ENDIF
      T2=0
      IST=.TRUE.
C X=0 (TRANSMISSION)
      IF(IALM) THEN
        DO 11 L=1,NL
        T1=T2
        T2=T1+TC(L)
        IF(UT.LE.T1) GOTO 13
        IF(UT.LT.T2) T2=UT
        IF(IST) THEN
          IST=.FALSE.
          E2=T1
        ENDIF
        E1=E2
        E2=T2
        SGLR=SGLR+W(L)*P(L)*(E2-E1)
 11     CONTINUE
       ELSE
C X<>0
      DO 21 L=1,NL
        T1=T2
        T2=T1+TC(L)
CC REFLECTION
        IF(AM1.LT.0) THEN
          IF(UT.GE.T2) GOTO 21
          IF(UT.GT.T1) T1=UT
CC TRANSMISSION
         ELSE
          IF(UT.LE.T1) GOTO 13
          IF(UT.LT.T2) T2=UT
        ENDIF
C
        IF(IST) THEN
          IST=.FALSE.
          E2=T1*X
          E2=EXPFN(E2)/X
        ENDIF
        E1=E2
        E2=T2*X
        E2=EXPFN(E2)/X
        SGLR=SGLR+W(L)*P(L)*(E2-E1)
 21   CONTINUE
      ENDIF
 13   E1=-UT/AM1
      SGLR=SGLR*EXPFN(E1)/ABS(AM1)
      RETURN
      END
      FUNCTION EXPFN(X)
C EXPONENTIAL FUNCTION WITH OVER/UNDER FLOW SETTING.
C--- HISTORY
C 89. 5. 4   CREATED BY T. NAKAJIMA.
C 90. 1.17   UPDATED WITH CPCON
C--- INPUT
C X        R         INDEPENDENT VARIABLE.
C--- OUTPUT
C EXPFN    F         EXP(X).
C                     IF X.LE. VMN THEN EXP(X) IS RESET AS   0.
C                     IF X.GE. VMX THEN EXP(X) IS RESET AS EXP(VMX).
C--- PARAMETERS
C SYSTEM SET THE -VMN- AND -VMX- BY THE FUNCTION R1MACH.
C
      SAVE INIT,VMN,VMX,EXPMN,EXPMX
      DIMENSION CCP(3)
      DATA INIT/1/
C
C SET VMN      R         ENDERFLOW LIMIT.
C     VMX      R         OVERFLOW LIMT.
      IF(INIT.GT.0) THEN
        INIT=0
        CALL CPCON(CCP)
        VMN=CCP(2)*0.8*2.3
        VMX=CCP(3)*0.8*2.3
        EXPMN=0
        EXPMX=EXP(VMX)
      ENDIF
C
      IF(X.LE.VMN) THEN
        EXPFN=EXPMN
        ELSE
        IF(X.GE.VMX) THEN
          EXPFN=EXPMX
          ELSE
          EXPFN=EXP(X)
        ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE CSPL1(N,X,Y,A,B,C,D)
C GETTING COEFFICIENTS OF NATURAL CUBIC SPLINE FITTING
C USE A LINEAR INTERPOLATION OUTSIDE THE MEANINGFUL RANGE OF X
C X-increasing order
C--- HISTORY
C 88. 1. 4  CREATED
C--- INPUT
C N      I     NBR OF DATA
C X    R(N)    INDEPENDENT VARIABLE DATA
C              X-increasing order
C Y    R(N)      DEPENDENT VARIABLE DATA
C--- OUTPUT
C A    R(N)    LAMBDA -> A   WHERE  Y=A+X*(B+X*(C+D*X))
C B    R(N)    D      -> B   I-TH FOR RANGE (X(I-1), X(I))
C C    R(N)    M         C   (A,B,C,D FOR I=1 ARE SAME AS THOSE FOR I=2)
C D    R(N)    M         D
C
C--- NOTES
C REF-1   P. F. DAVIS AND PHILIP RABINOWITZ (1984)
C         METHODS PF NUMERICAL INTEGRATION, SECOND EDITION
C         ACADEMIC PRESS, INC., PP612.
C$ENDI
C
      DIMENSION X(N),Y(N),A(N),B(N),C(N),D(N)
C
      GOTO (1,2,3,4), N
C N > 4
      A(1)=0.0
      A(N)=1.0
      B(1)=0.0
      B(N)=0.0
      H2=X(2)-X(1)
      S2=(Y(2)-Y(1))/H2
      DO 5 I=2,N-1
      H1=H2
      H2=X(I+1)-X(I)
      S1=S2
      S2=(Y(I+1)-Y(I))/H2
      A(I)=1.0/(1.0+H1/H2)
    5 B(I)=6.0*(S2-S1)/(H2+H1)
CC
      Q2=0.0
      U2=0.0
      DO 6 I=1,N
      Q1=Q2
      U1=U2
      P2=(1.0-A(I))*Q1+2.0
      Q2=-A(I)/P2
      U2=(B(I)-(1.0-A(I))*U1)/P2
      A(I)=Q2
    6 B(I)=U2
      C(N)=B(N)
      DO 7 I=N-1,1,-1
    7 C(I)=A(I)*C(I+1)+B(I)
CC
      X2=X(1)
      C2=C(1)/6.0
      DO 8 I=2,N
      H2=X(I)-X(I-1)
      X1=X2
      X2=X(I)
      C1=C2
      C2=C(I)/6.0
      P1=Y(I-1)/H2-C1*H2
      P2=Y(I  )/H2-C2*H2
      A(I)=    ( C1*X2**3-C2*X1**3)/H2+P1*X2-P2*X1
      B(I)=3.0*(-C1*X2**2+C2*X1**2)/H2-P1   +P2
      C(I)=3.0*( C1*X2   -C2*X1   )/H2
    8 D(I)=    (-C1      +C2      )/H2
      GOTO 11
C N=1
    1 A(1)=Y(1)
      B(1)=0.0
      C(1)=0.0
      D(1)=0.0
      RETURN
C N=2
    2 X1=X(1)
      X2=X(2)
      Z1=Y(1)/(X1-X2)
      Z2=Y(2)/(X2-X1)
      A(2)=-X2*Z1-X1*Z2
      B(2)=Z1+Z2
      C(2)=0.0
      D(2)=0.0
      GOTO 11
C N=3
    3 X1=X(1)
      X2=X(2)
      X3=X(3)
      Z1=Y(1)/(X1-X2)/(X1-X3)
      Z2=Y(2)/(X2-X3)/(X2-X1)
      Z3=Y(3)/(X3-X1)/(X3-X2)
      A(2)=X2*X3*Z1+X3*X1*Z2+X1*X2*Z3
      B(2)=-(X2+X3)*Z1-(X3+X1)*Z2-(X1+X2)*Z3
      C(2)=Z1+Z2+Z3
      D(2)=0.0
      A(3)=A(2)
      B(3)=B(2)
      C(3)=C(2)
      D(3)=D(2)
      GOTO 11
C N=4
    4 X1=X(1)
      X2=X(2)
      X3=X(3)
      X4=X(4)
      Z1=Y(1)/(X1-X2)/(X1-X3)/(X1-X4)
      Z2=Y(2)/(X2-X3)/(X2-X4)/(X2-X1)
      Z3=Y(3)/(X3-X4)/(X3-X1)/(X3-X2)
      Z4=Y(4)/(X4-X1)/(X4-X2)/(X4-X3)
      A(2)=-X2*X3*X4*Z1-X3*X4*X1*Z2-X4*X1*X2*Z3-X1*X2*X3*Z4
      B(2)=(X2*X3+X3*X4+X4*X2)*Z1+(X3*X4+X4*X1+X1*X3)*Z2
     &    +(X4*X1+X1*X2+X2*X4)*Z3+(X1*X2+X2*X3+X3*X1)*Z4
      C(2)=-(X2+X3+X4)*Z1-(X3+X4+X1)*Z2-(X4+X1+X2)*Z3-(X1+X2+X3)*Z4
      D(2)=Z1+Z2+Z3+Z4
      DO 10 I=3,4
      A(I)=A(2)
      B(I)=B(2)
      C(I)=C(2)
   10 D(I)=D(2)
C
   11 A(1)=A(2)
      B(1)=B(2)
      C(1)=C(2)
      D(1)=D(2)
      RETURN
      END
      FUNCTION CSPLI(X1,N,X,A,B,C,D)
C GET INTERPOLATED VALUE USING CUBIC SPLINE (PAIR WITH CSPL1)
C--- HISTORY
C 88. 6. 6  REGISTERED BY T. NAKAJIMA
C--- INPUT
C X1      R     INTERPOLATION POINT
C N       I     NO. OF INTERVALS + 1
C X     R(N)    DIVISION POINTS OF THE INTERVALS.
C A, B, C, D
C       R(N)    Y=A+X*(B+X*(C+X*D)))
C--- OUTPUT
C CSPLI  RF     INTERPOLATED VALUE
C$ENDI
      DIMENSION X(N),A(N),B(N),C(N),D(N)
      DO 7 J=1,N
      IF(X1.LE.X(J)) GOTO 8
    7 CONTINUE
      J=N
    8 CSPLI=A(J)+X1*(B(J)+X1*(C(J)+X1*D(J)))
      RETURN
      END
      SUBROUTINE EQ12(A,B,NI,J,NB1)
C  A(*) = B(*, J)
C--- HISTORY
C 88. 6. 6  REGISTERED BY T. NAKAJIMA
C--- INPUT
C B      R(NB1,*)    2-DIM ARRAY  B.
C NI        I        A(I) = B(I,J),   I=1,NI
C J         I
C--- OUTPUT
C A       R(NI)      1-DIM ARRAY   A=B(*,J)
C$ENDI
      DIMENSION A(NI),B(NB1,J)
      DO 1 I=1,NI
    1 A(I)=B(I,J)
      RETURN
      END
      SUBROUTINE EQ21(A,B,NI,J,NA1)
C  A(*,J) = B(*)
C--- HISTORY
C 88. 6. 6  REGISTERED BY T. NAKAJIMA
C B      R(NI)     SOURCE 1-DIM ARRAY  B.
C NI       I       A(I,J) = B(I),  I=1,NI
C J        I
C MA1      I       SIZE FOR DIMENSION A(NA1,*)
C--- OUTPUT
C A     R(NA1,*)   DESTINATION 2-DIM ARRAY A(*,J)=B
C$ENDI
      DIMENSION A(NA1,J),B(NI)
      DO 1 I=1,NI
    1 A(I,J)=B(I)
      RETURN
      END
      SUBROUTINE EQ32(A,B,NI,NJ,K,NA1,NA2,NA3,NB1,NB2)
C  A(*,*,K)=B(*,*)
C--- HISTORY
C 88. 5. 6   REGISTERED BY T. NAKAJIMA
C--- INPUT
C B      R(NB1,NB2)    SOURCE 2-DIM ARRAY
C NI,NJ,K   I          A(I,J,K)=B(I,J),   I=1,NI;  J=1,NJ
C NA1,NA2,NA3  I       DIM A(NA1,NA2,NA3)
C NB1,NB2      I       DIM B(NB1,NB2)
C--- OUTPUT
C A    R(NA1,NA2,NA3)  DESTINATION 3-DIM ARRAY A(*,*,K)=B(*,*)
C$ENDI
      DIMENSION A(NA1,NA2,NA3),B(NB1,NB2)
      DO 1 J=1,NJ
      DO 1 I=1,NI
    1 A(I,J,K)=B(I,J)
      RETURN
      END
      SUBROUTINE GRNDL3(FSOL,GALB,BGND,DPT,M,N1,AM,W,N0,AM0,R,T,ER,ET)
C--- HISTORY
C 95. 5.26 Generated from GRNDL introducing AM and W
C            and R=2*GALB... -> SUMWM*GALB...
C--- INPUT
C FSOL       R           SOLAR IRRADIANCE AT THE SYSTEM TOP.
C GALB       R           FLUX REFLECTIVITY OF THE SURFACE.
C BGND       R           THERMAL INTENSITY FROM THE SURFACE=(1-r)*PLANK.
C DPT        R           OPTICAL DEPTH AT THE SURFACE.
C M          I           FOURIER ORDER.
C N1         I           NUMBER OF QUADRATURE POINTS.
C AM      R(KNDM)        MU(I)
C W       R(KNDM)        W (I)
C N0         I           NUMBER OF SOLAR ANGLES.
C AM0     R(KNA0)        COS(SOLAR ZENITH ANGLE)
C--- OUTPUT
C R      R(KNDM,KNDM)    SCALED REFELCTION MATRIX.
C T      R(KNDM,KNDM)    SCALED TRANSMISSION MATRIX.
C ER     R(KNDM,KNA0)    SCALED SOURCE MATRIX FOR UPWARD RADIANCE.
C ET     R(KNDM,KNA0)    SCALED SOURCE MATRIX FOR DOWNWARD RADIANCE.
C
C LAMBERT SURFACE
      PARAMETER (PI=3.141592654)
      PARAMETER (KNDM  =16)
      PARAMETER (KNA0  =1)
      DIMENSION R(KNDM,KNDM),T(KNDM,KNDM),ER(KNDM,KNA0),ET(KNDM,KNA0)
     & ,AM0(KNA0),AM(KNDM),W(KNDM)
      SUMWM=0
      DO 3 I=1,N1
    3 SUMWM=SUMWM+AM(I)*W(I)
      IF(GALB.LE.0.0 .OR. M.GT.0) THEN
        CALL EQ20( R,0.0,N1,N1,KNDM,KNDM)
        CALL EQ20( T,0.0,N1,N1,KNDM,KNDM)
       ELSE
        DO 2 I=1,N1
        DO 2 J=1,N1
        T (I,J)=0
    2   R (I,J)=GALB*SQRT(AM(I)*W(I)*AM(J)*W(J))/SUMWM
      ENDIF
      IF(M.GT.0) THEN
        CALL EQ20(ER,0.0,N1,N0,KNDM,KNA0)
        CALL EQ20(ET,0.0,N1,N0,KNDM,KNA0)
       ELSE
        DO 1 J=1,N0
          TRNS0=EXP(-DPT/AM0(J))*FSOL
          X=GALB*AM0(J)*TRNS0/SUMWM+2*PI*BGND
          DO 1 I=1,N1
          ET(I,J)=0
    1     ER(I,J)=SQRT(AM(I)*W(I))*X
      ENDIF
      RETURN
      END
      SUBROUTINE PHAS2(M1,MMAX1,N1,N2,KNP,KN1,KN2,PT,PR,G
     & ,PL1,PL2)
C--- HISTORY
C 90. 1.20  GENERATED FROM PHASE (SAME AS PHASE BUT USING THE MOMENTS G)
C--- INPUT
C M1         I       M + 1    FOURIER ORDER + 1.
C MMAX1      I       MMAX+1   MAX FOURIER ORDER + 1.
C N1         I       NO. OF EMERGENT ZENITH ANGLES.
C N2         I       NO. OF INCIDENT ZENITH ANGLES.
C KNP        I       SIZE OF N1 FOR ARRAY PT AND PR.
C KN1        I       SIZE OF N1 FOR ARRAY PL1.
C KN2        I       SIZE OF N2 FOR ARRAY PL2.
C G     R(MMAX1)     LEGENDRE MOMENTS OF PHASE FUNCTION.
C                      G(1)=1
C PL1   R(KN1,MMAX1) LEGENDRE POLYNOMIALS FOR EMERGENT DIRECTION.
C                    (I,M1) = (EMERGENT DIRECTIONS, ORDER+1).
C PL2   R(KN2,MMAX1) LEGENDRE POLYNOMIALS FOR INCIDENT DIRECTION.
C--- OUTPUT
C PT    R(KNP,N2)    PHASE MATRIX FOR TRANSMISSION.
C PR    R(KNP,N2)    PHASE MATRIX FOR REFLECTION.
C$ENDI
C
      PARAMETER (PI=3.141592653)
      DIMENSION PT(KNP,N2),PR(KNP,N2),G(MMAX1)
     & ,PL1(KN1,MMAX1),PL2(KN2,MMAX1)
C
      CALL EQ20(PT,0.0,N1,N2,KNP,N2)
      CALL EQ20(PR,0.0,N1,N2,KNP,N2)
      IF(M1.LE.MMAX1) THEN
        DO 1 J=1,N2
        DO 1 I=1,N1
        SIGN=-1
        DO 1 K1=M1,MMAX1
        SIGN=-SIGN
        C4=(2*K1-1)/2.0*G(K1)*PL1(I,K1)*PL2(J,K1)
        PT(I,J)=PT(I,J)+C4
    1   PR(I,J)=PR(I,J)+SIGN*C4
      ENDIF
      RETURN
      END
      SUBROUTINE PLGND(M1,MMX1,NX,NX0,X,PL)
C NORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
C------ HISTORY
C 1987.03.04
C------ INPUT VARIABLES
C VARIABLE  TYPE    INTERPRETATION
C M1         I      ORDER OF THE FOURIER SERIES + 1
C MMX1       I      MAX ORDER OF M  + 1
C NX         I      NBR OF X
C NX0        I      DECLARED NBR OF NX
C X       R(NX)     XI, I=1,NX
C------ OUTPUT VARIABLES
C PL      R(NX0,    NORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS
C           MMX1)
C$ENDI
C------ GLOBAL PARAMETER
      PARAMETER (PI=3.141592653)
C------ VARIABLES FOR THE ROUTINE
      DIMENSION X(NX),PL(NX0,MMX1)
C------ NORMALIZED ASSOCIATED LEGENDRE POLYNOMIALS.
C  CC=P(M,L,X)*SQRT((L-M)|/(L+M)|), WHERE X= COS(MU).
      M=M1-1
      DO 15 L1=M1,MMX1
      L=L1-1
      IF(L-(M+1)) 16,17,18
   16 IF(M) 21,21,22
   21 DO 2 I=1,NX
    2 PL(I,L1)=1.0
      GOTO 15
   22 K=2
      ETA=1.0
      DO 1 J=1,M
      EPSI=1.0-1.0/REAL(K)
      ETA=ETA*EPSI
    1 K=K+2
      ETA=SQRT(ETA)
      DO 24 I=1,NX
      IF(X(I).GE.1.0) THEN
      PL(I,L1)=0.0
      ELSE
      EXPCC=(REAL(M)*LOG(ABS(1.0-X(I)**2)))/2.0
      IF(EXPCC.LE.-100.0) THEN
      PL(I,L1)=0.0
      ELSE
      PL(I,L1)=ETA*EXP(EXPCC)
      ENDIF
      ENDIF
   24 CONTINUE
      GO TO 15
   17 C0=SQRT(REAL(2*M+1))
      DO 26 I=1,NX
   26 PL(I,L1)=C0*PL(I,L1-1)*X(I)
      C0=SQRT(REAL((L-M)*(L+M)))
      GO TO 15
   18 C1=SQRT(REAL((L-M)*(L+M)))
      C2=REAL(2*L-1)
      DO 27 I=1,NX
   27 PL(I,L1)=(C2*X(I)*PL(I,L1-1)-C0*PL(I,L1-2))/C1
      C0=C1
   15 CONTINUE
      RETURN
      END
      SUBROUTINE TRN1(NSB,NDD,NA0,IUP,IDN,RE,TE,SER,SET,RUP,RDN,ERR)
C SOLVE THE RADIATIVE TRANSFER IN THE MULTI-SUBLAYER SYSTEM
C  BY THE ADDING METHOD.
C--- HISTORY
C 87. 3. 4 FOR INTENSITY CALCULATION.
C 89.10.30 RE-EDIT.
C 96.01.19   EQ22(S2U,ST,N2,NA0,KNDM,KNDM,KNDM,KNDM)
C          ->EQ22(S2U,ST,N2,NA0,KNDM,KNA0,KNDM,KNA0)
C--- INPUT
C NSB       I              NUMBER OF SUBLAYERS.
C NDD    R(NSB1)           NUMBER OF QUADRATURE POINTS AT INTERFACES.
C IUP    I(NSB)            ELEMENT NUMBER OF UPWELLING OPERATORS.
C IDN    I(NSB)            ELEMENT NUMBER OF DOWNGOING OPERATORS.
C NA0       I              NUMBER OF SOLAR ZENITH ANGLES.
C RE    R(KNDM,KNDM,KNLT)  REFLECTION MATRICES OF SUBLAYERS.
C TE    R(KNDM,KNDM,KNLT)  TRANSMISSION MATRICES OF SUBLAYERS.
C SER   R(KNDM,KNA0,KNSB)  UPWELLING SOURCE MATRICES.
C SET   R(KNDM,KNA0,KNSB)  DOWNGOING SOURCE MATRICES.
C--- OUTPUT
C RUP   R(KNDM,KNA0,KNSB1) UPWELLING INTERNAL INTENSITIES.
C RDN   R(KNDM,KNA0,KNSB1) DOWNGOING INTERNAL INTENSITIES.
C ERR      C*64            ERROR INDEX.
C--- PARAMETER
C KNA0      I              NUMBER OF SOLAR ZENITH ANGLES.
C KNDM      I              NUMBER OF QUADRATURE POINTS.
C KNSB      I              NUMBER OF SUBLAYERS.
C KNLT      I              TOTAL NUMBER OF ELEMENTARY OPERATORS
C                           TAKING PORALITY INTO ACCOUNT .GE. KNSB.
C--- AREAS FOR THIS ROUTINE
      PARAMETER (KNA0  =1)
      PARAMETER (KNDM  =16)
      PARAMETER (KNSB  =5)
      PARAMETER (KNLT  =5)
      PARAMETER (PI=3.141592653)
      PARAMETER (KNSB1=KNSB+1)
C--- AREAS FOR THE ROUTINE
      CHARACTER ERR*64
      DIMENSION NDD(KNSB1),IUP(KNSB),IDN(KNSB)
     &, RE(KNDM,KNDM,KNLT),  TE(KNDM,KNDM,KNLT)
     &,SER(KNDM,KNA0,KNSB), SET(KNDM,KNA0,KNSB)
     &,RUP(KNDM,KNA0,KNSB1),RDN(KNDM,KNA0,KNSB1)
C--- WORKING AREAS
      DIMENSION R1D(KNDM,KNDM),R1U(KNDM,KNDM),T1D(KNDM,KNDM)
     &,T1U(KNDM,KNDM),S1D(KNDM,KNA0),S1U(KNDM,KNA0)
     &,R2D(KNDM,KNDM),T2D(KNDM,KNDM),T2U(KNDM,KNDM)
     &,S2D(KNDM,KNA0),S2U(KNDM,KNA0),TU(KNDM,KNDM),RD(KNDM,KNDM)
     &,SD(KNDM,KNA0),SU(KNDM,KNA0),RL(KNDM,KNDM,KNSB),SL(KNDM,KNA0,KNSB)
     &,RT(KNDM,KNDM),ST(KNDM,KNA0)
C--- UPWARD ADDING
      NSB1=NSB+1
      ID=IDN(NSB)
      N2=NDD(NSB)
      N3=NDD(NSB1)
      CALL RP33(RL, RE,N2, N2,NSB, ID,KNDM,KNDM,KNSB,KNDM,KNDM,KNLT)
      CALL RP33(SL,SER,N2,NA0,NSB,NSB,KNDM,KNA0,KNSB,KNDM,KNA0,KNSB)
      IF(NSB.GE.2) THEN
        DO 52 L=NSB-1,1,-1
        LB=L
        L1=L+1
        N1=NDD(L)
        IU=IUP(L)
        ID=IDN(L)
        N2=NDD(L1)
        CALL EQ23(R1D, RE,N1, N1,ID,KNDM,KNDM,KNDM,KNDM,KNLT)
        CALL EQ23(R1U, RE,N2, N2,IU,KNDM,KNDM,KNDM,KNDM,KNLT)
        CALL EQ23(T1D, TE,N2, N1,ID,KNDM,KNDM,KNDM,KNDM,KNLT)
        CALL EQ23(T1U, TE,N1, N2,IU,KNDM,KNDM,KNDM,KNDM,KNLT)
        CALL EQ23(S1D,SET,N2,NA0,LB,KNDM,KNA0,KNDM,KNA0,KNSB)
        CALL EQ23(S1U,SER,N1,NA0,LB,KNDM,KNA0,KNDM,KNA0,KNSB)
        CALL EQ23(R2D, RL,N2, N2,L1,KNDM,KNDM,KNDM,KNDM,KNSB)
        CALL EQ23(S2U, SL,N2,NA0,L1,KNDM,KNA0,KNDM,KNA0,KNSB)
        CALL ADD(0,NA0,N1,N2,N3,R1D,R1U,T1D,T1U,S1D,S1U
     &   ,R2D,T2D,T2U,S2D,S2U,TU,RD,SD,SU,ERR)
        IF(ERR.NE.' ') THEN
          ERR='ERROR IN ADD FOR UPWARD ADDING (TRN1)'
          RETURN
        ENDIF
        CALL EQ32( RL, RD,N1, N1,LB,KNDM,KNDM,KNSB,KNDM,KNDM)
        CALL EQ32( SL, SU,N1,NA0,LB,KNDM,KNA0,KNSB,KNDM,KNA0)
   52   CONTINUE
      ENDIF
C FIELD AT THE TOP OF THE SYSTEM
      N1=NDD(1)
      CALL RP33(RUP, SL,N1,NA0,1,1,KNDM,KNA0,KNSB1,KNDM,KNA0,KNSB)
      CALL RP30(RDN,0.0,N1,NA0,1,  KNDM,KNA0,KNSB1)
C--- DOWNWARD ADDING
      IU=IUP(1)
      N1=NDD(1)
      N2=NDD(2)
      CALL EQ23(RT, RE,N2, N2,IU,KNDM,KNDM,KNDM,KNDM,KNLT)
      CALL EQ23(ST,SET,N2,NA0, 1,KNDM,KNA0,KNDM,KNA0,KNSB)
      IF(NSB.GE.2) THEN
        DO 26 L=2,NSB
        LB=L
        N2=NDD(L)
C INTERNAL FIELD
        CALL EQ23(S2U,SL,N2,NA0,LB,KNDM,KNA0,KNDM,KNA0,KNSB)
        CALL EQ23(R2D,RL,N2, N2,LB,KNDM,KNDM,KNDM,KNDM,KNSB)
        CALL AXB(SD,RT,S2U,N2,N2,NA0,KNDM,KNDM,KNDM)
        CALL AAPB(SD,ST,N2,NA0,KNDM,KNDM)
        CALL AXB(RD,RT,R2D,N2,N2,N2,KNDM,KNDM,KNDM)
        CALL MULTI(N2,RD,R1D,ERR)
        IF(ERR.NE.' ') THEN
          ERR='ERROR IN MULTI FOR INTERNAL FIELD (TRN1)'
          RETURN
        ENDIF
        CALL AXB(S1D,R1D,SD,N2,N2,NA0,KNDM,KNDM,KNDM)
        CALL AXB(S1U,R2D,S1D,N2,N2,NA0,KNDM,KNDM,KNDM)
        CALL AAPB(S1U,S2U,N2,NA0,KNDM,KNDM)
        CALL EQ32(RDN,S1D,N2,NA0,LB,KNDM,KNA0,KNSB1,KNDM,KNA0)
        CALL EQ32(RUP,S1U,N2,NA0,LB,KNDM,KNA0,KNSB1,KNDM,KNA0)
C UPSIDE DOWN DIRECTION FOR APPLICATION OF THE ROUTINE ADD.
        IU=IUP(L)
        ID=IDN(L)
        N3=NDD(L+1)
        CALL EQ23(R1D,RE,N3,N3,IU,KNDM,KNDM,KNDM,KNDM,KNLT)
        CALL EQ23(R1U,RE,N2,N2,ID,KNDM,KNDM,KNDM,KNDM,KNLT)
        CALL EQ23(T1D,TE,N2,N3,IU,KNDM,KNDM,KNDM,KNDM,KNLT)
        CALL EQ23(T1U,TE,N3,N2,ID,KNDM,KNDM,KNDM,KNDM,KNLT)
        CALL EQ22(R2D,RT,N2,N2,KNDM,KNDM,KNDM,KNDM)
        CALL EQ23(S1D,SER,N2,NA0,LB,KNDM,KNA0,KNDM,KNA0,KNSB)
        CALL EQ23(S1U,SET,N3,NA0,LB,KNDM,KNA0,KNDM,KNA0,KNSB)
C Bug9601       CALL EQ22(S2U,ST,N2,NA0,KNDM,KNDM,KNDM,KNDM)
        CALL EQ22(S2U,ST,N2,NA0,KNDM,KNA0,KNDM,KNA0)
        CALL ADD(0,NA0,N3,N2,N1,R1D,R1U,T1D,T1U,S1D,S1U
     &  ,R2D,T2D,T2U,S2D,S2U,TU,RT,SD,ST,ERR)
        IF(ERR.NE.' ') THEN
          ERR='ERROR IN ADD FOR DOWNWARD ADDING (TRN1)'
          RETURN
        ENDIF
   26   CONTINUE
      ENDIF
C FIELD AT THE BOTTOM OF THE SYSTEM.
      CALL EQ32(RDN,ST,N3,NA0,NSB1,KNDM,KNA0,KNSB1,KNDM,KNA0)
      CALL RP30(RUP,0.0,N3,NA0,NSB1,KNDM,KNA0,KNSB1)
C--- DOWNWARD ADDING END
      RETURN
      END
      SUBROUTINE AXB(C,A,B,NI,NK,NJ,NCI,NAI,NBI)
C  C = A*B
C--- HISTORY
C 88. 6. 6   REGISTERED BY T. NAKAJIMA
C--- INPUT
C A      R(NAI,*)     2-DIM ARRAY  A.
C B      R(NBI,*)     2-DIM ARRAY  B.
C NI, NK, NJ  I       C(I,J) = A(I,K)*B(K,J)
C                     I=1,NI; K=1,NK; J=1,NJ
C NCI      I          SIZE FOR C(NCI,*)
C NAI      I          SIZE FOR A(NAI,*)
C NBI      I          SIZE FOR B(NBI,*)
C--- OUTPUT
C C      R(NCI,*)     A*B.
C$ENDI
      DIMENSION A(NAI,NK),B(NBI,NJ),C(NCI,NJ)
      DO 1 I=1,NI
      DO 1 J=1,NJ
      S=0.0
      DO 2 K=1,NK
    2 S=S+A(I,K)*B(K,J)
    1 C(I,J)=S
      RETURN
      END
      SUBROUTINE    TNVSS2(N,A,DT,E,NN,IW,ERR)
C     INVERSION OF REAL MATRIX. SWEEP OUT, COMPLETE POSITIONING.
C--- HISTORY
C 71. 4.30 CREATED BY SAKATA MASATO
C 89.11.10 ADDED ERR
C 90. 1. 6 ERR*(*)
C--- INPUT
C N       I        DIMENSION OF THE MATRIX
C A     R(NN,N)    MATRIX
C NN      I        SIZE OF FIRST ARGUMENT
C E       R        CONVERGENCE CRITERION (0 IS OK)
C--- OUTPUT
C DT      R        DETERMINATION OF THE MATRIX
C IW    I(2*N)     WORKING AREA
C ERR    C*64      ERROR INDICATER. IF ' ' then no error.
C$ENDI
      CHARACTER ERR*(*)
      DIMENSION     A(NN,N)  ,IW(*)
      ERR=' '
      IF(N-1)    910,930,101
  101 IF(N.GT.NN)    GO TO  900
      EPS=0.0
      DT=1.0
      DO  100     K=1,N
      PIV=0.0
      DO  110       I=K,N
      DO  110       J=K,N
      IF(ABS(A(I,J)).LE.ABS(PIV))   GO TO  110
      IPIV=I
      JPIV=J
      PIV=A(I,J)
  110 CONTINUE
      DT=DT*PIV
      IF(ABS(PIV).LE.EPS)  GO TO 920
      IF(K.EQ.1)   EPS=ABS(PIV)*E
      IF(IPIV.EQ.K)      GO TO 130
      DT=-DT
      DO 120   J=1,N
      WORK=A(IPIV,J)
      A(IPIV,J)=A(K,J)
  120 A(K,J)=WORK
  130 IF(JPIV.EQ.K)      GO TO  150
      DT=-DT
      DO 140   I=1,N
      WORK=A(I,JPIV)
      A(I,JPIV)=A(I,K)
  140 A(I,K)=WORK
  150 IW(2*K-1)=IPIV
      AA=1.0/PIV
      IW(2*K)=JPIV
      DO 210   J=1,N
  210 A(K,J)=A(K,J)*AA
      DO 220  I=1,N
      IF(I.EQ.K)    GO TO  220
      AZ=A(I,K)
      IF(AZ.EQ.0.0)   GO TO  220
      DO 230   J=1,N
  230 A(I,J)=A(I,J)-A(K,J)*AZ
      A(I,K)=-AA*AZ
  220 CONTINUE
  100 A(K,K)=AA
      DO  400 KK=2,N
      K=N+1-KK
      IJ=IW(2*K)
      IF(IJ.EQ.K)   GO TO  420
      DO 410   J=1,N
      WORK=A(IJ,J)
      A(IJ,J)=A(K,J)
  410 A(K,J)=WORK
  420 IJ=IW(2*K-1)
      IF(IJ.EQ.K)   GO TO  400
      DO 430   I=1,N
      WORK=A(I,IJ)
      A(I,IJ)=A(I,K)
  430 A(I,K)=WORK
  400 CONTINUE
      RETURN
  910 ERR='ERROR IN TINVSS: N.LE.0'
      RETURN
  900 ERR='ERROR IN TINVSS: N.GT.NN'
      RETURN
  920 DT=0.0
      INDER=N-K+1
      NNN=K-1
      WRITE(ERR,1) NNN
    1 FORMAT('TINVSS: ILL CONDITIONED MATRIX WITH RANK ',I5)
      RETURN
  930 DT=A(1,1)
      K=1
      IF(DT.EQ.0.0)     GO TO  920
      A(1,1)=1.0/A(1,1)
      RETURN
      END
      SUBROUTINE SYMTRX(ROWIJ,M,ROOT,EIGV,NI,WK,IER)
C SOLVES EIGENFUNCTION PROBLEM FOR SYMMETRIC MATRIX
C--- HISTORY
C 89.12. 4 MODIFIED WITH CNCPU
C 90. 1. 6 CNCPU IS REPLACED BY PREC.
C--- INPUT
C M       I      ORDER OF ORIGINAL SYMMETRIC MATRIX
C NI      I      INITIAL DIMENSION OF -ROOT-, -EIGV- AND -WK-
C ROWIJ  R(*)    SYMMETRIC STORAGE MODE OF ORDER M*(M+1)/2
C--- OUTPUT
C EIGV   R(NI,M) EIGENVECTORS OF ORIGINAL SYMMETRIC MATRIX
C IER     I      INDEX FOR ROOT(J) FAILED TO CONVERGE (J=IER-128)
C ROOT   R(M)    EIGENVALUES OF ORIGINAL SYMMETRIC MATRIX
C ROWIJ          STORAGE OF HOUSEHOLDER REDUCTION ELEMENTS
C WK             WORK AREA
C$ENDI
      REAL ROWIJ(*),ROOT(*),WK(*),EIGV(NI,*),CCP(3)
C+++ ADD EPSCP
C     DATA  RDELP/ 1.1921E-07 /
      CALL CPCON(CCP)
      RDELP=CCP(1)*10
C+++
      IER = 0
      MP1 = M + 1
      MM = (M*MP1)/2 - 1
      MBEG = MM + 1- M
C
C+---------------------------------------------------------------------+
C|          LOOP-100 REDUCE -ROWIJ- (SYMMETRIC STORAGE MODE) TO A      |
C|          SYMMETRIC TRIDIAGONAL FORM BY HOUSEHOLDER METHOD           |
C|                      CF. WILKINSON, J.H., 1968,                     |
C|              THE ALGEBRAIC EIGENVALUE PROBLEM, PP 290-293.          |
C|          LOOP-30&40 AND 50 FORM ELEMENT OF A*U AND ELEMENT P        |
C+---------------------------------------------------------------------+
      DO 100 II=1,M
      I = MP1 - II
      L = I - 1
      H = 0.0
      SCALE = 0.0
      IF (L.LT.1) THEN
C|          SCALE ROW (ALGOL TOL THEN NOT NEEDED)
      WK(I) = 0.0
      GO TO 90
      END IF
      MK = MM
      DO 10 K=1,L
      SCALE = SCALE + ABS(ROWIJ(MK))
      MK = MK - 1
   10 CONTINUE
      IF (SCALE.LE.0.0) THEN
      WK(I) = 0.0
      GO TO 90
      END IF
C
      MK = MM
      DO 20 K = 1,L
      ROWIJ(MK) = ROWIJ(MK)/SCALE
      H = H + ROWIJ(MK)*ROWIJ(MK)
      MK = MK - 1
   20 CONTINUE
      WK(I) = SCALE*SCALE*H
      F = ROWIJ(MM)
      G = - SIGN(SQRT(H),F)
      WK(I) = SCALE*G
      H = H - F*G
      ROWIJ(MM) = F - G
      IF (L.GT.1) THEN
      F = 0.0
      JK1 = 1
      DO 50 J=1,L
      G = 0.0
      IK = MBEG + 1
      JK = JK1
      DO 30 K=1,J
      G = G + ROWIJ(JK)*ROWIJ(IK)
      JK = JK + 1
      IK = IK + 1
   30 CONTINUE
      JP1 = J + 1
      IF (L.GE.JP1) THEN
      JK = JK + J - 1
      DO 40 K=JP1,L
      G = G + ROWIJ(JK)*ROWIJ(IK)
      JK = JK + K
      IK = IK + 1
   40 CONTINUE
      END IF
      WK(J) = G/H
      F = F + WK(J)*ROWIJ(MBEG+J)
      JK1 = JK1 + J
   50 CONTINUE
      HH = F/(H+H)
C
      JK = 1
      DO 70 J=1,L
      F = ROWIJ(MBEG+J)
      G = WK(J) - HH*F
      WK(J) = G
      DO 60 K=1,J
      ROWIJ(JK) = ROWIJ(JK) - F*WK(K) - G*ROWIJ(MBEG+K)
      JK = JK + 1
   60 CONTINUE
   70 CONTINUE
      END IF
C
      DO 80 K=1,L
      ROWIJ(MBEG+K) = SCALE*ROWIJ(MBEG+K)
   80 CONTINUE
   90 ROOT(I) = ROWIJ(MBEG+I)
      ROWIJ(MBEG+I) = H*SCALE*SCALE
      MBEG = MBEG - I + 1
      MM = MM - I
  100 CONTINUE
C
C+---------------------------------------------------------------------+
C|          LOOP-210 COMPUTE EIGENVALUES AND EIGENVECTORS              |
C|          SETUP WORK AREA LOCATION EIGV TO THE IDENTITY MATRIX       |
C|          LOOP-140 FOR FINDING SMALL SUB-DIAGONAL ELEMENT            |
C|          LOOP-160 FOR CONVERGENCE OF EIGENVALUE J (MAX. 30 TIMES)   |
C|          LOOP-190 FOR QL TRANSFORMATION AND LOOP-180 FORM VECTORS   |
C+---------------------------------------------------------------------+
      DO 110 I=1,M-1
  110 WK(I) = WK(I+1)
      WK(M) = 0.0
      B = 0.0
      F = 0.0
      DO 130 I=1,M
      DO 120 J=1,M
  120 EIGV(I,J) = 0.0
      EIGV(I,I) = 1.0
  130 CONTINUE
C
      DO 210 L=1,M
      J = 0
      H = RDELP*(ABS(ROOT(L))+ABS(WK(L)))
      IF (B.LT.H) B = H
      DO 140 N=L,M
      K = N
      IF (ABS(WK(K)).LE.B) GO TO 150
  140 CONTINUE
  150 N = K
      IF (N.EQ.L) GO TO 200
C
  160 CONTINUE
      IF (J.EQ.30) THEN
      IER = 128 + L
      RETURN
      END IF
C
      J = J + 1
      L1 = L + 1
      G = ROOT(L)
      P = (ROOT(L1)-G)/(WK(L)+WK(L))
      R = ABS(P)
      IF (RDELP*ABS(P).LT.1.0) R = SQRT(P*P+1.0)
      ROOT(L) = WK(L)/(P+SIGN(R,P))
      H = G - ROOT(L)
      DO 170 I=L1,M
      ROOT(I) = ROOT(I) - H
  170 CONTINUE
      F = F + H
C
      P = ROOT(N)
      C = 1.0
      S = 0.0
      NN1 = N - 1
      NN1PL = NN1 + L
      IF (L.LE.NN1) THEN
      DO 190 II=L,NN1
      I = NN1PL - II
      G = C*WK(I)
      H = C*P
      IF (ABS(P).LT.ABS(WK(I))) THEN
      C = P/WK(I)
      R = SQRT(C*C+1.0)
      WK(I+1) = S*WK(I)*R
      S = 1.0/R
      C = C*S
      ELSE
      C = WK(I)/P
      R = SQRT(C*C+1.0)
      WK(I+1) = S*P*R
      S = C/R
      C = 1.0/R
      END IF
      P = C*ROOT(I) - S*G
      ROOT(I+1) = H + S*(C*G+S*ROOT(I))
      IF (NI.GE.M) THEN
      DO 180 K=1,M
      H = EIGV(K,I+1)
      EIGV(K,I+1) = S*EIGV(K,I) + C*H
      EIGV(K,I) = C*EIGV(K,I) - S*H
  180 CONTINUE
      END IF
  190 CONTINUE
      END IF
      WK(L) = S*P
      ROOT(L) = C*P
      IF (ABS(WK(L)).GT.B) GO TO 160
  200 ROOT(L) = ROOT(L) + F
  210 CONTINUE
C
C+---------------------------------------------------------------------+
C|          BACK TRANSFORM EIGENVECTORS OF THE ORIGINAL MATRIX FROM    |
C|          EIGENVECTORS 1 TO M OF THE SYMMETRIC TRIDIAGONAL MATRIX    |
C+---------------------------------------------------------------------+
      DO 250 I=2,M
      L = I - 1
      IA = (I*L)/2
      IF (ABS(ROWIJ(IA+I)).GT.0.0) THEN
      DO 240 J=1,M
      SUM = 0.0
      DO 220 K=1,L
      SUM = SUM + ROWIJ(IA+K)*EIGV(K,J)
  220 CONTINUE
      SUM = SUM/ROWIJ(IA+I)
      DO 230 K=1,L
      EIGV(K,J) = EIGV(K,J) - SUM*ROWIJ(IA+K)
  230 CONTINUE
  240 CONTINUE
      END IF
  250 CONTINUE
C
      RETURN
      END
      SUBROUTINE EQ20(A,B,NI,NJ,NA1,NA2)
C A(*,*) = B
C--- HISTORY
C B        R      SOURCE SCALER.
C NI       I      A(I,J)= B,    I=1,NI;   J=1,NJ
C NJ       I
C NA1      I      SIZE FOR A(NA1,NA2)
C NA2      I
C--- OUTPUT
C A   R(NA1,NA2)  DESTINATION 2-DIM ARRAY   A(*,*) = B
C$ENDI
      DIMENSION A(NA1,NA2)
      DO 1 J=1,NJ
      DO 1 I=1,NI
    1 A(I,J)=B
      RETURN
      END
      SUBROUTINE AAPB(A,B,NI,NJ,NAI,NBI)
C A=A+B
C--- HISTORY
C 88. 6. 6  REGISTERED BY T.NAKAJIMA
C--- INPUT
C A    R(NAI,*)    2-DIM ARRAY  A.
C B    R(NBI,*)    2-DIM ARRAY  B.
C NI      I        SUM (I,J) I=1,NI
C NJ      I        SUM (I,J) J=1,NI
C--- OUTPUT
C A                A+B
C$ENDI
      DIMENSION A(NAI,NJ),B(NBI,NJ)
      DO 1 J=1,NJ
      DO 1 I=1,NI
    1 A(I,J)=A(I,J)+B(I,J)
      RETURN
      END
      SUBROUTINE ADD(IT,NA0,N1,N2,N3,R1D,R1U,T1D,T1U,S1D,S1U
     & ,R2D,T2D,T2U,S2D,S2U,TU,RD,SD,SU,ERR)
C ADDING TWO LAYERS 1 AND 2.
C--- HISTORY
C 86.10.15  CHECK OK.
C 89.10.30  RE-EDIT.
C 90. 2. 2  ELEMINATE KND00
C--- INPUT
C IT         I      INDICATOR FOR CALCULATION OF TU AND SD.
C NA0        I      NO. OF SOLAR DIRECTIONS.
C N1,N2,N3   I       ----------------------------     MU1(I), I=1,N1
C R1D, R1U R(KNDM,   R1D, R1U, T1D, T1U, S1D, S1U    (LAYER-1)
C T1D, T1U   KNDM)   ----------------------------     MU2(I), I=1,N2
C R2D, R2U           R2D, R2U, T2D, T2U, S2D, S2U    (LAYER-2)
C T2D, T2U           ----------------------------     MU3(I),I=1,N3
C S1D, S1U R(KNDM,   SUFFIX U = UPGOING,   D = DOWNGOING INCIDENCES.
C S2D, S2U   KNA0)          R = REFLECTION,T = TRANSMISSION MATRICES.
C                           S = SOURCE MATRIX.
C--- OUTPUT
C RD       R(KNDM,   -----------------    MU1(I), I=1,N1
C TU         KNDM)   RD, TU, SD, SU       (LAYER 1+2)
C SD       R(KNDM,   -----------------    MU3(I), I=1,N3
C SU         KNA0)
C ERR       C*64     ERROR INDEX.
C--- PARAMETER
C KNA0       I       NUMBER OF SOLAR ZENITH ANGLES.
C KNDM       I       NUMBER OF QUADRATURE POINTS.
C--- AREAS FOR THIS ROUTINE
      PARAMETER (KNA0  =1)
      PARAMETER (KNDM  =16)
C
      CHARACTER ERR*64
      DIMENSION R1D(KNDM,KNDM),R1U(KNDM,KNDM),T1D(KNDM,KNDM)
     &,T1U(KNDM,KNDM),S1D(KNDM,KNA0),S1U(KNDM,KNA0)
     &,R2D(KNDM,KNDM),T2D(KNDM,KNDM),T2U(KNDM,KNDM)
     &,S2D(KNDM,KNA0),S2U(KNDM,KNA0)
     &,TU(KNDM,KNDM),RD(KNDM,KNDM),SD(KNDM,KNA0),SU(KNDM,KNA0)
C--- WORKING AREAS
      DIMENSION AA(KNDM,KNA0),BB(KNDM,KNDM),CC(KNDM,KNDM),DD(KNDM,KNDM)
C
      CALL AXB(AA,R2D,S1D,N2,N2,NA0,KNDM,KNDM,KNDM)
      CALL APB(SU,AA,S2U,N2,NA0,KNDM,KNDM,KNDM)
      CALL AXB(CC,R2D,R1U,N2,N2,N2,KNDM,KNDM,KNDM)
      CALL MULTI(N2,CC,BB,ERR)
      IF(ERR.NE.' ') THEN
        ERR='ERROR IN MULTI OF ADD'
        RETURN
      ENDIF
      CALL AXB(SD,BB,SU,N2,N2,NA0,KNDM,KNDM,KNDM)
      CALL AXB(SU,T1U,SD,N1,N2,NA0,KNDM,KNDM,KNDM)
      CALL AAPB(SU,S1U,N1,NA0,KNDM,KNDM)
      CALL AXB(CC,T1U,BB,N1,N2,N2,KNDM,KNDM,KNDM)
      CALL AXB(DD,R2D,T1D,N2,N2,N1,KNDM,KNDM,KNDM)
      CALL AXB(BB,CC,DD,N1,N2,N1,KNDM,KNDM,KNDM)
      CALL APB(RD,R1D,BB,N1,N1,KNDM,KNDM,KNDM)
      IF(IT.LE.0) RETURN
      CALL AXB(TU,CC,T2U,N1,N2,N3,KNDM,KNDM,KNDM)
      CALL AXB(AA,R1U,SD,N2,N2,NA0,KNDM,KNDM,KNDM)
      CALL AAPB(AA,S1D,N2,NA0,KNDM,KNDM)
      CALL AXB(SD,T2D,AA,N3,N2,NA0,KNDM,KNDM,KNDM)
      CALL AAPB(SD,S2D,N3,NA0,KNDM,KNDM)
      RETURN
      END
      SUBROUTINE EQ22(A,B,NI,NJ,NA1,NA2,NB1,NB2)
C  A = B
C--- HISTORY
C 88. 6. 6  REGISTERED BY T. NAKAJIMA
C--- INPUT
C B      R(NB1,NB2)     SOURCE 2-DIM ARRAY  B.
C NI       I           A(I,J) = B(I,J),  I=1,NI;  J=1,NJ
C NJ       I
C NA1      I           DIM  A(NA1,NA2)
C NA2      I
C NB1      I           DIM  B(NB1,NB2)
C--- OUTPUT
C A     R(NA1,NA2)     DESTINATION 2-DIM ARRAY  A = B.
C$ENDI
      DIMENSION A(NA1,NA2),B(NB1,NB2)
      DO 1 J=1,NJ
      DO 1 I=1,NI
    1 A(I,J)=B(I,J)
      RETURN
      END
      SUBROUTINE EQ23(A,B,NI,NJ,K,NA1,NA2,NB1,NB2,NB3)
C  A(*,*)= B(*,*,K)
C--- HISTORY
C 88. 6. 6  REGISTERED BY T. NAKAJIMA
C--- INPUT
C B      R(NB1,NB2,NB3)     SOURCE 3-DIM ARRAY B.
C NI         I              A(I,J)= B(I,J,K)
C NJ         I              I=1,NI;  J=1,NJ
C K          I
C NA1        I              DIM A(NA1, NA2)
C NA2        I
C NB1        I              DIM B(NB1,NB2,NB3)
C NB2, NB3   I
C--- OUTPUT
C A     R(NA1,NA2)          DESTINATION 2-DIM ARRAY A(*,*) = B(*,*,K)
C$ENDI
      DIMENSION A(NA1,NA2),B(NB1,NB2,NB3)
      DO 1 J=1,NJ
      DO 1 I=1,NI
    1 A(I,J)=B(I,J,K)
      RETURN
      END
      SUBROUTINE MULTI(ND,CD,CC,ERR)
C CALCULATION OF MULTIPLE REFLECTION BETWEEN TWO LAYERS.
C--- HISTORY
C 86.10.15  CHECK OK.
C 89.10.30  RE-EDIT.
C 94. 5. 7  ERR*64 -> ERR*(*)
C--- INPUT
C ND          I     NO. OF STREAMS AT THE INTERFACE BETWEEN TWO LAYERS.
C CD      R(KNDM,   R = R1 * R2.
C           KNDM)
C--- OUTPUT
C CC      R(KNDM,   CC = ( 1 - CD )**-1 = 1 + CD + CD**2 + CD**3 + ...
C           KNDM)
C ERR     C*64      ERROR INDEX.
C--- PARAMETER
C KNDM        I     NUMBER OF QUADRATURE POINTS.
C--- AREAS FOR THIS ROUTINE
      PARAMETER (KNDM  =16)
C
      CHARACTER ERR*(*)
      DIMENSION CD(KNDM,KNDM),CC(KNDM,KNDM)
C--- WORKING AREA
      PARAMETER (KNDM2=2*KNDM)
      DIMENSION IW(KNDM2)
C---
      DO 1 I=1,ND
      DO 2 J=1,ND
    2 CC(I,J)= -CD(I,J)
    1 CC(I,I)=1 - CD(I,I)
C INVERSION OF -CC-.
      EPS=0
      CALL TNVSS2(ND,CC,DT,EPS,KNDM,IW,ERR)
      IF(ERR.NE.' ') ERR='ERROR IN MULTI'
      RETURN
      END
      SUBROUTINE RP33(A,B,NI,NJ,K,L,NA1,NA2,NA3,NB1,NB2,NB3)
C A(*,*,K) = B(*,*,L)
C--- HISTORY
C 88. 5. 6  REGISTERED BY T. NAKAJIMA
C--- INPUT
C B      R(NB1,NB2,NB3)     SOURCE 3-DIM ARRAY  B.
C NI,NJ,K,L    I            A(I,J,K)=B(I,J,L), I=1,NI; J=1,NJ
C NA1,NA2,NA3  I            DIM A(NA1,NA2,NA3)
C NB1,NB2,NB3  I            DIM B(NB1,NB2,NB3)
C--- OUTPUT
C A      R(NA1,NA2,NA3)     DESTINATION 3-DIM ARRAY A(*,*,K)=B(*,*,L)
C$ENDI
      DIMENSION A(NA1,NA2,NA3),B(NB1,NB2,NB3)
      DO 1 J=1,NJ
      DO 1 I=1,NI
    1 A(I,J,K)=B(I,J,L)
      RETURN
      END
      SUBROUTINE RP30(A,B,NI,NJ,K,NA1,NA2,NA3)
C A(*,*,K) = B
C--- HISTORY
C 88. 5. 6   REGISTERED BY T. NAKAJIMA
C--- INPUT
C B            R       SOURCE SCALER
C NI,NJ,K      I       A(I,J,K)=B,  I=1,NI; J=1,NJ
C NA1,NA2,NA3  I       DIM A(NA1,NA2,NA3)
C--- OUTPUT
C A    R(NA1,NA2,NA3)  DESTINATION 3-DIM ARRAY A(*,*,K)=B
C$ENDI
      DIMENSION A(NA1,NA2,NA3)
      DO 1 J=1,NJ
      DO 1 I=1,NI
    1 A(I,J,K)=B
      RETURN
      END
      SUBROUTINE APB(C,A,B,NI,NJ,NCI,NAI,NBI)
C C=A+B
C--- HISTORY
C 88. 6. 6  REGISTERED BY T. NAKAJIMA
C--- INPUT
C  A       R(NAI,*)     2-DIM ARRAY A.
C  B       R(NBI,*)     2-DIM ARRAY B.
C NI         I          SUM (I,J) , I=1,NI
C NJ         I          SUM (I,J) , J=1,NJ
C NCI        I          SIZE FOR C(NCI,*)
C NAI        I          SIZE FOR A(NAI,*)
C NBI        I          SIZE FOR B(NBI,*)
C--- OUTPUT
C  C       R(NCI,*)     A+B
C$ENDI
      DIMENSION A(NAI,NJ),B(NBI,NJ),C(NCI,NJ)
      DO 1 J=1,NJ
      DO 1 I=1,NI
    1 C(I,J)=A(I,J)+B(I,J)
      RETURN
      END
