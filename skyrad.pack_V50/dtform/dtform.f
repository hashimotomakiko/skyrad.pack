C$OUTF dtform.f
C$LIBF /Users/makibo/Fortran/LIB/SKYLIB/A
C$LIBF /Users/makibo/Fortran/LIB/POMLIB/A
C$LIBF /Users/makibo/Fortran/LIB/CCDLIB/A
C$LIBF /Users/makibo/Fortran/LIB/ANGLIB/A
C$LIBF /Users/makibo/Fortran/LIB/MATLIB/A
C$LIBF /Users/makibo/Fortran/LIB/BASLIB/A
C$LIBF /Users/makibo/Fortran/LIB_Nakajima/LBR1/A
C$LIBF /Users/makibo/Fortran/LIB_Nakajima/LBM1/A
C
C --- Note
C   Skyradiometer data processing program
C   For PREDE instruments POM-01L/POM-01MKII
C
C  (1) Datafile Conversions for skyradiance analyses
C
C --- history
C  2001.02.14 Created by M.Yamano
C  2001.04.17 Modified to accept new file name(yyyymmdd.dat).
C  2001.07.13 Added output tag file(Tag/yymmdd.tag).
C  2001.08.03 Changed log file to single file(dtform.log).
C  2001.08.22 Time in tag file is always LT.
C  2001.09.05 Changed tag file format.
C  2001.10.08 Renewed.
C  2002.05.19 Renewed by M.Yamano for Ver.5.
C  2002.11.23 Modified.
C  2002.12.04 Renewed.
C  2003.03.06 Added procedure for ITYP=22.
C  2003.06.03 FG0 priority: FGN > FGN0 > FGP> FGP0 > FGF > FGF0 (sub AGCRCT)
C  2003.06.04 DPMX=10/ DFMX=2/ DNMX=0.5 (sub AGCRCT)
C  2003.11.18 Added processing for ITYP=30 (new format) data.
C  2004.01.14 Debugged.
C  2004.01.19 direct data only(NA=1) are available in sub DTLAND & DTLND2.
C  2004.12.06 Start on Mac.
C  2006.03.14 Changed control for S/N check in case of ITYP=30
C  2024.01.17 Edited for Ver5.0 by M.Hashimoto
C
C --- I/O files
C  dtform.par          I    processing control file
C  obs.para            I    observation parameters file(any name is OK)
C  ins.para            I    instrument parameters file(any name is OK)
C  CCD.para            I    CCD parameters file for POM-01MKII(any name is OK)
C  fname               I    processing dates file
C  DAT/yymmddnn.dat    I    measurement data file (ITYP=30)
C  or  yymmdd.DAT           measurement data file (ITYP<30)
C  or  yyyymmdd.dat
C  DT3/(yy)yymmdd.DT3  O    input file for analysis ver.3(if IDT3>0)
C  DT4/(yy)yymmdd.DT4  O    input file for analysis ver.4(if IDT4>0)
C  DT5/(yy)yymmdd.DT5  O    input file for analysis ver.5(if IDT5>0)
C  Lst/(yy)yymmdd.lst  O    measurement data list file(if ILST>0)
C  Tag/(yy)yymmdd.tag  O    tag file for available measurements
C  dtform.tag          O    merged tag file of all Tag/*.tag files
C  dtform.log          O    processing log file
C  (yy)yymmdd.fg       O    calculated Fg angle log file(if IDBG>0)
C  (yy)yymmdd.fi       O    calculated azimuth angle log file(if IDBG>0)
C  (yy)yymmdd.ht       O    calculated height angle log file(if IDBG>0)
C  (yy)yymmdd.fgM      O    calculated average Fg angle log file(if IDBG>0)
C  dtform.tmp          W    work file
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KNA   =50)
      PARAMETER (KDAY  =30)
C
      CHARACTER ERC*64
      CHARACTER ALINE*256,ALIN0*256,FNM*80,MRK*1
      CHARACTER DFIL*80,OFL1*80,OFL2*80,OFL3*80,LFIL*80,TFIL*80
     &         ,WFIL*80,WFL1*80,WFL2*80,WFL3*80,WFL4*80
C
      CHARACTER CMNT*40,SITE*20,CNTY*20,SRNO*20,SRNOO*20,SRNOD*7
      DIMENSION WLO(KNW),WL(KNW),SOLID(KNW),SA(KNA),THM(KNA),FIM(KNA)
     &         ,TH(KNW,KNA),FI(KNW,KNA),U(KNW,KNA),U2(KNW,KNA)
      DIMENSION FC(5)
C
      DATA SRNOD /'0000000'/
C
      ERC=' '
      IUI=1
      OPEN (IUI,FILE='dtform.par',STATUS='OLD')
      CALL RDDTFP(IUI,IDT3,IDT4,IDT5,ILST,IDBG,IANG
     &               ,CMNT,SITE,CNTY,ALNGSO,ALNGO,ALATO,ALT,NA0,SA
     &               ,SRNOO,ITYP,NWO,WLO,SOLID
     &               ,PT0,RL0,DHT,DFAI
     &               ,SH0,SV0,ASPN,CSA,DZN,FH0,FV0,ASPF,CO,DFI,FC,ERC)
      CLOSE (IUI)
      IF (ERC.NE.' ') GOTO 900
C
      IF ((ITYP.LT.20).OR.(ITYP.GE.30)) THEN
         IF (SA(1).GE.0.1) THEN
            ERC='No direct measurement !'
            GOTO 900
         ENDIF
      ENDIF
C
      IF ((ITYP.EQ.10).OR.(ITYP.GE.30)) THEN
         MRK='/'
      ELSE
         MRK=':'
      ENDIF
C
C --- processing
C
      IUO1=2
      IUO2=3
      IUO3=4
      IUL=8
      IUT=9
      IUW=10
      WFIL='dtform.tmp'
      IUTG=20
      OPEN (IUTG,FILE='dtform.tag',STATUS='UNKNOWN')
      IUR=21
      OPEN (IUR,FILE='dtform.log',STATUS='UNKNOWN')
      IUF=7
      OPEN (IUF,FILE='fname',STATUS='OLD')
      MTG=0
C
  110 READ(IUF,112,END=800,ERR=800) FNM
  112 FORMAT(A80)
      CALL GETPOS(FNM,'.',80,IL)
      IF (IL.EQ.1) GOTO 800
      IF (IL.GT.0) FNM=FNM(1:IL-1)
      CALL UTLSPC(FNM,NFNM)
C
      WRITE(IUR,115) FNM(1:NFNM)
  115 FORMAT(A8,'.log')
      TFIL='Tag/'//FNM(1:NFNM)//'.tag'
      OPEN (IUT,FILE=TFIL,STATUS='UNKNOWN')
      IF (ITYP.GE.30) THEN
         DFIL='DAT/'//FNM(1:NFNM)//'.dat'
      ELSE
         IF (NFNM.EQ.6) THEN
            DFIL='DAT/'//FNM(1:NFNM)//'.DAT'
         ELSE
            DFIL='DAT/'//FNM(1:NFNM)//'.dat'
         ENDIF
      ENDIF
      OPEN (IUI,FILE=DFIL,STATUS='OLD')
      IF (ITYP.GE.30) THEN
         CALL RDHD30(IUI,SRNO,ALNGS,ALNG,ALAT,NW,WL,ERC)
         IF (ERC.EQ.' ') THEN
            DLGS=ABS(ALNGSO-ALNGS)
            DLG=ABS(ALNGO-ALNG)
            DLT=ABS(ALATO-ALAT)
            IF ((DLGS.GE.0.5).OR.(DLG.GE.0.01).OR.(DLT.GE.0.01))
     &                     ERC='Longitude or latitude mismatch !'
            IF (NW.NE.NWO) ERC='Number of wavelength mismatch !'
            IF (SRNO.NE.SRNOO) THEN
               IF (SRNO(1:7).EQ.SRNOD) THEN
c                 WRITE(*,117) FNM(1:NFNM)
c 117             FORMAT('No check for S/N! in file (',A8,'.dat)')
               ELSE
                  ERC='S/N mismatch !'
               ENDIF
            ENDIF
            IW=0
            DO WHILE ((IW.LT.NW).AND.(ERC.EQ.' '))
              IW=IW+1
              WL0=WLO(IW)*1.E7
              IF (ABS(WL0-WL(IW)).GE.1.) ERC='Wavelength mismatch !'
            ENDDO
         ENDIF
         IF (ERC.EQ.' ') THEN
            DO IW=1,NW
              WL(IW)=WL(IW)*1.E-7
            ENDDO
         ELSE
            WRITE(*,118) FNM(1:NFNM),ERC
  118       FORMAT('Skip file (',A8,'.dat) : ',A64)
            ERC=' '
            GOTO 110
         ENDIF
      ELSE
         NW=NWO
         DO IW=1,NW
           WL(IW)=WLO(IW)
         ENDDO
         ALNGS=ALNGSO
         ALNG=ALNGO
         ALAT=ALATO
      ENDIF
C
      IF (IDT3.GT.0) THEN
         OFL1='DT3/'//FNM(1:NFNM)//'.DT3'
         OPEN (IUO1,FILE=OFL1,STATUS='UNKNOWN')
      ENDIF
      IF (IDT4.GT.0) THEN
         OFL2='DT4/'//FNM(1:NFNM)//'.DT4'
         OPEN (IUO2,FILE=OFL2,STATUS='UNKNOWN')
      ENDIF
      IF (IDT5.GT.0) THEN
         OFL3='DT5/'//FNM(1:NFNM)//'.DT5'
         OPEN (IUO3,FILE=OFL3,STATUS='UNKNOWN')
      ENDIF
      IF (ILST.GT.0) THEN
         LFIL='Lst/'//FNM(1:NFNM)//'.lst'
         OPEN (IUL,FILE=LFIL,STATUS='UNKNOWN')
      ENDIF
      IF (IDBG.GT.0) THEN
         WFL1=FNM(1:NFNM)//'.fg'
         OPEN (IUW+1,FILE=WFL1,STATUS='UNKNOWN')
         WFL2=FNM(1:NFNM)//'.fi'
         OPEN (IUW+2,FILE=WFL2,STATUS='UNKNOWN')
         WFL3=FNM(1:NFNM)//'.ht'
         OPEN (IUW+3,FILE=WFL3,STATUS='UNKNOWN')
         WFL4=FNM(1:NFNM)//'.fgM'
         OPEN (IUW+4,FILE=WFL4,STATUS='UNKNOWN')
         NG1=0
         NG2=0
      ENDIF
      NTG=0
      IFG=0
      NDT=0
C
  120 READ(IUI,121,ERR=200,END=200) ALINE
  121 FORMAT(A256)
      IF (ALINE(1:1).EQ.' ') ALINE=ALINE(2:256)
      IF (ALINE(1:10).EQ.'          ') GOTO 120
      IF (ALINE(3:3).EQ.MRK) THEN
         CLOSE(IUW)
C
         IF (IFG.GT.0) THEN
            OPEN (IUW,FILE=WFIL,STATUS='UNKNOWN')
            IF (ITYP.GE.30) THEN
               CALL DTLND2(IUW,ITYP,ALNGS,ALNG,ALAT,NW,WL,SOLID
     &                    ,NA0,NA,SA,IY,IM,ID,JH,JM,JS,TM
     &                    ,TH0,FAI,DST,THM,FIM,TH,FI,U,U2,ERC,IFLG)
            ELSEIF (ITYP.LT.20) THEN
               CALL DTLAND(IUW,ITYP,ALNGS,ALNG,ALAT,NW,WL,SOLID
     &                    ,NA0,NA,SA,IY,IM,ID,JH,JM,JS,TM
     &                    ,TH0,FAI,DST,THM,FIM,TH,FI,U,U2,ERC,IFLG)
            ELSE
               CALL DTSHIP(IUW,ITYP,IANG,IDBG,NG1,NG2,ALNGS,NW,WL,SOLID
     &                    ,IY,IM,ID,JH,JM,JS,TM,ALAT,ALNG,NA
     &                    ,TH0,FAI,DST,THM,FIM,TH,FI,U,U2
     &                    ,PT0,RL0,DHT,DFAI,SH0,SV0,ASPN,CSA,DZN
     &                    ,FH0,FV0,ASPF,CO,DFI,FC,ERC,IFLG)
            ENDIF
            CLOSE(IUW)
            WRITE(IUR,125) ERC
  125       FORMAT(A64)
            IF (IFLG.EQ.0) THEN
               NTG=NTG+1
               SAMX=SA(NA)
               CALL WTTAG0(IUT,NTG,IY,IM,ID,TM,ALNG,ALAT,TH0,SAMX)
               MTG=MTG+1
               IF (MTG.EQ.1) WRITE(IUTG,130)
  130          FORMAT('  TNo  No yyyy mm dd Hour   Long    Lat    '
     &               ,'Hs    SA(max)')
               WRITE(IUTG,140) MTG,NTG,IY,IM,ID,TM,ALNG,ALAT,90.0-TH0
     &                        ,SA(NA)
  140          FORMAT(I5,I4,I5,2I3,F6.2,F8.2,2F7.2,F6.1)
               IF (IDT3.GT.0)
     &            CALL WTDAT3(IUO1,IY,IM,ID,JH,JM,JS,TM,ALNGS,ALNG,ALAT
     &                       ,ALT,NA,TH0,FAI,DST,NW,WL,THM,FIM,U)
               IF (IDT4.GT.0)
     &            CALL WTDAT4(IUO2,IY,IM,ID,JH,JM,JS,TM,ALNGS,ALNG,ALAT
     &                       ,ALT,NA,TH0,FAI,DST,NW,WL,TH,FI,U)
               IF (IDT5.GT.0)
     &            CALL WTDAT4(IUO3,IY,IM,ID,JH,JM,JS,TM,ALNGS,ALNG,ALAT
     &                       ,ALT,NA,TH0,FAI,DST,NW,WL,TH,FI,U)
CO hashimoto20240116     &                       ,ALT,NA,TH0,FAI,DST,NW,WL,TH,FI,U2)
            ENDIF
         ENDIF
C
         NDT=NDT+1
         ALIN0=ALINE
         OPEN (IUW,FILE=WFIL,STATUS='UNKNOWN')
         WRITE(IUW,121) ALIN0
         IF (ILST.GT.0) WRITE(IUL,150) NDT,ALIN0
  150    FORMAT(I3,2X,A72)
         IFG=1
c        print*,ALIN0
      ELSE
         WRITE(IUW,121) ALINE
      ENDIF
      GOTO 120
C
  200 CLOSE(IUW)
      IF (IFG.GT.0) THEN
         OPEN (IUW,FILE=WFIL,STATUS='UNKNOWN')
         IF (ITYP.GE.30) THEN
            CALL DTLND2(IUW,ITYP,ALNGS,ALNG,ALAT,NW,WL,SOLID
     &                 ,NA0,NA,SA,IY,IM,ID,JH,JM,JS,TM
     &                 ,TH0,FAI,DST,THM,FIM,TH,FI,U,U2,ERC,IFLG)
         ELSEIF (ITYP.LT.20) THEN
            CALL DTLAND(IUW,ITYP,ALNGS,ALNG,ALAT,NW,WL,SOLID
     &                 ,NA0,NA,SA,IY,IM,ID,JH,JM,JS,TM
     &                 ,TH0,FAI,DST,THM,FIM,TH,FI,U,U2,ERC,IFLG)
         ELSE
            CALL DTSHIP(IUW,ITYP,IANG,IDBG,NG1,NG2,ALNGS,NW,WL,SOLID
     &                 ,IY,IM,ID,JH,JM,JS,TM,ALAT,ALNG,NA
     &                 ,TH0,FAI,DST,THM,FIM,TH,FI,U,U2
     &                 ,PT0,RL0,DHT,DFAI,SH0,SV0,ASPN,CSA,DZN
     &                 ,FH0,FV0,ASPF,CO,DFI,FC,ERC,IFLG)
         ENDIF
         CLOSE(IUW)
         WRITE(IUR,125) ERC
         IF (IFLG.EQ.0) THEN
            NTG=NTG+1
            SAMX=SA(NA)
            CALL WTTAG0(IUT,NTG,IY,IM,ID,TM,ALNG,ALAT,TH0,SAMX)
            MTG=MTG+1
            IF (MTG.EQ.1) WRITE(IUTG,130)
            WRITE(IUTG,140)MTG,NTG,IY,IM,ID,TM,ALNG,ALAT,90.0-TH0,SA(NA)
            IF (IDT3.GT.0)
     &         CALL WTDAT3(IUO1,IY,IM,ID,JH,JM,JS,TM,ALNGS,ALNG,ALAT
     &                    ,ALT,NA,TH0,FAI,DST,NW,WL,THM,FIM,U)
            IF (IDT4.GT.0)
     &         CALL WTDAT4(IUO2,IY,IM,ID,JH,JM,JS,TM,ALNGS,ALNG,ALAT
     &                    ,ALT,NA,TH0,FAI,DST,NW,WL,TH,FI,U)
            IF (IDT5.GT.0)
     &         CALL WTDAT4(IUO3,IY,IM,ID,JH,JM,JS,TM,ALNGS,ALNG,ALAT
     &                    ,ALT,NA,TH0,FAI,DST,NW,WL,TH,FI,U)
CO hashimoto20240116     &                    ,ALT,NA,TH0,FAI,DST,NW,WL,TH,FI,U2)
         ENDIF
      ENDIF
      CLOSE(IUI)
      CLOSE(IUL)
      CLOSE(IUO1)
      CLOSE(IUO2)
      CLOSE(IUO3)
      CLOSE(IUT)
      GOTO 110
C
  800 CLOSE (IUF)
      OPEN (IUW,FILE=WFIL,STATUS='UNKNOWN')
      CLOSE (IUW,STATUS='DELETE')
      GOTO 1000
C
  900 WRITE(*,*) 'ERROR CODE: ', ERC
 1000 STOP
      END
C
      SUBROUTINE DTLND2(IUI,ITYP,ALNGS,ALNG,ALAT,NW,WL,SOLID
     &                 ,NA0,NA,SA,IY,IM,ID,JH,JM,JS,TM
     &                 ,TH0,FAI,DST,THM,FIM,TH,FI,U,U2,ERC,IFLG)
      PARAMETER (KNW   =20)
      PARAMETER (KNA   =50)
C
      CHARACTER ERC*(*)
      CHARACTER ATIME*8,DMY*9
      DIMENSION WL(KNW),SOLID(KNW),SA(KNA),THM(KNA),FIM(KNA)
     &         ,TH(KNW,KNA),FI(KNW,KNA),U(KNW,KNA),U2(KNW,KNA)
      DIMENSION HTS(KNA)
      DIMENSION AP(KNA),HP(KNA)
C
      DT=51.0/3600.0
C
      ERC=' '
      IFLG=-9
      CALL RDDT30(IUI,NW,IY,IM,ID,JH,JM,JS,JFLG,NA,AP,HP,U,ERC)
      IF (ERC.NE.' ') GOTO 900
      IF (NA.GT.NA0) NA=NA0
C
      WRITE(DMY,10) 100+JH,100+JM,100+JS
   10 FORMAT(3I3)
      ATIME=DMY(2:3)//':'//DMY(5:6)//':'//DMY(8:9)
      ERC=ATIME//'  '
      IFLG=0
      IF (NA.LE.0) THEN
         ERC(11:25)='DTLND2: NO DATA'
         IFLG=-8
         GOTO 900
      ENDIF
C
      TM=JH+JM/60.0+JS/3600.0
      CALL SUNH(IY,IM,ID,TM,ALAT,ALNG,ALNGS,DT,DST,HGT,FAI)
      FAI=DEG(FAI,0)
      TH0=90.-HGT
C
      DAP=AP(1)-FAI
      DHP=HP(1)-HGT
c     print*, ATIME,DAP,DHP
      FIM0=AP(1)-DAP
      THM(1)=TH0
      FIM(1)=0.
      HTS(1)=HGT
      IF (NA.GT.1) THEN
        DO J=2,NA
          AP1=(AP(J)-DAP)-FIM0
          ZP1=90.-(HP(J)-DHP)
          SA1=SCAANG(AP1,ZP1,0.,TH0)
c         IF (ABS(SA1/SA(J)-1.).GT.0.05) print*,J,SA(J),SA1
          THM(J)=ZP1
          FIM(J)=AP1
          HTS(J)=HGT
        ENDDO
      ENDIF
C
      CALL INTENS(NW,SOLID,NA,HTS,U,U2,DST,IFLG)
      IF (IFLG.EQ.-2) IFLG=0
      IF (IFLG.NE.0) THEN
         IF (IFLG.EQ.-1) ERC(11:32)='DTLND2: NO DIRECT DATA'
         IF (IFLG.EQ.-2) ERC(11:33)='DTLND2: NO DIFFUSE DATA'
         GOTO 900
      ENDIF
C
      DO I=1,NW
        DO J=1,NA
          TH(I,J)=THM(J)
          FI(I,J)=FIM(J)
        ENDDO
      ENDDO
C
  900 RETURN
      END
C
      SUBROUTINE DTLAND(IUI,ITYP,ALNGS,ALNG,ALAT,NW,WL,SOLID
     &                 ,NA0,NA,SA,IY,IM,ID,JH,JM,JS,TM
     &                 ,TH0,FAI,DST,THM,FIM,TH,FI,U,U2,ERC,IFLG)
      PARAMETER (KNW   =20)
      PARAMETER (KNA   =50)
C
      CHARACTER ERC*(*)
      CHARACTER ALINE*256
      DIMENSION WL(KNW),SOLID(KNW),SA(KNA),THM(KNA),FIM(KNA)
     &         ,TH(KNW,KNA),FI(KNW,KNA),U(KNW,KNA),U2(KNW,KNA)
      DIMENSION HTS(KNA)
C
      DT=51.0/3600.0
C
      ERC=' '
      IFLG=-9
      IF (ITYP.EQ.11) THEN
         CALL RDDT11(IUI,ALINE,NW,NA,U,ERC)
         IST=1
      ELSE
         CALL RDDT10(IUI,ALINE,NW,NA,SA,U,ERC)
         IST=10
      ENDIF
      IF (ERC.NE.' ') GOTO 900
      IF (NA.GT.NA0) NA=NA0
C
      ERC=ALINE(IST:IST+7)//'  '
      IFLG=0
      IF (NA.LE.0) THEN
         ERC(11:25)='DTLAND: NO DATA'
         IFLG=-8
         GOTO 900
      ENDIF
C
      IF (ITYP.EQ.11) THEN
         CALL RDLB11(ALINE,IY,IM,ID,JH,JM,JS,JFLG)
      ELSE
         CALL RDLB10(ALINE,IY,IM,ID,JH,JM,JS,JFLG)
      ENDIF
      TM=JH+JM/60.0+JS/3600.0
      CALL SUNH(IY,IM,ID,TM,ALAT,ALNG,ALNGS,DT,DST,HGT,FAI)
      FAI=DEG(FAI,0)
      TH0=90.-HGT
C
      CALL CALANG(JFLG,TH0,NA,SA,THM,FIM,IDRCT,MXIDX)
      NA=MXIDX
      IF (IDRCT.EQ.1) THEN
         DO J=1,NA
           HTS(J)=HGT
         ENDDO
         CALL INTENS(NW,SOLID,NA,HTS,U,U2,DST,IFLG)
         IF (IFLG.EQ.-2) IFLG=0
         IF (IFLG.NE.0) THEN
            IF (IFLG.EQ.-1) ERC(11:32)='DTLAND: NO DIRECT DATA'
            IF (IFLG.EQ.-2) ERC(11:33)='DTLAND: NO DIFFUSE DATA'
            GOTO 900
         ENDIF
      ELSE
         IF (IDRCT.LT.1) ERC(11:39)='DTLAND: NO DIRECT MEASUERMENT'
         IF (IDRCT.GT.1) ERC(11:38)='DTLAND: ABNORMAL DATA FORMAT'
         IFLG=-7
         GOTO 900
      ENDIF
C
      DO I=1,NW
        DO J=1,NA
          TH(I,J)=THM(J)
          FI(I,J)=FIM(J)
        ENDDO
      ENDDO
C
  900 RETURN
      END
C
      SUBROUTINE DTSHIP(IUI,ITYP,IANG,IDBG,NG1,NG2,ALNGS,NW,WL,SOLID
     &                 ,IY,IM,ID,JH,JM,JS,TM,ALAT,ALNG,NA
     &                 ,TH0,FAI,DST,THM,FIM,TH,FI,U,U2
     &                 ,PT0,RL0,DHT,DFAI,SH0,SV0,ASPN,CSA,DZN
     &                 ,FH0,FV0,ASPF,CO,DFI,FC,ERC,IFLG)
      PARAMETER (KNW   =20)
      PARAMETER (KNA   =50)
C
      CHARACTER ERC*(*)
      CHARACTER ALINE*256
      DIMENSION WL(KNW),SOLID(KNW),SA(KNA),THM(KNA),FIM(KNA)
     &         ,TH(KNW,KNA),FI(KNW,KNA),U(KNW,KNA),U2(KNW,KNA)
      DIMENSION FC(5)
      DIMENSION AP(KNW,KNA),HP(KNW,KNA),KFLG(KNA)
     &         ,FH(KNW,KNA),FV(KNW,KNA),SH(KNW,KNA),SV(KNW,KNA)
     &         ,PT(KNW,KNA),RL(KNW,KNA)
      DIMENSION UTM(KNA),FIS(KNA),HTS(KNA),APC(KNW,KNA),HPC(KNW,KNA)
C
      DATA TISEC /12.0/
C
      DT=51.0/3600.0
C
      ERC=' '
      IFLG=-9
      IF (ITYP.EQ.20) THEN
         CALL RDDT20(IUI,2,ALINE,NW,NA,SA
     &                          ,U,AP,HP,FH,FV,SH,SV,PT,RL,MFLG,ERC)
      ELSE IF (ITYP.EQ.21) THEN
         CALL RDDT21(IUI,2,ALINE,NW,NA,SA,UTM
     &                          ,U,AP,HP,FH,FV,SH,SV,PT,RL,MFLG,ERC)
      ELSE
         CALL RDDT22(IUI,2,ALINE,NW,NA,SA,UTM
     &                          ,U,AP,HP,FH,FV,SH,SV,PT,RL,MFLG,ERC)
      ENDIF
      IF (ERC.NE.' ') GOTO 900
C
      ERC=ALINE(1:8)//'  '
      IFLG=0
      IF (NA.LE.0) THEN
         ERC(11:25)='DTSHIP: NO DATA'
         IFLG=-8
         GOTO 900
      ENDIF
C
      IF (ITYP.EQ.20) THEN
         CALL RDLB20(ALINE,IY,IM,ID,JH,JM,JS,ALAT,ALNG,JFLG)
         TM=JH+JM/60.0+JS/3600.0
         DO J=1,NA
           T=TM+TISEC*(J-1)/3600.0
           CALL SUNH(IY,IM,ID,T,ALAT,ALNG,ALNGS,DT,DST,HGT1,FAI1)
           FIS(J)=DEG(FAI1,0)
           HTS(J)=HGT1
         ENDDO
      ELSE
         CALL RDLB21(ALINE,IY,IM,ID,JH,JM,JS,ALAT,ALNG,JFLG
     &                    ,IYU,IMU,IDU,JHU,JMU,JSU)
         TM=JH+JM/60.0+JS/3600.0
         TM1=JHU+JMU/60.0+JSU/3600.0
         DO J=1,NA
           IF (UTM(J).LT.0) THEN
              T=TM1+TISEC*(J-1)/3600.0
              CALL SUNH(IYU,IMU,IDU,T,ALAT,ALNG,0.,DT,DST,HGT1,FAI1)
           ELSE
              T=UTM(J)
              CALL SUNH(IYU,IMU,IDU,T,ALAT,ALNG,0.,DT,DST,HGT1,FAI1)
           ENDIF
           FIS(J)=DEG(FAI1,0)
           HTS(J)=HGT1
         ENDDO
      ENDIF
      FAI=FIS(1)
      TH0=90.-HTS(1)
C
      CALL INTENS(NW,SOLID,NA,HTS,U,U2,DST,IFLG)
      IF (IFLG.NE.0) THEN
         IF (IFLG.EQ.-1) ERC(11:32)='DTSHIP: NO DIRECT DATA'
         IF (IFLG.EQ.-2) ERC(11:33)='DTSHIP: NO DIFFUSE DATA'
         GOTO 900
      ENDIF
C
      DO J=1,NA
        IF (HP(1,J).GT.90.) THEN
           KFLG(J)=-1
        ELSE
           KFLG(J)=1
        ENDIF
        DO I=1,NW
          AP(I,J)=DEG(AP(I,J),0)-DFAI
          HP(I,J)=HP(I,J)-DHT
          IF (KFLG(J).LT.0) THEN
             HP(I,J)=180.-HP(I,J)
             AP(I,J)=DEG(AP(I,J)+180.,0)
          ENDIF
          PT(I,J)=PT(I,J)-PT0
          RL(I,J)=RL(I,J)-RL0
        ENDDO
      ENDDO
C
      IF (IANG.GT.0) THEN
         CALL AGCRCT(IUI,IDBG,NG1,NG2,FIS,HTS,SA
     &           ,NW,NA,KFLG,AP,HP,FH,FV,SH,SV,PT,RL,APC,HPC,JFG
     &           ,SH0,SV0,ASPN,CSA,DZN,FH0,FV0,ASPF,CO,DFI,FC,ERC,IFLG)
         IF (IFLG.NE.0) GOTO 900
      ELSE
         APM=0.
         HPM=0.
         DO I=1,NW
           APM=APM+AP(I,1)
           HPM=HPM+HP(I,1)
         ENDDO
         APM=APM/FLOAT(NW)
         HPM=HPM/FLOAT(NW)
         JFG=0
         DAP=DEG(APM,JFG)-DEG(FIS(1),JFG)
         DO I=1,NW
           DHP=HPM-HTS(1)
           APC(I,1)=AP(I,1)-DAP
           HPC(I,1)=HP(I,1)-DHP
           DO J=2,NA
             APC(I,J)=AP(I,J)-DAP
             IF (ABS(APC(I,J)-APC(I,J-1)).GT.175.) DHP=-DHP
             HPC(I,J)=HP(I,J)-DHP
           ENDDO
         ENDDO
      ENDIF
C
      DO J=1,NA
        N=0
        THM(J)=0.
        FIM(J)=0.
        DO  I=1,NW
          IF (HPC(I,J).GT.0.) THEN
             TH(I,J)=90.-HPC(I,J)
             FI(I,J)=DEG(APC(I,J),JFG)-DEG(FIS(J),JFG)
             FI(I,J)=DEG(FI(I,J),JFG)
             N=N+1
             THM(J)=THM(J)+TH(I,J)
             FIM(J)=FIM(J)+FI(I,J)
          ELSE
             TH(I,J)=-99.
             FI(I,J)=-999.
          ENDIF
        ENDDO
        IF (N.GT.0) THEN
           THM(J)=THM(J)/FLOAT(N)
           FIM(J)=FIM(J)/FLOAT(N)
        ELSE
           THM(J)=-99.
           FIM(J)=-999.
        ENDIF
      ENDDO
C
  900 RETURN
      END
C
      SUBROUTINE INTENS(NW,SOLID,NA,HTS,U,U2,DST,IFLG)
      PARAMETER (KNW   =20)
      PARAMETER (KNA   =50)
C
      PARAMETER (PI=3.141592653590D0,RD=PI/180.0D0)
C
      DIMENSION SOLID(KNW),HTS(KNA),U(KNW,KNA),U2(KNW,KNA)
C
      N=0
      DO I = 1,NW
        IF (U(I,1).GT.0.) THEN
           U(I,1)=U(I,1)*DST**2
           N=N+1
         ELSE
           U(I,1)=0.
        ENDIF
        U2(I,1)=U(I,1)
      ENDDO
      IF (N.LE.0) THEN
         IFLG=-1
         GOTO 900
      ENDIF
C
      N=0
      IF (NA.GE.2) THEN
         DO J=2,NA
           EM=1/COS((90.-HTS(J))*RD)
           DO I = 1,NW
             IF (U(I,J).GT.0.) THEN
                U2(I,J)=U(I,J)*DST**2/SOLID(I)
                IF (U(I,1).GT.0.) THEN
                   U(I,J)=U2(I,J)/U(I,1)/EM
                   N=N+1
                ELSE
                   U(I,J)=0.
                ENDIF
             ELSE
                U(I,J)=0.
                U2(I,J)=U(I,J)
             ENDIF
           ENDDO
         ENDDO
      ENDIF
      IF (N.LE.0) IFLG=-2
C
  900 RETURN
      END
C
      SUBROUTINE AGCRCT(IUI,IDBG,NG1,NG2,FIS,HTS,SA
     &           ,NW,NA,KFLG,AP,HP,FH,FV,SH,SV,PT,RL,APC,HPC,JFG
     &           ,SH0,SV0,ASPN,CSA,DZN,FH0,FV0,ASPF,CO,DFI,FC,ERC,IFLG)
      PARAMETER (KNW   =20)
      PARAMETER (KNA   =50)
C
      PARAMETER (PI=3.141592653590D0,RD=PI/180.0D0)
C
      CHARACTER ERC*(*)
      DIMENSION FC(5)
      DIMENSION AP(KNW,KNA),HP(KNW,KNA),KFLG(KNA)
     &         ,FH(KNW,KNA),FV(KNW,KNA),SH(KNW,KNA),SV(KNW,KNA)
     &         ,PT(KNW,KNA),RL(KNW,KNA)
      DIMENSION FIS(KNA),HTS(KNA),SA(KNA),APC(KNW,KNA),HPC(KNW,KNA)
C
      DATA DPMX,DFMX,DNMX /10., 2., 0.5/
C
      J=1
      FAI=FIS(J)
      HGT=HTS(J)
      KFG=KFLG(J)
      SA0=SA(J)
      JFGG=-1
      JFG=-1
      N=0
      DO WHILE((JFG.LT.0).AND.(N.LT.NW))
        N=N+1
        IF (HP(N,J).GT.0.) THEN
           PT1=PT(N,J)
           RL1=RL(N,J)
           AP1=AP(N,J)
           HP1=HP(N,J)
           CALL GTSUNP(AP1,HP1,PFI,PHT,FGP,PT1,RL1,FAI,HGT,IFGP)
           IF (IFGP.GT.0) THEN
              IF (ABS(FGP).GE.135.) THEN
                 JFGG=1
              ELSE
                 JFGG=0
              ENDIF
              SAP=SCAANG(PFI,90.-PHT,FAI,90.-HGT)
              IF (ABS(SAP).LE.DPMX) JFG=JFGG
           ENDIF
        ENDIF
      ENDDO
      IF (JFG.LT.0) JFG=JFGG
      IF (JFG.LT.0) THEN
         ERC(11:32)='AGCRCT: NO FIRST ANGLE'
         IFLG=-3
         GOTO 900
      ENDIF
C
      CALL MESPOS(NW,NA,J,KFG,SA0,AP,HP,FH,FV,SH,SV,PT,RL,FAI,HGT
     &           ,SH0,SV0,ASPN,CSA,DZN,FH0,FV0,ASPF,CO,DFI,FC
     &           ,DPMX,DFMX,DNMX
     &           ,NP,AP0,HP0,FGP0,NF,AF0,HF0,FGF0,NN,AN0,HN0,FGN0,JFG)
      FG0=FGP0
      IF (IDBG.GT.0) THEN
         NG1=NG1+1
         WRITE(IUI+4,10) NG1,J,AF0,HF0,FGF0,AN0,HN0,FGN0,AP0,HP0,FGP0
   10    FORMAT(I5,I3,9F8.2)
      ENDIF
C
      DO I=1,NW
        IF (HP(I,J).GT.0.) THEN
           PT1=PT(I,J)
           RL1=RL(I,J)
           AP1=AP(I,J)
           HP1=HP(I,J)
           FH1=FH(I,J)
           FV1=FV(I,J)
           SH1=SH(I,J)
           SV1=SV(I,J)
           AP2=AP1
           IF (KFG.LT.0) AP2=DEG(AP1-180.,JFG)
           CALL PSFISH(FH1,FV1,AP2,FH0,FV0,ASPF,CO,DFI,FC,AF,HF)
           CALL GTSUNP(AF,HF,FFI,FHT,FGF,PT1,RL1,FAI,HGT,IFGF)
           CALL PSNARW(SH1,SV1,AP1,HP1,SH0,SV0,ASPN,CSA,DZN,SA1,AN,HN)
           CALL GTSUNP(AN,HN,SFI,SHT,FGN,PT1,RL1,FAI,HGT,IFGN)
           CALL GTSUNP(AP1,HP1,PFI,PHT,FGP,PT1,RL1,FAI,HGT,IFGP)
           IF (ABS(FGN).LE.360.) THEN
              FG0=FGN
           ELSE IF (ABS(FGN0).LE.360.) THEN
              FG0=FGN0
           ELSE IF (ABS(FGP).LE.360.) THEN
              FG0=FGP
           ELSE IF (ABS(FGP0).LE.360.) THEN
              FG0=FGP0
           ELSE IF (ABS(FGF).LE.360.) THEN
              FG0=FGF
           ELSE IF (ABS(FGF0).LE.360.) THEN
              FG0=FGF0
           ELSE
              FG0=-999.
           ENDIF
           IF (ABS(FG0).LE.360.) THEN
              CALL GTANGP(AP1,HP1,PFI,PHT,FG0,PT1,RL1,IFGP)
              CALL GTANGP(AF,HF,FFI,FHT,FG0,PT1,RL1,IFGF)
              CALL GTANGP(AN,HN,SFI,SHT,FG0,PT1,RL1,IFGN)
              IF (IDBG.GT.0) THEN
                 NG2=NG2+1
                 SAP=SCAANG(PFI,90.-PHT,FAI,90.-HGT)
                 SAF=SCAANG(FFI,90.-FHT,FAI,90.-HGT)
                 SAN=SCAANG(SFI,90.-SHT,FAI,90.-HGT)
c                write(*,10) NG2,J,FGF,FGF0,FGN,FGN0,FGP,FGP0,FG0
                 WRITE(IUI+1,10) NG2,J,AF,HF,FGF,AN,HN,FGN
     &                                          ,AP1,HP1,FG0
                 WRITE(IUI+2,10) NG2,J,SA0,SAP,SA1,SAF,SAN
     &                                          ,PFI,FFI,SFI,FAI
                 WRITE(IUI+3,10) NG2,J,SA0,SAP,SA1,SAF,SAN
     &                                          ,PHT,FHT,SHT,HGT
              ENDIF
              APC(I,J)=PFI
              HPC(I,J)=PHT
           ELSE
              APC(I,J)=-999.
              HPC(I,J)=-99.
           ENDIF
        ELSE
           APC(I,J)=-999.
           HPC(I,J)=-99.
        ENDIF
      ENDDO
C
      DO J=2,NA
        FAI1=FIS(J)
        HGT1=HTS(J)
        KFG=KFLG(J)
        SA0=SA(J)
        CALL MESPOS(NW,NA,J,KFG,SA0,AP,HP,FH,FV,SH,SV,PT,RL,FAI1,HGT1
     &        ,SH0,SV0,ASPN,CSA,DZN,FH0,FV0,ASPF,CO,DFI,FC
     &        ,DPMX,DFMX,DNMX
     &        ,NP,AP0,HP0,FGP0,NF,AF0,HF0,FGF0,NN,AN0,HN0,FGN0,JFG)
        IF (NN.GT.0) THEN
           FG0=FGN0
        ELSE
           FG0=FGF0
        ENDIF
        IF (IDBG.GT.0) THEN
           NG1=NG1+1
           WRITE(IUI+4,10) NG1,J,AF0,HF0,FGF0,AN0,HN0,FGN0,AP0,HP0,FGP0
        ENDIF
C
        DO I = 1,NW
          IF (HP(I,J).GT.0.) THEN
             PT1=PT(I,J)
             RL1=RL(I,J)
             AP1=AP(I,J)
             HP1=HP(I,J)
             FH1=FH(I,J)
             FV1=FV(I,J)
             SH1=SH(I,J)
             SV1=SV(I,J)
             AP2=AP1
             IF (KFG.LT.0) AP2=DEG(AP1-180.,JFG)
             CALL PSFISH(FH1,FV1,AP2,FH0,FV0,ASPF,CO,DFI,FC,AF,HF)
             CALL GTSUNP(AF,HF,FFI,FHT,FGF,PT1,RL1,FAI1,HGT1,IFGF)
             SAF=SCAANG(FFI,90.-FHT,FAI1,90.-HGT1)
             CALL PSNARW(SH1,SV1,AP1,HP1,SH0,SV0,ASPN,CSA,DZN,SA1,AN,HN)
             CALL GTSUNP(AN,HN,SFI,SHT,FGN,PT1,RL1,FAI,HGT,IFGN)
             SAN=SCAANG(SFI,90.-SHT,FAI1,90.-HGT1)
             IF (ABS(FGN).LE.360.) THEN
                FG0=FGN
             ELSE IF (ABS(FGN0).LE.360.) THEN
                FG0=FGN0
             ELSE IF (ABS(FGF).LE.360.) THEN
                FG0=FGF
             ELSE IF (ABS(FGF0).LE.360.) THEN
                FG0=FGF0
             ELSE
                FG0=-999.
             ENDIF
             IF (ABS(FG0).LE.360.) THEN
                CALL GTANGP(AP1,HP1,PFI,PHT,FG0,PT1,RL1,IFGP)
                CALL GTANGP(AF,HF,FFI,FHT,FG0,PT1,RL1,IFGF)
                CALL GTANGP(AN,HN,SFI,SHT,FG0,PT1,RL1,IFGN)
                SAP=SCAANG(PFI,90.-PHT,FAI1,90.-HGT1)
                SAF=SCAANG(FFI,90.-FHT,FAI1,90.-HGT1)
                SAN=SCAANG(SFI,90.-SHT,FAI1,90.-HGT1)
                IF (IDBG.GT.0) THEN
                   NG2=NG2+1
                   FGP=FGP0
c                  write(*,10) NG2,J,FGF,FGF0,FGN,FGN0,FGP,FGP0,FG0
                   WRITE(IUI+1,10) NG2,J,AF,HF,FGF,AN,HN,FGN
     &                                      ,AP1,HP1,FG0
                   WRITE(IUI+2,10) NG2,J,SA0,SAP,SA1,SAF,SAN
     &                                      ,PFI,FFI,SFI,FAI1
                   WRITE(IUI+3,10) NG2,J,SA0,SAP,SA1,SAF,SAN
     &                                      ,PHT,FHT,SHT,HGT1
                ENDIF
                APC(I,J)=PFI
                HPC(I,J)=PHT
             ELSE
                APC(I,J)=-999.
                HPC(I,J)=-99.
             ENDIF
          ELSE
             APC(I,J)=-999.
             HPC(I,J)=-99.
          ENDIF
        ENDDO
      ENDDO
C
  900 RETURN
      END
C
C ----------------------------------------------------------------------
      SUBROUTINE MESPOS(NW,NA,J,KFG,SA0,AP,HP,FH,FV,SH,SV,PT,RL,FAI,HGT
     &           ,SH0,SV0,ASPN,CSA,DZN,FH0,FV0,ASPF,CO,DFI,FC
     &           ,DPMX,DFMX,DNMX
     &           ,NP,AP0,HP0,FGP0,NF,AF0,HF0,FGF0,NN,AN0,HN0,FGN0,JFG)
      PARAMETER (KNW   =20)
      PARAMETER (KNA   =50)
C
      DIMENSION FC(5)
      DIMENSION AP(KNW,KNA),HP(KNW,KNA)
     &         ,FH(KNW,KNA),FV(KNW,KNA),SH(KNW,KNA),SV(KNW,KNA)
     &         ,PT(KNW,KNA),RL(KNW,KNA)
      DIMENSION APM(KNW),HPM(KNW),AFM(KNW),HFM(KNW),ANM(KNW),HNM(KNW)
     &         ,FGPM(KNW),FGFM(KNW),FGNM(KNW)
C
      NP=0
      NF=0
      NN=0
      DO I=1,NW
        IF (HP(I,J).GT.0.) THEN
           PT1=PT(I,J)
           RL1=RL(I,J)
           AP1=AP(I,J)
           HP1=HP(I,J)
           FH1=FH(I,J)
           FV1=FV(I,J)
           SH1=SH(I,J)
           SV1=SV(I,J)
           AP2=AP1
           IF (KFG.LT.0) AP2=DEG(AP1-180.,JFG)
           CALL PSFISH(FH1,FV1,AP2,FH0,FV0,ASPF,CO,DFI,FC,AF,HF)
           CALL GTSUNP(AF,HF,FFI,FHT,FGF,PT1,RL1,FAI,HGT,IFGF)
           CALL PSNARW(SH1,SV1,AP1,HP1,SH0,SV0,ASPN,CSA,DZN,SA1,AN,HN)
           CALL GTSUNP(AN,HN,SFI,SHT,FGN,PT1,RL1,FAI,HGT,IFGN)
           IF (J.EQ.1) THEN
              CALL GTSUNP(AP1,HP1,PFI,PHT,FGP,PT1,RL1,FAI,HGT,IFGP)
              IF (IFGP.GT.0) THEN
                 SAP=SCAANG(PFI,90.-PHT,FAI,90.-HGT)
                 IF (ABS(SAP).LE.DPMX) THEN
                    NP=NP+1
                    APM(NP)=PFI
                    HPM(NP)=PHT
                    FGPM(NP)=DEG(FGP,JFG)
                 ENDIF
              ENDIF
           ENDIF
           IF (IFGF.GT.0) THEN
              SAF=SCAANG(FFI,90.-FHT,FAI,90.-HGT)
              IF (ABS(SAF).LE.DFMX) THEN
                 NF=NF+1
                 AFM(NF)=FFI
                 HFM(NF)=FHT
                 FGFM(NF)=DEG(FGF,JFG)
              ENDIF
           ENDIF
           IF (IFGN.GT.0) THEN
              SAN=SCAANG(SFI,90.-SHT,FAI,90.-HGT)
              IF (ABS(SA1-SA0).LE.DNMX) THEN
                 NN=NN+1
                 ANM(NN)=SFI
                 HNM(NN)=SHT
                 FGNM(NN)=DEG(FGN,JFG)
              ENDIF
           ENDIF
        ENDIF
      ENDDO
C
      IF (NP.GT.0) THEN
         CALL CLMEAN(NP,APM,AP0,STD,RMS)
         CALL CLMEAN(NP,HPM,HP0,STD,RMS)
         CALL CLMEAN(NP,FGPM,FGP0,STD,RMS)
      ELSE
         HP0=-99.
         AP0=-999.
         FGP0=-999.
      ENDIF
      IF (NF.GT.0) THEN
         CALL CLMEAN(NF,AFM,AF0,STD,RMS)
         CALL CLMEAN(NF,HFM,HF0,STD,RMS)
         CALL CLMEAN(NF,FGFM,FGF0,STD,RMS)
      ELSE
         HF0=-99.
         AF0=-999.
         FGF0=-999.
      ENDIF
      IF (NN.GT.0) THEN
         CALL CLMEAN(NN,ANM,AN0,STD,RMS)
         CALL CLMEAN(NN,HNM,HN0,STD,RMS)
         CALL CLMEAN(NN,FGNM,FGN0,STD,RMS)
      ELSE
         HN0=-99.
         AN0=-999.
         FGN0=-999.
      ENDIF
      RETURN
      END
C
      SUBROUTINE GTSUNP(AP,HP,PFI,PHT,FG,PT,RL,FAI,HGT,IFGP)
C
      IF ((ABS(AP).LE.360.).AND.(HP.GT.0.)) THEN
         CALL GETXYZ(FAI,HGT,RX,RY,RZ)
         CALL GETXYZ(AP,HP,PXB,PYB,PZB)
         CALL ROTATU(-RL,PZB,PXB,PZ1,PX1)
         CALL ROTATU(-PT,PYB,PZ1,PY1,PZR)
         CALL CLDRFT(PX1,PY1,RX,RY,FG)
         CALL ROTATU( FG,PX1,PY1,PXR,PYR)
         CALL GETPOL(PXR,PYR,PZR,PFI,PHT)
         PFI=DEG(PFI,0)
         IFGP=1
      ELSE
         PFI=-999.
         PHT=-99.
         FG=-999.
         IFGP=-1
      ENDIF
      RETURN
      END
C
      SUBROUTINE GTANGP(AP,HP,PFI,PHT,FG,PT,RL,IFGP)
C
      IF ((ABS(AP).LE.360.).AND.(HP.GT.0.)) THEN
         CALL GETXYZ(AP,HP,PXB,PYB,PZB)
         CALL ROTATU(-RL,PZB,PXB,PZ1,PX1)
         CALL ROTATU(-PT,PYB,PZ1,PY1,PZR)
         CALL ROTATU( FG,PX1,PY1,PXR,PYR)
         CALL GETPOL(PXR,PYR,PZR,PFI,PHT)
         PFI=DEG(PFI,0)
         IFGP=1
      ELSE
         PFI=-999.
         PHT=-99.
         IFGP=-1
      ENDIF
      RETURN
      END
C
      SUBROUTINE RDDTFP(IU,IDT3,IDT4,IDT5,ILST,IDBG,IANG
     &                 ,CMNT,SITE,CNTY,ALNGS,ALNG,ALAT,ALT,NA,SA
     &                 ,SRNO,ITYP,NW,WL,SOLID
     &                 ,PT0,RL0,DHT,DFAI
     &                 ,SH0,SV0,ASPN,CSA,DZN,FH0,FV0,ASPF,CO,DFI,FC,ERC)
C
C   Read parameters from file('dtform.par')
C
C --- history
C   2001.01.26 Created by M.Yamano
C   2001.10.07 Renewed.
C   2002.05.18 Renewed for Ver5.
C   2002.11.20 Modified.
C   2002.12.02 Changed the order of parameters.
C   2003.03.06 Modified.
C
C --- Input
C IU            I     device No. for reading
C
C --- Output
C IDT3          I     creation option of data file for ver.3(1: create/0: not)
C IDT4          I     creation option of data file for ver.4(1: create/0: not)
C IDT5          I     creation option of data file for ver.5(1: create/0: not)
C ILST          I     creation option of list file(1: create/0: not)
C IDBG          I     creation option of debug files(1: create/0: not)
C IANG          I     angle correction option for POM-01MKII(1: correct/0: not)
C CMNT          C*40  project name
C SITE          C*20  observation site (or ship) name
C CNTY          C*20  country name
C ALNGS         R     standard longitude[deg] for the time zone
C ALNG          R     longitude of observation site[deg]
C ALAT          R     latitude of observation site[deg]
C ALT           R     height of observation site[m]
C NA            I     number of observation scattering angles
C SA            R(NA) scattering angles[deg]
C SRNO          C*20  instrument S/N
C ITYP          I     data file format type No.
C                       10,11   : PREDE - skyradiometer POM-01L
C                       20,21,22: PREDE - skyradiometer POM-01MKII
C NW            I     number of wavelengths
C WL            R(NW) wavelengths[cm]
C SOLID         R(NW) solid view angles[sr]
C PT0,RL0       R     pitching/rolling offset angle[deg]
C DHT,DFAI      R     height offset(=Hs-Hp)/azimuth offset(=Aps-As)[deg]
C SH0,SV0,ASPN  R     narrow CCD parameters
C CSA,DZ
C FH0,FV0,ASPF  R     fish CCD parameters
C CO,DFI,FC(5)
C ERC           C*64  ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KNA   =50)
      PARAMETER (KDAY  =30)
C
      CHARACTER ERC*(*)
      CHARACTER PFIL*80
      CHARACTER CMNT*(*),SITE*(*),CNTY*(*),SRNO*(*)
      DIMENSION SA(KNA),FC(5)
      DIMENSION WL(KNW),SOLID(KNW),TJLC(KDAY),F0D(KNW,KDAY)
C
      ERC=' '
      READ(IU,*,ERR=800,END=800) IDT3,IDT4,IDT5,ILST,IDBG
      READ(IU,*,ERR=800,END=800) IANG
C
C --- read 'obs.para' parameters
C
      READ(IU,*,ERR=800,END=800) PFIL
      CALL UTLSPC(PFIL,NL)
      IF (NL.LE.0) THEN
         ERC='RDDTFP: illegal obs.para filename'
         GOTO 900
      ENDIF
      OPEN (IU+1,FILE=PFIL(1:NL),STATUS='OLD')
      CALL RDOBSP(IU+1,CMNT,SITE,CNTY,ALNGS,ALNG,ALAT,ALT,NA,SA,ERC)
      CLOSE (IU+1)
      IF (ERC.NE.' ') GOTO 900
C
C --- read 'ins.para' parameters
C
      READ(IU,*,ERR=800,END=800) PFIL
      CALL UTLSPC(PFIL,NL)
      IF (NL.LE.0) THEN
         ERC='RDDTFP: illegal ins.para filename'
         GOTO 900
      ENDIF
      OPEN (IU+1,FILE=PFIL(1:NL),STATUS='OLD')
      CALL RDINSP(IU+1,SRNO,ITYP,NW,WL,SOLID,NDY,F0D,TJLC,ERC)
      CLOSE (IU+1)
      IF (ERC.NE.' ') GOTO 900
C
      IF (ITYP.LT.20) THEN
         IANG=0
         IDBG=0
      ELSE
         IF (IANG.EQ.0) IDBG=0
      ENDIF
C
C --- read 'CCD.para' parameters
C
      IF (IANG.GT.0) THEN
         READ(IU,*,ERR=800,END=800) PFIL
         CALL UTLSPC(PFIL,NL)
         IF (NL.LE.0) THEN
            ERC='RDDTFP: illegal CCD.para filename'
            GOTO 900
         ENDIF
         OPEN (IU+1,FILE=PFIL(1:NL),STATUS='OLD')
         CALL RDCCDP(IU+1,PT0,RL0,DHT,DFAI,SH0,SV0,ASPN,CSA,DZN
     &                   ,FH0,FV0,ASPF,CO,DFI,FC,ERC)
         CLOSE (IU+1)
      ELSE
         PT0=0.
         RL0=0.
         DHT=0.
         DFAI=0.
      ENDIF
      GOTO 900
C
  800 ERC='RDDTFP: Read Error !'
C
  900 RETURN
      END
C
      SUBROUTINE RDHD30(IU,SRNO,ALNGS,ALNG,ALAT,NW,WL,ERC)
C
C   Read file header information for PREDE - skyradiometer POM-01,-02.
C   Ver.4.10 data format(2003 - ):ITYP=30
C
C --- history
C   2003.11.05 Created by M.Yamano
C   2003.11.19 Added output argument ALNGS
C   2003.12.05 Debugged.
C
C --- Input
C IU         I         device No. for reading
C
C --- Output
C SRNO       C*20      instrument S/N
C ALNGS      R         standard longitude in degree for the time zone
C ALAT       R         latitude in degree
C ALNG       R         longitude in degree
C NW         I         number of wavelengths
C WL         R(NW)     wavelengths[cm]
C ERC        C*64      ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
C
      CHARACTER ERC*(*)
      CHARACTER ALIN1*128,ALIN2*128,OTHER*128
      CHARACTER SRNO*20
      DIMENSION WL(KNW)
C
      ERC=' '
      READ(IU,10,ERR=800,END=800) ALIN1
      READ(IU,10,ERR=800,END=800) ALIN2
   10 FORMAT(A128)
C
      DO I=1,2
        CALL GETPOS(ALIN1,',',128,IL)
        IF (IL.EQ.0) GOTO 800
        ALIN1=ALIN1(IL+1:128)
      ENDDO
      CALL GETPOS(ALIN1,',',128,IL)
      IF (IL.LE.1) GOTO 800
      SRNO=ALIN1(1:IL-1)
      OTHER=ALIN1(IL+1:128)
      READ(OTHER,*) ALNG,ALAT
c     print*,'-',SRNO,'-',ALNG,ALAT
C
      DO I=1,2
        CALL GETPOS(OTHER,',',128,IL)
        IF (IL.LE.1) GOTO 800
        OTHER=OTHER(IL+1:128)
      ENDDO
      READ(OTHER,20) IY,IT,ID,IH,IM,IS,LY,LT,LD,LH,LM,LS
   20 FORMAT(12(I2,1X))
c     write(*,20) IY,IT,ID,IH,IM,IS,LY,LT,LD,LH,LM,LS
      JFG=LY-IY
      IF (JFG.EQ.0) JFG=LT-IT
      IF (JFG.EQ.0) JFG=LD-ID
      IF ((ABS(JFG).GT.1).OR.(LS.NE.IS)) GOTO 810
      MM=ABS(LM-IM)
      IF ((MM.NE.0).AND.(MM.NE.30)) GOTO 810
      DTM=FLOAT(LH-IH)+FLOAT(LM-IM)/60.+FLOAT(JFG)*24.
      ALNGS=DTM*15.
C
      CALL GETPOS(ALIN2,',',128,IL)
      IF (IL.LE.1) GOTO 800
      READ(ALIN2,*) NW
      OTHER=ALIN2(IL+1:128)
      READ(OTHER,*,ERR=800,END=800) (WL(I),I=1,NW)
c     print*, (WL(I),I=1,NW)
      GOTO 900
C
  800 ERC='RDHD30: Read error.'
      GOTO 900
  810 ERC='RDHD30: Time mismatch !'
C
  900 RETURN
      END
C
      SUBROUTINE WTTAG0(IU,NTG,IY,IM,ID,TM,ALNG,ALAT,TH0,SAMX)
C
C   Write tag data into tag file(yymmdd.tag).
C
C --- history
C   2002.05.11 Created by M.Yamano
C
C --- Input
C IU          I         device No. for writing
C NTG         I         data No.
C IY,IM,ID    I         measured date IY/IM/ID
C TM          R         measured time in hour
C ALNG,ALAT   R         longitude,latitude in degree of measurement site
C TH0         R         solar zenith angle in degree
C SAMX        R         max, scattering angle in degree
C ----------------------------------------------------------------------
C
      IF (NTG.EQ.1) WRITE(IU,10)
   10 FORMAT('  No yyyy mm dd Hour   Long    Lat    Hs    SA(max)')
      WRITE(IU,20) NTG,IY,IM,ID,TM,ALNG,ALAT,90.0-TH0,SAMX
   20 FORMAT(I4,I5,2I3,F6.2,F8.2,2F7.2,F6.1)
      RETURN
      END
C
      SUBROUTINE WTDAT3(IU,IY,IM,ID,IHH,IMM,ISS,TM,ALNGS,ALNG,ALAT,ALT
     &                 ,NA,TH0,FI0,DST,NW,WL,TH,FI,AUR)
C
C   Write datafile(yymmdd.DT3) for analysis Ver.3.0.
C
C --- History
C   2001.02.23 Created by M.Yamano
C
C --- Input
C IU          I         device No. for writing
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
C TH,FI       R(NA)     zenith and azimuth angles of measured direction
C AUR         R(NW,NA)  intensity F(sca=0) or I/EM/F/(solid view angle)
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KNA   =50)
C
      DIMENSION WL(KNW),TH(KNA),FI(KNA),AUR(KNW,KNA)
C
      WRITE(IU,10) IY,IM,ID,IHH,IMM,ISS,TM
   10 FORMAT(I8,5I3,F9.4,' : IY IM ID IHH IMM ISS TM(hr)')
      WRITE(IU,20) ALNGS,ALNG,ALAT,ALT
   20 FORMAT(4F8.2,' : ALNGS ALNG ALAT ALT(m)')
      WRITE(IU,30) NA,TH0,FI0,DST
   30 FORMAT(I8,2F8.1,F8.4,' : NA TH0 FI0 DST')
C
      WRITE(IU,50) NW
   50 FORMAT(I4, ' : NW/ WL')
      WRITE(IU,60) (WL(IW),IW=1,NW)
   60 FORMAT(13X,(0P5E13.3))
      WRITE(IU,80)
   80 FORMAT('  TH     FI     F(SCA=0)  or R=U/F/M/SOLID')
      DO I=1,NA
        WRITE(IU,90) TH(I),FI(I),(AUR(IW,I),IW=1,NW)
   90   FORMAT(F6.1,F7.1,(1P5E13.4))
      ENDDO
      RETURN
      END
C
      SUBROUTINE WTDAT4(IU,IY,IM,ID,IHH,IMM,ISS,TM,ALNGS,ALNG,ALAT,ALT
     &                 ,NA,TH0,FI0,DST,NW,WL,TH,FI,AUR)
C
C   Write datafile(yymmdd.DT4) for analysis Ver.4.0.
C
C --- History
C   2001.01.16 Created by M.Yamano
C
C --- Input
C IU          I         device No. for writing
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
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KNA   =50)
C
      DIMENSION WL(KNW),TH(KNW,KNA),FI(KNW,KNA),AUR(KNW,KNA)
C
      WRITE(IU,10) IY,IM,ID,IHH,IMM,ISS,TM
   10 FORMAT(I8,5I3,F9.4,' : IY IM ID IHH IMM ISS TM(hr)')
      WRITE(IU,20) ALNGS,ALNG,ALAT,ALT
   20 FORMAT(4F8.2,' : ALNGS ALNG ALAT ALT(m)')
      WRITE(IU,30) NA,TH0,FI0,DST
   30 FORMAT(I8,2F8.1,F8.4,' : NA TH0 FI0 DST')
C
      WRITE(IU,50) NW
   50 FORMAT(I4, ' : NW/ WL')
      WRITE(IU,60) (WL(IW),IW=1,NW)
   60 FORMAT(3(14X,1PE12.3))
      WRITE(IU,80)
   80 FORMAT('   TH     FI    F (IF SCA=0)  or R=U/F/M/SOLID')
      DO I=1,NA
        WRITE(IU,90) (TH(IW,I),FI(IW,I),AUR(IW,I),IW=1,NW)
   90   FORMAT(3(0P2F7.1,1PE12.4))
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
      SUBROUTINE RDDT30(IU,NW,IY,IM,ID,JH,JM,JS,JFLG,NA,AP,HP,EM,ERC)
C
C   Read measured data for PREDE - skyradiometer POM-01,-02.
C   for datafile format Ver.4.10 (2003 - ):ITYP=30
C
C --- history
C   2003.11.18 Created by M.Yamano
C   2006.03.13 Debugged. (format(yy) -> 20yy)
C
C --- Input
C IU     I         device No. for reading
C NW     I         number of wavelengths
C
C --- Output
C IY,IM,ID   I         date(yyyy/mm/dd) - LT
C JH,JM,JS   I         time(hh:mm:ss)   - LT
C JFLG       I         flag of scanning pattern
C                        1: Horizontal scanning
C                        2: Vertical   scanning
C NA         I         number of scattering angles
C AP         R(NA)     azimuth angle of measurement direction
C HP         R(NA)     height angle of measurement direction
C EM         R(NW,NA)  intensity
C ERC        C*64      ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KNA   =50)
C
      CHARACTER ERC*(*)
      CHARACTER ALIN*256,ADATE*8,ATIME*8,HV*1
      DIMENSION AP(KNA),HP(KNA),EM(KNW,KNA)
C
      ERC=' '
      READ(IU,10,ERR=800,END=800) ALIN
   10 FORMAT(A256)
      IF (ALIN(1:1).EQ.' ') ALIN=ALIN(2:256)
      READ(ALIN,15) ADATE,ATIME,HV
   15 FORMAT(18X,A8,1X,A8,1X,A1)
c     print*,ADATE,ATIME,HV
      READ(ADATE,20) IY,IM,ID
      IY=IY+2000
      READ(ATIME,20) JH,JM,JS
   20 FORMAT(I2,1X,I2,1X,I2)
      IF (HV.EQ.'H') THEN
         JFLG=1
      ELSE
         JFLG=2
      ENDIF
C
      NA=0
   50 NA=NA+1
      READ(IU,10,ERR=100,END=100) ALIN
      ALIN=ALIN(19:256)
      READ(ALIN,*) AP(NA),HP(NA),(EM(IW,NA),IW=1,NW)
      IF (NA.LT.KNA) GOTO 50
      GOTO 900
C
  100 NA=NA-1
      IF (NA.LE.0) ERC='RDDT30: Data Read error !'
      GOTO 900
C
  800 ERC='RDDT30: Header Read error !'
C
  900 RETURN
      END
C
      SUBROUTINE SUNH(IY,IM,ID,TM,ALAT,ALNG,ALNGS,DT,DST,HGT,FAI)
C--- HISTORY
C 88.10.20 CREATED
C 92. 2. 8 ADD A1
C 94.11.26 Change the comment for FAI
C--- INPUT
C IY      I      YEAR YYYY
C IM      I      MONTH
C ID      I      DAY
C TM      R      TIME (HOUR, I.E., 10:30 = 10.5)
C ALAT    R      LATITUDE  (DEGREE) (+: NORTH, -: SOUTH)
C ALNG    R      LONGITUDE (DEGREE) (+: EAST,  -: WEST)
C ALNGS   R      LONGITUDE FOR LOCAL STANDARD TIME (DEGREE)
C DT      R      DELTA(T)  51SEC/3600
C--- OUTPUT
C DST     R      SUN-EARTH DISTANCE (ASTRONOMICAL UNIT)
C HGT     R      SOLAR HEIGHT (degree, 90 when the sun is at zenith)
C FAI     R      AZIMUTHAL ANGLE (DEGREES, South=0, East=-90)
C
C                           180 North
C            -90(270) East                 90 West
C                             0 South
C
C In order to get compus-angle, you have to add the following two statements.
C
C     FAI=FAI-180
C     IF(FAI.LT.0) FAI=FAI+360
C
C then you have
C
C                             0 North
C                  90 East                270 West
C                           180 South
C$ENDI
      PARAMETER (PI=3.141592653590D0,RD=PI/180.0D0)
      U=TM-ALNGS/15
      E=U+DT
      XI=INT((14-IM)/12)
      XJ=-XI+IY+4800
      XJ=INT(1461*XJ/4.0)+INT((12*XI+IM-2)*367/12.0)
     & -INT(INT(XJ/100.0+1)*3/4.0)+ID-2447095.5+E/24
      T=XJ/36525
      C=0.01675104-4.18D-5*T
      X=358.47583+35999.04975*T-1.5D-4*T**2
      E=0
      DO 1 I=1,5
    1 E=X+C/RD*SIN(E*RD)
      V=2*ATAN(SQRT((1+C)/(1-C))*TAN(E/2*RD))/RD
      DST=1-C*COS(E*RD)
      G=281.22083+1.71917331*T+4.53D-4*T**2
      W=T-1
      X=20*COS((32964*W+158)*RD)+18*COS((19*W+159)*RD)
     & +18*COS((445267*W+208)*RD)+15*COS((45038*W+254)*RD)
      X=V+G+(X+13*COS((22519*W+352)*RD))/1.0D4
      W=259.183275-1934.142008*T+2.078D-3*T**2
      O=23.452294-0.0130125*T
      X=X-(478686*SIN(W*RD)+569328/DST)/1.0D8
      O=O+2.5583D-3*COS(W*RD)
      B=ASIN(SIN(O*RD)*SIN(X*RD))/RD
      A=ASIN(COS(O*RD)*SIN(X*RD)/COS(B*RD))/RD
      W=COS(X*RD)/COS(B*RD)
      IF(W.LT.0.0) A=180-A
      A=AMOD(A,360.0)
      Q=18.64606556+2400.051262*T
      G=U+12+Q
      T=15*G+ALNG-A
C
      HGT=ASIN(SIN(ALAT*RD)*SIN(B*RD)
     &    +COS(ALAT*RD)*COS(B*RD)*COS(T*RD))/RD
      A1=SIN(T*RD)*COS(B*RD)/COS(HGT*RD)
      IF(ABS(A1).GT.1.0) A1=SIGN(1.0,A1)
      FAI=ASIN(A1)/RD
      W=(-COS(ALAT*RD)*SIN(B*RD)
     &    +SIN(ALAT*RD)*COS(B*RD)*COS(T*RD))/COS(HGT*RD)
      IF(W.LT.0.0) FAI=180-FAI
      FAI=AMOD(FAI,360.0)
      RETURN
      END
      FUNCTION SCAANG(AP1,ZP1,AP0,ZP0)
C
C    Calculation of scattering angle.
C
C --- history
C   2000.09.01 Created by M.Yamano
C   2001.09.29 Modified.
C   2007.04.25 Debugged.
C
C --- Input
C AP1,ZP1       R         azimuth/zenith angles of point-1
C AP0,ZP0       R         azimuth/zenith angles of point-2
C
C --- Output
C SCAANG        R         scattering angle in degree between point-1 and 2
C ----------------------------------------------------------------------
C
      PARAMETER (PI=3.141592653590D0,RD=PI/180.0D0)
C
      IF ((AP1.GE.-360.).AND.(AP0.GE.-360.)
c    &                  .AND.(ZP1.GT.0.).AND.(ZP0.GT.0.)) THEN
     &                  .AND.(ZP1.GE.0.).AND.(ZP0.GE.0.)) THEN
         DP1=DEG(AP1,0)
         DP0=DEG(AP0,0)
         IF ((ABS(DP1-DP0).LE.1.E-8).AND.(ABS(ZP1-ZP0).LE.1.E-8)) THEN
            SCAANG=0.0
         ELSE
            SCAANG=ACOS(COS(ZP1*RD)*COS(ZP0*RD)
     &            +SIN(ZP1*RD)*SIN(ZP0*RD)*COS((DP1-DP0)*RD))/RD
         ENDIF
      ELSE
         SCAANG=-99.
      ENDIF
      RETURN
      END
C
      FUNCTION DEG(ANG,IANG)
C
C   Conversion of angle range.
C
C --- history
C   2000.08.18 Created by M.Yamano
C   2007.04.26 Debugged.
C
C --- Input
C ANG           R         angle in degree( > -990 deg.)
C IANG          I         conversion option
C                         -1: -360 <  DEG <=   0
C                          0: -180 <  DEG <= 180
C                          1:    0 <= DEG <  360
C --- Output
C DEG           R         converted angle in degree
C ----------------------------------------------------------------------
C
      IF (ANG.GT.-990.) THEN
         AG=ANG
         DO WHILE(AG.LT.0.)
           AG=AG+360.
         ENDDO
         AG=AMOD(AG,360.)
         IF (((IANG.LT.0).AND.(AG.GT.0.)).OR.
     &       ((IANG.EQ.0).AND.(AG.GT.180.))) THEN
            DEG=AG-360.
         ELSE
            DEG=AG
         ENDIF
      ELSE
         DEG=ANG
      ENDIF
      RETURN
      END
C
      SUBROUTINE RDDT11(IU,ALIN,NW,NA,EM,ERC)
C
C   Read measured data for PREDE - skyradiometer POM-01L.
C   for new datafile format( 1999 - ):ITYP=11
C
C --- history
C   2001.02.14 Created by M.Yamano
C
C --- Input
C IU     I         device No. for reading
C NW     I         number of wavelengths
C
C --- Output
C ALIN   C*128     header
C NA     I         number of scattering angles
C EM     R(NW,NA)  intensity
C ERC    C*64      ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KNA   =50)
C
      CHARACTER ERC*(*)
      CHARACTER ALIN*128
      DIMENSION EM(KNW,KNA)
C
      ERC=' '
      READ(IU,10,ERR=800,END=800) ALIN
   10 FORMAT(A128)
      NA=0
   50 NA=NA+1
      READ(IU,*,ERR=100,END=100) (EM(IW,NA),IW=1,NW)
      IF (NA.LT.KNA) GOTO 50
      GOTO 900
C
  100 NA=NA-1
      IF (NA.LE.0) ERC='RDDT11: Data Read error !'
      GOTO 900
C
  800 ERC='RDDT11: Header Read error !'
C
  900 RETURN
      END
C
      SUBROUTINE RDLB11(ALIN,IY,IM,ID,JH,JM,JS,JFLG)
C
C   Read header information for PREDE - skyradiometer POM-01L.
C   for new datafile format( 1999 - ):ITYP=11
C
C --- history
C   2001.02.14 Created by M.Yamano
C   2005.10.05 Modified to be applied to both formats for year(yyyy or yy).
C                                        format(yy) -> 1990 - 2089
C
C --- Input
C ALIN       C*128  header
C
C --- Output
C IY,IM,ID   I      date(yyyy/mm/dd)
C JH,JM,JS   I      time(hh:mm:ss)
C JFLG       I      flag of scanning pattern
C                     1: Horizontal scanning
C                     2: Vertical   scanning
C ----------------------------------------------------------------------
C
      CHARACTER ALIN*128,ATIME*8,ADATE*10,HV*1
C
      READ(ALIN,10) ATIME,HV,ADATE
   10 FORMAT(A8,1X,A1,1X,A10)
c     print*,ADATE,ATIME
C
      CALL GETPOS(ADATE,'/',10,IL)
      IF (IL.EQ.5) THEN
         READ(ADATE,15) IY,IM,ID
   15    FORMAT(I4,1X,I2,1X,I2)
      ELSE
         READ(ADATE,16) IY,IM,ID
   16    FORMAT(I2,1X,I2,1X,I2)
         IF (IY.GE.90) THEN
            IY=IY+1900
         ELSE
            IY=IY+2000
         ENDIF
      ENDIF
      READ(ATIME,20) JH,JM,JS
   20 FORMAT(I2,1X,I2,1X,I2)
      IF (HV.EQ.'H') THEN
         JFLG=1
      ELSE
         JFLG=2
      ENDIF
      RETURN
      END
C
      SUBROUTINE RDDT10(IU,ALIN,NW,NA,SA,EM,ERC)
C
C   Read measured data for PREDE - skyradiometer POM-01L.
C   for old datafile format( - 1998 ):ITYP=10
C
C --- history
C   2001.02.14 Created by M.Yamano
C
C --- Input
C IU     I         device No. for reading
C NW     I         number of wavelengths
C
C --- Output
C ALIN   C*128     header
C NA     I         number of scattering angles
C SA     R(NA)     scattering angles
C EM     R(NW,NA)  intensity
C ERC    C*64      ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KNA   =50)
C
      CHARACTER ERC*(*)
      CHARACTER ALIN*128
      DIMENSION EM(KNW,KNA),SA(KNA)
C
      ERC=' '
      READ(IU,10,ERR=800,END=800) ALIN
   10 FORMAT(A128)
      NA=0
   50 NA=NA+1
      READ(IU,*,ERR=100,END=100) SA(NA),(EM(IW,NA),IW=1,NW)
      IFG=0
      IW=0
      DO WHILE ((IW.LT.NW).AND.(IFG.EQ.0))
        IW=IW+1
        IF (EM(IW,NA).GT.0.) IFG=IW
      ENDDO
      IF (IFG.EQ.0) NA=NA-1
      IF (NA.LT.KNA) GOTO 50
      GOTO 900
C
  100 NA=NA-1
      IF (NA.LE.0) ERC='RDDT10: Data Read error !'
      GOTO 900
C
  800 ERC='RDDT10: Header Read error !'
C
  900 RETURN
      END
C
      SUBROUTINE RDLB10(ALIN,IY,IM,ID,JH,JM,JS,JFLG)
C
C   Read header information for PREDE - skyradiometer POM-01L.
C   for old datafile format( - 1998 ):ITYP=10
C
C --- history
C   2001.02.14 Created by M.Yamano
C
C --- Input
C ALIN       C*128  header
C
C --- Output
C IY,IM,ID   I      date(yyyy/mm/dd)
C JH,JM,JS   I      time(hh:mm:ss)
C JFLG       I      flag of scanning pattern
C                     1: Horizontal scanning
C                     2: Vertical   scanning
C ----------------------------------------------------------------------
C
      CHARACTER ALIN*128,ATIME*8,ADATE*8,HV*1
C
      READ(ALIN,10) ADATE,ATIME,HV
   10 FORMAT(A8,1X,A8,1X,A1)
c     print*,ADATE,ATIME
      READ(ADATE,20) IY,IM,ID
      READ(ATIME,20) JH,JM,JS
   20 FORMAT(I2,1X,I2,1X,I2)
      IF (IY.GE.90) THEN
         IY=IY+1900
      ELSE
         IY=IY+2000
      ENDIF
      IF (HV.EQ.'H') THEN
         JFLG=1
      ELSE
         JFLG=2
      ENDIF
      RETURN
      END
C
      SUBROUTINE CALANG(JFLG,TH0,NA,SCA,TH,FI,IDRCT,MXIDX)
C
C   Calculation of zenith/azimuth angles from scattering angle.
C
C --- history
C   2000.09.22 Created by M.Yamano
C
C --- Input
C JFLG   I       flag of scanning pattern
C                  1: almucantar scanning around the nadir clockwisely.
C                 -1: almucantar scanning around the nadir anti-clockwisely.
C                  2: plane scaning in a plane torward the nadir.
C                 -2: plane scaning in a plane backward to the nadir.
C TH0    R       solar zenith angle in degree
C NA     I       number of scattering angles
C SCA    R(NA)   scattering angles in degree
C
C --- Output
C TH     R(NA)   zenith angles
C FI     R(NA)   azimuth angles relative to the sun
C IDRCT  I       no. of SCA array where SCA(IDRCT)=0.
C MXIDX  I       no. of maximum scattering angle in SCA array
C
C ----------------------------------------------------------------------
C
      PARAMETER (PI=3.141592653590D0,RD=PI/180.0D0)
C
      DIMENSION SCA(NA),TH(NA),FI(NA)
C
      C0=COS(TH0*RD)
      S0=SIN(TH0*RD)
      IDRCT=0
      MXIDX=-1
      DO I=1,NA
        IF (SCA(I).EQ.0.) IDRCT=I
        IF (ABS(JFLG).EQ.1) THEN
           TH(I)=TH0
           IF (I.EQ.IDRCT) THEN
              FI(I)=0.
           ELSE
              IF (SCA(I).LE.TH0*2.) THEN
                 CS=COS(SCA(I)*RD)
                 CF=(CS-C0*C0)/(S0*S0)
                 FI(I)=ACOS(CF)/RD
              ELSE
                 FI(I)=180.
                 IF (MXIDX.LT.0) MXIDX=I-1
              ENDIF
           ENDIF
           IF (JFLG.LT.0) FI(I)=-FI(I)
        ELSE
           FI(I)=0.
           IF (I.EQ.IDRCT) THEN
              TH(I)=TH0
           ELSE
              IF (JFLG.GT.0) THEN
                 TH(I)=TH0-SCA(I)
                 IF (TH(I).LT.0.) THEN
                    TH(I)=ABS(TH(I))
                    IF (TH(I).GE.90.) THEN
                       TH(I)=90.
                       IF (MXIDX.LT.0) MXIDX=I-1
                    ENDIF
                    FI(I)=FI(I)+180.
                 ENDIF
              ELSE
                 TH(I)=TH0+SCA(I)
                 IF (TH(I).GE.90.) THEN
                    TH(I)=90.
                    IF (MXIDX.LT.0) MXIDX=I-1
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
      ENDDO
      IF (MXIDX.LT.0) MXIDX=NA
      RETURN
      END
C
      SUBROUTINE RDDT20(IU,IRFG,ALIN,NW,NA,SA
     &                         ,EM,AP,HP,FH,FV,SH,SV,PT,RL,MFLG,ERC)
C
C   Read measured data for PREDE - skyradiometer POM-01MKII.
C   1st version data format( - 1999 ):ITYP=20
C
C --- history
C   2001.02.14 Created by M.Yamano
C
C --- Input
C IU     I         device No. for reading
C IRFG   I         read option
C                   0: read only header
C                   1: read (header + direct) data
C                   2: read full data(header + direct + diffuse)
C NW     I         number of wavelengths
C
C --- Output
C ALIN   C*128     header
C NA     I         number of scattering angles
C SA     R(NA)     scattering angles
C EM     R(NW,NA)  intensity
C AP     R(NW,NA)  setting azimuth angle
C HP     R(NW,NA)  setting height angle
C FH     R(NW,NA)  horizontal value of fish CCD
C FV     R(NW,NA)  vertical value of fish CCD
C SH     R(NW,NA)  horizontal value of narrow CCD
C SV     R(NW,NA)  vertical value of narrow CCD
C PT     R(NW,NA)  angle of ship pitching in degree
C RL     R(NW,NA)  angle of ship rolling in degree
C MFLG   I         data flag (if MFLG=0 then NORMAL)
C ERC    C*64      ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KNA   =50)
C
      CHARACTER ERC*(*)
      CHARACTER ALIN*128
      DIMENSION SA(KNA),EM(KNW,KNA),AP(KNW,KNA),HP(KNW,KNA)
     &                 ,FH(KNW,KNA),FV(KNW,KNA),SH(KNW,KNA),SV(KNW,KNA)
     &                 ,PT(KNW,KNA),RL(KNW,KNA)
C
      ERC=' '
      MFLG=-9
      NA=0
      READ(IU,10,ERR=800,END=800) ALIN
   10 FORMAT(A128)
      MFLG=0
      IF (IRFG.LT.1) GOTO 900
C
      NA=1
      SA(NA)=0.
      I=0
      DO WHILE((I.LT.NW).AND.(MFLG.EQ.0))
        I=I+1
        READ(IU,*,ERR=30,END=30) NO,EM1,AP1,HP1,FH1,FV1,SH1,SV1,PT1,RL1
        IF (NO.EQ.I) GOTO 40
   30   MFLG=-8
        ERC='RDDT20: Direct data Read error !'
        GOTO 50
C
   40   EM(I,NA)=EM1
        AP(I,NA)=AP1
        HP(I,NA)=HP1
        FH(I,NA)=FH1
        FV(I,NA)=FV1
        SH(I,NA)=SH1
        SV(I,NA)=SV1
        PT(I,NA)=PT1
        RL(I,NA)=RL1
   50 ENDDO
      IF (MFLG.LT.0) GOTO 900
      IF (IRFG.LT.2) GOTO 900
C
   60 READ(IU,*,ERR=900,END=900) SA1
      NA=NA+1
      SA(NA)=SA1
      I=0
      DO WHILE((I.LT.NW).AND.(MFLG.EQ.0))
        I=I+1
        READ(IU,*,ERR=70,END=70) NO,EM1,AP1,HP1,FH1,FV1,SH1,SV1,PT1,RL1
        IF (NO.EQ.I) GOTO 80
   70   MFLG=-7
        ERC='RDDT20: Diffuse data Read error !'
        GOTO 90
C
   80   EM(I,NA)=EM1
        AP(I,NA)=AP1
        HP(I,NA)=HP1
        FH(I,NA)=FH1
        FV(I,NA)=FV1
        SH(I,NA)=SH1
        SV(I,NA)=SV1
        PT(I,NA)=PT1
        RL(I,NA)=RL1
   90 ENDDO
      IF (MFLG.EQ.0) GOTO 60
      NA=NA-1
      GOTO 900
C
  800 ERC='RDDT20: Header Read error !'
C
  900 RETURN
      END
C
      SUBROUTINE RDLB20(ALIN,IY,IM,ID,JH,JM,JS,ALAT,ALNG,JFLG)
C
C   Read header information for PREDE - skyradiometer POM-01MKII.
C   1st version data format( - 1999 ):ITYP=20
C
C --- history
C   2001.02.14 Created by M.Yamano
C
C --- Input
C ALIN       C*128  header
C
C --- Output
C IY,IM,ID   I      date(yyyy/mm/dd)
C JH,JM,JS   I      time(hh:mm:ss)
C ALAT       R      latitude in degree
C ALNG       R      longitude in degree
C JFLG       I      flag of scanning pattern
C                     1: Horizontal scanning
C                     2: Vertical   scanning
C ----------------------------------------------------------------------
C
      CHARACTER ALIN*128,OTHER*100,ATIME*8,ADATE*8,HV*1
C
      NJ=0
      N=0
      DO WHILE (NJ.EQ.0)
        N=N+1
        IF (ALIN(N:N).EQ.'/') NJ=N
      ENDDO
      ADATE=ALIN(NJ-2:NJ+5)
      READ(ALIN,10) ATIME,HV,OTHER
   10 FORMAT(A8,1X,A1,5X,A100)
      READ(OTHER,*) ALNG,ALAT
c     print*,ADATE,ATIME,ALNG,ALAT
      READ(ADATE,20) IY,IM,ID
      READ(ATIME,20) JH,JM,JS
   20 FORMAT(I2,1X,I2,1X,I2)
      IF (IY.GE.98) THEN
         IY=IY+1900
      ELSE
         IY=IY+2000
      ENDIF
      IF (HV.EQ.'H') THEN
         JFLG=1
      ELSE
         JFLG=2
      ENDIF
      RETURN
      END
C
      SUBROUTINE RDDT21(IU,IRFG,ALIN,NW,NA,SA,TM
     &                         ,EM,AP,HP,FH,FV,SH,SV,PT,RL,MFLG,ERC)
C
C   Read measured data for PREDE - skyradiometer POM-01MKII.
C   2nd version data format( 2000 - ):ITYP=21
C
C --- history
C   2001.02.14 Created by M.Yamano
C   2003.03.06 Cut ITYP=22 format.
C
C --- Input
C IU     I         device No. for reading
C IRFG   I         read option
C                   0: read only header
C                   1: read (header + direct) data
C                   2: read full data(header + direct + diffuse)
C NW     I         number of wavelengths
C
C --- Output
C ALIN   C*128     header
C NA     I         number of scattering angles
C SA     R(NA)     scattering angles
C TM     R(NA)     time(UT) in second
C                    if ITYP=21(no time data), set TM=-1
C EM     R(NW,NA)  intensity
C AP     R(NW,NA)  setting azimuth angle
C HP     R(NW,NA)  setting height angle
C FH     R(NW,NA)  horizontal value of fish CCD
C FV     R(NW,NA)  vertical value of fish CCD
C SH     R(NW,NA)  horizontal value of narrow CCD
C SV     R(NW,NA)  vertical value of narrow CCD
C PT     R(NW,NA)  angle of ship pitching in degree
C RL     R(NW,NA)  angle of ship rolling in degree
C MFLG   I         data flag (if MFLG=0 then NORMAL)
C ERC    C*64      ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KNA   =50)
C
      CHARACTER ERC*(*)
      CHARACTER ALIN*128,BLIN*15,ATIME*8
      DIMENSION SA(KNA),TM(KNA),EM(KNW,KNA),AP(KNW,KNA),HP(KNW,KNA)
     &                 ,FH(KNW,KNA),FV(KNW,KNA),SH(KNW,KNA),SV(KNW,KNA)
     &                 ,PT(KNW,KNA),RL(KNW,KNA)
C
      ERC=' '
      MFLG=-9
      NA=0
      READ(IU,10,ERR=800,END=800) ALIN
   10 FORMAT(A128)
      MFLG=0
      IF (IRFG.LT.1) GOTO 900
C
      READ(IU,20,ERR=900,END=900) BLIN
   20 FORMAT(A15)
      NJ=0
      N=0
      DO WHILE((NJ.EQ.0).AND.(N.LT.15))
        N=N+1
        IF (BLIN(N:N).EQ.':') NJ=N
      ENDDO
      IF (NJ.GT.2) THEN
         ATIME=BLIN(NJ-2:NJ+5)
         READ(ATIME,25) JH,JM,JS
   25    FORMAT(I2,1X,I2,1X,I2)
         TM1=JH+JM/60.0+JS/3600.0
      ELSE
         TM1=-1
      ENDIF
      READ(BLIN,*,ERR=900,END=900) SA1
C
      NA=NA+1
      SA(NA)=SA1
      TM(NA)=TM1
      I=0
      DO WHILE((I.LT.NW).AND.(MFLG.EQ.0))
        I=I+1
        READ(IU,*,ERR=30,END=30) NO,EM1,AP1,HP1,FH1,FV1,SH1,SV1,PT1,RL1
        IF (NO.EQ.I) GOTO 40
   30   MFLG=-8
        ERC='RDDT21: Direct data Read error !'
        GOTO 50
C
   40   EM(I,NA)=EM1
        AP(I,NA)=AP1
        HP(I,NA)=HP1
        FH(I,NA)=FH1
        FV(I,NA)=FV1
        SH(I,NA)=SH1
        SV(I,NA)=SV1
        PT(I,NA)=PT1
        RL(I,NA)=RL1
   50 ENDDO
      IF (MFLG.LT.0) GOTO 900
      IF (IRFG.LT.2) GOTO 900
C
   60 READ(IU,20,ERR=900,END=900) BLIN
      NJ=0
      N=0
      DO WHILE((NJ.EQ.0).AND.(N.LT.15))
        N=N+1
        IF (BLIN(N:N).EQ.':') NJ=N
      ENDDO
      IF (NJ.GT.2) THEN
         ATIME=BLIN(NJ-2:NJ+5)
         READ(ATIME,25) JH,JM,JS
         TM1=JH+JM/60.0+JS/3600.0
      ELSE
         TM1=-1
      ENDIF
      READ(BLIN,*,ERR=900,END=900) SA1
C
      NA=NA+1
      SA(NA)=SA1
      TM(NA)=TM1
      I=0
      DO WHILE((I.LT.NW).AND.(MFLG.EQ.0))
        I=I+1
        READ(IU,*,ERR=70,END=70) NO,EM1,AP1,HP1,FH1,FV1,SH1,SV1,PT1,RL1
        IF (NO.EQ.I) GOTO 80
   70   MFLG=-7
        ERC='RDDT21: Diffuse data Read error !'
        GOTO 90
C
   80   EM(I,NA)=EM1
        AP(I,NA)=AP1
        HP(I,NA)=HP1
        FH(I,NA)=FH1
        FV(I,NA)=FV1
        SH(I,NA)=SH1
        SV(I,NA)=SV1
        PT(I,NA)=PT1
        RL(I,NA)=RL1
   90 ENDDO
      IF (MFLG.EQ.0) GOTO 60
      NA=NA-1
      GOTO 900
C
  800 ERC='RDDT21: Header Read error !'
C
  900 RETURN
      END
C
      SUBROUTINE RDLB21(ALIN,IY,IM,ID,JH,JM,JS,ALAT,ALNG,JFLG
     &                      ,IYU,IMU,IDU,JHU,JMU,JSU)
C
C   Read header information for PREDE - skyradiometer POM-01MKII.
C   2nd/3rd version data format ( 2000 - ):ITYP=21,22
C
C --- history
C   2001.02.14 Created by M.Yamano
C
C --- Input
C ALIN       C*128  header
C
C --- Output
C IY ,IM ,ID    I      date(yyyy/mm/dd) - LT
C JH ,JM ,JS    I      time(hh:mm:ss)   - LT
C ALAT          R      latitude in degree
C ALNG          R      longitude in degree
C JFLG          I      flag of scanning pattern
C                       1: Horizontal scanning
C                       2: Vertical   scanning
C IYU,IMU,IDU   I      date(yyyy/mm/dd) - UT
C JHU,JMU,JSU   I      time(hh:mm:ss)   - UT
C ----------------------------------------------------------------------
C
      CHARACTER ALIN*128,OTHER*100,ATIME*8,ADATE*8,HV*1
     &                            ,UTIME*8,UDATE*8
C
      NJ=0
      N=0
      DO WHILE (NJ.EQ.0)
        N=N+1
        IF (ALIN(N:N).EQ.'/') NJ=N
      ENDDO
      ADATE=ALIN(NJ-2:NJ+5)
      UDATE=ALIN(NJ+7:NJ+14)
      UTIME=ALIN(NJ+16:NJ+23)
      READ(ALIN,10) ATIME,HV,OTHER
   10 FORMAT(A8,1X,A1,5X,A100)
      READ(OTHER,*) ALNG,ALAT
c     print*,ADATE,ATIME,ALNG,ALAT
c     print*,UDATE,UTIME
      READ(ADATE,20) IY,IM,ID
      READ(ATIME,20) JH,JM,JS
      READ(UDATE,20) IYU,IMU,IDU
      READ(UTIME,20) JHU,JMU,JSU
   20 FORMAT(I2,1X,I2,1X,I2)
      IF (IY.GE.98) THEN
         IY=IY+1900
      ELSE
         IY=IY+2000
      ENDIF
      IF (IYU.GE.98) THEN
         IYU=IYU+1900
      ELSE
         IYU=IYU+2000
      ENDIF
      IF (HV.EQ.'H') THEN
         JFLG=1
      ELSE
         JFLG=2
      ENDIF
      RETURN
      END
C
      SUBROUTINE RDDT22(IU,IRFG,ALIN,NW,NA,SA,TM
     &                         ,EM,AP,HP,FH,FV,SH,SV,PT,RL,MFLG,ERC)
C
C   Read measured data for PREDE - skyradiometer POM-01MKII.
C   3rd version data format( 2000 - ):ITYP=22
C
C --- history
C   2003.03.06 Created by M.Yamano
C
C --- Input
C IU     I         device No. for reading
C IRFG   I         read option
C                   0: read only header
C                   1: read (header + direct) data
C                   2: read full data(header + direct + diffuse)
C NW     I         number of wavelengths
C
C --- Output
C ALIN   C*128     header
C NA     I         number of scattering angles
C SA     R(NA)     scattering angles
C TM     R(NA)     time(UT) in second
C                    if ITYP=21(no time data), set TM=-1
C EM     R(NW,NA)  intensity
C AP     R(NW,NA)  setting azimuth angle
C HP     R(NW,NA)  setting height angle
C FH     R(NW,NA)  horizontal value of fish CCD
C FV     R(NW,NA)  vertical value of fish CCD
C SH     R(NW,NA)  horizontal value of narrow CCD
C SV     R(NW,NA)  vertical value of narrow CCD
C PT     R(NW,NA)  angle of ship pitching in degree
C RL     R(NW,NA)  angle of ship rolling in degree
C MFLG   I         data flag (if MFLG=0 then NORMAL)
C ERC    C*64      ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KNA   =50)
C
      CHARACTER ERC*(*)
      CHARACTER ALIN*128,ATIME*8
      DIMENSION SA(KNA),TM(KNA),EM(KNW,KNA),AP(KNW,KNA),HP(KNW,KNA)
     &                 ,FH(KNW,KNA),FV(KNW,KNA),SH(KNW,KNA),SV(KNW,KNA)
     &                 ,PT(KNW,KNA),RL(KNW,KNA)
C
      ERC=' '
      MFLG=-9
      NA=0
      READ(IU,10,ERR=800,END=800) ALIN
   10 FORMAT(A128)
      MFLG=0
      IF (IRFG.LT.1) GOTO 900
C
      READ(IU,*,ERR=900,END=900) SA1
C
      NA=NA+1
      SA(NA)=SA1
      NT=0
      UTM=0.
      I=0
      DO WHILE((I.LT.NW).AND.(MFLG.EQ.0))
        I=I+1
        READ(IU,*,ERR=30,END=30) NO,EM1,AP1,HP1,ATIME
     &                          ,FH1,FV1,SH1,SV1,PT1,RL1
        IF (NO.EQ.I) GOTO 40
   30   MFLG=-8
        ERC='RDDT22: Direct data Read error !'
        GOTO 50
C
   40   READ(ATIME,45) IH,IM,IS
   45   FORMAT(I2,1X,I2,1X,I2)
        NT=NT+1
        UTM=UTM+IH+IM/60.0+IS/3600.0
        EM(I,NA)=EM1
        AP(I,NA)=AP1
        HP(I,NA)=HP1
        FH(I,NA)=FH1
        FV(I,NA)=FV1
        SH(I,NA)=SH1
        SV(I,NA)=SV1
        PT(I,NA)=PT1
        RL(I,NA)=RL1
   50 ENDDO
      IF (NT.GT.0) THEN
         TM(NA)=UTM/FLOAT(NT)
      ELSE
         TM(NA)=-1
      ENDIF
      IF (MFLG.LT.0) GOTO 900
      IF (IRFG.LT.2) GOTO 900
C
   60 READ(IU,*,ERR=900,END=900) SA1
C
      NA=NA+1
      SA(NA)=SA1
      NT=0
      UTM=0.
      I=0
      DO WHILE((I.LT.NW).AND.(MFLG.EQ.0))
        I=I+1
        READ(IU,*,ERR=70,END=70) NO,EM1,AP1,HP1,ATIME
     &                          ,FH1,FV1,SH1,SV1,PT1,RL1
        IF (NO.EQ.I) GOTO 80
   70   MFLG=-7
        ERC='RDDT22: Diffuse data Read error !'
        GOTO 90
C
   80   READ(ATIME,45) IH,IM,IS
        NT=NT+1
        UTM=UTM+IH+IM/60.0+IS/3600.0
        EM(I,NA)=EM1
        AP(I,NA)=AP1
        HP(I,NA)=HP1
        FH(I,NA)=FH1
        FV(I,NA)=FV1
        SH(I,NA)=SH1
        SV(I,NA)=SV1
        PT(I,NA)=PT1
        RL(I,NA)=RL1
   90 ENDDO
      IF (NT.GT.0) THEN
         TM(NA)=UTM/FLOAT(NT)
      ELSE
         TM(NA)=-1
      ENDIF
      IF (MFLG.EQ.0) GOTO 60
      NA=NA-1
      GOTO 900
C
  800 ERC='RDDT22: Header Read error !'
C
  900 RETURN
      END
C
      SUBROUTINE PSFISH(FH,FV,AP,FH0,FV0,ASPF,CO,DFI,FC,FI,HT)
C
C   Get solar direction from detected values of fish CCD.
C
C --- history
C   2000.08.25 Created by M.Yamano
C   2001.10.06 Renewed.
C   2002.10.02 Modified.
C   2002.11.21 Modified.
C   2002 11.28 Debugged.
C
C --- Input
C FH,FV         R         detected values of fish CCD
C AP            R         azimuth angle driven by motor
C FH0,FV0       R         fish CCD parameters
C CO,ASPF,FC(5) R
C DFI           R
C
C --- Output
C FI,HT         R         azimuth/height angles of solar direction
C ----------------------------------------------------------------------
C
      PARAMETER (PI=3.141592653590D0,RD=PI/180.0D0)
C
      DIMENSION FC(5)
C
      DATA DCCD/0.2/
C
      IF ((FH.LT.DCCD).OR.(FV.LT.DCCD)) THEN
         FI=-999.
         HT=-99.
      ELSE
         X=-(FV-FV0)/ASPF
         Y=FH-FH0
         R=SQRT(X*X+Y*Y)
         HT=HGT4(CO*R,FC)
         FI=ATAN2(X,Y)/RD-DFI+AP
         FI=DEG(FI,0)
      ENDIF
      RETURN
      END
C
      SUBROUTINE PSNARW(SH,SV,AP,HP,SH0,SV0,ASPN,CSA,DZ,SA,FI,HT)
C
C   Get solar direction from detected values of narrow CCD.
C
C --- history
C   2000.08.18 Created by M.Yamano
C   2001.10.09 Renewed.
C   2002.11.01 Modified.
C   2002.11.21 Modified.
C
C --- Input
C SH,SV         R         detected value of narrow CCD
C AP,HP         R         azimuth/height angles driven by motor
C SH0,SV0       R         narrow CCD parameters
C CSA,ASPN,DZ   R
C
C --- Output
C SA            R         scattering angle in degree
C FI,HT         R         azimuth/height angles of solar direction
C ----------------------------------------------------------------------
C
      PARAMETER (PI=3.141592653590D0,RD=PI/180.0D0)
C
      DATA DCCD,SAMX,DDG /0.2, 15.0, 0.15/
C
      IF ((SH.LT.DCCD).OR.(SV.LT.DCCD)) THEN
         SA=-99.
         FI=-999.
         HT=-99.
      ELSE
         X1=-(SH-SH0)
         Y1=-(SV-SV0)/ASPN
         CALL ROTATU(-DZ,X1,Y1,X2,Y2)
         R=SQRT(X2*X2+Y2*Y2)
         SA=ATAN(R/CSA)/RD
         IF (SA.GT.SAMX) THEN
            SA=-99.
            FI=-999.
            HT=-99.
         ELSE
            DF= ATAN(X2/CSA)/RD
            DH=-ATAN(Y2/CSA)/RD
            CALL STNARW(DF,DH,AP,HP,FI,HT)
            IF (HP.LE.90.) THEN
               HP1=HP
               AP1=AP
            ELSE
               HP1=180.-HP
               AP1=DEG(AP+180.,0)
            ENDIF
            SAC=SCAANG(AP1,90.-HP1,FI,90.-HT)
            IF (ABS(SA-SAC).GT.DDG) THEN
               WRITE(*,10) 'PSNARW :',SA,SAC,SA-SAC,AP,HP,FI,HT
   10          FORMAT(A8,3F9.3,4F9.2)
               SA=-99.
               FI=-999.
               HT=-99.
            ENDIF
         ENDIF
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
      SUBROUTINE GETXYZ(FI,HT,X,Y,Z)
C
C   Get coordinates(X,Y,Z) from (FI,HT).
C
C --- history
C   2000.08.18 Created by M.Yamano
C
C --- Input
C FI,HT         R         azimuth/height angles in degree
C
C --- Output
C X,Y,Z         R         (X,Y,Z) coordinates value
C ----------------------------------------------------------------------
C
      PARAMETER (PI=3.141592653590D0,RD=PI/180.0D0)
C
      X=COS(HT*RD)*SIN(FI*RD)
      Y=COS(HT*RD)*COS(FI*RD)
      Z=SIN(HT*RD)
      RETURN
      END
C
      SUBROUTINE GETPOL(X,Y,Z,FI,HT)
C
C   Get polar coordinate value(FI,HT) from (X,Y,Z).
C
C --- history
C   2000.08.18 Created by M.Yamano
C
C --- Input
C X,Y,Z         R         (X,Y,Z) coordinates value
C
C --- Output
C FI,HT         R         azimuth/height angles in polar coordinate
C ----------------------------------------------------------------------
C
      PARAMETER (PI=3.141592653590D0,RD=PI/180.0D0)
C
      IF (ABS(Z).GE.1) THEN
         Z=SGN(Z)*1.
         FI=0.
      ELSE
         FI=ATAN2(X,Y)/RD
         IF (FI.LT.0.) FI=FI+360.
      ENDIF
      HT=ASIN(Z)/RD
      RETURN
      END
C
      SUBROUTINE ROTATU(ANG,RY,RZ,WY,WZ)
C
C   Rotation of coordinates.
C
C --- history
C   2000.08.18 Created by M.Yamano
C
C --- Input
C ANG           R         rotational angle in degree
C RY,RZ         R         coordinates before rotation
C
C --- Output
C WY,WZ         R         coordinates after rotation
C ----------------------------------------------------------------------
C
      PARAMETER (PI=3.141592653590D0,RD=PI/180.0D0)
C
      WY= COS(ANG*RD)*RY+SIN(ANG*RD)*RZ
      WZ=-SIN(ANG*RD)*RY+COS(ANG*RD)*RZ
      RETURN
      END
C
      SUBROUTINE CLDRFT(CX,CY,RX,RY,FG)
C
C   Get rotational angle of coordinates.
C
C --- history
C   2000.08.18 Created by M.Yamano
C   2003.04.02 Changed FG: [##] -> [##.#] deg.
C
C --- Input
C CX,CY         R         coordinate values before roration
C RX,RY         R         coordinate values after roration
C
C --- Output
C FG            R         rorational angle in degree
C ----------------------------------------------------------------------
C
      PARAMETER (PI=3.141592653590D0,RD=PI/180.0D0)
C
      IR=0
      RRMIN=10000.
      DO I = 1,36
        FG=FLOAT(I)*10.
        CALL ROTATU(FG,CX,CY,X1,Y1)
        RR=(X1-RX)*(X1-RX)+(Y1-RY)*(Y1-RY)
        IF (RR.LE.RRMIN) THEN
           IR=I*10
           RRMIN=RR
        ENDIF
      ENDDO
      IL=IR
      DO I = IL-10,IL+10
        FG=FLOAT(I)
        CALL ROTATU(FG,CX,CY,X1,Y1)
        RR=(X1-RX)*(X1-RX)+(Y1-RY)*(Y1-RY)
        IF (RR.LE.RRMIN) THEN
           IR=I
           RRMIN=RR
        ENDIF
      ENDDO
      FG0=FLOAT(IR)
      FGR=FG0
      DO I = 1,21
        FG=FG0+0.1*FLOAT(I-11)
        CALL ROTATU(FG,CX,CY,X1,Y1)
        RR=(X1-RX)*(X1-RX)+(Y1-RY)*(Y1-RY)
        IF (RR.LE.RRMIN) THEN
           FGR=FG
           RRMIN=RR
        ENDIF
      ENDDO
      FG=DEG(FGR,0)
      RETURN
      END
C
      SUBROUTINE RDOBSP(IU,CMNT,SITE,CNTY
     &                 ,ALNGS,ALNG,ALAT,ALT,NA,SA,ERC)
C
C   Read observational parameters.(from 'obs.para' file)
C
C --- History
C   2001.02.14 Created by M.Yamano
C   2002.05.11 Renewed by M.Yamano for Ver.5.
C   2002.11.20 Modified.
C
C --- Input
C IU     I       device No. for reading
C
C --- Output
C CMNT   C*40    project name
C SITE   C*20    observation site (or ship) name
C CNTY   C*20    country name
C ALNGS  R       standard longitude in degree for the time zone
C ALNG   R       longitude of observational site in degree
C ALAT   R       latitude of observational site in degree
C ALT    R       height of observational site in m.
C NA     I       number of observational scattering angles
C SA     R(NA)   scattering angles in degree
C ERC    C*64    ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNA   =50)
C
      CHARACTER ERC*(*)
      CHARACTER CMNT*(*),SITE*(*),CNTY*(*)
      DIMENSION SA(KNA)
C
      ERC=' '
      READ(IU,*,ERR=800,END=800)
      READ(IU,*,ERR=800,END=800) CMNT
      READ(IU,*,ERR=800,END=800) SITE
      READ(IU,*,ERR=800,END=800) CNTY
C
      READ(IU,*,ERR=800,END=800) ALNGS
      READ(IU,*,ERR=800,END=800) ALNG,ALAT,ALT
      READ(IU,*,ERR=800,END=800) NA
      IF (NA.GT.KNA) THEN
         ERC='RDOBSP: NA.GT.KNA !'
         GOTO 900
      ENDIF
      READ(IU,*,ERR=800,END=800) (SA(I),I=1,NA)
      GOTO 900
C
  800 ERC='RDOBSP: Read Error !'
C
  900 RETURN
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
      SUBROUTINE RDCCDP(IU,PT0,RL0,DHT,DFAI,SH0,SV0,ASPN,CSA,DZ
     &                    ,FH0,FV0,ASPF,CO,DFI,FC,ERC)
C
C   Read CCD calibration parameters.(from 'CCD.para' file)
C
C --- history
C   2000.11.29 Created by M.Yamano
C   2001.09.29 Renewed.
C   2002.11.20 Changed format.
C
C --- Input
C IU            I         device No. for reading
C
C --- Output
C PT0,RL0       R         pitching/rolling angle offset
C DHT,DFAI      R         height offset(=Hps-Hs)/azimuth offset(=Aps-As)
C SH0,SV0,ASPN  R         narrow CCD parameters
C CSA,DZ
C FH0,FV0,ASPF  R         fish CCD parameters
C CO,DFI,FC(5)
C ERC           C*64      ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
C
      CHARACTER ERC*(*)
      DIMENSION FC(5)
C
      ERC=' '
      READ(IU,*,ERR=800,END=800)
      READ(IU,*,ERR=800,END=800) PT0,RL0
      READ(IU,*,ERR=800,END=800) DHT,DFAI
      READ(IU,*,ERR=800,END=800) SH0,SV0,ASPN,CSA,DZ
      READ(IU,*,ERR=800,END=800) FH0,FV0,ASPF,CO,DFI
      READ(IU,*,ERR=800,END=800) (FC(I),I=1,5)
      GOTO 900
C
  800 ERC='RDCCDP: Read error.'
C
  900 RETURN
      END
C
      FUNCTION HGT4(R,CR)
C
C   Calculate height of detected point in fish CCD.
C
C --- history
C   2000.08.18 Created by M.Yamano
C
C --- Input
C R             R         radius of detected point in fish CCD
C CR            R(5)      coefficients of 4-order fitting function
C
C --- Output
C HGT4          R         height in degree of detected point in fish CCD
C ----------------------------------------------------------------------
C
      DIMENSION CR(5)
C
      HGT4=CR(5)*R**4+CR(4)*R**3+CR(3)*R**2+CR(2)*R+CR(1)
      RETURN
      END
C
      SUBROUTINE STNARW(FDG,HDG,AS,HS,AI,HI)
C
C   Get sensor direction of narrow CCD.
C
C --- history
C   2000.09.11 Created by M.Yamano
C
C --- Input
C FDG           R         azimuthal deviation from solar direction
C HDG           R         zenithal deviation from solar direction
C AS            R         azimuth angle of solar direction
C HS            R         height angle of solar direction
C
C --- Output
C AI            R         azimuth angle of sensor direction
C HI            R         height angle of sensor direction
C ----------------------------------------------------------------------
C
      CALL GETXYZ(FDG,HDG,XM,YM,ZM)
      CALL ROTATU(-HS,YM,ZM,Y1,Z1)
      CALL ROTATU( AS,XM,Y1,X2,Y2)
      CALL GETPOL(X2,Y2,Z1,AI,HI)
      RETURN
      END
C
       FUNCTION SGN(X)
C
C   Get plus/minus sign of value X.
C
C --- history
C   2000.08.18 Created by M.Yamano
C
C --- Input
C X             R         variable
C
C --- Output
C SGN           R          1: plus value
C                          0: zero value
C                         -1: minus value
C ----------------------------------------------------------------------
C
      SGN=0.
      IF (X.GT.0.) SGN=1.
      IF (X.LT.0.) SGN=-1.
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
