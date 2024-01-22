C$OUTF solid3m.f
C$LIBF /Users/makibo/Fortran/LIB/SKYLIB/A
C$LIBF /Users/makibo/Fortran/LIB/POMLIB/A
C$LIBF /Users/makibo/Fortran/LIB/BASLIB/A
C$LIBF /Users/makibo/Fortran/LIB_Nakajima/LBC1/A
C$LIBF /Users/makibo/Fortran/LIB_Nakajima/LBM1/A
C$LIBF /Users/makibo/Fortran/LIB_Nakajima/LBR1/A
C
C --- Note
C   Skyradiometer data processing program
C   For PREDE instruments POM-01L/POM-01MKII
C
C   Calculation of solid view angles.
C
C   * Define the scanning is from down to up, from left to right
C   * Get z=f(th)
C   * (x, y+vx)=(x0,y0)+Rot(ang) (r cos th, br sin th)
C   * data are taken with solar movement correction.
C
C --- History
C 95. 2.17 Created from solid3
C 96. 3. 5 Kazuma Aoki
C          Created from solid3b
C          A(IZ,IA)-DMIN
C 96. 3.10 Created from solid3b
C          Kazuma Aoki
C 96. 9.26 Created from solid3m.f
C          Kazuma Aoki
C *******************************************************************
C 01.01.26 Overhauled by M.Yamano for skyrad.pack.
C *******************************************************************
C 01.07.06 Renewed by M.Yamano
C 01.08.17 Modified by M.Yamano
C 01.10.06 'SVA.out' is added.
C 02.02.19 Modified by M.Yamano
C 03.03.20 Renewed.
C 03.03.26 Changed file name -> (yy)yymmdd.V??
C 03.11.05 Added processing for ITYP=30 (new format) data.
C 03.11.18 Modified.
C 04.12.06 Start on Mac.
C 06.03.14 Changed control for S/N check in case of ITYP=30
C
C --- I/O files
C  solid.par           I    processing parameters file
C  ins.para            I    instrument parameters file(any name is OK)
C  fname               I    processing dates file((yy)yymmdd or DISK?.DAT)
C  DSK/yymmddnn.V??    I    measured disk-scan data file (ITYP=30)
C  or  (yy)yymmdd.V??       measured disk-scan data file (ITYP<30)
C  or  DISK?.DAT
C  or  disk?.dat
C  disk??.out          O    solid view angles results
C  SVA.out             O    merged file of all disk?.out files
C  disk??.log          O    processing log file
C  disk??.map          O    disk-scan map data(if IMAP>0)
C  disk??.fnc          O    response function data(if IFNC>0)
C  disk??.lst          O    measured data listfile(if ILST>0)
C  solid3m.tmp         W    work file
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KDAY  =50)
      PARAMETER (KN    =21)
      PARAMETER (KCH   =47)
      PARAMETER (KN2   =201)
      PARAMETER (KNN=KN*KN)
C
      PARAMETER (PI=3.141592653589793,RAD=PI/180.0)
C
      CHARACTER ERC*64
      CHARACTER ALINE*255,ALIN0*255,FNM*80,CCH*3,MRK*1
      CHARACTER DFIL*80,OFIL*80,RFIL*80,SFIL*80,CFIL*80,LFIL*80,WFIL*80
      CHARACTER SRNO*20,SRNOO*20,SRNOD*7
      DIMENSION WL(KNW),WLO(KNW)
      DIMENSION DX(KN,KN),DY(KN,KN),AM(KN,KN),RM(KN,KN),ZM(KN,KN)
     &         ,RRR(KN2),ZZZ(KN2)
C
      DATA DR,AERR /0.02, 0.2/
      DATA SRNOD /'0000000'/
C
      ERC=' '
      IUI=1
      OPEN (IUI,FILE='solid.par',STATUS='OLD')
      CALL RDSLDP(IUI,IMAP,IFNC,ILST,SRNOO,ITYP,NZ,NA,DA,NCH,WLO,ERC)
      CLOSE (IUI)
      IF (ERC.NE.' ') GOTO 900
C
      IF (ITYP.LT.20) THEN
         IRK=5
      ELSE
         IRK=3
      ENDIF
      NDR=2.5/DR+1
      DO I=1,NDR
        RRR(I)=(I-1)*DR
      ENDDO
C
C --- processing
C
      IUO=2
      IUF=3
      IUR=4
      IUS=7
      IUC=8
      IUL=9
      IUW=10
      WFIL='solid3m.tmp'
      IUT=12
      OPEN (IUT,FILE='SVA.out',STATUS='UNKNOWN')
      DO ICH=1,NCH
        WRITE(CCH,10) ICH+100
   10   FORMAT(I3)
        IF (ICH.LT.10) THEN
           IST=3
        ELSE
           IST=2
        ENDIF
        IFIRST=0
        OPEN (IUF,FILE='fname',STATUS='OLD')
        NDT=0
        NAV=0
        SVAA=0.
C
   20   READ(IUF,25,END=200,ERR=200) FNM
   25   FORMAT(A80)
        CALL GETPOS(FNM,'.',80,IL)
        IF (IL.EQ.1) GOTO 200
        IF (IL.GT.0) FNM=FNM(1:IL-1)
        CALL UTLSPC(FNM,NFNM)
        IF (ITYP.GE.30) THEN
           DFIL='DSK/'//FNM(1:NFNM)//'.V'//CCH(2:3)
        ELSE
           IF (FNM(1:4).EQ.'DISK') THEN
              DFIL='DSK/DISK'//CCH(IST:3)//'.DAT'
           ELSEIF (FNM(1:4).EQ.'disk') THEN
              DFIL='DSK/disk'//CCH(IST:3)//'.dat'
           ELSE
              DFIL='DSK/'//FNM(1:NFNM)//'.V'//CCH(2:3)
           ENDIF
        ENDIF
        IFG=0
        OPEN (IUI,FILE=DFIL,STATUS='OLD',ERR=100)
        IF (ITYP.GE.30) THEN
           CALL RDHD30(IUI,SRNO,ALNGS,ALNG,ALAT,NW,WL,ERC)
           IF (ERC.EQ.' ') THEN
              WL0=WLO(ICH)*1.E7
              IF (NCH.EQ.NW) THEN
                 IF (ABS(WL0-WL(ICH)).GE.1.) ERC='Wavelength mismatch !'
              ELSE
                 NFG=0
                 N=0
                 DO WHILE ((NFG.EQ.0).AND.(N.LT.NW))
                   N=N+1
                   IF (ABS(WL0-WL(N)).LT.1.) NFG=1
                 ENDDO
                 IF (NFG.EQ.0) ERC='Wavelength mismatch !'
              ENDIF
c             IF (NCH.NE.NW) ERC='Number of wavelength mismatch !'
              IF (SRNO.NE.SRNOO) THEN
                 IF (SRNO(1:7).EQ.SRNOD) THEN
c                   WRITE(*,27) FNM(1:NFNM),CCH(2:3)
c  27               FORMAT('No check for S/N! in file (',A8,'.V',A2,')')
                 ELSE
                    ERC='S/N mismatch !'
                 ENDIF
              ENDIF
           ENDIF
           IF (ERC.NE.' ') THEN
              WRITE(*,28) FNM(1:NFNM),CCH(2:3),ERC
   28         FORMAT('Skip file (',A8,'.V',A2,') : ',A64)
              ERC=' '
              GOTO 20
           ENDIF
        ENDIF
C
   30   READ(IUI,35,ERR=100,END=100) ALINE
   35   FORMAT(A255)
        IF (ALINE(1:1).EQ.' ') ALINE=ALINE(2:255)
        IF (ALINE(1:10).EQ.'          ') GOTO 30
        IF (ALINE(IRK:IRK).EQ.'/') THEN
           CLOSE(IUW)
C
           IF (IFG.GT.0) THEN
              OPEN (IUW,FILE=WFIL,STATUS='UNKNOWN')
              CALL PRCDSK(IUW,ITYP,IY,IM,ID,TIM
     &                   ,NZ,NA,DA,DX,DY,AM,BX1,BY1,RAA,BET,SVA,RM,ZM
     &                   ,NDR,RRR,ZZZ,EROR,ERC)
              CLOSE(IUW)
              WRITE(IUR,40) IY,IM,ID,TIM,ERC
   40         FORMAT(I4,2I3,F6.2,1X,A63)
              IF (ERC.EQ.' ') THEN
                 IF (EROR.LT.AERR) THEN
                    NAV=NAV+1
                    SVAA=SVAA+SVA
                    MRK='*'
                 ELSE
                    MRK=' '
                 ENDIF
                 CALL WTDSKP(IUO,NDT,IY,IM,ID,TIM,BX1,BY1,RAA,BET,SVA
     &                      ,EROR,MRK)
                 CALL WTDSKP(IUT,NDT,IY,IM,ID,TIM,BX1,BY1,RAA,BET,SVA
     &                      ,EROR,MRK)
                 IF (IMAP.GT.0)
     &              CALL WTDSKD(IUS,NDT,IY,IM,ID,TIM,BX1,BY1,RAA,BET,SVA
     &                         ,EROR,MRK,NZ,NA,DX,DY,RM,AM,ZM,ZZZ(1))
                 IF (IFNC.GT.0)
     &              CALL WTDSKF(IUC,NDT,IY,IM,ID,TIM,BX1,BY1,RAA,BET,SVA
     &                         ,EROR,MRK,NDR,RRR,ZZZ)
              ENDIF
           ENDIF
C
           NDT=NDT+1
           ALIN0=ALINE
           OPEN (IUW,FILE=WFIL,STATUS='UNKNOWN')
           WRITE(IUW,35) ALIN0
           IF (IFIRST.EQ.0) THEN
              RFIL='disk'//CCH(2:3)//'.log'
              OPEN (IUR,FILE=RFIL,STATUS='UNKNOWN')
              OFIL='disk'//CCH(2:3)//'.out'
              OPEN (IUO,FILE=OFIL,STATUS='UNKNOWN')
              WRITE(IUO,43) ICH,WLO(ICH)*1E4
              WRITE(IUT,43) ICH,WLO(ICH)*1E4
   43         FORMAT('CH = ',I2,'   WL = ',F7.3
     &             //' No YYYY MM DD   HR      X0      Y0'
     &              ,'    Lambda   Beta    SVA        error')
              IF (IMAP.GT.0) THEN
                 SFIL='disk'//CCH(2:3)//'.map'
                 OPEN (IUS,FILE=SFIL,STATUS='UNKNOWN')
              ENDIF
              IF (IFNC.GT.0) THEN
                 CFIL='disk'//CCH(2:3)//'.fnc'
                 OPEN (IUC,FILE=CFIL,STATUS='UNKNOWN')
              ENDIF
              IF (ILST.GT.0) THEN
                 LFIL='disk'//CCH(2:3)//'.lst'
                 OPEN (IUL,FILE=LFIL,STATUS='UNKNOWN')
              ENDIF
              IFIRST=1
           ENDIF
           IF (ILST.GT.0) WRITE(IUL,45) NDT,ALIN0
   45      FORMAT(I3,2X,A72)
           IFG=1
c          print*,ALIN0
        ELSE
           WRITE(IUW,35) ALINE
        ENDIF
        GOTO 30
C
  100   CLOSE(IUW)
        IF (IFG.GT.0) THEN
           OPEN (IUW,FILE=WFIL,STATUS='UNKNOWN')
           CALL PRCDSK(IUW,ITYP,IY,IM,ID,TIM
     &                ,NZ,NA,DA,DX,DY,AM,BX1,BY1,RAA,BET,SVA,RM,ZM
     &                ,NDR,RRR,ZZZ,EROR,ERC)
           CLOSE(IUW)
           WRITE(IUR,40) IY,IM,ID,TIM,ERC
           IF (ERC.EQ.' ') THEN
              IF (EROR.LT.AERR) THEN
                 NAV=NAV+1
                 SVAA=SVAA+SVA
                 MRK='*'
              ELSE
                 MRK=' '
              ENDIF
              CALL WTDSKP(IUO,NDT,IY,IM,ID,TIM,BX1,BY1,RAA,BET,SVA
     &                   ,EROR,MRK)
              CALL WTDSKP(IUT,NDT,IY,IM,ID,TIM,BX1,BY1,RAA,BET,SVA
     &                   ,EROR,MRK)
              IF (IMAP.GT.0)
     &           CALL WTDSKD(IUS,NDT,IY,IM,ID,TIM,BX1,BY1,RAA,BET,SVA
     &                      ,EROR,MRK,NZ,NA,DX,DY,RM,AM,ZM,ZZZ(1))
              IF (IFNC.GT.0)
     &           CALL WTDSKF(IUC,NDT,IY,IM,ID,TIM,BX1,BY1,RAA,BET,SVA
     &                      ,EROR,MRK,NDR,RRR,ZZZ)
           ENDIF
        ENDIF
        CLOSE(IUI)
        GOTO 20
C
  200   CLOSE(IUF)
        IF (IFIRST.GT.0) THEN
           IF (NAV.GT.0) THEN
              SVAA=SVAA/FLOAT(NAV)
              WRITE(IUO,210) AERR,SVAA,NAV
              WRITE(IUT,210) AERR,SVAA,NAV
  210         FORMAT(/33X,'ave.( error <',0PF4.1,' ) =',1PE10.3
     &                ,9X,'[',I3,']')
           ELSE
              WRITE(IUO,220) AERR,NAV
              WRITE(IUT,220) AERR,NAV
  220         FORMAT(/33X,'ave.( error <',0PF4.1,' ) = *********'
     &                ,9X,'[',I3,']')
           ENDIF
        ENDIF
        CLOSE(IUO)
        CLOSE(IUR)
        CLOSE(IUS)
        CLOSE(IUC)
        CLOSE(IUL)
      ENDDO
C
      CLOSE(IUT)
      OPEN (IUW,FILE=WFIL,STATUS='UNKNOWN')
      CLOSE (IUW,STATUS='DELETE')
C
  900 IF (ERC.NE.' ') WRITE(*,*) 'ERROR CODE: ', ERC
C
 1000 STOP
      END
C
      SUBROUTINE PRCDSK(IU,ITYP,IY,IM,ID,TIM
     &                 ,NZ,NA,DA,DX,DY,AM,BX1,BY1,RAA,BET,SVA,RM,ZM
     &                 ,NDR,RRR,ZZZ,EROR,ERC)
      PARAMETER (KN    =21)
      PARAMETER (KN2   =201)
      PARAMETER (KNN=KN*KN)
C
      PARAMETER (PI=3.141592653589793,RAD=PI/180.0)
C
      CHARACTER ERC*(*)
      DIMENSION AO(KN,KN)
      DIMENSION DX(KN,KN),DY(KN,KN),AM(KN,KN),RM(KN,KN),ZM(KN,KN)
     &         ,RRR(KN2),ZZZ(KN2)
C
      ERC=' '
      CALL RDDISK(IU,ITYP,IY,IM,ID,JH,JM,JS,NZ,NA,AO,ERC)
      IF (ERC.NE.' ') GOTO 900
C
      TIM=JH+JM/60.+JS/3600.
C
      DMIN=AO(1,1)
      DO IA=1,NA
        DO IZ=1,NZ
          IF (AO(IZ,IA).LT.DMIN) DMIN=AO(IZ,IA)
        ENDDO
      ENDDO
      DO IA=1,NA
        DO IZ=1,NZ
          AM(IZ,IA)=AO(IZ,IA)-DMIN
          DX(IZ,IA)=DA*(IA-(NA+1)/2)
          DY(IZ,IA)=-DA*(IZ-(NZ+1)/2)
        ENDDO
      ENDDO
      CALL CALSVA(NZ,NA,DX,DY,AM,BX1,BY1,RAA,BET,SVA,RM,ZM
     &                 ,NDR,RRR,ZZZ,EROR,ERC)
C
  900 RETURN
      END
C
      SUBROUTINE WTDSKP(IU,NDT,IY,IM,ID,TIM,BX1,BY1,RAA,BET,SVA
     &                 ,EROR,MRK)
C
      CHARACTER MRK*1
C
      WRITE(IU,10) NDT,IY,IM,ID,TIM,BX1,BY1,RAA,BET,SVA,EROR,MRK
   10 FORMAT(I3,I5,2I3,F7.3,2F8.3,F8.2,F8.4,1P2E11.3,1X,A1)
      RETURN
      END
C
      SUBROUTINE WTDSKD(IU,NDT,IY,IM,ID,TIM,BX1,BY1,RAA,BET,SVA
     &                    ,EROR,MRK,NZ,NA,DX,DY,RM,AM,ZM,ZZZ1)
      PARAMETER (KN    =21)
C
      CHARACTER MRK*1
      DIMENSION DX(KN,KN),DY(KN,KN),AM(KN,KN),RM(KN,KN),ZM(KN,KN)
C
      CALL WTDSKP(IU,NDT,IY,IM,ID,TIM,BX1,BY1,RAA,BET,SVA,EROR,MRK)
      DO IA=1,NA
        DO IZ=1,NZ
          WRITE(IU,10) DX(IZ,IA)-BX1,DY(IZ,IA)-BY1
     &                ,RM(IZ,IA),AM(IZ,IA)/ZZZ1,ZM(IZ,IA)
   10     FORMAT(2F8.3,1P3E12.4)
        ENDDO
      ENDDO
      RETURN
      END
C
      SUBROUTINE WTDSKF(IU,NDT,IY,IM,ID,TIM,BX1,BY1,RAA,BET,SVA
     &                    ,EROR,MRK,NDR,RRR,ZZZ)
      PARAMETER (KN2   =201)
C
      CHARACTER MRK*1
      DIMENSION RRR(KN2),ZZZ(KN2)
C
      CALL WTDSKP(IU,NDT,IY,IM,ID,TIM,BX1,BY1,RAA,BET,SVA,EROR,MRK)
      DO I=1,NDR
        WRITE(IU,10) RRR(I),ZZZ(I)/ZZZ(1)
   10   FORMAT(F6.2,1PE11.3)
      ENDDO
      RETURN
      END
C
      SUBROUTINE RDSLDP(IU,IMAP,IFNC,ILST,SRNO,ITYP,NZ,NA,DA,NW,WL,ERC)
C
C   Read parameters from file('solid.par')
C
C --- history
C   2001.02.16 Created by M.Yamano
C   2001.05.23 Renewed. DA is added and output options are changed
C   2002.02.19 Modified.
C   2002.11.20 Modified.
C   2003.01.08 Renewd.
C   2003.11.05 Added SRNO to output arguments.
C
C --- Input
C IU       I       device No. for reading
C
C --- Output
C IMAP     I       disk-scan map data output option
C IFNC     I       respose function outout option
C ILST     I       measured data list output option
C SRNO     C*20    instrument S/N
C ITYP     I       PREDE - skyradiometer instrument type
C                    10,11   : POM-01L
C                    20,21,22: POM-01MKII
C NZ       I       number of grids in zenithal direction
C NA       I       number of grids in azimuthal direction
C DA       I       angle in degree between neighboring grids
C NW       I       number of wavelengths
C WL       R(NW)   wavelengths in cm
C ERC      C*64    ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KNW   =20)
      PARAMETER (KDAY  =50)
C
      CHARACTER ERC*(*)
      CHARACTER PFIL*80
      CHARACTER SRNO*20
      DIMENSION WL(KNW),SOLID(KNW),F0D(KNW,KDAY)
C
      ERC=' '
      READ(IU,*,ERR=800,END=800) IMAP,IFNC,ILST
      READ(IU,*,ERR=800,END=800) NZ,NA,DA
C
C --- read 'ins.para' parameters
C
      READ(IU,*,ERR=800,END=800) PFIL
      CALL UTLSPC(PFIL,NL)
      IF (NL.LE.0) THEN
         ERC='RDSLDP: illegal ins.para filename'
         GOTO 900
      ENDIF
      OPEN (IU+1,FILE=PFIL(1:NL),STATUS='OLD')
      CALL RDINSP(IU+1,SRNO,ITYP,NW,WL,SOLID,NDY,F0D,TJLC,ERC)
      CLOSE (IU+1)
      GOTO 900
C
  800 ERC='RDSLDP: Read Error !'
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
      SUBROUTINE RDDISK(IU,ITYP,IY,IM,ID,JH,JM,JS,NZ,NA,DAT,ERC)
C
C   Read disk-scan data for PREDE - skyradiometer.
C
C --- history
C   2001.02.16 Created by M.Yamano
C   2003.11.05 Added processing for ITYP=30 (new format) data.
C
C --- Input
C IU         I           device No. for reading
C ITYP       I           PREDE - skyradiometer instrument type
C                          10,11   : POM-01L
C                          20,21,22: POM-01MKII
C                          30      : POM-02 (Ver.4.10)
C NZ         I           number of grids in zenithal direction
C NA         I           number of grids in azimuthal direction
C
C --- Output
C IY,IM,ID   I           date(yyyy/mm/dd) - LT
C JH,JM,JS   I           time(hh:mm:ss)   - LT
C DAT        R(NZ,NA)    disk-scan data array
C ERC        C*64        ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KN    =21)
C
      CHARACTER ERC*(*)
      CHARACTER ALIN*128,ADATE*10,ATIME*8
      DIMENSION DAT(KN,KN)
C
      ERC=' '
      READ(IU,10,ERR=800,END=800) ALIN
   10 FORMAT(A128)
      IF (ALIN(1:1).EQ.' ') ALIN=ALIN(2:128)
      IF (ITYP.GE.30) ALIN=ALIN(19:128)
      IF (ITYP.LT.20) THEN
         READ(ALIN,20) ADATE,ATIME
   20    FORMAT(A10,1X,A8)
         READ(ADATE,25) IY,IM,ID
   25    FORMAT(I4,1X,I2,1X,I2)
      ELSE
         READ(ALIN,30) ADATE,ATIME
   30    FORMAT(A8,1X,A8)
         READ(ADATE,40) IY,IM,ID
         IF (IY.GE.98) THEN
            IY=IY+1900
         ELSE
            IY=IY+2000
         ENDIF
      ENDIF
c     print*,ADATE,ATIME
      READ(ATIME,40) JH,JM,JS
   40 FORMAT(I2,1X,I2,1X,I2)
C
      IF (ITYP.GE.30) THEN
         DO IA=1,NA
           READ(IU,*,ERR=800,END=800) DD,(DAT(IZ,IA),IZ=1,NZ)
         ENDDO
      ELSE
         DO IA=1,NA
           READ(IU,*,ERR=800,END=800) (DAT(IZ,IA),IZ=1,NZ)
         ENDDO
      ENDIF
      GOTO 900
C
  800 ERC='RDDISK: Read error !'
  900 RETURN
      END
C
      SUBROUTINE CALSVA(NZ,NA,DX,DY,AM,BX1,BY1,RAA,BET,SVA,RM,ZM
     &                 ,NDR,RRR,ZZZ,EROR,ERC)
C
C   Calculation of solid view angle.
C
C --- history
C   2001.07.06 Created from 'solid3m' by M.Yamano
C
C --- Input
C NZ         I           number of grids in zenithal direction
C NA         I           number of grids in azimuthal direction
C DX         R(NZ,NA)    azimuthal angles in degree
C                                  relative to the solar disk center
C DY         R(NZ,NA)    zenithal angles in degree
C                                  relative to the solar disk center
C AM         R(NZ,NA)    measured signals at (DX,DY)
C NDR        I           number of grids in integrating distance
C RRR        R(NDR)      integrating distance in degree
C                                  from the solar disk center
C --- Output
C BX1,BY1    R           origin of the ellipse-signals
C RAA        R           direction of major axis of the ellipse in degree
C BET        R           ratio of minor axis to major (< 1)
C SVA        R           solid view angle
C RM         R(NZ,NA)    distance data of grids in degree
C ZM         R(NZ,NA)    calculated signals at (DX,DY)
C ZZZ        R(NDR)      calculated signals at RRR
C EROR       R           RMSD of ABS(ZM-AM)
C ERC        C*64        ERROR CODE. IF ' ' THEN NORMAL.
C ----------------------------------------------------------------------
      PARAMETER (KN    =21)
      PARAMETER (KN2   =201)
      PARAMETER (KNN=KN*KN)
C
      PARAMETER (PI=3.141592653589793,RAD=PI/180.0)
C
      CHARACTER ERC*(*)
      DIMENSION DX(KN,KN),DY(KN,KN),AM(KN,KN),RM(KN,KN),ZM(KN,KN)
     &         ,RR(KNN),ZZ(KNN),IRR(KNN),RR1(KNN),ZZ1(KNN)
     &         ,RRR(KN2),ZZZ(KN2)
     &         ,ZZ0(KN2),CA(KN2),CB(KN2),CC(KN2),CD(KN2)
C
      DATA NFLT,SERR /9, 0.005/
C
      ERC=' '
      B0=0.
      BX1=0.
      BY1=0.
      DO IA=1,NA
        DO IZ=1,NZ
          B0=B0+AM(IZ,IA)
          BX1=BX1+DX(IZ,IA)*AM(IZ,IA)
          BY1=BY1+DY(IZ,IA)*AM(IZ,IA)
        ENDDO
      ENDDO
      BX1=BX1/B0
      BY1=BY1/B0
      BX2=0.
      BY2=0.
      BXY=0.
      DO IA=1,NA
        DO IZ=1,NZ
          DX1=DX(IZ,IA)-BX1
          DY1=DY(IZ,IA)-BY1
          BX2=BX2+DX1**2*AM(IZ,IA)
          BXY=BXY+DX1*DY1*AM(IZ,IA)
          BY2=BY2+DY1**2*AM(IZ,IA)
        ENDDO
      ENDDO
      BX2=BX2/B0
      BXY=BXY/B0
      BY2=BY2/B0
      B2P=2*(BX2+BY2)
      XX=2*(BX2-BY2)
      YY=4*BXY
      RA=ATAN2(YY,XX)/2
      IF (ABS(XX).GT.ABS(YY)) THEN
         B2M=XX/COS(2*RA)
      ELSE
         B2M=YY/SIN(2*RA)
      ENDIF
      BET=SQRT((1-B2M/B2P)/(1+B2M/B2P))
      RAA=RA/RAD
C
      CS=COS(RA)
      SN=SIN(RA)
      NR=0
      DO IA=1,NA
        DO IZ=1,NZ
          X=DX(IZ,IA)-BX1
          Y=DY(IZ,IA)-BY1
C
C normalized to the major-axis
          XP= CS*X+SN*Y
          YP=(-SN*X+CS*Y)/BET
          RM(IZ,IA)=SQRT(XP**2+YP**2)
          NR=NR+1
          RR(NR)=RM(IZ,IA)
          ZZ(NR)=AM(IZ,IA)
        ENDDO
      ENDDO
      CALL SORT1(NR,RR,IRR)
      DO I=1,NR
        RR1(I)=RR(IRR(I))
        ZZ1(I)=ZZ(IRR(I))
      ENDDO
      CALL SMTHF1(NR,RR1,ZZ1,NDR,RRR,ZZZ,NFLT,ERC)
      IF (ERC.NE.' ') GOTO 900
C
C Integration in the direction of major-axis
C
      DO I=1,NDR
        ZZ0(I)=ZZZ(I)/ZZZ(1)*2*PI*RRR(I)
      ENDDO
      CALL CSPL1(NDR,RRR,ZZ0,CA,CB,CC,CD)
      CALL CSPIN(NDR,RRR,CA,CB,CC,CD,SVA)
C
      S=0.
      N=0
      DO IA=1,NA
        DO IZ=1,NZ
          X=RM(IZ,IA)
          ZM(IZ,IA)=CSPLI(X,NDR,RRR,CA,CB,CC,CD)/2/PI/X
          IF (ZM(IZ,IA).GE.SERR) THEN
             N=N+1
             Y=AM(IZ,IA)/ZZZ(1)
             S=S+(Y/ZM(IZ,IA)-1.)**2
          ENDIF
        ENDDO
      ENDDO
      EROR=SQRT(S/FLOAT(N))
C
C Radian, re-scaling by the minor-axis
      SVA=BET*SVA*RAD**2
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
      PARAMETER (KDAY  =50)
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
      SUBROUTINE CSPIN(N,X,A,B,C,D,S)
C GETTING INTEGRATION with
C COEFFICIENTS OF NATURAL CUBIC SPLINE FITTING
C USE with CSPL1
C--- HISTORY
C 94. 8. 9  CREATED
C--- INPUT
C N      I     NBR OF DATA
C X    R(N)    INDEPENDENT VARIABLE DATA
C              X-increasing order
C A    R(N)    LAMBDA -> A   WHERE  Y=A+X*(B+X*(C+D*X))
C B    R(N)    D      -> B   I-TH FOR RANGE (X(I-1), X(I))
C C    R(N)    M         C   (A,B,C,D FOR I=1 ARE SAME AS THOSE FOR I=2)
C D    R(N)    M         D
C  These values come from CSPL1
C--- Output
C S     R      Integ(X(1),X(N)) Y(X) dx
C--- NOTES
C REF-1   P. F. DAVIS AND PHILIP RABINOWITZ (1984)
C         METHODS PF NUMERICAL INTEGRATION, SECOND EDITION
C         ACADEMIC PRESS, INC., PP612.
C$ENDI
C
      DIMENSION X(N),A(N),B(N),C(N),D(N)
C
      S=0
      IF(N.LE.1) return
      DO I=2,N
        X1=X(I-1)
        X2=X(I)
        S=S+A(I)*(X2-X1)+B(I)*(X2**2-X1**2)/2+C(I)*(X2**3-X1**3)/3
     &     +D(I)*(X2**4-X1**4)/4
      END DO
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
      SUBROUTINE SMTHF1(N,X,Y,NS,XS,YS,NFIL,ERR)
C Function smoothing by Gauss-filter method
C (X, Y) -smoothing-> (X1, Y1)
C X and X1 should be in increasing order
C (This is for computational efficiency)
C--- history
C 94. 8.13 Created by T. Nakajima
C--- Input
C N     I      Number of X-values
C X    R(N)    X (Should be X(1) <= X(2) <=...)
C Y    R(N)    Y
C NS    I      Number of X1-values
C XS   R(N1)   X1 (Should be X1(1) <= X1(2) <=...)
C NFIL  I      Number of X-values for using filtering
C              1 <= NFIL <= NS
C--- Output
C YS   R(N1)   Smoothed Y-values at X1
C ERR  C*64    Error code.  If ' ' then no-error.
C$ENDI
C
      character ERR*(*)
      dimension X(N),Y(N),XS(NS),YS(NS)
      ERR=' '
      NF=MIN(NFIL,N)
      J0=1
      DO I=1,NS
        S1=0
        S0=0
    7   IF(X(J0).LT.XS(I).AND.J0.LT.N) THEN
          J0=J0+1
          GOTO 7
         ELSE
          J1=J0-NF/2
          IF(J1.LE.0) J1=1
          J2=J1-1+NF
          IF(J2.GT.N) THEN
            J2=N
            J1=N-NFIL+1
          ENDIF
          SGM=0
          DO J=J1,J2
            SGM=SGM+(XS(I)-X(J))**2
          END DO
          SGM=SQRT(SGM/(J2-J1+1))
          DO J=J1,J2
            F=EXP(-((XS(I)-X(J))/SGM)**2)
            S0=S0+F
            S1=S1+F*Y(J)
          END DO
          YS(I)=S1/S0
        ENDIF
      END DO
      RETURN
      END
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
