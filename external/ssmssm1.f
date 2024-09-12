CDECK  ID>, SSMSSM.
      SUBROUTINE SSMSSM1(XMG,XMU,XMHA,XTANB,XMQ1,XMDR,XMUR,
     $XML1,XMER,XMQ2,XMSR,XMCR,XML2,XMMR,XMQ3,XMBR,XMTR,
     $XML3,XMLR,XAT,XAB,XAL,XM1,XM2,XMT,
     $XMH, XMHH, XMHC, XALPHA,  IALLOW,IMODEL)
C-----------------------------------------------------------------------
C
C     Calculate MSSM masses and decays using parameters:
C       XM1    = U(1) mass
C              > 1e19: use scaling from XMG
C       XM2    = SU(2) mass
C              > 1e19: use scaling from XMG
C       XMG    = gluino mass
C       XMQ1,...  = 1st gen. su(2) soft squark mass,...
C       XMTL   = m(stop-left)
C       XMTR   = m(stop-right)
C       XMBR   = m(sbot-right)
C       XML1   = left selectron mass
C       XMER   = right selectron mass
C       XMN1   = 1st ge. sneutrino mass
C       XTANB  = v/v' = ratio of vev's
C       XMU    = -2*m_1 = SUSY Higgs mass
C       XMHA   = m(pseudo-scalar-Higgs)
C       XMT    = m(top)
C       XAT    = stop trilinear coupling
C       XAB    = sbottom trilinear coupling
C       XAL    = stau trilinear coupling
C       IALLOW = 0 for valid point, 1 otherwise
C       IMODEL = 1 for SUGRA or MSSM, 2 for GMSB
C
C     Program outline:
C     SSMSSM:  Initialize standard model parameters in /SSSM/ and 
C              SUSY parameters in /SSPAR/.
C     SSMASS:  Calculate dependent SUSY masses and mixings.
C     SSTPBF:  Calculate top decays; save in /SSMODE/.
C     SSSTBF:  Calculate stop decays; save in /SSMODE/.
C     SSGLBF:  Calcualte gluino decays; save in /SSMODE/.
C     SSQKBF:  Calculate squark decays; save in /SSMODE/.
C     SSWZBF:  Calculate gaugino decays; save in /SSMODE/.
C     SSHIBF:  Calculate Higgs decays; save in /SSMODE/.
C
C     Notes: 
C  1) All particle ID codes are defined with symbolic names in 
C     /SSTYPE/, making it easy to change them.
C
C  2) /SSMODE/ contains the parent, the daughters, the width, and
C     the branching ratio for each mode. Decay modes for a given parent
C     need not be adjacent, so they must be sorted at the end.
C
C  3) Some of Baer's original routines used single precision and others
C     double precision. To accomodate this, the variable names used in
C     /SSSM/ and /SSPAR/ have all been changed to longer, more 
C     descriptive ones.
C
C  4) All routines have been strongly typed.
C
C     Source: H. Baer, et al.
C     Modified: F. Paige, Aug. 1992
C-----------------------------------------------------------------------
      IMPLICIT NONE
      COMMON/SSLUN/LOUT
      INTEGER LOUT
      SAVE /SSLUN/
C          MXSS         =  maximum number of modes
C          NSSMOD       = number of modes
C          ISSMOD       = initial particle
C          JSSMOD       = final particles
C          GSSMOD       = width
C          BSSMOD       = branching ratio
C          MSSMOD       = decay matrix element pointer
C          LSSMOD       = logical flag used internally by SSME3
      INTEGER MXSS
      PARAMETER (MXSS=1000)
      COMMON/SSMODE/NSSMOD,ISSMOD(MXSS),JSSMOD(5,MXSS),GSSMOD(MXSS)
     $,BSSMOD(MXSS),MSSMOD(MXSS),LSSMOD
      INTEGER NSSMOD,ISSMOD,JSSMOD,MSSMOD
      REAL GSSMOD,BSSMOD
      LOGICAL LSSMOD
      SAVE /SSMODE/
C          Standard model parameters
C          AMUP,...,AMTP        = quark masses
C          AME,AMMU,AMTAU       = lepton masses
C          AMW,AMZ              = W,Z masses
C          GAMW,GAMZ            = W,Z widths
C          ALFAEM,SN2THW,ALFA3  = SM couplings
C          ALQCD4               = 4 flavor lambda
      COMMON/SSSM/AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      REAL AMUP,AMDN,AMST,AMCH,AMBT,AMTP,AME,AMMU,AMTAU
     $,AMW,AMZ,GAMW,GAMZ,ALFAEM,SN2THW,ALFA2,ALFA3,ALQCD4
      SAVE /SSSM/
C          SUSY parameters
C          AMGLSS               = gluino mass
C          AMULSS               = up-left squark mass
C          AMELSS               = left-selectron mass
C          AMERSS               = right-slepton mass
C          AMNiSS               = sneutrino mass for generation i
C          TWOM1                = Higgsino mass = - mu
C          RV2V1                = ratio v2/v1 of vev's
C          AMTLSS,AMTRSS        = left,right stop masses
C          AMT1SS,AMT2SS        = light,heavy stop masses
C          AMBLSS,AMBRSS        = left,right sbottom masses
C          AMB1SS,AMB2SS        = light,heavy sbottom masses
C          AMLLSS,AMLRSS        = left,right stau masses
C          AML1SS,AML2SS        = light,heavy stau masses
C          AMZiSS               = signed mass of Zi
C          ZMIXSS               = Zi mixing matrix
C          AMWiSS               = signed Wi mass
C          GAMMAL,GAMMAR        = Wi left, right mixing angles
C          AMHL,AMHH,AMHA       = neutral Higgs h0, H0, A0 masses
C          AMHC                 = charged Higgs H+ mass
C          ALFAH                = Higgs mixing angle
C          AAT                  = stop trilinear term
C          THETAT               = stop mixing angle
C          AAB                  = sbottom trilinear term
C          THETAB               = sbottom mixing angle
C          AAL                  = stau trilinear term
C          THETAL               = stau mixing angle
C          AMGVSS               = gravitino mass
      COMMON/SSPAR/AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS(4,4)
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS
      REAL AMGLSS,AMULSS,AMURSS,AMDLSS,AMDRSS,AMSLSS
     $,AMSRSS,AMCLSS,AMCRSS,AMBLSS,AMBRSS,AMB1SS,AMB2SS
     $,AMTLSS,AMTRSS,AMT1SS,AMT2SS,AMELSS,AMERSS,AMMLSS,AMMRSS
     $,AMLLSS,AMLRSS,AML1SS,AML2SS,AMN1SS,AMN2SS,AMN3SS
     $,TWOM1,RV2V1,AMZ1SS,AMZ2SS,AMZ3SS,AMZ4SS,ZMIXSS
     $,AMW1SS,AMW2SS
     $,GAMMAL,GAMMAR,AMHL,AMHH,AMHA,AMHC,ALFAH,AAT,THETAT
     $,AAB,THETAB,AAL,THETAL,AMGVSS
      REAL AMZISS(4)
      EQUIVALENCE (AMZISS(1),AMZ1SS)
      SAVE /SSPAR/
C
C          Data for SUSY 3-body matrix elements. There is a double 
C          pointer structure, first to modes, and then to poles that
C          make up the matrix element for that mode:
C          MELEM=-I in /DKYTAB/ points to the mode information:
C            J1SS3(I) = start of pole list for this mode
C            J2SS3(I) = end of pole list for this mode
C            WTSS3(I) = maximum weight for this mode
C          J1SS3<J<J2SS3 points to the corresponding poles:
C            KSS3(J)    = pole type
C            AMSS3(J)   = pole mass
C            ZISS3(2,J) = initial couplings
C            ZFSS3(2,J) = final couplings
C          For gaugino -> gaugino f fbar, the pole types are
C            KSS3=1: spin-1 pole in f-fbar channel
C            KSS3=2: spin-0 pole in gaugino-f channel
C            KSS3=3: spin-0 pole in gaugino-fbar channel
C            KSS3=4: spin-0 pole in f-fbar channel
C          The two couplings are the coefficients of 1,gamma_5 or of
C          gamma_mu,gamma_mu*gamma_5. 
C
      INTEGER MXMSS3,MXPSS3
      PARAMETER (MXMSS3=1000)
      PARAMETER (MXPSS3=2000)
      COMMON/DKYSS3/NMSS3,NPSS3,
     $J1SS3(MXMSS3),J2SS3(MXMSS3),WTSS3(MXMSS3),
     $KSS3(MXPSS3),AMSS3(MXPSS3),ZISS3(2,MXPSS3),ZFSS3(2,MXPSS3)
      INTEGER NMSS3,NPSS3,KSS3,J1SS3,J2SS3
      REAL WTSS3,AMSS3
      COMPLEX ZISS3,ZFSS3
C
      REAL XR21,PI,SR2
      REAL XMG,XMU,XMHA,XTANB,XMQ1,XMDR,XMUR,XML1,XMER,XMQ2,XMSR,
     $XMCR,XML2,XMMR,XMQ3,XMBR,XMTR,XML3,XMLR,XAT,XAB,XAL,XM1,XM2,
     $XMT,MU1,MU2,BETA,COS2B, XMH, XMHH, XMHC, XALPHA
      INTEGER IALLOW,MHLNEG,MHCNEG,IMODEL,I
      COMMON/ISAWDTH/ IWIDTH
      REAL IWIDTH(32)
      INTEGER ISUPL,ISDNL,ISSTL,ISCHL,ISBT1,ISTP1
      INTEGER ISNEL,ISEL,ISNML,ISMUL,ISNTL,ISTAU1
      INTEGER ISUPR,ISDNR,ISSTR,ISCHR,ISBT2,ISTP2
      INTEGER ISNER,ISER,ISNMR,ISMUR,ISNTR,ISTAU2
      INTEGER ISZ1,ISZ2,ISZ3,ISZ4,ISW1,ISW2,ISGL
      INTEGER ISHL,ISHH,ISHA,ISHC
      INTEGER ISGRAV
      PARAMETER (ISUPL=21,ISDNL=22,ISSTL=23,ISCHL=24,ISBT1=25,ISTP1=26)
      PARAMETER (ISNEL=31,ISEL=32,ISNML=33,ISMUL=34,ISNTL=35,ISTAU1=36)
      PARAMETER (ISUPR=41,ISDNR=42,ISSTR=43,ISCHR=44,ISBT2=45,ISTP2=46)
      PARAMETER (ISNER=51,ISER=52,ISNMR=53,ISMUR=54,ISNTR=55,ISTAU2=56)
      PARAMETER (ISGL=29)
      PARAMETER (ISZ1=30,ISZ2=40,ISZ3=50,ISZ4=60,ISW1=39,ISW2=49)
      PARAMETER (ISHL=82,ISHH=83,ISHA=84,ISHC=86)
      PARAMETER (ISGRAV=91)
C
      NSSMOD=0
C
C          Standard model and SUSY parameters
C
      IALLOW=0
      XR21=1./XTANB
      PI=4.*ATAN(1.)
      SR2=SQRT(2.)
      AMDN=0.0099
      AMUP=0.0056
      AMST=0.199
      AMCH=1.35
      AMBT=5.0
      AMTP=XMT
      AME=0.511E-3
      AMMU=0.105
      AMTAU=1.777
      AMW=80.0
      AMZ=91.17
      GAMW=2.12
      GAMZ=2.487
      ALFAEM=1./128.
      SN2THW=0.232
      ALFA2=ALFAEM/SN2THW
      BETA=ATAN(XTANB)
      COS2B=COS(2*BETA)
C
C          SU(2) and U(1) gaugino masses are reset in SSMASS if
C          they are > 1e19.
C
      MU2=XM2
      MU1=XM1
C          Set 2nd gen soft terms equal to 1st gen. soft terms 
c          unless previously set by user.
      IF (XMQ2.GE.1.E19) THEN
        XMQ2=XMQ1
        XMSR=XMDR
        XMCR=XMUR
        XML2=XML1
        XMMR=XMER
      END IF
C
C          The results can be quite sensitive to the choice of the
C          4-flavor QCD scale ALQCD4 and the expression for the QCD
C          coupling ALFA3. Select among the following lines:
C
      ALQCD4=0.177
      ALFA3=0.12
C
C          Calculate simple masses; other masses via SSMASS
      AMGLSS=XMG
      AMULSS=SQRT(XMQ1**2+AMUP**2+(.5-2.*SN2THW/3.)*AMZ**2*COS2B)
      AMURSS=SQRT(XMUR**2+AMUP**2+2./3.*SN2THW*AMZ**2*COS2B)
      AMDLSS=SQRT(XMQ1**2+AMDN**2+(-.5+SN2THW/3.)*AMZ**2*COS2B)
      AMDRSS=SQRT(XMDR**2+AMDN**2-1./3.*SN2THW*AMZ**2*COS2B)
      AMCLSS=SQRT(XMQ2**2+AMCH**2+(.5-2.*SN2THW/3.)*AMZ**2*COS2B)
      AMCRSS=SQRT(XMCR**2+AMCH**2+2./3.*SN2THW*AMZ**2*COS2B)
      AMSLSS=SQRT(XMQ2**2+AMST**2+(-.5+SN2THW/3.)*AMZ**2*COS2B)
      AMSRSS=SQRT(XMSR**2+AMST**2-1./3.*SN2THW*AMZ**2*COS2B)
      AMELSS=SQRT(XML1**2+AME**2-(.5-SN2THW)*AMZ**2*COS2B)
      AMERSS=SQRT(XMER**2+AME**2-SN2THW*AMZ**2*COS2B)
      AMMLSS=SQRT(XML2**2+AMMU**2-(.5-SN2THW)*AMZ**2*COS2B)
      AMMRSS=SQRT(XMMR**2+AMMU**2-SN2THW*AMZ**2*COS2B)
      AMN1SS=SQRT(XML1**2+.5*AMZ**2*COS2B)
      AMN2SS=SQRT(XML2**2+.5*AMZ**2*COS2B)
      AMN3SS=SQRT(XML3**2+.5*AMZ**2*COS2B)
      AMTLSS=XMQ3
      AMTRSS=XMTR
      AMBLSS=XMQ3
      AMBRSS=XMBR
      AMLLSS=XML3
      AMLRSS=XMLR
      AMHA=XMHA
      AAT=XAT
      AAB=XAB
      AAL=XAL
      TWOM1=-XMU
      RV2V1=XR21
      
      amgvss=1.E20
C
C          Calculate mass eigenstates and check Z1SS = LSP
C
      CALL SSMASS(MU1,MU2,IALLOW,1,MHLNEG,MHCNEG,IMODEL)
      IF (MHLNEG.EQ.1.OR.MHCNEG.EQ.1) IALLOW=10
C     IF(IALLOW.NE.0) RETURN
C
C          Initialize counters for matrix elements
C          Calculate decay widths and branching rations
C
      write(*,*) 'ISASUSY values MHL=',amhl, ' MHH=', amhh

      write(*,*) 'value of alpha: ', ALFAH, ' changed to ',
     + -xalpha
     
     
      ALFAH=-XALPHA
      AMHL=XMH
      AMHH=XMHH
      AMHC=XMHC
      NMSS3=0
      NPSS3=0
      CALL SSTPBF
      CALL SSGLBF
      CALL SSQKBF
      CALL SSSTBF
      CALL SSLPBF
      CALL SSWZBF
      CALL SSHIBF
C

      DO 120 I=1,32
 120  IWIDTH(I)=0

C gaugino
      
      DO 121 I=1,NSSMOD
 121  IF(ISSMOD(I).EQ.ISGL) IWIDTH(1)=IWIDTH(1)+GSSMOD(I)
      
      DO 122 I=1,NSSMOD
 122  IF(ISSMOD(I).EQ.ISZ1) IWIDTH(2)=IWIDTH(2)+GSSMOD(I)
      
      DO 123 I=1,NSSMOD
 123  IF(ISSMOD(I).EQ.ISZ2) IWIDTH(3)=IWIDTH(3)+GSSMOD(I)
      
      DO 124 I=1,NSSMOD
 124  IF(ISSMOD(I).EQ.ISZ3) IWIDTH(4)=IWIDTH(4)+GSSMOD(I)
      
      DO 125 I=1,NSSMOD
 125  IF(ISSMOD(I).EQ.ISZ4) IWIDTH(5)=IWIDTH(5)+GSSMOD(I)
            
      DO 126 I=1,NSSMOD
 126  IF(ISSMOD(I).EQ.ISW1) IWIDTH(6)=IWIDTH(6)+GSSMOD(I)
      
      DO 127 I=1,NSSMOD
 127  IF(ISSMOD(I).EQ.ISW2) IWIDTH(7)=IWIDTH(7)+GSSMOD(I)
      
C higgses
      
      DO 128 I=1,NSSMOD
 128  IF(ISSMOD(I).EQ.ISHL) IWIDTH(8)=IWIDTH(3)+GSSMOD(I)
      
      DO 129 I=1,NSSMOD
 129  IF(ISSMOD(I).EQ.ISHH) IWIDTH(9)=IWIDTH(4)+GSSMOD(I)
      
      DO 130 I=1,NSSMOD
 130  IF(ISSMOD(I).EQ.ISHA) IWIDTH(10)=IWIDTH(10)+GSSMOD(I)
      
      DO 131 I=1,NSSMOD
 131  IF(ISSMOD(I).EQ.ISHC) IWIDTH(11)=IWIDTH(11)+GSSMOD(I)

C  sleptons

      DO 132 I=1,NSSMOD
 132  IF(ISSMOD(I).EQ.ISEL) IWIDTH(12)=IWIDTH(12)+GSSMOD(I)
      
      DO 133 I=1,NSSMOD
 133  IF(ISSMOD(I).EQ.ISER) IWIDTH(13)=IWIDTH(13)+GSSMOD(I)

      DO 134 I=1,NSSMOD
 134  IF(ISSMOD(I).EQ.ISNEL) IWIDTH(14)=IWIDTH(14)+GSSMOD(I)
      
      DO 135 I=1,NSSMOD
 135  IF(ISSMOD(I).EQ.ISMUL) IWIDTH(15)=IWIDTH(15)+GSSMOD(I)

      DO 136 I=1,NSSMOD
 136  IF(ISSMOD(I).EQ.ISMUR) IWIDTH(16)=IWIDTH(16)+GSSMOD(I)
      
      DO 137 I=1,NSSMOD
 137  IF(ISSMOD(I).EQ.ISNML) IWIDTH(17)=IWIDTH(17)+GSSMOD(I)

      DO 138 I=1,NSSMOD
 138  IF(ISSMOD(I).EQ.ISTAU1) IWIDTH(18)=IWIDTH(18)+GSSMOD(I)
      
      DO 139 I=1,NSSMOD
 139  IF(ISSMOD(I).EQ.ISTAU2) IWIDTH(19)=IWIDTH(19)+GSSMOD(I)

      DO 140 I=1,NSSMOD
 140  IF(ISSMOD(I).EQ.ISNTL) IWIDTH(20)=IWIDTH(20)+GSSMOD(I)
      
C  squarks
      
      DO 141 I=1,NSSMOD
 141  IF(ISSMOD(I).EQ.ISUPL) IWIDTH(21)=IWIDTH(21)+GSSMOD(I)

      DO 142 I=1,NSSMOD
 142  IF(ISSMOD(I).EQ.ISUPR) IWIDTH(22)=IWIDTH(22)+GSSMOD(I)

      DO 143 I=1,NSSMOD
 143  IF(ISSMOD(I).EQ.ISDNL) IWIDTH(23)=IWIDTH(23)+GSSMOD(I)

      DO 144 I=1,NSSMOD
 144  IF(ISSMOD(I).EQ.ISDNR) IWIDTH(24)=IWIDTH(24)+GSSMOD(I)

      DO 145 I=1,NSSMOD
 145  IF(ISSMOD(I).EQ.ISSTL) IWIDTH(25)=IWIDTH(25)+GSSMOD(I)

      DO 146 I=1,NSSMOD
 146  IF(ISSMOD(I).EQ.ISSTR) IWIDTH(26)=IWIDTH(26)+GSSMOD(I)

      DO 147 I=1,NSSMOD
 147  IF(ISSMOD(I).EQ.ISCHL) IWIDTH(27)=IWIDTH(27)+GSSMOD(I)

      DO 148 I=1,NSSMOD
 148  IF(ISSMOD(I).EQ.ISCHR) IWIDTH(28)=IWIDTH(28)+GSSMOD(I)

      DO 149 I=1,NSSMOD
 149  IF(ISSMOD(I).EQ.ISBT1) IWIDTH(29)=IWIDTH(29)+GSSMOD(I)

      DO 150 I=1,NSSMOD
 150  IF(ISSMOD(I).EQ.ISBT2) IWIDTH(30)=IWIDTH(30)+GSSMOD(I)
      
      DO 151 I=1,NSSMOD
 151  IF(ISSMOD(I).EQ.ISTP1) IWIDTH(31)=IWIDTH(31)+GSSMOD(I)
      
      DO 152 I=1,NSSMOD
 152  IF(ISSMOD(I).EQ.ISTP2) IWIDTH(32)=IWIDTH(32)+GSSMOD(I)
      
      write(*,*) ' ----------------------- ISASUSY masses and widthes',
     +           ' -----------------------'      
      write(*,800) 'Z1', AMZ1SS, IWIDTH(2),
     +             'Z2', AMZ2SS, IWIDTH(3)
      write(*,800) 'Z3', AMZ3SS, IWIDTH(4),
     +             'Z4', AmZ4SS, IWIDTH(5)
      write(*,800) 'W1', AMW1SS, IWIDTH(6),
     +             'W2', AMW2SS, IWIDTH(7)
      write(*,800) 'HL', AMHL,   IWIDTH(8),
     +             'HH', AMHH,   IWIDTH(9)
      write(*,800) 'HA', AMHA,   IWIDTH(10),
     +             'HC', AMHC,   IWIDTH(11)
      write(*,800) 'EL', AMELSS, IWIDTH(12),
     +             'ER', AMERSS, IWIDTH(13)
      write(*,800) 'N1', AMN1SS, IWIDTH(14),
     +             'ML', AMMLSS, IWIDTH(15)
      write(*,800) 'MR', AMMRSS, IWIDTH(16),
     +             'N2', AMN2SS, IWIDTH(17)
      write(*,800) 'T1', AML1SS, IWIDTH(18),
     +             'T2', AML2SS, IWIDTH(19)
      write(*,800) 'N3', AMN3SS, IWIDTH(20),
     +             'GL', AMGLSS, IWIDTH(1)
      write(*,800) 'UL', AMULSS, IWIDTH(21),
     +             'UR', AMURSS, IWIDTH(22)
      write(*,800) 'DL', AMDLSS, IWIDTH(23),
     +             'DR', AMDRSS, IWIDTH(24)
      write(*,800) 'SL', AMSLSS, IWIDTH(25),
     +             'SR', AMSRSS, IWIDTH(26)
      write(*,800) 'CL', AMCLSS, IWIDTH(27),
     +             'CR', AMCRSS, IWIDTH(28)
      write(*,800) 'B1', AMB1SS, IWIDTH(29),
     +             'B2', AMB2SS, IWIDTH(30)
      write(*,800) 'T1', AMT1SS, IWIDTH(31),
     +             'T2', AMT2SS, IWIDTH(32)
     
 800  FORMAT(2X,A4,3X,F7.1,3X,E9.3,6X,A4,3X,F7.1,3X,E9.3)
 
      
C
c      WRITE(*,*) NSSMOD
c      do 223 I=1,NSSMOD
c      write(*,*) label(ISSMOD(I)),' --> ',label(JSSMOD(1,I)),' ',
c     + label(JSSMOD(2,I)),' ',label(JSSMOD(3,I)),
c     + '   ',GSSMOD(I),BSSMOD(I)
c 223   CONTINUE


      RETURN
      END
      
      REAL FUNCTION SSVALUE(ind)
      INTEGER ind
      COMMON/ISAWDTH/ IWIDTH
      REAL IWIDTH(32)
      
      SSVALUE=IWIDTH(ind)
      RETURN
      END
