C
C Copyright (C) 2007-2009, CompHEP Collaboration
C Author: Alexander Sherstnev 
C ------------------------------------------------------
C
      subroutine makelhapdfsilent
      CHARACTER*20 LHAPARM(20)
      REAL*8 LHAVALUE(20)
      INTEGER lhasilent
      COMMON/LHACONTROL/LHAPARM,LHAVALUE
      common/lhasilent/lhasilent

      LHAPARM(19)='SILENT'
      call setlhaparm('SILENT')
      lhasilent=1
      end

      function lhapdfqcdlam (set,mem)
      real*8 lhapdfqcdlam
      integer set,mem
#ifdef LHAPDF
      call GetLam5(mem,lhapdfqcdlam)
#else
c     dummy value...
      lhapdfqcdlam=0.0*real(set)*real(mem)
#endif
      return
      end

      function lhapdfqcdorder ()
      integer lhapdfqcdorder
#ifdef LHAPDF
      call GetOrderAs(lhapdfqcdorder)
#else
c     dummy value...
      lhapdfqcdorder=0
#endif
      return
      end

      function lhapdfalphas (Q)
      real*8 lhapdfalphas
      real Q
#ifdef LHAPDF
      lhapdfalphas = alphasPDF(Q)
#else
c     dummy value...
      lhapdfalphas=0.0*Q
#endif
      return
      end

      function locevolvepdf (X,Q,lPDF)
      integer I, locevolvepdf
      real*8 X,Q,lPDF(0:12),lF(-6:6)

      DO I=0,12 
        lPDF(I)=0.0
      ENDDO
#ifdef LHAPDF
      call evolvePDF(X,Q,lF)
      DO I=0,12 
        lPDF(I)=lF(I-6)
      ENDDO
#endif
      locevolvepdf=1
      return
      end
