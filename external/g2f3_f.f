      REAL*8 FUNCTION g2f3(EE_p,MZ_p,GF_p,PS_p,MA_p,gt_p,rg_p)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

      real*8 EE_p,MZ_p,GF_p,PS_p,MA_p,gt_p,rg_p
      real*8 EE,MZ,GF,PS,MA,gt,rg

      COMMON/det_vars/EE,MZ,GF,PS,MA,gt,rg
      EXTERNAL det_eq

      EE=EE_p
      MZ=MZ_p
      GF=GF_p
      PS=PS_p
      MA=MA_p
      gt=gt_p
      rg=rg_p
      x1=0.5
      x2=1.
      xacc=1E-6
      g2f3=zriddr(det_eq,x1,x2,xacc)
      RETURN
      END


      REAL*8 FUNCTION det_eq(g2)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)

      real*8 EE_p,MZ_p,GF_p,PS_p,MA_p,gt_p,rg_p
      real*8 EE,MZ,GF,PS,MA,gt,rg

      COMMON/det_vars/EE,MZ,GF,PS,MA,gt,rg

      g1=1/Sqrt(EE**(-2) - g2**(-2) - 2/gt**2)
      v=
     > Sqrt((2*Sqrt(2.)*3.14159)/GF + 
     > MA**2*((16*3.14159)/gt**2 - PS) - 
     > (4*Sqrt(2*3.14159)*Sqrt(GF**4*MA**4*(8*3.14159 - 
     > gt**2*PS)))/(GF**2*gt**2))/(2.*Sqrt(3.14159))

      r3= 
     > (8*GF**3*MA**4*Sqrt(3.14159)*PS)/
     > (8*Sqrt(2.)*GF**2*MA**2*3.14159**1.5 + 
     > 4*GF**3*MA**4*Sqrt(3.14159)*
     > PS + 4*3.14159*Sqrt(GF**4*MA**4*(8*3.14159 - gt**2*PS)) - 
     > Sqrt(2.)*GF*MA**2*PS*Sqrt(GF**4*MA**4*(8*3.14159 - gt**2*PS)))

      det_eq=
     > (MZ**2*(gt**2*(32*gt**2*MZ**2*(-MA**2 + MZ**2)*(-2*MA**2 + 
     > 2*MZ**2 + gt**2*(-1 + r3)*v**2) + 
     > g2**2*(32*gt**2*MZ**4*(-1 + r3)*v**2 + 2*gt**4*MZ**2*(4 + 
     > r3*(-8 + 3*r3))*v**4 - gt**6*(-1 + r3)*r3**2*v**6 + 
     > 16*MA**4*(4*MZ**2 - gt**2*v**2) + 2*MA**2*(-32*MZ**4 + 
     > 8*gt**2*MZ**2*(4 - 3*r3)*v**2 + gt**4*(-4 + r3*(4 + 
     > r3))*v**4))) + g1**2*(2*g2**2*(2*MA**2 - gt**2*(-1 + r3)*v**2)*
     > (16*MA**2*MZ**2 - 8*gt**2*(MA**2 + MZ**2*(-1 + r3))*v**2 + 
     > gt**4*r3**2*v**4) + gt**2*(32*gt**2*MZ**4*(-1 + r3)*v**2 + 
     > 2*gt**4*MZ**2*(4 + r3*(-8 + 3*r3))*v**4 - gt**6*(-1 + 
     > r3)*r3**2*v**6 + 16*MA**4*(4*MZ**2 - gt**2*v**2) + 
     > 2*MA**2*(-32*MZ**4 + 8*gt**2*MZ**2*(4 - 3*r3)*v**2 + 
     > gt**4*(-4 + r3*(4 + r3))*v**4)))))/(64.*gt**4*MA**8)
       RETURN
       END


      REAL*8 FUNCTION zriddr(func,x1,x2,xacc)
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
      INTEGER MAXIT
CAB      REAL zriddr,x1,x2,xacc,func,UNUSED
      PARAMETER (MAXIT=60,UNUSED=-1.11E30)
      EXTERNAL func
CU    USES func
      INTEGER j
CAB      REAL fh,fl,fm,fnew,s,xh,xl,xm,xnew
      fl=func(x1)
      fh=func(x2)
      if((fl.gt.0..and.fh.lt.0.).or.(fl.lt.0..and.fh.gt.0.))then
        xl=x1
        xh=x2
        zriddr=UNUSED
        do 11 j=1,MAXIT
          xm=0.5*(xl+xh)
          fm=func(xm)
          s=sqrt(fm**2-fl*fh)
          if(s.eq.0.)return
          xnew=xm+(xm-xl)*(sign(1d0,fl-fh)*fm/s)
          if (abs(xnew-zriddr).le.xacc) return
          zriddr=xnew
          fnew=func(zriddr)
          if (fnew.eq.0.) return
          if(sign(fm,fnew).ne.fm) then
            xl=xm
            fl=fm
            xh=zriddr
            fh=fnew
          else if(sign(fl,fnew).ne.fl) then
            xh=zriddr
            fh=fnew
          else if(sign(fh,fnew).ne.fh) then
            xl=zriddr
            fl=fnew
          else
cAB            pause 'never get here in zriddr'
          endif
          if(abs(xh-xl).le.xacc) return
11      continue
cAB        pause 'zriddr exceed maximum iterations'
      else if (fl.eq.0.) then
        zriddr=x1
      else if (fh.eq.0.) then
        zriddr=x2
      else
cAB        pause 'root must be bracketed in zriddr'
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<0(9p#3.
