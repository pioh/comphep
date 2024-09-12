*
* Copyright (C) 2003-2009, CompHEP Collaboration
* Copyright (C) 2003, Slava Bunichev
*
* Form procedures for CompHEP ver 0.92

*----------- PROCEDURES ------------------


************ Vector calculation functions ***************

#procedure SubstOptimal()
  id '$MaxDpMom'= '$MaxDpSub';  
#endprocedure

#procedure DefineF()
  .sort
  id F(s?,t?)=d_(s,t);
  id F(z?,s?)=z(s);
  .sort
#endprocedure

#procedure RepG(x,ind)
  .sort
  skip; nskip Fl'x';
  id g_(ln,j?)=g_('ind',j);
  id g_(ln,5_)=g_('ind',5_);
  id gi_(ln)=gi_('ind');
  .sort
#endprocedure

#procedure RepIND(x,ind,mom)
  .sort
  skip; nskip Vrt'x';
  id z?('ind')=z.moment;
  id d_('ind',s?)=moment(s);
  id g_(j?,'ind')=g_(j,moment);
  id moment = 'mom';
  .sort
#endprocedure

#procedure G5ravno0(x)
  .sort
  skip; nskip Vrt'x';
  id g_(ln,5_)=0;
  .sort
#endprocedure

#procedure Iravno0(x,y)
  .sort
  skip; nskip Vrt1r,Vrt2r;
  id g_(j?,5_)=0;
  .sort
  L Vrt'x'=Vrt'x'-Vrt1r;
  L Vrt'y'=Vrt'y'-Vrt2r;
  .sort
    #endprocedure



***************** Denominator Functions ****************************

#procedure DefineDenominator()
  .sort
  L NumDenominator=denominator;
  .sort
  skip; nskip NumDenominator;
  id propDen(mom?,mas?,0)=1;
  id propDen(mom?,mas?,width?)=dum_(mom.mom-mas^2)-i_*mas*width;
  id i_=0;
  .sort
  id propDen(mom?,mas?,0)=dum_(mom.mom-mas^2);
  id propDen(mom?,mas?,width?)=dum_(dum_(mom.mom-mas^2)^2+(mas*width)^2);
  .sort
  skip; nskip NumDenominator,denominator;

  argument;
    #call Substitution()
    #call SubstOptimal()
    #call DotProduct()
    argument;
      #call Substitution()
      #call SubstOptimal()
      #call DotProduct()
    endargument;
  endargument;
  .sort
#endprocedure



********************** Gorner Functions ****************************

#procedure OtherVar()
  .sort 
  Ab MZ,MW,Mt,Mb,Mc,Mu,Md,Ms,Mtop,Me,MH,Mm,SW,EE,GG,Cw,CW,
  Vus,Vcd,Vtb,Vcb,Vud,Vub,Vts,Vtd,Vcs,wtop,wZ,wW,wH,
  s12,s23,s13,c12,c23,c13,Sqrt2;
  .sort 
  collect dum_;
  .sort
#endprocedure


#procedure GornerInfo(DPnum)
  .sort
*
  #$steps=0;
  if(match(dp(?x)) > $steps) $steps=count_(dp,1);
  .sort
  #$steps=$steps-1;    
*
  #do i=0,27
    id dp('i')=d'i';
  #enddo
*
  #$maxcount=0;
  #do i=0,'DPnum'
     if(match(d'i') > $maxcount) $maxcount=count_(d'i',1);
     .sort    
  #enddo
*    
  #do k=1,'$maxcount'
     #do i=0,'DPnum'  
        id once d'i'=cc'i'c'k';
     #enddo
  #enddo
*
  .sort
#endprocedure



#procedure Min(x,y)
  .sort
  #$min=10000000000000000;
  #do s=1,'y'
     #do k=0,'x'
       #$j=0;
       if(match(cc'k'c's')>0) $j=$j+1;
       .sort
       #if(('$j'<='$min')&&('$j'>0))
         #$ind1='k';
         #$ind2='s';
         #$min=$j;
       #endif  
     #enddo
  #enddo
  .sort
#endprocedure



#procedure Gorner(DPnum)
  .sort
  #call GornerInfo('DPnum')
*
  #do t=0,'$steps'    
      #do l=0,'$maxcount'
      #do i=0,'DPnum'
        #call Min('DPnum','$maxcount')
        if(match(ANT(x?))==0) id once cc'$ind1'c'$ind2'=ANT('$ind1'); 
        .sort
        id once cc'$ind1'c'$ind2'=H(cc'$ind1'c'$ind2'); 
        .sort
      #enddo
      #enddo 
*       
      Ab ANT,dum_;
      .sort
      collect dum_;
      .sort
*   
      argument;
       id ANT(x?)=dp(x);
      endargument;
      .sort
*
      id H(x?)=x;
      .sort
  #enddo
  .sort
#endprocedure





**************** Output Functions **************************

#procedure DotProduct()
  #$kk=0;
  #do i=2,8
    #do j=1,'i'-1
      id p'j'.p'i'=dp('$kk');
      #$kk=$kk+1;
    #enddo
  #enddo
#endprocedure


#procedure OutputTest(xx,xxx)
  .sort
  #write <./results/t'xx'_'xxx'.prc> " #\procedure FromForm'xx'p'xxx'() "
  #write <./results/t'xx'_'xxx'.prc> " .sort"
  #write <./results/t'xx'_'xxx'.prc>  " L totFact= %e",totFactor
  #write <./results/t'xx'_'xxx'.prc>  " L numerat= %e",numerator
  #write <./results/t'xx'_'xxx'.prc> " .sort"
  #write <./results/t'xx'_'xxx'.prc> " id FromForm()= totFact*numerat;"
  #write <./results/t'xx'_'xxx'.prc> " .sort"
  #write <./results/t'xx'_'xxx'.prc> " #\endprocedure"
  .sort
#endprocedure


#procedure OutputC(x)
  .sort
  format c;
  #write <./results/f'x'.c> " #include<math.h>"
  #write <./results/f'x'.c> " #define real double "
  #write <./results/f'x'.c> " extern real va[35]; "
  #write <./results/f'x'.c> " #include\"out_ext.h\" "
  #write <./results/f'x'.c> " #include\"out_int.h\" "
  #write <./results/f'x'.c> " #include\"var_def.h\" "
  #write <./results/f'x'.c> " FNN F'x'; "
  #write <./results/f'x'.c> " real F'x'(void)"
  #write <./results/f'x'.c> " {"
  #write <./results/f'x'.c> " real FACTOR,RNUM,NUMDENOM,DENOM,result;"
  #write <./results/f'x'.c>  " FACTOR= %e",totFactor(FACTOR)
  #write <./results/f'x'.c>  " RNUM= %e",numerator(RNUM)
  #write <./results/f'x'.c>  " NUMDENOM= %e",NumDenominator(NUMDENOM)
  #write <./results/f'x'.c>  " DENOM= %e",denominator(DENOM)
  #write <./results/f'x'.c> " result=FACTOR*RNUM*(NUMDENOM/DENOM);"
  #write <./results/f'x'.c> " if(result>Fmax) Fmax=result; else if(result<-Fmax) Fmax=-result;"
  #write <./results/f'x'.c> "if(color_weights)"
  #write <./results/f'x'.c> "{"
  #write <./results/f'x'.c> "color_weights[0] += result*(1)/(1);"
  #write <./results/f'x'.c> "}"
  #write <./results/f'x'.c> " return result;"
  #write <./results/f'x'.c> " }"
  .sort
#endprocedure

