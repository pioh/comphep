%    ==============================
%    *  CompHEP version 4.5.2rc10    *
%    ==============================
inParticles:={"u","b"}$
outParticles:={"d","t"}$
%
parameters:={ EE=>3.134500E-01,SW=>4.807600E-01,s12=>2.229000E-01,s13=>
 0.000000E+00,MZ=>9.118760E+01,Mtop=>1.725000E+02,Mb=>4.850000E+00,FL1=>
 1.000000E+00,FR1=>0.000000E+00,FL2=>0.000000E+00,FR2=>0.000000E+00}$
%
substitutions:={ FFR2=>-FR2/2/MZ/CW,FFL2=>-FL2/2/MZ/CW,MW=>MZ*CW,Vud=>c12*
 c13,c13=>sqrt(1-s13^2),c12=>sqrt(1-s12^2),CW=>sqrt(1-SW^2)}$


vector p1,p2,p3,p4,p5,p6$

let p4 = +p1+p2-p3;
mass p1  = 0; Mshell p1;
mass p2  = Mb; Mshell p2;
mass p3  = 0; Mshell p3;
let p2.p3 = -1*(Mtop^2-0^2-Mb^2-0^2-2*p1.p2+2*p1.p3)/2;

vector !=p_,!=q_$
operator propDen$
for all p_,q_,m,w let propDen(0*p_+q_,m,w)=propDen(q_,m,w)$
for all p_,m,w such that ordp(p_,-p_) let propDen(p_,m,w)=propDen(-p_,m,w);$

initSum();

DiagrNumber:="1_1"$

%                     u     d    !  d     u                          
%                   ==>==@==>====!==>==@==>==                        
%                     P1 |  P3   !  P3 |  P1                         
%                      W+|       !   W+|                             
%                     b  |  t    !  t  |  b                          
%                   ==>==@==>====!==>==@==>==                        
%                     P2    P4   !  P4    P2                         
totFactor:=(Vud^2*EE^4)/(SW^4)$
numerator:=8*p3.p4*p2.p3*p1.p3*FFL2^2+4*p3.p4*p1.p3*FFL2*FL1*Mb+p3.p4*p1.p2*
 FL1^2+p1.p4*p2.p3*FR1^2+8*p1.p4*p1.p3*p1.p2*FFR2^2-4*p1.p4*p1.p3*FFR2*FR1*
 Mb-4*p2.p3*p1.p3*FFL2*FR1*Mtop+8*p1.p3^2*FFR2*FFL2*Mb*Mtop+4*p1.p3*p1.p2*
 FFR2*FL1*Mtop-p1.p3*FR1*FL1*Mb*Mtop$
denominator:=propDen(-p1+p3,MW,0)^2$


addToSum()$
finishSum();
End$
