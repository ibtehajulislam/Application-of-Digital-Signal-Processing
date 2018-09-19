function [y,zz]=ALCLAS(xx,dd,FG,pp,mm,w1,w2,u1,u2)
% This program evaluate the line search proposed by Prof. A. Antoniou.
% To be used for prog. FIROPT.m
%--------------------------------------------------------------------------
% Reference: A. Antoniou, Digital Filters, Analysis, Design, and Applications
% 2nd Ed., McGrow Hill, To appear on November 1992.
% written by : Esam Abdel-Raheem, 12/7/92
%--------------------------------------------------------------------------
% xx: initial point
% dd=-s*g: search direction
% F:Function
% G:Gradient
% pp:pnorm
% mm : No. of function evaluations.
% output: y=a0(alpha) which is close to a* (the true min. alpha)
%         zz=final no. of function evaluations.
%tolerances and limits
ro=0.1;sigma=0.7;taw=0.1;exi=0.75;Mmax=600;epsilon1=10e-6;epsilon2=10e-10;

%initialize algorithm
m=mm;
eval(['[f0,g0]=' FG '(xx,pp,u1,u2,w1,w2);']);f00=f0;deltaf0=f0;
%eval(['g0=' G '(xx,pp,u1,u2,w1,w2);']);
m=m+2;
aL=0;aU=10e99;
fL=f0;
eval(['[fs,gs]=' FG '(xx+aL*dd,pp,u1,u2,w1,w2);']);fLd=gs'*dd;

%estimate a0
if abs(fLd)> epsilon2,
a0=-2*deltaf0/fLd;
else a0=1;
end
if (a0 <=0 | a0>1),
a0=1;
end
noviol=1;
while noviol,
delta=a0*dd;
eval(['[f0,g0]=' FG '(xx+delta,pp,u1,u2,w1,w2);']);m=m+1;

%Interpolation
%---------------
s1=(f0>(fL+ro*(a0-aL)*fLd));
s2=(abs(fL-f0)>epsilon2);
s3=(m<Mmax);
if s1&s2&s3,
if a0<aU, aU=a0;
end
a0u=aL+((a0-aL)^2*fLd)/(2*(fL-f0+(a0-aL)*fLd));
a0uL=aL+taw*(aU-aL);
if a0u<a0uL, a0u=a0uL;
end
a0uU=aU-taw*(aU-aL);
if a0u>a0uU, a0u=a0uU;
end
[a0u a0uL a0uU];
a0=a0u;
else

eval(['[fss,gss]=' FG '(xx+a0*dd,pp,u1,u2,w1,w2);']);f0d=gss'*dd;
m=m+1;

%Extrapolation
%---------------

s1=(f0d<sigma*fLd);
s2=(abs(fL-f0)>epsilon2);
s3=(m<Mmax);

if s1&s2&s3,
deltaa=(a0-aL)*f0d/(fLd-f0d);
if deltaa<=0, a0u=2*a0;
else a0u=a0+deltaa;
end
a0uU=a0+exi*(aU-a0);
if a0u>a0uU,
a0u=a0uU;
end
aL=a0;a0=a0u;fL=f0;fLd=f0d;
else
     noviol=0;
   end
  end
end
y=a0;
z=y*dd;
zz=m;
