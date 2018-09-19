function [xstar,mm]=QNclas(x0,n,FG,pp,w1,w2,u1,u2)
% this is a quasi-Newton algorithm for prog. QMFclassic.m 
% Reference: A. Antoniou, Digital Filters, Analysis, Design, and Applications
% 2nd Ed., McGraw Hill, To appear on November 1992.
% written by : Esam Abdel-Raheem, 28/6/92
%--------------------------------------------------------------------------
% x0:the starting point
% n: no. of variables
% F:objective fn.
% G:the Gradient
% The outputs from the program is : xstar,the optimum point, and mm,
% the final count of function evaluations.
% tolerances
ro=0.1;sigma=0.7;taw=0.1;exi=0.75;Mmax=600;epsilon1=10e-6;epsilon2=10e-10;

%Initialize Algorithm
%---------------------
k=0;m=0;
ni=1;
eval(['[f(1),g]=' FG '(x0,pp,u1,u2,w1,w2);']);
x=x0;
X=x;
%eval(['g=' G '(x0,pp,u1,u2,w1,w2);']);
%g=gra(x);
S=eye(n);
while k<1000,
%'I am here'
k;
d=-S*g;
[alpha,mm]=ALSCLAS(x,d,FG,pp,m,w1,w2,u1,u2);
m=mm;
delta=alpha*d;
x=x+delta;        %the new value of x
%check termination criteria and output results
%------------------------------------------------
X=[X x];
ni=ni+1;
eval(['[f(ni),gf]=' FG '(x,pp,u1,u2,w1,w2);']);
deltaf0=(f(ni-1)-f(ni));
if ((norm(delta)<epsilon1) &(abs(deltaf0)<epsilon1)|(m>Mmax)),
xstar=x
k
eval(['[fxstar,gxstar]=' FG '(xstar,pp,u1,u2,w1,w2);']);
m;
fxstar
break
end

%prepare for the next iteration
%-------------------------------
eval(['[f1,g1]=' FG '(x,pp,u1,u2,w1,w2);']);
%g1=gra(x);
gama=g1-g;
g=g1;
D=(delta)'*gama;
if D<=0,
S=eye(n);
else
d1=delta;
g11=gama;g12=g11*g11';d2=g11'*S*g11;d3=g11'*d1;
S=S+(d1*d1')/(d3)-(S*g12*S)/(d2); %DFB
%u1=gama'*S*gama;u2=gama'*delta;u3=delta*delta';u4=gama'*delta;
%u5=delta*gama'*S;u6=S*gama*delta';u7=gama'*delta;
%S=S+(1+u1./u2)*(u3./u4)-((u5+u6)./u7); %BFGS
end
k=k+1;
end
