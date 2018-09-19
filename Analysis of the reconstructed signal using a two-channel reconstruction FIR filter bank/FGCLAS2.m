% FGCLAS2.m : a file to create the fn. and grad. for QMFclassic.m

function [errf,gra]=FGCLAS2(x,pp,u1,u2,w1,w2);
ALFA =1;
M=length(x);
for i=1:length(w1);
u=w1(i);
j=1:M;
co=cos((j-.5)*u);copi=cos((j-.5)*(u+pi));
Q1=co'*co+copi'*copi;
z1(i)=(x'*Q1*x);
s=2.*Q1*x;
G(:,i)=s;
end

for i=1: length(w2); 
    v=w2(i);
    k=1:M;
    co1=cos((k-.5)*v);
    Q2=co1'*co1;
    z2(i)=(x'*Q2*x);
    s2=2.*Q2*x;
    G2(:,i)=s2;
    
    
end
error2=sum(z2);
error1=z1'-u1';

gra1=2.*G*error1;

% Similarly, find the function and gradient related to stopband
% error,error2 and gra2. This has been left for you to try

errf=error1'*error1+ALFA.*error2;

gra2=sum(G2'); % you can get G2 similar to G
gra2=gra2';
gra=(ALFA).*gra2+gra1;


