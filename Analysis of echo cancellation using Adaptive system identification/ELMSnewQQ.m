function [x,h,Er]=ELMSnewQQ(x,d,h,mu,px) 
% This function performs adaptive normalized LMS algorithm 
% x--the input signal
% d--the desired response
% h--the coefficients of the adaptive filter
% mu-the convergence parameter

h=h(:);
% px-a dummy vector
N=length(h);
n=length(x);
Er=[];
for k=1:n;
	px=[x(k) px(1: N-1)];
	
	y(k)=px*h;
	E=(d(k)-y(k));
    ES = sign(E);
    ER=PO2Error(E);%Error with power of two quantizarion.
	Er=[Er E];
	h=h+mu*ER*px';%Put here, ES for observing sign error, ER for error with PO2 and E for error.
    
end