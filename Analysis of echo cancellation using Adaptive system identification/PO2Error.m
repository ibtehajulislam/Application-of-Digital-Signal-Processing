function [ Q ] = PO2Error (E)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
for B=10
    tau=0;
end

if abs(E)>= 1
    Q = sign(E);
     
elseif power(2, -B +1)< abs(E) && abs(E)<1
    Q = (sign(E) * 2^floor(log2(abs(E))));
    
else
    Q = tau * sign(E);
end
