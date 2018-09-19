clear all
clf
% QMF_ibt.m: A program to design classic 2-ch QMF filter banks.
% Source: Johnson, A filter family in use for QMFB
format long

ws=.6*pi;wp=.4*pi;N=34;
ALFA=1; % a weight constant used in the function and gradient in FGCLAS.m
% If you donot know how to use it, consider it as 1
save ALFA ALFA
%consider a starting point as the impulse response from an
%FIR filter designed using window technique with the same
%edge frequencies
x0=fir1(N-1,.5);
x0=2.*x0(N/2+1:N); % note that we are designing a linear phase filter and 
%hence we only need half the coefficients to proceed in forming the 
%frequency response

x0=x0';

M=N/2; % This is the length of the vector
n=length(x0); % basically M=n
%n =number of elements in the unknown vector (half the coefficients)
% to minimize the cost function the function and the gradient for the 
%cost function to be used with the Quasi Newton Algorithm
FG='FGCLAS2'
%FGCLAS.m is a function to provide the function and the gradient for the
%quasi Newton algorithm
                                  
%---------------------------------------------
%form the ideal amp. vector and the freq. vector
% form w1 and w2, the frequency points in the band 0-pi/2 and ws-pi,
% respectively..the number of frequency points in the total band (0 pi)
%would be equal to 6-8 the filter length...put w1 and w2 in the form of an
%array.
%provide u1 and u2, the desired mag. res. for the channel and the stopband
%regions. Note that u1 and u2 are 2 arrays with the same dimensions of w1
%and w2, respectively
L1=125; %No. of points in (0-pi/2) interval.
Lm1=L1-1;
L2=100; %no. of points in the stop band.
Lm2=L2-1;
q1=0:Lm1;
w1=(pi/(2*Lm1))*q1;
q2=0:Lm2;
w2=ws+((pi-ws)/Lm2)*q2;
u1=ones(1,L1);
u2=zeros(1,L2);

pp=2;
%-------------------------------------------------------
%the quasi Newton optimization algorithm
%_____________________________________________
[xstar,m]=QNCLAS(x0,n,FG,pp,w1,w2,u1,u2)
xstar
m
% m is the number of iterations to arrive at the solution
%form the impulse response out of the optimal solution, xstar
b1=xstar';
b=.5.*[fliplr(b1) b1]; 

%impulse response of the low pass filter in the analysis bank
[hh,ii]=freqz(b,1,2048);
plot(ii/(2*pi),20*log10(abs(hh)))
grid
ylabel("gain, db")
xlabel(" normalized frequency (LPF in the Analysis Bank)")
%xlabel(" Impulse response of the low pass filter")
pause
h0w=b;

%___________________________%

%impulse responses of the HPF in the analysis bank
for i=1:length(h0w) 
    h1w(i)=(-1)^(i-1)*h0w(i);
end
c=h1w;
[hh1,i2]=freqz(c,1,2048);
plot(i2/(2*pi),20*log10(abs(hh1)))
grid
ylabel("gain, db")
xlabel(" normalized frequency (HPF in the Analysis Bank)")
pause
%___________________________%


%impulse responses of the LPF in the synthesis bank
for u1=1:length(h1w);
    f0w(u1)=(-1)^(u1-1)*h1w(u1);
end
%d=f0w;
%[hh2,i3]=freqz(d,1,2048);
%plot(i3/(2*pi),20*log10(abs(hh2)))
%grid
%ylabel("gain, db")
%xlabel(" normalized frequency (LPF in the Synthesis Bank)")
%pause
%___________________________%


%impulse responses of the HPF in the synthesis bank
for u1=1:length(h0w);
    f1w(u1)=-(-1)^(u1-1)*h0w(u1);
end
%e=f1w;
%[hh3,i4]=freqz(e,1,2048);
%plot(i4/(2*pi),20*log10(abs(hh3)))
%grid
%ylabel("gain, db")
%xlabel(" normalized frequency (HPF in the Synthesis Bank)")
%pause

%Impulse response in the Analysis bank
%AB=h0w+c;
%tAB=conv(h0w,h1w);
%[hhAB,iAB]=freqz(tAB,1,2048);
%plot(iAB/(2*pi),20*log10(abs(hhAB)))
%grid
%ylabel("gain, db")
%xlabel(" normalized frequency (Analysis Bank Response)")
%pause
%___________________________%

% Channel Response for Low pass filter
%cLPF=conv(h0w,f0w);
%[hLPF,iLPF]=freqz(cLPF,1,2048);
%plot(iLPF/(2*pi),20*log10(abs(hLPF)))
%grid
%ylabel("gain, db")
%xlabel(" normalized frequency (Channel Response for Low pass filter)")
%pause
%___________________________%

% Channel Response for Low pass filter
%cHPF=conv(h1w,f1w);
%[hHPF,iHPF]=freqz(cHPF,1,2048);
%plot(iHPF/(2*pi),20*log10(abs(hHPF)))
%grid
%ylabel("gain, db")
%xlabel(" normalized frequency (Channel Response for High pass filter)")
%pause
%___________________________%

%Test for reconstruction
x=0:100;
x=.01*[zeros(1,(N-1)),x];
m=length(x);
t=conv(h0w,f0w)+conv(h1w,f1w);
L=t;
opww=filter(t,1,x);
Xs=x(1:length(x)-(N-1));
aha1=opww(N:length(opww));
errw=Xs-aha1;
mse=errw*errw';
msew=mse/length(Xs);
DB1=10*log10((Xs*Xs')/mse);
[H1n,ii]=freqz(t,1,2048);
LW1=(20*log10(abs(H1n)));
i=ii/(2*pi);

plot(i/(2*pi),LW1)
grid
ylabel("gain, db")
xlabel(" normalized frequency (Channel Response)")
%the impulse response of the lowpass filter
% Let h0w=b as the impulse response of the LPF in the analysis bank
%Find the impulse responses of the HPF in the analysis bank, h1w, LPF in
%the synthesis bank, f0w, and the hpf of the synthesis bank, f1w
% Draw the analysis bank response and the channel response.
% apply a 100 point ramp function to the analysis bank and get the output
% of the synthesis bank... find either the Signal to reconstruction noise
% ratio of the MSE of the reconstructed signal by comparing the two signals



% but note that the reconstructed signal has a delay N-1 samples from the
% original signal (due to the filter bank system) and hence it would be
% useful to pad the ramp with number of zeros equal to that delay

