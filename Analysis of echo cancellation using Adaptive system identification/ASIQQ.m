clear
% ASI.m: a program for Adaptive System Identification
load hyb2

n=1000000; % no. of sample points
N=98; % filter length
h=zeros(1,N); % initial filter coefficients
%h=[zeros(1,(N/2)-1) .5 .5 zeros(1,(N/2)-1)];
px=zeros(1,N); % dummy vector
mu=0.0002; % convergence factor

x=sqrt(12).*(rand(1,n)-.5); %input as a random signal
%[x,fs] = audioread('rec_1.wav'); % input as a voice signal
v =(1/10).*x;%noise

%The unknown system
%------------------
%b=fir1(20,.2);
b=hyb2;
N=98;
d=filter(b,1,x);%echo
new=d+v;%echo+noise
[y,h,Er]= ELMSnewQQ(x,new,h,mu,px);

% check the result:
%------------------
stem(h)
hold
stem(b,'*')
hold
gtext('o the adaptive filter')
gtext('* the unknown system')
xlabel('sample number, n')
ylabel('Amplitude')

pause
'press any key to continue';
plot(20*log10(abs(Er)));
xlabel('sample number, n')
ylabel('Error, dB')
title('learning curve')