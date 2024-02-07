clear;
clc;
close all;
%%
% Mohammad Javad Amin 401211193
% Problem 2

%% Read audio file

[y,Fs] = audioread('s1.wav');
s=y(:,1);                          % Extrac 1 channel form audio
s = s - mean(s);                   % Normalize s
s = s/std(s);

%% Noise

n0=randn(size(s));
b=rand(1,5);
a=1;
n=filter(b,a,n0);
disp("please wait until all souds played")

%% Noisy signal

x=s+n;
sound(x,Fs)

%%  FIR Filter with unit variance noise
% N :length of filter
% M : length of input signal
% alpha : mu tilde
% e : errors
% w : weights of filter
M = length(s);
N =10;

alpha = 1;
[~,e]=NLMS(n0,x,N,alpha,M);
pause(14)
sound(e,Fs)

%% FIR Filter with a noise of variance 10

n0 = sqrt(10).*randn(size(s));
n=filter(b,a,n0);

pause(14)
x=s+n;
sound(x,Fs)

[~,e]=NLMS(n0,x,N,alpha,M);
pause(14)
sound(e,Fs)

%% IIR filter

a= [1,0.5];
b=[1,-0.9];

n=filter(b,a,n0);


x=s+n;
[w,e]=NLMS(n0,x,N,alpha,M);
pause(14)
sound(e,Fs)

figure(1)
plot(x)
title('noisy signal')
figure(2)
plot(e)
title('out signal');

%% NLMS algorithms

function[w,e]=NLMS(inputs,d,N,alpha,M)
% e : error
% u_temp : because LMS run when the first sample arrive, we put M-1 zeros in beging of inputs, if whe don't put this zeros we must wait to m sample arrive
    u_temp=[zeros(1,N-1),inputs'];   
    e=zeros(1,M);
    w=zeros(1,N);
    for i=N:M
        u=u_temp(i:-1:i-N+1);
        y=dot(w,u);
        e(i-N+1)=d(i-N+1)-y;
        w =  w + (alpha/(norm(u)^2))*e(i-N+1)*u;
    end
    w=w';
end
                  