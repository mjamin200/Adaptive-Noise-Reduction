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
disp("please wait until all sounds play")
%% Noisy signal

x=s+n;
sound(x,Fs)

%%  FIR Filter with unit variance noise
% N :length of filter
% M : length of input signal
% alpha : learning rate
% e : errors
% w : weights of filter
M = length(s);
N = 10;

[~,e]=RLS(n0,x,N,M);
pause(8)
sound(e,Fs)

%% FIR Filter with a noise of variance 10

n0 = sqrt(10).*randn(size(s));
n=filter(b,a,n0);

pause(12)
x=s+n;
sound(x,Fs)

[~,e]=RLS(n0,x,N,M);
pause(8)
sound(e,Fs)

%% IIR filter

a= [1,0.5];
b=[1,-0.9];

n=filter(b,a,n0);


x=s+n;
[w,e]=RLS(n0,x,N,M);
pause(8)
sound(e,Fs)

figure(1)
plot(x)
title('noisy signal')
figure(2)
plot(e)
title('out signal');
disp("The difference between LMS (or VSLMS) and RLS algorithm is that the LMS is faster than RLS but the output quality of RLS is better than LMS" + ...
    "and the number of iterations the algotithm need  to converge in the RLS is less than LMS ")

%% RLS algorithm

function[w,z]=RLS(inputs,d,N,M)
% z : error
% N :length of filter
% M : length of input signal
    z=zeros(1,M-N+1);
    w=zeros(1,N);
    lambda=0.999;
    delta= 1e-10;
    
    p=delta*eye(N);

    for i=N:M-1
        u=(inputs(i:-1:i-N+1))';
        y=dot(w,u);
        z(i-N+1)=d(i)-y;
        k=(p*u')/(lambda+u*p*u');
        w=w+k'*conj(z(i-N+1));
        p=(p -k*conj(u)*p)/lambda;
        
    end
    
end

                  