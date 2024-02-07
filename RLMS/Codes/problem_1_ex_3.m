clear;
clc;
close all;
%%
% Mohammad Javad Amin 401211193
% Problem 1 , exercise 3

%% definition
% d : desired signal
% N :length of filter
% M : length of input signal
% e : errors
% w : weights of filter
% l : noise amplitude
% d_t : corrupted desired signal

a=[1,0.5];
b=[1,-0.9];         % impulse response
inputs=randn(1,300);
d=filter(b,a,inputs);    
M=length(inputs);

%% part a

N = 4;
k=5;
m_error=zeros(1,M-N+1);

for i=1:k
    [w,cost]=RLS(inputs,d,N,M);
    m_error=m_error+cost;
end
m_error=m_error/5;

disp("weights for N=4 :");
disp(w');


figure
plot(m_error);
title('squared error, N=4 ');
xlabel('Number of iterations');

%% part b

N = [2,3,5,7,10];
for i=N
   
    m_error=zeros(1,M-i+1);

    for g=1:k
        [w,cost]=RLS(inputs,d,i,M);
        m_error=m_error+cost;
    end 
    m_error=m_error/5;

    disp(['weights forand N=',num2str(i),':']);
    disp(w');
       
    figure
    plot(m_error);
    title(['squared error, N=',num2str(i)]);
    xlabel('Number of iterations');
end

%% part c
l = [0.1,0.3,1];
N = 4;
v = randn(1,300);
k=5;
m_error=zeros(1,M-N+1);


for g=l
    d_t=d+g*v;

    for i=1:k
        [w,cost]=RLS(inputs,d_t,N,M);
        m_error=m_error+cost;
    end
    m_error=m_error/5;
    
    disp(['weights for N=4 and l=',num2str(g), ': ']);
    disp(w')
    
    figure
    plot(m_error);
    title(['squared error, N=4 and l=', num2str(g),' : ']);
    xlabel('Number of iterations');
end
disp(" The LMS is more quicker than RLS algorithm but the error in RLS is much better than LMS and" + ...
    " the number of iterations the algotithm need  to converge in the RLS is less than LMS ");
%% RLS algorithms

function[w,cost,J_min,J_inf]=RLS(inputs,d,N,M)
% z : error
% N :length of filter
% M : length of input signal

    z=zeros(1,M-N+1);
    w=zeros(1,N);
    lambda=0.5;
    delta= 1e-10;
    
    p=delta*eye(N);
   
    for i=N:M-1
        u=inputs(i:-1:i-N+1);
        y=dot(w,u);
        z(i-N+1)=d(i)-y;
        k=(p*u')/(lambda+u*p*u');
        w=w+k'*conj(z(i-N+1));
        p=(p -k*conj(u)*p)/lambda;
        
    end
    cost=z.^2;
    J_min=min(z);
    J_inf=sum(z(M-N-19:M-N))/20;

end

