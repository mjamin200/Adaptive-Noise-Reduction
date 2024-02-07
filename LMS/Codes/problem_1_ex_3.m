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
% alpha : learning rate
% e : errors
% w : weights of filter
% p : power of input signal
% l : noise amplitude
% d_t : corrupted desired signal

a=[1,0.5];
b=[1,-0.9];         % impulse response
inputs=randn(1,300);
d=filter(b,a,inputs);    
M=length(inputs);

%% part a

alpha = 0.3;
N = 4;
k=5;
m_error=zeros(1,M);

% calulate mu max for N=4
p= inputs*inputs'/M;
 alpha_max=2/(2*N*p);
 disp('mu max for N=4 and is :');
 disp(alpha_max);

for i=1:k
    [w,cost]=LMS(inputs,d,N,alpha,M);
    m_error=m_error+cost;
end
m_error=m_error/5;

disp("weights for mu=0.3 and N=4 :");
disp(w');
disp('if mu is begger than u_max may be LMS algorithm not converged')

figure
plot(m_error);
title('Mean squared error, N=4 , mu=0.3');
xlabel('Number of iterations');

%% part b

N = [2,3,5,7,10];
for i=N
    alpha_max=2/(3*i*p);
    disp(['mu max for N=',num2str(i),' and mu= 0.3 is :']);
    disp(alpha_max);
    
    m_error=zeros(1,M);

    for g=1:k
        [w,cost]=LMS(inputs,d,i,alpha,M);
        m_error=m_error+cost;
    end 
    m_error=m_error/5;

    disp(['weights for mu=0.3 and N=',num2str(i),':']);
    disp(w');
    disp('if mu is begger than u_max may be LMS algorithm not converged')
    
    figure
    plot(m_error);
    title(['Mean squared error, N=',num2str(i),' , mu=0.3']);
    xlabel('Number of iterations');
end

%% part c
l = [0.1,0.3,1];
N = 4;
v = randn(1,300);
alpha = 0.3;
k=5;

% calulate mu max for N=4
p= inputs*inputs'/M;
alpha_max=2/(2*N*p);
disp('mu max for N=4 and is :');
disp(alpha_max);
disp('if mu is begger than u_max may be LMS algorithm not converged')
for g=l
    m_error=zeros(1,M);
    d_t=d+g*v;
    
    
    for i=1:k
        [w,cost]=LMS(inputs,d_t,N,alpha,M);
        m_error=m_error+cost;
    end
    m_error=m_error/5;
    
    disp(['weights for mu=0.3 , N=4 and l=',num2str(g), ': ']);
    disp(w')
    
    figure
    plot(m_error);
    title(['Mean squared error, N=4 , mu=0.3 and l=', num2str(g),' : ']);
    xlabel('Number of iterations');
end

%% LMS algorithms

function[w,cost,J_min,J_inf]=LMS(inputs,d,N,alpha,M)
% e : error
% u_temp : because LMS run when the first sample arrive, we put M-1 zeros in beging of inputs, if whe don't put this zeros we must wait to m sample arrive
    u_temp=[zeros(1,N-1),inputs];   
    e=zeros(1,M);
    w=zeros(1,N);
    for i=N:M
        u=u_temp(i:-1:i-N+1);
        y=dot(w,u);
        e(i-N+1)=d(i-N+1)-y;
        w =  w + alpha*e(i-N+1)*u;
    end
    cost=e.^2;
    J_min=min(cost);
    J_inf=sum(cost(M-19:M))/20;

end


