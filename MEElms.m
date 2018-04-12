function [W, e] = MEElms(u, d, m);
% Maximum number of time step that can be predicted
avg=0;
N = length(u);
mu=0.1;
lambda=0.05
sigma=4.5;
L=100;
V=zeros(m,1);
z=mean(d)
u=vertcat(zeros(L,1),u)
err=zeros(L,1)
a=[1 -1];
b=[1 0 0 0 0 0 0 0 0 -1];
w_opt=filter(b,a,[1,zeros(1,m-1)]);
w_opt=(w_opt)'
WSNR=0
WSNR_final=0

% Intializatize weight matrix and associated parameters for LMS predictor
w = zeros(m,1);
W = [];

for n=201:N-m+1
x=u(n:n+m-1)
x_flip=flipud(x)
y=d(n+m-1-100)
if n == 103
    display(n);
end
% Predict next sample and error
xp = w'*x_flip+z;
e= (y-xp);%/(x_flip'*x_flip+0.01);
%nmse(n)=(e(n)*e(n))/(x_flip'*x_flip+0.01);
% Adapt weight matrix ans step size
%err=vertcat(err,e);
err(n)=e;
gaussian=zeros(m,1);
for i=n-L:n-1;
    x_i=u(i:i+m-1);
    error=err(i)-err(n);
    x_diff=flipud(x-x_i);
    new_error=error*error*x_diff;
    gaussian=gaussian+normpdf(error,avg,sigma)*error*x_diff;    
end    
V=(1-lambda)*V+(lambda/(sigma*sigma*L))*gaussian;
w = w - mu * V;


W(1:m,n-L+1)=w;
end; % for n t
WSNR=10*log10(((w_opt)'*w_opt)/((w_opt-w)'*(w_opt-w)));

str=sprintf('MEE LMS weight tracks for filter order %d',m);
figure
plot(W')
%plot(nmse);
title(str);
xlabel('Iterations')
ylabel('Weights')
%plot(e.^2)
end