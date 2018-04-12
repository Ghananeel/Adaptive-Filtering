function [ W,xp ] = rls( u,d,m)
lambda=0.5;
N = length(u);
delta=100;
W = [];
w = zeros(m,1);
P=eye(m)*delta; %reset filter variable between monte carlo runs
for n=1:N-m,
    x=u(n:n+m-1); 
    x_flip=flipud(x);
    y=d(n+m-1);
    xp = w'*x_flip; % Predicted Output
    e(n) = (y-xp);%(x_flip'*x_flip+0.01); % Compute error
    nmse(n)=(e(n)*e(n))/(x_flip'*x_flip+0.01); % Compute Normalized Mean-Sq Error
    kappa=lambda^(-1)*P*x_flip/(1+lambda^(-1)*x_flip'*P*x_flip); % update kappa as perRLS
    w=w+(kappa'*nmse(n))'; % update weights
    P=lambda^(-1)*P-lambda^(-1)*x_flip'*kappa*P; %update as per R
    W(1:m,n)=w;

end
str=sprintf('RLS learning curve for filter order %d',m);

figure
%plot(nmse)

plot(W');

title(str);
xlabel('Iterations')
ylabel('NMSE')
end