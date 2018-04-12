function y = Kernel_lms(u,d,m)
% Learning Rate for KLMS
lr_k = .3;
%   lr_k = .6;

N_tr=length(u);

%init
e_k = zeros(N_tr,1);        % Error
y = zeros(N_tr,1);          % Predicted Output
mse_te_k = zeros(N_tr,1);   % Mean Squared Error

T=d(1:N_tr-m+1);

for k=1:N_tr-m+1
    X(:,k)=u(k:k+m-1);
end

% n=1 init
e_k(1) = T(1);
y(1) = 0;

% start KLMS
for n=2:N_tr-m+1
    % training
    ii = 1:n-1;
    % Compute the predicted output in Gaussian Kernel Space
    y(n) = lr_k*e_k(ii)'*(exp(-sum((X(:,n)*ones(1,n-1)-X(:,ii)).^2)))';
    % Compute Error between true and predicted output
    e_k(n) = T(n) - y(n);
    % Compute Normalized Mean-Squared Error
    nmse(n)=(e_k(n)*e_k(n))/(u(:,n)*u(:,n)'+0.01);
    
end
%str=sprintf('LMS NMSE for filter order %d',m);
% Learning Curve
plot(nmse);
title('Learning for KLMS');
xlabel('Iterations')
ylabel('NMSE')

end
