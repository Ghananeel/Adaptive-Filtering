function [y_predicted,W, e] = lms(u, d, m);
    % Maximum number of time step that can be predicted
    N = length(u);
    y_predicted=[];
    step = 0.1;
    % Intializatize weight matrix and associated parameters for LMS predictor
    w = zeros(m,1);
    W = [];
    x=[];
    %Transforming the input
    for k=1:N
        x(:,k)=u(k:k+m-1);
    end

    % LMS algortihm
    for n=1:N
        % Predict next sample and error
        y_predicted(n) = w'*x(:,n);
        e(n) = (y-y_predicted(i));
        w = w + step * e(n) * x(:,n);
        W(1:m,n)=w;  
    end; % for n t
    
    mse = e'e/N;
    %str=sprintf('LMS NMSE for filter order %d',m);
    %plot(nmse);
    %title(str);
    %xlabel('Iterations')
    %ylabel('NMSE')
    %W=W';
    %plot(W(:,2));
    %figure
    %axis([])

end
