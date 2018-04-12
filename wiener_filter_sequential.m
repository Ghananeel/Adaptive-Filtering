%m=order w_size=window size x=input y=output
function [nmse,WSNR_final]=wiener_filter_sequential(x,y,m,w_size)
n=length(x);
i=0;
j=0;
w=zeros(m,1);
k=0
b=[1 0 0 0 0 0 0 0 0 -1];           % Optimal Filter Coefficients for b
a=[1 -1];                           % Optimal Filter Coefficients for a
w_opt=filter(b,a,[1,zeros(1,m-1)]); % Optimal Filter Coefficients
w_opt=(w_opt)'
nmse=0
Lamda=0.01;
if w_size==10000
    win_start=[1];
else
    win_start=randperm(10000-w_size,100);
end
k=length(win_start);
WSNR=0
WSNR_final=0



weight_array=zeros(ceil(length(x)/w_size),m);

% Weiner Filter

for i=1:w_size:length(x)-w_size
    w_window=0
    X=x(i:i+w_size-1);
    
    Y=y(i:i+w_size-1);
    L=length(X);
    

    
    WSNR_new=0
    y_in=Y(m:L);

    for j=1:L-m+1
    x_in(j,1:m)=X(j:j+m-1);
    end
    

    % Wiener Filter main computation steps
    I=eye(m);
    R=(((x_in'*x_in)+Lamda*I)^-1); % Compute Autocorrelation inverse -->(R+lambda*I)^-1
    P=x_in'*y_in                   % Compute Cross-Correlation
    w_window=R*P;                  % Filter Coefficients obtained by Wiener Filter    
    weights_array(i,:)=w_window;
    
    %Compute WSNR between actual coefficients and those obtained by Wiener
    %Filtering
    WSNR_new=10*log10(((w_opt)'*w_opt)/((w_opt-w_window)'*(w_opt-w_window)));
    WSNR=WSNR+WSNR_new
    %Compute NMSE obtained by wiener filtering
    nmse_new=sum((y_in-((x_in)*w_window)).^2/(X'*X))/length(y_in)
    nmse=nmse_new+nmse
         
end

WSNR_final=WSNR/k
nmse=nmse/k
w=w./k


        
  
    