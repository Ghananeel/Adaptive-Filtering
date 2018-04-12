function [wn,nmse,wsnr]=weiner(ip,desired,order,windowsize,lamda,wo)

L=size(desired,1);
wn=zeros(order,1);
c=0;
while(c<1000)
    
        i=floor(rand(1)*(L-windowsize)+1);
        R=ip(i:i+windowsize-1,:)'*ip(i:i+windowsize-1,:)+lamda*eye(order);
        w=(R^-1)*ip(i:i+windowsize-1,:)'*desired(i:i+windowsize-1);
%         e=desired(i:i+windowsize-1)-(ip(i:i+windowsize-1,:)*w);
        nmse(c+1)=nmsecal(ip(i:i+windowsize-1,:),desired(i:i+windowsize-1),w,lamda);
        wsnr(c+1)=wsnrCompute(wo,w);
    wn=wn+w;
    c=c+1;
end
wn=wn/c;

% % nmse=0;

% L=size(desired,1);
% wn=zeros(order,1);
% c=0;
% for i=1:windowsize:L
%     
%     if L-i+1<windowsize
% %          R=checkSingular(ip(i:L,:)'*ip(i:L,:),lamda);
%         R=ip(i:L,:)'*ip(i:L,:)+lamda*eye(order);
%          w=(R^-1)*ip(i:L,:)'*desired(i:L);
%     else
% %         R=checkSingular(ip(i:windowsize-1,:)'*ip(i:windowsize-1,:),lamda);
%         R=ip(i:windowsize-1,:)'*ip(i:windowsize-1,:)+lamda*eye(order);
%         w=(R^-1)*ip(i:windowsize-1,:)'*desired(i:windowsize-1);
%     end
%     wn=wn+w;
%     c=c+1;
% end
% wn=wn/c;

end

function nmsec=nmsecal(inp,ds,we,lamda)
nmsec=0;
for j=1:length(ds)
  e=ds(j)-inp(j,:)*we;
  nmsec=nmsec+e^2/(inp(j,:)*inp(j,:)'+lamda);
end
nmsec=nmsec/length(ds);
end
