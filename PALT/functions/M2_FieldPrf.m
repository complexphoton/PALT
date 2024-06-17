function [E,P,D,B]=M2_FieldPrf(x,D0,Px0,w0,wr,w_ab,rp,rl,Green_full)
[K,Nx0] = size(Px0);K0=(K-1)/2;
x0=linspace(0,1,Nx0);
dw = wr*(-K0:K0)'; w = dw + w0;

Nx=length(x); 
E=zeros(K,Nx); D=zeros(K,Nx); P=zeros(K,Nx);

for k =1:K
G = Green_full(x,x0,w(k)+w_ab);
E(k,:) = trapz(x0,G.*Px0(k,:),2).';
end

B=1i*(E(:,2:end)-E(:,1:end-1))./(x(2:end)-x(1:end-1))./(w+w_ab);

Nx1=find(x>=0,1,'first');Nx2=find(x<=1,1,'last');

Rp = diag(rp./(w+1i*rp));
Rl = diag(rl./(dw+1i*rl));Rl(K0+1,K0+1)=-1i;

% D=zeros(K,Nx2-Nx1+1); %comment this line
% P=zeros(K,Nx2-Nx1+1);

for n = Nx1:Nx2
En = spdiags((flipud(E(:,n)).').*ones(K,K),-K0:K0,K,K);
Ep = fliplr(En);
d0 = sparse(K0+1,1,D0(x(n)),K,1);
D(:,n) = (eye(K,K)-0.5*Rl*((En'*Rp*En)-Ep*(Ep*Rp)'))\d0; 
DD = spdiags((flipud(D(:,n)).').*ones(K,K),-K0:K0,K,K);
P(:,n) = Rp*DD*E(:,n); 
end



