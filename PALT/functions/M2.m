function [y,Eg,P,NR]=M2(In,D,w_ab,rp,rl,Green)
[K,N] = size(In);N=N/2;
x=linspace(0,1,N);
K0 = (K+1)/2; K1 = K0+1; 
M = K0-1; 
w0 = In(K0,1); dw = In(K1,1)*(-M:M)';
w = dw + w0;

I0 = In(K0,2); 
In(K0,1)=0; In(K1,1) = 0;In(K0,2)=1; 
E = 1i*In(:,1:2:end) + In(:,2:2:end);

P = zeros(K,N); Eg = zeros(K,N); 

Rp = diag(rp./(w+1i*rp)); Rn=rot90(Rp,2);
Rl = diag(rl./(dw+1i*rl));Rl(K0,K0)=-1i;

NR=zeros(K,N); %comment this line

for n = 1:N
En = spdiags((flipud(E(:,n)).').*ones(K,K),-M:M,K,K);
% Ep = fliplr(En);
DD = sparse(K0,1,D(n),K,1);
NN = (eye(K,K)-0.5*I0*Rl*((En'*Rp*En)-En*(En*Rn)'))\DD; 
% Nn = spdiags((flipud(NN).').*ones(K,K),-M:M,K,K);
P(:,n) = Rp*En*NN; 
NR(:,n)= NN; %comment this line
end

for k =1:K
G = Green(x,x,w(k)+w_ab);
Eg(k,:) = trapz(x,G.*P(k,:),2).';
end

Y = Eg-E;
y = [real(Y),imag(Y)];


