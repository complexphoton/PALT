function [y,E1]=M1_E1_SALT(w1,D,w0,E0,w_ab,rp,Green)
M=length(E0);
xm = linspace(0,1,M);

Rp0 = rp./(1i*rp+w0); 
Rp1 = rp./(1i*rp+w1); 

N0=D./(1+abs(Rp0*E0).^2);
G = Green(xm,xm,w1+w_ab); G(:,[1,end])=G(:,[1,end])/2;

dx=xm(2)-xm(1);

T = Rp1*G.*N0*dx;

[E1,y]=eigs(T-eye(M),1,'smallestabs');