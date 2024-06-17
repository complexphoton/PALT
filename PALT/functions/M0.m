function [y,E0]=M0(w,x,D,w_ab,rp,Green)
dx = x(2)-x(1);
Rp = rp./(1i*rp+w); 
chi = D(x)*Rp;
G=Green(x,x,w+w_ab);G(:,[1,end])=G(:,[1,end])/2;
T=G.*chi*dx;

[E0,y]=eigs(T-eye(length(T)),1,'smallestabs');