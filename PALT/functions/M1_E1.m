function [y,E1]=M1_E1(In,Dp,w0,E0,w_ab,rp,rl,Green)
dw = In(1)+1i*In(2);
M=length(E0);
xm = linspace(0,1,M);
dW=[-dw',0,dw]; w = w0+dW;
Rp = rp./(1i*rp+w); Rl=rl/(1i*rl+dw);

ch1=0.5*Rl*(Rp(3)-Rp(2)')./(1-0.5*Rl*(Rp(3)-Rp(1)')*abs(E0).^2);
ch2=0.5*Rl*(Rp(2)-Rp(1)')./(1-0.5*Rl*(Rp(3)-Rp(1)')*abs(E0).^2);

D0=Dp./(1+abs(Rp(2)*E0).^2);
c11 = Rp(3)*D0.*(1+ch1.*abs(E0).^2); c11([1,end])=0.5*c11([1,end]);
c12 = Rp(3)*D0.*ch2.*E0.^2; c12([1,end])=0.5*c12([1,end]);
c21 = Rp(1)'*D0.*ch1.*conj(E0).^2; c21([1,end])=0.5*c21([1,end]);
c22 = Rp(1)'*D0.*(1+ch2.*abs(E0).^2); c22([1,end])=0.5*c22([1,end]);

G1 = Green(xm,xm,w(3)+w_ab);
G2 = Green(xm,xm,w(1)+w_ab); G2=conj(G2);

dx=xm(2)-xm(1);

C11 = G1.*c11*dx;
C12 = G1.*c12*dx;
C21 = G2.*c21*dx;
C22 = G2.*c22*dx;

T = [C11,C12;C21,C22];

% T = C11+C12*((eye(M,M)-C22)\C21)

[E1,yc]=eigs(T-eye(length(T)),1,'smallestabs');
y = [real(yc),imag(yc)];
