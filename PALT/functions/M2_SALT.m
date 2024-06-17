function [y,E0,E1,P0,P1]=M2_SALT(In,D,w_ab,rp,Green)
w0 = In(1,1);I0=In(1,2);
w1 = In(2,1);I1=In(2,2);
In(1,1)=0; In(1,2)=1;
In(2,1)=0; In(2,2)=1;
E0=1i*In(1,1:2:end)+In(1,2:2:end);
E1=1i*In(2,1:2:end)+In(2,2:2:end);

x=linspace(0,1,length(E0));
Rp0 = rp./(1i*rp+w0); 
Rp1 = rp./(1i*rp+w1); 
N = D(x)./(1+I0*abs(Rp0*E0).^2+I1*abs(Rp1*E1).^2);

P0 = Rp0*N.*E0; 
P1 = Rp1*N.*E1;

G0 = Green(x,x,w0+w_ab);
G1 = Green(x,x,w1+w_ab);

Eg0 = trapz(x,G0.*P0,2).';
Eg1 = trapz(x,G1.*P1,2).';

Y=[Eg0-E0;Eg1-E1];
y=[real(Y),imag(Y)];