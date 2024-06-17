function [y,E0,D0]=M1_E0(In,Dp,w_ab,rp,Green)
w = In(1);I0=In(2);
In(1)=0; In(2)=1;
E0=1i*In(1:2:end)+In(2:2:end);
x=linspace(0,1,length(E0));
Rp = rp./(1i*rp+w); 
D0=Dp./(1+I0*abs(Rp*E0).^2);
P = Rp*D0.*E0;
G = Green(x,x,w+w_ab);
Eg = trapz(x,G.*P,2).';

Y=Eg-E0;
y=[real(Y),imag(Y)];

