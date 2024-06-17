function E=M1_FieldPrf_SALT(x,Ex0,D0,w,w_ab,rp,Green_full)
Nx0=length(D0);
x0=linspace(0,1,Nx0);
Rp=rp/(1i*rp+w);
G = Green_full(x,x0,w+w_ab);
E = trapz(x0,Rp*G.*(D0.*Ex0),2).';


