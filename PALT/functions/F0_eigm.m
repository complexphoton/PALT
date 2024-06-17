function [v,E]=F0_eigm(x,w,cn,d)
cn=[cn,1];
Ln = length(d);
X = [0,cumsum(d)]; 
Lx = zeros(1,Ln+1);Lx(1)=sum(x<=0);
v = zeros(2,Ln+1);
v(:,1)=[1;-1i*w];
for n = 1:Ln
    phase=w*cn(n)*d(n);
    M=[cos(phase),sin(phase)/(cn(n)*w);
       -sin(phase)*cn(n)*w, cos(phase)];
    v(:,n+1)=M*v(:,n);
    Lx(n+1)=sum(x<=X(n+1));
end
Lx=[Lx,length(x)];
E=zeros(1,length(x));
for n=1:Ln+1
    xn=x(Lx(n)+1:Lx(n+1))-X(n);
    E(Lx(n)+1:Lx(n+1))=v(1,n)*cos(w*cn(n)*xn)...
                      +v(2,n)*sin(w*cn(n)*xn)/(cn(n)*w);
end
E(1:Lx(1))=exp(-1i*w*x(1:Lx(1)));
