function f=F0_eigf(w,er,sigma,d,D0,w_ab,rp)

% w = w(1,:)+1i*w(2,:);

Ln = length(d);
v1 = 1; v2 = -1i; 

for n=1:Ln
    g = D0(n)*rp./(w-w_ab+1i*rp);
    cn=sqrt(er(n)*(1+1i*sigma(n)./w)+g);
    [M11,M12,M21,M22]=F0_trM(w,cn,d(n));
    x1 = M11.*v1 + M12.*v2;
    v2 = M21.*v1 + M22.*v2;
    v1 = x1;
end

f = -1i.*v1+v2;

% f = [real(f);imag(f)];

