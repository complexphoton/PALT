function [M11,M12,M21,M22]=F0_trM(w,n,l)
M11=cos(n.*w*l); 
M12=sin(n.*w*l)./n;
M21=-n.*sin(n.*w*l);
M22=cos(n.*w*l);