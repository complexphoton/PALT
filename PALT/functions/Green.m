function [G,detZ,MR,ML] = Green(x,x0,w,na,L,R)
dL=L(1,:); nL=L(2,:).*sqrt(1+1i*L(3,:)/w);
dR=R(1,:); nR=R(2,:).*sqrt(1+1i*R(3,:)/w);

% vR = [-1i,1];
% for m=length(dR):-1:1
%     p=w*dR(m)*nR(m);
%     vR=vR*[cos(p),sin(p)/nR(m);-sin(p)*nR(m),cos(p)];
% end
% 
% vL = [1i,1];
% for m=1:length(dL)
%     p=w*dL(m)*nL(m);
%     vL=vL*[cos(p),-sin(p)/nL(m);sin(p)*nL(m),cos(p)];
% end
% 
% pR=w*(1-x0)*na;
% pL=w*x0*na;
% 
% Z11 = vR(1)*cos(pR)-vR(2)*sin(pR)*na;
% Z12 = vR(2)*cos(pR)+vR(1)*sin(pR)/na;
% Z21 = vL(1)*cos(pL)+vL(2)*sin(pL)*na;
% Z22 = vL(2)*cos(pL)-vL(1)*sin(pL)/na;
% 
% detZ = Z11.*Z22-Z21.*Z12;
%--------------------------------------------------------
MR = eye(2,2);
for m=1:length(dR)
    p=w*dR(m)*nR(m);
    MR=[cos(p),sin(p)/nR(m);-sin(p)*nR(m),cos(p)]*MR;
end

ML = eye(2,2);
for m=1:length(dL)
    p=w*dL(m)*nL(m);
    ML=[cos(p),sin(p)/nL(m);-sin(p)*nL(m),cos(p)]*ML;
end

vR=[-1i,1]*MR;
vL=ML*[1;-1i];
detZ = vR*[cos(w*na),sin(w*na)/na;-sin(w*na)*na,cos(w*na)]*vL;

pR=w*(1-x0)*na;
pL=w*x0*na;

Z11 = vR(1)*cos(pR)-vR(2)*sin(pR)*na;
Z12 = vR(2)*cos(pR)+vR(1)*sin(pR)/na;
Z22 = vL(1)*cos(pL)+vL(2)*sin(pL)/na;

G0 = Z22.*Z12./detZ;
dG = -Z22.*Z11./detZ;

p_x=(x*ones(1,length(x0))-ones(length(x),1)*x0);
P=w*p_x*na;

G=cos(P).*G0+sin(P).*(dG+(p_x<0))/na;
G=G*w;
