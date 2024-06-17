function [G,detZ] = Green_full(x,x0,w,na,L,R)
dL=L(1,:); nL=L(2,:).*sqrt(1+1i*L(3,:)/w);
dR=R(1,:); nR=R(2,:).*sqrt(1+1i*R(3,:)/w);
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

%middle
xM=x((x>=0)&(x<=1));
p_x=(xM*ones(1,length(x0))-ones(length(xM),1)*x0);
P=w*p_x*na;

GM=cos(P).*G0+sin(P).*(dG+(p_x<0))/na;

%left
xL=x(x<0);GL=zeros(length(xL),length(x0));
VL = [cos(pL).*G0-sin(pL).*(dG+1)/na;
       cos(pL).*(dG+1)+sin(pL).*G0*na];

XL=[0,cumsum(dL)]-sum(dL);

for m=length(dL):-1:1
    p=w*dL(m)*nL(m);
    px=w*nL(m)*(xL((xL>=XL(m))&(xL<XL(m+1)))-XL(m+1));
    GL((xL>=XL(m))&(xL<XL(m+1)),:)=[cos(px),sin(px)/nL(m)]*VL;
    VL=[cos(p),-sin(p)/nL(m);sin(p)*nL(m),cos(p)]*VL;
end
px=w*(xL(xL<XL(1))-XL(1));
GL(xL<XL(1),:)=[cos(px),sin(px)]*VL;

%Right
xR=x(x>1);GR=zeros(length(xR),length(x0));
VR = [cos(pR).*G0+sin(pR).*dG/na;
       cos(pR).*dG-sin(pR).*G0*na];

XR=[0,cumsum(dR)]+1;

for m=1:length(dR)
    p=w*dR(m)*nR(m);
    px=w*nR(m)*(xR((xR>XR(m))&(xR<=XR(m+1)))-XR(m));
    GR((xR>XR(m))&(xR<=XR(m+1)),:)=[cos(px),sin(px)/nR(m)]*VR;
    VR=[cos(p),sin(p)/nR(m);-sin(p)*nR(m),cos(p)]*VR;
end
px=w*(xR(xR>XR(end))-XR(end));
GR(xR>XR(end),:)=[cos(px),sin(px)]*VR;

G=[GL;GM;GR]*w;