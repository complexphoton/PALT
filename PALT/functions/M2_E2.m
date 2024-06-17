function [y,E1]=M2_E2(wf,w0,wr,E,P,N,w_ab,rp,rl,Green)
% wf = In(1)+1i*In(2);
[K,M]=size(E); k=(K-1)/2;
x = linspace(0,1,M); dx=x(2)-x(1);

W=wf+[-k:k]*wr;
R_p=rp./(w0+W+1i*rp); R_p = diag(R_p);
R_n=rp./(w0-W-1i*rp); R_n = diag(R_n);
Rl =rl./(W+1i*rl);    Rl = diag(Rl);

C11 = zeros(K,K,M); C12 = zeros(K,K,M);
C21 = zeros(K,K,M); C22 = zeros(K,K,M);

for m=1:M
    Eb = spdiags((flipud(E(:,m)).').*ones(K,K),-k:k,K,K);
    Pb = spdiags((flipud(P(:,m)).').*ones(K,K),-k:k,K,K);
    dd = spdiags((flipud(N(:,m)).').*ones(K,K),-k:k,K,K);
    Nf = eye(K,K)-0.5*Rl*(Eb'*R_p*Eb-Eb*R_n*Eb');
    X_p = 0.5*(Nf\Rl*(Eb'*R_p*dd-Pb'));
    X_n = 0.5*(Nf\Rl*(Pb-Eb*R_n*dd));
    C11(:,:,m) = dd + Eb*X_p; 
    C12(:,:,m) = Eb*X_n;
    C21(:,:,m) = Eb'*X_p;
    C22(:,:,m) = dd + Eb'*X_n; 
end

H11 = zeros(K*M,K*M); H12 = zeros(K*M,K*M);
H21 = zeros(K*M,K*M); H22 = zeros(K*M,K*M);

for m=1:K
    G_p = Green(x,x,W(m)+w0+w_ab);
    G_n = Green(x,x,W(m)-w0-w_ab);
    for n=1:K
        c11= reshape(C11(m,n,:),[1,M]); c11([1,end])=0.5*c11([1,end]);
        c12= reshape(C12(m,n,:),[1,M]); c12([1,end])=0.5*c12([1,end]);
        c21= reshape(C21(m,n,:),[1,M]); c21([1,end])=0.5*c21([1,end]);
        c22= reshape(C22(m,n,:),[1,M]); c22([1,end])=0.5*c22([1,end]);
        m1 = (m-1)*M+1;m2=m*M;
        n1 = (n-1)*M+1;n2=n*M;
        H11(m1:m2,n1:n2) = R_p(m,m)*c11.*G_p;
        H12(m1:m2,n1:n2) = R_p(m,m)*c12.*G_p;
        H21(m1:m2,n1:n2) = R_n(m,m)*c21.*G_n;
        H22(m1:m2,n1:n2) = R_n(m,m)*c22.*G_n;
    end
end

H = [H11,H12;H21,H22]*dx;
[E1,y]=eigs(H-eye(2*K*M),1,'smallestabs');


