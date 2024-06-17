clc;clear; close all;
addpath('functions');

load('data\Par_gain');
load('data\Par_Passive'); 
w_ab=Par_gain(1);na=Par_gain(2);rp=Par_gain(3); rl=Par_gain(4);
Gn=@(x,x0,w)Green(x',x0,w,na,Lft,Rt);
x0=linspace(0,1,201); Mx=length(x0);

load('data\Pump');
load('data\Dth1');
load('data\Dth1_w');Dth1_w0=Dth1_w(1);Dth1_w1=Dth1_w(2);%D0_w2=D0_w(4);%D0_w3=D0_w(4);

%% Calculating the lasing mode above the first lasing thrshold
Dmax=linspace(Dth1,0.1,200);
LM=zeros(length(Dmax),2*Mx);

[y0,Dth1_E0]=M0(Dth1_w0,x0,@(x)Pump(Dth1,x),w_ab,rp,Gn);
Dth1_E0=Dth1_E0/Dth1_E0(1);
In=zeros(1,2*Mx); In(1:2:end)=imag(Dth1_E0);In(2:2:end)=real(Dth1_E0);In(1)=real(Dth1_w0);In(2)=0;

for n=1:length(Dmax)
    Dp = Pump(Dmax(n),x0);
    Gn=@(x,x0,w)Green(x',x0,w,na,Lft,Rt);
    LM(n,:)=fsolve(@(In)M1_E0(In,Dp,w_ab,rp,Gn),In);
    In = LM(n,:);
    n
end
 
w0=LM(:,1);I0=LM(:,2);
E0=[ones(length(Dmax),1),1i*LM(:,3:2:end)+LM(:,4:2:end)].*sqrt(I0);

%% Tracking the cavity modes under SIA (SALT)
w1_SALT=zeros(1,length(Dmax));
In=Dth1_w1;w1_SALT(1)=In;
for n=2:length(Dmax)
    D=Pump(Dmax(n),x0);
    Gn=@(x,x0,w)Green(x',x0,w,na,Lft,Rt);
    w1_SALT(n)=fsolve(@(w)M1_E1_SALT(w,D,w0(n),E0(n,:),w_ab,rp,Gn),w1_SALT(n-1));
    disp(['SALT',num2str(n)])
end

%% PALT stability analysis

%==========================eigenvalue 1=========================
w1_PALT=zeros(1,length(Dmax));
In=[real(Dth1_w1)-w0(1),imag(Dth1_w1)];
for n=1:length(w1_PALT)
    D=Pump(Dmax(n),x0);
    Gn=@(x,x0,w)Green(x',x0,w,na,Lft,Rt);
    Out=fsolve(@(In)M1_E1(In,D,w0(n),E0(n,:),w_ab,rp,rl,Gn),In);
    w1_PALT(n)=Out(1)+1i*Out(2)+w0(n);
    In = Out;
    disp(['PALT',num2str(n)])
end

%==========================eigenvalue 2=========================
w2_PALT=zeros(1,length(Dmax));
In=[0.00102,0.0002];
for n=length(w2_PALT):-1:1
    D=Pump(Dmax(n),x0);
    Gn=@(x,x0,w)Green(x',x0,w,na,Lft,Rt);
    Out=fsolve(@(In)M1_E1(In,D,w0(n),E0(n,:),w_ab,rp,rl,Gn),In);
    w2_PALT(n)=Out(1)+1i*Out(2)+w0(n);
    In = Out;
    disp(['PALT',num2str(n)])
end

%% Plotting the eigen values vs Pumping 
figure(1) 
clf;
subplot(2,1,1)
plot(Dmax,w0,'r');hold on;
plot(Dmax,real(w1_SALT)','b--');hold on;
plot(Dmax,real(w1_PALT)','b-');hold on;
plot(Dmax,real(w2_PALT)','color','#F7931E');hold on;
ylabel('Re(\omega)');

subplot(2,1,2)
P0=plot(Dmax,imag(w0),'r-');hold on;
S1=plot(Dmax,imag(w1_SALT)','b--');hold on;
P1=plot(Dmax,imag(w1_PALT)','b-');hold on;
P2=plot(Dmax,imag(w2_PALT)','color','#F7931E');hold on;
ylabel('Im(\omega)');
xlabel('Dmax');

legend([P0,P1,P2,S1],'Lasing mode','PALT','PALT','SALT');
%% Plotting the lasing-field intensity vs Pumping
figure(2)
clf;
plot(Dmax,I0,'r');hold on;
set(gca,'XColor','k','fontsize',14);

%% Determining the second threshold
Dth2_n=find(imag(w2_PALT)<0, 1 ,'last')+1;
Dth2 = Dmax(Dth2_n);
Dth2_LM = LM(Dth2_n,:);
Dth2_wd = w2_PALT(Dth2_n)-w0(Dth2_n);

% save the results for Dth2_Dth3.m
save('data\Dth2','Dth2');
save('data\Dth2_LM','Dth2_LM');
save('data\Dth2_wd','Dth2_wd');