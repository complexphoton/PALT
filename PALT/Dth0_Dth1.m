clc;clear; close all;
addpath('functions');
%% Physical Parameters of the laser system
%=============active cavity==============
lambda=0.82; %(um) central lasing wavelength
L=2.05;%(um) active cavity length
w_ab=2*pi*L/lambda; % gain center
rp=(10*L/3e2); %normalized dephasing rate, aka, gain bandwidth 
rl=rp*1e-4; %normalized relaxation rate of population inversion 
na=3.4; %refractive index
%==============passive cavity============
np=3.67; %refractive index
s=0.0025; %material loss
L2=1.34/L; %passive cavity length
%===============DBR======================
n_1 = 3.0; %index of material A
n_2 = 1.5; %index of material B
DBR_pl = 3; %number of DBR1 periods  
DBR_pm = 6; %number of DBR3 periods  
DBR_pr = 3; %number of DBR2 periods 
dm1 = 0.08/L;dm2=0.13/L; %layer thickness of DBR3
dl1=0.07/L;dl2=0.14/L; %layer thickness of DBR1
dr1=0.07/L;dr2=0.14/L;%layer thickness of DBR2
%=============Pumping===============
Pump=@(Dmax,x) 0.5*Dmax*(1-cos(2*pi*x)); %Pumping profile Dp(x)
Dmax=0.0124; %The peak of Pumping strength
%===========DBR setup====================
DBR_e = n_1^2;
DBR_er1 = real(DBR_e);
DBR_er2= n_2^2;

d_left = repmat([dl1 dl2],1,DBR_pl);
d_mid  = [dm2,repmat([dm1 dm2],1,DBR_pm)];
d_right  = repmat([dr2 dr1],1,DBR_pr);
er_left = repmat([DBR_er1 DBR_er2],1,DBR_pl);
er_mid  = [DBR_er2,repmat([DBR_er1 DBR_er2],1,DBR_pm)];
er_right  = repmat([DBR_er2 DBR_er1],1,DBR_pr);

sigma = w_ab*imag(DBR_e)/real(DBR_e);
s_left = repmat([sigma 0],1,DBR_pl);
s_mid  = [0,repmat([sigma 0],1,DBR_pm)];
s_right  = repmat([0 sigma],1,DBR_pr);


d =[d_left,1,d_mid,L2,d_right];
er=[er_left,na^2,er_mid,np^2,er_right];
si=[s_left,0,s_mid,s,s_right];

%======================Display the setup==============================
X = cumsum(d);
x=zeros(1,2*length(d));
dmin = -0.5; dmax = X(end)+0.5;%dmax = 3.5;
% % 
erx=zeros(1,2*length(d));
x(1:2:end)=[0,X(1:end-1)];x(2:2:end)=X(1:end);
erx(1:2:end)=er;erx(2:2:end)=er;

figure(1)
clf;
subplot(2,1,1)
plot([dmin, 0, x, x(end), dmax],sqrt([1,1,erx,1,1]),'k');
xlabel('x/L'),ylabel('index')
axis([dmin dmax 0 4]);

sx=zeros(1,2*length(d));
sx(1:2:end)=si;sx(2:2:end)=si;
% % % 
% % % 
subplot(2,1,2)
yyaxis left
x00=linspace(0,1,1000);
area([x00+X(2*DBR_pl)],[Pump(Dmax,x00)],'basevalue',0,'LineStyle','none','FaceColor','#F7931E','EdgeColor','none');ylim([0 0.25]);hold on;
ylim([0 0.015]); 
ylabel('D_p'); 
yyaxis right
nl=2*(DBR_pl+DBR_pm)+2;Nr=nl+1;
area([X(nl),X(Nr)] ,[s,s]*3e2/L*np^2,'basevalue',0,'LineStyle','none','FaceColor','#0071BC','EdgeColor','none');
ylim([0 8]);
xlim([dmin dmax]);
ylabel('absorption'); 

%% Solving eigen frequencies and determining the first lasing threshold

% ===================Solver for Homogeneous Pumping======================
% DO NOT use this solver if the Pump function is x-dependent
% Dc=zeros(1,length(d));
% Dc(length(d_left)+1)=Dmax/2;
% Dc(length(d_left)+1)=Dmax/2;

%Initial guess
% w0 = [-0.0175,-0.0189-0.00054i];  
%To get a good initial guess of eigen values, scan the function value y=F0_eigf(...) on the complex frequency
%plane and look for minimum function values of M0.m 

% w0=fsolve(@(w)F0_eigf(w+w_ab,er,si,d,Dc,w_ab,rp),w0);
% 
% figure(2)
% clf;
% plot(real(w0),imag(w0),'ro');
%=======================================================================

%====================Solver for arbitrary NON-ZERO pumping==============
% For zero pumping, use the above solver instead 

Lft = [d_left;sqrt(er_left);s_left];
Rt = [[d_mid,L2,d_right];sqrt([er_mid,np^2,er_right]);[s_mid,s,s_right]];
Par_gain=[w_ab,na,rp,rl];

w_ab=Par_gain(1);na=Par_gain(2);rp=Par_gain(3);
Gn=@(x,x0,w)Green(x',x0,w,na,Lft,Rt);
x0=linspace(0,1,201); Mx=length(x0);

%----------------initial guess-----------------
%To get a good initial guess of eigen values, scan the function value of
% y=M0(w,...) on the complex frequency plane.
% figure(2)
% clf;
% Nr=50;Ni=50;
% re_w = linspace(-0.02,-0.015,Nr); im_w=linspace(-1e-3,1e-3,Ni);
% f = zeros(Ni,Nr);
% for nr=1:Nr
%     for ni=1:Ni
%         f(ni,nr) = M0(re_w(nr)+1i*im_w(ni),x0,@(x)Pump(Dmax,x),w_ab,rp,Gn);
%     end
%     nr
% end
% y = min(abs(f),[],'all')./abs(f);
% surf(re_w,im_w,log10(y));hold on;shading interp; axis tight;view(2)

w0 = [-0.0175,-0.0189-0.00054i]; %initial guess 
%----------------------------------------------

% Dth1_Solver solves for the first threshold, Dth1, using Dmax as the initial guess. Dth1 is
% required as an input for the other scripts. But if you just want the eigenvalues below Dth1,
% simply disable line 136 and line 151-155.
Dmax=fsolve(@(Dmax)Dth1_Solver(w0(1),x0,Dmax,w_ab,rp,Gn,Pump),Dmax); 

Dth1=Dmax;
% Determing the eigen frequencies at the first threshold 
Dth1_w=w0;
[y1,Dth1_w(1)]=Dth1_Solver(w0(1),x0,Dth1,w_ab,rp,Gn,Pump);
[y2,Dth1_w(2)]=Dth1_Solver(w0(2),x0,Dth1,w_ab,rp,Gn,Pump);

figure(2)
plot3(real(Dth1_w),imag(Dth1_w),[1 1],'ro');view(2);
xlabel('Re(\omega)-\omega_a_b (c/L)');
ylabel('Im(\omega) (c/L)');
%=======================================================================

%===============Save the results for Dth1_Dth2.mat======================
save('data\Par_gain','Par_gain');
save('data\Par_Passive','Lft','Rt');
save('data\Dth1','Dth1');
save('data\Dth1_w','Dth1_w');
save('data\Pump','Pump');

%%
function [y,w]=Dth1_Solver(w0,x0,Dmax,w_ab,rp,Gn,Pump)
         Dp=@(x) Pump(Dmax,x);  
         options = optimset('TolFun',1e-9,'TolX',1e-9,'display','on');
         w=fsolve(@(w)M0(w,x0,Dp,w_ab,rp,Gn),w0,options);
         y=imag(w);
end

