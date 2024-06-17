clc;clear; close all;
addpath('functions');
load('data\Par_gain');
load('data\Par_Passive'); [~,PsLayer]=max(Rt(1,:)); %m = layer number of the passive cavity
w_ab=Par_gain(1);na=Par_gain(2);rp=Par_gain(3); rl=Par_gain(4);
load('data\Pump');

%% Enable only one of the following two initialization methods:

%% Initial guess for Dmax equal to or slightly larger than the seond threshold (must pre-run Dth1_Dth2.m) 

% load('data\Dth2');
% load('data\Dth2_LM');
% load('data\Dth2_wd');
% 
% K = length(Dth2_LM)/2; x =linspace(0,1,K);
% Dth2_w0=Dth2_LM(1);
% Dth2_I0=Dth2_LM(2);
% Dth2_E0=[ones(length(Dth2_I0),1),1i*Dth2_LM(3:2:end)+Dth2_LM(4:2:end)].*sqrt(Dth2_I0);
% Gn=@(x,x0,w)Green(x',x0,w,na,Lft,Rt);
% [y,Dth2_E1]=M1_E1([real(Dth2_wd),imag(Dth2_wd)],Pump(Dth2,x),Dth2_w0,Dth2_E0,w_ab,rp,rl,Gn);%y
% Dth2_E1 = Dth2_E1/Dth2_E1(1)*0.1;
% Dth2_Ep = Dth2_E1(1:K); %max(abs(D1_Ep))
% Dth2_En = conj(Dth2_E1(K+1:end)); %max(abs(D1_En))
% 
% M=3; % Set the total number=(2*M+1) of the comb teeth
% In = zeros(2*M+1,K*2);
% In(M+1,:) = Dth2_LM;
% In(M,1:2:end) = imag(Dth2_En); In(M,2:2:end) = real(Dth2_En);
% In(M+2,1:2:end) = imag(Dth2_Ep); In(M+2,2:2:end) = real(Dth2_Ep);
% In(M+2,1) = real(Dth2_wd);
% 
% s=0.0025; % set the absorption
% Dmax = Dth2; % set the pumping strength

%% Using an existing PALT solution above the second threshold as the initial guess
load('data\COMB_s0.0041D0.2M7.mat'); In=Out; K=size(In,2)/2; 

s=0.0041; % set the absorption
Dmax = 0.21; % set the pumping strength

%% Solving PALT equations for the comb
x=linspace(0,1,K);
D=Pump(Dmax,x);
Rt(3,PsLayer)=s;
Gn=@(x,x0,w)Green(x',x0,w,na,Lft,Rt);

options = optimoptions('fsolve','Display','iter');
Out=fsolve(@(In)M2(In,D,w_ab,rp,rl,Gn),In,options);
[MM,K]=size(Out);M=(MM-1)/2;
save(['data\COMB_s',num2str(Rt(3,PsLayer)),'D',num2str(Dmax),'M',num2str(M),'.mat'],'Out');

w0 = Out(M+1,1); dw = Out(M+2,1); I0 = Out(M+1,2);
Out(M+1,1)=0;Out(M+2,1)=0;Out(M+1,2)=1;
E = 1i*Out(:,1:2:end)+Out(:,2:2:end);
I= abs(E(:,1)).^2*I0;

figure(1)
clf;
semilogy((dw*[-M:M]+w0+w_ab),I,'ro');
ylim([1e-4 1e2]);
yticks([10.^(-4:2)]);
title(['Dmax=',num2str(Dmax),', $$\sigma/\varepsilon_0 =',num2str(Rt(3,PsLayer)),' (n^2c/L)$$'],'interpreter','latex');
xlabel('\omega (c/L)'); ylabel('|E_m(x=0)|^2');

