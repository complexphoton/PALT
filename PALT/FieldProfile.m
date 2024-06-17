clc;clear; close all;
addpath('functions');
load('data\Par_gain');
load('data\Par_Passive'); [~,PsLayer]=max(Rt(1,:)); %m = layer number of the passive cavity
w_ab=Par_gain(1);na=Par_gain(2);rp=Par_gain(3); rl=Par_gain(4);
load('data\Pump');

Dmax=0.2;s=0.0041;
load(['data\COMB_s',num2str(s),'D',num2str(Dmax),'M7.mat']);

[MM,K]=size(Out);M=(MM-1)/2;K=K/2;
w0 = Out(M+1,1); wr = Out(M+2,1); I0 = Out(M+1,2);
Rt(3,PsLayer)=s;Gn=@(x,x0,w)Green(x',x0,w,na,Lft,Rt);
D0x=Pump(Dmax,linspace(0,1,K));
[y,Ex0,Px0,Dx0]=M2(Out,D0x,w_ab,rp,rl,Gn);
disp(norm(y)); %The imported comb solution is valid iff norm(y) is approximately 0.
Ex0=Ex0.*sqrt(I0);Px0=Px0.*sqrt(I0);

x=linspace(-sum(Lft(1,:))-0.1,sum(Rt(1,:))+1.1,2000); %set the x coordinates to plot the field
Gn_full=@(x,x0,w)Green_full(x',x0,w,na,Lft,Rt);
[E,P,D,B]=M2_FieldPrf(x,@(x)Pump(Dmax,x),Px0,w0,wr,w_ab,rp,rl,Gn_full);

figure(1)
clf;
subplot(3,1,1)
plot(x,abs(E(M,:)).^2);
ylabel('$$|E_{-1}|^2$$','interpreter','latex');
subplot(3,1,2)
plot(x,abs(E(M+1,:)).^2);
ylabel('$$|E_{0}|^2$$','interpreter','latex');
subplot(3,1,3)
plot(x,abs(E(M+2,:)).^2);
ylabel('$$|E_{1}|^2$$','interpreter','latex');
xlabel('x/L');