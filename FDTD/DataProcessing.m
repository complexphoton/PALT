clc;clear;close all;
set(0, 'DefaultLineLineWidth', 1);
set(groot,'defaultAxesFontName','Times New Roman')
set(groot,'defaultAxesFontSize',16)
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

%% Reading the time-domain data
load(['E.mat']);dt=data(end);Et=data(1:end-1); N0=length(Et);

% Sampling the data points to accelerate FFT (But be aware of Nyquist Theorem!)
Sampling = 10;
Et=Et(1:Sampling:end);dt=Sampling*dt;
t = (0:(length(Et)-1))'*dt;

%Plot the field
figure(1)
clf;
plot(t/3e5,Et,'color',[1 1 1]*0,'linewidth',1);hold on;
% ylim([-20 20]);yticks([-10:10:10]);
h=gca;
set(gca,'XColor','k','fontsize',16);
h.TickLength = [0.02, 0.01];
xlabel('time t (ps)');
ylabel('E(t)');
%% FFT
% Cut out the transient data before the limit cycle
t_cut=50*3e5;n_cut=ceil(t_cut/dt);
Et=Et(n_cut:end);
t = (0:(length(Et)-1))'*dt; Tw=t(end);

%Padding the time-domain data to improve the spectral resolution
Padding=4*length(Et);
t=dt*(0:Padding-1); 
T=t(end);

%Envelope the time-domain data with Hann window (Hann smoothing)
w=2*pi/T*[0:(length(t)-1)];
Ew=2*dt*fft(Et.*hann(length(Et)),Padding)/Tw;
I_sim = abs(Ew).^2;

%=====================Plots=====================
figure(2)
clf;
P_FDTD=semilogy(w*3e5/2/pi,I_sim,'k','linewidth',0.50);hold on;

%Adding the PALT solution to the plot
L=2.05;%um, Active cavity thickness
f = 3e2/L/2/pi; % Normalizing factor of frequency
load('PALT');
P_PALT=scatter(w*f,abs(E).^2,'filled');
xlim([365.25 365.6]);
ylim([1e-4 1e2]); 
ylabel('$$|E_{\omega}|^2$$','interpreter','latex');
xlabel('Frequency $$f$$ (THz)','interpreter','latex');
legend([P_FDTD, P_PALT],'FDTD','PALT');