% Continuous-time control 
clc
clear 
addpath('functions/');

%% Discretizes system
A = [0 1 0; 0 0 1; 0 -200 -30];
B = [0; 0; 4000];
C = [1 0 0];
D = 0;
ss_c = ss(A,B,C,D);
Ts = 1/50;
ss_d = c2d(ss_c, Ts);
Ad = ss_d.A;
Bd = ss_d.B;
Cd = ss_d.C;
Dd = ss_d.D;

%% Pole placement
Cb = [Bd, Ad*Bd, Ad^2*Bd]; % controllable

damping_ratio = 0.69;
theta = (pi-acos(damping_ratio));
mod = -10;
real = mod*cos(theta);
img = mod*sin(theta);
pole1 = real + img*1j;
pole2 = real - img*1j;
pCont=[-5.268 -2.9133+27.9300i -2.9133-27.9300i]; %poles in s-domain
pDisc=exp(pCont.*Ts); % poles in z-domain
kCont=place(A,B,pCont);
kDisc=place(Ad,Bd,pDisc);
sysclCont = ss(A-B*kCont,B,C,D);
step(sysclCont)

%{
kDis = place(Ad,Bd,[0.6703 0.8187 0.9802]);
sys1 = ss(Ad-Bd,Bd,Cd,Dd);
Kdc = dcgain(sys1);
L = 1/Kdc;
sys2 = ss(Ad-Bd*kDis,Bd*L,Cd,Dd*L);
sysclDis = ss(Ad-Bd*kDis, Bd*kDis, Cd, Dd);
t = 0:0.01:10;
r = [ones(1,length(t)); 
    zeros(1,length(t)); 
    zeros(1,length(t))];
x0 = zeros(size(sysclDis.A, 1), 1);
[y, t, x] = lsim(sysclDis, r, t, x0);
plot(t, y);
%}