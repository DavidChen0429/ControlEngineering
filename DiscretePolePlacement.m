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

damping_ratio = 0.45;
theta = (pi-acos(damping_ratio));
mod = -5;
real = mod*cos(theta);
img = mod*sin(theta);
pole1 = real + img*1j;
pole2 = real - img*1j;
pCont1=[-20 -2.9133+27.9300i -2.9133-27.9300i]; %poles in s-domain
pCont2=[-17.8337 -12.3574+34.7369i -12.3574-34.7369i]; %poles in s-domain
pCont3=[-35 -22.2900+44.8028i -22.2900-44.8028i]; %poles in s-domain
pCont4=[-10 -2.9133+27.9300i -2.9133-27.9300i];
pDisc1=exp(pCont1.*Ts); % poles in z-domain
pDisc2=exp(pCont2.*Ts); % poles in z-domain
pDisc3=exp(pCont3.*Ts); % poles in z-domain
pDisc4=exp(pCont4.*Ts); % poles in z-domain
kCont1=place(A,B,pCont1);
kCont2=place(A,B,pCont2);
kCont3=place(A,B,pCont3);
kCont4=place(A,B,pCont4);
kDisc1=place(Ad,Bd,pDisc1);
kDisc2=place(Ad,Bd,pDisc2);
kDisc3=place(Ad,Bd,pDisc3);
kDisc4=place(Ad,Bd,pDisc4);
kDiscF=place(Ad,Bd,[0.01 0.01+0.1j 0.01-0.1j]);
kContF=place(A,B,log([0.01 0.01+0.1j 0.01-0.1j])/Ts);
sysclCont1 = ss(A-B*kCont1,B*kCont1,C,D);
sysclCont2 = ss(A-B*kCont2,B*kCont2,C,D);
sysclCont3 = ss(A-B*kCont3,B*kCont3,C,D);
sysclCont4 = ss(A-B*kCont4,B*kCont4,C,D);
sysclContF = ss(A-B*kContF,B*kContF,C,D);
%{
sysclDisc1 = ss(Ad-Bd*kDisc1,Bd*kDisc1,Cd,Dd);
sysclDisc2 = ss(Ad-Bd*kDisc2,Bd*kDisc2,Cd,Dd);
sysclDisc3 = ss(Ad-Bd*kDisc3,Bd*kDisc3,Cd,Dd);
sysclDisc4 = ss(Ad-Bd*kDisc4,Bd*kDisc4,Cd,Dd);
%}
t = 0:0.01:2;
r = [ones(1,length(t)); 
    zeros(1,length(t)); 
    zeros(1,length(t))];
x0 = zeros(size(sysclCont1.A, 1), 1);
% Continuous
%% Close the loop in code
sysclDiscBuf = ss(Ad-Bd*kDiscF,Bd,Cd,Dd);
Kdc = dcgain(sysclDiscBuf);
sysclDiscF = ss(Ad-Bd*kDiscF,Bd*(1/Kdc),Cd,Dd);
step(sysclDiscF)

%% Visualize the step response
[y1, t, x1] = lsim(sysclCont1, r, t, x0);
[y2, t, x2] = lsim(sysclCont2, r, t, x0);
[y3, t, x3] = lsim(sysclCont3, r, t, x0);
[y4, t, x4] = lsim(sysclCont4, r, t, x0);
[yF, t, xF] = lsim(sysclContF, r, t, x0);
figure;
plot(t, y1);
hold on;
plot(t, y2);
plot(t, y3);
plot(t, y4);
plot(t, yF);
hold off;
legend;
stepinfo(sysclCont1(1))
stepinfo(sysclCont2(1))
stepinfo(sysclCont3(1))
stepinfo(sysclCont4(1))
stepinfo(sysclContF(1))

%{
% Discrete
[y1d, t, x1d] = lsim(sysclDisc1, r, t, x0);
[y2d, t, x2d] = lsim(sysclDisc2, r, t, x0);
[y3d, t, x3d] = lsim(sysclDisc3, r, t, x0);
[y4d, t, x4d] = lsim(sysclDisc4, r, t, x0);
figure;
plot(t, y1d);
hold on;
plot(t, y2d);
plot(t, y3d);
plot(t, y4d);
hold off;
legend;
%}

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