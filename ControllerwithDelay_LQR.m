% Continuous-time control 
clc
clear 
addpath('functions/');

%% System
% Transfer function 
Ts = 0.02;
A = [0 1 0; 0 0 1; 0 -200 -30];
B = [0; 0; 4000];
C = [1 0 0];
D = 0;
ss_c = ss(A,B,C,D);
ss_d = c2d(ss_c, Ts);
Ad = ss_d.A;
Bd = ss_d.B;
Cd = ss_d.C;
Dd = ss_d.D;

%% Augument State Space System
% reference tracking
Aaug = [Ad Bd;zeros(1,3) 0];
Baug = [0;0;0;1];
%Cbm = [Baug Aaug*Baug Aaug^2*Baug Aaug^3*Baug]  % rank = 4
Caug = [1 0 0 0];
Daug = 0;
ss_daug = ss(Aaug,Baug,Caug,Daug);

%% LQR control
% without delay
clc
R = 0.1;
v = [1000 0.001 0.001];
Q = diag(v);
[Klqr_old,S,P] = dlqr(Ad,Bd,Q,R);

% without delay without redesign
R2 = 0.1;
v2 = [1000 0.001 0.001 1];
Q2 = diag(v2);
[Klqr_old2,S2,P2] = dlqr(Aaug,Baug,Q2,R2);
Klqr_old2
eig(Aaug-Baug*Klqr_old2)

R3 = 100000;
v3 = [.001 0.0001 0.00001 .001];
Q3 = diag(v3);
[Klqr_new,S3,P3] = dlqr(Aaug,Baug,Q3,R3);
Klqr_new
eig(Aaug-Baug*Klqr_new)
Klqr_new=place(Aaug,Baug,[-0.4 -0.5069+0.2i -0.5069-0.2i -0.1])
eig(Aaug-Baug*Klqr_new)