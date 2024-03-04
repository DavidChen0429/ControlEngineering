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

%% LQR control
R1 = 1;
v1 = [1 1 1];
Q1 = diag(v1);
[K1,S1,P1] = dlqr(Ad,Bd,Q1,R1);

R2 = 1;
v2 = [1000 1 1];
Q2 = diag(v2);
[K2,S2,P2] = dlqr(Ad,Bd,Q2,R2);

R3 = 1;
v3 = [1000 0.01 0.01];
Q3 = diag(v3);
[K3,S3,P3] = dlqr(Ad,Bd,Q3,R3);

R4 = 1;
v4 = [1000 0.001 0.001];
Q4 = diag(v4);
[K4,S4,P4] = dlqr(Ad,Bd,Q4,R4);

R5 = 0.1;
v5 = [1000 0.001 0.001];
Q5 = diag(v5);
[K5,S5,P5] = dlqr(Ad,Bd,Q5,R5);

R6 = 10;
v6 = [1000 0.001 0.001];
Q6 = diag(v6);
[K6,S6,P6] = dlqr(Ad,Bd,Q6,R6);
