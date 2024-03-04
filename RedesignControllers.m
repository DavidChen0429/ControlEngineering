% Continuous-time control 
clc
clear 
addpath('functions/');

%% System
% Transfer function 
syms s;
numerator = 4000;
denominator = [1, 30, 200, 0];
L = tf(numerator, denominator); 

% State space 
A = [0 1 0; 0 0 1; 0 -200 -30];
B = [0; 0; 4000];
C = [1 0 0];
D = 0;
ss_c = ss(A,B,C,D);
Ts = 0.02;
ss_d = c2d(ss_c, Ts);
Ad = ss_d.A;
Bd = ss_d.B;
Cd = ss_d.C;
Dd = ss_d.D;

%% Original Controllers
% ========================== Discrete PID
Tf = 0.01;
C_crt = pid(0.49, 0, 0.08, Tf); %p, i, d, tf
C_cdr = pid(2.9, 15, 0.3, Tf); %p, i, d, tf
C_drt_tustin = c2d(C_crt, Ts, 'matched');
C_ddr_tustin = c2d(C_cdr, Ts, 'matched');
%{
PD_ControlledSys = L * C_crt;
PD_ControlledCLSys = feedback(PD_ControlledSys, 1);
figure();
PlotSysAll(PD_ControlledSys)
%}

% ========================== Pole Placement 
kDiscF=place(Ad,Bd,[0.01 0.01+0.1j 0.01-0.1j]);
p_cont = log([0.01 0.01+0.1j 0.01-0.1j])*50;
kContF=place(A,B,p_cont);

% ========================== LQRoriginal 
R = 0.1;
Q = diag([1000 0.001 0.001]);
[Klqr,S,P] = dlqr(Ad,Bd,Q,R);

%% Redesigned Controllers
% ============================ PID
C_crt2 = pid(0.14, 0, 0.015, Tf); %p, i, d, tf
C_cdr2 = pid(2.9, 15, 0.3, Tf); %p, i, d, tf
C_drt_tustin2 = c2d(C_crt2, Ts, 'matched');
C_ddr_tustin2 = c2d(C_cdr2, Ts, 'matched');
%{
PD_ControlledSys = L * C_crt2;
PD_ControlledCLSys = feedback(PD_ControlledSys, 1);
figure();
PlotSysAll(PD_ControlledSys)
%}

% ============================ Pole Placement
kDiscF_redesign=place(Ad,Bd,[0.01 0.9+0.085j 0.9-0.085j]);
p_cont = log([0.01 0.9+0.085j 0.9-0.085j])*50;
kContF_redesign=place(A,B,p_cont);

% ============================ LQR
R2 = 750;
Q2 = diag([1000 0.001 0.001]);
[Klqr2,S2,P2] = dlqr(Ad,Bd,Q2,R2);