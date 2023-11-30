% Continuous-time control 
clc
clear 
addpath('functions/');

%% System
% Transfer function 
Ts = 0.02;
syms s;
numerator = 4000;
denominator = [1, 30, 200, 0];
L = tf(numerator, denominator); 
Ld = c2d(L, Ts);
% State space 
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

%% PID 
Tf = 1/1000;
C_cdr = pid(20, 100, 2, Tf); %p, i, d, tf
C_ddr_tustin = c2d(C_cdr, Ts, 'tustin');
% Disturbance Rejection
C = pid(2.9, 15, .3, 0.01);
L_new = tf(numerator, [1, 30, 200, 0, 0]); 
C2 = pid(0.9, 0.7, 0.3, 0.01);
c2d(C2, Ts, 'tustin')
% Reference Tracking
Crt1 = pid(0.49, 0, 0.08, 0.01);
Crt2 = pid(0.4, 0, 0.05, 0.01);
PD_rt1= L_new * Crt1;
PD_clrt1 = feedback(PD_rt1, 1);
PD_rt2= L_new * Crt2;
PD_clrt2 = feedback(PD_rt2, 1);
stepinfo(PD_clrt1)
stepinfo(PD_clrt2)

[y1, t1] = DisRejectVisual(L, C);
[y3, t3] = DisRejectVisual(L_new, C);
[y2, t2] = DisRejectVisual(L_new, C2);
figure()
title("Disturbance Rejection")
stairs(t1,y1);
hold on
grid on
stairs(t2,y2)
hold off 
grid off
legend('Original','Time Delay Redesign')

