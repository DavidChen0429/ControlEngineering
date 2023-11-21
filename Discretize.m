% Continuous-time control 
clc
clear 
addpath('functions/');

%% Original system dynamic analysis
syms s;
numerator = 4000;
denominator = [1, 30, 200, 0];
L = tf(numerator, denominator); 

%% Discretizes system
[A,B,C,D] = tf2ss(numerator, denominator);
ss_c = ss(A,B,C,D);
Ts = 1/50;
ss_d = c2d(ss_c, Ts);
Ad = ss_d.A;
Bd = ss_d.B;
Cd = ss_d.C;
Dd = ss_d.D;

%% Discretizes controller

% continuous reference tracking
Kp = 0.49;
Ki = 0;
Kd = 0.08;
Tf = 0.01;
C_crt = pid(Kp, Ki, Kd, Tf);
% continuous disturbance rejection 
Kp = 2.9;
Ki = 15;
Kd = 0.3;
Tf = 0.01;
C_cdr = pid(Kp, Ki, Kd, Tf);

Ts = 0.02;
% discrete reference tracking (bode compare, step compare)
C_drt_zoh = c2d(C_crt, Ts, 'zoh');
C_drt_tustin = c2d(C_crt, Ts, 'tustin');
C_drt_matched = c2d(C_crt, Ts, 'matched');
% discrete disturbance rejection (bode compare, step compare)
C_ddr_zoh = c2d(C_cdr, Ts, 'zoh');
C_ddr_tustin = c2d(C_cdr, Ts, 'tustin');
C_ddr_matched = c2d(C_cdr, Ts, 'matched');