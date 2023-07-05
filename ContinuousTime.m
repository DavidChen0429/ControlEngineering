%% Continuous-time control 
clc
clear 
addpath('functions/');

%% Original system dynamic analysis
syms s;
numerator = 4000;
denominator = [1, 30, 200, 0];
L = tf(numerator, denominator); 
CL_sys = feedback(L, 1);
PlotSysAll(L);

%% Higher order system estimation
L_approximate = tf(4000, [1, 10, 0]);   % high order approximation
CL_sys_approximate = feedback(L_approximate, 1);

%% PID design for reference tracking
% Design PID controller (transfer function version)
% P controller
Kp = 0.15;
Ki = 0;
Kd = 0;
Tf = 0.01;
P = pid(Kp, Ki, Kd, Tf);

P_ControlledSys = L * P;
P_ControlledCLSys = feedback(P_ControlledSys, 1);
figure();
stepinfo(P_ControlledCLSys)
PlotSysAll(P_ControlledSys)

% PI controller
Kp = 0.15;
Ki = 0.1;
Kd = 0;
Tf = 0.01;
PI = pid(Kp, Ki, Kd, Tf);

PI_ControlledSys = L * PI;
PI_ControlledCLSys = feedback(PI_ControlledSys, 1);
figure();
stepinfo(PI_ControlledCLSys)
PlotSysAll(PI_ControlledSys)

% PD controller
Kp = 0.48;
Ki = 0;
Kd = 0.08;
Tf = 0.01;
PD = pid(Kp, Ki, Kd, Tf);

PD_ControlledSys = L * PD;
PD_ControlledCLSys = feedback(PD_ControlledSys, 1);
figure();
stepinfo(PD_ControlledCLSys)
PlotSysAll(PD_ControlledSys)
