%% Continuous-time control 
clc
clear 
addpath('functions/');

%% 
syms s;
numerator = 4000;
denominator = [1, 30, 200, 0];
L = tf(numerator, denominator); 
L_approximate = tf(4000, [1, 10, 0]);   % high order approximation
PlotSysAll(L);

% Created close loop system
CL_sys = feedback(L, 1);
CL_sys_approximate = feedback(L_approximate, 1);

% % Design PID controller (transfer function version) 
% Kp = 2;
% Ki = 0;
% Kd = 2;
% Tf = 0.01;
% PID_controller = pid(Kp,Ki, Kd, Tf);
% 
% % PID controlled system
% new_open_loop_sys = open_loop_sys * PID_controller;
% Controlled_system = feedback(new_open_loop_sys, 1);
% plot_system_all(Controlled_system)