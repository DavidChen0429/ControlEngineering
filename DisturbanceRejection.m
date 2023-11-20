%% Continuous-time control 
clc
clear 
addpath('functions/');

%% Original system dynamic analysis
syms s;
numerator = 4000;
denominator = [1, 30, 200, 0];
G = tf(numerator, denominator); 

Kp = 0.49;
Ki = 1;
Kd = 0.08;
%Kp = 2.9;
%Ki = 15;
%Kd = .3;
Tf = 0.01;
C = pid(Kp, Ki, Kd, Tf);

DisRejectVisual(G, C)