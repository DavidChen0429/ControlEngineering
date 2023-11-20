% Continuous-time control 
clc
clear 
addpath('functions/');

%% Original system dynamic analysis
syms s;
numerator = 4000;
denominator = [1, 30, 200, 0];
L = tf(numerator, denominator); 

%% Discritize
[A,B,C,D] = tf2ss(numerator, denominator);
ss_c = ss(A,B,C,D);
Ts = 1/50;
ss_d = c2d(ss_c, Ts);