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
Aaug = [Ad Bd;zeros(1,3) 0];
Baug = [0;0;0;1];
%Cbm = [Baug Aaug*Baug Aaug^2*Baug Aaug^3*Baug]  % rank = 4
Caug = [1 0 0 0];