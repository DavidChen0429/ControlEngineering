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

% disturbance rejection
Adaug = [Ad zeros(3,1);-Cd 1];
Bdaug = [Bd; 0];
Cdaug = [1 0 0 0];

Aaug2 = [Ad zeros(3,1) Bd; -Cd 1 0;zeros(1,5)];
Baug2 = [zeros(3,1);0;1];
Caug2 = [1 0 0 0 0];

%% Pole placement 
% reference tracking
kDisc_old=place(Ad,Bd,[0.5 0.4+0.5j 0.4-0.5j]);
kDisc_old2=place(Aaug,Baug,[0.5 0.4+0.5j 0.4-0.5j 1.1]);          % old control
kDisc_new=place(Aaug,Baug,[0.5 0.4+0.5j 0.4-0.5j 0.6]);
kDisc_new2=place(Aaug,Baug,[0.5 0.4+0.5j 0.4-0.5j 0.1]);
kDisc_new3=place(Aaug,Baug,[0.2 0.2+0.1j 0.2-0.1j 0.1]);

% dr without time delay
kDisc_dr=place(Adaug,Bdaug,[0.2 0.2+0.1j 0.2-0.1j 0.8]); 

% dr with time delay
kDisc_dr2=place(Aaug2,Baug2,[0.2 0.2+0.1j 0.2-0.1j 0.8 1]); % old control
kDisc_dr_delay=place(Aaug2,Baug2,[0.2 0.2+0.1j 0.2-0.1j 0.8 0.1]); 
kDisc_dr_delay2=place(Aaug2,Baug2,[0.2 0.2+0.1j 0.2-0.1j 0.12 0.1]); 