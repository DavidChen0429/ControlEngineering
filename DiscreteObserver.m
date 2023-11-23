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

%% Pole placement
%Cb = [Cd; Cd*Ad; Cd*Ad^2]; % observability
kDisc=place(Ad,Bd,[0.5 0.4+0.3j 0.4-0.3j]);

%% Create observer
KDiscObe=acker(Ad',Cd',[0.2; 0.3+0.2j; 0.3-0.2j])'; %3*1
%KDiscObe=acker(Ad',Cd',[0.1; 0.1+0.1j; 0.1-0.1j])'; %3*1
Aobs = Ad-KDiscObe*Cd;
Bobs = [Bd KDiscObe];
Cobs = eye(3);
Dobs = zeros(3,3);