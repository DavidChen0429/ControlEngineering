clear all
%% Try to understand the lead and lag compensator 
s = tf('s');
G = 1/(s^2+2*s+5); 
H = feedback(G, -1);

%% Before adding compensator
figure(1);
subplot(3,1,1)
step(H)
subplot(3,1,2)
rlocus(H)
subplot(3,1,3)
bode(H)

%% Adding a lag compensator
C1 = (s+1)/(s+2);
H1 = feedback(C1*G, -1);
figure(2);
subplot(3,1,1)
step(H1)
subplot(3,1,2)
rlocus(H1)
subplot(3,1,3)
bode(H1)
% The root locus moves right, results in less stability, slow down settling 
% time

%% Add a lead compensator
C2 = (s+2)/(s+1);
H2 = feedback(C2*G, -1);
figure(3);
subplot(3,1,1)
step(H2)
subplot(3,1,2)
rlocus(H2)
subplot(3,1,3)
bode(H2)
% The root locus move left, results in more stability, speed up settling
% time

%% Visualization comparsion
figure(4)
step(H)
hold on
step(H1)
step(H2)
legend('original', 'lag', 'lead')