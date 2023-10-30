clear all
%% Try to understand the lead and lag compensator 
s = tf('s');
G = 1/(s^2+2*s+5); 
H = feedback(G, -1);

%% Before adding compensator
figure('Name', 'Original system');
subplot(3,1,1)
step(H)
subplot(3,1,2)
rlocus(H)
subplot(3,1,3)
bode(G)

%% Adding a lag compensator
C1 = (s+1)/(s+2);
H1 = feedback(C1*G, -1);
figure('Name', 'Lag compensator');
subplot(3,1,1)
step(H1)
subplot(3,1,2)
rlocus(H1)
subplot(3,1,3)
bode(C1*G)
% The root locus moves right, results in less stability, slow down settling 
% time

%% Add a lead compensator
C2 = (s+2)/(s+1);
H2 = feedback(C2*G, -1);
figure('Name', 'Lead compensator');
subplot(3,1,1)
step(H2)
subplot(3,1,2)
rlocus(H2)
subplot(3,1,3)
bode(C2*G)
% The root locus move left, results in more stability, speed up settling
% time

%% Add a pole
%{
C3 = (1)/(s+1);
H3 = feedback(C3*G, -1);
figure('Name', 'Add a pole');
subplot(3,1,1)
step(H3)
subplot(3,1,2)
rlocus(H3)
subplot(3,1,3)
bode(C3*G)
%}

%% Visualization comparsion
figure('Name', 'Compare step response')
step(H)
hold on
step(H1)
step(H2)
legend('original', 'lag', 'lead')