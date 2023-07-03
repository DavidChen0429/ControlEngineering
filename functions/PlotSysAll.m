%% Plotting function for system
function [Gain_Margin, Phase_Margin, Wcg, Wcp] = PlotSysAll(sys)
% This function is defined for visualizing general behavior of a system
% Noted that the given parameter is open-loop transfer function 'L'
% Step respond
subplot(2,2,1);
step(feedback(sys, 1));                               

% Bode diagram with margin
subplot(2,2,2);
margin(sys);
grid on;
[Gain_Margin,Phase_Margin, Wcg, Wcp] = margin(sys);

% Nyquist plot 
subplot(2,2,4)
nyquist(sys);
grid on;

% Root locus 
subplot(2,2,3)
rlocus(sys);
grid on; 
end