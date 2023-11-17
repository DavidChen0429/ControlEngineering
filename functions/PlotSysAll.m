%% Plotting function for system
function [Gain_Margin, Phase_Margin, Wcg, Wcp] = PlotSysAll(sys)
% This function is defined for visualizing general behavior of a system
% Noted that the given parameter is open-loop transfer function 'L'
% Step respond
subplot(2,2,1);
ol_sys = sys;
cl_sys = feedback(sys, 1);
step(cl_sys);                               

% Bode diagram with margin
subplot(2,2,2);
margin(ol_sys);
grid on;
[Gain_Margin,Phase_Margin, Wcg, Wcp] = margin(sys);

% Nyquist plot 
subplot(2,2,4)
nyquist(sys);
grid on;

% Root locus 
subplot(2,2,3)
rlocus(cl_sys);
grid on; 
end