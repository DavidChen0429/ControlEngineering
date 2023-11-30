function [y, t] = DisRejectVisual(G ,C)
figure()
% time domain response
subplot(1,2,1);
cl_sys = G/(1+C*G);
t = 0:0.1:6;  % Time vector
[y, t] = step(cl_sys, t);
%step(cl_sys)
stairs(t, y);

% freq domain response with margin
subplot(1,2,2);
margin(G*C);
hold on;
bode(G);
bode(C);
hold off;
legend('G*C', 'G', 'C')
end