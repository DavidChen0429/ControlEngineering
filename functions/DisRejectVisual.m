function [] = DisRejectVisual(G ,C)
% time domain response
subplot(1,2,1);
cl_sys = G/(1+C*G);
step(cl_sys)

% freq domain response with margin
subplot(1,2,2);
margin(G*C);
hold on;
bode(G);
bode(C);
hold off;
legend('G*C', 'G', 'C')
end