%% Q5 Pole placement method

b = 4000;
a = [1 30 200 0];

[A1, B1, C1, D1] = tf2ss(b, a);
sys_CT = ss(A1, B1, C1, D1);

%since the biplane has to react relatively fast I start with a sampling
%time of 0.1[s] however after testing better results are achieved with
%0.05[s].
ts = 0.122/10; % sampling time determined by 10 samples in the rise time
sys_DT = c2d(sys_CT, ts, 'zoh'); %zoh good method for unit step, 
% impulse good method for (output) disturbance rejection 'impulse response'

%since all poles are in the RHP we want to shift them or make them obselete
p_con_h = [0.1, 0.5+0.175i, 0.5-0.175i]; %[0.1, 0.75+0.175i, 0.75-0.175i]; slower with large dominant poles
p_con_l = [0.1, 0.09+0.3i, 0.09-0.3i]; % [0.1, 0.75+0.175i, 0.75-0.175i]; faster but a higher gain is required for
% the poles close to the origin
p_rel_h = [0.15, 0.1, 0.6]; %[0.15, 0.1, 0.75]; slower with large dominant pole
p_rel_l = [0.11, 0.1, 0.2]; %[0.11, 0.1, 0.2]; faster but a higher gain is required for
% the poles close to the origin

%Defing the gains for both complex and real pole sets with and without
%dominant poles
K_con_h = place(sys_DT.A,sys_DT.B,p_con_h);
K_con_l = place(sys_DT.A,sys_DT.B,p_con_l);
K_rel_h = place(sys_DT.A,sys_DT.B,p_rel_h);
K_rel_l = place(sys_DT.A,sys_DT.B,p_rel_l);
%[K,prec] = place(sys_DT.A,sys_DT.B,p);
    
Acl_con_h = sys_DT.A-sys_DT.B*K_con_h;
Acl_con_l = sys_DT.A-sys_DT.B*K_con_l;
Acl_rel_h = sys_DT.A-sys_DT.B*K_rel_h;
Acl_rel_l = sys_DT.A-sys_DT.B*K_rel_l;

syscl_con_h = ss(Acl_con_h,sys_DT.B,sys_DT.C,sys_DT.D,ts);
syscl_con_l = ss(Acl_con_l,sys_DT.B,sys_DT.C,sys_DT.D,ts);
syscl_rel_h = ss(Acl_rel_h,sys_DT.B,sys_DT.C,sys_DT.D,ts);
syscl_rel_l = ss(Acl_rel_l,sys_DT.B,sys_DT.C,sys_DT.D,ts);

Kdc_con_h = dcgain(syscl_con_h);
Kdc_con_l = dcgain(syscl_con_l);
Kdc_rel_h = dcgain(syscl_rel_h);
Kdc_rel_l = dcgain(syscl_rel_l);

syscl_scaled_con_h = ss(Acl_con_h,(1/Kdc_con_h)*sys_DT.B,sys_DT.C,sys_DT.D,ts);
syscl_scaled_con_l = ss(Acl_con_l,(1/Kdc_con_l)*sys_DT.B,sys_DT.C,sys_DT.D,ts);
syscl_scaled_rel_h = ss(Acl_rel_h,(1/Kdc_rel_h)*sys_DT.B,sys_DT.C,sys_DT.D,ts);
syscl_scaled_rel_l = ss(Acl_rel_l,(1/Kdc_rel_l)*sys_DT.B,sys_DT.C,sys_DT.D,ts);

figure()
step(syscl_scaled_con_h), grid on, hold on
step(syscl_scaled_rel_h)
step(syscl_scaled_con_l)
step(syscl_scaled_rel_l)
legend('Dominant complex conjugate poles','Real poles with 1 dominant pole','Complex conjugate poles with a fast dominant pole','Real poles with a fast dominant pole')
hold off

stepinfo(syscl_scaled_con_h)
stepinfo(syscl_scaled_con_l)
stepinfo(syscl_scaled_rel_h)
stepinfo(syscl_scaled_rel_l)