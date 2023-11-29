%% Setup
clc; clear;
set(groot,'defaultAxesXgrid','on')
set(groot,'defaultAxesYgrid','on')

s = tf([1 0],[1]);
G = 5/(s*(0.5*s+1)*(0.2*s+1));
G = zpk(G) %just to check
%figure()
%step(G/(1+G))
%% Question 1
q1.Kc = 1;
q1.Kp = 1;
q1.Ki = 0;
q1.Tau1 = 0.02; %0.02
q1.Tau2 = 0.01; %0.01
q1.Kd1 = (-q1.Kp*q1.Tau1^2+0.7*q1.Tau1-0.1)/(q1.Tau1-q1.Tau2);
q1.Kd2 = 0.7-q1.Kp*q1.Tau1-q1.Kp*q1.Tau2-q1.Kd1;
q1.C = PIDD(q1.Kc,q1.Kp,q1.Ki,q1.Kd1,q1.Kd2,q1.Tau1,q1.Tau2);    
%
q1.L = G*q1.C;
figure()
step(q1.L/(1+q1.L),2)

q1.stepinf = stepinfo(q1.L/(1+q1.L),"SettlingTimeThreshold",0.01);
q1.stepinf.Overshoot
q1.stepinf.SettlingTime
%
figure()
rlocus(q1.L);
%
q1.Kc = 3.41; 
q1.Cc = PIDD(q1.Kc,q1.Kp,q1.Ki,q1.Kd1,q1.Kd2,q1.Tau1,q1.Tau2);
q1.Lc=G*q1.Cc;
figure()
step(q1.Lc/(1+q1.Lc));
grid on;
q1.stepinf = stepinfo(q1.Lc/(1+q1.Lc),"SettlingTimeThreshold",0.01);
fprintf("Closed loop system results in an overshoot of %f \n and a settling time of %f\n",q1.stepinf.Overshoot,q1.stepinf.SettlingTime)
%
figure()
step(q1.Cc/(1+q1.Lc));
ylabel("Voltage")
grid on
%% Analysis q1
figure();
margin(q1.Lc);
grid on;
%% Question 2
%disturbance step with same contoller
%figure();
%step(G/(1+q1.Cc*G)); 

q2.Tau1 = q1.Tau1;
q2.Tau2 = q1.Tau2;
q2.Ki=1;
q2.Kp=0.7-q2.Tau1*q2.Ki-q2.Tau2*q2.Ki;
q2.Kd2 = (0.1+q2.Tau1*q2.Kp-q2.Tau1*q2.Tau2*q2.Ki-q2.Tau2*q2.Kp-q2.Tau1*q2.Kp)/(1-q2.Tau1/q2.Tau2);
%q2.Kd2=-q2.Tau1*q2.Kp-q2.Tau2/q2.Tau1*(0.7-q2.Tau1*q2.Tau2*q2.Ki-q2.Tau1*q2.Kp-q2.Tau2*q2.Kp)/(1-q2.Tau2/q2.Tau1)
q2.Kd1=-(q2.Tau1/q2.Tau2)*q2.Kd2-q2.Tau1*q2.Kp;
%q2.Kd1=0.7-q2.Kd2-q2.Tau1*q2.Tau2*q2.Ki-q2.Tau1*q2.Kp-q2.Tau2*q2.Kp;
q2.Kd1=(0.1-q2.Kp*q2.Tau1-q2.Kp*q2.Tau2-q2.Ki*q2.Tau1*q2.Tau2+q2.Kp*q2.Tau2)/(1-q2.Tau2/q2.Tau1);
q2.Kd2=0.1-q2.Kp*q2.Tau1-q2.Kp*q2.Tau2-q2.Ki*q2.Tau1*q2.Tau2-q2.Kd1;
q2.C=PIDD(1,q2.Kp,q2.Ki,q2.Kd1,q2.Kd2,q2.Tau1,q2.Tau2)

q2.C=PIDD(1,q2.Kp,q2.Ki,q2.Kd1,q2.Kd2,q2.Tau1,q2.Tau2)

minreal(G*q2.C)
step(minreal(G/(1+G*q2.C)))
%%
figure
hold on
for i=1:5
    q2.C=PIDD(1,1,i,1,0,0.001,0)
    step(G/(1+G*q2.C))
end
legend("Ki=1","Ki=2","Ki=3","Ki=4","Ki=5")

figure
hold on
for i=1:1:5
    q2.C=PIDD(1,1,1,i,0,0.001,0)
    step(G/(1+G*q2.C))
end
legend("Kd=1","Kd=2","Kd=3","Kd=4","Kd=5")

figure
hold on
for i=2:2:10
    q2.C=PIDD(i,1,1,1,0,0.001,0)
    step(G/(1+G*q2.C))
end
legend("Kc=2","Kc=4","Kc=6","Kc=8","Kc=10")

figure
margin(minreal(G*q2.C))
%%
figure()
hold on
q2.C=PIDD(100,1,1,1,0,0.01,0)
step(G/(1+G*q2.C))
q2.C=PIDD(100,1,10,1,0,0.01,0)
step(G/(1+G*q2.C))
q2.C=PIDD(10,1,1,1,0,0.01,0)
step(G/(1+G*q2.C))
q2.C=PIDD(10,1,10,1,0,0.01,0)
step(G/(1+G*q2.C))
%
figure()
hold on
q2.C=PIDD(100,1,1,1,0,0.01,0)
step(1/(1+G*q2.C))
q2.C=PIDD(100,1,10,1,0,0.01,0)
step(1/(1+G*q2.C))
q2.C=PIDD(10,1,1,1,0,0.01,0)
step(1/(1+G*q2.C))
q2.C=PIDD(10,1,10,1,0,0.01,0)
step(1/(1+G*q2.C))
%% PI
q2.C=PIDD(1,1,1,0,0,0,0)
figure()
margin(G*q2.C)
grid on
figure()
step(G/(1+G*q2.C))
grid on
%% PID
q2.C=PIDD(1,1,1,1,0,0.01,0)
figure()
margin(G*q2.C)
grid on
figure()
step(G/(1+G*q2.C))
grid on
stepinfo(G/(1+G*q2.C),"SettlingTimeThreshold",0.1)
%% collective gain
figure()
hold on
%q2.C=PIDD(100,1,2,0.5,0,0.001,0);
q2.C=PIDD(10,1,1,1,0,0.01,0);
step(G/(1+G*q2.C))
grid on
figure()
margin(G*q2.C)
grid on
stepinfo(G/(1+G*q2.C),"SettlingTimeThreshold",0.1)
%% final
figure()
hold on
q2.C=PIDD(1,10,20,10,0,0.01,0);
step(G/(1+G*q2.C))
grid on
figure()
margin(G*q2.C)
grid on
stepinfo(G/(1+G*q2.C),"SettlingTimeThreshold",0.1)
max(step(G/(1+G*q2.C)))
% for simulink verification
[num den]=tfdata(q2.C,'v');
[Gnum Gden]=tfdata(G,'v');
%% Question 3
q3.A = [0 1 0; 0 0 1; 0 -10 -7];
q3.B = [0; 0; 50];
q3.C = [1 0 0];
q3.D = [];
q3.sysc = ss(q3.A,q3.B,q3.C,q3.D);
% discretize
q3.sysd = c2d(q3.sysc,0.008);
%% Question 4 set point
% load discrete system
q4.sysd=q3.sysd;
% discretize controller
q4.Con1zoh=c2d(q1.Cc,0.008);
q4.Con1tust=c2d(q1.Cc,0.008,"Tustin");
% create controlled system
q4.Lzoh1=q4.Con1zoh*q4.sysd;
q4.Ltust1=q4.Con1tust*q4.sysd;
figure()
hold on;
step(q4.Lzoh1/(1+q4.Lzoh1));
stepinfo(q4.Lzoh1/(1+q4.Lzoh1),"SettlingTimeThreshold",0.01)
step(q4.Ltust1/(1+q4.Ltust1));
stepinfo(q4.Ltust1/(1+q4.Ltust1),"SettlingTimeThreshold",0.01)
step(q1.Lc/(1+q1.Lc));
grid on
legend("DT - ZOH","DT - Tustin","CT")

%% Question 4 disturbance
q4.h=0.001;
q4.sysd2=c2d(q3.sysc,q4.h)
% discretize controller
q4.Con2zoh=c2d(q2.C,q4.h);
q4.Con2tust=c2d(q2.C,q4.h,"Tustin");
%create controlled system
q4.Lzoh2=q4.Con2zoh*q4.sysd2;
q4.Ltust2=q4.Con2tust*q4.sysd2;
figure()
hold on;
step(q4.sysd2/(1+q4.Lzoh2));
stepinfo(q4.sysd2/(1+q4.Lzoh2),"SettlingTimeThreshold",0.01)
step(q4.sysd2/(1+q4.Ltust2));
stepinfo(q4.sysd2/(1+q4.Ltust2),"SettlingTimeThreshold",0.01)
step(G/(1+q2.C*G));
grid on
legend("DT - ZOH","DT - Tustin","CT")
%%
q4.h=0.008;
q4.sysd2=c2d(q3.sysc,q4.h)
% discretize controller
q2.C4=PIDD(1,1,1,1,0,0.01,0)
q4.Con2zoh=c2d(q2.C4,q4.h);
q4.Con2tust=c2d(q2.C4,q4.h,"Tustin");
%create controlled system
q4.Lzoh2=q4.Con2zoh*q4.sysd2;
q4.Ltust2=q4.Con2tust*q4.sysd2;
figure()
hold on;
step(q4.sysd2/(1+q4.Lzoh2));
stepinfo(q4.sysd2/(1+q4.Lzoh2),"SettlingTimeThreshold",0.01)
step(q4.sysd2/(1+q4.Ltust2));
stepinfo(q4.sysd2/(1+q4.Ltust2),"SettlingTimeThreshold",0.01)
step(G/(1+q2.C4*G));
grid on
legend("DT - ZOH","DT - Tustin","CT")
%% Question 5
q5.h = 0.008;
pole(minreal(q1.L))
exp(pole(minreal(q1.L))*q5.h)
q5.sysd=c2d(q3.sysc,q5.h);
%
q5.zeta=abs(log(0.05))/sqrt(pi^2+log(0.05)^2);
figure()
hold on
for i=25:2:35
    q5.w=i;
    q5.G=1/(s^2+2*q5.w*q5.zeta*s+q5.w^2)
    q5.G1=1/((s^2+2*q5.w*q5.zeta*s+q5.w^2)*(0.5*s+1));
    q5.G2=1/((s^2+2*q5.w*q5.zeta*s+q5.w^2)*(0.2*s+1));
    q5.G3=1/((s^2+2*q5.w*q5.zeta*s+q5.w^2)*(s));
    step(1/dcgain(feedback(q5.G,1))*feedback(q5.G,1));
%     step(1/dcgain(feedback(q5.G1,1))*feedback(q5.G1,1));
%     step(1/dcgain(feedback(q5.G2,1))*feedback(q5.G2,1));
%     step(1/dcgain(feedback(q5.G3,1))*feedback(q5.G3,1));
end
legend("5","10","15","20","25","30");
%%
figure
hold on
q5.w=20;
q5.G=1/(s^2+2*q5.w*q5.zeta*s+q5.w^2);
step(1/dcgain(feedback(q5.G,1))*feedback(q5.G,1))
q5.w=25;
q5.G=1/(s^2+2*q5.w*q5.zeta*s+q5.w^2);
step(1/dcgain(feedback(q5.G,1))*feedback(q5.G,1))
q5.w=30;
q5.G=1/(s^2+2*q5.w*q5.zeta*s+q5.w^2);
step(1/dcgain(feedback(q5.G,1))*feedback(q5.G,1))
step(q1.Lc/(1+q1.Lc))
grid on
legend("w=20","w=25","w=30","PDD reference")
%%
q5.w = 25;
q5.G=1/(s^2+2*q5.w*q5.zeta*s+q5.w^2);
% polynoom z^2+pol1z+pol2
q5.pol1=-2*exp(-q5.zeta*q5.w*q5.h)*cos(q5.w*q5.h*sqrt(1-q5.zeta^2));
q5.pol2=exp(-2*q5.zeta*q5.w*q5.h);
q5.p1=-0.5*q5.pol1+0.5*sqrt(q5.pol1^2-4*q5.pol2);
q5.p2=-0.5*q5.pol1-0.5*sqrt(q5.pol1^2-4*q5.pol2);
q5.p=[q5.p1,q5.p2,min(pole(q5.sysd))]
%q5.p=[q5.p1,q5.p2,0.75]
q5.p=[q5.p1,q5.p2,-0.75]
q5.K = place(q5.sysd.A,q5.sysd.B,q5.p);
q5.sysdfb=q5.sysd;
q5.sysdfb.A = q5.sysdfb.A-q5.sysd.B*q5.K; 
figure
hold on
step(1/dcgain(q5.sysdfb)*q5.sysdfb)
step(1/dcgain(feedback(q5.G,1))*feedback(q5.G,1))
step(q1.Lc/(1+q1.Lc))
legend("DT - Poleplacement","CT - 2^{nd} order reference","CT- Reference Q1","Location","southeast")
stepinfo(1/dcgain(q5.sysdfb)*q5.sysdfb)
grid on
%%
figure()
hold on
pzmap(q5.sysd)
pzmap(q5.sysdfb)
xlim([-1.1 1.1])
ylim([-1.1 1.1])
legend("Original Poles/ Zeros","New Poles/Peros")
%% Question 6 Observer
% poles of observer faster than controller
q6.h=q5.h;
q6.sysd=c2d(q3.sysc,q6.h);
q6.pC = q5.p;
q6.Kc = place(q6.sysd.A,q6.sysd.B,q6.pC);%Controller gain
q6.pO = q5.p/5;%[]; %Observer poles
q6.Lo = place(q6.sysd.A',q6.sysd.C',q6.pC).'; %Observer gain
q6.fwd = 1/dcgain(ss(q6.sysd.A-q6.sysd.B*q6.Kc,q6.sysd.B,q6.sysd.C,q6.sysd.D,q6.h));

% with disturbance estimator
q6.Ldis=200; %integrator gain
q6.Acl = [q6.sysd.A, -q6.sysd.B*q6.Kc -q6.sysd.B; 
    q6.Lo*q6.sysd.C, q6.sysd.A-q6.Lo*q6.sysd.C-q6.sysd.B*q6.Kc zeros(3,1); 
    q6.Ldis*q6.sysd.C -q6.Ldis*q6.sysd.C 1];
q6.Bcl = [q6.sysd.B q6.sysd.B*q6.fwd; zeros(size(q6.sysd.B)) q6.sysd.B*q6.fwd; 0 0];
q6.Ccl =  [1 0 0 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 0 0 1];
syscl=ss(q6.Acl,q6.Bcl,q6.Ccl,[],q6.h);

q6.T = 0:q6.h:1;
q6.x0=[0.5 0.5 0.5 0 0 0 0]';
q6.d = zeros(1,length(q6.T));
q6.r = ones(1,length(q6.T));
q6.u=[q6.d; q6.r];

%figue
%lsim(syscl,q6.u,q6.T,q6.x0)
[Yobs,Tobs,Xobs] = lsim(syscl,q6.u,q6.T,q6.x0);
[Ypp,Tpp,Xpp] = lsim(1/dcgain(q5.sysdfb)*q5.sysdfb,q6.u(2,:),q6.T,q6.x0(1:3)*dcgain(q5.sysdfb));
%
figure()
hold on
stairs(Tobs,Yobs(:,1))
stairs(Tpp,Ypp)
xlabel("Time (s)")
ylabel("Output")
legend("Observer","Pole placement")
figure()
hold on
stairs(Tobs,Yobs(:,1)-Yobs(:,2))
stairs(Tobs,Ypp(:,1)-Yobs(:,1))
ylabel("Error")
xlabel("Time (s)")
legend("$y-\hat{y}$","Poleplacement-Observer","Interpreter","latex")
%% disturbance rejection
q6.T = 0:q6.h:2;
q6.x0=[0 0 0 0 0 0 0]';
q6.d = ones(1,length(q6.T))*10;
q6.r = ones(1,length(q6.T));
q6.d(1:50) = 0;
q6.d(end-100:end)=0;
q6.u=[q6.d; q6.r];
[Yobs,Tobs,Xobs] = lsim(syscl,q6.u,q6.T,q6.x0);
figure()
subplot(2,1,1)
hold on
stairs(Tobs,Yobs(:,1))
stairs(Tobs,q6.u(2,:),"--")
ylabel("Ouput")
legend("Output","Reference","Location","southeast","Interpreter","latex")
subplot(2,1,2)
hold on
stairs(Tobs,Yobs(:,3))
stairs(Tobs,q6.u(1,:),"--")
ylabel("Disturbance")
xlabel("Time (s)")
legend("$\hat{d}$","d","Interpreter","latex")

%% Question 7 LQ control
q7.h=0.008;
q7.sysd=c2d(q3.sysc,q7.h);
%
figure()
grid on
hold on
q7.R=0.0001;
%
q7.Q=diag([1 1 1]);
[q7.X,q7.K,q7.L] = idare(q7.sysd.A,q7.sysd.B,q7.Q,q7.R)
q7.sysdcl = ss(q7.sysd.A-q7.sysd.B*q7.K,q7.sysd.B,q7.sysd.C,[],q7.h);
q7.sysdcl = 1/dcgain(q7.sysdcl)*q7.sysdcl;
step(q7.sysdcl)
%
q7.Q=diag([10 1 1]);
[q7.X,q7.K,q7.L] = idare(q7.sysd.A,q7.sysd.B,q7.Q,q7.R)
q7.sysdcl = ss(q7.sysd.A-q7.sysd.B*q7.K,q7.sysd.B,q7.sysd.C,[],q7.h);
q7.sysdcl = 1/dcgain(q7.sysdcl)*q7.sysdcl;
step(q7.sysdcl)
%
q7.Q=diag([1 10 1]);
[q7.X,q7.K,q7.L] = idare(q7.sysd.A,q7.sysd.B,q7.Q,q7.R)
q7.sysdcl = ss(q7.sysd.A-q7.sysd.B*q7.K,q7.sysd.B,q7.sysd.C,[],q7.h);
q7.sysdcl = 1/dcgain(q7.sysdcl)*q7.sysdcl;
step(q7.sysdcl)
%
q7.Q=diag([1 1 10]);
[q7.X,q7.K,q7.L] = idare(q7.sysd.A,q7.sysd.B,q7.Q,q7.R)
q7.sysdcl = ss(q7.sysd.A-q7.sysd.B*q7.K,q7.sysd.B,q7.sysd.C,[],q7.h);
q7.sysdcl = 1/dcgain(q7.sysdcl)*q7.sysdcl;
step(q7.sysdcl)
legend('$Q\_{11}=1,Q\_{22}=1,Q\_{33}=1$',"$Q\_{11}=10,Q\_{22}=1,Q\_{33}=1$",'$Q\_{11}=1,Q\_{22}=10,Q\_{33}=1$',"$Q\_{11}=1,Q\_{22}=1,Q\_{33}=10$",'Location',"southeast",'Interpreter','latex')
%%
figure
hold on
q7.R=1;
for q=3:7
    q7.Q=diag([10^q 1 1]);
    [q7.X,q7.K,q7.L] = idare(q7.sysd.A,q7.sysd.B,q7.Q,q7.R)
    q7.sysdcl = ss(q7.sysd.A-q7.sysd.B*q7.K,q7.sysd.B,q7.sysd.C,[],q7.h);
    q7.sysdcl = 1/dcgain(q7.sysdcl)*q7.sysdcl;
    step(q7.sysdcl)
end
legend('$Q\_{11}=10^3$','$Q\_{11}=10^4$','$Q\_{11}=10^5$','$Q\_{11}=10^6$','$Q\_{11}=10^7$','Location',"southeast",'Interpreter','latex')
%%
figure
hold on
q7.Q=diag([10^7 1 1]);
for r=1:5
    q7.R=10^(2-r);
    [q7.X,q7.K,q7.L] = idare(q7.sysd.A,q7.sysd.B,q7.Q,q7.R)
    q7.sysdcl = ss(q7.sysd.A-q7.sysd.B*q7.K,q7.sysd.B,q7.sysd.C,[],q7.h);
    q7.sysdcl = 1/dcgain(q7.sysdcl)*q7.sysdcl;
    step(q7.sysdcl)
end
legend('$R=10^2$','$R=10^1$','$R=10^0$','$R=10^1$','$R=10^2$','Location',"southeast",'Interpreter','latex')
%% Question 7 final
q7.Q=diag([10^7 1 1]);
q7.R=10^(0);
[q7.X,q7.K,q7.L] = idare(q7.sysd.A,q7.sysd.B,q7.Q,q7.R)
q7.sysdcl = ss(q7.sysd.A-q7.sysd.B*q7.K,q7.sysd.B,q7.sysd.C,[],q7.h);
q7.fwd=1/dcgain(q7.sysdcl);
q7.sysdcl = 1/dcgain(q7.sysdcl)*q7.sysdcl;
figure
step(q7.sysdcl)
stepinfo(q7.sysdcl,'SettlingTimeThreshold',0.01)
%% Question 8 input limit so far
% Question 1 reference tracking CT
figure("Name","q1 reference tracking");
hold on
step(q1.Cc/(1+q1.Lc));
plot([0,1],[1,1],'--r')
plot([0,1],[-1,-1],'--r')
%xlabel("Time (seconds)")
ylabel("Input")
grid on
%%
% Question 2 disturbance rejection CT
figure("Name","q2 disturbance rejection");
hold on
step(1/(1+q2.C*G));
plot([0,10],[1,1],'--r')
plot([0,10],[-1,-1],'--r')
ylim([-1.1 1.1])
%xlabel("Time (seconds)")
ylabel("Input")
grid on
%%
% Question 4 reference tracking DT (tustin)
figure("Name","q4 reference tracking");
hold on
step(q4.Con1tust/(1+q4.Ltust1));
plot([0,1],[1,1],'--r')
plot([0,1],[-1,-1],'--r')
%xlabel("Time (seconds)")
ylabel("Input")
grid on
%% Question 4 disturbance rejection DT
figure("Name","q4 disturbance rejection");
hold on
step(1/(1+q4.Ltust2));
plot([0,10],[1,1],'--r')
plot([0,10],[-1,-1],'--r')
ylim([-1.1 1.1])
%xlabel("Time (seconds)")
ylabel("Input")
grid on
%% Question 5 pole placement
figure("Name","q5 pole placement-tracking");
hold on
[Y,T,X] = step(1/dcgain(q5.sysdfb)*q5.sysdfb);
stairs(T,1/dcgain(q5.sysdfb)-1/dcgain(q5.sysdfb)*q5.K*X')
plot([0,1],[1,1],'--r')
plot([0,1],[-1,-1],'--r')
xlabel("Time (seconds)")
ylabel("Input")
grid on
% Question 6 Observer
% Question 7 LQ
figure("Name","q7 LQ tracking");
hold on
[Y,T,X] = step(q7.sysdcl);
stairs(T,q7.fwd-q7.K*X'*q7.fwd)
plot([0,1],[1,1],'--r')
plot([0,1],[-1,-1],'--r')
xlabel("Time (seconds)")
ylabel("Input")
grid on
%% Question 9 input constraint PID
% reference tracking
q9.Kc = 1;
q9.Kp = 1;
q9.Ki = 0;
q9.Tau1 = 0.02; %0.02
q9.Tau2 = 0.01; %0.01
q9.Kd1 = (-q9.Kp*q9.Tau1^2+0.7*q9.Tau1-0.1)/(q9.Tau1-q9.Tau2);
q9.Kd2 = 0.7-q9.Kp*q9.Tau1-q9.Kp*q9.Tau2-q9.Kd1;
q9.C = PIDD(q9.Kc,q9.Kp,q9.Ki,q9.Kd1,q9.Kd2,q9.Tau1,q9.Tau2);  

figure
step(q9.C*G/(1+q9.C*G))
grid on

figure
hold on
step(q9.C/(1+q9.C*G))
plot([0,1],[1,1],'--r')
plot([0,1],[-1,-1],'--r')
grid on

%%
q9.Kc=0.23
q9.Tau1 = 0.2; %0.02
q9.Tau2 = 0.19; %0.01
q9.Kd1 = (-q9.Kp*q9.Tau1^2+0.7*q9.Tau1-0.1)/(q9.Tau1-q9.Tau2);
q9.Kd2 = 0.7-q9.Kp*q9.Tau1-q9.Kp*q9.Tau2-q9.Kd1;
q9.C = PIDD(q9.Kc,q9.Kp,q9.Ki,q9.Kd1,q9.Kd2,q9.Tau1,q9.Tau2);  
%
stepinfo(q9.C*G/(1+q9.C*G))

figure
step(q9.C*G/(1+q9.C*G))
grid on

figure
hold on
step(q9.C/(1+q9.C*G))
plot([0,10],[1,1],'--r')
plot([0,10],[-1,-1],'--r')
ylim([-1.1 1.1])
grid on
%%
q9.h = 0.1
q9.Cd=c2d(q9.C,q9.h,'Tustin');
q9.G=c2d(G,q9.h)
stepinfo(q9.Cd*q9.G/(1+q9.Cd*q9.G),"SettlingTimeThreshold",0.01)
figure();
step(q9.Cd*q9.G/(1+q9.Cd*q9.G))
grid on
figure
hold on
step(q9.Cd/(1+q9.Cd*q9.G))
plot([0,10],[1,1],'--r')
plot([0,10],[-1,-1],'--r')
ylim([-1.1 1.1])
grid on
%% distrubance rejection
%% Question 10 input constraint pole placement
q10.zeta=q5.zeta

figure
hold on
q10.w=0.5;
q10.G=1/(s^2+2*q10.w*q10.zeta*s+q10.w^2);
step(1/dcgain(feedback(q10.G,1))*feedback(q10.G,1))
q10.w=1;
q10.G=1/(s^2+2*q10.w*q10.zeta*s+q10.w^2);
step(1/dcgain(feedback(q10.G,1))*feedback(q10.G,1))
q10.w=2;
q10.G=1/(s^2+2*q10.w*q10.zeta*s+q10.w^2);
step(1/dcgain(feedback(q10.G,1))*feedback(q10.G,1))
step(q9.C*G/(1+q9.C*G))
grid on
legend("w=0.5","w=1","w=2","PDD reference")
%%
q10.h=0.1;
q10.sysd=c2d(q3.sysc,q10.h);
rank(ctrb(q10.sysd.A,q10.sysd.B));

q10.w = 2;
q10.G=1/(s^2+2*q10.w*q10.zeta*s+q10.w^2);
% polynoom z^2+pol1z+pol2
q10.pol1=-2*exp(-q10.zeta*q10.w*q10.h)*cos(q10.w*q10.h*sqrt(1-q10.zeta^2));
q10.pol2=exp(-2*q10.zeta*q10.w*q10.h);
q10.p1=-0.5*q10.pol1+0.5*sqrt(q10.pol1^2-4*q10.pol2);
q10.p2=-0.5*q10.pol1-0.5*sqrt(q10.pol1^2-4*q10.pol2);
q10.p=[q10.p1,q10.p2,min(pole(q10.sysd))];

%q10.p=[q10.p1,q10.p2,-0.75]

q10.K = place(q10.sysd.A,q10.sysd.B,q10.p);
q10.sysdfb=q10.sysd;
q10.sysdfb.A = q10.sysdfb.A-q10.sysd.B*q10.K; 
figure
hold on
step(1/dcgain(q10.sysdfb)*q10.sysdfb)
step(1/dcgain(feedback(q10.G,1))*feedback(q10.G,1))
step(q9.C*G/(1+q9.C*G))
legend("DT - Poleplacement","CT - 2^{nd} order reference","CT- Reference Q9","Location","southeast")
stepinfo(1/dcgain(q10.sysdfb)*q10.sysdfb)
grid on
%% input
figure("Name","q10 pole placement-tracking");
hold on
[Y,T,X] = step(1/dcgain(q10.sysdfb)*q10.sysdfb);
stairs(T,1/dcgain(q10.sysdfb)-1/dcgain(q10.sysdfb)*q10.K*X')
plot([0,max(T)],[1,1],'--r')
plot([0,max(T)],[-1,-1],'--r')
ylim([-1.1 1.1])
xlabel("Time (seconds)")
ylabel("Input")
grid on

%% Question 11 input constrained LQ controller - tracking
q11.sysd = q7.sysd;
q11.h = q7.h;
figure("Name","Reference Q7");
subplot(2,1,1)
step(q7.sysdcl)
subplot(2,1,2)
[Y,T,X] = step(q7.sysdcl);
plot([0,1],[1,1],'--r')
plot([0,1],[-1,-1],'--r')
plot(T,q7.fwd-q7.K*X'*q7.fwd)

%%
figure()
hold on
grid on
q11.R=10^0;
for i=1:5
    q11.Q=diag([10^(5-i) 1 1]);
    [q11.X,q11.K,q11.L] = idare(q11.sysd.A,q11.sysd.B,q11.Q,q11.R)
    q11.sysdcl = ss(q11.sysd.A-q11.sysd.B*q11.K,q11.sysd.B,q11.sysd.C,[],q11.h);
    q11.fwd=1/dcgain(q11.sysdcl);
    q11.sysdcl = 1/dcgain(q11.sysdcl)*q11.sysdcl;
    step(q11.sysdcl)
end
legend("$Q\_{11}=10^4$","$Q\_{11}=10^3$","$Q\_{11}=10^2$","$Q\_{11}=10^1$","$Q\_{11}=10^0$",'Location',"southeast",'Interpreter','latex')
%%
figure()
grid on
hold on
for i=1:5
    q11.Q=diag([10^(5-i) 1 1]);
    [q11.X,q11.K,q11.L] = idare(q11.sysd.A,q11.sysd.B,q11.Q,q11.R)
    q11.sysdcl = ss(q11.sysd.A-q11.sysd.B*q11.K,q11.sysd.B,q11.sysd.C,[],q11.h);
    q11.fwd=1/dcgain(q11.sysdcl);
    q11.sysdcl = 1/dcgain(q11.sysdcl)*q11.sysdcl;
    [Y,T,X] = step(q11.sysdcl);
    plot(T,q11.fwd-q11.K*X'*q11.fwd)
end
xlim([0 0.5])
plot([0,2],[1,1],'--r')
plot([0,2],[-1,-1],'--r')
legend("$Q\_{11}=10^4$","$Q\_{11}=10^3$","$Q\_{11}=10^2$","$Q\_{11}=10^1$","$Q\_{11}=10^0$",'Location',"northeast",'Interpreter','latex')
%% final
q11.h=0.1;
q11.sysd=c2d(q3.sysc,q11.h);
q11.R=65;
q11.Q=diag([10^2 1 1]);
[q11.X,q11.K,q11.L] = idare(q11.sysd.A,q11.sysd.B,q11.Q,q11.R);
q11.sysdcl = ss(q11.sysd.A-q11.sysd.B*q11.K,q11.sysd.B,q11.sysd.C,[],q11.h);
q11.fwd=1/dcgain(q11.sysdcl);
q11.sysdcl = 1/dcgain(q11.sysdcl)*q11.sysdcl;
[Y,T,X] = step(q11.sysdcl);

figure()
subplot(2,1,1)
grid on
hold on
stairs(T,Y)
ylabel("Output")
subplot(2,1,2)
grid on
hold on;
stairs(T,q11.fwd-q11.K*X'*q11.fwd)
plot([0,max(T)],[1,1],'--r')
plot([0,max(T)],[-1,-1],'--r')
ylabel("Input")
xlabel("Time (seconds)")
stepinfo(q11.sysdcl,'SettlingTimeThreshold',0.01)


%% Question 12 steady state error

%% Question 13 time delay
% CT PDD tracking
q13.h=0.1;
q13.tau=q13.h; %q13.h
q13.G=exp(-q13.tau*s)*G;
figure();
hold on
step(G*q9.C/(1+G*q9.C));
step(q13.G*q9.C/(1+q13.G*q9.C));
grid on
legend("$\tau=0$","$\tau=0.1$","interpreter","latex")
% DT PDD tracking
q13.sysd=c2d(q13.G,q13.h);
figure
hold on
step(q9.G*q9.Cd/(1+q9.G*q9.Cd));
step(q13.sysd*q9.Cd/(1+q13.sysd*q9.Cd));
grid on
legend("$\tau=0$","$\tau=0.1$","interpreter","latex")
stepinfo(step(q13.sysd*q9.Cd/(1+q13.sysd*q9.Cd)),"SettlingTimeThreshold",0.01)

%% redesign pdd
q13.PDD = PIDD(0.189,q9.Kp,q9.Ki,q9.Kd1,q9.Kd2,q9.Tau1,q9.Tau2)
q13.PDD = c2d(q13.PDD,q13.h,'Tustin')
figure
hold on
step(q9.G*q9.Cd/(1+q9.G*q9.Cd));
step(q13.sysd*q13.PDD/(1+q13.sysd*q13.PDD))
grid on
legend("$\tau=0$","$\tau=0.1$, mod. controller","interpreter","latex")
stepinfo(q13.sysd*q13.PDD/(1+q13.sysd*q13.PDD),"SettlingTimeThreshold",0.01)

%%
q13.sysaugA = [q4.sysd.A q4.sysd.B; zeros(1,4)];
q13.sysaugB = [zeros(3,1);1];
q13.sysaugC = [q4.sysd.C 0];
q13.ss = ss(q13.sysaugA,q13.sysaugB,q13.sysaugC,[],q13.tau)
%%
% CT PID disturbance rejection
q13.tau=0.008;
q13.h=q13.tau;
q13.G=exp(-q13.tau*s)*G;
q13.sysd=c2d(q13.G,q13.h);

figure();
hold on
step(minreal(G/(1+G*q2.C)));
step(q13.G/(1+q13.G*q2.C));
grid on
legend("$\tau=0$","$\tau=0.008$","interpreter","latex")
%%
q13.tau=0.008;
q13.h=q13.tau;
q13.G=exp(-q13.tau*s)*G;
q13.sysd=c2d(q13.G,q13.h);

figure
hold on
step(q4.sysd/(1+q4.sysd*q4.Con2tust));
step(q13.sysd/(1+q13.sysd*q4.Con2tust));
grid on
legend("$\tau=0$","$\tau=0.008$","interpreter","latex")

%%
% pole placement
q13.h=q10.h;
q13.tau=q13.h;
q13.sysd=c2d(ss(q13.G),q13.tau);
q13.sysd.A=q13.sysd.A-q13.sysd.B*q10.K;
eig(q13.sysd.A)
figure
hold on
step(1/dcgain(q10.sysdfb)*q10.sysdfb)
step(1/dcgain(q13.sysd)*q13.sysd);
grid on
legend("$\tau=0$","$\tau=0.1$","interpreter","latex")
% pole replacement
q13.sysd=c2d(ss(q13.G),q13.h);
q13.K=place(q13.sysd.A,q13.sysd.B,q10.p)
q13.sysd.A=q13.sysd.A-q13.sysd.B*q13.K;
figure
hold on
step(1/dcgain(q10.sysdfb)*q10.sysdfb)
step(1/dcgain(q13.sysd)*q13.sysd);
grid on
legend("$\tau=0$","$\tau=0.1$","interpreter","latex")
%%
[Y,T,X]=step(1/dcgain(q13.sysd)*q13.sysd);
figure
grid on
hold on;
stairs(T,1/dcgain(q13.sysd)-q13.K*X'*(1/dcgain(q13.sysd)))
plot([0,max(T)],[1,1],'--r')
plot([0,max(T)],[-1,-1],'--r')
%% Observer delay
% q13.h=q6.h;
% q13.tau=q13.h;
% q13.sysd=c2d(ss(q13.G),q13.tau);
% q13.Kc=q13.K;
% q13.Lo=q6.Lo;
% q13.fwd = 1/dcgain(ss(q13.sysd.A-q13.sysd.B*q13.Kc,q13.sysd.B,q13.sysd.C,q13.sysd.D,q13.h));
% 
% q13.Acl = [q13.sysd.A, -q13.sysd.B*q13.Kc -q13.sysd.B; q6.Lo*q13.sysd.C, q13.sysd.A-q6.Lo*q13.sysd.C-q13.sysd.B*q13.Kc zeros(3,1); q6.Ldis*q13.sysd.C -q6.Ldis*q13.sysd.C 1];
% q13.Bcl = [q13.sysd.B q13.sysd.B*q13.fwd; zeros(size(q13.sysd.B)) q13.sysd.B*q13.fwd; 0 0];
% q13.Ccl =  [1 0 0 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 0 0 1];
% q13.syscl=ss(q13.Acl,q13.Bcl,q13.Ccl,[],q13.tau);
% 
% q6.u(1,:)=zeros(1,length(q6.T))
% 
% [Yobs,Tobs,Xobs] = lsim(q13.syscl,q6.u,q6.T,q6.x0);
% figure()
% subplot(2,1,1)
% hold on
% stairs(Tobs,Yobs(:,1))
% stairs(Tobs,q6.u(2,:),"--")
% ylabel("Ouput")
% legend("Output","Reference","Location","southeast","Interpreter","latex")
% subplot(2,1,2)
% hold on
% stairs(Tobs,Yobs(:,3))
% stairs(Tobs,q6.u(1,:),"--")
% ylabel("Disturbance")
% xlabel("Time (s)")
% legend("$\hat{d}$","d","Interpreter","latex")
%%
q13.h=q11.h;
q13.tau=q13.tau;
q13.sysd=c2d(ss(q13.G),q13.tau);
q13.sysd.A=q13.sysd.A-q13.sysd.B*q11.K;

figure
hold on
step(1/dcgain(q11.sysdcl)*q11.sysdcl)
step(1/dcgain(q13.sysd)*q13.sysd);
grid on
legend("$\tau=0$","$\tau=0.1$","interpreter","latex")

q13.sysd=c2d(ss(q13.G),q13.tau);
q13.Q=diag([2000,1,1]);
q13.R=65;
[X q13.K L]=idare(q13.sysd.A,q13.sysd.B,q13.Q,q11.R);
q13.sysd.A=q13.sysd.A-q13.sysd.B*q13.K;
figure
hold on
step(1/dcgain(q11.sysdcl)*q11.sysdcl)
step(1/dcgain(q13.sysd)*q13.sysd);
grid on
legend("$\tau=0$","$\tau=0.1$","interpreter","latex")

[Y,T,X]=step(1/dcgain(q13.sysd)*q13.sysd);
figure
grid on
hold on;
stairs(T,1/dcgain(q13.sysd)-q13.K*X'*(1/dcgain(q13.sysd)))
plot([0,max(T)],[1,1],'--r')
plot([0,max(T)],[-1,-1],'--r')
%% 
close all