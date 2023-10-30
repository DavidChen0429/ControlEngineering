syms s;
numerator = 4000;
denominator = [1, 30, 200, 0];
L = tf(numerator, denominator); 

Kp = 0.49;
Ki = 0;
Kd = 0.08;
Tf = 0.01;
PD = pid(Kp, Ki, Kd, Tf);

P_ControlledSys = L * PD;
P_ControlledCLSys = feedback(P_ControlledSys, 1);
stepinfo(P_ControlledCLSys)
subplot(2,1,1)
step(P_ControlledCLSys);
subplot(2,1,2)
margin(P_ControlledSys)
pole(P_ControlledCLSys)