function C = PIDD(Kc,Kp,Ki,Kd1,Kd2,Tau1,Tau2)
    s = tf([1 0],[1]);
    C = Kc*(Kp+Ki/s+Kd1*s/(Tau1*s+1)+Kd2*s/(Tau2*s+1));
end