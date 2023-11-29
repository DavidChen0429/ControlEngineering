function [M,Tr,Mp,Ts] = ordersys(omega0,gamma)
    k=1;
    s = tf([1 0],[1]);
    Tr=1/omega0*exp(acos(gamma)/tan(acos(gamma)));
    Mp=exp(-pi*gamma/sqrt(1-gamma^2));
    Ts=4.6/(omega0*gamma);
    M = (k*omega0^2)/(s^2+2*gamma*omega0*s+omega0^2);
end

