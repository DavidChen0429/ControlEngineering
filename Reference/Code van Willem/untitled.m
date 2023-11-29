function [M,Tr,Mp] = ordersys(omega0,gamma)
    Tr=1/omega0*exp(acos(gamma)/tan(acos(gamma)));
    Mp=exp(-pi*gamma/sqrt(1-gamma^2))
    M = (k*omega0^2)/(s^2+2*gamma*omega0*s+omega0^2);
end

