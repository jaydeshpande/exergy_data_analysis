function [hw,hsw]=calculate_spenthalpy(T,ws) %Returns seawater and pure water enthalpies w.r.t. temperature T and concentration of salt ws 
b1=-2.348e04;
b2=3.152e05;
b3=2.803e06;
b4=-1.446e07;
b5=7.826e03;
b6=-4.417e01;
b7=2.139e-01;
b8=-1.991e04;
b9=2.778e04;
b10=9.728e01;
hw=141.355+(4202.070*T)-(0.535*T*T)+(0.004*T*T*T); % enthalpy of pure water
hsw=hw-(ws*(b1 + (b2*ws) + (b3*ws*ws) + (b4*ws*ws*ws) + (b5*T) + (b6*T*T) + (b7*T*T*T) + (b8*ws*T) + (b9*ws*ws*T) + (b10*ws*T*T))); % enthalpy of saline water
end