function [rhosw,rhow,spvol]=calculate_spvol(T,ws) %Returns seawater and pure water densities and specific volume w.r.t. temperature T and concentration of salt ws 
a1=8.020e02;
a2=-2.001;
a3=1.667e-02;
a4=-3.060e-05;
a5=-1.613e-05;
rhow=(9.999e02)+((2.034e-2)*T)-((6.162e-03)*T*T)+((2.261e-05)*T*T*T)-((4.657e-08)*T*T*T*T);
rhosw=rhow+(ws*(a1 + (a2*T) + (a3*T*T) + (a4*T*T*T) + (a5*ws*T*T)));
spvol=1/(rhosw);
end