function [sw,ssw]=calculate_spentropy(T,ws) %Returns seawater and pure water entropies w.r.t. temperature T and concentration of salt ws 
c1=-4.231e02;
c2=1.463e04;
c3=-9.880e04;
c4=3.095e05;
c5=2.562e01;
c6=-1.443e-01;
c7=5.879e-04;
c8=-6.111e01;
c9=8.041e01;
c10=3.035e-01;
sw=0.1543 + (15.383*T) - ((2.996e-02)*T*T) + ((8.193e-05)*T*T*T) - ((1.370e-07)*T*T*T*T); % entropy of pure water 
ssw=sw-(ws*(c1 + (c2*ws) + (c3*ws*ws) + (c4*ws*ws*ws) + (c5*T) + (c6*T*T) + (c7*T*T*T) + (c8*ws*T) + (c9*ws*ws*T) + (c10*ws*T*T))); % entropy of saline water 
end