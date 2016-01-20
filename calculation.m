clc;
clear all;
close all;
T0=25;
data=xlsread('dstar.xlsx');
[k,b]=size(data);
%--------------------------
% Get temperature, velocity values from the flux excel file
T=25;
ws=35e-03;
[rhosw,rhow,spvol]=calculate_spvol(T,ws);
rhoswd=rhosw;
rhowd=rhow;
spvold=spvol;
[hw,hsw]=calculate_spenthalpy(T,ws);
hwd=hw;
hswd=hsw;

[sw,ssw]=calculate_spentropy(T,ws);
swd=sw;
sswd=ssw;
[muw,mus]=calculate_potential(T,ws,hsw,ssw);
muwd=muw;
musd=mus;
% Thermodynamic Properties of Restricted Dead State
ws=50e-03;
[rhosw,rhow,spvol]=calculate_spvol(T,ws);
rhoswr=rhosw;
rhowr=rhow;
spvolr=spvol;
[hw,hsw]=calculate_spenthalpy(T,ws);
hwr=hw;
hswr=hsw;
[sw,ssw]=calculate_spentropy(T,ws);
swr=sw;
sswr=ssw;
[muw,mus]=calculate_potential(T,ws,hsw,ssw);
muwr=muw;
musr=mus;
rr=zeros(k,1);
ex_destroyed=zeros(k,1);
ex_eff=zeros(k,1);
%--------------------------
for n=1:1:k
%---------------------------
% Thermodynamic Properties of Saline Water at Inlet 
T1=data(n,6)-273.15;
T=T1;
[rhosw,rhow,spvol]=calculate_spvol(T,ws);
rhosw_saline_in=rhosw;
rhow_saline_in=rhow;
spvol_saline_in=spvol;
[hw,hsw]=calculate_spenthalpy(T,ws);
hw_saline_in=hw;
hsw_saline_in=hsw;
[sw,ssw]=calculate_spentropy(T,ws);
sw_saline_in=sw;
ssw_saline_in=ssw;
[muw,mus]=calculate_potential(T,ws,hsw,ssw);
muw_saline_in=muw;
mus_saline_in=mus;
%--------------------------
% Thermodynamic Properties of Saline Water at Outlet
T2=data(n,8)-273.15;
T=T2;
[rhosw,rhow,spvol]=calculate_spvol(T,ws);
rhosw_saline_out=rhosw;
rhow_saline_out=rhow;
spvol_saline_out=spvol;
[hw,hsw]=calculate_spenthalpy(T,ws);
hw_saline_out=hw;
hsw_saline_out=hsw;
[sw,ssw]=calculate_spentropy(T,ws);
sw_saline_out=sw;
ssw_saline_out=ssw;
[muw,mus]=calculate_potential(T,ws,hsw,ssw);
muw_saline_out=muw;
mus_saline_out=mus;
%-------------------------
% Thermodynamic Properties of Permeate at Inlet
T3=data(n,5)-273.15;
T=T3;
[rhosw,rhow,spvol]=calculate_spvol(T,ws);
rhosw_permeate_in=rhosw;
rhow_permeate_in=rhow;
spvol_permeate_in=spvol;
[hw,hsw]=calculate_spenthalpy(T,ws);
hw_permeate_in=hw;
hsw_permeate_in=hsw;
[sw,ssw]=calculate_spentropy(T,ws);
sw_permeate_in=sw;
ssw_permeate_in=ssw;
[muw,mus]=calculate_potential(T,ws,hsw,ssw);
muw_permeate_in=muw;
mus_permeate_in=mus;
%-------------------------
% Thermodynamic Properties of Permeate at Outlet
T4=data(n,7)-273.15;
T=T4;
[rhosw,rhow,spvol]=calculate_spvol(T,ws);
rhosw_permeate_out=rhosw;
rhow_permeate_out=rhow;
spvol_permeate_out=spvol;
[hw,hsw]=calculate_spenthalpy(T,ws);
hw_permeate_out=hw;
hsw_permeate_out=hsw;
[sw,ssw]=calculate_spentropy(T,ws);
sw_permeate_out=sw;
ssw_permeate_out=ssw;
[muw,mus]=calculate_potential(T,ws,hsw,ssw);
muw_permeate_out=muw;
mus_permeate_out=mus;
%--------------------------
% Property Calculation Complete
%--------------------------
% Calculate Physical Parameters: Surface Areas, Volumes, Mass Flow Rates, Recovery Ratio: 
% v_saline_in=data(n,1);
% v_permeate_in=data(n,2);
mass_flux_permeate=4*(data(n,10)+data(n,12));
mass_flux_saline_inlet=4*data(n,10);
mass_flux_saline_outlet=-4*data(n,12);
mass_flux_permeate_inlet=4*data(n,9);
mass_flux_permeate_outlet=-4*data(n,11);
rr(n)=100*mass_flux_permeate/mass_flux_saline_inlet;
% fprintf('\n Mass Flux Saline Inlet= %d @ %d %cC @ %d m/s',mass_flux_saline_inlet,T1,char(176),v_saline_in);
% fprintf('\n Mass Flux Permeate Inlet= %d @ %d %cC @ %d m/s',mass_flux_permeate_inlet,T3,char(176),v_permeate_in);
% fprintf('\n Mass Flux Permeate= %d',mass_flux_permeate);
%-------------------------
% Physical Parameter Calculation Complete
%-------------------------
% Calculate flow exergies and exergy destroyed across MD
ex_sal_in=((hsw_saline_in-hswr)-(298*(ssw_saline_in-sswr))+ (0.05*(musr-musd)) + (0.95*(muwr-muwd)))*mass_flux_saline_inlet;
ex_permeate_in=((hw_permeate_in-hwr)-(298*(sw_permeate_in-swr))+ (0.05*(musr-musd)) + (0.95*(muwr-muwd)))*mass_flux_permeate_inlet;
ex_sal_out=((hsw_saline_out-hswr)-(298*(ssw_saline_out-sswr))+ (0.05*(musr-musd)) + (0.95*(muwr-muwd)))*mass_flux_saline_outlet;
ex_permeate_out=((hw_permeate_out-hswr)-(298*(sw_permeate_out-sswr))+ (0.05*(musr-musd)) + (0.95*(muwr-muwd)))*(mass_flux_permeate_inlet+mass_flux_permeate);
ex_destroyed(n)=ex_sal_in+ex_permeate_in-ex_sal_out-ex_permeate_out;
ex_eff(n)=100*(1-(ex_destroyed(n)/(ex_sal_in+ex_permeate_in)));
end
fprintf('\n Recovery Ratio= %d',rr);
fprintf('\n Exergetic Efficiency = %d',ex_eff);
fprintf('\n Exergy Destroyed Across MD Module = %d',ex_destroyed);
xlswrite('dstar_exergy.xlsx',rr,1)
xlswrite('dstar_exergy.xlsx',ex_eff,2)
xlswrite('dstar_exergy.xlsx',ex_destroyed,3)
