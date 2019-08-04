function [DATA] = func(out,CR,supar,OF,Pc)

%% specify oxidizer & fuel conditions inputs
ox1   = 'O2(L)';             %% primary oxidizer
ox2   = '';             %% secondary oxidizer
fu1   = 'H2(L)';   %% primary fuel
fu2   = '';            %% secondary fuel
ox1chem   = '';         %% optional primary oxidizer chemical formula if required (captilize chemical symbols)
ox2chem   = '';         %% optional secondary oxidizer chemical formula if required (captilize chemical symbols)
fu1chem   = '';         %% optional primary fuel chemical formula if requied (captilize chemical symbols)
fu2chem   = '';         %% optional secondary fuel chemical formula if requied (captilize chemical symbols)

ox1wt = 100;               %% wt fraction of primary oxid in total oxid [1]
ox2wt = 0;              %% wt fraction of secondary oxid in total oxid [1]
fu1wt = 100;            %% wt fraction of primary oxid in total oxid [1]
fu2wt = 0;              %% wt fraction of secondary oxid in total oxid [1]
ox1T  = 90;%90.17;               %% optional input of primary oxid temperature which enthalpy is evaluated [degR]
ox2T  = 0;              %% optional input of secondary oxid temperature which enthalpy is evaluated [degR]
fu1T  = 20;%20.27;               %% optional input of primary fuel temperature which enthalpy is evaluated [degR]
fu2T  = 0;              %% optional input of secondary fuel temperature which enthalpy is evaluated [degR]
ox1H  = 0;              %% optional input of primary oxid enthalpy of formation [cal/mol]
ox2H  = 0;              %% optional input of secondary oxid enthalpy of formation [cal/mol]
fu1H  = 0;              %% optional input of primary fuel enthalpy of formation [cal/kg]
fu2H  = 0;              %% optional input of secondary fuel enthalpy of formation [cal/mol]

Pc    = [Pc*14.7/101325];             %% Chamber pressure [psia]
% OF    = [5.0];             %% mixture ratio [wt oxid/wt fuel]
Phi   = [];            %% equivalence ratio
PR    = [];             %% pressure ratios [Pc/Pe]
subar = [];             %% subsonic area ratios [A/At]
% supar = [];             %% supersonic area ratios [A/At]
% CR    = [];              %% chamber contraction ratio [Ac/At]
flow  = 'eq';           %% flow type [eq or fz]
% out   = 'cp'; %%'aeat p t';  %maximum eight output parameters in one call to CEA
%outlab= ['Area ratio' 'Temperature' 'Pressure' 'xH2O' 'xCO2' 'xCO'];
% out   = ''; %maximum eight output parameters in one call to CEA
                        % typical outputs
                        % p: pressure
                        % t: temperature
                        % rho: density
                        % h: enthalpy
                        % g: gibbs energy
                        % s: entropy
                        % m: molecular weight (per definition eq(2.3a) cea manual sp-1)
                        % mw: molecular weight (per definition eq(2.4a) cea
                        % manual sp-2
                        % cp: specific heat
                        % gam: gamma(s)
                        % son: sonic velocity
                        % pip: pressure ratio, Pinj/Pe
                        % mach: Mach number
                        % ae: area ratio
                        % cf: coefficient of thrust, Cf
                        % ivac: vacuum specific impulse, Ivac
                        % isp: specific impulse, Isp

inpgen_rocket(ox1,ox2,ox1wt,ox2wt,ox1T,ox2T,ox1chem,ox2chem,ox1H,ox2H,fu1,fu2,fu1wt,fu2wt,fu1T,fu2T,fu1chem,fu2chem,fu1H,fu2H,Pc,OF,Phi,PR,subar,supar,CR,flow,out);

%% execute CEA600
system('CEA600.exe');

%% extract output from plot file (.plt)
DATA = load('Detn.plt');


close all;

end

