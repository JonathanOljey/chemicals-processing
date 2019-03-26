clear all
close all
clc
% Glycolic Acid Production
syms x Vm k

%% Don't touch
CO = 1;         % mols of CO
P = 900;        % Pressure (atm)
a = 1.505;      % a Vanderwaals constant (L^2*atm/mol^2)
b = 0.03985;    % B Vanderwaals constant (L/mol)
R = 0.08205736; % R gas constant
T = 200+273;    % Temperature (K)
Vmon = double(solve(R*T == (P+a/(Vm^2)*(Vm-b))))
Vmon = Vmon(1);     % Molar volume of Carbon Monoxide (L/mol)
% Assumption: Given molar values are already in concentration form
CCo = 1;    % Concentration of Carbon Monoxide mol/L
CForma = 1; % Concentration of Formaldehyde mol/L
Ch2o = 6;   %Concentration of Water mol/L
X = 0.755;
tmin = 60;
V = 1;
Na = CCo*V;
% Molar Ratios
A = CCo/CCo;
B = CForma/CCo;
C = Ch2o/CCo;
k = double(solve(tmin ==...
    int(Na/(k*V*CCo*(1-x/A)*CForma*(1-x/B)*Ch2o*(1-x/C)),x,0,X))) % min^-1

%%  Adjustable
NForma = 1; % Formaldehyde mol
NCo = NForma;    % Carbon Monoxide mol
Nh2o = 6*NCo;   % Water mol
X = 0.755;
tmin = 60;
Na = NForma;
% Molar Ratios
A = NCo/NCo;
B = NForma/NCo;
C = Nh2o/NCo;
Vm = double(solve(tmin ==...
    int(Na/(k*Vm*NCo/Vm*(1-x/A)*NForma/Vm*(1-x/B)*Nh2o/Vm*(1-x/C)),x,0,X))) % L