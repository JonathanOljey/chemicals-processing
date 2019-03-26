clear all
close all
clc

syms  Vm x

%% PLA Production
%% Adjustable
k = 0.0270
X = 0.9;             % Conversion (Molar Basis)
Na = 39000;          % Mols of GAO in (mol)
Ca = 1;              % Assumed Concentration (mol/L)
n = 1;
time = 2;            % Time (hrs) in the reactor at constant temp
tmin = time * 60;    % Residence Time (min)
Vm1 = double(solve(tmin ...
    == int(Na/(Vm*k*(Ca*(1-x))^n),x,0,X))) % Volume of the reactor (L)

k = 0.012
X = 0.9;             % Conversion (Molar Basis)
Na = 39000;          % Mols of GAO in (mol)
n = 1.15;
time = 2;            % Time (hrs) in the reactor at constant temp
tmin = time * 60;    % Residence Time (min)
Vm2 = double(solve(tmin ...
    == int(Na/(Vm*k*(Na/Vm*(1-x))^n),x,0,X))) % Volume of the reactor (L)




