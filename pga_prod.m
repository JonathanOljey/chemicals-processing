clear all
close all
clc

syms  Vm x

%% PGA Production
%% Adjustable
% Bulk Single Order (Methanol initiator)
k = 0.5171
X = 0.9;             % Conversion (Molar Basis)
Na = 39000;          % Mols of glycolide in (mol)
Ca = 1;              % Assumed Concentration (mol/L)
Cat = Na/100;        % Mols of Catalyst in
Init = Cat/2;        % Mols of initiator in
n = 1;
time = 2;            % Time (hrs) in the reactor at constant temp
tmin = time * 60;    % Residence Time (min)
Vm1 = double(solve(tmin ...
    == int(Na/(Vm*k*(Ca*(1-x))^n),x,0,X))) % Volume of the reactor (L)
