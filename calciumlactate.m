clear all
close all
clc

%% Lactide Production Kinectics (Don't Touch)
% Based on Example 14 in PLA Feedstock Patent
syms x k Vm           % Converts to equialent mol stannous octoate (kg)
temp = 160;           % Temperature (C)
mLA = 90.08;          % Atomic Mass of Lactic Acid (kg/kmol)
mCC = 100.0869;       % Atomic Mass of Limestone (kg/kmol)
mass = 1260000;       % Mass of the feed (g)
massCC = 26000;       % Mass of Limestone (g)
purity = 0.016;       % Purity of the feed (wt %)
massLA = purity*mLA;  % mass of Lactic Acid (kg) 
V = 3000;             % Volume of the reactor (L)
X = 0.99;             % Conversion
Na = mass*purity/mLA; % kMols of Lactic Acid in (mol)
Nb = massCC/mCC;      % kMols of Lactic Acid in (mol)
Ca = Na/V;            % Concentraton Lactic Acid(mol/L)
Cb = Nb/V;            % Concentraton Limestone(mol/L)
time = 2;             % Time (hrs) in the reactor 
tmin = time * 60;     % Residence Time (min)
k = double(solve(tmin...
    ==int(Na/(V*k*Cb*(1-x)*(Ca*(1-x))^2),x,0,X))) % kmol^-1 min^-1

%% Determines Volume (Adjustable)
X = 0.99;             % Conversion
Na = 1;               % Mols of Lactic Acid in (kmol)
Nb = Na/2;            % Mols of  in (kmol)
time = 2;             % Time (hrs) in the reactor 
tmin = time * 60;     % Residence Time (min)
Vm = double(solve(tmin...
    ==int(Na/(Vm*k*Nb/Vm*(1-x)*(Na/Vm*(1-x))^2),x,0,X))) % L


