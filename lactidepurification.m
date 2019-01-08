clear all
close all
clc

%% Lactide Purification Kinectics
% Based on Example 18 in Patent 2012/0149920
syms x k Vm

temp = 2;             % Temperature (C)
mLACT = 144.13;       % Atomic Mass of Lactide (kg/kmol)
mass = 47;            % Mass of the feed (g)
V = 0.1;              % Volume of the reactor (L)
X = 0.94;             % Conversion
Na = mass/mLACT;      % mols of Lactide
Nethyl = 3/5*Na;      % mols of ethyl acetate
Ca = Na/V;            % Concentraton (mol/L)
time = 4;             % Total time (hrs) purifying in the cooler
tmin = time * 60;     % Residence Time (min)
k = double(solve(tmin ...
    == int(Na/(V*k*Ca*(1-x)),x,0,X))) % mol^-3 min^-1

%% Adjustable
X = 0.94;             % Conversion
Na = 1;               % mols of Lactide
Nethyl = 3/5*Na;      % mols of ethyl acetate
time = 4;             % Total time (hrs) purifying in the cooler
tmin = time * 60;     % Residence Time (min)
Vm = double(solve(tmin ...
    == int(Na/(Vm*k*Ca*(1-x)),x,0,X))) % mol^-3 min^-1