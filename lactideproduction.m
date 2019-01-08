clear all
close all
clc

%% Lactide Production Kinectics (Do not touch)
% Based on Example 17 in PLA Feedstock Patent
syms x k Vm
catweight = 340;      % Stannous Octoate Catalyst Weight (g)
mLA = 90.08;          % Atomic Mass of Lactic Acid (kg/kmol)
mass = 42500;         % Mass of the feed (g)
purity = 0.447;       % Purity of the feed (wt %)
V = 50;               % Volume of the reactor (L)
X = 0.87              % Conversion
Na = mass*purity/mLA; % kMols of Lactic Acid in (mol)
Ca = Na/V;            % Concentraton (mol/L)
time = 8;             % Time (hrs) in the reactor before distillation
                      % (does not include water removal)
tmin = time * 60;     % Residence Time (min)
k = double(solve(tmin ...
    == int(Na/(V*k*(Ca*(1-x))^2),x,0,X))) % min^-1

%% Adjustable
X = 0.87              % Conversion
Na = 1;               % Mols of Lactic Acid in (mol)
Mass_cat = Na/mLA*340/42500 % Mass of Catalyst (g)
Ca = Na/V;            % Concentraton (mol/L)
time = 8;             % Time (hrs) in the reactor before distillation
                      % (does not include water removal)
tmin = time * 60;     % Residence Time (min)
Vm = double(solve(tmin ...
    == int(Na/(Vm*k*(Na/Vm*(1-x))^2),x,0,X))) % min^-1

