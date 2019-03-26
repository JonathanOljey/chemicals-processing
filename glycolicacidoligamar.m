clear all
close all
clc

%% GAO Production Kinectics (Don't touch)
% Based on Example 1 in Glycolide Patent
syms x k Vm

temp = 230;           % Temperature (C)
mGA = 76.05;          % Atomic Mass of Glycolic Acid (kg/kmol)
mass = 1000;          % Mass of the feed (g)
purity = 0.7;         % Purity of the feed (wt %)
V = 1;                % Volume of the reactor (L)
mGAO = 480;           % Mass of Oligomar (g)
X = mGAO/(mass*purity)% Conversion
Na = mass*purity/mGA; % Mols of Lactic Acid in (mol)
Ca = Na/V;            % Concentraton (mol/L)
time = 7;             % Time (hrs) in the reactor at constant temp
                      % (includes water removal)
tmin = time * 60;     % Residence Time (min)
k = double(solve(tmin ...
    == int(Na/(V*k*Ca*(1-x)),x,0,X))) % min^-1


%% Adjustable
X =  X;           % Conversion
Na = 32064;           % Mols of Glycolic Acid in (mol)
time = 7;             % Time (hrs) in the reactor at constant temp
                      % (includes water removal)
tmin = time * 60;     % Residence Time (min)

Vm = double(solve(tmin ...
    == int(Na/(Vm*k*Ca*(1-x)),x,0,X))) % min^-1


