clear all
close all
clc

%% Glycolide Production Kinectics (Don't Touch)
% Based on Example 1 in Glycolide Patent
syms x k Vm

temp = 230;           % Temperature (C)
mGA = 76.05;          % Atomic Mass of Glycolic Acid (g/mol)
mGly = 116.072;       % Atomic Mass of Glycolide (g/mol)
mass = 160;           % Mass of the feed (g)
purity = 0.7;         % Purity of the feed (wt %)
V = 0.5;              % Volume of the reactor (L)
Xmass = 0.963;        % Conversion (wt%)
X = (Xmass/mGly/...   % Conversion (Molar Basis)
    ((Xmass/mGly)+(1-Xmass)/mGA)) 
Na = mass*purity/mGA; % Mols of GAO in (mol)
Ca = Na/V;            % Concentraton (mol/L)
time = 10;            % Time (hrs) in the reactor at constant temp
tmin = time * 60;     % Residence Time (min)
k = double(solve(tmin ...
    == int(Na/(V*k*(Ca*(1-x))^2),x,0,X))) % min^-1

%% Adjustable
X = 0.99;             % Conversion (Molar Basis)
Na = 2;               % Mols of GAO in (mol)
massGAO = Na/mGA
massGAO = Na/mGA
mTEGDB = 100/160*massGAO            % Mass of TEG-DB (g)
mOTEG = 89/160*massGAO              % Mass of OTEG (g)
mOct = 0.17/160*massGAO             % Mass of Stannous (g)
time = 10;            % Time (hrs) in the reactor at constant temp
tmin = time * 60;     % Residence Time (min)
Vm = double(solve(tmin ...
    == int(Na/(Vm*k*(Ca*(1-x))^2),x,0,X))) % Volume of the reactor (L)


