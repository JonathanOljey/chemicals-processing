clear all 
close all 
clc

syms Go x k Vm tcal W
% All concentrations are based on the volume of the reactor unless
% otherwised stated.
%% Calcium Lactate Production Kinectics (Don't Touch)
% Based on Example 14 in PLA Feedstock Patent

temp = 160;           % Temperature (C)
mLA = 90.08;          % Atomic Mass of Lactic Acid (g/mol)
mCC = 100.0869;       % Atomic Mass of Limestone (g/mol)
mass = 1260000;       % Mass of the feed (g)
massCC = 26000;       % Mass of Limestone (g)
massh2O = 156000;     % Mass of Water (g)
purity = 0.016;       % Purity of the feed (wt %)
massLA = purity*mLA;  % mass of Lactic Acid (g) 
V = 3000;             % Volume of the reactor (L)
X = 0.99;             % Conversion
Na = mass*purity/mLA; % mols of Lactic Acid in (mol)
Nb = massCC/mCC;      % mols of Lactic Acid in (mol)
Ca = Na/V;            % Concentraton Lactic Acid(mol/L)
Cb = Nb/V;            % Concentraton Limestone(mol/L)
time = 2;             % Time (hrs) in the reactor 
tmin = time * 60;     % Residence Time (min)
k_lactate = double(solve(tmin...
    ==int(Na/(V*k*Cb*(1-x)*(Ca*(1-x))^2),x,0,X))) % L mol^-1 min^-1

%% Calcium Lactate Production Requirements (Adjustable)
time = 24;                     % Total Reaction Time (hr)
N_l_hour = 9870.56;            % Calcium Lactate requirement (mol/hr)
N_lactate = N_l_hour*time;     % Calcium Lactate requirement (mol)
N_carbonate = N_lactate/X;     % Calcium Carbonate requirement (mol)
N_lactic_acid = N_carbonate*2; % Lactic Acid requirement (mol)
X = 0.99;                      % Conversion

%% Lactic Acid Production Kinectics

%% Adjustable Parameters
rho_G = 1540;                  % Density of Glucose (g/L)
D = 1200;                      % Density of the glucose water solution (g/L) 
scale = .1;                    % Scaling Factor
Go =  1.2*rho_G*scale;         % Initial Concentration of Glucose (g/L)
                               % 20% Added Content
Xo = Go*0.001;                 % Initial Biomass Concentration (g/L)
number = 2;                    % Number of reactors
tf = 24;                       % Duration (hrs)

%% Set Parameters
prod = N_lactic_acid;        % Production requirement (mol)
Dmol =  prod/number;         % Desired amount of Lactic acid 
                             % per reactor(mol)
mwl = 90.08;                 % mol. weight of lactic acid (g/gmol)
mwg = 180.16;                % mol. weight of glucose (g/mol)
Desiredmass =  Dmol*mwl;     % Desired amount of Lactic acid 
                             % per reactor (g)
mg = 0.057;                  % m Glucose(g) /m Biomass (g)
K   = 1.171/0.012*Xo;        % Maximum Biomass (g/L)
umx = 0.687;                 % Specific Maximum Growth Rate (g/L*h)
Yxg = 0.139;                 % Yield for Biomass formation (g biomass/g glucose)
Ygl = 1.601;                 % Yield for glucose consumption per lactic acid production
                             % (g glucose/g lactic acid)
ti = 0;                      % Starting Time (hrs)



%% Final Concentration

c = log(K/Xo-1);                        % Constant
vmx = K*umx/4;                          % max growth rate (g/L*h)
lamdax = (c-2)/umx;                     % Lag (h) 
tspan = 2*tf+1;
tsp = linspace(ti,tf,tspan);
Xf = K/(1+exp(2+4*vmx/K*(lamdax-tf)))   % Biomass (g/L)
Gf =(Go+Xo/Yxg-1/Yxg*K/(1+(K/Xo-1)*exp(-4*vmx*tf/K))...
        -(mg*K^2)/(4*vmx)*log((Xo*exp(4*vmx/K*tf)+K)/K))  
                                        % Final Glucose Concentration (g/L)
Lf =  -Xo/(Yxg*Ygl)+1/(Yxg*Ygl)*K/(1+(K/Xo-1)*exp(-4*vmx*tf/K))...
        +(mg*K^2)/(4*vmx)*log((Xo*exp(4*vmx/K*tf)+K)/K)     
                                        % Final Lactic Acid Concentration (g/L)

V_l = Desiredmass/Lf                    % Volume of reactor (L)
mass_g = Go*V_l                         % Mass of Glucose Required (g)
mass_b = Xo*V_l                         % Mass of Biomass Required (g)
V_g = mass_g/rho_G                      % Volume of glucose (L)
W_g = double(vpasolve(D==(mass_g+1000*W)/(V_g+W))) 
                                        % Voume of Water in glucose 
                                        % water solution (L)
Vm3 = V_l/1000                          % Required Reactor Volume per reactor (m^3)
M_cal_water_t = (N_carbonate*mCC/massCC*massh2O); 
                                        % Total Water requirement for Calcium Carbonate
                                        % (g)
M_cal_water = M_cal_water_t/number*1/1000% Per reactor Water requirement glucos
                                        % (kg)


%% Projected Concentration
for i = 1:tspan
    Xspan(i) = K/(1+exp(c-umx*tsp(i)));  % remaining Biomass (g/L)
    Gspan(i) =Go+Xo/Yxg-1/Yxg*K/(1+(K/Xo-1)*exp(-4*vmx*tsp(i)/K))...
        -(mg*K^2)/(4*vmx)*log((Xo*exp(4*vmx/K*tsp(i))+K)/K);  % Final Glucose Concentration (g/L)
    Lspan(i) = -Xo/(Yxg*Ygl)+1/(Yxg*Ygl)*K/(1+(K/Xo-1)*exp(-4*vmx*tsp(i)/K))...
        +(mg*K^2)/(4*vmx)*log((Xo*exp(4*vmx/K*tsp(i))+K)/K); % Lactic Acid Concentration (g/L)

end

plot(tsp,Xspan,tsp,Gspan,tsp,Lspan)
legend('Bacteria Concentration','Glucose Concentration','Lactic Acid Concentration')
xlabel('Time(hrs)')
ylabel('Concentration(g/L)')


%% Calcium Carbonate Calculation

tcal = double(solve(tcal...
    ==int(N_carbonate/(V_l*k_lactate*(N_carbonate/V_l)*(1-x)*(N_lactic_acid/V_l*(1-x))^2),x,0,X))) % L mol^-1 min^-1