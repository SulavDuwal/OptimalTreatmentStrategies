% Code Courtesy From Stefanie Winkelmann
clear all
%New Paramaters
k_MtoH     = 0.1308;
k_LtoM     = 0.2027;
k_HtoM_max = 0.0662;
k_MtoL     = 0.0113;
mu_H1      = 0.1739; 
mu_M1      = 0.1739;
mu_H2      = 0.0254;
mu_M2      = 0.0254;
eta1       = 0.9894;
eta2       = 0.9825;
SaturationParam = 0;
%
NrOfStatesPerStrain = 4; %[0 L M H]
NrOfStrains = 4; % [HR M1 M2 Wt]
% define parameters

%South Africa
% Country Dates
GDP                             = 6620; % 2013
LossOfProductivityH             = 0.4;         % 30 - 43% according to Sendi et al (2007) PLoS One 
LossOfProductivityM             = 0.1;         % ~10% according to Sendi et al (2007) PLoS One 
life_expectancy_overall         = 49.3; % South Africa
age_at_detection                = 25;           % http://www.aids.org.za/page/hiv-prevalence-south-africa
residual_life_expectancy_StateH = 5;  % y ??
residual_life_expectancy_StateM = 10; % y ?
prevalence = 0.2;

k_info  = 200; %USd
% action 1 (no treatment)
c_action1 = 0; % Costs for treatment per day
% action (treatment 1)
c_action2 = 0.3;  % Costs for treatment per day UsD
% action (treatment 2) 
c_action3 = 1.08; % Costs for treatment per day
lambda=5.4/(365*100);%[1/day] ~5.4% inflation/year
%% lifetime and death rate computations

ResidualLifeExpectancy_overall = (life_expectancy_overall-age_at_detection);% years
ResidualLifeExpectancy_healthy = (ResidualLifeExpectancy_overall - prevalence*residual_life_expectancy_StateH)/(1-prevalence);% years

deathrate_healthy = 1/(ResidualLifeExpectancy_healthy*365);  % [1/days]
deathrate_StateM  = 1/(residual_life_expectancy_StateM*365); % [1/days]
deathrate_StateH  = 1/(residual_life_expectancy_StateH*365); % [1/days]

PatientDeathRates =  [deathrate_healthy deathrate_StateM deathrate_StateH];% [1/days] not in publications/analysis from data sources

%%
% optimization algorithm 
t_stepsize=1;%days
N=1000;

c_action=[c_action1 c_action2 c_action3]; % for cost splitting in Nettocosts.m

% Start calculation

c_state=MakeCostfunction(GDP,LossOfProductivityH,LossOfProductivityM);

c(:,1)=c_state+c_action1;
c(:,2)=c_state+c_action2;
c(:,3)=c_state+c_action3;
c(1,:)=c_state(1);% Beim Tod kosten die Medikamente nichts

%% Create Generator Matrix
G = MakeGenerator(k_HtoM_max, k_MtoL, k_MtoH, k_LtoM, 0, 0, mu_H1, mu_M1, mu_H2, mu_M2, NrOfStatesPerStrain,NrOfStrains,SaturationParam,PatientDeathRates);
L(:,:,1)=full(G);
G = MakeGenerator(k_HtoM_max, k_MtoL, k_MtoH, k_LtoM, eta1, 0, mu_H1, mu_M1, mu_H2, mu_M2, NrOfStatesPerStrain,NrOfStrains,SaturationParam,PatientDeathRates);%  treatment #1
L(:,:,2)=full(G);
G = MakeGenerator(k_HtoM_max, k_MtoL, k_MtoH, k_LtoM, 0, eta2, mu_H1, mu_M1, mu_H2, mu_M2, NrOfStatesPerStrain,NrOfStrains,SaturationParam,PatientDeathRates);%  treatment #2
L(:,:,3)=full(G);

%% Function to Optimize the Algorithm
[V,action,t_opt] = DisInfoPolicyFunction(L,c,k_info,lambda,t_stepsize,N);
% Ausgabe
for hr = 0:3
    for m2 = 0:3
        for m1 = 0:3
            for wt = 0:3
                CurrentState = [hr m2 m1 wt];
                StateIndex = CurrentState(1)*4^3 + CurrentState(2)*4^2 + CurrentState(3)*4^1 + CurrentState(4)*4^0 + 1;
                State(StateIndex,:) = CurrentState;
                Strategy(StateIndex,:) = [State(StateIndex,:) action(StateIndex) t_opt(StateIndex) V(StateIndex)];
            end
        end
    end
end

[A,v1,v2,v3,V_info, V_action, V_state, V_total] = Nettocosts(L,lambda,t_opt,action,c_action,c_state,c,k_info);

% Results
c = '';
dlmwrite(strcat('./Results/Strategy',c,'.txt'),Strategy);
% 2) Parameter
Params = [k_MtoH ;k_LtoM;k_HtoM_max;k_MtoL;mu_H1;mu_M1;mu_H2;mu_M2;eta1;eta2;SaturationParam; GDP;LossOfProductivityH;LossOfProductivityM;k_info;c_action1;c_action2;c_action3;lambda;t_stepsize;N;ResidualLifeExpectancy_healthy;age_at_detection;residual_life_expectancy_StateH;residual_life_expectancy_StateM];
dlmwrite(strcat('./Results/Parameters',c,'.txt'),Params);

% 3) Kosten-Splitting
dlmwrite(strcat('./Results/CostSplitting',c,'.txt'),[V_info, V_action, V_state, V_total]);
