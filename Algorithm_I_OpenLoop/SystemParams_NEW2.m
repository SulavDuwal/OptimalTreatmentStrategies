%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Params
% The script assigns the value to variables :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% *System Dynamic Parameters*
% Growth rates
k_MtoH     =  0.1308; % k_m,phi
k_LtoM     =  0.2027; % k_l,phi
% % Decline rates
k_HtoM_max =  0.0662; % delta_h
k_MtoL     =  0.0113; % delta_m
% Propensity of mutation
mu_H1           = 0.1739; % mu_R1,phi
mu_M1           = 0.1739; % mu_R1,phi
mu_H2           = 0.0254; % mu_R2,phi
mu_M2           = 0.0254; % mu_R2,phi
% Efficacy of the treatments
eta1            = 9.894e-01; % Treatment 1 - eta(a1,{W1,R2})
eta2            = 9.825e-01; % Treatment 2 - eta(a2,{W1,R1})
SaturationParam = 0;
%% *Param for States*
NrOfStatesPerStrain = 4; %[0 L M H]
NrOfStrains = 4; % [HR M1 M2 Wt]
%% *Country Specific parameters
%% *South Africa*
% Country Data
SouthAfrica.GDP                             = 6620; %8300; % GDP in the considered country in USD
SouthAfrica.LossOfProductivityH             = 0.4;  % 30 - 43% according to Sendi et al (2007) PLoS One 
SouthAfrica.LossOfProductivityM             = 0.1;  % ~10% according to Sendi et al (2007) PLoS One 
SouthAfrica.life_expectancy_overall         = 49.3; %South Africa
SouthAfrica.age_at_detection                = 25; %http://www.aids.org.za/page/hiv-prevalence-south-africa
SouthAfrica.residual_life_expectancy_StateH = 5;% y ??
SouthAfrica.residual_life_expectancy_StateM = 10; % y ?
SouthAfrica.prevalence                      = 0.2;

% action 1 (no treatment)
SouthAfrica.c_action1                       = 0; % Costs for treatment per day

% action (treatment 1)
SouthAfrica.c_action2                       = 0.3;% Costs for treatment per day UsD

% action (treatment 2) 
SouthAfrica.c_action3                       = 1.08; % Costs for treatment per day 

% Inflation
SouthAfrica.lambda                          = 5.4/(365*100);%[1/day] ~5.4% inflation/year
%%
if strcmp(country,'SouthAfrica')
    GDP                             = SouthAfrica.GDP;
    LossOfProductivityH             = SouthAfrica.LossOfProductivityH;
    LossOfProductivityM             = SouthAfrica.LossOfProductivityM;
    life_expectancy_overall         = SouthAfrica.life_expectancy_overall;
    age_at_detection                = SouthAfrica.age_at_detection;
    residual_life_expectancy_StateH = SouthAfrica.residual_life_expectancy_StateH;
    residual_life_expectancy_StateM = SouthAfrica.residual_life_expectancy_StateM;
    prevalence                      = SouthAfrica.prevalence;
    %action 1 (no treatment)
    c_action1                       = SouthAfrica.c_action1; 
    %action (treatment 1)
    c_action2                       = SouthAfrica.c_action2;
    %action (treatment 2)
    c_action3                       = SouthAfrica.c_action3;
    % Discount factor based on inflation
    lambda                          = SouthAfrica.lambda;
end
%% *lifetime and death rate computations*
ResidualLifeExpectancy_overall = (life_expectancy_overall-age_at_detection);% years
ResidualLifeExpectancy_healthy = (ResidualLifeExpectancy_overall - prevalence*residual_life_expectancy_StateH)/(1-prevalence);% years

deathrate_healthy = 1/(ResidualLifeExpectancy_healthy*365);  % [1/days]
deathrate_StateM  = 1/(residual_life_expectancy_StateM*365); % [1/days]
deathrate_StateH  = 1/(residual_life_expectancy_StateH*365); % [1/days]

PatientDeathRates =  [deathrate_healthy deathrate_StateM deathrate_StateH]; % [1/days] not in publications/analysis from data sources
%%
% Number of treatment options
treatments = 2;
% Number of states
states     = 256;
% treatment cost of different options
ct         = [c_action2 c_action3];
% Terminal_Cost_types
Terminal_Costs = {'San',...   % No terminal cost
                  'Beta',...  
                  'Futur_NoDrug',... % Continue without drug after the horizon
                  'Futur_u1',... % Continue with u1 after the horizon
                  'Futur_u2'}; % Continue with u2 after the horizons
%
Terminal_Cost_Type = Terminal_Costs{terminalcost_index};

if any(terminalcost_index == [3,4,5])
    t_inf = varargin{1};
end
% eof