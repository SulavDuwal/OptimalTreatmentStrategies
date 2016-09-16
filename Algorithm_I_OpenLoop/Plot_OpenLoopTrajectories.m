%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author  : Sulav
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear memory
clear all
% close graphics
%close all
%% Initialization :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Country - Germany or SouthAfrica
country                 = 'SouthAfrica';
% Total Number of intervals
k                       = 5000;
% Length of intervals [days]
tau                     = 1;
% Initial State Vector
P0                      = zeros(4^4,1);
P0(4)                   = 1; % [HR R2 R1 wt] = [0 0 0 h]
% Terminal Cost Vector
% Required to be set but does not play a role in plotting or the
% computation :
terminalcost_index      = 3; 
varargin                = cell({999999999999}); % dummy setting
% Transmission rate per 100 person-years - estimated from Attia-2009
TransmissionRate        = [0.20 1.85 13.18]./36500; %% Per Person, per day
% 
timepoints              = (0:k)*tau; %% [days]
%% *Loading System Params*
% System Dynamics params ....
run('SystemParams_NEW2')
%% *Calculating other parameters*
% Generators L
% Transfer operator A
%% *Generate Generators*
G_null = MakeGenerator(k_HtoM_max, k_MtoL, k_MtoH, k_LtoM, 0, 0, mu_H1, mu_M1, mu_H2, mu_M2, NrOfStatesPerStrain,NrOfStrains,SaturationParam,PatientDeathRates);
% L{1}=full(G)';
% treatment 1
G1 = MakeGenerator(k_HtoM_max, k_MtoL, k_MtoH, k_LtoM, eta1, 0, mu_H1, mu_M1, mu_H2, mu_M2, NrOfStatesPerStrain,NrOfStrains,SaturationParam,PatientDeathRates);%  treatment #1
L{1}=full(G1)';
% treatment 2
G2 = MakeGenerator(k_HtoM_max, k_MtoL, k_MtoH, k_LtoM, 0, eta2, mu_H1, mu_M1, mu_H2, mu_M2, NrOfStatesPerStrain,NrOfStrains,SaturationParam,PatientDeathRates);%  treatment #2
L{2}=full(G2)';
%% *Transfer Operator, Cost vector of treatment and states* 
% Transfer matrices from generators : 
A = cell(1,2);
for i = 1:treatments
    A{i} = (expm(tau*L{i}));
end
%% Define the control :
u_opt     = horzcat(1*ones(1,k-14),2*ones(1,k-14));
% Trajectory of optimal state trajectory - P :
P_opt          = zeros(states,k+1);
P_opt(:,1)     = P0; %initialize recursion
%% State Vector Propagation :
for i = 1:k
    P_opt(:,i+1) = A{u_opt(i)}*P_opt(:,i);  
end
%% Classify into H,M,L
times = (0:k)*tau;
eval_vector = [64 16 4 1];
% species = {'H','R2','R1','W'};
L_prob_mat = [];
M_prob_mat = [];
H_prob_mat = [];
for i = 0:3
    for j = 0:3
        for ii = 0:3
            for kk = 0:3
                tempvar1 = [i j ii kk];
                tempvar2 = eval_vector*tempvar1'+1;
                tempvar3 = max(tempvar1);
                if tempvar3 == 1
                    L_prob_mat = [L_prob_mat tempvar2];
                elseif tempvar3 == 2
                    M_prob_mat = [M_prob_mat tempvar2];
                elseif tempvar3 == 3
                    H_prob_mat = [H_prob_mat tempvar2];
                end
            end
        end
    end
end
%%
%% Optimal
Opt.prob_H = sum(P_opt(H_prob_mat,:));
Opt.prob_M = sum(P_opt(M_prob_mat,:));
Opt.prob_L = sum(P_opt(L_prob_mat,:));

Opt.ProbDeath    = P_opt(1,1:end);
Opt.ProbSurvival = 1-P_opt(1,1:end);
%%
figure(1)
hold on
% area([1 3],[1 1],'FaceColor',[1 1 1],'EdgeColor','none')
area([14 k],[1 1],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none')
plot(Opt.prob_H,'r--')
plot(Opt.prob_M,'m--')
plot(Opt.prob_L,'b-')
plot(Opt.ProbDeath,'k:')
hold off
xlabel('days after treatment initiation')
ylabel('probabilities')
set(gca,'XScale','log')
set(gca,'XTick',          [0 14 100 1000 5000]+1)
set(gca,'XTickLabel',(tau*[0 14 100 1000 5000]))
title('Pro-active strategy')
xlim([1 k+1])
MakePlot_Nice(1)
% eof