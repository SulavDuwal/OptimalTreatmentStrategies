%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author  : Sulav Duwal
% Purpose : A wrapper script to execute dynamic programming algorithms
%           with backward propagation in order to find the optimal 
%           treatment strategy against HIV-1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear memory
clear all
% close graphics
close all
%% Initialization :
% Country
country_val_array  = {'SouthAfrica'};
% Parameters
k_val_array        = 15;  % Number of intervals  
tau_val_array      = 2;   % Interval length [days]
P_index_array      = 4;   % HIV status for eg. 4 is [WT == H,0,0,0]
%% Which algorithms to run :
% Brute Force Algorithm - Naive     = 1
% BDA Parallel                      = 9
Flags_index_array  = 9;  % Default is BDA Parallel (5) 
% 'Backward dynamic programming algorithm' require linear programming 
% solver.
% IBM CPLEX Suite for solver linear programming is recommended.
%
% The default is Matlab's linear programming solver 'linprog'
%% Save Results ?
% Path to Folder to save results
results_folder     = 'SavedVariables/';
% Flag to save results
Flag.SaveResults   = true;
%% Local Path in your computer to IBM ( linear programming solver ) - cplex
% Uncomment, if you have CPLEX installed.
% addpath('Defince your local path /IBM_Cplex_files/cplex/matlab')
%% If you access to matlab parallel toolbox - set the number of slaves
% if matlabpool('size') > 0
%     matlabpool close;
% end
% matlabpool 2;
%% Terminal Cost Type
% No Terminal Cost                            = 1
% Terminal Cost Scaled by beta                = 2
% Large Time Horizon - Treatment Discontinued = 3 ( Scenario 0  ) 
% Large Time Horizon - Treatment u1 continued = 4 ( Scenario I  ) 
% Large Time Horizon - Treatment u2 continued = 5 ( Scenario II ) 
for iii = [3,4,5] %%iii = [3,4,5]
    terminalcost_index = iii;
    if any(terminalcost_index == [3,4,5])
        t_inf = 3e5; % Select value above 300,000 [days]
    end
    
    for i = 1:length(country_val_array) % country
        for j = 1 % k and tau
            % Assigning value from array
            country_val = country_val_array{i};
            k_val       = k_val_array(j);
            tau_val     = tau_val_array(j);
            P_index     = P_index_array(1);
            Flags_index = Flags_index_array(1);
            % Executing the algorithms
            if any(terminalcost_index == [3,4,5])
                %display('xyz')
                [u,J,candidates_size,minJArray,maxJArray,exitflags,exec_time] = ...
                    func_execute(country_val,k_val,tau_val,P_index,Flags_index,terminalcost_index,t_inf);
            else
                [u,J,candidates_size,minJArray,maxJArray,exitflags,exec_time] = ...
                    func_execute(country_val,k_val,tau_val,P_index,Flags_index,terminalcost_index);
            end
            %J(end-1)
            % Common suffix for saving files
            common_suffix = strcat(country_val,'_k',num2str(k_val),'_tau',num2str(tau_val),'_P',num2str(P_index),'_Algo',num2str(Flags_index),'_Sce',num2str(terminalcost_index-3));
            %% Displaying result in console
            display('=======================================================================================')
            display('Settings  :')
            display('---------------------------------------------------------------------------------------')
            display(sprintf('Country = %s ; #Interval = %d ; Interval Length = %d ; Terminal Cost Type = %d', ...
                country_val,k_val,tau_val,terminalcost_index))
            display('---------------------------------------------------------------------------------------')
            display('Results   :')
            display('---------------------------------------------------------------------------------------')
            display(sprintf('                      %s','Optimal Switching Sequence'))
            display(sprintf('%d     ',u));
            display('---------------------------------------------------------------------------------------')
            display(sprintf('Total Cost = %d', J(end)));
            display('=======================================================================================')
        end
    end
    %% Saving variables
    if Flag.SaveResults
        csvwrite(strcat(results_folder,'u_',common_suffix,'.csv'),u);
        csvwrite(strcat(results_folder,'J',common_suffix,'.csv'),J);
        csvwrite(strcat(results_folder,'cand_size_',common_suffix,'.csv'),candidates_size);
        csvwrite(strcat(results_folder,'exec_time_',common_suffix,'.csv'),exec_time);
    end
end
% eof