function [u,J,candidates_size,minJArray,maxJArray,exitflags,exec_time] = func_execute(country_val,k_val,tau_val,P_index,Flags_index,terminalcost_index,varargin)
%% *Algorithm Control Params*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Country - SouthAfrica
country                 = country_val;
% Total Number of intervals
k                       = k_val;
% Length of intervals [days]
tau                     = tau_val;
% Initial State Vector
P0                      = zeros(4^4,1);
P0(P_index)             = 1; % [HR R2 R1 wt] = [0 0 0 h]
% Flags - Which algorithms to run
Exec_Flag               = false(1,10); 
Exec_Flag(Flags_index)  = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% *Loading System Params*
% System Dynamics params ....
run('SystemParams_NEW2')
%% *Calculating other parameters*
% Generators L
% Transfer operator A
% State Cost Vectors in each interval
% Treatment Cost Vectors in each interval
% A_max A_min
% Upper- and lower bound of Co-state Vectors
% State and Co-State Trajectories For Path without switches
% Initial Upperbound of J*
run('CalculationParams')
%C_max
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization of values to be returned :
u               = [];
J               = [];
candidates_size = [];
minJArray       = [];
maxJArray       = [];
exitflags       = [];
exec_time       = NaN;
%% *Find Optimal Policy*
% Brute force algorithm
if k < 16 && Exec_Flag(1)
    tic
    [u1,~,~,~,~,~,J1,~] = runBruteForceF(A,P0,c,k,cs_integral,ct_integral,treatments,states);
    time_BF = toc;
%   display(time_BF);
%   display(u1);
    u = u1;
    J = J1; 
    exec_time = time_BF;
elseif k > 16 && Exec_Flag(1)
    display('Execution Aborted - k larger than allowed !!!')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% *Backward - box constraint parallel*
if Exec_Flag(9)
    tic
    [u_box,~,~,J_box,~,~,~,~,Q_size_box,minJArray_box,maxJArray_box,~,~,exitflags_box] = ...
       runBDA_Parallel(A,P0,c,k,cs_integral,ct_integral,treatments,states,A_max_series,A_min_series,maxJ,P11,P22,J11,J22,DeathProb,h_R1_index,h_R1,h_R2_index,h_R2);
    time_BB = toc;
%   display(time_BB);
%   display(u_box);
%   Storing value for return
    u = u_box;
    J = J_box;
    candidates_size = Q_size_box;
    minJArray = minJArray_box;
    maxJArray = maxJArray_box;
    exitflags = exitflags_box;
    exec_time = time_BB;
end
%eof