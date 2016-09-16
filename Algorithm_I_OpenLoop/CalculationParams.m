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
%% *Cost of action & state*
c_states = MakeCostfunction(GDP,LossOfProductivityH,LossOfProductivityM)';
% Expected State Cost vector  
Id = eye(states);
L_lambdaId{1} = L{1} - lambda*Id;
L_lambdaId{2} = L{2} - lambda*Id;
LexpmL_Id_tau{1} = L_lambdaId{1}\(expm(L_lambdaId{1}*tau)-Id); 
LexpmL_Id_tau{2} = L_lambdaId{2}\(expm(L_lambdaId{2}*tau)-Id);
cs_integral = NaN(states,1,k,treatments);
for tr = 1:treatments
    for t = 1:k    
        cs_integral(:,:,t,tr) = exp(-lambda*(t-1)*tau)*c_states'*LexpmL_Id_tau{tr};
    end
end

% Expected Treatment Cost vector
deathpenalty = ones(states,1);
deathpenalty(1) = 0;
ct_integral = NaN(states,1,k,treatments);
for tr = 1:treatments
    for t = 1:k
        ct_integral(:,:,t,tr) = exp(-lambda*(t-1)*tau)*ct(tr)*deathpenalty'*LexpmL_Id_tau{tr};
    end
end

% terminal cost vector
switch Terminal_Cost_Type
    case 'San' % No Terminal cost
        c_s = zeros(states,1);
        c   = zeros(states,1);
    case 'Beta' % Terminal cost scaled by beta
        c_s = MakeTerminalCostfunction4(c_states,tau*k,50000,lambda)';
        c = c_s;
    case 'Futur_NoDrug' % Treatment discontinued
        c_a = zeros(states,1);
        c_s  = MakeTerminalCostfunction3(G_null,c_states,tau*k,t_inf-tau*k,lambda)';
        c = c_a + c_s;
    case 'Futur_u1' % Treatment u1 continued
        c_a = MakeTerminalCostfunction3(G1,ct(1)*deathpenalty,tau*k,t_inf-tau*k,lambda)';
        c_s  = MakeTerminalCostfunction3(G1,c_states,tau*k,t_inf-tau*k,lambda)';
        c = c_a + c_s;
    case 'Futur_u2' % Treatment u2 continued
        c_a = MakeTerminalCostfunction3(G2,ct(2)*deathpenalty,tau*k,t_inf-tau*k,lambda)';
        c_s  = MakeTerminalCostfunction3(G2,c_states,tau*k,t_inf-tau*k,lambda)';
        c = c_a + c_s;
end
%% *A max and A min*
A_max = NaN(states,states);
A_min = NaN(states,states);
filter_A = A{1} >= A{2};
% Amax
A_max(filter_A)  = A{1}(filter_A);
A_max(~filter_A) = A{2}(~filter_A);
% Amin
A_min(filter_A)  = A{2}(filter_A);
A_min(~filter_A) = A{1}(~filter_A);
A_max_series = cell(1,k);
A_min_series = cell(1,k);
for i = 1:k
    % k-i+1
    A_max_series{i} = A_max^(i);
    A_min_series{i} = A_min^(i);
end
%% Treatment without switches
% treatment 1
P11 = NaN(states,k+1);
P11(:,1) = P0;
R11 = zeros(1,k+1);
C11 = NaN(states,k+1);
C11(:,1) = c;
u11 = ones(1,k);
for i = 1:k
    P11(:,i+1) = A{u11(i)}*P11(:,i);
    R11(i+1)   = cs_integral(:,:,i,u11(i))'*P11(:,i) + ct_integral(:,:,i,u11(i))'*P11(:,i);
    C11(:,i+1) = ( C11(:,i)'*A{u11(i)} + cs_integral(:,:,k-i+1,u11(i))' + ct_integral(:,:,k-i+1,u11(i))' )';
end
J11 = cumsum(R11); % Cummulative running cost
J11(end) = J11(end) + c'*P11(:,k+1); % Terminal cost
clear u11 
% treatment 2
P22 = NaN(states,k+1);
P22(:,1) = P0;
R22 = zeros(1,k+1);
C22 = NaN(states,k+1);
C22(:,1) = c;
u22 = 2*ones(1,k);
for i = 1:k
    P22(:,i+1) = A{u22(i)}*P22(:,i);
    R22(i+1)   = cs_integral(:,:,i,u22(i))'*P22(:,i) + ct_integral(:,:,i,u22(i))'*P22(:,i);
    C22(:,i+1) = (C22(:,i)'*A{u22(i)} + cs_integral(:,:,k-i+1,u22(i))' + ct_integral(:,:,k-i+1,u22(i))' )';
end
J22 = cumsum(R22); % Cummulative running cost
J22(end) = J22(end) + c'*P22(:,k+1); % Terminal cost
clear u22
%% Upperbound of J*
maxJ = min(J11(end),J22(end));
%% *Maximum Death rate probability*
DeathProb = zeros(1,k+1);
L_null = full(G_null)';
A_null = expm(tau*L_null);
for i = 1:k
    tempVar = A_null^k*P0;
    DeathProb(i+1) = tempVar(1);
end
%% *Maximum H-R1 and H-R2 rate*
times = (0:k)*tau;
eval_vector = [64 16 4 1];
% species = {'H','R2','R1','W'};
h_R1_index = [];
h_R2_index = [];
for i = 0:3
    for j = 0:3
        for ii = 0:3
            for kk = 0:3
                tempvar1 = [i j ii kk];
                tempvar2 = eval_vector*tempvar1'+1;
                if ii == 3
                h_R1_index = [h_R1_index tempvar2];
                end
                if j == 3
                h_R2_index = [h_R2_index tempvar2];
                end
            end
        end
    end
end
h_R1 = sum(P11(h_R1_index,end));
h_R2 = sum(P22(h_R2_index,end));
%% Upper- and lowerbound of co-state vectors
C_st_max = NaN(states,k);
C_st_min = NaN(states,k);
for t = 1:k
    tempVector1 = ct_integral(:,:,t,1) + cs_integral(:,:,t,1);
    tempVector2 = ct_integral(:,:,t,2) + cs_integral(:,:,t,2);
    filterVect_temp1 = tempVector1 >= tempVector2;
    filterVect_temp2 = tempVector2 >= tempVector1;
    
    C_st_max(filterVect_temp1,t) = tempVector1(filterVect_temp1);
    C_st_max(filterVect_temp2,t) = tempVector2(filterVect_temp2);
    
    C_st_min(filterVect_temp2,t) = tempVector1(filterVect_temp2);
    C_st_min(filterVect_temp1,t) = tempVector2(filterVect_temp1);
end

C_max = NaN(states,k+1);
C_min = NaN(states,k+1);
C_max(:,k+1) = c;
C_min(:,k+1) = c;
for i = k:-1:1
        C_max(:,i) = A_max'*C_max(:,i+1) + C_st_max(:,i);
        C_min(:,i) = A_min'*C_min(:,i+1) + C_st_min(:,i);
end

% eof
%% Upper and Lower bound of R
Rmin = zeros(k,1);
Rmax = zeros(k,1);
Rmin(2) = C_st_min(:,1)'*P0;
Rmin(2) = C_st_min(:,1)'*P0;
for i = 3:k
    Rmin(i+1) = Rmin(i) + C_st_min(:,i)'*A_min_series{i-1}*P0;
    Rmax(i+1) = Rmax(i) + C_st_min(:,i)'*A_max_series{i-1}*P0;
end

%%
max_A1 = max(A{1}');
max_A2 = max(A{2}');
P_max = max(vertcat(max_A1,max_A2))';
