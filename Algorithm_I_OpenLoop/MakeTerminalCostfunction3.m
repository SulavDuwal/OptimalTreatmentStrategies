function [c]=MakeTerminalCostfunction3(G_null,c_states,t_start,t_end,lambda)
% Reassigning the variables
L = full(G_null)';
% Total time in future 
t = t_end - t_start;
% Calculation of Terminal cost vector
Id = eye(size(L,1));
L_lambdaId = L - lambda*Id;
LexpmL_Id_tau = L_lambdaId\(expm(L_lambdaId*t)-Id);
c = exp(-lambda*t_start)*c_states'*LexpmL_Id_tau;   
