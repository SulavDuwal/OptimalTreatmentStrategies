% Compute Information cost , action cost and state cost for given optimal
% strategy

% gegeben:
% [V,action,t_opt] = DisInfoPolicyFunction(L,c,k_info,lambda,t_stepsize,N)
% c_action: action cost
% c_state: state cost

function [A,v1,v2,v3,V_info, V_action, V_state, V_total] = Nettocosts(L,lambda,t_opt,action,c_action,c_state,c,k_info)

size_S = size(L,1); % number of states

Id = eye(size_S);

%t_opt(find(t_opt == 2000)) = 1e6;
%t_opt(find(t_opt == 2000)) = 10;

for s = 1:size_S
    e_tau(s) = exp(-lambda*t_opt(s));
    eL = expm(L(:,:,action(s))*t_opt(s));
    P(s,:) = eL(s,:);
    cc(s) = c_action(action(s));
    cc(1) = 0;
    A_state=inv(L(:,:,action(s))-lambda*Id)*(expm((L(:,:,action(s))-lambda*Id)*t_opt(s))-Id);
    A(s,:) = A_state(s,:);
    CC = A_state*c(:,action(s));
    C(s) = CC(s);
    c_a(s) = c_action(action(s))*(1/lambda*(1-e_tau(s))-A_state(s,1)); % for HIV-Model: weight costs
end

cc(1) = 0; % If a person dies, action cost is zero.

D = diag(e_tau);

INV = inv(Id-D*P);

V_info = k_info*INV*e_tau';

% For HIV-Model: 
V_action = INV*c_a';
V_state = INV*A*c_state';
V_total = INV*(C'+k_info*e_tau');

% Test:
v1=1/lambda*(Id-D)*cc';
v2=A*c_state';
v3=C';


