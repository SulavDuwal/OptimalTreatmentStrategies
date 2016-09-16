function [V,action,t_opt] = DisInfoPolicyFunction(L,c,k_info,lambda,t_stepsize,N)

% L(:,:,a) = Generator for action a 
% c(:,a) = cost function for action a (costs per unit of time)
% k_info = information costs (absolute value)
% lambda: exp(-lambda) = discount factor 
% t_stepsize = stepsize for the time
% N = maximal number of timesteps (N*t_stepsize = maximal time)

size_S = size(L,1); % number of states
size_A = size(L,3); % number of actions

Id = eye(size_S);

%% Pre-initialization
t = NaN(1,N);
C = NaN(size_S,N,size_A);
expL = NaN(size_S,size_S,size_A);
action = NaN(1,size_S);
ind_t  = NaN(1,size_S);
K_info = NaN(1,size_S);

% define new cost function:
for n = 1:N
    t(n) = n*t_stepsize;
    for a = 1:size_A
        C(:,n,a) = (L(:,:,a)-lambda*Id)^(-1)*(expm((L(:,:,a)-lambda*Id)*t(n))-Id)*c(:,a);
    end
end

for a = 1:size_A
    expL(:,:,a) = expm(L(:,:,a));
end

% initialise strategy:
for s = 1:size_S
    action(s) = 1;
    ind_t(s)  = N;
    K_info(s) = k_info;
end

test = 0;
l = 0;

while test < size_S+1 && l < 100
    dis = zeros(size_S,size_S);
    c_strat = NaN(1,size_S);
    for s = 1:size_S
        dis(s,s) = exp(-lambda*t(ind_t(s)));
        c_strat(s) = C(s,ind_t(s),action(s));
        Q = expL(:,:,action(s))^t(ind_t(s));
        P(s,:) = Q(s,:);
    end
    V = inv(Id-dis*P)*(c_strat'+dis*K_info'); 
    D = NaN(size_S,N,size_A);
    for a = 1:size_A
        for n = 1:N
            D(:,n,a) = C(:,n,a)+exp(-lambda*t(n))*(K_info'+expL(:,:,a)^t(n)*V);
        end
    end
    test=1;
    for s = 1:size_S
        m = NaN(1,size_A);
        num = NaN(1,size_A);
        for a = 1:size_A
            [mini,nmini] = min(D(s,:,a));
            m(a) = mini;
            num(a) = nmini;
        end
        [~,nmini] = min(m);
        if nmini == action(s) && num(nmini) == ind_t(s)
            test = test+1;
        else
            action(s) = nmini;
            ind_t(s) = num(nmini);
        end
    end
    l = l+1;
end


t_opt = ind_t*t_stepsize;   
t_opt = t_opt';
action = action';
    