function [u,V,allu,allV,B,P,J,X] = runBruteForceF(A,P0,c,k,cs_integral,ct_integral,treatments,states)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose
% Finds the optimal control history for cost function with running and 
% terminal cost :
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% close all graphics 
% close all;
% clear memory 
% clear all
%##########################################################################
% *Display controls*
%%
% Initialization of struct containing
X{1} = P0; % States
L{1} = 0; % Running cost (L = short for Lagrange type)
cumL{1} = 0; % Cummulative running cost 
% Calculating states and running cost:
for i = 1:k % Timesteps
    %% All possible child nodes from given parent nodes
    X{i+1} = NaN(states,treatments^i);
    L{i+1} = NaN(1,treatments^i);
    cumL{i+1} = NaN(1,treatments^i);
    for j = 1:treatments % treatments
        if j == 1
            X{i+1}(1:states,1:2^(i-1)) = A{j}*X{i};
            L{i+1}(1,1:2^(i-1)) = cs_integral(:,:,i,j)'*X{i} + ct_integral(:,:,i,j)'*X{i};
            cumL{i+1}(1,1:2^(i-1)) = L{i+1}(1,1:2^(i-1)) + cumL{i};
        else
            X{i+1}(1:states,2^(i-1)+1:2^i) = A{j}*X{i}; 
            L{i+1}(1,2^(i-1)+1:2^i) = cs_integral(:,:,i,j)'*X{i} + ct_integral(:,:,i,j)'*X{i};
            cumL{i+1}(1,2^(i-1)+1:2^i) = L{i+1}(1,2^(i-1)+1:2^i) + cumL{i};
        end
    end
end
% Terminal cost (M = short for oMega or Mayer type)
M = c'*X{k+1}; 
% Total cost = M + L (B = short for Bolza type)
B = M + cumL{k+1};
%% Backtracing
allu = NaN(length(B),k); % control history
tracebackIndices = NaN(length(B),k);
for j = 1:length(B)
    li = j;
    u = NaN(1,k); 
    for i = k:-1:1
        tracebackIndices(j,i) = li;
        temp1 = treatments^(i-1);
        temp2 = li/temp1;
        if temp2 > 1
            temp3 = 2;
        else
            temp3 = 1;
        end
        u(i) = temp3;
        li = mod(li,temp1);
        if li == 0
            li = temp1;
        end
    end
    allu(j,:) = u;
end
%% Value Function plot
allV = NaN(length(B),k+1);
for j = 1:length(B)
    allV(j,k+1) = M(tracebackIndices(j,k));
    for i = k:-1:1
        allV(j,i) = allV(j,i+1) + L{i+1}(tracebackIndices(j,i));
    end
end
%% Backtracing
u = NaN(1,k); % control history
tracebackIndices = NaN(1,k);
[~,li] = min(B);
for i = k:-1:1
    tracebackIndices(i) = li;
    temp1 = treatments^(i-1);
    temp2 = li/temp1;
    if temp2 > 1
        temp3 = 2;
    else
        temp3 = 1;
    end
    u(i) = temp3;
    li = mod(li,temp1);
    if li == 0
        li = temp1;
    end
end
%% Value Function plot
% V = NaN(1,k+1);
% V(k+1) = M(tracebackIndices(k));
P = NaN(states,k+1);
P(:,1) = P0;
for i = k:-1:1
%     V(i) = V(i+1) + L{i+1}(tracebackIndices(i));
    P(:,i+1) = X{i+1}(:,tracebackIndices(i)); 
end
%% n.b 11-10-2012
[~,id] = min(allV(:,1));
V = allV(id,:);
u = allu(id,:);
J = NaN(1,k+1);
tempVar = 0;
for i = 1:k
    tempVar = tempVar + cs_integral(:,:,i,u(i))'*P(:,i) + ct_integral(:,:,i,u(i))'*P(:,i);
    J(i) = tempVar;
end
J(k+1) = J(k) + c'*P(:,k+1);
% c_array(:,1) = c;
% for i = 1:k
%     c_array(:,i+1) = c_array(:,i)'*A{u(k-i+1)}+q{u(k-i+1)}';
%     tempV(i) = c_array(:,i)'*P(:,k-i+2);
% end
% eof