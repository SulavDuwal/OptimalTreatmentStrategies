function [u,V,Q,J,P,SJ,TrJ,Resi,Q_size,minJArray,maxJArray,fvalArray,ThresholdArray,exitflags] = runBDA_Parallel(A,P0,c,k,cs_integral,ct_integral,treatments,states,A_max_series,A_min_series,maxJ,P11,P22,J11,J22,DeathProb,h_R1_index,h_R1,h_R2_index,h_R2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author  : Sulav
% Purpose : Dynamic Programming Algorithms for Bolza Type of cost
% return - Control history (Optimal Switching Sequence) and other details 
% 1) BDA-III parallel *BDA-III parallel*
% 2) With bounds on optimal cost and state vectors
% 3) For LP - default package Linprog !
% 4) Parallel execution of redundancy tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% *Backward solution*
%% Setting option for linear programming solver
% Uncomment you use cplex
% options = cplexoptimset('Display','off'); % if cplex is used
options = optimset('Display','off','LargeScale','on'); % if linprog of matlab is used
% Initializations
Q = c;    % time zero
% For tracing back optimal sequence
const_size = 20;
sizeTraceback = floor((k-1)/const_size)+1;
traceback = zeros(sizeTraceback,1);
% For storage
exitflags = [];
Resi = [];
Q_size = NaN(1,k);
minJArray = NaN(1,k);
maxJArray = NaN(1,k);
fvalArray = [];
ThresholdArray = [];
%R_min = 0;
residual_vector = [ones(states,1); 0; 0];
h_R1_vector = zeros(states+2,1);
h_R1_vector(h_R1_index) = 1;
h_R2_vector = zeros(states+2,1);
h_R2_vector(h_R2_index) = 1;
% Relaxation factor
lax_factor = 1e-3;
maxJ = maxJ + lax_factor*maxJ;
for i = 1:k % timesteps
    %display(i);
    %% All possible child nodes from given parent nodes
    numPreviousLeaves = size(Q,2);
    numLeaves = treatments*numPreviousLeaves;
    Q_new = zeros(states,numLeaves);
    T_new = zeros(states+2,numLeaves);
    traceback_new = repmat(traceback,1,treatments);
    idx = floor((i-1)/const_size)+1;
    modi = mod(i-1,const_size);
    for j = 1:treatments % treatments
        if j == 1
            interval = 1:numPreviousLeaves;
        else
            interval = (j-1)*numPreviousLeaves+(1:numPreviousLeaves);
        end
         Q_new(:,interval) = A{j}'*Q + ...   % Backward Propagation of previous cost vector 
                                 repmat(ct_integral(:,:,k-i+1,j),1,numPreviousLeaves) + ... % State cost vector
                                 repmat(cs_integral(:,:,k-i+1,j),1,numPreviousLeaves);      % Treatment cost vector
         traceback_new(idx,interval) = treatments^modi*(j-1) + traceback_new(idx,interval); 
    end
    % c'x + R - V* >= 0
    % (c' 1 -1)*(x R V*)
    T_new(1:states,:) = Q_new; %  c  
    T_new(states+1,:) = 1;     %  1  
    T_new(states+2,:) = -1;    % -1
    %% Bounds for state vectors
    if i < k
        Pmax = A_max_series{k-i}*P0;
        Pmax = Pmax + lax_factor*Pmax;
        Pmax = min(Pmax,1);
        Pmin = A_min_series{k-i}*P0;
        Pmin = Pmin - lax_factor*Pmin;
        Pmin(Pmin < 0) = 0;
        %R_min = max(min(C_st_min(:,k-i+1)'*Pmin),R_min); 
    else 
        Pmax = P0+lax_factor*P0;
        Pmin = P0-lax_factor*P0;
        flag = Pmin < 0;
        Pmin(flag) = 0;
        %R_min = max(min(C_st_min(:,k-i+1)'*Pmin),R_min);
    end
    Pmax(1) = DeathProb(k-i+1);
    %% Dynamically estimation of lowerbound of J*
    minJ = min(Q_new'*Pmin);
    minJ = minJ - lax_factor * minJ;
    minJArray(i) = minJ;
    lowerbound = [Pmin;0;minJ];
    upperbound = [Pmax;maxJ;maxJ];
    %% Trimming redundant nodes (using LP) - parallel execution of tests !!
    flag_nonredundant = true(1,numLeaves);
    exitflag_arr = NaN(1,numLeaves);
    h_R1 = sum(P11(h_R1_index,i+1));
    h_R2 = sum(P22(h_R2_index,i+1));
    parfor par_i = 1:numLeaves
        indices = true(1,numLeaves);
        indices(1,par_i) = false;
        if ~isempty(indices)
            %% Redundancy Test
            % Objective function LP
            f = T_new(:,par_i);
            % Constraints for LP
            AA = vertcat((-1*T_new(:,indices))',h_R1_vector',h_R2_vector'); 
            b = vertcat(zeros(sum(indices),1),h_R1,h_R2); 
            %% Solving linear programming
            % Uncomment this line if you are using CPLEX
            % [~,fval,exitflag] = cplexlp(f,AA,b,residual_vector',1,lowerbound,upperbound,[],options);
            [~,fval,exitflag] = linprog(f,AA,b,residual_vector',1,lowerbound,upperbound,[],options);
           
            exitflag_arr(1,par_i) = exitflag;
            if fval > 0
                % If the column (vector) of the matrix, is recognized as 
                % redundant, then delete column and update matrix
                flag_nonredundant(1,par_i) = false;
            else
                % The optimal co-state trajectory may pass through column
            end
        end
    end
    exitflags = [exitflags exitflag_arr];
    T_new = T_new(:,flag_nonredundant);
    Q = Q_new(:,flag_nonredundant);
    traceback = traceback_new(:,flag_nonredundant);
    sizeCurQ = size(Q,2);
    %display(sizeCurQ);
    display(sprintf('Iteration Nr %3d - Cand Size %4d',i,sizeCurQ));
    Q_size(i) = sizeCurQ;
    if sizeCurQ > 25000
        error('Nr. Candidates exceeded 25e3')
    end
    %% Dynamically estimating upper bound of J*
    maxJArray(i) = maxJ;
    tempOneJ = min(Q'*P11(:,k-i+1) + J11(k-i+1));
    tempOneJ = tempOneJ + lax_factor*tempOneJ;
    tempTwoJ = min(Q'*P22(:,k-i+1) + J22(k-i+1));
    tempTwoJ = tempTwoJ + lax_factor*tempTwoJ;
    maxJ = min([tempOneJ,tempTwoJ,maxJ]);
%     display(length(indices));
    clear Q_temp t_temp traceback_new
end
%% *Traceback 2*
P = zeros(states,k+1);
P(:,1) = P0; %initialize recursion

% Optimal cost and index of optimal sequence
[V,Index] = min(Q'*P(:,1));
u_str = '';
for z = 1:(sizeTraceback)
    path_index = traceback(z,Index);
    str_size = const_size;
    if z == sizeTraceback
        str_size = mod(k,const_size);
        if str_size == 0
           str_size = const_size; 
        end
    end
    u_str = [dec2base(path_index,treatments,str_size) u_str];
end
% Optimal switching sequences :
u = NaN(1,length(u_str));
for i = 1:length(u_str)
    u(i) = str2double(u_str(i)) + 1;
end
% Trajectory of optimal state trajectory - P :
for i = 1:k
    P(:,i+1) = A{u(i)}*P(:,i); 
end

% Cost splitting :
SJ     = zeros(1,k+1); % State cost
TrJ    = zeros(1,k+1); % Treatment cost
ST     = zeros(1,k+1);   % State + Treatment cost
for i = 1:k
    SJ(i)  = cs_integral(:,:,i,u(i))'*P(:,i);
    TrJ(i) = ct_integral(:,:,i,u(i))'*P(:,i);
    ST(i)  = SJ(i)    + TrJ(i); 
end
J = cumsum(ST); % Cummulative sum of running cost
J(end) = J(end-1) + c'*P(:,k+1); % Adding terminal cost
% display('Backward')
% eof