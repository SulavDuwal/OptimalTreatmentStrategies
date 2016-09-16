function G= MakeGenerator(k_HtoM_max, k_MtoL, k_MtoH, k_LtoM, eta1, eta2, mu_H1, mu_M1, mu_H2, mu_M2, NrOfStatesPerStrain,NrOfStrains,SaturationParam,PatientDeathRates)
%Take care with the condition of G. i.e. machine precision!!


%deathrate_healthy = 1/((life_expectancy_healthy-age_at_detection)*365);% [1/days]

%                    l -> death     m -> death      h -> death
%PatientDeathRates =  [deathrate_healthy log(2)/(residual_life_expectancy_StateM*365) log(2)/(residual_life_expectancy_StateH*365)];% not in publications/analysis from data sources

%                    l -> m     m -> h     
GrowthRates =       [k_LtoM     k_MtoH];

GrowthRateMatrix = zeros(2,NrOfStrains);
for i =1:NrOfStrains
    if i == 4 % wt -> all treatments efficient
        GrowthRateMatrix(1,i)= (1-eta1)*(1-eta2)*GrowthRates(1);
        GrowthRateMatrix(2,i)= (1-eta1)*(1-eta2)*GrowthRates(2);
    elseif i == 2  % m2 -> treatment#1 efficient; treatment#2 inefficient
        GrowthRateMatrix(1,i)= (1-eta1)*GrowthRates(1);
        GrowthRateMatrix(2,i)= (1-eta1)*GrowthRates(2);
    elseif i == 3  % m1 ->  treatment#1 inefficient; treatment#2 efficient
        GrowthRateMatrix(1,i)= (1-eta2)*GrowthRates(1);
        GrowthRateMatrix(2,i)= (1-eta2)*GrowthRates(2);
    else % HR -> all treatments inefficient
        GrowthRateMatrix(1,i)= GrowthRates(1);
        GrowthRateMatrix(2,i)= GrowthRates(2);
    end
end
%                    l -> 0    m -> l     h -> m    
DeclineRates =      [0         k_MtoL     k_HtoM_max]; 

DeclineRateMatrix = zeros(3,NrOfStrains);
for i =1:NrOfStrains
    DeclineRateMatrix(2,i)=DeclineRates(2);
    if i == 4 % wt -> all treatments efficient
        DeclineRateMatrix(3,i)=DeclineRates(3)/(SaturationParam*(1-eta1)*(1-eta2)+1); % DeclineRates affected by treatment 
    elseif i == 2  % m2 -> treatment#1 efficient; treatment#2 inefficient
        DeclineRateMatrix(3,i)=DeclineRates(3)/(SaturationParam*(1-eta1)+1); % DeclineRates affected by treatment 
    elseif i == 3  % m1 -> treatment#1 inefficient; treatment#2 efficient
        DeclineRateMatrix(3,i)=DeclineRates(3)/(SaturationParam*(1-eta2)+1); % DeclineRates affected by treatment 
    else % HR -> all treatments inefficient
        DeclineRateMatrix(3,i)=DeclineRates(3)/(SaturationParam+1); % DeclineRates affected by treatment 
    end
end

MutationRateMatrixH = zeros(NrOfStrains,NrOfStrains);
MutationRateMatrixM = zeros(NrOfStrains,NrOfStrains);
% From source to destination. i.e. MutationRateMatrixH(1,2) = Mutation from
% HR to M2
for i=1:NrOfStrains
    if i == 4 % wt -> all treatments efficient
        MutationRateMatrixH(i,2)= mu_H2*(1-eta1)*(1-eta2); % wt -> M2
        MutationRateMatrixM(i,2)= mu_M2*(1-eta1)*(1-eta2); % wt -> M2
        
        MutationRateMatrixH(i,3)= mu_H1*(1-eta1)*(1-eta2); % wt -> M1
        MutationRateMatrixM(i,3)= mu_M1*(1-eta1)*(1-eta2); % wt -> M1

    elseif i == 2  % m2 -> treatment#1 efficient; treatment#2 inefficient
        MutationRateMatrixH(i,4)= mu_H2*(1-eta1); % M2 -> wt
        MutationRateMatrixM(i,4)= mu_M2*(1-eta1); % M2 -> wt
        
        MutationRateMatrixH(i,1)= mu_H1*(1-eta1); % M2 -> HR (HR = M2 + M1)
        MutationRateMatrixM(i,1)= mu_M1*(1-eta1); % M2 -> HR
        
    elseif i == 3  % m1 -> treatment#1 inefficient; treatment#2 efficient
        MutationRateMatrixH(i,4)= mu_H1*(1-eta2); % M1 -> wt
        MutationRateMatrixM(i,4)= mu_M1*(1-eta2); % M1 -> wt
        
        MutationRateMatrixH(i,1)= mu_H2*(1-eta2); % M1 -> HR (HR = M2 + M1)
        MutationRateMatrixM(i,1)= mu_M2*(1-eta2); % M1 -> HR
    else % HR -> all treatments inefficient
        MutationRateMatrixH(i,3)= mu_H2; % HR -> M1
        MutationRateMatrixM(i,3)= mu_M2; % HR -> M1
        
        MutationRateMatrixH(i,2)= mu_H1; % HR -> M2 (HR = M2 + M1)
        MutationRateMatrixM(i,2)= mu_M1; % HR -> M2
    end
end


%                      HR   m2      m1      wt
MutationalNeighbours =[2     4       4       2;
                       3     1       1       3];

G = sparse(NrOfStatesPerStrain^NrOfStrains,NrOfStatesPerStrain^NrOfStrains);

% index of state = hr*4^3 + m2 *4^2 + m1*4^1 + wt*4^0 + 1 ; counting from 1...4^4
% states are counted the following way:
% [0 0 0 0] = death of patient; StateIndex 1
% [0 0 0 1] = StateIndex 2; hr = 0, m2 = 0, m1 = 0, wt = 1;
%  :    :
% [0 0 0 3] = [HR = 0; M2 = 0; M1 = 0; S = h]; StateIndex 4; hr = 0, m2 = 0, m1 = 0, wt = 3;
% [0 0 1 0] = [HR = 0; M2 = 0; M1 = ell; S = 0]; StateIndex 5;  hr = 0, m2 = 0, m1 = 1, wt = 0;

for hr = 0: NrOfStatesPerStrain-1 % highly resistant
    for m2 = 0:NrOfStatesPerStrain-1 % mutant 2
        for m1 = 0: NrOfStatesPerStrain-1 % mutant1
            for wt = 0: NrOfStatesPerStrain-1 % wild type (susceptible)
                CurrentState = [hr m2 m1 wt];
                %CurrentState = [m2 m1 wt];
                %CurrentStateIndex =  CurrentState(1)*NrOfStatesPerStrain^2 + CurrentState(2)*NrOfStatesPerStrain^1 + CurrentState(3)*NrOfStatesPerStrain^0 + 1; 
                CurrentStateIndex = CurrentState(1)*NrOfStatesPerStrain^3 + CurrentState(2)*NrOfStatesPerStrain^2 + CurrentState(3)*NrOfStatesPerStrain^1 + CurrentState(4)*NrOfStatesPerStrain^0 + 1; 
                
                %virus growth
                RelevantStrains = find(CurrentState > 0 & CurrentState < 3);
                if ~(isempty(RelevantStrains))
                    for i = 1:length(RelevantStrains);
                        K = zeros(1,NrOfStrains);
                        K(RelevantStrains(i))= 1;
                        %DestinationStateIndex = (K(1) + m2) *NrOfStatesPerStrain^2 + (K(2) + m1)*NrOfStatesPerStrain^1 + (K(3) + wt)*NrOfStatesPerStrain^0 + 1;
                        DestinationStateIndex = (K(1) + hr)*NrOfStatesPerStrain^3 + (K(2) + m2) *NrOfStatesPerStrain^2 + (K(3) + m1)*NrOfStatesPerStrain^1 + (K(4) + wt)*NrOfStatesPerStrain^0 + 1; 
                        %G(CurrentStateIndex,DestinationStateIndex) = G(CurrentStateIndex,DestinationStateIndex) + GrowthRateMatrix(CurrentState(RelevantStrains(i)),RelevantStrains(i));
                        G(CurrentStateIndex,DestinationStateIndex) = GrowthRateMatrix(CurrentState(RelevantStrains(i)),RelevantStrains(i));
                    end
                end
                %virus decline
                RelevantStrains = find(CurrentState > 1);
                for i = 1:length(RelevantStrains)
                    K = zeros(1,NrOfStrains);
                    K(RelevantStrains(i))= -1;
                    DestinationStateIndex = (K(1) + hr)*NrOfStatesPerStrain^3 + (K(2) + m2) *NrOfStatesPerStrain^2 + (K(3) + m1)*NrOfStatesPerStrain^1 + (K(4) + wt)*NrOfStatesPerStrain^0 + 1; 
                    %DestinationStateIndex = (K(1) + m2) *NrOfStatesPerStrain^2 + (K(2) + m1)*NrOfStatesPerStrain^1 + (K(3) + wt)*NrOfStatesPerStrain^0 + 1; 
                    %G(CurrentStateIndex,DestinationStateIndex) = G(CurrentStateIndex,DestinationStateIndex) + DeclineRateMatrix(CurrentState(RelevantStrains(i)),RelevantStrains(i));
                    G(CurrentStateIndex,DestinationStateIndex) = DeclineRateMatrix(CurrentState(RelevantStrains(i)),RelevantStrains(i));
                end
                
                %virus mutation from 'H' States
                RelevantStrains = find(CurrentState == 3); % mutation only from 'h' states
                PossibleDestination = find(CurrentState == 0); % makes only sense into '0' states
                if ~isempty(PossibleDestination) && ~isempty(RelevantStrains) 
                    for i = 1:length(RelevantStrains)
                        Neighbours = MutationalNeighbours(:,RelevantStrains(i));
                        c = intersect(PossibleDestination,Neighbours');
                        if ~isempty(c)
                            for j = 1:length(c)
                                K = CurrentState;
                                K(c(j)) = 1;
                                DestinationStateIndex = K(1)*NrOfStatesPerStrain^3 + K(2) *NrOfStatesPerStrain^2 + K(3)*NrOfStatesPerStrain^1 + K(4)*NrOfStatesPerStrain^0 + 1; 
                                %DestinationStateIndex =  K(1) *NrOfStatesPerStrain^2 + K(2)*NrOfStatesPerStrain^1 + K(3)*NrOfStatesPerStrain^0 + 1; 
                                %if DestinationStateIndex ~= CurrentStateIndex
                                %G(CurrentStateIndex,DestinationStateIndex) = G(CurrentStateIndex,DestinationStateIndex) + MutationRateMatrixH(RelevantStrains(i));
                                G(CurrentStateIndex,DestinationStateIndex) = MutationRateMatrixH(RelevantStrains(i),c(j));
                                %end
                            end
                        end
                    end
                end
                
                %virus mutation from 'M' States
                RelevantStrains = find(CurrentState == 2); % mutation only from 'h' states
                PossibleDestination = find(CurrentState == 0); % makes only sense into '0' states
                if ~isempty(PossibleDestination) && ~isempty(RelevantStrains) 
                    for i = 1:length(RelevantStrains)
                        Neighbours = MutationalNeighbours(:,RelevantStrains(i));
                        c = intersect(PossibleDestination,Neighbours');
                        if ~isempty(c)
                            for j = 1:length(c)
                                K = CurrentState;
                                K(c(j)) = 1;
                                DestinationStateIndex = K(1)*NrOfStatesPerStrain^3 + K(2) *NrOfStatesPerStrain^2 + K(3)*NrOfStatesPerStrain^1 + K(4)*NrOfStatesPerStrain^0 + 1; 
                                %DestinationStateIndex =  K(1) *NrOfStatesPerStrain^2 + K(2)*NrOfStatesPerStrain^1 + K(3)*NrOfStatesPerStrain^0 + 1; 
                                %if DestinationStateIndex ~= CurrentStateIndex
                                %G(CurrentStateIndex,DestinationStateIndex) = G(CurrentStateIndex,DestinationStateIndex) + MutationRateMatrixM(RelevantStrains(i));
                                G(CurrentStateIndex,DestinationStateIndex) = MutationRateMatrixM(RelevantStrains(i),c(j));
                                %end
                            end
                        end
                    end
                end                

                %Patient death
                deathrateIndex =  max([hr m2 m1 wt]);
                if deathrateIndex
                    G(CurrentStateIndex,1) = PatientDeathRates(deathrateIndex);
                end
            end
        end
    end
end


%%SET Diagonal element - row sum

for i = 1:NrOfStatesPerStrain^NrOfStrains
    G(i,i) = -sum(G(i,:));
end
