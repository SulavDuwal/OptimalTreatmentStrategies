function [c]=MakeCostfunction(GDP,LossOfProductivityH,LossOfProductivityM)

% Create a cost function c

CostGDPH=LossOfProductivityH*GDP/365; % Cost per day because of loss of productivity
CostGDPM=LossOfProductivityM*GDP/365; % Cost per day because of loss of productivity


NrOfStatesPerStrain = 4;


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
                CurrentStateIndex = CurrentState(1)*4^3 + CurrentState(2)*4^2 + CurrentState(3)*4^1 + CurrentState(4)*4^0 + 1; 
                if any(CurrentState>=3)==1
                    c(CurrentStateIndex)=CostGDPH;
                elseif any(CurrentState>=2)==1
                    c(CurrentStateIndex)=CostGDPM;
                else
                    c(CurrentStateIndex)=0;
                end                
            end
        end
    end
end

c(1)=GDP/365;