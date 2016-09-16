function [c] = MakeTerminalCostfunction4(c_states,t_start,t_end,lambda)
c = c_states'*(exp(-lambda*t_start)-exp(-lambda*t_end))/lambda;
