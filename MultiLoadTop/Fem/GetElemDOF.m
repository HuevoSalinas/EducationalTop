function [eDOF] = GetElemDOF(ELEM,DpN)

% Initial parameters
Ne = size(ELEM,1); % Number of elments
NpE = size(ELEM,2); % Nodes per element

% Degrees of freedom 
eDOF = nan(Ne,DpN*NpE);
for i=1:DpN
    eDOF(:,i:DpN:end) = (DpN*ELEM) - (DpN-i);
end 

return