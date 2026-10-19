function [L] = ParseDisplacements(DISP,Nn,DpN)

Ndof = DpN * Nn;
if ~isempty(DISP)
    if max(DISP(:,2)) > Nn
        error('Node number > Nn is used on target displacements')
    elseif min(DISP(:,2)) < 1
        error('Invalid node number is used on target displacements')
    elseif any(isnan(DISP(:,1)))
        error('Node number cannot be NaN')
    end
    
    Nt = max(DISP(:, 1));
    L = sparse(Nt, Ndof);
    for i = 1:DpN
        L = L + sparse(DISP(:,1), DpN* DISP(:,2) -(DpN-i), DISP(:,i+2),Nt, Ndof) ;
    end
else
    L = sparse(1, Ndof);
end 

return