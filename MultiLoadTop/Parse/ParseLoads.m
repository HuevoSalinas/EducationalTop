function [f] = ParseLoads(LOAD,Nn, DpN)

Ndof = DpN * Nn;

if ~isempty(LOAD)
    if max(LOAD(:,1)) > Nn
        error('Node number > Nn is being loaded')
    elseif min(LOAD(:,1)) < 1
        error('Invalid node number loaded')
    elseif any(isnan(LOAD(:,1)))
        error('Node number cannot be NaN')
    end
    
    [i,j] = find(~isnan(LOAD(:,2:end)));
    val = LOAD(sub2ind(size(LOAD),i,j+1));
    if isrow(j), j = j'; end
    myDOF = DpN*(LOAD(i,1)-1) + j;
    f = accumarray(myDOF,val,[Ndof 1]);
else
    f = zeros(Ndof, 1);
end

return