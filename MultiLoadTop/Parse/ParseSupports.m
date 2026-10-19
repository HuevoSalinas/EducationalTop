function [suppDOF,suppVAL] = ParseSupports(SUPP,Nn, DpN)

if ~isempty(SUPP)
    if max(SUPP(:,1)) > Nn
        error('Node number > Nn in the supports list')
    elseif min(SUPP(:,1)) < 1
        error('Invalid node number in the supports list')
    elseif any(isnan(SUPP(:,1)))
        error('Node number cannot be NaN')
    end
    
    [i,j] = find(~isnan(SUPP(:,2:end)));
    suppVAL = SUPP(sub2ind(size(SUPP),i,j+1));
    if isrow(j), j = j'; suppVAL = suppVAL'; end
    suppDOF = DpN*(SUPP(i,1)-1) + j;
    
    [suppDOF,order] = sort(suppDOF);
    suppVAL = suppVAL(order);
else
    suppDOF = [];
    suppVAL = [];
end 