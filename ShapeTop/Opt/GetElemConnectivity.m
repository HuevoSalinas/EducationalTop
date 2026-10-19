function [opt] = GetElemConnectivity(ELEM,opt)

if ~isfield(opt,'shapeElemNeigh'), opt.shapeElemNeigh = 'node'; end
Ne = size(ELEM,1);
owner = repmat((1:Ne)',4,1);
edges = [ELEM(:,[1 2]); ELEM(:,[2 3]); ELEM(:,[3 4]); ELEM(:,[4 1])];
[~,~,edgeID] = unique(sort(edges,2),'rows');

switch lower(opt.shapeElemNeigh)
    case 'node'
        B = sparse(owner,ELEM(:),1,Ne,max(ELEM(:)));
    case 'edge'
        B = sparse(owner,edgeID,1,Ne,max(edgeID));
    otherwise
        error('Neighbor type must be node or edge.')
end
opt.L = B*B' > 0;
opt.L(1:Ne+1:end) = false;

edgeCount = accumarray(edgeID,1);
opt.elemBORDER = unique(owner(edgeCount(edgeID)==1));

end
