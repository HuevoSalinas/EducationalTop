function [CM,V] = GetElementCM2(NODE, ELEM, MATprop, Nfunc, dNfunc, intOrder)

% Initial settings
[GP_N,GP_w,GP_xi] = GetGaussQuad(intOrder);
N_GP = GetShapeData(Nfunc,GP_xi);
dN_GP = GetShapeData(dNfunc,GP_xi);
Ne = size(ELEM,1);
if isempty(MATprop)
    MATprop = struct('t',ones(Ne,1));
end

% Volume and center of mass integration
CM = zeros(Ne,size(NODE,2));
V = zeros(Ne,1);
for i=1:Ne
    XYZ = NODE(ELEM(i,:),:);
    for j=1:GP_N
        N = N_GP(:,:,j);
        dN = dN_GP(:,:,j);
        xhat = N * XYZ;
        J = dN * XYZ;
        detJ = det(J);
        w = GP_w(j)*detJ*MATprop.t(i);
        V(i, :) = V(i, :) + w;
        CM(i, :) = CM(i, :)  + w*xhat;
    end
end
CM = CM ./ repmat(V,1,size(NODE,2));

end