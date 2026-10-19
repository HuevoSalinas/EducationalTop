function [CM,V] = GetElementCM2(NODE, ELEM, MATprop, Nfunc, dNfunc, intOrder)

% Initial settings
[GP_N,GP_w,GP_xi] = GetGaussQuad(intOrder);
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
        N = Nfunc(GP_xi(:,j));
        dN = dNfunc(GP_xi(:,j));
        xhat = N * XYZ;
        J = dN * XYZ;
        detJ = det(J);
        V(i, :) = V(i, :) + GP_w(j) * detJ * MATprop.t(i);
        CM(i, :) = CM(i, :)  + GP_w(j) * (xhat)  * detJ * MATprop.t(i);
    end
end
CM = CM ./ repmat(V,1,size(NODE,2));

return 