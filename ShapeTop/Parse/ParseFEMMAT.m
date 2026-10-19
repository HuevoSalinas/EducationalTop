function [mat,fem] = ParseFEMMAT(mat,fem,NODE,ELEM,SUPP,LOAD,DISP)

Nn = size(NODE,1);
Ne = size(ELEM,1);
props = {'E', 'nu', 't'};
default = {1, 0.3, 1};
for i = 1:size(props,2)
    if ~isfield(mat,props{i}), mat.(props{i}) = default{i}; end
    if isscalar(mat.(props{i})), mat.(props{i}) = repmat(mat.(props{i}),Ne,1); end
end
fem.DpN = 2;
Nfunc = @ShapeFuncQ4;
dNfunc = @ShapeDerivQ4;
intOrder = 2;
[fem.CM,fem.ve] = GetElementCM2(NODE,ELEM,mat,Nfunc,dNfunc,intOrder);
fem.f = ParseLoads(LOAD, Nn, fem.DpN);
[fem.suppDOF,fem.suppVAL] = ParseSupports(SUPP, Nn, fem.DpN);
fem.L = ParseDisplacements(DISP, Nn, fem.DpN);

end
