function [mat,fem] = ParseMultiLoadFEMMAT(mat,fem,NODE,ELEM,SUPP,LOAD,DISP)

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

[fem.suppDOF,fem.suppVAL] = ParseSupports(SUPP,Nn,fem.DpN);

f = cell(1, size(LOAD, 2));
for i = 1: size(LOAD, 2)
    f{i} = ParseLoads(LOAD{i}, Nn, fem.DpN);
end
fem.f = cell2mat(f);

L = cell(1, size(DISP, 2));
for i = 1: size(DISP, 2)
    L{i} = cell(1, size(DISP{i}, 2));
    for j = 1: size(DISP{i}, 2)
        L{i}{j} = ParseDisplacements(DISP{i}{j}, Nn, fem.DpN);
    end
end
fem.L = L;

end