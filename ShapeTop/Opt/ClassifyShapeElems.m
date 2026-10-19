function [opt] = ClassifyShapeElems(x,opt)

if opt.current_mu>=0.5, error('Shape optimization mu must be positive and <0.5.'), end
if ~isfield(opt,'shapeIsBorderTrans'), opt.shapeIsBorderTrans = true; end

% Classify the elements in the void, solid, and neighboring these two phases.
flagVoid  = double(x <= 0.5 - opt.current_mu);
flagSolid = double(x >= 0.5 + opt.current_mu);
flagVoidNeigh  = logical(opt.L*flagVoid);
flagSolidNeigh = logical(opt.L*flagSolid);

% Classify void elements surrounded by void: pure void elems
opt.flagPureVoid  = and(flagVoid,and(flagVoidNeigh,not(flagSolidNeigh)));

% Classify solid elements surrounded by solid: pure solid elems
opt.flagPureSolid = and(flagSolid,and(flagSolidNeigh,not(flagVoidNeigh)));
if opt.shapeIsBorderTrans
    % Do not consider pure solid those that are in the domain's border
    % Classify elements in the boundary
    flagBorder = false(numel(x),1);
    flagBorder(opt.elemBORDER) = true;
    opt.flagPureSolid = and(opt.flagPureSolid,not(flagBorder));
end

end