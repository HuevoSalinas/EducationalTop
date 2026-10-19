function [MaxIter, hist_J, hist_V, hist_w, x, rho, opt, fem] = ...
    MultiLoadTOPSolve(NODE, ELEM, SUPP, LOAD, fem, opt, x, optoutput)

% Selection of the multiple load topology optimization method
switch lower(opt.combinationMethod)
    case 'pareto'
        [MaxIter, hist_J, hist_V, hist_w, x, rho, opt, fem] = ...
                ParetoTOPSolve(NODE, ELEM, SUPP, LOAD, fem, opt, x, optoutput);
    case 'pareto-search'
        [MaxIter, hist_J, hist_V, hist_w, x, rho, opt, fem] = ...
                ParetoSearchTOPSolve(NODE, ELEM, SUPP, LOAD, fem, opt, x, optoutput);
    case 'arbitrary'
        WupdateFunc = @ArbitraryWeightUpdate; % Function for weight updating
        [MaxIter, hist_J, hist_V, hist_w, x, rho, opt, fem] = ...
                VariableWeightTOPSolve(NODE, ELEM, SUPP, LOAD, fem, opt, x, ...
                                      WupdateFunc, optoutput);
    case 'p-norm'
        WupdateFunc = @PNormWeightUpdate;     % Function for weight updating
        [MaxIter, hist_J, hist_V, hist_w, x, rho, opt, fem] = ...
                VariableWeightTOPSolve(NODE, ELEM, SUPP, LOAD, fem, opt, x, ...
                                       WupdateFunc, optoutput);
    case 'ks'
        WupdateFunc = @KSWeightUpdate;        % Function for weight updating
        [MaxIter, hist_J, hist_V, hist_w, x, rho, opt, fem] = ...
                VariableWeightTOPSolve(NODE, ELEM, SUPP, LOAD, fem, opt, x, ...
                                       WupdateFunc, optoutput);
    case 'bounded-dual'
        WupdateFunc = @DualWeightUpdate;      % Function for weight updating
        [MaxIter, hist_J, hist_V, hist_w, x, rho, opt, fem] = ...
                VariableWeightTOPSolve(NODE, ELEM, SUPP, LOAD, fem, opt, x, ...
                                       WupdateFunc, optoutput);
    case 'bound-constrained'
        [MaxIter, hist_J, hist_V, x, rho, opt, fem] = ...
                BoundConstrainedTOPSolve(NODE, ELEM, SUPP, LOAD, fem, opt, x, optoutput);
        hist_w = [];
    otherwise
        error('Unknown combination method')
 
end



return 