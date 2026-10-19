function [opt, w] = KSWeightUpdate(opt)

eta = opt.norm;
if ~isfield(opt,'J0')
    opt.J0 = max(opt.MJ);
end
    
fKS = (1/eta) * log(sum(exp(eta * ((opt.MJ)/opt.J0 -1))));
w = exp(eta * (opt.MJ/opt.J0-1) )'/( opt.J0 * exp(eta * fKS));
opt.J = fKS;
end
