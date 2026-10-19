function [opt, w] = PNormWeightUpdate(opt)

p = opt.norm;
pnorm = norm(opt.MJ, p);
w = (opt.MJ/ pnorm)'.^(p-1);
opt.J = pnorm;

end
