function [opt, w] = ArbitraryWeightUpdate(opt)

w = opt.w;
opt.J = opt.MJ * w;
return
