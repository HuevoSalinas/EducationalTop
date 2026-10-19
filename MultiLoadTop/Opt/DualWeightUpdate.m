function [opt, w] = DualWeightUpdate(opt)

No = size(opt.MJ, 2);

if ~isfield(opt,'dual')
    w = ones(No, 1)/No;
    f0 =  opt.MJ * w;
    low = [];
    upp = [];
    opt.dual = struct('w', w, 'wold1', w, 'wold2', w, 'f0', f0/1000, ...
                     'low', low, 'upp', upp, 'iter', 1);
end
iter = opt.dual.iter;
low = opt.dual.low;
upp = opt.dual.upp;
w = opt.dual.w;
opt.J = opt.MJ * w;
f0 = opt.dual.f0;
wold1 = opt.dual.wold1;
wold2 = opt.dual.wold2;
dw = opt.movew;
w_lb = max(zeros(No, 1), w - dw);
w_ub = min(ones(No, 1), w + dw);
f = -opt.J/abs(f0);
df = -opt.MJ'/abs(f0);
df2 = zeros(No, 1);
g = [sum(w) - 1; 1 - sum(w)];
dg = [ones(1, No); -ones(1, No)];
dg2 = zeros(2, No);
a0 = 1;
a = zeros(2, 1);
c = 1e3 * ones(2, 1);
d = zeros(2, 1);

warning off all
[w,~ ,~, ~, ~, ~, ~, ~, ~,low,upp] = ...
    mmasub(2, No, iter, w, w_lb, w_ub, wold1, wold2, ...
    f, df, df2, g, dg, dg2, low, upp, a0,a,c,d);
warning on all
opt.dual.wold2 = opt.dual.wold1;
opt.dual.wold1 = opt.dual.w;
opt.dual.w = w;
opt.dual.low = low;
opt.dual.upp = upp;
opt.dual.iter = iter + 1;

end
