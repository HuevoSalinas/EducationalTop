function [x,opt] = ShapeTOPcontinuation(x,fem,opt,iter)
narginchk(3,4)

rebuildH = false;
remap = false;
maxIter = opt.contin(end,1);
update = false;

if (nargin==3 || iter==1)
    rebuildH = true;
    ind = 1;
    update = true;
    if opt.contin(1,1)~=0, fprintf('WARNING - continuation matrix should start at step 0.\n'), end 
elseif ((iter<maxIter) && any(iter==opt.contin(:,1)))
    ind = find(iter==opt.contin(:,1),1,'last');
    update = true;
    if ~isequal(opt.contin(ind-1,3),opt.contin(ind,3))
        remap = true;   
    end
    if ~isequal(opt.contin(ind-1,5),opt.contin(ind,5))
        rebuildH = true;
    end
end

if update
    opt.current_penal = opt.contin(ind,2);
    opt.current_mu = opt.contin(ind,3);
    opt.current_move = opt.contin(ind,4);
    opt.current_rmin = opt.contin(ind,5);
    fprintf('STEP %g:\tpenal = %g\tmu = %g\tmove = %g\trmin = %g\n',...
        ind,opt.current_penal,opt.current_mu,opt.current_move,opt.current_rmin)
end

if rebuildH
    fprintf('--- Updating the filter matrix: rmin = %g and exponent s = %g\n',...
            opt.current_rmin,opt.filterS);
    
    opt.H = GetFilterH(fem.CM,opt.current_rmin,opt.filterS,opt.active);
end

if remap 
    x = SigmoidRemap(x, opt.contin(ind-1,3), opt.contin(ind,3));
    x(opt.pass_solid) = 1;
    x(opt.pass_void) = opt.xmin;
end

end

function [x] = SigmoidRemap(x, mu_old, mu_new)

eps = 1e-5;
n = 0.5;
x_old = 0.5 + mu_old;
x_new = 0.5 + mu_new;

func = @(x,penal)((tanh(penal*n)+tanh(penal*(x-n)))./(tanh(penal*n)+tanh(penal*(1-n))));
func_p = @(penal)(func(x_old, penal) - x_new - eps);
p = fzero(func_p, 1);
x = func(x,p);

end