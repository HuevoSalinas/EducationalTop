function [opt] = TOPcontinuation(fem,opt,iter)
narginchk(2,3)

rebuildH = false;
maxIter = opt.contin(end,1);
update = false;

if (nargin==2 || iter==1)
    rebuildH = true;
    ind = 1;
    update = true;
    if opt.contin(1,1)~=0, fprintf('WARNING - continuation matrix should start at step 0.\n'), end
elseif ((iter<maxIter) && any(iter==opt.contin(:,1)))
    ind = find(iter==opt.contin(:,1),1,'last');
    update = true;
    if ~isequal(opt.contin(ind-1,4),opt.contin(ind,4))
        rebuildH = true;
    end
end

if update
    opt.current_penal = opt.contin(ind,2);
    opt.current_move = opt.contin(ind,3);
    opt.current_rmin = opt.contin(ind,4);
    fprintf('STEP %g:\tpenal = %g\tmove = %g\trmin = %g\n',...
        ind,opt.current_penal,opt.current_move,opt.current_rmin)
end

if rebuildH
    fprintf('--- Updating the filter matrix: rmin = %g and exponent s = %g\n',opt.current_rmin,opt.filterS);
    opt.H = GetFilterH(fem.CM,opt.current_rmin,opt.filterS,opt.active);
end
end
