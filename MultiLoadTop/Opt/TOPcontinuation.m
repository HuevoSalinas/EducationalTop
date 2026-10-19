function [opt] = TOPcontinuation(fem,opt,iter)
narginchk(2,3)

rebuildH = false;
maxIter = opt.contin(end,1);

if (nargin==2 || iter==1)
    rebuildH = true;
    opt.current_penal = opt.contin(1,2);
    opt.current_move = opt.contin(1,3);
    opt.current_rmin = opt.contin(1,4);
    if opt.contin(1,1)~=0, fprintf('WARNING - continuation matrix should start at step 0.\n'), end
    fprintf('STEP %g:\tpenal = %g\tmove = %g\trmin = %g\n',...
        1,opt.current_penal,opt.current_move,opt.current_rmin)
    
elseif ((iter<maxIter) && any(iter==opt.contin(:,1)))
    ind = find(iter==opt.contin(:,1),1,'last');
    opt.current_penal = opt.contin(ind,2);
    opt.current_move = opt.contin(ind,3);
    if ~isequal(opt.contin(ind-1,4),opt.contin(ind,4))
        rebuildH = true;
        opt.current_rmin = opt.contin(ind,4);
    end

    fprintf('STEP %g:\tpenal = %g\tmove = %g\trmin = %g\n',...
        ind,opt.current_penal,opt.current_move,opt.current_rmin)
end

% H computed using trees.
if rebuildH
    fprintf('--- Updating the filter matrix: rmin = %g and exponent s = %g\n',opt.current_rmin,opt.filterS);
    opt.H = GetFilterH(fem.CM,opt.current_rmin,opt.filterS,opt.active);
end
return