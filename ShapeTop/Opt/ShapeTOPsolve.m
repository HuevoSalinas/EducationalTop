function [MaxIter, hist_J, hist_V, x, rho, opt, fem] = ...
    ShapeTOPsolve(NODE, ELEM, fem, opt, x, optoutput)

% Initial parameters
MaxIter = opt.contin(end,1);
hist_J = nan(MaxIter+1,1);
hist_V = nan(MaxIter+1,1);
Ne = size(ELEM,1);
Ndof = size(fem.f,1);

% Create output dir if needed
[folder,name] = fileparts(optoutput.storefileprefix);
if ~exist(folder, 'dir')
    mkdir(folder);
end

% Delete old output files
if (optoutput.storefinal || optoutput.storeiters)
    pattern = fullfile(folder, [name '*.mat']);
    files = dir(pattern);
    if ~isempty(files)
        delete(pattern);
    end
end

% Get interpolation functions 
[func_E,func_dEdx] = GetTOPinterpFunc(opt.matmodel);

% Volume constraint
opt.va = sum(fem.ve(opt.active));
opt.ve = fem.ve;
opt.v_func = @(x, opt)(sum((opt.H(opt.active,:) *x).*opt.ve(opt.active))/opt.va);

% Start iterations
fprintf('=== ITERATIONS BEGIN... ===\n')

for iter=1:MaxIter
    % Continuation update 
    [x,opt] = ShapeTOPcontinuation(x,fem,opt,iter);
    rho = opt.H * x;
    E = opt.Emin + func_E(rho,opt.current_penal)*(opt.E0-opt.Emin);

    % Classify current topology into phases
    opt = ClassifyShapeElems(x,opt);

    % Plot densities
    if optoutput.displayflag
       optoutput = ...
           PlotCurrentDensity(NODE, ELEM, rho, optoutput, iter);
       optoutput = ...
           PlotCurrentPhase(NODE, ELEM, x, opt, optoutput, iter);
    end

    % FEM Analysis
    ke = pagemtimes(reshape(E,[1 1 Ne]),fem.k0);
    fem.K = sparse(fem.indI(:),fem.indJ(:),ke(:));
    fem = FEMsolve(fem);

    % Objective function and sensitivity analysis
    dE_drho = func_dEdx(rho,opt.current_penal)*(opt.E0-opt.Emin);
    switch opt.objective
        case 'compliance'
            U = permute(fem.u(fem.edof),[2 1 3]);
            UT = fem.u(fem.edof);
            K = fem.k0;
            J0 = squeeze(pagemtimes(pagemtimes(UT,K),U));
            opt.J = sum(E.*J0);
            opt.dJ = dE_drho.*-J0;
        case 'displacement'
            Kf = fem.K(fem.freeDOF, fem.freeDOF);
            ut = fem.L' * (fem.L * fem.u);
            uf = ut(fem.freeDOF);
            lambda = zeros(Ndof, 1);
            lambda(fem.freeDOF) = 2 * (Kf \ uf);
            U = permute(fem.u(fem.edof),[2 1 3]);
            LT = lambda(fem.edof);
            K = fem.k0;
            J0 = squeeze(pagemtimes(pagemtimes(LT,K),U));
            opt.J = fem.u' * fem.L' * fem.L * fem.u;
            opt.dJ = dE_drho.*-J0;
        otherwise
           error('Unknown objective. Cannot solve.')
    end

    opt.v = opt.v_func(x, opt);
    opt.dv = fem.ve /opt.va;
    opt.dJ(~opt.active) = 0; 
    
    opt.dJ = opt.H' * opt.dJ;
    opt.dv = opt.H' * opt.dv;
  
    % Add history
    hist_J(iter) = opt.J;
    hist_V(iter) = opt.v;

    % Save iter data
    if optoutput.storeiters
        filename = sprintf('%s%03.0f.mat',optoutput.storefileprefix,iter-1);
        J = opt.J;
        V = opt.v;
        save(filename,'rho','x','J','V');
    end

    % Design variables update
    opt.xLB = max(opt.xmin,max(opt.flagPureSolid*(0.5+opt.current_mu),x-opt.current_move));
    opt.xUB = min(max(opt.flagPureVoid*(0.5-opt.current_mu),1-opt.flagPureVoid),x+opt.current_move);
    [x,opt] = TOPupdateDV(x,opt);
    
    fprintf(' It.:%03.0f Obj.:%g Vol.:%g\n',iter,opt.J,opt.v);
end
rho = opt.H * x;
opt = ClassifyShapeElems(x,opt);

% Plot densities
if optoutput.displayflag
    optoutput = ...
        PlotCurrentDensity(NODE, ELEM, rho, optoutput, MaxIter + 1);
    optoutput = ...
        PlotCurrentPhase(NODE, ELEM, x, opt, optoutput, MaxIter + 1);
end

% Save history data
if (optoutput.storefinal || optoutput.storeiters)
    J = NaN; 
    V = NaN; 
    hist_J(end) = NaN;
    hist_V(end) = NaN;
    filename = sprintf('%s%03.0f.mat',optoutput.storefileprefix,MaxIter);
    save(filename,'rho','x','J','V');
    filename = sprintf('%sMESH.mat',optoutput.storefileprefix);
    exportdata = optoutput.exportdata;
    exportdata.hist_J = hist_J;
    exportdata.hist_V = hist_V;
    save(filename,'-struct','exportdata');
end

end