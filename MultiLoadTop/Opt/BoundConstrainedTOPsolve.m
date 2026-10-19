function [MaxIter, hist_J, hist_V, x, rho, opt, fem] = ...
    BoundConstrainedTOPsolve(NODE, ELEM, SUPP, LOAD, fem, opt, x, optoutput)

% Initial parameters
Ne = size(ELEM,1);           % Number of elements
Nl = size(fem.f, 2);         % Number of load cases
Na = sum(opt.active);        % Number of active elements
Ndof = size(fem.f,1);        % Number of degrees of freedom 
MaxIter = opt.contin(end,1); % Number of maximum iterations
switch opt.objective
    case 'compliance'
        No = Nl;             % Number of objectives for compliance objective
    case 'displacement'
        No_list = zeros(Nl, 1);
        for i = 1:Nl
            No_list(i) = size(fem.L{i},2);
        end
        No = sum(No_list);   % Number of objectives for displacement objectives
    otherwise
        error('Unknown objective. Cannot solve.')
end
hist_J = nan(MaxIter+1, No);        % Objective history
hist_V = nan(MaxIter+1, 1);         % Volume history 
hist_x = nan(MaxIter+1, size(x,1)); % Design variables history 
hist_gamma = nan(MaxIter+1, 1 );    % Design variables history 
low = [];
upp = [];

% Get interpolation functions 
[func_E,func_dEdx] = GetTOPinterpFunc(opt.matmodel);

% Create output dir if needed
[folder,name] = fileparts(optoutput.storefileprefix);
if ~exist(folder, 'dir')
    mkdir(folder);
end

% Delete old output files
if (optoutput.storefinal || optoutput.storeiters) % Delete previous data
    pattern = fullfile(folder, [name '*.mat']);
    files = dir(pattern);
    if ~isempty(files)
        delete(pattern);
    end
end

opt.va = sum(fem.ve(opt.active)); % Volume of active elements
opt.ve = fem.ve;                  % Element volume
opt.v_func = @(x, opt)(sum((opt.H(opt.active,:) *x).*opt.ve(opt.active))/opt.va); % Volume function

for iter=1:MaxIter
    % Continuation update (update filter matrix if needed)
    opt = TOPcontinuation(fem,opt,iter);
    rho = opt.H * x; % Filtered density variables
    E = opt.Emin + func_E(rho,opt.current_penal)*(opt.E0-opt.Emin); % Penalized modulus

    % Plot densities
    if optoutput.displayflag
       [optoutput] = ...
           PlotCurrentDensity(NODE, ELEM, rho, optoutput, iter);
        drawnow
    end
    
    % FEM Analysis
    ke = pagemtimes(reshape(E,[1 1 Ne]),fem.k0);   % Penalized stiffness matrices
    fem.K = sparse(fem.indI(:),fem.indJ(:),ke(:)); % Assembled stiffness matrix 
    fem = FEMsolve(fem);           % Displacement computation

    % Objective function and sensitivity analysis
    dE_drho = func_dEdx(rho,opt.current_penal)*(opt.E0-opt.Emin); % Elasticity modulus sensitivity
    switch opt.objective
        case 'compliance'
            J = zeros(1, No);
            dJ = zeros(Ne, No);
            K = fem.k0;
            for j = 1:No
                uj = fem.u(:, j);
                U = permute(uj(fem.edof),[2 1 3]);
                UT = uj(fem.edof);
                J0 = squeeze(pagemtimes(pagemtimes(UT,K),U));
                J(:, j) =  sum(E.*J0);
                dJ(:, j) = dE_drho.*-J0;
            end
            opt.MJ = J;
            opt.MdJ = dJ;
        case 'displacement'
            J = zeros(1, No);
            dJ = zeros(Ne, No);
            Kf = fem.K(fem.freeDOF, fem.freeDOF);
            [L,~,s] = chol(Kf,'lower','vector');
            K = fem.k0;
            for j = 1:Nl
                uj = fem.u(:, j);
                U = permute(uj(fem.edof),[2 1 3]);
                for k = 1: No_list(j)
                    ut = fem.L{j}{k}' * (fem.L{j}{k} * uj);
                    uf = ut(fem.freeDOF);
                    lambda = zeros(Ndof, 1);
                    lambda(fem.freeDOF(s), :) = 2 * (L'\(L\uf(s,:)));  % Adjoint method
                    LT = lambda(fem.edof);
                    J0 = squeeze(pagemtimes(pagemtimes(LT,K),U));
                    ind = k + sum(No_list(1:j-1));
                    J(:, ind) = uj' * fem.L{j}{k}' * fem.L{j}{k} * uj;
                    dJ(:, ind) =  dE_drho.*-J0;
                end
            end
            opt.MJ = J;
            opt.MdJ = dJ;
        otherwise
           error('Unknown objective. Cannot solve.')
    end
    
    opt.v = opt.v_func(x, opt);       % Current volume 
    opt.dv = fem.ve /opt.va;          % Volume sensitivity
    opt.MdJ(~opt.active, :) = 0;      % Null sensitivity for passive elements

    opt.MdJ = opt.H' * opt.MdJ;       % Objectives sensitivities considering filter
    opt.dv = opt.H' * opt.dv;         % Volume sensitivities considering filter

    % Get gamma
    if iter == 1
        gamma = max(opt.MJ);          % Initial gamma value 
        gamma_ub = max(opt.MJ);       % gamma upper bound
    end

    % Add history
    hist_J(iter, :) = opt.MJ;
    hist_V(iter, :) = opt.v;
    hist_gamma(iter, :) = gamma;
    hist_x(iter, :) = x;

    % Variables update
    x_lb = max(zeros(Na, 1), x(opt.active,:) - opt.current_move); % Design variables lower bound
    x_ub = min(ones(Na, 1), x(opt.active,:) + opt.current_move);  % Design variables upper bound
    
    rparam = 10;
    f = gamma/gamma_ub*rparam;               % Normalized objective
    df =  [zeros(Na, 1); 1]/gamma_ub*rparam; % Normalized objective sensitivity
    df2 = zeros(Na + 1, 1);                  % Normalized objective second derivative   
    g = [(opt.MJ' - gamma)/gamma_ub;
        opt.v_func(x, opt) - opt.volfrac];   % Normalized constraints
    dg = [[opt.MdJ(opt.active, :); -ones(1,No)]'/gamma_ub; 
          [opt.dv(opt.active, :); 0]'];      % Normalized constraints sensitivity
    dg2 = zeros(No + 1, Na + 1);             % Normalized constraints second derivative   
    a0 = 1;                    % mma param a0
    a = zeros(No + 1, 1);      % mma param a
    c = 1e3 * ones(No + 1, 1); % mma param c
    d = zeros(No + 1, 1);      % mma param d
   
    if iter == 1
        xold1 = x(opt.active);
        xold2 = x(opt.active);
        gammaold1 = gamma;
        gammaold2 = gamma;
    elseif iter == 2
        xold2 = x(opt.active);
        xold1 = hist_x(iter - 1, opt.active)';
        gammaold2 = gamma;
        gammaold1 = hist_gamma(iter - 1, :);
    else
        xold2 = hist_x(iter - 2, opt.active)';
        xold1 = hist_x(iter - 1, opt.active)';
        gammaold2 = hist_gamma(iter - 2, :);
        gammaold1 = hist_gamma(iter - 1, :);
    end
    warning off all
    [var,~ ,~, ~, ~, ~, ~, ~, ~,low,upp] = ...
        mmasub(No + 1, Na + 1, iter, [x(opt.active); gamma], [x_lb; 0 ], [x_ub; gamma_ub], ...
               [xold1 ; gammaold1], [xold2 ; gammaold2], ...
        f, df, df2, g, dg, dg2, low, upp, a0,a,c,d); % Variables updating
    warning on all
    x = zeros(size(x));
    x(opt.active) = var(1:Na, :);
    x(opt.pass_solid) = 1;
    x(opt.pass_void) = opt.xmin;
    gamma = var(Na + 1, :);

    % Save iter data
    if optoutput.storeiters
        filename = sprintf('%s%03.0f.mat',optoutput.storefileprefix,iter-1);
        J = opt.MJ;
        save(filename,'rho','J');
    end
    fprintf(' It.:%03.0f Obj.:%g Vol.:%g\n',iter,max(opt.MJ),opt.v);
end
rho = opt.H * x;

% Plot densities
if optoutput.displayflag
    [optoutput] = ...
        PlotCurrentDensity(NODE, ELEM, rho, optoutput, MaxIter + 1);
    drawnow
end

% Save history data
if (optoutput.storefinal || optoutput.storeiters)
    J = NaN;
    hist_J(end) = NaN;
    hist_V(end) = NaN;
    filename = sprintf('%s%03.0f.mat',optoutput.storefileprefix,MaxIter);
    save(filename,'rho','J');
    filename = sprintf('%sMESH.mat',optoutput.storefileprefix);
    save(filename,'NODE','ELEM','SUPP','LOAD','hist_J','hist_V');
end

end
