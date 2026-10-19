function [MaxIter, hist_J, hist_V, hist_w, x, rho, opt, fem] = ...
    VariableWeightTOPsolve(NODE, ELEM, SUPP, LOAD, fem, opt, x, WupdateFunc, optoutput)

% Initial parameters
Ne = size(ELEM,1);           % Number of elements
Nl = size(fem.f, 2);         % Number of load cases 
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
     
hist_J = nan(MaxIter+1, No);  % Objective history
hist_V = nan(MaxIter+1, 1);   % Volume history 
hist_w = nan(MaxIter+1, No);  % Weights history

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

% START ITERATION
fprintf('=== ITERATIONS BEGIN... ===\n')
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
       [optoutput] = PlotCurrentDensity(NODE, ELEM, rho, optoutput, iter);
        drawnow
    end

    % FEM Analysis
    ke = pagemtimes(reshape(E,[1 1 Ne]),fem.k0);   % Penalized stiffness matrices
    fem.K = sparse(fem.indI(:),fem.indJ(:),ke(:)); % Assembled stiffness matrix 
    fem = FEMsolve(fem);           % Displacement computation

    % Objective function 
    switch opt.objective
        case 'compliance'
            J = zeros(1, No);
            for j = 1:No
                uj = fem.u(:, j);
                J(:, j) =  uj' * fem.K * uj;
            end
            opt.MJ = J; 
        case 'displacement'
            J = zeros(1, No);
            for j = 1:Nl
                uj = fem.u(:, j);
                for k = 1: No_list(j)
                    ind = k + sum(No_list(1:j-1));
                    J(:,ind) = uj' * fem.L{j}{k}' * fem.L{j}{k} * uj;
                end
            end
            opt.MJ = J; 
        otherwise
           error('Unknown objective. Cannot solve.')
    end

    % Weights update
    [opt, w] = WupdateFunc(opt);

    % Objective function sensitivity
    dE_drho = func_dEdx(rho,opt.current_penal)*(opt.E0-opt.Emin); % Elasticity modulus sensitivity
    switch opt.objective
        case 'compliance'
            dJ = zeros(Ne, No);
            K = fem.k0;
            for j = 1:No
                uj = fem.u(:, j);
                U = permute(uj(fem.edof),[2 1 3]);
                UT = uj(fem.edof);
                J0 = squeeze(pagemtimes(pagemtimes(UT,K),U));
                dJ(:, j) =  dE_drho.*-J0;
            end
            opt.dJ = dJ * w;
        case 'displacement'
            opt.dJ = zeros(Ne, 1);
            Kf = fem.K(fem.freeDOF, fem.freeDOF);
            [L,~,s] = chol(Kf,'lower','vector');
            K = fem.k0;
            for j = 1:Nl
                uj = fem.u(:, j);
                ut = zeros(Ndof, 1);
                for k = 1: No_list(j)
                    ind = k + sum(No_list(1:j-1));
                    ut = ut + w(ind) * fem.L{j}{k}' * (fem.L{j}{k} * uj);
                end
                uf = ut(fem.freeDOF);
                lambda = zeros(Ndof, 1);
                lambda(fem.freeDOF(s), :) = 2 * (L'\(L\uf(s,:))); % Adjoint method
                U = permute(uj(fem.edof),[2 1 3]);
                LT = lambda(fem.edof);
                J0 = squeeze(pagemtimes(pagemtimes(LT,K),U));
                opt.dJ = opt.dJ - dE_drho .* J0;
            end
        otherwise
           error('Unknown objective. Cannot solve.')
    end
     
    opt.v = opt.v_func(x, opt);       % Current volume 
    opt.dv = fem.ve /opt.va;          % Volume sensitivity
    opt.dJ(~opt.active) = 0;      % Null sensitivity for passive elements

    opt.dJ = opt.H' * opt.dJ;       % Objectives sensitivities considering filter
    opt.dv = opt.H' * opt.dv;         % Volume sensitivities considering filter


    % Design variables update
    opt.xLB = max(opt.xmin,x-opt.current_move); % Design variables lower bound
    opt.xUB = min(1,x+opt.current_move);        % Design variables upper bound
    [x,opt] = TOPupdateDV(x,opt);               % Design update

     % Add history
    hist_J(iter, :) = opt.MJ;
    hist_V(iter, :) = opt.v;
    hist_w(iter, :) = w;

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
    [optoutput] = PlotCurrentDensity(NODE, ELEM, rho, optoutput, MaxIter + 1);
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
    save(filename,'NODE','ELEM','SUPP','LOAD','hist_J','hist_V','hist_w');
end

end
