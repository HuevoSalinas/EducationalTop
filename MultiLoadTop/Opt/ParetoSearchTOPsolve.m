function [MaxIter, hist_J_comb, hist_V_comb, w_comb, x_comb, rho, opt, fem] = ...
    ParetoSearchTOPsolve(NODE, ELEM, SUPP, LOAD, fem, opt, x, optoutput)

% Initial parameters
Ne = size(ELEM,1);              % Number of elements
Nl = size(fem.f, 2);            % Number of load cases
switch opt.objective            % Number of objectives
    case 'compliance'
        No = Nl;
    case 'displacement'
        No_list = zeros(Nl, 1); % Number of objectives for compliance objective
        for i = 1:Nl
            No_list(i) = size(fem.L{i},2);
        end
        No = sum(No_list);      % Number of objectives for displacement objectives
    otherwise
        error('Unknown objective. Cannot solve.')
end
MaxIter = opt.contin(end,1);
Nw = opt.Nw;
x0 = x;

% Create result variables
x_comb = nan(Nw, Ne );                % Final design variables for each weight combination 
J_comb = nan(Nw, No );                % Final objective for each weight combination 
w_comb = nan(Nw, No );                % Weight combinations 
V_comb = nan(Nw, 1 );                 % Final volume for each weight combination 
hist_J_comb = nan(Nw, MaxIter+1, No); % Objective history for each weight combination 
hist_V_comb = nan(Nw, MaxIter+1, 1 ); % Volume history for each weight combination 

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

opt.w = ones(No, 1)/No;
WupdateFunc = @ArbitraryWeightUpdate; % Function for weight updating
for e = 1:Nw
    x = x0;
    fprintf(' It.:%03.0f Weight\n',e);
   
    [MaxIter, hist_J, hist_V, ~, x, rho, opt, fem] = ...
            VariableWeightTOPsolve(NODE, ELEM, SUPP, LOAD, fem, opt, x, ...
                                  WupdateFunc, optoutput); % Variable weight optimization
    % Save weight data
    J_comb(e, :) = opt.MJ;
    w_comb(e, :) = opt.w';
    V_comb(e, :) = opt.v;
    x_comb(e, :) = x;
    hist_J_comb(e, :, :) = hist_J;
    hist_V_comb(e, :, :) = hist_V;

    % Update weight
    [opt, w] = DualWeightUpdate(opt); % Dual update for weight combinations
    opt.w = w;

end

% Plot densities
if optoutput.displayflag
    [optoutput] = ...
        PlotCurrentDensity(NODE, ELEM, rho, optoutput, 1);
     [optoutput] = ...
        PlotCurrentDensity(NODE, ELEM, rho, optoutput, MaxIter + 1);
    drawnow
end

% Save history data
if (optoutput.storefinal || optoutput.storeiters)
    filename = sprintf('%s.mat',optoutput.storefileprefix);
    save(filename,'x_comb','J_comb','V_comb','w_comb');
    filename = sprintf('%sMESH.mat',optoutput.storefileprefix);
    save(filename,'NODE','ELEM','SUPP','LOAD','hist_J_comb','hist_V_comb','w_comb');
end

end
