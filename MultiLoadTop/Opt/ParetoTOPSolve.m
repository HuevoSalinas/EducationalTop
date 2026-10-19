function [MaxIter, hist_J_comb, hist_V_comb, w_comb, x_comb, rho, opt, fem] = ...
    ParetoTOPSolve(NODE, ELEM, SUPP, LOAD, fem, opt, x, optoutput)

% Initial parameters
Ne = size(ELEM,1);         % Number of elements
Nl = size(fem.f, 2);       % Number of load cases
switch opt.objective       % Number of objectives
    case 'compliance'
        No = Nl;           % Number of objectives for compliance objective
    case 'displacement'
        No_list = zeros(Nl, 1);
        for i = 1:Nl
            No_list(i) = size(fem.L{i},2);
        end
        No = sum(No_list); % Number of objectives for displacement objectives
    otherwise
        error('Unknown objective. Cannot solve.')
end
MaxIter = opt.contin(end,1);
dw = opt.movew;
x0 = x;

% Get weight combinations
w = 0:dw:1;
n = No;
nw = size(w,2);
w_comb = nan( nw^(n), n);
tol = 1e-5 * dw;
for i = 1:n
    for j = 1:n-1
        if j == 1
            if j < i
                w_column = repmat(w, 1, nw);
            else
                w_column = reshape(repmat(w, nw, 1), 1, nw^2);
            end
        else
            if j < i
                w_column = repmat(w_column,1, nw);
            else
                w_column = reshape(repmat(w_column,nw, 1), 1, nw^(j + 1));
            end
        end
    end
    w_comb(:, i) = w_column;
end
w_comb = w_comb(abs(sum(w_comb,2) - 1)<tol, :);    % Weight combinations with total weight equal to 1

% Create result variables
x_comb = nan(size(w_comb, 1), Ne );                % Final design variables for each weight combination 
J_comb = nan(size(w_comb, 1), No );                % Final objective for each weight combination 
V_comb = nan(size(w_comb, 1), 1 );                 % Final volume for each weight combination 
hist_J_comb = nan(size(w_comb, 1), MaxIter+1, No); % Objective history for each weight combination 
hist_V_comb = nan(size(w_comb, 1), MaxIter+1, 1 ); % Volume history for each weight combination 

% Create output dir if needed
folder = extractBefore(optoutput.storefileprefix,strfind(optoutput.storefileprefix, filesep));
if ~exist(folder, 'dir')
    mkdir(folder);
end

% Delete old output files
if (optoutput.storefinal || optoutput.storeiters) % Delete previous data
    pattern = fullfile(optoutput.storefileprefix, '*.mat');
    files = dir(pattern);
    if ~isempty(files)
        delete(pattern);
    end
end

WupdateFunc = @ArbitraryWeightUpdate; % Function for weight updating
for e = 1:size(w_comb, 1)
    opt.w = w_comb(e, :)';
    x = x0;
    fprintf(' It.:%03.0f Weight\n',e);
   
    [MaxIter, hist_J, hist_V, ~, x, ~, opt, fem] = ...
            VariableWeightTOPSolve(NODE, ELEM, SUPP, LOAD, fem, opt, x, ...
                                  WupdateFunc, optoutput); % Variable weight optimization
    % Save weight data
    J_comb(e, :) = opt.MJ;
    V_comb(e, :) = opt.v;
    x_comb(e, :) = x;
    hist_J_comb(e, :, :) = hist_J;
    hist_V_comb(e, :, :) = hist_V;

end
% Update for best solution
[~, ind] = min(max(hist_J_comb(:,MaxIter,:), [], 3));
x = x_comb(ind, :)';
rho = opt.H * x;

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
    J = NaN; 
    filename = sprintf('%s%03.0f.mat',optoutput.storefileprefix,MaxIter);
    save(filename,'rho','J');
    filename = sprintf('%sMESH.mat',optoutput.storefileprefix);
    save(filename,'NODE','ELEM','SUPP','LOAD','hist_J_comb','hist_V_comb','w_comb');
end

return 
