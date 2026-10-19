% MULTI-LOAD TOPOLOGY OPTIMIZATION 
%
% This script solves a multi-load topology optimization problem.
% The objective is to minimize a combination of objectives for multiple load cases.
%
% Main steps:
% 1) Define problem
% 2) Assemble FEM matrices
% 3) Run optimization loop
% 4) Post-process results

clc, clear, close all
%% === Add paths ==========================================================
addpath('./Path');
AddLocalPaths() % Load subfolders

%% == Examples ============================================================
problemID = 1;   % Problem identification number
elemType = 'Q4'; % Finite element type
constLaw = 'PS'; % Constitutive law for the element

[NODE, ...       % Nodes coodinates [Nn x 2] 
 ELEM, ...       % Elements connectivity [Ne x 4]
 SUPP, ...       % Supported nodes and imposed displacements [node, dx, dy]
 LOAD, ...       % Loaded nodes [node, fx, fy]
 DISP, ...       % Control points [id, node, dx, dy]
 pass_solid, ... % Passive solid elements
 pass_void,...   % Passive void elements
 rmin] ...       % Filter radius 
    = LoadMultiLoadOpt2DPlaneProblem(problemID); % Load problem 
% prob 1: Beam
% prob 2: Bridge
% prob 3: Building 

%% == Parameters ==========================================================
% --- User defined FEM parameters -----------------------------------------
mat = struct('E',129, ...             % Elasticity modulus of solid material  
             'nu',0.3, ...            % Poisson ratio of solid material 
             't', 1, ...              % Element thickness 
             'constLaw',constLaw);    % Element constitutive law
fem = struct('intType','S',...        % Element integration type
             'elemType',elemType,...  % Finite element type
             'solveMethod','direct'); % Solve method for Ku=f

% --- User defined optimization parameters ----------------------------------
opt = struct(...
 'E0',1, ...                            % Upper bound for elasticity modulus
 'Emin',1e-9,...                        % Lower bound for elasticity modulus
 'volfrac',0.3, ...                     % Maximum volume fraction 
 'xmin',1e-3, ...                       % Lower bound for design variables
 'filterS',1.5,...                      % Filter parameter
 'eta',0.3, ...                         % Damping parameter for oc update method
 'matmodel','SIMP', ...                 % Material interpolation model
 'objective','compliance', ...          % Objective function 
 'combinationMethod','bounded-dual',... % Objectives combination method 
 'movew', 0.2,...                       % Weight move (for bounded-dual & pareto methods)
 'Nw',5, ...                            % Number of weights (for pareto-search method)
 'norm',3, ...                          % p or eta (for p-norm & ks methods)
 'w',[0.4 0.1 0.5]');                       % Initial weights (for arbitrary method)

% Objectives combination method options:
% pareto
% pareto-search
% arbitrary
% p-norm
% ks
% bounded-dual
% bound-constrained

%             iter P  move rmin
opt.contin = [  0  3.  0.25  rmin;
              100 NaN NaN  NaN];
opt.pass_solid = pass_solid;
opt.pass_void = pass_void;
opt.active = (~opt.pass_solid & ~opt.pass_void);

% --- User defined output parameters --------------------------------------
optoutput = struct('displayflag', true,'storeiters', true,'storefinal', true, ...
                   'storefileprefix',fullfile('output', 'topout'));
dispFactor = 1;

%% == Domain Properties ===================================================
% --- Number of Nodes and Elements ----------------------------------------
Nn = size(NODE,1); Ne = size(ELEM,1);

% --- FEM data ------------------------------------------------------------
[mat, fem] = ParseMultiLoadFEMMAT(mat, fem, NODE, ELEM, SUPP, LOAD, DISP); % Inputs parsing
dNfunc = @ShapeDerivQ4; % Derivative of the corresponding shape function 
switch fem.intType
    case 'F', intOrder = 2; % 2x2 rule
    case 'S', intOrder = [2 1]; % 2x2 rule for volumetric & 1x1 for shear
    case 'R', intOrder = 1; % 1x1 rule
    otherwise, error('intType can only be F, S, or R.')
end
[fem.k0,fem.indI] = intK_continuum2(NODE,ELEM,mat,mat.constLaw,...
                                    dNfunc,intOrder); % Stiffness matrix integration
fem.indJ = permute(fem.indI,[2 1 3]);
fem.edof = permute(GetElemDOF(ELEM,fem.DpN),[3 2 1]); % Element degrees of freedom

% --- Plot Domain ---------------------------------------------------------
PlotMultiLoadDomain(NODE,ELEM,SUPP,LOAD); % Domain plotting 

% --- Initial guess -------------------------------------------------------
x = repmat(opt.volfrac,Ne,1); % Initial uniform distribution 
x(opt.pass_solid) = 1;        % Desity 1 for passive solid elements
x(opt.pass_void) = opt.xmin;  % Density xmin for passive void elements

%% == Optimization ========================================================
% --- Topology optimization -----------------------------------------------
[MaxIter, hist_J, hist_V, hist_w, x, rho, opt, fem] = ...
    MultiLoadTOPSolve(NODE, ELEM, SUPP, LOAD, fem, opt, x, optoutput);

% --- Optimization results ------------------------------------------------
PlotMultiLoadResults(opt.combinationMethod, MaxIter, hist_J, hist_V, hist_w)

%% === Remove paths =======================================================
RemoveLocalPaths() % Remove subfolders
rmpath('./Path');