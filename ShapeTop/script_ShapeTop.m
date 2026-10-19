% SHAPE TOPOLOGY OPTIMIZATION
%
% This script solves a 2D shape optimization problem with Q4 elements.
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
problemID = 1;
elemType = 'Q4';
constLaw = 'PS';

[NODE, ELEM, SUPP, LOAD, DISP, x0, volfrac, pass_solid, pass_void, rmin] ...
    = LoadShapeOpt2DProblem(problemID);
% prob 1: MBB with circular holes
% prob 2: Beam height
% prob 3: Vierendeel bridge
% prob 4: L shaped beam
% prob 5: Plane seat bracket
% prob 6: Full MBB beam

%% == Parameters ==========================================================
% --- User defined FEM parameters -----------------------------------------
mat = struct('E',129,'nu',0.3,'t',1,'constLaw',constLaw);
fem = struct('intType','S','elemType',elemType,'solveMethod','direct');

% --- User defined optimization parameters --------------------------------
opt = struct('E0', 1, 'Emin', 1e-9, 'xmin', 1e-3, 'volfrac', volfrac,...
             'filterS', 1.5, 'eta', 0.5,...
             'matmodel','SIMP', 'objective', 'compliance',...
             'shapeElemNeigh','edge','shapeIsBorderTrans', true);
%             iter P   mu  move rmin
opt.contin = [  0  3.5  0.4 0.2  rmin;
              200  NaN  NaN NaN  NaN];
opt.pass_solid = pass_solid;
opt.pass_void = pass_void;
opt.active = (~opt.pass_solid & ~opt.pass_void);

% --- User defined output parameters --------------------------------------
exportdata = struct('NODE', NODE, 'ELEM', ELEM, 'SUPP', SUPP, 'LOAD', LOAD);
optoutput = struct('displayflag', true,'storeiters', true,'storefinal', true, ...
                   'storefileprefix',fullfile('output', 'shapeout'),...
                   'exportdata', exportdata);
dispFactor = 0.1;

%% == Domain Properties ===================================================
% --- Number of Nodes and Elements ----------------------------------------
Ne = size(ELEM,1);

% --- Plot Domain ---------------------------------------------------------
PlotDomain(NODE,ELEM,SUPP,LOAD);

% --- FEM data ------------------------------------------------------------
[mat, fem] = ParseFEMMAT(mat, fem, NODE, ELEM, SUPP, LOAD, DISP);
dNfunc = @ShapeDerivQ4;
switch fem.intType
    case 'F', intOrder = 2;
    case 'S', intOrder = [2 1];
    case 'R', intOrder = 1;
    otherwise, error('intType can only be F, S, or R.')
end
[fem.k0,fem.indI] = intK_continuum2(NODE,ELEM,mat,mat.constLaw,...
                                    dNfunc,intOrder);
fem.indJ = permute(fem.indI,[2 1 3]);
fem.edof = permute(GetElemDOF(ELEM,fem.DpN),[3 2 1]);

%% == Optimization ========================================================
% --- Initial guess -------------------------------------------------------
x = x0;
x(opt.pass_solid) = 1;
x(opt.pass_void) = opt.xmin;

% --- Connectivity --------------------------------------------------------
opt = GetElemConnectivity(ELEM,opt);

% --- Shape optimization --------------------------------------------------
[MaxIter, hist_J, hist_V, x, rho, opt, fem] = ...
    ShapeTOPsolve(NODE, ELEM, fem, opt, x, optoutput);

% --- Optimization results ------------------------------------------------
PlotResults(MaxIter, hist_J, hist_V)

%% == Initial Distribution Analysis =======================================
% --- FEM Analysis --------------------------------------------------------
func_E = GetTOPinterpFunc(opt.matmodel);
rho0 = opt.H * x0;
E = opt.Emin + func_E(rho0,opt.current_penal)*(opt.E0-opt.Emin);
ke = pagemtimes(reshape(E,[1 1 Ne]),fem.k0);
fem.K = sparse(fem.indI(:),fem.indJ(:),ke(:));
fem = FEMsolve(fem);
D0 = fem.D;

%% == Final Distribution Analysis =========================================
% --- FEM Analysis --------------------------------------------------------
func_E = GetTOPinterpFunc(opt.matmodel);
E = opt.Emin + func_E(rho,opt.current_penal)*(opt.E0-opt.Emin);
ke = pagemtimes(reshape(E,[1 1 Ne]),fem.k0);
fem.K = sparse(fem.indI(:),fem.indJ(:),ke(:));
fem = FEMsolve(fem);
D = fem.D;

%% == Results =============================================================
% --- Plot Displacements -------------------------------------------------- 
PlotDisplacementsDiffTOP(NODE, ELEM, D0, rho0, D, rho, dispFactor);

%% === Remove paths =======================================================
RemoveLocalPaths() % Remove subfolders
rmpath('./Path');