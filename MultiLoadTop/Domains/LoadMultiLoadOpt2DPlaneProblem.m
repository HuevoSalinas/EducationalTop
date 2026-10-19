function [NODE, ELEM, SUPP, LOAD, DISP, pass_solid, pass_void, rmin] ...
    = LoadMultiLoadOpt2DPlaneProblem(problemID)

% Load problem
switch problemID
    case 1
        [NODE, ELEM, SUPP, LOAD, DISP, pass_solid, pass_void, rmin] = ...
            BeamProblem();
    case 2
        [NODE, ELEM, SUPP, LOAD, DISP, pass_solid, pass_void, rmin] = ...
            BridgeProblem();
    case 3
        [NODE, ELEM, SUPP, LOAD, DISP, pass_solid, pass_void, rmin] = ...
            BuildingProblem();

    otherwise
        error('Unknown problem. Cannot load.')
end

end

function [NODE, ELEM, SUPP, LOAD, DISP, pass_solid, pass_void, rmin] ...
    = BeamProblem() 

nx = 200;
ny = 80;
L = 10;
H = 4;
NODE = zeros((nx +1)*(ny + 1), 2);
ELEM = zeros(nx * ny, 4);

%% Nodes
for i = 1:nx+1
    for j = 1:ny+1
        NODE(i + (j - 1) *(nx + 1), 1) = (i - 1) * L / nx;
        NODE(i + (j - 1) *(nx + 1), 2) = (j - 1) * H / ny;
    end
end

%% Elements
for i = 1:nx
    for j = 1:ny
        ELEM(i + (j - 1) *(nx), 1) = i + (j - 1) * (nx + 1);
        ELEM(i + (j - 1) *(nx), 2) = (i + 1)+ (j - 1) * (nx + 1);
        ELEM(i + (j - 1) *(nx), 3) = (i + 1)+ (j) * (nx + 1);
        ELEM(i + (j - 1) *(nx), 4) = (i)+ (j) * (nx + 1);
    end
end

%% Supports
tol = 1e-3;
ind = find(abs(NODE(:,1) - min(NODE(:,1))) < tol);
SUPP = [ind repmat([0 0],length(ind),1)];

%% Loads
LOAD1 = [(nx + 1)*(ny) + nx/2 + 1, 0, -2];
LOAD2 = [nx + 1, 0, -1];
LOAD3 = [(nx + 1)*(ny + 1), 0, -1];
LOAD = cell(1, 3);
LOAD{1} = LOAD1;
LOAD{2} = LOAD2;
LOAD{3} = LOAD3;

%% Displacement
DISP1 = [1, (nx + 1)*(ny) + nx/2 + 1, 0, 1];
DISP2 = [1, nx + 1, 0, 1];
DISP3 = [1, (nx + 1)*(ny + 1), 0, 1];
DISP = cell(1, 3);
DISP{1} = {DISP1};
DISP{2} = {DISP2};
DISP{3} = {DISP3};

%% Passive elements 
pass_solid = false(nx*ny, 1); 
pass_void = false(nx*ny, 1); 

%% Filter radius
rmin = 0.15;

end

function [NODE, ELEM, SUPP, LOAD, DISP, pass_solid, pass_void, rmin] = ...
    BridgeProblem() 

nx = 320;
ny = 80;
L = 40;
Lt = 4;
H = 8;
h = 0.1;
NODE = zeros((nx +1)*(ny + 1), 2);
ELEM = zeros(nx * ny, 4);

%% Nodes
for i = 1:nx+1
    for j = 1:ny+1
        NODE(i + (j - 1) *(nx + 1), 1) = (i - 1) * L / nx;
        NODE(i + (j - 1) *(nx + 1), 2) = (j - 1) * H / ny;
    end
end

%% Elements
for i = 1:nx
    for j = 1:ny
        ELEM(i + (j - 1) *(nx), 1) = i + (j - 1) * (nx + 1);
        ELEM(i + (j - 1) *(nx), 2) = (i + 1)+ (j - 1) * (nx + 1);
        ELEM(i + (j - 1) *(nx), 3) = (i + 1)+ (j) * (nx + 1);
        ELEM(i + (j - 1) *(nx), 4) = (i)+ (j) * (nx + 1);
    end
end

%% Supports
tol = 1e-3;
ind_1 = find(abs(NODE(:,2) - min(NODE(:,2)))< tol & abs(NODE(:,1) - L/8)< tol);
ind_2 = find(abs(NODE(:,2) - min(NODE(:,2)))< tol & abs(NODE(:,1) - 7*L/8)< tol);
ind_3 = find(abs(NODE(:,2) - H/2)< tol & abs(NODE(:,1) - min(NODE(:,1)))< tol);
ind_4 = find(abs(NODE(:,2) - H/2)< tol & abs(NODE(:,1) - max(NODE(:,1)))< tol);
SUPP = [ ind_1 0 0;
         ind_2 0 0;
         ind_3 0 0;
         ind_4 0 0]; 

%% Loads
n = 8;
LOAD = cell(1, 4*n);
for i = 1:n
    ind_1 = find(abs(NODE(:,2) - H/2) < tol & abs(NODE(:,1) - (i - 1)*L/2/n - L/2) < tol);
    ind_2 = find(abs(NODE(:,2) - H/2) < tol & abs(NODE(:,1) - (i - 1)*L/2/n - 3*Lt/4 - L/2) < tol);
    ind_3 = find(abs(NODE(:,2) - H/2) < tol & abs(NODE(:,1) - (i - 1)*L/2/n - Lt - L/2) < tol);
    LOAD_i = [ind_1 repmat([0 -4],length(ind_1),1);
              ind_2 repmat([0 -4],length(ind_2),1);
              ind_3 repmat([0 -1],length(ind_3),1)];
    LOAD{4*i} = LOAD_i;
    ind_1 = find(abs(NODE(:,2) - H/2) < tol & abs(NODE(:,1) - (i - 1)*L/2/n - L/2) < tol);
    ind_2 = find(abs(NODE(:,2) - H/2) < tol & abs(NODE(:,1) - (i - 1)*L/2/n + 3*Lt/4 - L/2) < tol);
    ind_3 = find(abs(NODE(:,2) - H/2) < tol & abs(NODE(:,1) - (i - 1)*L/2/n + Lt - L/2) < tol);
    LOAD_i = [ind_1 repmat([0 -4],length(ind_1),1);
              ind_2 repmat([0 -4],length(ind_2),1);
              ind_3 repmat([0 -1],length(ind_3),1)];
    LOAD{4*i - 1} = LOAD_i;
    ind_1 = find(abs(NODE(:,2) - H/2) < tol & abs(NODE(:,1) - (L - (i - 1)*L/2/n)+ L/2) < tol);
    ind_2 = find(abs(NODE(:,2) - H/2) < tol & abs(NODE(:,1) - (L - (i - 1)*L/2/n - 3*Lt/4)+ L/2) < tol);
    ind_3 = find(abs(NODE(:,2) - H/2) < tol & abs(NODE(:,1) - (L - (i - 1)*L/2/n - Lt)+ L/2) < tol);
    LOAD_i = [ind_1 repmat([0 -4],length(ind_1),1);
              ind_2 repmat([0 -4],length(ind_2),1);
              ind_3 repmat([0 -1],length(ind_3),1)];
    LOAD{4*i - 2} = LOAD_i;
    ind_1 = find(abs(NODE(:,2) - H/2) < tol & abs(NODE(:,1) - (L - (i - 1)*L/2/n)+ L/2) < tol);
    ind_2 = find(abs(NODE(:,2) - H/2) < tol & abs(NODE(:,1) - (L - (i - 1)*L/2/n + 3*Lt/4)+ L/2) < tol);
    ind_3 = find(abs(NODE(:,2) - H/2) < tol & abs(NODE(:,1) - (L - (i - 1)*L/2/n + Lt)+ L/2) < tol);
    LOAD_i = [ind_1 repmat([0 -4],length(ind_1),1);
              ind_2 repmat([0 -4],length(ind_2),1);
              ind_3 repmat([0 -1],length(ind_3),1)];
    LOAD{4*i - 3} = LOAD_i;
end

%% Displacement 
DISP = cell(1, 4*n);
ind = find(abs(NODE(:,2) - H/2)< tol);
for i = 1:n
    DISP{4*i} = {[(1:length(ind))' ind repmat([0 1],length(ind),1)]};
    DISP{4*i-1} = {[ (1:length(ind))' ind repmat([0 1],length(ind),1)]};
    DISP{4*i-2} = {[ (1:length(ind))' ind repmat([0 1],length(ind),1)]};
    DISP{4*i-3} = {[ (1:length(ind))' ind repmat([0 1],length(ind),1)]};
end

%% Passive elements 
Nfunc = @ShapeFuncQ4;
dNfunc = @ShapeDerivQ4;
intOrder = 1;
CM = GetElementCM2(NODE,ELEM,[],Nfunc,dNfunc,intOrder);
pass_solid = CM(:,2) < H/2 + h &  CM(:,2) > H/2 - h ; 
pass_void = false(nx*ny, 1); 

%% Filter radius
rmin = 0.25;

end

function [NODE, ELEM, SUPP, LOAD, DISP, pass_solid, pass_void, rmin] = ...
    BuildingProblem() 

nx = 80;
ny = 200;
L = 8;
H = 20;
d = 0.2;
NODE = zeros((nx +1)*(ny + 1), 2);
ELEM = zeros(nx * ny, 4);

%% Nodes
for i = 1:nx+1
    for j = 1:ny+1
        NODE(i + (j - 1) *(nx + 1), 1) = (i - 1) * L / nx;
        NODE(i + (j - 1) *(nx + 1), 2) = (j - 1) * H / ny;
    end
end

%% Elements
for i = 1:nx
    for j = 1:ny
        ELEM(i + (j - 1) *(nx), 1) = i + (j - 1) * (nx + 1);
        ELEM(i + (j - 1) *(nx), 2) = (i + 1)+ (j - 1) * (nx + 1);
        ELEM(i + (j - 1) *(nx), 3) = (i + 1)+ (j) * (nx + 1);
        ELEM(i + (j - 1) *(nx), 4) = (i)+ (j) * (nx + 1);
    end
end

%% Supports
tol = 1e-3;
ind_supp = find(abs(NODE(:,2) - min(NODE(:,2)))< tol);
SUPP = [ ind_supp    zeros(size(ind_supp,1),1)   zeros(size(ind_supp,1),1)]; 

%% Loads
LOAD = cell(1, 2);
ind_p = find((abs(NODE(:, 2) - 1*H/10) < tol & abs(NODE(:, 1)) < tol)| ...
           (abs(NODE(:, 2) - 2*H/10) < tol & abs(NODE(:, 1)) < tol)| ...
           (abs(NODE(:, 2) - 3*H/10) < tol & abs(NODE(:, 1)) < tol)| ...
           (abs(NODE(:, 2) - 4*H/10) < tol & abs(NODE(:, 1)) < tol)| ...
           (abs(NODE(:, 2) - 5*H/10) < tol & abs(NODE(:, 1)) < tol)| ...
           (abs(NODE(:, 2) - 6*H/10) < tol & abs(NODE(:, 1)) < tol)| ...
           (abs(NODE(:, 2) - 7*H/10) < tol & abs(NODE(:, 1)) < tol)| ...
           (abs(NODE(:, 2) - 8*H/10) < tol & abs(NODE(:, 1)) < tol)| ...
           (abs(NODE(:, 2) - 9*H/10) < tol & abs(NODE(:, 1)) < tol)| ...
           (abs(NODE(:, 2) - 10*H/10) < tol & abs(NODE(:, 1)) < tol));
ind_n = find((abs(NODE(:, 2) - 1*H/10) < tol & abs(NODE(:, 1)- L) < tol)| ...
           (abs(NODE(:, 2) - 2*H/10) < tol & abs(NODE(:, 1)- L) < tol)| ...
           (abs(NODE(:, 2) - 3*H/10) < tol & abs(NODE(:, 1)- L) < tol)| ...
           (abs(NODE(:, 2) - 4*H/10) < tol & abs(NODE(:, 1)- L) < tol)| ...
           (abs(NODE(:, 2) - 5*H/10) < tol & abs(NODE(:, 1)- L) < tol)| ...
           (abs(NODE(:, 2) - 6*H/10) < tol & abs(NODE(:, 1)- L) < tol)| ...
           (abs(NODE(:, 2) - 7*H/10) < tol & abs(NODE(:, 1)- L) < tol)| ...
           (abs(NODE(:, 2) - 8*H/10) < tol & abs(NODE(:, 1)- L) < tol)| ...
           (abs(NODE(:, 2) - 9*H/10) < tol & abs(NODE(:, 1)- L) < tol)| ...
           (abs(NODE(:, 2) - 10*H/10) < tol & abs(NODE(:, 1)- L) < tol));

LOAD{1} = [ind_p repmat([1/10 0],length(ind_p),1)];
LOAD{2} = [ind_n repmat([-1/10 0],length(ind_n),1)];

%% Displacement
ind = find((abs(NODE(:, 2) - 1*H/10) < tol & abs(NODE(:, 1)- L/2) < tol)| ...
           (abs(NODE(:, 2) - 2*H/10) < tol & abs(NODE(:, 1)- L/2) < tol)| ...
           (abs(NODE(:, 2) - 3*H/10) < tol & abs(NODE(:, 1)- L/2) < tol)| ...
           (abs(NODE(:, 2) - 4*H/10) < tol & abs(NODE(:, 1)- L/2) < tol)| ...
           (abs(NODE(:, 2) - 5*H/10) < tol & abs(NODE(:, 1)- L/2) < tol)| ...
           (abs(NODE(:, 2) - 6*H/10) < tol & abs(NODE(:, 1)- L/2) < tol)| ...
           (abs(NODE(:, 2) - 7*H/10) < tol & abs(NODE(:, 1)- L/2) < tol)| ...
           (abs(NODE(:, 2) - 8*H/10) < tol & abs(NODE(:, 1)- L/2) < tol)| ...
           (abs(NODE(:, 2) - 9*H/10) < tol & abs(NODE(:, 1)- L/2) < tol)| ...
           (abs(NODE(:, 2) - 10*H/10) < tol & abs(NODE(:, 1)- L/2) < tol));
DISP = cell(1, 2);
DISP{1} = cell(1, 10);
DISP{2} = cell(1, 10);
for i = 1:10
    hi = H/10;
    if i ==1
        DISP{1}{i} = [1 ind(i) 1/hi 0];
        DISP{2}{i} = [1 ind(i) 1/hi 0];
    else
        DISP{1}{i} = [1 ind(i) 1/hi 0;
                     1 ind(i-1) -1/hi 0];
        DISP{2}{i} = [1 ind(i) 1/hi 0;
                     1 ind(i-1) -1/hi 0];
    end

end

%% Passive elements 
Nfunc = @ShapeFuncQ4;
dNfunc = @ShapeDerivQ4;
intOrder = 1;
CM = GetElementCM2(NODE,ELEM,[],Nfunc,dNfunc,intOrder);
pass_solid = CM(:,2)  < H/10 + tol & CM(:,2)  > H/10 -d - tol | ...
             CM(:,2)  < 2*H/10 + tol & CM(:,2)  > 2*H/10 -d - tol | ...
             CM(:,2)  < 3*H/10 + tol & CM(:,2)  > 3*H/10 -d - tol | ...
             CM(:,2)  < 4*H/10 + tol & CM(:,2)  > 4*H/10 -d - tol | ...
             CM(:,2)  < 5*H/10 + tol & CM(:,2)  > 5*H/10 -d - tol | ...
             CM(:,2)  < 6*H/10 + tol & CM(:,2)  > 6*H/10 -d - tol | ...
             CM(:,2)  < 7*H/10 + tol & CM(:,2)  > 7*H/10 -d - tol | ...
             CM(:,2)  < 8*H/10 + tol & CM(:,2)  > 8*H/10 -d - tol | ...
             CM(:,2)  < 9*H/10 + tol & CM(:,2)  > 9*H/10 -d - tol | ...
             CM(:,2)  < 10*H/10 + tol & CM(:,2)  > 10*H/10 -d - tol |...
             CM(:,1) - 0 < d + tol | ...
             CM(:,1)  >L - d - tol;
pass_void = false(nx*ny, 1); 

%% Filter radius
rmin = 0.2;

end
