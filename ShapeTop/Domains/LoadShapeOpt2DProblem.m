function [NODE, ELEM, SUPP, LOAD, DISP, x, volfrac, pass_solid, pass_void, rmin] ...
    = LoadShapeOpt2DProblem(problemID,refine)
narginchk(1,2)

if (nargin<2 || isempty(refine))
    switch problemID
        case 1
            refine = 2;
        case 2
            refine = 4;
        case 3
            refine = 6;
        case 4
            refine = 4;
        case 5
            refine = 1/8;
        case 6
            refine = 2;
    end
end

% Load problem
switch problemID
    case 1
        [NODE, ELEM, SUPP, LOAD, DISP, x, pass_solid, pass_void, rmin] = ...
            LoadMBBholesProblem(refine);
    case 2
        [NODE, ELEM, SUPP, LOAD, DISP, x, pass_solid, pass_void, rmin] = ...
            LoadBeamHeightProblem(refine);
    case 3
        [NODE, ELEM, SUPP, LOAD, DISP, x, pass_solid, pass_void, rmin] = ...
            LoadMBBvierendeelProblem(refine);
    case 4
        [NODE, ELEM, SUPP, LOAD, DISP, x, pass_solid, pass_void, rmin] = ...
            LoadLShapedBeamProblem(refine);
    case 5
        [NODE, ELEM, SUPP, LOAD, DISP, x, pass_solid, pass_void, rmin] = ...
            LoadBracketProblem(refine);
    case 6
        [NODE, ELEM, SUPP, LOAD, DISP, x, pass_solid, pass_void, rmin] = ...
            LoadMBBFullProblem(refine);
    otherwise
        error('Unknown problem. Cannot load.')
end

if ~iscolumn(x)
    fprintf('problemID = %g is not returning x in a column: FIX IT!\n',problemID)
    x = x(:);
end

% Get volfrac
[~,V] = GetElementCM2(NODE,ELEM,[],@ShapeFuncQ4,@ShapeDerivQ4,2);
active = (~pass_solid & ~pass_void);
volfrac = sum(x(active).*V(active))/sum(V(active));

end

%% Problem ID 1
function [NODE, ELEM, SUPP, LOAD, DISP, x, pass_solid, pass_void, rmin] = ...
    LoadMBBholesProblem(refine)

Lx = 3;
Ly = 1;
B = 0.75;
Nx = round(75*refine);
Ny = round(25*refine);

[NODE,ELEM,~,~,CM] = BuildOrthoMesh_embedded(Lx,Ly,Nx,Ny);
Ne = size(ELEM,1);

SUPP = nan(Ny+2,3);
SUPP(1:end-1,1) = 1:Ny+1;
SUPP(1:end-1,2) = 0;
SUPP(end,1) = Nx*(Ny+1)+1;
SUPP(end,3) = 0;

LOAD = [Ny+1 0 -1];
DISP = [1 Ny+1 0 1];

% Cut the circular holes
x = ones(Ne,1);
for i=1:3
    C = [-0.5+i 0.5];
    dist = vecnorm(CM-repmat(C,Ne,1),2,2);
    bool_dist = dist<=B/2;
    x(bool_dist) = 0;
end

ind_Node = find(NODE(:,1)<0.1);
pass_solid = any(ismember(ELEM,ind_Node),2);
pass_void = false(Ne,1);
rmin = 0.05;
end

%% Problem ID 2
function [NODE, ELEM, SUPP, LOAD, DISP, x, pass_solid, pass_void, rmin] = ...
    LoadBeamHeightProblem(refine)

L = 1;
H = 0.2;
Nx = round(25*refine);
Ny = round(10*refine);

[NODE,ELEM,Xs,~,CM] = BuildOrthoMesh_embedded(L,H,Nx,Ny);
Nn = size(NODE,1);
Ne = size(ELEM,1);

% Compute the loads and supports
dx = Xs(2)-Xs(1);
f_ref = -1e6*((Xs/L).^4).*(1-Xs/L);
ind = (Ny+1):(Ny+1):Nn;
LOAD = [ind' zeros(Nx+1,1) dx*f_ref'];
LOAD([1 end],3) = -(dx/2)*f_ref([1 end]);

ind = 1:(Ny+1);
SUPP = [ind' zeros(Ny+1,2)]; 
DISP = [1 Nn 0 1];

% Apply the top solid slab
ind = find(CM(:,2)>0.9*H);
pass_solid = false(Ne,1);
pass_solid(ind) = true;
pass_void = false(Ne,1);
ind = find(CM(:,2)>0.3*H);
x = zeros(Ne,1);
x(ind) = 1;
rmin = 0.02;
end

%% Problem ID 3
function [NODE, ELEM, SUPP, LOAD, DISP, x, pass_solid, pass_void, rmin] = ...
    LoadMBBvierendeelProblem(refine)

Nx = round(30*refine);
Ny = round(5*refine);
Lx = Nx;
Ly = Ny;

[NODE,ELEM,~,~,CM] = BuildOrthoMesh_embedded(Lx,Ly,Nx,Ny);
Nn = size(NODE,1);
Ne = size(ELEM,1);

SUPP = [  1    0  0;
        Nn-Ny NaN 0];
LOAD = [(round(Nx/2)+1)*(Ny+1) 0 -1];
DISP = [1 (round(Nx/2)+1)*(Ny+1) 0 1];

% Cut the square holes
x = zeros(Ne,1);
ind = find(CM(:,2)<1*refine | CM(:,2)>4*refine);
x(ind) = 1;
for i=1:4
    hangX = 4.5+7*(i-1);
    ind = find(CM(:,1)>(hangX-0.5)*refine & CM(:,1)<(hangX+0.5)*refine);
    x(ind) = 1;
end
pass_solid = false(Ne,1);
pass_void = false(Ne,1);
rmin = 2;
end

%% Problem ID 4
function [NODE, ELEM, SUPP, LOAD, DISP, x, pass_solid, pass_void, rmin] = ...
    LoadLShapedBeamProblem(refine)

Nx = round(30*refine);
Ny = Nx;
Lx = 2;
Ly = 2;

[NODE,ELEM,~,~,CM] = BuildOrthoMesh_embedded(Lx,Ly,Nx,Ny);
Nn = size(NODE,1);
Ne = size(ELEM,1);

m = 1/35*Lx; % element half-width
% segments defined as [x1 x2 y1 y2 half-width]
S = [1/9*Lx 1/9*Lx 2/5*Lx-m    Lx    m;
     1/3*Lx 1/3*Lx   0+m       Lx    m;
     1/9*Lx 1/3*Lx 2/5*Lx-m   0+m    m;
     1/9*Lx   Lx   2/5*Lx-m 2/5*Lx-m m;
     1/3*Lx 4/5*Lx   0+m      0+m    m;
     4/5*Lx  Lx      0+m    2/5*Lx-m m;
     1/3*Lx 4/5*Lx 2/5*Lx-m   0+m    m];

% Determine the initial "x" by calculating the distance to all segments
D = Inf(Nn,1);
for i=1:size(S,1)
    thisD = DistanceToSegment(S(i,[1 3]),S(i,[2 4]),NODE);
    ind = find(thisD>S(i,5));
    thisD(ind) = NaN;
    D = min([D thisD],[],2);
end
aux = isfinite(D);
D = D - m;
D(~aux) = 0;
D = double(D>=0);

% Element density is the average of the nodes
x = nan(Ne,1);
for i=1:Ne
    x(i) = mean(D(ELEM(i,:)));
end
x = 1-x;

ind = (Ny+1):(Ny+1):Nn;
SUPP = [ind' zeros(length(ind),2)];

dist = vecnorm(NODE-repmat([Lx 2/5*Lx],Nn,1),2,2);
[~,ind] = min(dist);
LOAD = [ind 0 -1];
DISP = [1 ind 0 1];

pass_solid = false(Ne,1);
pass_void = false(Ne,1);
ind = find(CM(:,1)>2/5*Lx & CM(:,2)>2/5*Lx);
pass_void(ind) = true;
x(pass_void) = 0;
rmin = 0.05;
end

%% Problem ID 5
function [NODE, ELEM, SUPP, LOAD, DISP, x, pass_solid, pass_void, rmin] = ...
    LoadBracketProblem(refine)

cutoff = 0.5;
topload = [3 -1];

[NODE,ELEM,x,Nx,Ny] = ReadImageTOP('Plane_Seat_Bracket.png',refine);
Nn = size(NODE,1);
Ne = Nx*Ny;

ElemSolid = (x > cutoff);
ind = ELEM(ElemSolid,:);
NodeUse = accumarray(ind(:),1,[Nn 1]);
NodeUse = (NodeUse>0);

ind1 = find(NODE(:,2)==0 & NodeUse);
ind2 = find(NODE(:,2)==Ny & NodeUse);

SUPP = [ind1 zeros(length(ind1),2)];
LOAD = [ind2 repmat(topload,length(ind2),1)];
DISP = [(1:length(ind2))' ind2 repmat([1 0],length(ind2),1)];
pass_solid = false(Ne,1);
pass_void = false(Ne,1); 
rmin = 2;
end

%% Problem ID 6
function [NODE, ELEM, SUPP, LOAD, DISP, x, pass_solid, pass_void, rmin] = ...
    LoadMBBFullProblem(refine)

B = 0.75;
Nx = round(150*refine);
Ny = round(25*refine);
Lx = 6;
Ly = 1;

[NODE,ELEM,~,~,CM] = BuildOrthoMesh_embedded(Lx,Ly,Nx,Ny);
% NODE(:,1) = NODE(:,1) - Lx/2; % Center the domain in X
% CM(:,1) = CM(:,1) - Lx/2;
Nn = size(NODE,1);
Ne = size(ELEM,1);

SUPP = nan(3,3);
BM = round(Nx/2)*(Ny+1)+1; % bottom-middle node
BR = Nn-Ny; % bottom-right node
SUPP = [ 1 NaN  0 ;
        BM  0  NaN;
        BR NaN  0 ];
LOAD = [BM+Ny 0 -1];
DISP = [1 BM+Ny 0 1];

% Cut the circular holes
x = ones(Ne,1);
for i=1:6
    % C = [-3.5+i 0.5];
    C = [-0.5+i 0.5];
    dist = vecnorm(CM-repmat(C,Ne,1),2,2);
    ind = find(dist<=B/2);
    x(ind) = 0;
end

pass_solid = false(Ne,1);
pass_void = false(Ne,1);
rmin = 0.05;
end

%% Additional embedded functions
function [normW] = DistanceToSegment(A,B,X)

% Segment director
D = B - A;
L = norm(D);
D = D / L;

% Determine proyection to segment and distance to it
Np = size(X,1);
V = X - repmat(A,Np,1);
normP = V * D';
P = normP * D;
W = V - P;
normW = vecnorm(W,2,2);

% Distance to points not proyecting to segment is the distance to endpoints
indA = find(normP<0);
indB = find(normP>L);

normW(indA) = vecnorm(V(indA,:),2,2);
V2 = X(indB,:) - repmat(B,length(indB),1);
normW(indB) = vecnorm(V2,2,2);
end

function [NODE,ELEM,Xs,Ys,CM] = BuildOrthoMesh_embedded(Lx,Ly,Nx,Ny)

Xs = linspace(0,Lx,Nx+1);
Ys = linspace(0,Ly,Ny+1);
[X,Y] = meshgrid(Xs,Ys);
NODE = [X(:) Y(:)];
% Nn = size(NODE,1);
Ne = Nx*Ny;

ELEM = nan(Ne,4);
for i=1:Nx
    aux = (i-1)*(Ny+1) + [0 Ny+1 Ny+2 1];
    ELEM(Ny*(i-1)+(1:Ny),:) = repmat(aux,Ny,1) + repmat((1:Ny)',1,4);
end

if (nargout==5)
    CM = nan(Ne,2);
    for i=1:Ne
        CM(i,:) = mean(NODE(ELEM(i,:),:),1);
    end
end
end

function [NODE,ELEM,A,Nx,Ny] = ReadImageTOP(filename,scale)

% Read the file
A = imread(filename);
A = imresize(A,scale);
A = 1 - max(double(flipud(A)),[],3) / 255;
Nx = size(A,2);
Ny = size(A,1);
A = A(:);

[X,Y] = meshgrid(0:Nx,0:Ny);
NODE = [X(:) Y(:)];
% Nn = size(NODE,1);
Ne = Nx*Ny;

ELEM = nan(Ne,4);
for i=1:Nx
    aux = (i-1)*(Ny+1) + [0 Ny+1 Ny+2 1];
    ELEM(Ny*(i-1)+(1:Ny),:) = repmat(aux,Ny,1) + repmat((1:Ny)',1,4);
end
end