function [H] = GetFilterH(CM,rmin,expS,active)
narginchk(2,4)

% Initial parameters 
Ne = size(CM,1); % Number of elements
if rmin <= 0
    H = speye(Ne);
    return
end
indI = cell(Ne,1);
indJ = cell(Ne,1);
indK = cell(Ne,1);

if (nargin<3 || isempty(expS))
    expS = 1;
end
if (nargin<4 || isempty(active))
    active = true(Ne,1);
end

% Distance between elements for filter
for i=1:Ne
    if active(i)
        dist = vecnorm(CM-repmat(CM(i,:),Ne,1),2,2);
        neighbors = find(dist<rmin & active);
        kernel = (1-dist(neighbors)/rmin).^expS;
        
        indI{i} = i*ones(length(neighbors),1);
        indJ{i} = neighbors;
        indK{i} = kernel/sum(kernel);
    else % Passive element
        indI{i} = i;
        indJ{i} = i;
        indK{i} = 1;
    end
end
H = sparse(vertcat(indI{:}),vertcat(indJ{:}),vertcat(indK{:}),Ne,Ne); % filter matrix 
end
