function [fem] = FEMsolve(fem)

if ~isfield(fem,'solveMethod')
    fem.solveMethod = 'direct';
end

% Initial parameters
Ndof = size(fem.f,1); % number of degrees of freedom
Nf = size(fem.f, 2); % nomber of load cases 
fem.freeDOF = setdiff(1:Ndof,fem.suppDOF); % free dregrees of freedom index

% Get displacements 
fem.u = zeros(size(fem.f)); % displacement vector 
fem.u(fem.suppDOF,:) = repmat(fem.suppVAL, 1,Nf); % displacements in supports
b = fem.f(fem.freeDOF,:) ...
  - repmat(fem.K(fem.freeDOF,fem.suppDOF)*fem.suppVAL, 1,Nf); % forces produced by support displacements
switch lower(fem.solveMethod)
    case 'direct'
        fem.K = (fem.K+fem.K')/2; % Enforce the symmetry of the matrix.
        fem.u(fem.freeDOF,:) = fem.K(fem.freeDOF,fem.freeDOF) \ b;
    otherwise
        error('Linear solve method is not supported')
end

fem.r = fem.K(fem.suppDOF,:) * fem.u - fem.f(fem.suppDOF,:);
Nn = Ndof / fem.DpN;
fem.D = permute(reshape(fem.u,fem.DpN,Nn,Nf),[2,1,3]);

end