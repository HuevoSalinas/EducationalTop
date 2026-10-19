function [Bv,Bs,detJ] = GetBmatrixContinuum(dN,XYZ)

% Derivates
J = dN * XYZ; % Jacoboian
detJ = det(J); % Jacobian determinant
dNdx = J \ dN; % Derivates in real coordinates

% B matrix
Bv = zeros(2,2*size(dNdx,2)); % volumetric component of B matrix 
Bv(1,1:2:end) = dNdx(1,:);
Bv(2,2:2:end) = dNdx(2,:);
Bs = zeros(1,2*size(dNdx,2)); % shear component of B matrix
Bs(1,1:2:end) = dNdx(2,:);
Bs(1,2:2:end) = dNdx(1,:);
  
return