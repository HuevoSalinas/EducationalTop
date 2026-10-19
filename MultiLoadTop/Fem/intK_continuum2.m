function [k,indI] = intK_continuum2(NODE,ELEM,MATprop,constLaw,dNfunc,intOrder)

% Initial parameters
DpN = 2; % Degrees of freedom per node
Ne = size(ELEM,1); % Number of elements
NpE = size(ELEM,2); % Number of nordes per element
edof = GetElemDOF(ELEM,DpN); % element degrees of fredom 

% Stiffness matrix integration
k = zeros(DpN*NpE,DpN*NpE,Ne);
indI = nan(DpN*NpE,DpN*NpE,Ne);
switch numel(intOrder)
    case 1 % equal integration for volume and shear 
        [GP_N,GP_w,GP_xi] = GetGaussQuad(intOrder);
        dN_GP = GetShapeData(dNfunc,GP_xi);
        for i=1:Ne
            XYZ = NODE(ELEM(i,:),:);
            [Dv,Ds] = GetDmatrixContinuum(MATprop.E(i),MATprop.nu(i),constLaw);
            for j=1:GP_N
                dN = dN_GP(:,:,j);
                [Bv,Bs,detJ] = GetBmatrixContinuum(dN,XYZ);
                k(:,:,i) = k(:,:,i) + GP_w(j) * (Bv'*Dv*Bv+Bs'*Ds*Bs) * MATprop.t(i) * detJ;
            end
            
            indI(:,:,i) = repmat(edof(i,:),size(indI,1),1);
        end
    case 2 % diferent integration for volume and shear
        [FGP_N,FGP_w,FGP_xi] = GetGaussQuad(intOrder(1));
        dN_FGP = GetShapeData(dNfunc,FGP_xi);
        [RGP_N,RGP_w,RGP_xi] = GetGaussQuad(intOrder(2));
        dN_RGP = GetShapeData(dNfunc,RGP_xi);
        for i=1:Ne
            XYZ = NODE(ELEM(i,:),:);
            [Dv,Ds] = GetDmatrixContinuum(MATprop.E(i),MATprop.nu(i),constLaw);
            for j=1:FGP_N % volume integration
                dN = dN_FGP(:,:,j);
                [Bv,~,detJ] = GetBmatrixContinuum(dN,XYZ);
                k(:,:,i) = k(:,:,i) + FGP_w(j) * (Bv'*Dv*Bv) * MATprop.t(i) * detJ;
            end
            for j=1:RGP_N % shear integration
                dN = dN_RGP(:,:,j);
                [~,Bs,detJ] = GetBmatrixContinuum(dN,XYZ);
                k(:,:,i) = k(:,:,i) + RGP_w(j) * (Bs'*Ds*Bs) * MATprop.t(i) * detJ;
            end
            indI(:,:,i) = repmat(edof(i,:),size(indI,1),1);
        end
    otherwise
        error('intOrder must have length 1 (full or reduced integration) or 2 (selective).')
end
end