function [Dv,Ds] = GetDmatrixContinuum(E,nu,constLaw)

% D matrix
switch upper(constLaw)
    case 'PS'
        Dv = E/(1-nu^2)*[1  nu;
                         nu 1 ]; % volumetric component
        Ds = E/(2*(1+nu)); % shear component
    case 'PE'
        Dv = E/(1+nu)/(1-2*nu)*[1-nu  nu ;
                                 nu  1-nu]; % volumetric component
        Ds = E/(2*(1+nu)); % shear component
    otherwise
        error('Constitutive law is unknown')
end
end