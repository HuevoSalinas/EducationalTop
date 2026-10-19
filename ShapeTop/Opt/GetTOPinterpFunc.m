function [E,dEdx] = GetTOPinterpFunc(matModel)

switch upper(matModel)
    case 'SIMP'
        E = @(x,penal)(x.^penal);
        dEdx = @(x,penal)(penal*x.^(penal-1));
    case 'RAMP'
        E = @(x,penal)(x./(1+penal*(1-x)));
        dEdx = @(x,penal)((penal+1)./(penal*(1-x)+1).^2);
    otherwise
        error('Material interpolation scheme not supported.')
end
end