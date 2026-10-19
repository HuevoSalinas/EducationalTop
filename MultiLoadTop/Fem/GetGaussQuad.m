function [GP_N,GP_w,GP_xi] = GetGaussQuad(order)

% Gauss points and weights
switch order
    case 1
        GP_xi = [0 0]';
        GP_w  = 4;
    case 2
        GP_xi = [-1 -1;-1 1;1 -1; 1 1]' / sqrt(3);
        GP_w  = [ 1 1 1 1];
  otherwise
        error('GetGaussQuad: Rule order n=%g not available',order)
end
GP_N = length(GP_w);

end
