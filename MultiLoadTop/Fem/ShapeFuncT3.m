function [N] = ShapeFuncT3(xi)

N = zeros(1,3);
N(1)= 1-xi(1)-xi(2);
N(2)= xi(1);
N(3)= xi(2);
return