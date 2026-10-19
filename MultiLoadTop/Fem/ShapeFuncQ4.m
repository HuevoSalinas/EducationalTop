function [N] = ShapeFuncQ4(xi)

N=zeros(1,4);
N(1)= 1/4*(1-xi(1))*(1-xi(2));
N(2)= 1/4*(1+xi(1))*(1-xi(2));
N(3)= 1/4*(1+xi(1))*(1+xi(2));
N(4)= 1/4*(1-xi(1))*(1+xi(2));
return