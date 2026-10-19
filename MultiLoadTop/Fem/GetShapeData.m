function [values] = GetShapeData(func,GP_xi)
GP_N = size(GP_xi,2);
value = func(GP_xi(:,1));
values = zeros(size(value,1),size(value,2),GP_N);
values(:,:,1) = value;
for j=2:GP_N
    values(:,:,j) = func(GP_xi(:,j));
end
end
